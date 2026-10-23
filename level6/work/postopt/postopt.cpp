// Post-optimizer for level 6 programs: shortens valid programs by local edits that keep them valid.
// Edits tried repeatedly until no improvement:
//   * drop everything after the tick of the count-th acceptance
//   * delete one line
//   * merge two consecutive lines (each emitter has an instruction in at most one of them)
// Every edit is checked by exact simulation (same rules as sim6.py).
// Usage: postopt <input file> <candidate output file> <result file>
#include <bits/stdc++.h>
#include <omp.h>
using namespace std;

static const int DY[6] = {0, 1, 1, 0, -1, -1}, DX[6] = {2, 1, -1, -2, -1, 1};
constexpr int W = 7;                                   // 7*64 >= 397 interior cells (r <= 12)
typedef array<uint64_t, W> Bits;

struct Geo {
    int r, n;                                          // n = perimeter length 6r
    vector<vector<array<Bits, 1>>> ray;               // ray[p][e]: interior cells hit from perimeter p
    map<pair<int,int>, int> idx;                       // interior cell -> bit index
};

Geo makeGeo(int r) {
    Geo g; g.r = r;
    set<pair<int,int>> interior, border;
    for (int i = 0; i <= 2 * r; i++) {
        int d = abs(i - r), m = 2 * r - d;
        for (int j = 0; j <= m; j++) ((d == r || j == 0 || j == m) ? border : interior).insert({i, d + 2 * j});
    }
    int k = 0; for (auto& c : interior) g.idx[c] = k++;
    vector<pair<int,int>> per;
    for (int q = 0; q < r; q++) per.push_back({r - q, q});
    for (int q = 0; q < r; q++) per.push_back({0, r + 2 * q});
    for (int q = 0; q < r; q++) per.push_back({q, 3 * r + q});
    for (int q = 0; q < r; q++) per.push_back({r + q, 4 * r - q});
    for (int q = 0; q < r; q++) per.push_back({2 * r, 3 * r - 2 * q});
    for (int q = 0; q < r; q++) per.push_back({2 * r - q, r - q});
    g.n = per.size();
    g.ray.assign(g.n, vector<array<Bits, 1>>(6));
    for (int p = 0; p < g.n; p++) for (int e = 0; e < 6; e++) {
        Bits b{}; int y = per[p].first + DY[e], x = per[p].second + DX[e];
        while (interior.count({y, x}) || border.count({y, x})) {
            auto it = g.idx.find({y, x});
            if (it != g.idx.end()) b[it->second >> 6] |= 1ULL << (it->second & 63);
            y += DY[e]; x += DX[e];
        }
        g.ray[p][e][0] = b;
    }
    return g;
}

struct Case { int r, cnt; Bits P{}, K{}; };

struct State {
    int pos[6]; bool act[6]; Bits m{}; int acc = 0;
    bool operator==(const State& o) const {
        for (int e = 0; e < 6; e++) if (pos[e] != o.pos[e] || act[e] != o.act[e]) return false;
        return m == o.m && acc == o.acc;
    }
};

typedef array<char, 6> Line;

State initState(const Geo& g) { State s; for (int e = 0; e < 6; e++) { s.pos[e] = e * g.r; s.act[e] = false; } return s; }

inline void step(const Geo& g, const Case& c, State& s, const Line& L) {
    for (int e = 0; e < 6; e++) {
        if (L[e] == 'F') s.pos[e] = (s.pos[e] + 1) % g.n;
        else if (L[e] == 'B') s.pos[e] = (s.pos[e] + g.n - 1) % g.n;
        else if (L[e] == 'T') s.act[e] = !s.act[e];
    }
    Bits one{}, two{}, three{};                       // bit-sliced saturating counter (>=3 sticky)
    for (int e = 0; e < 6; e++) if (s.act[e]) {
        const Bits& b = g.ray[s.pos[e]][e][0];
        for (int w = 0; w < W; w++) {
            uint64_t x = b[w];
            uint64_t c1 = one[w] & x;
            one[w] ^= x;
            uint64_t c2 = two[w] & c1;
            two[w] ^= c1;
            three[w] |= c2;
        }
    }
    bool ok = true;
    for (int w = 0; w < W; w++) {
        uint64_t exact2 = two[w] & ~one[w] & ~three[w];
        uint64_t ge3 = three[w] | (two[w] & one[w]);
        s.m[w] = (s.m[w] | exact2) & ~ge3;
        if ((s.m[w] & c.P[w]) != c.P[w] || (s.m[w] & c.K[w])) ok = false;
    }
    if (ok) { for (int w = 0; w < W; w++) s.m[w] &= ~c.P[w]; s.acc++; }
}

// simulate prog from state s starting at index from; returns tick count at which acc reaches cnt (+1), or -1
int finishTick(const Geo& g, const Case& c, State s, const vector<Line>& prog, int from) {
    if (s.acc >= c.cnt) return from;
    for (int i = from; i < (int)prog.size(); i++) { step(g, c, s, prog[i]); if (s.acc >= c.cnt) return i + 1; }
    return -1;
}

vector<Line> optimize(const Geo& g, const Case& c, vector<Line> prog) {
    // prefix states: st[i] = state before line i
    auto prefix = [&](const vector<Line>& p) { vector<State> st(p.size() + 1); st[0] = initState(g);
        for (size_t i = 0; i < p.size(); i++) { st[i + 1] = st[i]; step(g, c, st[i + 1], p[i]); } return st; };
    {   int t = finishTick(g, c, initState(g), prog, 0); if (t < 0) return prog; prog.resize(t); }
    bool improved = true;
    while (improved) {
        improved = false;
        vector<State> st = prefix(prog);
        // merges first (cheap and common), from the end backwards
        for (int i = (int)prog.size() - 2; i >= 0; i--) {
            if (i + 1 >= (int)prog.size()) continue;
            Line a = prog[i], b = prog[i + 1], m; bool comp = true;
            for (int e = 0; e < 6; e++) { if (a[e] != '-' && b[e] != '-') { comp = false; break; } m[e] = a[e] != '-' ? a[e] : b[e]; }
            if (!comp) continue;
            vector<Line> q(prog.begin(), prog.begin() + i); q.push_back(m); q.insert(q.end(), prog.begin() + i + 2, prog.end());
            State s = st[i]; step(g, c, s, m);
            int t;
            if (s == st[i + 2]) t = (int)prog.size() - 1;     // resynchronized: rest is identical
            else t = finishTick(g, c, s, q, i + 1);
            if (t >= 0) { q.resize(t); prog = q; st = prefix(prog); improved = true; }
        }
        for (int i = (int)prog.size() - 1; i >= 0; i--) {
            if (i >= (int)prog.size()) continue;
            vector<Line> q(prog.begin(), prog.begin() + i); q.insert(q.end(), prog.begin() + i + 1, prog.end());
            int t = finishTick(g, c, st[i], q, i);
            if (t >= 0) { q.resize(t); prog = q; st = prefix(prog); improved = true; }
        }
    }
    return prog;
}

int main(int argc, char** argv) {
    if (argc < 4) { fprintf(stderr, "usage: postopt in cand out\n"); return 1; }
    ifstream in(argv[1]); int n; in >> n;
    vector<Case> cases(n); map<int, Geo> geos;
    for (auto& c : cases) {
        int rows, cols; in >> c.r >> rows >> cols >> c.cnt;
        vector<string> gr(rows); for (auto& s : gr) in >> s;
        if (!geos.count(c.r)) geos[c.r] = makeGeo(c.r);
        Geo& g = geos[c.r];
        set<pair<int,int>> P;
        for (int y = 0; y < rows; y++) for (int x = 0; x < cols; x++) if (gr[y][x] == 'M') P.insert({y, x});
        for (auto& p : P) { int b = g.idx.at(p); c.P[b >> 6] |= 1ULL << (b & 63); }
        int ny[6] = {0, 0, 1, 1, -1, -1}, nx[6] = {2, -2, 1, -1, 1, -1};
        for (auto& p : P) for (int d = 0; d < 6; d++) {
            pair<int,int> q = {p.first + ny[d], p.second + nx[d]};
            auto it = g.idx.find(q);
            if (it != g.idx.end() && !P.count(q)) c.K[it->second >> 6] |= 1ULL << (it->second & 63);
        }
    }
    vector<vector<Line>> progs; { ifstream cf(argv[2]); string line; vector<Line> cur;
        while (getline(cf, line)) {
            if (!line.empty() && line.back() == '\r') line.pop_back();
            bool blank = line.find_first_not_of(" \t") == string::npos;
            if (blank) { if (!cur.empty()) { progs.push_back(cur); cur.clear(); } continue; }
            istringstream ss(line); Line L; string tok; int k = 0;
            while (ss >> tok && k < 6) L[k++] = tok[0];
            while (k < 6) L[k++] = '-';
            cur.push_back(L);
        }
        if (!cur.empty()) progs.push_back(cur); }
    if ((int)progs.size() != n) { fprintf(stderr, "expected %d programs, got %zu\n", n, progs.size()); return 1; }
    long long before = 0, after = 0;
    vector<vector<Line>> res(n);
    #pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < n; i++) res[i] = optimize(geos.at(cases[i].r), cases[i], progs[i]);
    ofstream out(argv[3], ios::binary);
    for (int i = 0; i < n; i++) {
        before += progs[i].size(); after += res[i].size();
        if (finishTick(geos.at(cases[i].r), cases[i], initState(geos.at(cases[i].r)), res[i], 0) < 0) res[i] = progs[i];  // safety
        if (i) out << "\n";
        for (auto& L : res[i]) { for (int e = 0; e < 6; e++) { if (e) out << ' '; out << L[e]; } out << "\n"; }
    }
    fprintf(stderr, "lines: %lld -> %lld\n", before, after);
}
