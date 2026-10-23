// shots: backward beam planner for CCC level 6.
// Backward view: going from the acceptance tick backwards, every tick "decides" the cells it touches
// (2 hits -> matter, 3+ -> empty). A build is valid iff no undecided K cell is ever created and no
// undecided P cell is ever destroyed, and all P cells are decided (created) by the build's first tick.
#include <bits/stdc++.h>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace std;
typedef unsigned long long u64;

static const int DY[6] = {0, 1, 1, 0, -1, -1}, DX[6] = {2, 1, -1, -2, -1, 1};

struct Geo {
    int r, n;
    vector<pair<int,int>> per;
    set<pair<int,int>> interior;
    vector<array<vector<pair<int,int>>,6>> ray;
};
Geo makeGeo(int r) {
    Geo g; g.r = r; set<pair<int,int>> border;
    for (int i = 0; i <= 2 * r; i++) { int d = abs(i - r), m = 2 * r - d;
        for (int j = 0; j <= m; j++) { pair<int,int> c = {i, d + 2 * j};
            if (d == r || j == 0 || j == m) border.insert(c); else g.interior.insert(c); } }
    for (int k = 0; k < r; k++) g.per.push_back({r - k, k});
    for (int k = 0; k < r; k++) g.per.push_back({0, r + 2 * k});
    for (int k = 0; k < r; k++) g.per.push_back({k, 3 * r + k});
    for (int k = 0; k < r; k++) g.per.push_back({r + k, 4 * r - k});
    for (int k = 0; k < r; k++) g.per.push_back({2 * r, 3 * r - 2 * k});
    for (int k = 0; k < r; k++) g.per.push_back({2 * r - k, r - k});
    g.n = 6 * r;
    for (auto [y, x] : g.per) { array<vector<pair<int,int>>,6> row;
        for (int e = 0; e < 6; e++) { int cy = y + DY[e], cx = x + DX[e];
            while (g.interior.count({cy, cx}) || border.count({cy, cx})) {
                if (g.interior.count({cy, cx})) row[e].push_back({cy, cx});
                cy += DY[e]; cx += DX[e]; } }
        g.ray.push_back(row); }
    return g;
}

struct Pat { int r, cnt; vector<string> g; };
vector<Pat> parseInput(const string& path) {
    ifstream in(path); int n; in >> n; vector<Pat> v(n);
    for (auto& p : v) { int rows, cols; in >> p.r >> rows >> cols >> p.cnt; p.g.resize(rows); for (auto& s : p.g) in >> s; }
    return v;
}

// ---------------- exact forward simulator (mirror of sim6.py) ----------------
// returns tick index (1-based count of lines) at which the cnt-th acceptance happens, or -1
int simulate(const Geo& G, const Pat& P, const vector<string>& prog, int need, int* accOut = nullptr) {
    int r = G.r, n = G.n;
    map<pair<int,int>, int> id; int ni = 0; for (auto c : G.interior) id[c] = ni++;
    vector<char> isP(ni, 0), isK(ni, 0);
    int ny[6] = {0, 0, 1, 1, -1, -1}, nx[6] = {2, -2, 1, -1, 1, -1};
    for (int y = 0; y < (int)P.g.size(); y++) for (int x = 0; x < (int)P.g[y].size(); x++) if (P.g[y][x] == 'M') isP[id.at({y, x})] = 1;
    for (int y = 0; y < (int)P.g.size(); y++) for (int x = 0; x < (int)P.g[y].size(); x++) if (P.g[y][x] == 'M')
        for (int d = 0; d < 6; d++) { auto it = id.find({y + ny[d], x + nx[d]}); if (it != id.end() && !isP[it->second]) isK[it->second] = 1; }
    vector<vector<int>> rayId(n * 6);
    for (int p = 0; p < n; p++) for (int e = 0; e < 6; e++) for (auto c : G.ray[p][e]) rayId[p * 6 + e].push_back(id.at(c));
    int pos[6], act[6] = {0, 0, 0, 0, 0, 0}; for (int e = 0; e < 6; e++) pos[e] = e * r;
    vector<char> matter(ni, 0); vector<int> hit(ni, 0);
    int acc = 0, t = 0;
    for (auto& line : prog) {
        t++;
        for (int e = 0; e < 6; e++) { char c = line[e];
            if (c == 'F') pos[e] = (pos[e] + 1) % n; else if (c == 'B') pos[e] = (pos[e] - 1 + n) % n; else if (c == 'T') act[e] ^= 1; }
        fill(hit.begin(), hit.end(), 0);
        for (int e = 0; e < 6; e++) if (act[e]) for (int c : rayId[pos[e] * 6 + e]) hit[c]++;
        for (int i = 0; i < ni; i++) { if (hit[i] == 2) matter[i] = 1; else if (hit[i] >= 3) matter[i] = 0; }
        bool ok = true;
        for (int i = 0; i < ni && ok; i++) { if (isP[i] && !matter[i]) ok = false; if (isK[i] && matter[i]) ok = false; }
        if (ok) { for (int i = 0; i < ni; i++) if (isP[i]) matter[i] = 0; acc++; if (acc >= need) { if (accOut) *accOut = acc; return t; } }
    }
    if (accOut) *accOut = acc;
    return -1;
}

// ---------------- configuration ----------------
struct Cfg { int8_t k[6]; uint8_t a; };
static inline bool cfgEq(const Cfg& x, const Cfg& y) { for (int e = 0; e < 6; e++) if (x.k[e] != y.k[e]) return false; return x.a == y.a; }
static const int SGN[6] = {-1, 1, 1, 1, -1, -1};

static inline u64 mix64(u64 x) { x += 0x9e3779b97f4a7c15ULL; x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL; x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL; return x ^ (x >> 31); }

string stepInstr(const Cfg& A, const Cfg& B) {
    string s(6, '-');
    for (int e = 0; e < 6; e++) {
        int da = (A.a >> e) & 1, db = (B.a >> e) & 1;
        if (A.k[e] == B.k[e] && da == db) s[e] = '-';
        else if (A.k[e] == B.k[e]) s[e] = 'T';
        else if (B.k[e] == A.k[e] + 1 && da == db) s[e] = 'F';
        else if (B.k[e] == A.k[e] - 1 && da == db) s[e] = 'B';
        else { fprintf(stderr, "bad step\n"); exit(1); }
    }
    return s;
}

struct Params { int beam = 400; int finW = 100; int maxDepth = 400; int rootN = 400; u64 seed = 1; int maxAct = 1; double closW = 1.0; };

template<int NW>
struct Solver {
    typedef array<u64, NW> BS;
    int r, R2, np, nk, ncell, NS;
    BS Pm, Km;
    vector<BS> mask;                  // mask[e*R2 + k+r]
    vector<int> noffT, nliveT; vector<char> directT;
    Cfg I;

    static inline bool zero(const BS& a) { for (int i = 0; i < NW; i++) if (a[i]) return false; return true; }
    static inline int pc(const BS& a) { int s = 0; for (int i = 0; i < NW; i++) s += __builtin_popcountll(a[i]); return s; }

    inline int sidx(const Cfg& c, int e) const { return (c.k[e] + r) * 2 + ((c.a >> e) & 1); }
    inline bool offState(int s) const { int k = s / 2 - r, a = s & 1; return !a || abs(k) >= r; }

    void init(int r_, const vector<pair<int,int>>& Pc, const vector<pair<int,int>>& Kc, const Geo& G) {
        r = r_; R2 = 2 * r + 1; np = Pc.size(); nk = Kc.size(); ncell = np + nk; NS = 2 * R2;
        map<pair<int,int>, int> id; for (int i = 0; i < np; i++) id[Pc[i]] = i; for (int i = 0; i < nk; i++) id[Kc[i]] = np + i;
        Pm.fill(0); Km.fill(0);
        for (int i = 0; i < np; i++) Pm[i >> 6] |= 1ULL << (i & 63);
        for (int i = np; i < ncell; i++) Km[i >> 6] |= 1ULL << (i & 63);
        mask.assign(6 * R2, BS{});
        for (int e = 0; e < 6; e++) for (int k = -r; k <= r; k++) {
            int p = ((e * r + k) % G.n + G.n) % G.n; BS m{}; m.fill(0);
            for (auto c : G.ray[p][e]) { auto it = id.find(c); if (it != id.end()) m[it->second >> 6] |= 1ULL << (it->second & 63); }
            mask[e * R2 + k + r] = m; }
        for (int e = 0; e < 6; e++) I.k[e] = 0; I.a = 0;
        // closure tables
        noffT.assign(NS * NS, 1 << 20); nliveT.assign(NS * NS, 1 << 20); directT.assign(NS * NS, 0);
        auto nbrs = [&](int s) { vector<int> v; int k = s / 2 - r, a = s & 1; v.push_back(s);
            if (k + 1 <= r) v.push_back((k + 1 + r) * 2 + a); if (k - 1 >= -r) v.push_back((k - 1 + r) * 2 + a); v.push_back(s ^ 1); return v; };
        // dist (any states)
        vector<vector<int>> dist(NS, vector<int>(NS, 1 << 20));
        for (int s = 0; s < NS; s++) { deque<int> q; dist[s][s] = 0; q.push_back(s);
            while (!q.empty()) { int x = q.front(); q.pop_front(); for (int y : nbrs(x)) if (dist[s][y] > dist[s][x] + 1) { dist[s][y] = dist[s][x] + 1; q.push_back(y); } } }
        // off-only dist
        vector<vector<int>> od(NS, vector<int>(NS, 1 << 20));
        for (int s = 0; s < NS; s++) if (offState(s)) { deque<int> q; od[s][s] = 0; q.push_back(s);
            while (!q.empty()) { int x = q.front(); q.pop_front(); for (int y : nbrs(x)) if (offState(y) && od[s][y] > od[s][x] + 1) { od[s][y] = od[s][x] + 1; q.push_back(y); } } }
        for (int s = 0; s < NS; s++) for (int t = 0; t < NS; t++) {
            auto ns = nbrs(s); directT[s * NS + t] = find(ns.begin(), ns.end(), t) != ns.end();
            nliveT[s * NS + t] = max(1, dist[s][t] - 1);
            int best = 1 << 20;
            for (int i1 : nbrs(s)) if (offState(i1)) for (int in : nbrs(t)) if (offState(in)) best = min(best, 1 + od[i1][in]);
            noffT[s * NS + t] = best;
        }
    }

    // effect of a config
    void effect(const Cfg& c, BS& S, BS& D) const {
        BS c1{}, c2{}, c3{}; c1.fill(0); c2.fill(0); c3.fill(0);
        for (int e = 0; e < 6; e++) if ((c.a >> e) & 1) { const BS& m = mask[e * R2 + c.k[e] + r];
            for (int i = 0; i < NW; i++) { u64 n3 = c3[i] | (c2[i] & m[i]); u64 n2 = (c2[i] & ~m[i]) | (c1[i] & m[i]); u64 n1 = (c1[i] & ~m[i]) | (~(c1[i] | c2[i] | c3[i]) & m[i]); c3[i] = n3; c2[i] = n2; c1[i] = n1; } }
        S = c2; D = c3;
    }

    // closure: number of transit ticks strictly between X and C (forward X -> C); liveE = emitter allowed live
    int closure(const Cfg& X, const Cfg& C, int* liveE = nullptr) const {
        bool dir = true; int m1 = -1, m2 = -1, i1 = -1;
        for (int e = 0; e < 6; e++) { int s = sidx(X, e), t = sidx(C, e);
            if (!directT[s * NS + t]) dir = false;
            int v = noffT[s * NS + t]; if (v > m1) { m2 = m1; m1 = v; i1 = e; } else if (v > m2) m2 = v; }
        if (liveE) *liveE = -1;
        if (dir) return 0;
        int best = m1;
        int alt = max(nliveT[sidx(X, i1) * NS + sidx(C, i1)], m2);
        if (alt < best) { best = alt; if (liveE) *liveE = i1; }
        return max(best, 1);
    }

    // build explicit transit configs (n intermediates), returns list of configs strictly between X and C
    vector<Cfg> transit(const Cfg& X, const Cfg& C, int n, int liveE) const {
        vector<Cfg> out(n, X);
        for (int e = 0; e < 6; e++) {
            int s = sidx(X, e), t = sidx(C, e); bool live = (e == liveE);
            // time-expanded BFS: layers 0..n+1
            vector<vector<int>> par(n + 2, vector<int>(NS, -2));
            par[0][s] = -1;
            for (int L = 0; L <= n; L++) for (int x = 0; x < NS; x++) if (par[L][x] != -2) {
                int k = x / 2 - r, a = x & 1; int cand[4] = {x, (k + 1 <= r) ? (k + 1 + r) * 2 + a : -1, (k - 1 >= -r) ? (k - 1 + r) * 2 + a : -1, x ^ 1};
                for (int y : cand) { if (y < 0) continue; if (L + 1 <= n && !live && !offState(y)) continue; if (L + 1 == n + 1 && y != t) continue;
                    if (par[L + 1][y] == -2) par[L + 1][y] = x; } }
            if (par[n + 1][t] == -2) { fprintf(stderr, "transit fail\n"); exit(1); }
            int x = t;
            for (int L = n + 1; L >= 1; L--) { int px = par[L][x]; if (L - 1 >= 1) { out[L - 2].k[e] = px / 2 - r; if (px & 1) out[L - 2].a |= 1 << e; else out[L - 2].a &= ~(1 << e); } x = px; }
        }
        return out;
    }

