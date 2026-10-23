// Common geometry / parsing for level 6 cycle search.
#pragma once
#include <bits/stdc++.h>
using namespace std;

static const int AXIS[6] = {0, 1, 2, 0, 1, 2};      // 0=H (h=v-u), 1=U (u), 2=V (v)
static const int SGN[6] = {-1, +1, +1, +1, -1, -1};  // line = SGN[e]*k

struct Pattern {
    int r = 0, cnt = 0, rows = 0, cols = 0;
    vector<string> g;
    vector<array<int, 2>> cells;  // (u,v): P cells first, then K cells
    int nP = 0, nK = 0, N = 0;
    vector<array<int, 3>> lineOf;  // per cell: line index (c + r - 1) for axes H,U,V
    int lb1 = 0;                   // 1 + max over P of min(|u|,|v|,|h|)
};

static inline bool interior(int r, int u, int v) {
    return abs(u) <= r - 1 && abs(v) <= r - 1 && abs(v - u) <= r - 1;
}

static vector<Pattern> readInput(const string& path) {
    ifstream in(path);
    int n;
    in >> n;
    vector<Pattern> res(n);
    for (int ci = 0; ci < n; ci++) {
        Pattern& p = res[ci];
        in >> p.r >> p.rows >> p.cols >> p.cnt;
        p.g.resize(p.rows);
        for (auto& s : p.g) in >> s;
        int r = p.r;
        set<pair<int, int>> P, K;
        for (int y = 0; y < p.rows; y++)
            for (int x = 0; x < p.cols; x++)
                if (p.g[y][x] == 'M') P.insert({(x - y - r) / 2, (x + y - 3 * r) / 2});
        int du[6] = {1, -1, 0, 0, 1, -1}, dv[6] = {1, -1, 1, -1, 0, 0};
        for (auto [u, v] : P)
            for (int d = 0; d < 6; d++) {
                int a = u + du[d], b = v + dv[d];
                if (interior(r, a, b) && !P.count({a, b})) K.insert({a, b});
            }
        for (auto [u, v] : P) p.cells.push_back({u, v});
        for (auto [u, v] : K) p.cells.push_back({u, v});
        p.nP = P.size();
        p.nK = K.size();
        p.N = p.nP + p.nK;
        p.lb1 = 0;
        for (auto& c : p.cells) {
            int u = c[0], v = c[1];
            p.lineOf.push_back({v - u + r - 1, u + r - 1, v + r - 1});
        }
        for (int i = 0; i < p.nP; i++) {
            int u = p.cells[i][0], v = p.cells[i][1];
            p.lb1 = max(p.lb1, 1 + min({abs(u), abs(v), abs(v - u)}));
        }
    }
    return res;
}

// Raw emitter state: offset k in [-r, r] (|k|==r is the dead zone), active flag.
struct ES {
    int8_t k;
    int8_t a;
    bool operator==(const ES& o) const { return k == o.k && a == o.a; }
};
using Cfg = array<ES, 6>;

static inline bool adjES(const ES& x, const ES& y) {
    if (x.a == y.a) return abs(x.k - y.k) <= 1;
    return x.k == y.k;
}
static inline bool adjCfg(const Cfg& x, const Cfg& y) {
    for (int e = 0; e < 6; e++)
        if (!adjES(x[e], y[e])) return false;
    return true;
}
static inline char instr(int e, const ES& x, const ES& y) {
    (void)e;
    if (x.a != y.a) return 'T';
    if (y.k == x.k + 1) return 'F';
    if (y.k == x.k - 1) return 'B';
    return '-';
}
// effective line index (c + r - 1) or -1
static inline int effLine(int r, int e, const ES& s) {
    if (!s.a || abs(s.k) >= r) return -1;
    return SGN[e] * s.k + r - 1;
}

// Convert a config sequence (starting from the initial config, not included) into program lines.
static vector<string> toProgram(const vector<Cfg>& seq) {
    Cfg cur;
    for (int e = 0; e < 6; e++) cur[e] = {0, 0};
    vector<string> out;
    for (auto& c : seq) {
        if (!adjCfg(cur, c)) {
            fprintf(stderr, "toProgram: non-adjacent configs\n");
            exit(1);
        }
        string s;
        for (int e = 0; e < 6; e++) {
            if (e) s += ' ';
            s += instr(e, cur[e], c[e]);
        }
        out.push_back(s);
        cur = c;
    }
    return out;
}

// Exact simulation over P u K cells (cells beyond are free and irrelevant). Returns #accepts.
struct Sim {
    const Pattern* p;
    int r;
    vector<vector<int>> lineCells[3];  // per axis, per line index: cell ids
    Sim(const Pattern& pat) : p(&pat), r(pat.r) {
        for (int a = 0; a < 3; a++) lineCells[a].assign(2 * r - 1, {});
        for (int i = 0; i < pat.N; i++)
            for (int a = 0; a < 3; a++) lineCells[a][pat.lineOf[i][a]].push_back(i);
    }
    // accept ticks (1-based) for the program given as config sequence
    vector<int> run(const vector<Cfg>& seq) const {
        vector<char> m(p->N, 0);
        vector<int> acc;
        vector<int> hit(p->N);
        for (size_t t = 0; t < seq.size(); t++) {
            fill(hit.begin(), hit.end(), 0);
            for (int e = 0; e < 6; e++) {
                int li = effLine(r, e, seq[t][e]);
                if (li < 0) continue;
                for (int c : lineCells[AXIS[e]][li]) hit[c]++;
            }
            for (int i = 0; i < p->N; i++) {
                if (hit[i] == 2) m[i] = 1;
                else if (hit[i] >= 3) m[i] = 0;
            }
            bool ok = true;
            for (int i = 0; i < p->nP && ok; i++) ok = m[i];
            for (int i = p->nP; i < p->N && ok; i++) ok = !m[i];
            if (ok) {
                acc.push_back((int)t + 1);
                for (int i = 0; i < p->nP; i++) m[i] = 0;
            }
        }
        return acc;
    }
};
