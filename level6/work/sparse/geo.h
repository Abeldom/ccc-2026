// Common infrastructure: geometry (exact port of sim6.py), bitsets, parsing, exact simulation.
#pragma once
#include <bits/stdc++.h>
#include <cassert>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace std;
typedef unsigned long long u64;

constexpr int NW = 7;  // 448 bits >= 397 cells (r=12)
struct Bits {
    u64 w[NW];
};
static inline Bits bz() { Bits b; memset(b.w, 0, sizeof b.w); return b; }
static inline void bset(Bits& b, int i) { b.w[i >> 6] |= 1ULL << (i & 63); }
static inline void bclr(Bits& b, int i) { b.w[i >> 6] &= ~(1ULL << (i & 63)); }
static inline bool bget(const Bits& b, int i) { return (b.w[i >> 6] >> (i & 63)) & 1; }
static inline Bits operator&(const Bits& a, const Bits& b) { Bits c; for (int i = 0; i < NW; i++) c.w[i] = a.w[i] & b.w[i]; return c; }
static inline Bits operator|(const Bits& a, const Bits& b) { Bits c; for (int i = 0; i < NW; i++) c.w[i] = a.w[i] | b.w[i]; return c; }
static inline Bits operator^(const Bits& a, const Bits& b) { Bits c; for (int i = 0; i < NW; i++) c.w[i] = a.w[i] ^ b.w[i]; return c; }
static inline Bits operator~(const Bits& a) { Bits c; for (int i = 0; i < NW; i++) c.w[i] = ~a.w[i]; return c; }
static inline Bits andnot(const Bits& a, const Bits& b) { Bits c; for (int i = 0; i < NW; i++) c.w[i] = a.w[i] & ~b.w[i]; return c; }
static inline int bpop(const Bits& a) { int s = 0; for (int i = 0; i < NW; i++) s += __builtin_popcountll(a.w[i]); return s; }
static inline bool bnone(const Bits& a) { u64 s = 0; for (int i = 0; i < NW; i++) s |= a.w[i]; return s == 0; }
static inline bool beq(const Bits& a, const Bits& b) { for (int i = 0; i < NW; i++) if (a.w[i] != b.w[i]) return false; return true; }
template <class F> static inline void bforeach(const Bits& a, F f) {
    for (int i = 0; i < NW; i++) { u64 x = a.w[i]; while (x) { int b = __builtin_ctzll(x); f(i * 64 + b); x &= x - 1; } }
}

static const int AX[6] = {0, 1, 2, 0, 1, 2};   // H U V H U V
static const int SG[6] = {-1, 1, 1, 1, -1, -1};

struct Geo {
    int r = 0, n6 = 0, nc = 0;
    vector<int> cy, cx, L[3];         // per cell: char coords and lines (H,U,V)
    vector<array<Bits, 6>> ray;       // [perimeter index][emitter]
    vector<Bits> line[3];             // [axis][c + r], c in [-r, r]; |c|=r empty
    Bits all;
    map<pair<int, int>, int> idx;
    int cellAt(int h, int u) const {  // returns -1 if not interior; v = h+u
        auto it = lid.find({h, u});
        return it == lid.end() ? -1 : it->second;
    }
    map<pair<int, int>, int> lid;     // (h,u) -> cell
};

static Geo makeGeo(int r) {
    Geo g; g.r = r; g.n6 = 6 * r;
    set<pair<int, int>> interior, border;
    for (int i = 0; i <= 2 * r; i++) {
        int d = abs(i - r), m = 2 * r - d;
        for (int j = 0; j <= m; j++) {
            pair<int, int> c = {i, d + 2 * j};
            if (d == r || j == 0 || j == m) border.insert(c); else interior.insert(c);
        }
    }
    vector<pair<int, int>> per;
    for (int k = 0; k < r; k++) per.push_back({r - k, k});
    for (int k = 0; k < r; k++) per.push_back({0, r + 2 * k});
    for (int k = 0; k < r; k++) per.push_back({k, 3 * r + k});
    for (int k = 0; k < r; k++) per.push_back({r + k, 4 * r - k});
    for (int k = 0; k < r; k++) per.push_back({2 * r, 3 * r - 2 * k});
    for (int k = 0; k < r; k++) per.push_back({2 * r - k, r - k});
    assert((int)per.size() == 6 * r);
    for (auto& c : interior) {
        int id = g.nc++;
        g.idx[c] = id; g.cy.push_back(c.first); g.cx.push_back(c.second);
        int y = c.first, x = c.second;
        int u = (x - y - r) / 2, v = (x + y - 3 * r) / 2, h = y - r;
        g.L[0].push_back(h); g.L[1].push_back(u); g.L[2].push_back(v);
        g.lid[{h, u}] = id;
    }
    assert(g.nc <= NW * 64);
    g.all = bz(); for (int i = 0; i < g.nc; i++) bset(g.all, i);
    set<pair<int, int>> allc = interior; for (auto& c : border) allc.insert(c);
    const int DY[6] = {0, 1, 1, 0, -1, -1}, DX[6] = {2, 1, -1, -2, -1, 1};
    g.ray.resize(per.size());
    for (size_t p = 0; p < per.size(); p++)
        for (int e = 0; e < 6; e++) {
            Bits b = bz();
            int cy = per[p].first + DY[e], cx = per[p].second + DX[e];
            while (allc.count({cy, cx})) {
                auto it = g.idx.find({cy, cx});
                if (it != g.idx.end()) bset(b, it->second);
                cy += DY[e]; cx += DX[e];
            }
            g.ray[p][e] = b;
        }
    for (int a = 0; a < 3; a++) {
        g.line[a].assign(2 * r + 1, bz());
        for (int i = 0; i < g.nc; i++) {
            int c = g.L[a][i];
            if (abs(c) <= r - 1) bset(g.line[a][c + r], i);
        }
    }
    // verify line formulas against ray tracing
    for (int e = 0; e < 6; e++)
        for (int k = -r; k <= r; k++) {
            int p = ((e * r + k) % g.n6 + g.n6) % g.n6;
            const Bits& m1 = g.ray[p][e];
            Bits m2 = (abs(k) >= r) ? bz() : g.line[AX[e]][SG[e] * k + r];
            if (!beq(m1, m2)) { fprintf(stderr, "geometry mismatch r=%d e=%d k=%d\n", r, e, k); exit(1); }
        }
    return g;
}

struct Pat {
    int r, cnt;
    vector<string> grid;
    Bits P, K;
    int np = 0, nk = 0;
};

typedef array<char, 6> Ins;
typedef vector<Ins> Prog;

static vector<Pat> parseInput(const string& path) {
    ifstream f(path);
    int n; f >> n;
    vector<Pat> v;
    for (int i = 0; i < n; i++) {
        Pat p; int rows, cols; f >> p.r >> rows >> cols >> p.cnt;
        p.grid.resize(rows);
        for (int j = 0; j < rows; j++) f >> p.grid[j];
        v.push_back(p);
    }
    return v;
}

static void initPat(Pat& p, const Geo& g) {
    p.P = bz(); p.K = bz();
    for (int y = 0; y < (int)p.grid.size(); y++)
        for (int x = 0; x < (int)p.grid[y].size(); x++)
            if (p.grid[y][x] == 'M') bset(p.P, g.idx.at({y, x}));
    const int NY[6] = {0, 0, 1, 1, -1, -1}, NX[6] = {2, -2, 1, -1, 1, -1};
    bforeach(p.P, [&](int c) {
        for (int d = 0; d < 6; d++) {
            auto it = g.idx.find({g.cy[c] + NY[d], g.cx[c] + NX[d]});
            if (it != g.idx.end() && !bget(p.P, it->second)) bset(p.K, it->second);
        }
    });
    p.np = bpop(p.P); p.nk = bpop(p.K);
}

static vector<Prog> parseOutput(const string& path) {
    ifstream f(path);
    vector<Prog> progs; Prog cur; string line;
    while (getline(f, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        bool blank = line.find_first_not_of(" \t") == string::npos;
        if (blank) { if (!cur.empty()) { progs.push_back(cur); cur.clear(); } continue; }
        istringstream ss(line); Ins in; string t;
        for (int e = 0; e < 6; e++) { ss >> t; in[e] = t[0]; }
        cur.push_back(in);
    }
    if (!cur.empty()) progs.push_back(cur);
    return progs;
}

static void writeOutput(const string& path, const vector<Prog>& progs) {
    FILE* f = fopen(path.c_str(), "wb");
    for (size_t i = 0; i < progs.size(); i++) {
        if (i) fputc('\n', f);
        for (auto& in : progs[i]) fprintf(f, "%c %c %c %c %c %c\n", in[0], in[1], in[2], in[3], in[4], in[5]);
    }
    fclose(f);
}

// exact simulation (port of sim6.simulate)
static int simulate(const Geo& g, const Pat& pt, const Prog& prog) {
    int pos[6], act[6] = {0, 0, 0, 0, 0, 0};
    for (int e = 0; e < 6; e++) pos[e] = e * g.r;
    Bits matter = bz(); int acc = 0;
    for (auto& in : prog) {
        for (int e = 0; e < 6; e++) {
            if (in[e] == 'F') pos[e] = (pos[e] + 1) % g.n6;
            else if (in[e] == 'B') pos[e] = (pos[e] - 1 + g.n6) % g.n6;
            else if (in[e] == 'T') act[e] ^= 1;
        }
        Bits a = bz(), b = bz(), c = bz();
        for (int e = 0; e < 6; e++) if (act[e]) {
            const Bits& m = g.ray[pos[e]][e];
            for (int i = 0; i < NW; i++) { c.w[i] |= b.w[i] & m.w[i]; b.w[i] |= a.w[i] & m.w[i]; a.w[i] |= m.w[i]; }
        }
        for (int i = 0; i < NW; i++) matter.w[i] = (matter.w[i] | (b.w[i] & ~c.w[i])) & ~c.w[i];
        if (bnone(andnot(pt.P, matter)) && bnone(pt.K & matter)) { matter = andnot(matter, pt.P); acc++; }
    }
    return acc;
}

