// core3.h: geometry, exact simulator, pattern context and generalized backward beam (shots3).
#pragma once
#include <bits/stdc++.h>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace std;
typedef unsigned long long u64;

static const int DY[6] = {0, 1, 1, 0, -1, -1}, DX[6] = {2, 1, -1, -2, -1, 1};
static const int SGN[6] = {-1, 1, 1, 1, -1, -1};

static inline u64 mix64(u64 x) { x += 0x9e3779b97f4a7c15ULL; x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL; x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL; return x ^ (x >> 31); }

struct Geo {
    int r, n;
    vector<pair<int,int>> per;
    set<pair<int,int>> interior;
    vector<array<vector<pair<int,int>>,6>> ray;
};
inline Geo makeGeo(int r) {
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
inline vector<Pat> parseInput(const string& path) {
    ifstream in(path); int n; in >> n; vector<Pat> v(n);
    for (auto& p : v) { int rows, cols; in >> p.r >> rows >> cols >> p.cnt; p.g.resize(rows); for (auto& s : p.g) in >> s; }
    return v;
}
// programs as lists of 6-char instruction strings
inline vector<vector<string>> readProgs(const string& path) {
    ifstream in(path); string line; vector<vector<string>> res; vector<string> cur;
    while (getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        string s; for (char c : line) if (c == 'F' || c == 'B' || c == 'T' || c == '-') s += c;
        if (s.empty()) { if (!cur.empty()) { res.push_back(cur); cur.clear(); } continue; }
        if (s.size() != 6) { fprintf(stderr, "bad line %s\n", line.c_str()); exit(1); }
        cur.push_back(s);
    }
    if (!cur.empty()) res.push_back(cur);
    return res;
}
inline void writeProgs(const string& path, const vector<vector<string>>& progs) {
    string tmp = path + ".tmp";
    { ofstream f(tmp, ios::binary);
      for (size_t i = 0; i < progs.size(); i++) {
          for (auto& s : progs[i]) f << s[0] << ' ' << s[1] << ' ' << s[2] << ' ' << s[3] << ' ' << s[4] << ' ' << s[5] << '\n';
          if (i + 1 < progs.size()) f << '\n';
      } }
    remove(path.c_str()); rename(tmp.c_str(), path.c_str());
}

// exact forward simulator (mirror of sim6.py). Returns tick (1-based) of the need-th acceptance or -1.
inline int simulateExact(const Geo& G, const Pat& P, const vector<string>& prog, int need, int* accOut = nullptr) {
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

struct Cfg { int8_t k[6]; uint8_t a; };
static inline bool cfgEq(const Cfg& x, const Cfg& y) { for (int e = 0; e < 6; e++) if (x.k[e] != y.k[e]) return false; return x.a == y.a; }

struct BParams { int W = 400, finW = 100, rootN = 400, maxAct = 1, maxDepth = 400; u64 seed = 1; double closW = 0.0; int sumMode = 0;
    double alpha = 0, beta = 0, gamma = 0.5, mult = 3.0, cloW = 0; };

template<int NW>
struct Ctx {
    typedef array<u64, NW> BS;
    int r = 0, M = 0, np = 0, nk = 0, ncell = 0, NS = 0;
    BS Pm{}, Km{}, All{};
    vector<BS> mask;                       // mask[e*M + k+3r], k in [-3r,3r)
    vector<int> noffT, nliveT; vector<char> directT;
    vector<int> cellLine[3];               // per cell: H, U, V line index
    Cfg I;
    const Geo* G = nullptr; const Pat* P = nullptr;

    static inline bool zero(const BS& a) { for (int i = 0; i < NW; i++) if (a[i]) return false; return true; }
    static inline int pc(const BS& a) { int s = 0; for (int i = 0; i < NW; i++) s += __builtin_popcountll(a[i]); return s; }
    inline int ki(int k) const { return k + 3 * r; }
    inline int wrapk(int k) const { if (k >= 3 * r) k -= M; if (k < -3 * r) k += M; return k; }
    inline int sidx(const Cfg& c, int e) const { return ki(c.k[e]) * 2 + ((c.a >> e) & 1); }
    inline bool offS(int s) const { int k = s / 2 - 3 * r; return !(s & 1) || abs(k) >= r; }
    inline const BS& em(int e, int k) const { return mask[e * M + ki(k)]; }

    void init(const Pat& P_, const Geo& G_) {
        P = &P_; G = &G_; r = P_.r; M = 6 * r; NS = 2 * M;
        vector<pair<int,int>> Pc, Kc; set<pair<int,int>> Ps;
        for (int y = 0; y < (int)P_.g.size(); y++) for (int x = 0; x < (int)P_.g[y].size(); x++) if (P_.g[y][x] == 'M') { Pc.push_back({y, x}); Ps.insert({y, x}); }
        int ny[6] = {0, 0, 1, 1, -1, -1}, nx[6] = {2, -2, 1, -1, 1, -1}; set<pair<int,int>> Ks;
        for (auto [y, x] : Pc) for (int d = 0; d < 6; d++) { pair<int,int> c = {y + ny[d], x + nx[d]}; if (G_.interior.count(c) && !Ps.count(c)) Ks.insert(c); }
        Kc.assign(Ks.begin(), Ks.end());
        np = Pc.size(); nk = Kc.size(); ncell = np + nk;
        if (ncell > 64 * NW) { fprintf(stderr, "NW too small\n"); exit(1); }
        map<pair<int,int>, int> id; for (int i = 0; i < np; i++) id[Pc[i]] = i; for (int i = 0; i < nk; i++) id[Kc[i]] = np + i;
        for (int c = 0; c < 3; c++) cellLine[c].assign(ncell, 0);
        for (auto& [yx, i] : id) { int y = yx.first, x = yx.second; cellLine[0][i] = y - r; cellLine[1][i] = (x - y - r) / 2; cellLine[2][i] = (x + y - 3 * r) / 2; }
        Pm.fill(0); Km.fill(0);
        for (int i = 0; i < np; i++) Pm[i >> 6] |= 1ULL << (i & 63);
        for (int i = np; i < ncell; i++) Km[i >> 6] |= 1ULL << (i & 63);
        for (int i = 0; i < NW; i++) All[i] = Pm[i] | Km[i];
        mask.assign(6 * M, BS{});
        for (int e = 0; e < 6; e++) for (int k = -3 * r; k < 3 * r; k++) {
            int p = ((e * r + k) % G_.n + G_.n) % G_.n; BS m{}; m.fill(0);
            for (auto c : G_.ray[p][e]) { auto it = id.find(c); if (it != id.end()) m[it->second >> 6] |= 1ULL << (it->second & 63); }
            mask[e * M + ki(k)] = m; }
        for (int e = 0; e < 6; e++) I.k[e] = 0; I.a = 0;
        // closure tables
        auto nbrs = [&](int s, int* v) { int k = s / 2 - 3 * r, a = s & 1; int c = 0; v[c++] = s;
            v[c++] = ki(wrapk(k + 1)) * 2 + a; v[c++] = ki(wrapk(k - 1)) * 2 + a; v[c++] = s ^ 1; return c; };
        vector<vector<int>> dist(NS, vector<int>(NS, 1 << 20)), od(NS, vector<int>(NS, 1 << 20));
        for (int s = 0; s < NS; s++) { deque<int> q; dist[s][s] = 0; q.push_back(s);
            while (!q.empty()) { int x = q.front(); q.pop_front(); int v[4]; int c = nbrs(x, v); for (int j = 0; j < c; j++) { int y = v[j]; if (dist[s][y] > dist[s][x] + 1) { dist[s][y] = dist[s][x] + 1; q.push_back(y); } } } }
        for (int s = 0; s < NS; s++) if (offS(s)) { deque<int> q; od[s][s] = 0; q.push_back(s);
            while (!q.empty()) { int x = q.front(); q.pop_front(); int v[4]; int c = nbrs(x, v); for (int j = 0; j < c; j++) { int y = v[j]; if (offS(y) && od[s][y] > od[s][x] + 1) { od[s][y] = od[s][x] + 1; q.push_back(y); } } } }
        noffT.assign(NS * NS, 1 << 20); nliveT.assign(NS * NS, 1 << 20); directT.assign(NS * NS, 0);
        for (int s = 0; s < NS; s++) for (int t = 0; t < NS; t++) {
            int v[4]; int c = nbrs(s, v); bool d = false; for (int j = 0; j < c; j++) if (v[j] == t) d = true;
            directT[s * NS + t] = d;
            nliveT[s * NS + t] = max(1, dist[s][t] - 1);
            int best = 1 << 20; int w[4]; int cw = nbrs(t, w);
            for (int j = 0; j < c; j++) if (offS(v[j])) for (int jj = 0; jj < cw; jj++) if (offS(w[jj])) best = min(best, 1 + od[v[j]][w[jj]]);
            noffT[s * NS + t] = best;
        }
    }

    inline void effect(const Cfg& c, BS& S, BS& D) const {
        BS c1{}, c2{}, c3{}; c1.fill(0); c2.fill(0); c3.fill(0);
        for (int e = 0; e < 6; e++) if ((c.a >> e) & 1) { const BS& m = em(e, c.k[e]);
            for (int i = 0; i < NW; i++) { u64 n3 = c3[i] | (c2[i] & m[i]); u64 n2 = (c2[i] & ~m[i]) | (c1[i] & m[i]); u64 n1 = (c1[i] & ~m[i]) | (~(c1[i] | c2[i] | c3[i]) & m[i]); c3[i] = n3; c2[i] = n2; c1[i] = n1; } }
        S = c2; D = c3;
    }
    // number of transit ticks strictly between X and C; emitters in freeMask only need to end up off
    inline int closure(const Cfg& X, const Cfg& C, uint8_t freeMask = 0, int* liveE = nullptr) const {
        bool dir = true; int m1 = 0, m2 = 0, i1 = -1;
        for (int e = 0; e < 6; e++) { if ((freeMask >> e) & 1) continue; int s = sidx(X, e), t = sidx(C, e);
            if (!directT[s * NS + t]) dir = false;
            int v = noffT[s * NS + t]; if (v > m1) { m2 = m1; m1 = v; i1 = e; } else if (v > m2) m2 = v; }
        if (liveE) *liveE = -1;
        if (dir) return 0;
        int best = m1;
        if (i1 >= 0) { int alt = max(nliveT[sidx(X, i1) * NS + sidx(C, i1)], m2); if (alt < best) { best = alt; if (liveE) *liveE = i1; } }
        return max(best, 1);
    }
    // explicit transit: n intermediate configs strictly between X and C (all effect-free)
    bool transit(const Cfg& X, const Cfg& C, int n, int liveE, vector<Cfg>& out) const {
        out.assign(n, X);
        for (int e = 0; e < 6; e++) {
            int s = sidx(X, e), t = sidx(C, e); bool live = (e == liveE);
            vector<vector<int>> par(n + 2, vector<int>(NS, -2));
            par[0][s] = -1;
            for (int L = 0; L <= n; L++) for (int x = 0; x < NS; x++) if (par[L][x] != -2) {
                int k = x / 2 - 3 * r, a = x & 1; int cand[4] = {x, ki(wrapk(k + 1)) * 2 + a, ki(wrapk(k - 1)) * 2 + a, x ^ 1};
                for (int y : cand) { if (L + 1 <= n && !live && !offS(y)) continue; if (L + 1 == n + 1 && y != t) continue;
                    if (par[L + 1][y] == -2) par[L + 1][y] = x; } }
            if (par[n + 1][t] == -2) return false;
            int x = t;
            for (int L = n + 1; L >= 1; L--) { int px = par[L][x]; if (L - 1 >= 1) { out[L - 2].k[e] = px / 2 - 3 * r; if (px & 1) out[L - 2].a |= 1 << e; else out[L - 2].a &= ~(1 << e); } x = px; }
        }
        return true;
    }
    // instruction between consecutive configs ("" if not adjacent)
    string stepInstr(const Cfg& A, const Cfg& B) const {
        string s(6, '-');
        for (int e = 0; e < 6; e++) {
            int da = (A.a >> e) & 1, db = (B.a >> e) & 1;
            if (A.k[e] == B.k[e] && da == db) s[e] = '-';
            else if (A.k[e] == B.k[e]) s[e] = 'T';
            else if (B.k[e] == wrapk(A.k[e] + 1) && da == db) s[e] = 'F';
            else if (B.k[e] == wrapk(A.k[e] - 1) && da == db) s[e] = 'B';
            else return "";
        }
        return s;
    }
    vector<Cfg> parseCfgs(const vector<string>& prog) const {
        vector<Cfg> c(prog.size() + 1); c[0] = I;
        for (size_t t = 0; t < prog.size(); t++) { Cfg x = c[t];
            for (int e = 0; e < 6; e++) { char ch = prog[t][e];
                if (ch == 'F') x.k[e] = wrapk(x.k[e] + 1); else if (ch == 'B') x.k[e] = wrapk(x.k[e] - 1); else if (ch == 'T') x.a ^= 1 << e; }
            c[t + 1] = x; }
        return c;
    }
    // fast forward sim on P∪K; F[t] = state after tick t; acc = ticks of acceptances
    void forward(const vector<Cfg>& cf, vector<BS>& F, vector<int>& acc) const {
        F.assign(cf.size(), BS{}); acc.clear(); BS m{}; m.fill(0); F[0] = m;
        for (size_t t = 1; t < cf.size(); t++) {
            BS S, D; effect(cf[t], S, D);
            for (int i = 0; i < NW; i++) m[i] = (m[i] | S[i]) & ~D[i];
            bool ok = true; for (int i = 0; i < NW; i++) if ((m[i] & Pm[i]) != Pm[i] || (m[i] & Km[i])) ok = false;
            if (ok) { acc.push_back(t); for (int i = 0; i < NW; i++) m[i] &= ~Pm[i]; }
            F[t] = m;
        }
    }
};

// ---------------------------------------------------------------------------------------------
// Generalized backward beam.
// Plans a window of configs D_1..D_L (forward order) that follows a forward state (X config, F matter)
// reached after a transit, and precedes either the detection tick (free root) or a fixed config Cnext
// (D_L must be adjacent to Cnext) with suffix-undecided cell set U0.
// Cost = closure(X, D_1) + L (+ extra for nextMode).
template<int NW>
struct Spec {
    typedef array<u64, NW> BS;
    Cfg X; BS F; BS U0; bool freeRoot = true; Cfg Cnext; int bound = INT_MAX;
    bool libMode = false;                     // ignore transit from X (pure build length), for library generation
    int nextMode = 0; Cfg W1; int cnt = 1;   // 1: + closure(Y, W1); 2: cycle: + (cnt-1)*(L + closure(Y, D_1)) and X = I
};

template<int NW>
struct Beam {
    typedef array<u64, NW> BS;
    const Ctx<NW>& C;
    int r, M;
    explicit Beam(const Ctx<NW>& c) : C(c), r(c.r), M(c.M) {}
    static inline int pc(const BS& a) { int s = 0; for (int i = 0; i < NW; i++) s += __builtin_popcountll(a[i]); return s; }

    struct Node { Cfg c; uint8_t flo; int8_t yk[6]; BS U; int par; int16_t need, nK, clo; };
    struct EOpt { int8_t k; uint8_t a, flo, act, idx; };
    struct AxOpt { int8_t ki, kj; uint8_t ai, aj, fi, fj, act, ii, ij; BS m1, m2; };
    struct Res { int total = INT_MAX; vector<Cfg> seq; Cfg Y; uint8_t freeM = 0; };

    void eopts(const Node& nd, int e, vector<EOpt>& v) const {
        v.clear();
        if ((nd.flo >> e) & 1) {
            v.push_back({nd.c.k[e], 0, 1, 0, 0});
            for (int k = -(r - 1); k <= r - 1; k++) v.push_back({(int8_t)k, 1, 0, 1, (uint8_t)(k + r)});
        } else {
            int k = nd.c.k[e], a = (nd.c.a >> e) & 1;
            v.push_back({(int8_t)k, (uint8_t)a, 0, 0, 0});
            v.push_back({(int8_t)C.wrapk(k + 1), (uint8_t)a, 0, 0, 1});
            v.push_back({(int8_t)C.wrapk(k - 1), (uint8_t)a, 0, 0, 2});
            v.push_back({(int8_t)k, (uint8_t)(a ^ 1), 0, 0, 3});
        }
    }
    EOpt eoptByIdx(const Node& nd, int e, int idx) const {
        if ((nd.flo >> e) & 1) { if (idx == 0) return {nd.c.k[e], 0, 1, 0, 0}; return {(int8_t)(idx - r), 1, 0, 1, (uint8_t)idx}; }
        int k = nd.c.k[e], a = (nd.c.a >> e) & 1;
        if (idx == 1) k = C.wrapk(k + 1); else if (idx == 2) k = C.wrapk(k - 1); else if (idx == 3) a ^= 1;
        return {(int8_t)k, (uint8_t)a, 0, 0, (uint8_t)idx};
    }
    inline void needOf(const BS& U, const BS& F, BS& out) const {
        for (int i = 0; i < NW; i++) out[i] = (U[i] & C.Pm[i] & ~F[i]) | (U[i] & C.Km[i] & F[i]);
    }
    Node child(const Node& nd, int combo, int pi, const Spec<NW>& sp) const {
        Node ch = nd; ch.par = pi;
        for (int e = 0; e < 6; e++) { int idx = (combo >> (5 * e)) & 31; EOpt o = eoptByIdx(nd, e, idx);
            ch.c.k[e] = o.k; if (o.a) ch.c.a |= 1 << e; else ch.c.a &= ~(1 << e);
            if (o.flo) ch.flo |= 1 << e; else ch.flo &= ~(1 << e);
            if (o.act) ch.yk[e] = o.k; }
        BS S, D; C.effect(ch.c, S, D);
        for (int i = 0; i < NW; i++) ch.U[i] = nd.U[i] & ~(S[i] | D[i]);
        BS nd2; needOf(ch.U, sp.F, nd2); ch.need = pc(nd2);
        BS t2; for (int i = 0; i < NW; i++) t2[i] = ch.U[i] & C.Km[i]; ch.nK = pc(t2);
        ch.clo = C.closure(sp.X, ch.c, ch.flo);
        return ch;
    }
    double wt[64];
    // emitter distance to line l (ticks until it can be active on that line)
    inline int edist(const Node& n, int e, int l) const {
        if ((n.flo >> e) & 1) return 1;
        int k = l * SGN[e]; int d = abs(n.c.k[e] - k); if (d > 3 * r) d = M - d;
        return ((n.c.a >> e) & 1) ? d : d + 1;
    }
    // sum over still-needed P cells of gamma^(ticks until the cell could be created)
    double potential(const Node& n, const BS& F) const {
        int L = 2 * r - 1;
        int mn[3][32], dp[3][32];
        for (int a = 0; a < 3; a++) for (int l = 0; l < L; l++) { int d1 = edist(n, a, l - (r - 1)), d2 = edist(n, a + 3, l - (r - 1)); mn[a][l] = min(d1, d2); dp[a][l] = max(d1, d2); }
        double s = 0;
        for (int w = 0; w < NW; w++) { u64 b = n.U[w] & C.Pm[w] & ~F[w];
            while (b) { int i = w * 64 + __builtin_ctzll(b); b &= b - 1;
                int h = C.cellLine[0][i] + r - 1, u = C.cellLine[1][i] + r - 1, v = C.cellLine[2][i] + r - 1;
                int x = mn[0][h], y = mn[1][u], z = mn[2][v];
                int med = max(min(x, y), min(max(x, y), z));
                int t = min(med, min(dp[0][h], min(dp[1][u], dp[2][v])));
                s += wt[min(max(t, 1), 63)]; } }
        return s;
    }
    u64 hashNode(const Node& n) const {
        u64 h = 0x12345;
        for (int a = 0; a < 3; a++) {
            int l1 = ((n.flo >> a) & 1) ? 999 : n.c.k[a] * 2 + ((n.c.a >> a) & 1);
            int l2 = ((n.flo >> (a + 3)) & 1) ? 999 : n.c.k[a + 3] * 2 + ((n.c.a >> (a + 3)) & 1);
            h = mix64(h ^ (u64)(l1 + 1000) ^ ((u64)(l2 + 1000) << 16));
        }
        for (int i = 0; i < NW; i++) h = mix64(h ^ n.U[i]);
        return h;
    }
    // end config of the window given the earliest node (floating emitters resolved)
    Cfg endCfg(const Node& fin, const Cfg& root, const Cfg& X) const {
        Cfg y = root;
        for (int e = 0; e < 6; e++) if (!((root.a >> e) & 1)) {
            if ((fin.flo >> e) & 1) y.k[e] = X.k[e]; else y.k[e] = fin.yk[e];
        }
        return y;
    }
    // concretize: forward-order configs of the window
    vector<Cfg> concretize(const vector<Node>& path, const Cfg& X) const {
        const Node& fin = path.back();
        vector<Cfg> seq;
        for (int i = (int)path.size() - 1; i >= 0; i--) {
            const Node& nd = path[i]; Cfg c = nd.c;
            for (int e = 0; e < 6; e++) if ((nd.flo >> e) & 1) { c.k[e] = ((fin.flo >> e) & 1) ? X.k[e] : fin.yk[e]; c.a &= ~(1 << e); }
            seq.push_back(c);
        }
        return seq;
    }

    Res run(const Spec<NW>& sp, const BParams& prm, vector<Res>* lib = nullptr, int libMax = 0, int libSlack = 2) {
        Res res;
        vector<array<int,4>> fins; int minFinDepth = INT_MAX;   // (depth, pi, combo, tot) of finished candidates for the library
        for (int t = 0; t < 64; t++) wt[t] = pow(prm.gamma, t - 1);
        BS UP0, UK0; for (int i = 0; i < NW; i++) { UP0[i] = sp.U0[i] & C.Pm[i]; UK0[i] = sp.U0[i] & C.Km[i]; }
        BS need0; needOf(sp.U0, sp.F, need0);
        int needInit = pc(need0);
        int wL = (sp.nextMode == 2) ? sp.cnt : 1;
        int bestTotal = sp.bound;     // only strictly better than bound
        // ---------------- roots
        vector<Node> roots; vector<Cfg> rootCfg;
        auto scoreCfg = [&](const BS& S, const BS& D, bool& ok) {
            ok = true; int sc = 0;
            for (int i = 0; i < NW; i++) {
                if ((S[i] & UK0[i]) | (D[i] & UP0[i])) { ok = false; return 0; }
                sc += 64 * __builtin_popcountll((S[i] | D[i]) & need0[i]) + __builtin_popcountll(D[i] & UK0[i] & ~need0[i]);
            }
            return sc;
        };
        auto mkRoot = [&](const Cfg& c, uint8_t flo) {
            Node nd; nd.c = c; nd.par = (int)rootCfg.size(); nd.flo = flo;
            for (int e = 0; e < 6; e++) nd.yk[e] = c.k[e];
            BS S, D; C.effect(c, S, D);
            for (int i = 0; i < NW; i++) nd.U[i] = sp.U0[i] & ~(S[i] | D[i]);
            BS nn; needOf(nd.U, sp.F, nn); nd.need = pc(nn);
            BS t2; for (int i = 0; i < NW; i++) t2[i] = nd.U[i] & C.Km[i]; nd.nK = pc(t2);
            nd.clo = C.closure(sp.X, c, flo);
            rootCfg.push_back(c); roots.push_back(nd);
        };
        typedef pair<long long, int> KI;
        if (sp.freeRoot) {
            // per-axis line options restricted to lines that pass through needed cells
            vector<int> useful[3];
            for (int ax = 0; ax < 3; ax++) { vector<char> u(2 * r - 1, 0);
                for (int i = 0; i < C.ncell; i++) if ((need0[i >> 6] >> (i & 63)) & 1) u[C.cellLine[ax][i] + r - 1] = 1;
                for (int l = 0; l < 2 * r - 1; l++) if (u[l]) useful[ax].push_back(l - (r - 1)); }
            vector<pair<int,int>> opts[3];
            for (int ax = 0; ax < 3; ax++) { opts[ax].push_back({99, 99});
                for (size_t i = 0; i < useful[ax].size(); i++) { int a = useful[ax][i]; opts[ax].push_back({a, 99}); opts[ax].push_back({a, a});
                    for (size_t j = i + 1; j < useful[ax].size(); j++) opts[ax].push_back({a, useful[ax][j]}); } }
            vector<array<BS,2>> om[3];
            for (int ax = 0; ax < 3; ax++) for (auto [l1, l2] : opts[ax]) {
                BS a{}, b{}, x1, x2; a.fill(0); b.fill(0);
                if (l1 != 99) a = C.em(ax, l1 * SGN[ax]);
                if (l2 != 99) b = C.em(ax + 3, l2 * SGN[ax + 3]);
                for (int i = 0; i < NW; i++) { x1[i] = a[i] ^ b[i]; x2[i] = a[i] & b[i]; }
                om[ax].push_back({x1, x2});
            }
            priority_queue<tuple<long long,int,int,int>, vector<tuple<long long,int,int,int>>, greater<tuple<long long,int,int,int>>> rh;
            int no0 = opts[0].size(), no1 = opts[1].size(), no2 = opts[2].size();
            for (int oa = 0; oa < no0; oa++) for (int ob = 0; ob < no1; ob++) {
                const BS &x1 = om[0][oa][0], &x2 = om[0][oa][1], &y1 = om[1][ob][0], &y2 = om[1][ob][1];
                BS ab0, ab1, ab2, ab3; bool bad = false;
                for (int i = 0; i < NW; i++) { u64 xa = x1[i] | x2[i], ya = y1[i] | y2[i];
                    ab1[i] = (x1[i] & ~ya) | (y1[i] & ~xa); ab2[i] = (x2[i] & ~ya) | (y2[i] & ~xa) | (x1[i] & y1[i]);
                    ab3[i] = (xa | ya) & ~ab1[i] & ~ab2[i]; ab0[i] = ~(xa | ya);
                    if (ab3[i] & UP0[i]) bad = true; }
                if (bad) continue;
                for (int oc = 0; oc < no2; oc++) {
                    const BS &z1 = om[2][oc][0], &z2 = om[2][oc][1];
                    bool ok = true; int sc = 0;
                    for (int i = 0; i < NW && ok; i++) { u64 za = z1[i] | z2[i];
                        u64 S = (ab2[i] & ~za) | (ab1[i] & z1[i]) | (ab0[i] & z2[i]);
                        u64 D = ab3[i] | (ab2[i] & za) | (ab1[i] & z2[i]);
                        if ((S & UK0[i]) | (D & UP0[i])) ok = false;
                        else sc += 64 * __builtin_popcountll((S | D) & need0[i]) + __builtin_popcountll(D & UK0[i] & ~need0[i]); }
                    if (!ok || sc < 64) continue;
                    long long key = (long long)sc * 1000003 + (long long)(mix64(prm.seed * 7 + oa * 1000000 + ob * 1000 + oc) % 1000003);
                    if ((int)rh.size() < prm.rootN) rh.push({key, oa, ob, oc});
                    else if (key > get<0>(rh.top())) { rh.pop(); rh.push({key, oa, ob, oc}); }
                }
            }
            while (!rh.empty()) { auto [key, oa, ob, oc] = rh.top(); rh.pop();
                Cfg c; c.a = 0; int oo[3] = {oa, ob, oc};
                for (int ax = 0; ax < 3; ax++) { auto [l1, l2] = opts[ax][oo[ax]];
                    int e1 = ax, e2 = ax + 3;
                    // assign lines to the two emitters of the axis: choose the cheaper assignment w.r.t. X
                    int ka = (l1 != 99) ? l1 * SGN[e1] : sp.X.k[e1], kb = (l2 != 99) ? l2 * SGN[e2] : sp.X.k[e2];
                    int ka2 = (l2 != 99) ? l2 * SGN[e1] : sp.X.k[e1], kb2 = (l1 != 99) ? l1 * SGN[e2] : sp.X.k[e2];
                    int c1 = abs(ka - sp.X.k[e1]) + abs(kb - sp.X.k[e2]), c2 = abs(ka2 - sp.X.k[e1]) + abs(kb2 - sp.X.k[e2]);
                    bool a1 = l1 != 99, a2 = l2 != 99;
                    if (c2 < c1) { c.k[e1] = ka2; c.k[e2] = kb2; if (a2) c.a |= 1 << e1; if (a1) c.a |= 1 << e2; }
                    else { c.k[e1] = ka; c.k[e2] = kb; if (a1) c.a |= 1 << e1; if (a2) c.a |= 1 << e2; }
                }
                mkRoot(c, (~c.a) & 63);
            }
        } else {
            // predecessors of Cnext
            int opt[6][4][2]; int no[6];
            for (int e = 0; e < 6; e++) { int k = sp.Cnext.k[e], a = (sp.Cnext.a >> e) & 1; no[e] = 0;
                opt[e][no[e]][0] = k; opt[e][no[e]++][1] = a;
                opt[e][no[e]][0] = C.wrapk(k - 1); opt[e][no[e]++][1] = a;
                opt[e][no[e]][0] = C.wrapk(k + 1); opt[e][no[e]++][1] = a;
                opt[e][no[e]][0] = k; opt[e][no[e]++][1] = a ^ 1; }
            vector<KI> cand;
            vector<Cfg> cfgs;
            for (int combo = 0; combo < 4096; combo++) {
                Cfg c; c.a = 0;
                for (int e = 0; e < 6; e++) { int j = (combo >> (2 * e)) & 3; c.k[e] = opt[e][j][0]; if (opt[e][j][1]) c.a |= 1 << e; }
                BS S, D; C.effect(c, S, D); bool ok; int sc = scoreCfg(S, D, ok);
                if (!ok) continue;
                int clo = C.closure(sp.X, c, 0);
                if (1 + clo >= bestTotal) continue;
                cand.push_back({-(long long)sc * 1000003 - (long long)(mix64(prm.seed * 13 + combo) % 1000003), (int)cfgs.size()}); cfgs.push_back(c);
            }
            sort(cand.begin(), cand.end());
            for (int i = 0; i < (int)cand.size() && i < prm.rootN; i++) mkRoot(cfgs[cand[i].second], 0);
        }
        if (roots.empty()) return res;
        // ---------------- finishing evaluation
        auto finishCost = [&](const Node& nd, int depthTicks, const Cfg& root) -> int {
            int clo = sp.libMode ? 0 : nd.clo;
            int tot = clo + wL * depthTicks;
            if (sp.nextMode == 1) {
                Cfg y = endCfg(nd, root, sp.X);
                tot += C.closure(y, sp.W1, 0);
            } else if (sp.nextMode == 2) {
                Cfg y = endCfg(nd, root, sp.X);
                // earliest concrete config: floating emitters at X (off)
                Cfg d1 = nd.c; for (int e = 0; e < 6; e++) if ((nd.flo >> e) & 1) { d1.k[e] = sp.X.k[e]; d1.a &= ~(1 << e); }
                tot += (sp.cnt - 1) * C.closure(y, d1, 0);
            }
            return tot;
        };
        vector<vector<Node>> layers; layers.reserve(prm.maxDepth + 2);
        layers.push_back(roots);
        int bestDepth = -1, bestPar = -1, bestCombo = -1;
        for (int i = 0; i < (int)layers[0].size(); i++) if (layers[0][i].need == 0) {
            int tot = finishCost(layers[0][i], 1, rootCfg[layers[0][i].par]);
            if (lib) { minFinDepth = 0; fins.push_back({0, i, -1, tot}); }
            if (tot < bestTotal) { bestTotal = tot; bestDepth = 0; bestPar = i; bestCombo = -1; } }
        // root index tracking: par of root nodes = root id; children par = index in previous layer
        vector<vector<int>> rootOf; rootOf.push_back(vector<int>(roots.size())); for (int i = 0; i < (int)roots.size(); i++) rootOf[0][i] = i;
        vector<AxOpt> axo[3]; vector<EOpt> e1, e2;
        for (int depth = 1; depth < prm.maxDepth; depth++) {
            int slackL = lib ? libSlack : 0;
            if (wL * (depth + 1) >= bestTotal + slackL) break;
            vector<Node>& cur = layers[depth - 1];
            if (cur.empty()) break;
            int mnp = INT_MAX; for (auto& nd : cur) mnp = min(mnp, (int)nd.need);
            if (getenv("BDBG")) fprintf(stderr, "      depth %d layer %d minNeed %d best %d\n", depth, (int)cur.size(), mnp, bestTotal);
            double rate = max(1.0, (double)(needInit - mnp) / depth);
            int rate16 = (int)(rate * 16 * prm.closW);
            typedef tuple<u64,int,int> CT;
            priority_queue<CT> h1, h2;
            size_t cap1 = (size_t)(prm.W * (prm.alpha > 0 ? prm.mult : 2.0)), cap2 = prm.finW * 2;
            for (int pi = 0; pi < (int)cur.size(); pi++) {
                const Node& nd = cur[pi];
                BS UP, UK; for (int i = 0; i < NW; i++) { UP[i] = nd.U[i] & C.Pm[i]; UK[i] = nd.U[i] & C.Km[i]; }
                BS NP; needOf(nd.U, sp.F, NP);
                const Cfg& rootC = rootCfg[rootOf[depth - 1][pi]];
                for (int a = 0; a < 3; a++) {
                    axo[a].clear(); eopts(nd, a, e1); eopts(nd, a + 3, e2);
                    for (auto& oi : e1) for (auto& oj : e2) {
                        if (oi.act + oj.act > prm.maxAct) continue;
                        AxOpt o; o.ki = oi.k; o.kj = oj.k; o.ai = oi.a; o.aj = oj.a; o.fi = oi.flo; o.fj = oj.flo; o.act = oi.act + oj.act; o.ii = oi.idx; o.ij = oj.idx;
                        BS A{}, B{}; A.fill(0); B.fill(0);
                        if (oi.a) A = C.em(a, oi.k); if (oj.a) B = C.em(a + 3, oj.k);
                        for (int i = 0; i < NW; i++) { o.m1[i] = A[i] ^ B[i]; o.m2[i] = A[i] & B[i]; }
                        axo[a].push_back(o);
                    }
                }
                for (auto& X1 : axo[0]) {
                for (auto& Y1 : axo[1]) { if (X1.act + Y1.act > prm.maxAct) continue;
                    BS ab0, ab1, ab2, ab3; bool bad = false;
                    for (int i = 0; i < NW; i++) { u64 xa = X1.m1[i] | X1.m2[i], ya = Y1.m1[i] | Y1.m2[i];
                        ab1[i] = (X1.m1[i] & ~ya) | (Y1.m1[i] & ~xa); ab2[i] = (X1.m2[i] & ~ya) | (Y1.m2[i] & ~xa) | (X1.m1[i] & Y1.m1[i]);
                        ab3[i] = (xa | ya) & ~ab1[i] & ~ab2[i]; ab0[i] = ~(xa | ya);
                        if (ab3[i] & UP[i]) bad = true; }
                    if (bad) continue;
                    for (auto& Z1 : axo[2]) { if (X1.act + Y1.act + Z1.act > prm.maxAct) continue;
                        bool ok = true; int nP = 0, nK = 0;
                        for (int i = 0; i < NW; i++) { u64 za = Z1.m1[i] | Z1.m2[i];
                            u64 S = (ab2[i] & ~za) | (ab1[i] & Z1.m1[i]) | (ab0[i] & Z1.m2[i]);
                            u64 D = ab3[i] | (ab2[i] & za) | (ab1[i] & Z1.m2[i]);
                            if ((S & UK[i]) | (D & UP[i])) { ok = false; break; }
                            u64 dec = S | D;
                            nP += __builtin_popcountll(NP[i] & ~dec); nK += __builtin_popcountll(UK[i] & ~D); }
                        if (!ok) continue;
                        int combo = X1.ii | (Y1.ii << 5) | (Z1.ii << 10) | (X1.ij << 15) | (Y1.ij << 20) | (Z1.ij << 25);
                        Cfg cc; cc.k[0] = X1.ki; cc.k[3] = X1.kj; cc.k[1] = Y1.ki; cc.k[4] = Y1.kj; cc.k[2] = Z1.ki; cc.k[5] = Z1.kj;
                        cc.a = X1.ai | (Y1.ai << 1) | (Z1.ai << 2) | (X1.aj << 3) | (Y1.aj << 4) | (Z1.aj << 5);
                        uint8_t fl = X1.fi | (Y1.fi << 1) | (Z1.fi << 2) | (X1.fj << 3) | (Y1.fj << 4) | (Z1.fj << 5);
                        int clo = C.closure(sp.X, cc, fl);
                        u64 noise = mix64(prm.seed ^ ((u64)depth << 40) ^ ((u64)pi << 30) ^ (u64)combo);
                        if (nP > 0) {
                            if (wL * (depth + 1) + (sp.libMode ? 1 : max(clo, 1)) >= bestTotal + slackL) continue;
                            u64 est;
                            if (prm.sumMode) est = (u64)nP * 16 + (u64)rate16 * clo;
                            else est = max((u64)nP * 16, (u64)rate16 * clo);
                            u64 key = (est << 40) | ((u64)min(nK, 65535) << 24) | (noise & 0xFFFFFF);
                            if (h1.size() < cap1) h1.push({key, pi, combo}); else if (key < get<0>(h1.top())) { h1.pop(); h1.push({key, pi, combo}); }
                        } else {
                            Node cl; cl.c = cc; cl.flo = fl; cl.clo = clo;
                            for (int e = 0; e < 6; e++) cl.yk[e] = nd.yk[e];
                            for (int e = 0; e < 6; e++) if (((nd.flo >> e) & 1) && !((fl >> e) & 1)) cl.yk[e] = cc.k[e];
                            int tot = finishCost(cl, depth + 1, rootC);
                            if (lib && depth <= minFinDepth + libSlack && fins.size() < 200000) { minFinDepth = min(minFinDepth, depth); fins.push_back({depth, pi, combo, tot}); }
                            if (tot < bestTotal) { bestTotal = tot; bestDepth = depth; bestPar = pi; bestCombo = combo; }
                            if (sp.nextMode != 0) {
                                u64 key = ((u64)(tot) << 32) | (noise & 0xFFFFFFFF);
                                if (h2.size() < cap2) h2.push({key, pi, combo}); else if (key < get<0>(h2.top())) { h2.pop(); h2.push({key, pi, combo}); }
                            }
                        }
                    } } }
            }
            layers.emplace_back(); rootOf.emplace_back();
            vector<Node>& nxt = layers.back(); vector<int>& nro = rootOf.back();
            vector<Node>& curL = layers[depth - 1];
            unordered_set<u64> seen;
            auto take = [&](priority_queue<CT>& h, int lim, bool refine) {
                vector<CT> v; while (!h.empty()) { v.push_back(h.top()); h.pop(); }
                sort(v.begin(), v.end());
                if (!refine || prm.alpha <= 0) {
                    int got = 0;
                    for (auto& [key, pi, combo] : v) { if (got >= lim) break;
                        Node ch = child(curL[pi], combo, pi, sp);
                        u64 hh = hashNode(ch); if (!seen.insert(hh).second) continue;
                        nxt.push_back(ch); nro.push_back(rootOf[depth - 1][pi]); got++; }
                    return;
                }
                vector<Node> chs; chs.reserve(v.size()); vector<pair<double,int>> sc; sc.reserve(v.size());
                for (auto& [key, pi, combo] : v) {
                    Node ch = child(curL[pi], combo, pi, sp);
                    double pot = potential(ch, sp.F);
                    double s = ch.need - prm.alpha * pot + prm.beta * ch.nK + prm.cloW * ch.clo + (double)(key & 0xFFFFFF) * 1e-9;
                    sc.push_back({s, (int)chs.size()}); chs.push_back(ch);
                }
                sort(sc.begin(), sc.end());
                int got = 0;
                for (auto& [s, ci] : sc) { if (got >= lim) break;
                    Node& ch = chs[ci];
                    u64 hh = hashNode(ch); if (!seen.insert(hh).second) continue;
                    nxt.push_back(ch); nro.push_back(rootOf[depth - 1][ch.par]); got++; }
            };
            take(h2, prm.finW, false); take(h1, prm.W, true);
        }
        auto rebuild = [&](int bD, int bP, int bC, Res& out) {
            vector<Node> path; int rootId;
            if (bC < 0) { path.push_back(layers[0][bP]); rootId = rootOf[0][bP]; }
            else {
                Node last = child(layers[bD - 1][bP], bC, bP, sp);
                path.push_back(last); rootId = rootOf[bD - 1][bP];
                int p = last.par;
                for (int li = bD - 1; li >= 0; li--) { const Node& nd = layers[li][p]; path.push_back(nd); p = nd.par; }
            }
            reverse(path.begin(), path.end());   // path[0] = root ... path.back() = earliest
            out.seq = concretize(path, sp.X);
            out.Y = endCfg(path.back(), rootCfg[rootId], sp.X);
            out.freeM = path.back().flo;
        };
        if (lib && !fins.empty()) {
            sort(fins.begin(), fins.end(), [](const array<int,4>& a, const array<int,4>& b) { return a[0] != b[0] ? a[0] < b[0] : a[3] < b[3]; });
            unordered_set<u64> seenL; int got = 0;
            for (auto& f : fins) { if (got >= libMax) break; if (f[0] > minFinDepth + libSlack) break;
                Res e; rebuild(f[0], f[1], f[2], e); e.total = f[3];
                u64 h = 0x777; auto hc = [&](const Cfg& c) { for (int q = 0; q < 6; q++) h = mix64(h ^ (u64)(c.k[q] + 100) ^ ((u64)q << 20)); h = mix64(h ^ c.a); };
                hc(e.seq[0]); hc(e.Y); h = mix64(h ^ e.seq.size() ^ ((u64)e.freeM << 30));
                if (!seenL.insert(h).second) continue;
                lib->push_back(e); got++; }
        }
        if (bestDepth < 0) return res;
        rebuild(bestDepth, bestPar, bestCombo, res);
        res.total = bestTotal;
        return res;
    }
};
