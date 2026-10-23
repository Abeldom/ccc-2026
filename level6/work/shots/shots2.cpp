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

    // closure ignoring emitters in freeMask (they only need to be off, always satisfiable directly)
    int closureF(const Cfg& X, const Cfg& C, uint8_t freeMask) const {
        bool dir = true; int m1 = 0, m2 = 0, i1 = -1;
        for (int e = 0; e < 6; e++) { if ((freeMask >> e) & 1) continue; int s = sidx(X, e), t = sidx(C, e);
            if (!directT[s * NS + t]) dir = false;
            int v = noffT[s * NS + t]; if (v > m1) { m2 = m1; m1 = v; i1 = e; } else if (v > m2) m2 = v; }
        if (dir) return 0;
        int best = m1;
        if (i1 >= 0) { int alt = max(nliveT[sidx(X, i1) * NS + sidx(C, i1)], m2); if (alt < best) best = alt; }
        return max(best, 1);
    }

    struct Node { Cfg c; uint8_t flo; int8_t yk[6]; BS U; int root; int par; int nP, nK; };
    struct EOpt { int8_t k; uint8_t a, flo, act, idx; };
    struct AxOpt { int8_t ki, kj; uint8_t ai, aj, fi, fj, act, ii, ij; BS m1, m2; };

    void eopts(const Node& nd, int e, vector<EOpt>& v) const {
        v.clear();
        if ((nd.flo >> e) & 1) {
            v.push_back({nd.c.k[e], 0, 1, 0, 0});
            for (int k = -(r - 1); k <= r - 1; k++) v.push_back({(int8_t)k, 1, 0, 1, (uint8_t)(k + r)});
        } else {
            int k = nd.c.k[e], a = (nd.c.a >> e) & 1;
            v.push_back({(int8_t)k, (uint8_t)a, 0, 0, 0});
            if (k + 1 <= r) v.push_back({(int8_t)(k + 1), (uint8_t)a, 0, 0, 1});
            if (k - 1 >= -r) v.push_back({(int8_t)(k - 1), (uint8_t)a, 0, 0, 2});
            v.push_back({(int8_t)k, (uint8_t)(a ^ 1), 0, 0, 3});
        }
    }
    // apply option index idx to emitter e of node -> (k,a,flo)
    EOpt eoptByIdx(const Node& nd, int e, int idx) const {
        if ((nd.flo >> e) & 1) { if (idx == 0) return {nd.c.k[e], 0, 1, 0, 0}; return {(int8_t)(idx - r), 1, 0, 1, (uint8_t)idx}; }
        int k = nd.c.k[e], a = (nd.c.a >> e) & 1;
        if (idx == 1) k++; else if (idx == 2) k--; else if (idx == 3) a ^= 1;
        return {(int8_t)k, (uint8_t)a, 0, 0, (uint8_t)idx};
    }

    Cfg yOf(const Node& nd, const vector<Cfg>& roots) const {
        Cfg y = roots[nd.root];
        for (int e = 0; e < 6; e++) if (!((y.a >> e) & 1) && !((nd.flo >> e) & 1)) y.k[e] = nd.yk[e];
        return y;
    }

    struct Result { int total = INT_MAX; vector<Node> path; vector<Cfg> roots; };

    Node child(const Node& nd, int combo, int pi) const {
        Node ch = nd; ch.par = pi;
        for (int e = 0; e < 6; e++) { int idx = (combo >> (5 * e)) & 31; EOpt o = eoptByIdx(nd, e, idx);
            ch.c.k[e] = o.k; if (o.a) ch.c.a |= 1 << e; else ch.c.a &= ~(1 << e);
            if (o.flo) ch.flo |= 1 << e; else ch.flo &= ~(1 << e);
            if (o.act) ch.yk[e] = o.k; }
        BS S, D; effect(ch.c, S, D);
        for (int i = 0; i < NW; i++) ch.U[i] = nd.U[i] & ~(S[i] | D[i]);
        BS t1, t2; for (int i = 0; i < NW; i++) { t1[i] = ch.U[i] & Pm[i]; t2[i] = ch.U[i] & Km[i]; }
        ch.nP = pc(t1); ch.nK = pc(t2);
        return ch;
    }
    u64 hashNode(const Node& n) const {
        u64 h = 0x12345;
        for (int a = 0; a < 3; a++) {
            int l1 = ((n.flo >> a) & 1) ? 999 : n.c.k[a] * SGN[a] * 2 + ((n.c.a >> a) & 1);
            int l2 = ((n.flo >> (a + 3)) & 1) ? 999 : n.c.k[a + 3] * SGN[a + 3] * 2 + ((n.c.a >> (a + 3)) & 1);
            if (l1 > l2) swap(l1, l2);
            h = mix64(h ^ (u64)(l1 + 100) ^ ((u64)(l2 + 100) << 16));
        }
        for (int i = 0; i < NW; i++) h = mix64(h ^ n.U[i]);
        return h;
    }

    // mode 0: chain from X: objective = closure(X,C1) + L
    // mode 1: cycle: objective = closure(I,C1) + cnt*L + (cnt-1)*closure(Y,C1)
    Result beam(const Cfg& X, int mode, int cnt, const Params& prm) {
        int wL = (mode == 0) ? 1 : cnt, wY = (mode == 0) ? 0 : cnt - 1;
        const Cfg& X0 = (mode == 0) ? X : I;
        Result res;
        BS All{}; for (int i = 0; i < NW; i++) All[i] = Pm[i] | Km[i];
        // ---- roots
        vector<pair<int,int>> opts; opts.push_back({99, 99});
        for (int a = -(r - 1); a <= r - 1; a++) { opts.push_back({a, 99}); opts.push_back({a, a}); for (int b = a + 1; b <= r - 1; b++) opts.push_back({a, b}); }
        int no = opts.size();
        vector<array<BS,2>> om(3 * no);
        for (int ax = 0; ax < 3; ax++) for (int o = 0; o < no; o++) {
            BS x1{}, x2{}; auto [l1, l2] = opts[o];
            BS a{}, b{}; a.fill(0); b.fill(0);
            if (l1 != 99) a = mask[ax * R2 + l1 * SGN[ax] + r];
            if (l2 != 99) b = mask[(ax + 3) * R2 + l2 * SGN[ax + 3] + r];
            for (int i = 0; i < NW; i++) { x1[i] = a[i] ^ b[i]; x2[i] = a[i] & b[i]; }
            om[ax * no + o] = {x1, x2};
        }
        typedef tuple<long long,int,int,int> RT;
        priority_queue<RT, vector<RT>, greater<RT>> rh;
        for (int oa = 0; oa < no; oa++) for (int ob = 0; ob < no; ob++) {
            const BS &x1 = om[oa][0], &x2 = om[oa][1], &y1 = om[no + ob][0], &y2 = om[no + ob][1];
            BS ab0, ab1, ab2, ab3;
            for (int i = 0; i < NW; i++) { u64 xa = x1[i] | x2[i], ya = y1[i] | y2[i];
                ab1[i] = (x1[i] & ~ya) | (y1[i] & ~xa); ab2[i] = (x2[i] & ~ya) | (y2[i] & ~xa) | (x1[i] & y1[i]);
                ab3[i] = (xa | ya) & ~ab1[i] & ~ab2[i]; ab0[i] = ~(xa | ya); }
            for (int oc = 0; oc < no; oc++) {
                const BS &z1 = om[2 * no + oc][0], &z2 = om[2 * no + oc][1];
                bool ok = true; int sc = 0;
                for (int i = 0; i < NW && ok; i++) { u64 za = z1[i] | z2[i];
                    u64 S = (ab2[i] & ~za) | (ab1[i] & z1[i]) | (ab0[i] & z2[i]);
                    u64 D = ab3[i] | (ab2[i] & za) | (ab1[i] & z2[i]);
                    if ((S & Km[i]) | (D & Pm[i])) ok = false; else sc += __builtin_popcountll(S & Pm[i]) * 64 + __builtin_popcountll(D & Km[i]); }
                if (!ok || sc < 64) continue;
                long long key = (long long)sc * 1000003 + (long long)(mix64(prm.seed * 7 + oa * 1000000 + ob * 1000 + oc) % 1000003);
                if ((int)rh.size() < prm.rootN) rh.push({key, oa, ob, oc});
                else if (key > get<0>(rh.top())) { rh.pop(); rh.push({key, oa, ob, oc}); }
            }
        }
        vector<Cfg>& roots = res.roots;
        vector<vector<Node>> layers(1);
        layers.reserve(prm.maxDepth + 4);
        while (!rh.empty()) { auto [key, oa, ob, oc] = rh.top(); rh.pop();
            Cfg c; c.a = 0; int oo[3] = {oa, ob, oc};
            for (int ax = 0; ax < 3; ax++) { auto [l1, l2] = opts[oo[ax]];
                c.k[ax] = (l1 != 99) ? l1 * SGN[ax] : 0; if (l1 != 99) c.a |= 1 << ax;
                c.k[ax + 3] = (l2 != 99) ? l2 * SGN[ax + 3] : 0; if (l2 != 99) c.a |= 1 << (ax + 3); }
            BS S, D; effect(c, S, D);
            Node nd; nd.c = c; nd.root = roots.size(); nd.par = -1; nd.flo = (~c.a) & 63;
            for (int e = 0; e < 6; e++) nd.yk[e] = c.k[e];
            for (int i = 0; i < NW; i++) nd.U[i] = All[i] & ~(S[i] | D[i]);
            BS t1, t2; for (int i = 0; i < NW; i++) { t1[i] = nd.U[i] & Pm[i]; t2[i] = nd.U[i] & Km[i]; }
            nd.nP = pc(t1); nd.nK = pc(t2);
            roots.push_back(c); layers[0].push_back(nd); }
        if (roots.empty()) return res;
        int bestDepth = -1, bestPar = -1, bestCombo = -1;
        auto clYof = [&](const Node& nd) { return closureF(yOf(nd, roots), nd.c, nd.flo); };
        auto evalFinish = [&](const Node& nd) { return closureF(X0, nd.c, nd.flo) + wY * ((wY > 0) ? clYof(nd) : 0); };
        for (int i = 0; i < (int)layers[0].size(); i++) if (layers[0][i].nP == 0) {
            int tot = wL * 1 + evalFinish(layers[0][i]);
            if (tot < res.total) { res.total = tot; bestDepth = 0; bestPar = i; bestCombo = -1; } }
        vector<AxOpt> axo[3]; vector<EOpt> e1, e2;
        int bestNPprev = np;
        for (int depth = 1; depth < prm.maxDepth; depth++) {
            if (wL * (depth + 1) >= res.total) break;
            vector<Node>& cur = layers[depth - 1];
            if (cur.empty()) break;
            int mnp = INT_MAX; for (auto& nd : cur) mnp = min(mnp, nd.nP);
            bestNPprev = mnp;
            double rate = max(1.0, (double)(np - bestNPprev) / depth);
            int rate16 = (int)(rate * 16 * prm.closW);
            typedef tuple<u64,int,int> CT;
            priority_queue<CT> h1, h2;
            size_t cap1 = prm.beam * 2, cap2 = prm.finW * 2;
            for (int pi = 0; pi < (int)cur.size(); pi++) {
                const Node& nd = cur[pi];
                BS UP, UK; for (int i = 0; i < NW; i++) { UP[i] = nd.U[i] & Pm[i]; UK[i] = nd.U[i] & Km[i]; }
                Cfg Ynd; if (mode == 1) Ynd = roots[nd.root];
                for (int a = 0; a < 3; a++) {
                    axo[a].clear(); eopts(nd, a, e1); eopts(nd, a + 3, e2);
                    for (auto& oi : e1) for (auto& oj : e2) {
                        if (oi.act + oj.act > prm.maxAct) continue;
                        AxOpt o; o.ki = oi.k; o.kj = oj.k; o.ai = oi.a; o.aj = oj.a; o.fi = oi.flo; o.fj = oj.flo; o.act = oi.act + oj.act; o.ii = oi.idx; o.ij = oj.idx;
                        BS A{}, B{}; A.fill(0); B.fill(0);
                        if (oi.a) A = mask[a * R2 + oi.k + r]; if (oj.a) B = mask[(a + 3) * R2 + oj.k + r];
                        for (int i = 0; i < NW; i++) { o.m1[i] = A[i] ^ B[i]; o.m2[i] = A[i] & B[i]; }
                        axo[a].push_back(o);
                    }
                }
                for (auto& X1 : axo[0]) {
                for (auto& Y1 : axo[1]) { if (X1.act + Y1.act > prm.maxAct) continue;
                    BS ab0, ab1, ab2, ab3;
                    for (int i = 0; i < NW; i++) { u64 xa = X1.m1[i] | X1.m2[i], ya = Y1.m1[i] | Y1.m2[i];
                        ab1[i] = (X1.m1[i] & ~ya) | (Y1.m1[i] & ~xa); ab2[i] = (X1.m2[i] & ~ya) | (Y1.m2[i] & ~xa) | (X1.m1[i] & Y1.m1[i]);
                        ab3[i] = (xa | ya) & ~ab1[i] & ~ab2[i]; ab0[i] = ~(xa | ya); }
                    for (auto& Z1 : axo[2]) { if (X1.act + Y1.act + Z1.act > prm.maxAct) continue;
                        bool ok = true; int nP = 0, nK = 0;
                        for (int i = 0; i < NW; i++) { u64 za = Z1.m1[i] | Z1.m2[i];
                            u64 S = (ab2[i] & ~za) | (ab1[i] & Z1.m1[i]) | (ab0[i] & Z1.m2[i]);
                            u64 D = ab3[i] | (ab2[i] & za) | (ab1[i] & Z1.m2[i]);
                            if ((S & UK[i]) | (D & UP[i])) { ok = false; break; }
                            nP += __builtin_popcountll(UP[i] & ~S); nK += __builtin_popcountll(UK[i] & ~D); }
                        if (!ok) continue;
                        int combo = X1.ii | (Y1.ii << 5) | (Z1.ii << 10) | (X1.ij << 15) | (Y1.ij << 20) | (Z1.ij << 25);
                        // child light-weight node for closure
                        Node cl; cl.c.k[0] = X1.ki; cl.c.k[3] = X1.kj; cl.c.k[1] = Y1.ki; cl.c.k[4] = Y1.kj; cl.c.k[2] = Z1.ki; cl.c.k[5] = Z1.kj;
                        cl.c.a = X1.ai | (Y1.ai << 1) | (Z1.ai << 2) | (X1.aj << 3) | (Y1.aj << 4) | (Z1.aj << 5);
                        cl.flo = X1.fi | (Y1.fi << 1) | (Z1.fi << 2) | (X1.fj << 3) | (Y1.fj << 4) | (Z1.fj << 5);
                        cl.root = nd.root;
                        for (int e = 0; e < 6; e++) cl.yk[e] = nd.yk[e];
                        if (X1.act | Y1.act | Z1.act) { for (int e = 0; e < 6; e++) if (((nd.flo >> e) & 1) && !((cl.flo >> e) & 1)) cl.yk[e] = cl.c.k[e]; }
                        u64 noise = mix64(prm.seed ^ ((u64)depth << 40) ^ ((u64)pi << 30) ^ (u64)combo);
                        if (nP > 0) {
                            u64 est = (u64)nP * 16;
                            if (rate16 > 0) {
                                int clo = (mode == 1 && cnt > 1) ? clYof(cl) : closureF(X0, cl.c, cl.flo);
                                est = max(est, (u64)rate16 * clo);
                            }
                            u64 key = (est << 40) | ((u64)min(nK, 65535) << 24) | (noise & 0xFFFFFF);
                            if (h1.size() < cap1) h1.push({key, pi, combo}); else if (key < get<0>(h1.top())) { h1.pop(); h1.push({key, pi, combo}); }
                        } else {
                            int rem = evalFinish(cl);
                            int tot = wL * (depth + 1) + rem;
                            if (tot < res.total) { res.total = tot; bestDepth = depth; bestPar = pi; bestCombo = combo; }
                            u64 key = ((u64)rem << 32) | (noise & 0xFFFFFFFF);
                            if (h2.size() < cap2) h2.push({key, pi, combo}); else if (key < get<0>(h2.top())) { h2.pop(); h2.push({key, pi, combo}); }
                        }
                    } } }
            }
            layers.emplace_back();
            vector<Node>& nxt = layers.back();
            unordered_set<u64> seen;
            auto take = [&](priority_queue<CT>& h, int lim) {
                vector<CT> v; while (!h.empty()) { v.push_back(h.top()); h.pop(); }
                sort(v.begin(), v.end());
                int got = 0;
                for (auto& [key, pi, combo] : v) { if (got >= lim) break;
                    Node ch = child(cur[pi], combo, pi);
                    u64 hh = hashNode(ch); if (!seen.insert(hh).second) continue;
                    nxt.push_back(ch); got++; } };
            take(h2, prm.finW); take(h1, prm.beam);
        }
        if (bestDepth < 0) return res;
        Node last;
        if (bestCombo < 0) last = layers[bestDepth][bestPar];
        else last = child(layers[bestDepth - 1][bestPar], bestCombo, bestPar);
        res.path.push_back(last);
        int p = last.par;
        for (int li = bestDepth - 1; li >= 0; li--) { const Node& nd = layers[li][p]; res.path.push_back(nd); p = nd.par; }
        return res;
    }

    // concrete config sequence (forward order) from a result path, given the build's start config X
    vector<Cfg> concretize(const Result& res, const Cfg& X) const {
        const Node& fin = res.path[0];
        vector<Cfg> seq;
        for (auto& nd : res.path) {
            Cfg c = nd.c;
            for (int e = 0; e < 6; e++) if ((nd.flo >> e) & 1) {
                c.k[e] = ((fin.flo >> e) & 1) ? X.k[e] : fin.yk[e];
                c.a &= ~(1 << e);
            }
            seq.push_back(c);
        }
        return seq;
    }

    // produce program lines for: transit X->seq[0], then seq
    void emit(const Cfg& X, const vector<Cfg>& seq, vector<string>& prog, Cfg& endCfg) const {
        int liveE; int n = closure(X, seq[0], &liveE);
        vector<Cfg> tr = (n > 0) ? transit(X, seq[0], n, liveE) : vector<Cfg>();
        Cfg cur = X;
        for (auto& c : tr) { prog.push_back(stepInstr(cur, c)); cur = c; }
        for (auto& c : seq) { prog.push_back(stepInstr(cur, c)); cur = c; }
        endCfg = cur;
    }
};

// state of emitters after running prog from start
Cfg cfgAfter(const vector<string>& prog, int r, int upto) {
    Cfg c; for (int e = 0; e < 6; e++) c.k[e] = 0; c.a = 0;
    for (int t = 0; t < upto; t++) for (int e = 0; e < 6; e++) { char ch = prog[t][e];
        if (ch == 'F') c.k[e]++; else if (ch == 'B') c.k[e]--; else if (ch == 'T') c.a ^= 1 << e; }
    (void)r; return c;
}

template<int NW>
vector<string> solveT(const Pat& P, const Geo& G, const Params& prm, int modeMask, string& info) {
    int r = P.r;
    vector<pair<int,int>> Pc, Kc; set<pair<int,int>> Ps;
    for (int y = 0; y < (int)P.g.size(); y++) for (int x = 0; x < (int)P.g[y].size(); x++) if (P.g[y][x] == 'M') { Pc.push_back({y, x}); Ps.insert({y, x}); }
    int ny[6] = {0, 0, 1, 1, -1, -1}, nx[6] = {2, -2, 1, -1, 1, -1}; set<pair<int,int>> Ks;
    for (auto [y, x] : Pc) for (int d = 0; d < 6; d++) { pair<int,int> c = {y + ny[d], x + nx[d]}; if (G.interior.count(c) && !Ps.count(c)) Ks.insert(c); }
    Kc.assign(Ks.begin(), Ks.end());
    Solver<NW> S; S.init(r, Pc, Kc, G);
    vector<string> best; int bestLen = INT_MAX; string bestInfo;
    // mode 1: cycle
    if (modeMask & 2) {
        auto res = S.beam(S.I, 1, P.cnt, prm);
        if (res.total < INT_MAX) {
            vector<Cfg> seq = S.concretize(res, S.I);
            vector<string> prog; Cfg cur = S.I;
            S.emit(cur, seq, prog, cur);
            for (int b = 1; b < P.cnt; b++) S.emit(cur, seq, prog, cur);
            int t = simulate(G, P, prog, P.cnt);
            if (t > 0) { prog.resize(t); if ((int)prog.size() < bestLen) { bestLen = prog.size(); best = prog; bestInfo = "cycle L=" + to_string(seq.size()) + " pred=" + to_string(res.total); } }
            else bestInfo = "cyclefail";
        }
    }
    if ((modeMask & 1) && !(P.cnt == 1 && (modeMask & 2))) {
        vector<string> prog; Cfg cur = S.I; bool ok = true; int acc = 0;
        while (acc < P.cnt) {
            Params p2 = prm; p2.seed = prm.seed + acc * 1000 + 17;
            auto res = S.beam(cur, 0, 1, p2);
            if (res.total == INT_MAX) { ok = false; break; }
            vector<Cfg> seq = S.concretize(res, cur);
            size_t before = prog.size();
            Cfg endc; S.emit(cur, seq, prog, endc);
            int a2 = 0; int t = simulate(G, P, prog, acc + 1, &a2);
            if (t < 0) { ok = false; bestInfo += " chainfail"; break; }
            if (t < (int)before + 1) { ok = false; break; }
            prog.resize(t); acc++;
            cur = cfgAfter(prog, r, prog.size());
            if ((int)prog.size() >= bestLen) { ok = false; break; }
        }
        if (ok && (int)prog.size() < bestLen) { bestLen = prog.size(); best = prog; bestInfo = "chain"; }
    }
    info = bestInfo;
    return best;
}

int main(int argc, char** argv) {
    if (argc < 3) { fprintf(stderr, "usage: shots in out [beam finW rootN modeMask seed maxAct closW only]\n"); return 1; }
    Params prm; int modeMask = 3; int only = -1;
    if (argc > 3) prm.beam = atoi(argv[3]);
    if (argc > 4) prm.finW = atoi(argv[4]);
    if (argc > 5) prm.rootN = atoi(argv[5]);
    if (argc > 6) modeMask = atoi(argv[6]);
    if (argc > 7) prm.seed = atoll(argv[7]);
    if (argc > 8) prm.maxAct = atoi(argv[8]);
    if (argc > 9) prm.closW = atof(argv[9]);
    if (argc > 10) only = atoi(argv[10]);
    auto pats = parseInput(argv[1]);
    int N = pats.size();
    map<int, Geo> geos; for (auto& p : pats) if (!geos.count(p.r)) geos[p.r] = makeGeo(p.r);
    vector<vector<string>> out(N); vector<string> infos(N);
    vector<double> tm(N);
    #pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < N; i++) {
        if (only >= 0 && i != only) continue;
        auto t0 = chrono::steady_clock::now();
        const Pat& P = pats[i]; const Geo& G = geos.at(P.r);
        int ncell = 3 * P.r * (P.r - 1) + 1;
        int nw = (ncell + 63) / 64;
        vector<string> prog; string info;
        switch (nw) {
            case 1: prog = solveT<1>(P, G, prm, modeMask, info); break;
            case 2: prog = solveT<2>(P, G, prm, modeMask, info); break;
            case 3: prog = solveT<3>(P, G, prm, modeMask, info); break;
            case 4: prog = solveT<4>(P, G, prm, modeMask, info); break;
            case 5: prog = solveT<5>(P, G, prm, modeMask, info); break;
            case 6: prog = solveT<6>(P, G, prm, modeMask, info); break;
            case 7: prog = solveT<7>(P, G, prm, modeMask, info); break;
            default: prog = solveT<8>(P, G, prm, modeMask, info); break;
        }
        out[i] = prog; infos[i] = info;
        tm[i] = chrono::duration<double>(chrono::steady_clock::now() - t0).count();
        #pragma omp critical
        { fprintf(stderr, "case %d r=%d cnt=%d -> %d lines (%s) %.1fs\n", i, P.r, P.cnt, (int)prog.size(), info.c_str(), tm[i]); }
    }
    ofstream f(argv[2], ios::binary);
    long long tot = 0;
    for (int i = 0; i < N; i++) {
        if (out[i].empty()) out[i].push_back("- - - - - -");
        for (auto& s : out[i]) { f << s[0] << ' ' << s[1] << ' ' << s[2] << ' ' << s[3] << ' ' << s[4] << ' ' << s[5] << '\n'; }
        tot += out[i].size();
        if (i + 1 < N) f << '\n';
    }
    fprintf(stderr, "total %lld\n", tot);
}
