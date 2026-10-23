// SA over config sequences (generalised): arbitrary start config, optional end adjacency, extra
// adjacency edges (cycle closure), segments that each must build P from canonical and accept at their end.
// Supports warm starts (run from the current state) and tick deletion with adjacency repair.
#pragma once
#include "common.h"
typedef unsigned long long u64;

template <int W>
struct BS {
    u64 w[W];
    static BS zero() { BS b; for (int i = 0; i < W; i++) b.w[i] = 0; return b; }
    BS operator|(const BS& o) const { BS b; for (int i = 0; i < W; i++) b.w[i] = w[i] | o.w[i]; return b; }
    BS operator&(const BS& o) const { BS b; for (int i = 0; i < W; i++) b.w[i] = w[i] & o.w[i]; return b; }
    BS operator~() const { BS b; for (int i = 0; i < W; i++) b.w[i] = ~w[i]; return b; }
    BS& operator|=(const BS& o) { for (int i = 0; i < W; i++) w[i] |= o.w[i]; return *this; }
    BS& operator&=(const BS& o) { for (int i = 0; i < W; i++) w[i] &= o.w[i]; return *this; }
    int pc() const { int s = 0; for (int i = 0; i < W; i++) s += __builtin_popcountll(w[i]); return s; }
    bool any() const { u64 s = 0; for (int i = 0; i < W; i++) s |= w[i]; return s != 0; }
    void set(int i) { w[i >> 6] |= 1ULL << (i & 63); }
};

struct RNG {
    u64 s;
    RNG(u64 seed) : s(seed * 0x9E3779B97F4A7C15ULL + 12345) { if (!s) s = 1; }
    u64 next() { s ^= s << 7; s ^= s >> 9; return s; }
    int rint(int n) { return (int)((next() >> 11) % (u64)n); }
    double rd() { return (next() >> 11) * (1.0 / 9007199254740992.0); }
};

struct EGraph {
    int r, S;
    vector<vector<int>> nb;
    vector<vector<char>> adj;
    EGraph(int r_) : r(r_) {
        S = 2 * (2 * r + 1);
        nb.assign(S, {});
        adj.assign(S, vector<char>(S, 0));
        for (int s = 0; s < S; s++)
            for (int t = 0; t < S; t++)
                if (adjES(st(s), st(t))) { nb[s].push_back(t); adj[s][t] = 1; }
    }
    ES st(int s) const { return ES{(int8_t)(s / 2 - r), (int8_t)(s % 2)}; }
    int id(ES x) const { return (x.k + r) * 2 + x.a; }
};

static Cfg initCfg() { Cfg c; for (int e = 0; e < 6; e++) c[e] = {0, 0}; return c; }

struct Problem {
    int T = 0;
    vector<pair<int, int>> extra;  // extra adjacency edges (a, b)
    bool fixStart = true;          // tick 0 adjacent to startCfg
    Cfg startCfg = initCfg();
    bool fixEnd = false;           // tick T-1 adjacent to endCfg
    Cfg endCfg = initCfg();
    vector<pair<int, int>> segs;   // [start, end] inclusive
};

// cycle with segment lengths lens (closure edge), no start fix
static Problem makeCycle(const vector<int>& lens) {
    Problem pr;
    int t = 0;
    for (int l : lens) { pr.segs.push_back({t, t + l - 1}); t += l; }
    pr.T = t;
    pr.fixStart = false;
    if (t >= 2) pr.extra.push_back({0, t - 1});
    return pr;
}
// open chain of segments from start config
static Problem makeChain(const vector<int>& lens, const Cfg& start) {
    Problem pr;
    int t = 0;
    for (int l : lens) { pr.segs.push_back({t, t + l - 1}); t += l; }
    pr.T = t;
    pr.fixStart = true;
    pr.startCfg = start;
    return pr;
}

template <int W>
struct SA {
    const Pattern& p;
    int r, NL;
    EGraph G;
    BS<W> Pm, Km, ALL;
    vector<BS<W>> lm[3];
    Problem pr;
    vector<array<int, 6>> s;
    vector<BS<W>> eq2, ge3;
    vector<vector<int>> extraNb;
    vector<int> segOf;  // per tick: segment index
    vector<int> segC;
    int total = 0;
    RNG rng;
    int startId[6], endId[6];

    SA(const Pattern& pat, u64 seed) : p(pat), r(pat.r), NL(2 * pat.r - 1), G(pat.r), rng(seed) {
        Pm = Km = ALL = BS<W>::zero();
        for (int i = 0; i < p.N; i++) { ALL.set(i); if (i < p.nP) Pm.set(i); else Km.set(i); }
        for (int a = 0; a < 3; a++) {
            lm[a].assign(NL, BS<W>::zero());
            for (int i = 0; i < p.N; i++) lm[a][p.lineOf[i][a]].set(i);
        }
    }
    inline int lineOf(int e, int sid) const {
        int k = sid / 2 - r, a = sid & 1;
        if (!a || abs(k) >= r) return -1;
        return SGN[e] * k + r - 1;
    }
    void tickMasks(int t) {
        BS<W> one[3], two[3];
        for (int a = 0; a < 3; a++) {
            int l1 = lineOf(a, s[t][a]), l2 = lineOf(a + 3, s[t][a + 3]);
            one[a] = two[a] = BS<W>::zero();
            if (l1 >= 0 && l1 == l2) two[a] = lm[a][l1];
            else {
                if (l1 >= 0) one[a] |= lm[a][l1];
                if (l2 >= 0) one[a] |= lm[a][l2];
            }
        }
        BS<W> aa = one[0] | two[0], bb = one[1] | two[1], cc = one[2] | two[2];
        eq2[t] = (two[0] & ~bb & ~cc) | (two[1] & ~aa & ~cc) | (two[2] & ~aa & ~bb) | (one[0] & one[1] & ~cc) |
                 (one[0] & one[2] & ~bb) | (one[1] & one[2] & ~aa);
        ge3[t] = (two[0] & (bb | cc)) | (two[1] & (aa | cc)) | (two[2] & (aa | bb)) | (one[0] & one[1] & one[2]);
    }
    int segCost(int j) const {
        int a = pr.segs[j].first, b = pr.segs[j].second;
        BS<W> und = ALL, bad = BS<W>::zero();
        for (int t = b; t >= a; t--) {
            BS<W> x = und & eq2[t];
            bad |= x & Km;
            und &= ~x;
            BS<W> y = und & ge3[t];
            bad |= y & Pm;
            und &= ~y;
        }
        bad |= und & Pm;
        return bad.pc();
    }
    void fullEval() {
        eq2.assign(pr.T, BS<W>::zero());
        ge3.assign(pr.T, BS<W>::zero());
        for (int t = 0; t < pr.T; t++) tickMasks(t);
        segC.assign(pr.segs.size(), 0);
        total = 0;
        for (size_t j = 0; j < pr.segs.size(); j++) { segC[j] = segCost(j); total += segC[j]; }
    }
    bool okState(int e, int t, int sid, int lo, int hi) const {
        if (t == 0 && pr.fixStart && !G.adj[startId[e]][sid]) return false;
        if (t == pr.T - 1 && pr.fixEnd && !G.adj[endId[e]][sid]) return false;
        if (t - 1 >= 0 && (t - 1 < lo || t - 1 > hi) && !G.adj[s[t - 1][e]][sid]) return false;
        if (t + 1 < pr.T && (t + 1 < lo || t + 1 > hi) && !G.adj[s[t + 1][e]][sid]) return false;
        for (int o : extraNb[t])
            if ((o < lo || o > hi) && !G.adj[s[o][e]][sid]) return false;
        return true;
    }
    bool resample(int e, int lo, int hi, vector<int>& out) {
        int len = hi - lo + 1, S = G.S;
        static thread_local vector<double> cnt;
        cnt.assign((size_t)len * S, 0.0);
        static thread_local vector<char> ok;
        ok.assign((size_t)len * S, 0);
        for (int i = 0; i < len; i++)
            for (int x = 0; x < S; x++) ok[i * S + x] = okState(e, lo + i, x, lo, hi);
        for (int x = 0; x < S; x++) cnt[(len - 1) * S + x] = ok[(len - 1) * S + x] ? 1.0 : 0.0;
        for (int i = len - 2; i >= 0; i--)
            for (int x = 0; x < S; x++) {
                if (!ok[i * S + x]) continue;
                double c = 0;
                for (int y : G.nb[x]) c += cnt[(i + 1) * S + y];
                cnt[i * S + x] = c;
            }
        double tot = 0;
        for (int x = 0; x < S; x++) tot += cnt[x];
        if (tot <= 0) return false;
        out.resize(len);
        double z = rng.rd() * tot;
        int cur = -1;
        for (int x = 0; x < S; x++) { if (cnt[x] <= 0) continue; z -= cnt[x]; cur = x; if (z < 0) break; }
        out[0] = cur;
        for (int i = 1; i < len; i++) {
            double c = 0;
            for (int y : G.nb[cur]) c += cnt[i * S + y];
            double z2 = rng.rd() * c;
            int nxt = -1;
            for (int y : G.nb[cur]) { if (cnt[i * S + y] <= 0) continue; z2 -= cnt[i * S + y]; nxt = y; if (z2 < 0) break; }
            cur = nxt;
            out[i] = cur;
        }
        for (auto [a, b] : pr.extra)
            if (a >= lo && a <= hi && b >= lo && b <= hi && !G.adj[out[a - lo]][out[b - lo]]) return false;
        return true;
    }
    bool validAll() const {
        for (int e = 0; e < 6; e++) {
            if (pr.fixStart && !G.adj[startId[e]][s[0][e]]) return false;
            if (pr.fixEnd && !G.adj[endId[e]][s[pr.T - 1][e]]) return false;
            for (int t = 0; t + 1 < pr.T; t++) if (!G.adj[s[t][e]][s[t + 1][e]]) return false;
            for (auto [a, b] : pr.extra) if (!G.adj[s[a][e]][s[b][e]]) return false;
        }
        return true;
    }
    void setup(const Problem& prob) {
        pr = prob;
        s.assign(pr.T, {});
        extraNb.assign(pr.T, {});
        for (auto [a, b] : pr.extra) { extraNb[a].push_back(b); extraNb[b].push_back(a); }
        segOf.assign(pr.T, -1);
        for (size_t j = 0; j < pr.segs.size(); j++)
            for (int t = pr.segs[j].first; t <= pr.segs[j].second; t++) segOf[t] = j;
        for (int e = 0; e < 6; e++) { startId[e] = G.id(pr.startCfg[e]); endId[e] = G.id(pr.endCfg[e]); }
    }
    bool initRandom() {
        int base = pr.fixStart ? -1 : G.id(ES{0, 0});
        for (int t = 0; t < pr.T; t++) for (int e = 0; e < 6; e++) s[t][e] = base >= 0 ? base : startId[e];
        vector<int> out;
        bool all = true;
        for (int e = 0; e < 6; e++) {
            bool ok = false;
            for (int tries = 0; tries < 30 && !ok; tries++)
                if (resample(e, 0, pr.T - 1, out)) { for (int t = 0; t < pr.T; t++) s[t][e] = out[t]; ok = true; }
            all &= ok;
        }
        fullEval();
        return all;
    }
    void loadCfgs(const vector<Cfg>& v) {
        for (int t = 0; t < pr.T; t++) for (int e = 0; e < 6; e++) s[t][e] = G.id(v[t][e]);
        fullEval();
    }
    vector<array<int, 6>> best;
    // anneal from the current state; returns best cost (0 = feasible, stored in best)
    int anneal(long long iters, double T0, double T1end) {
        int cur = total;
        int bestc = cur;
        best = s;
        if (cur == 0 && validAll()) return 0;
        vector<int> out;
        vector<array<int, 6>> saved;
        vector<int> nc;
        double lt0 = log(T0), lt1 = log(T1end);
        int nseg = pr.segs.size();
        for (long long it = 0; it < iters; it++) {
            double temp = exp(lt0 + (lt1 - lt0) * it / iters);
            int ne = (rng.rint(4) == 0) ? 2 : 1;
            int len = 1 + rng.rint(min(pr.T, 4));
            if (rng.rint(10) == 0) len = pr.T;
            int lo = rng.rint(pr.T - len + 1), hi = lo + len - 1;
            int es[2];
            es[0] = rng.rint(6);
            if (ne == 2) {
                if (rng.rint(2)) es[1] = (es[0] + 3) % 6;
                else { do es[1] = rng.rint(6); while (es[1] == es[0]); }
            }
            saved.assign(s.begin() + lo, s.begin() + hi + 1);
            bool okm = true;
            for (int q = 0; q < ne && okm; q++) {
                if (!resample(es[q], lo, hi, out)) { okm = false; break; }
                for (int t = lo; t <= hi; t++) s[t][es[q]] = out[t - lo];
            }
            if (!okm) { for (int t = lo; t <= hi; t++) s[t] = saved[t - lo]; continue; }
            for (int t = lo; t <= hi; t++) tickMasks(t);
            int j0 = segOf[lo], j1 = segOf[hi];
            if (j0 < 0) j0 = 0;
            if (j1 < 0) j1 = nseg - 1;
            int nt = total;
            nc.resize(j1 - j0 + 1);
            for (int j = j0; j <= j1; j++) { nc[j - j0] = segCost(j); nt += nc[j - j0] - segC[j]; }
            if (nt <= cur || rng.rd() < exp((cur - nt) / temp)) {
                cur = total = nt;
                for (int j = j0; j <= j1; j++) segC[j] = nc[j - j0];
                if (cur < bestc) { bestc = cur; best = s; if (cur == 0) return validAll() ? 0 : 1; }
            } else {
                for (int t = lo; t <= hi; t++) { s[t] = saved[t - lo]; tickMasks(t); }
            }
        }
        return bestc;
    }
    int run(long long iters, double T0, double T1end) {
        if (!initRandom()) return INT_MAX / 2;
        return anneal(iters, T0, T1end);
    }
    vector<Cfg> bestCfgs() const {
        vector<Cfg> v(pr.T);
        for (int t = 0; t < pr.T; t++) for (int e = 0; e < 6; e++) v[t][e] = G.st(best[t][e]);
        return v;
    }
    vector<Cfg> curCfgs() const {
        vector<Cfg> v(pr.T);
        for (int t = 0; t < pr.T; t++) for (int e = 0; e < 6; e++) v[t][e] = G.st(s[t][e]);
        return v;
    }
    // repair adjacency violations by resampling windows around bad edges; returns success
    bool repair() {
        vector<int> out;
        for (int pass = 0; pass < 6; pass++) {
            bool any = false;
            for (int e = 0; e < 6; e++) {
                vector<int> badT;
                if (pr.fixStart && !G.adj[startId[e]][s[0][e]]) badT.push_back(0);
                if (pr.fixEnd && !G.adj[endId[e]][s[pr.T - 1][e]]) badT.push_back(pr.T - 1);
                for (int t = 0; t + 1 < pr.T; t++) if (!G.adj[s[t][e]][s[t + 1][e]]) badT.push_back(t);
                for (auto [a, b] : pr.extra) if (!G.adj[s[a][e]][s[b][e]]) { badT.push_back(a); badT.push_back(b); }
                for (int d : badT) {
                    any = true;
                    bool fixed = false;
                    for (int w = 1; w <= pr.T && !fixed; w++) {
                        int lo = max(0, d - w), hi = min(pr.T - 1, d + w);
                        if (resample(e, lo, hi, out)) { for (int t = lo; t <= hi; t++) s[t][e] = out[t - lo]; fixed = true; }
                        if (lo == 0 && hi == pr.T - 1) break;
                    }
                    if (!fixed) return false;
                }
            }
            if (!any) break;
        }
        fullEval();
        return validAll();
    }
};
