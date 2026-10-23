// Backward beam search for one build.
// A build is a tick sequence starting right after config X (end of the previous build / initial config) and
// ending with an acceptance. Only the LAST >=2-hit event of each cell matters, so we construct the build
// backwards from the acceptance tick: state = (config, determined set D). Going one tick earlier with config
// C', undetermined K cells must not be created (eq2), undetermined P cells must not be destroyed (ge3);
// newly determined: eq2 & P, ge3 & K. Done when D contains P; then transit ticks from X to the earliest config
// (emitters inactive while moving) complete the build.
// Modes: target X fixed (chain / first build), or cycle mode (X = the build's own acceptance config).
#pragma once
#include "sa3.h"

template <int W>
struct BBeam {
    const Pattern& p;
    int r, NL, S, nP, nK;
    BS<W> Pm, Km;
    vector<BS<W>> pairOne[3], pairTwo[3];  // per axis per (l1+1)*(NL+1)+(l2+1)
    vector<int> lineSt[6];                  // per emitter per state: line index or -1
    vector<vector<int>> nb;
    vector<vector<int>> dI;  // inactive-transit distance between emitter states
    EGraph G;

    BBeam(const Pattern& pat) : p(pat), r(pat.r), NL(2 * pat.r - 1), G(pat.r) {
        S = G.S;
        nP = p.nP; nK = p.nK;
        Pm = Km = BS<W>::zero();
        for (int i = 0; i < p.N; i++) { if (i < p.nP) Pm.set(i); else Km.set(i); }
        vector<BS<W>> lm[3];
        for (int a = 0; a < 3; a++) {
            lm[a].assign(NL, BS<W>::zero());
            for (int i = 0; i < p.N; i++) lm[a][p.lineOf[i][a]].set(i);
            int Q = (NL + 1) * (NL + 1);
            pairOne[a].assign(Q, BS<W>::zero());
            pairTwo[a].assign(Q, BS<W>::zero());
            for (int l1 = -1; l1 < NL; l1++)
                for (int l2 = -1; l2 < NL; l2++) {
                    int q = (l1 + 1) * (NL + 1) + (l2 + 1);
                    if (l1 >= 0 && l1 == l2) pairTwo[a][q] = lm[a][l1];
                    else {
                        if (l1 >= 0) pairOne[a][q] |= lm[a][l1];
                        if (l2 >= 0) pairOne[a][q] |= lm[a][l2];
                    }
                }
        }
        for (int e = 0; e < 6; e++) {
            lineSt[e].assign(S, -1);
            for (int s = 0; s < S; s++) {
                ES x = G.st(s);
                if (x.a && abs(x.k) < r) lineSt[e][s] = SGN[e] * x.k + r - 1;
            }
        }
        nb = G.nb;
        dI.assign(S, vector<int>(S, 0));
        for (int s = 0; s < S; s++)
            for (int t = 0; t < S; t++) dI[s][t] = transitDist(G.st(s), G.st(t));
    }
    static int transitDist(ES x, ES y) {
        if (x == y) return 0;
        if (x.k == y.k) return 1;
        return abs(x.k - y.k) + x.a + y.a;
    }
    inline int pairIdx(int a, int s1, int s2) const {
        return (lineSt[a][s1] + 1) * (NL + 1) + (lineSt[a + 3][s2] + 1);
    }
    inline void masks(int qa, int qb, int qc, BS<W>& eq2, BS<W>& ge3) const {
        const BS<W>&o0 = pairOne[0][qa], &t0 = pairTwo[0][qa];
        const BS<W>&o1 = pairOne[1][qb], &t1 = pairTwo[1][qb];
        const BS<W>&o2 = pairOne[2][qc], &t2 = pairTwo[2][qc];
        for (int w = 0; w < W; w++) {
            u64 a1 = o0.w[w], a2 = t0.w[w], b1 = o1.w[w], b2 = t1.w[w], c1 = o2.w[w], c2 = t2.w[w];
            u64 aa = a1 | a2, bb = b1 | b2, cc = c1 | c2;
            eq2.w[w] = (a2 & ~bb & ~cc) | (b2 & ~aa & ~cc) | (c2 & ~aa & ~bb) | (a1 & b1 & ~cc) | (a1 & c1 & ~bb) | (b1 & c1 & ~aa);
            ge3.w[w] = (a2 & (bb | cc)) | (b2 & (aa | cc)) | (c2 & (aa | bb)) | (a1 & b1 & c1);
        }
    }
    int cfgDist(const array<uint8_t, 6>& x, const array<uint8_t, 6>& y) const {
        int d = 0;
        for (int e = 0; e < 6; e++) d = max(d, dI[x[e]][y[e]]);
        return d;
    }

    struct Node {
        array<uint8_t, 6> c;
        BS<W> D;
        int parent;
        int root;  // index in layer 0 (for cycle mode)
        int sP, sK;
    };
    struct Opt {
        int B = 2000;
        int maxDepth = 30;
        double wK = 0.3;   // weight of determined K cells (relative to P = 1)
        double wD = 0.5;   // weight of transit distance
        bool cycle = false;
        array<uint8_t, 6> X;  // target config (chain mode)
        vector<array<uint8_t, 6>> layer0;  // optional explicit acceptance configs (empty = enumerate)
        u64 seed = 1;
    };
    struct Result {
        int len = INT_MAX;
        vector<Cfg> build;  // configs of the build (transit + path), in forward order
    };

    double score(const Node& n, int dist, const Opt& o) const { return n.sP + o.wK * n.sK - o.wD * dist; }

    // transit configs strictly between X and C (emitters inactive while moving), length max(dI)-1;
    // returns the config list of length T-1 where T = max(1, max dI)
    vector<Cfg> transit(const Cfg& X, const Cfg& C) const {
        int T = 0;
        for (int e = 0; e < 6; e++) T = max(T, transitDist(X[e], C[e]));
        vector<Cfg> v(max(0, T - 1));
        for (int e = 0; e < 6; e++) {
            ES x = X[e], y = C[e];
            // path: if needs moving: toggle off (if active), move inactive, wait, toggle on at the last tick.
            vector<ES> path;  // states after each of the T ticks (last = y)
            if (x == y) { for (int t = 0; t < T; t++) path.push_back(x); }
            else if (x.k == y.k) { for (int t = 0; t < T - 1; t++) path.push_back(x); path.push_back(y); }
            else {
                ES cur = x;
                if (cur.a) { cur.a = 0; path.push_back(cur); }
                while (cur.k != y.k) { cur.k += (y.k > cur.k) ? 1 : -1; path.push_back(cur); }
                while ((int)path.size() < T - (y.a ? 1 : 0)) path.push_back(cur);
                if (y.a) path.push_back(y);
            }
            for (int t = 0; t + 1 < T; t++) v[t][e] = path[t];
        }
        return v;
    }

    Result search(const Opt& o) {
        Result res;
        RNG rng(o.seed);
        vector<vector<Node>> layers;
        // ---- layer 0: valid acceptance configs ----
        vector<Node> L0;
        {
            vector<pair<double, Node>> cand;
            auto consider0 = [&](const array<uint8_t, 6>& c) {
                BS<W> eq2, ge3;
                masks(pairIdx(0, c[0], c[3]), pairIdx(1, c[1], c[4]), pairIdx(2, c[2], c[5]), eq2, ge3);
                if ((eq2 & Km).any() || (ge3 & Pm).any()) return;
                Node n;
                n.c = c;
                n.D = (eq2 & Pm) | (ge3 & Km);
                n.parent = -1;
                n.root = -1;
                n.sP = (n.D & Pm).pc();
                n.sK = (n.D & Km).pc();
                if (n.sP == 0) return;
                int dist = o.cycle ? 0 : cfgDist(o.X, c);
                cand.push_back({score(n, dist, o) + 1e-3 * rng.rd(), n});
            };
            if (!o.layer0.empty()) {
                for (auto& c : o.layer0) consider0(c);
            } else {
                // enumerate effective line pairs per axis; inactive emitters parked at X's k (or 0 in cycle mode)
                vector<vector<uint8_t>> opts(6);
                for (int e = 0; e < 6; e++) {
                    int kpark = o.cycle ? 0 : G.st(o.X[e]).k;
                    opts[e].push_back(G.id(ES{(int8_t)kpark, 0}));
                    for (int k = -(r - 1); k <= r - 1; k++) opts[e].push_back(G.id(ES{(int8_t)k, 1}));
                }
                array<uint8_t, 6> c;
                for (int a0 : opts[0]) for (int a3 : opts[3]) {
                    c[0] = a0; c[3] = a3;
                    for (int a1 : opts[1]) for (int a4 : opts[4]) {
                        c[1] = a1; c[4] = a4;
                        for (int a2 : opts[2]) for (int a5 : opts[5]) {
                            c[2] = a2; c[5] = a5;
                            consider0(c);
                        }
                    }
                }
            }
            sort(cand.begin(), cand.end(), [](auto& x, auto& y) { return x.first > y.first; });
            // dedup by (D) to keep diversity: at most 4 configs per identical D
            unordered_map<u64, int> seenD;
            for (auto& [sc, n] : cand) {
                if ((int)L0.size() >= o.B) break;
                u64 h = 1469598103934665603ULL;
                for (int w = 0; w < W; w++) h = (h ^ n.D.w[w]) * 1099511628211ULL;
                if (++seenD[h] > 4) continue;
                n.root = L0.size();
                L0.push_back(n);
            }
        }
        if (L0.empty()) return res;
        layers.push_back(L0);
        auto finish = [&](int layer, int idx) {
            // reconstruct path
            vector<array<uint8_t, 6>> path;
            int li = layer, ni = idx;
            while (li >= 0) { path.push_back(layers[li][ni].c); ni = layers[li][ni].parent; li--; }
            // path[0] = earliest config, path.back() = acceptance config
            Cfg Cfirst, Xc;
            for (int e = 0; e < 6; e++) Cfirst[e] = G.st(path[0][e]);
            const auto& rootc = path.back();
            for (int e = 0; e < 6; e++) Xc[e] = G.st(o.cycle ? rootc[e] : o.X[e]);
            vector<Cfg> b = transit(Xc, Cfirst);
            for (auto& c : path) { Cfg cc; for (int e = 0; e < 6; e++) cc[e] = G.st(c[e]); b.push_back(cc); }
            if ((int)b.size() < res.len) { res.len = b.size(); res.build = b; }
        };
        for (int i = 0; i < (int)layers[0].size(); i++)
            if (layers[0][i].sP == nP) finish(0, i);
        // ---- backward layers ----
        struct Child { double sc; int parent; array<uint8_t, 6> c; u64 h; };
        vector<Child> kids;
        for (int depth = 1; depth <= o.maxDepth; depth++) {
            if (depth + 1 >= res.len) break;
            const vector<Node>& prev = layers.back();
            kids.clear();
            double thr = -1e18;
            unordered_set<u64> seen;
            seen.reserve(o.B * 8);
            for (int pi = 0; pi < (int)prev.size(); pi++) {
                const Node& n = prev[pi];
                BS<W> U = ~n.D;
                BS<W> KU = Km & U, PU = Pm & U;
                const array<uint8_t, 6>& rootc = layers[0][n.root].c;
                // neighbour pairs per axis
                struct AP { uint8_t s1, s2; int q; };
                AP ap[3][16];
                int na[3];
                for (int a = 0; a < 3; a++) {
                    na[a] = 0;
                    for (int s1 : nb[n.c[a]])
                        for (int s2 : nb[n.c[a + 3]]) ap[a][na[a]++] = {(uint8_t)s1, (uint8_t)s2, pairIdx(a, s1, s2)};
                }
                for (int i0 = 0; i0 < na[0]; i0++)
                    for (int i1 = 0; i1 < na[1]; i1++)
                        for (int i2 = 0; i2 < na[2]; i2++) {
                            BS<W> eq2, ge3;
                            masks(ap[0][i0].q, ap[1][i1].q, ap[2][i2].q, eq2, ge3);
                            if ((eq2 & KU).any() || (ge3 & PU).any()) continue;
                            array<uint8_t, 6> c;
                            c[0] = ap[0][i0].s1; c[3] = ap[0][i0].s2;
                            c[1] = ap[1][i1].s1; c[4] = ap[1][i1].s2;
                            c[2] = ap[2][i2].s1; c[5] = ap[2][i2].s2;
                            BS<W> D2 = n.D | (eq2 & Pm) | (ge3 & Km);
                            int sP = (D2 & Pm).pc(), sK = (D2 & Km).pc();
                            int dist = cfgDist(o.cycle ? rootc : o.X, c);
                            double sc = sP + o.wK * sK - o.wD * dist + 1e-3 * rng.rd();
                            if (sc < thr) continue;
                            u64 h = 1469598103934665603ULL;
                            for (int w = 0; w < W; w++) h = (h ^ D2.w[w]) * 1099511628211ULL;
                            for (int e = 0; e < 6; e++) h = (h ^ c[e]) * 1099511628211ULL;
                            if (o.cycle) h = (h ^ (u64)n.root) * 1099511628211ULL;
                            if (!seen.insert(h).second) continue;
                            kids.push_back({sc, pi, c, h});
                            if ((int)kids.size() >= 4 * o.B) {
                                nth_element(kids.begin(), kids.begin() + o.B, kids.end(), [](const Child& x, const Child& y) { return x.sc > y.sc; });
                                kids.resize(o.B);
                                thr = kids[0].sc;
                                for (auto& k : kids) thr = min(thr, k.sc);
                            }
                        }
            }
            if (kids.empty()) break;
            if ((int)kids.size() > o.B) {
                nth_element(kids.begin(), kids.begin() + o.B, kids.end(), [](const Child& x, const Child& y) { return x.sc > y.sc; });
                kids.resize(o.B);
            }
            vector<Node> cur;
            cur.reserve(kids.size());
            for (auto& k : kids) {
                const Node& n = prev[k.parent];
                BS<W> eq2, ge3;
                masks(pairIdx(0, k.c[0], k.c[3]), pairIdx(1, k.c[1], k.c[4]), pairIdx(2, k.c[2], k.c[5]), eq2, ge3);
                Node m;
                m.c = k.c;
                m.D = n.D | (eq2 & Pm) | (ge3 & Km);
                m.parent = k.parent;
                m.root = n.root;
                m.sP = (m.D & Pm).pc();
                m.sK = (m.D & Km).pc();
                cur.push_back(m);
            }
            layers.push_back(move(cur));
            for (int i = 0; i < (int)layers.back().size(); i++)
                if (layers.back()[i].sP == nP) {
                    const Node& m = layers.back()[i];
                    int dist = cfgDist(o.cycle ? layers[0][m.root].c : o.X, m.c);
                    int len = depth + max(dist, 1);
                    if (len < res.len) finish(depth, i);
                }
        }
        return res;
    }
};
