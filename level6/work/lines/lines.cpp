// lines: line-cover planner. Pairs (the two emitters of an axis) fire whole lines; the other emitters
// cut the skirt cells on those lines. A plan = per build a list of fires (axis, line, tick); every cell
// gets a "finalize" tick (last fire through it). Emitter trajectories are then optimised by an
// alternating exact DP over axis pairs (each DP optimises the two emitters of one axis given the
// other four) under an additive violation model; the result is checked by exact simulation and
// remaining defects are fixed by a cleanup beam search.
// usage: lines <in> <best_out> <res_dir> [opts]
#include "common.h"

static int VERB = 0;
constexpr int NONE = -100;

struct ES { signed char k; unsigned char a; };
typedef array<ES, 6> Conf;
typedef vector<Conf> Sched;  // S[0] = conf before tick 1, S[t] = conf during tick t

static inline int lineOf(const Geo& g, int e, const ES& s) {
    if (!s.a || s.k >= g.r || s.k <= -g.r) return NONE;
    return SG[e] * s.k;
}

static Prog toProg(const Sched& S, int from = 1, int to = -1) {
    if (to < 0) to = (int)S.size() - 1;
    Prog p;
    for (int t = from; t <= to; t++) {
        Ins in;
        for (int e = 0; e < 6; e++) {
            const ES &a = S[t - 1][e], &b = S[t][e];
            if (a.a != b.a) { if (a.k != b.k) { fprintf(stderr, "bad sched\n"); exit(1); } in[e] = 'T'; }
            else if (b.k == a.k + 1) in[e] = 'F';
            else if (b.k == a.k - 1) in[e] = 'B';
            else { if (a.k != b.k) { fprintf(stderr, "bad sched2\n"); exit(1); } in[e] = '-'; }
        }
        p.push_back(in);
    }
    return p;
}

static Conf startConf() { Conf c; for (int e = 0; e < 6; e++) { c[e].k = 0; c[e].a = 0; } return c; }

static inline void confMasks(const Geo& g, const Conf& c, const Bits* m[6], const Bits& Z) {
    for (int e = 0; e < 6; e++) { int l = lineOf(g, e, c[e]); m[e] = l == NONE ? &Z : &g.line[AX[e]][l + g.r]; }
}

// exact simulation of ticks t0..t1 of schedule; appends acceptance ticks
static void simRange(const Geo& g, const Pat& pt, const Sched& S, int t0, int t1, Bits& board, vector<int>& accT) {
    const Bits Z = bz();
    for (int t = t0; t <= t1; t++) {
        const Bits* m[6]; confMasks(g, S[t], m, Z);
        tickBoard(board, m);
        if (bnone(andnot(pt.P, board)) && bnone(pt.K & board)) { board = andnot(board, pt.P); accT.push_back(t); }
    }
}

// ---------------------------------------------------------------------------------------------
struct Fire { int ax, line, t; };

struct Model {
    const Geo& g; const Pat& pt;
    int t0 = 1, T = 0;                 // optimise ticks t0..T (S[t0-1] fixed)
    vector<int> bOf;                   // tick -> build index
    vector<vector<short>> fx;          // [build][cell] finalize tick or -1
    vector<unsigned char> typ;         // 0 free, 1 P, 2 K
    vector<int> rel;                   // relevant cells
    vector<int> lcells[3];             // flattened: [axis][(c+r)] -> list in lc
    vector<vector<int>> lc[3];
    int WP = 10, WK = 10;
    vector<int> wt;                    // per-cell violation weight
    Model(const Geo& g_, const Pat& p_) : g(g_), pt(p_) {
        typ.assign(g.nc, 0);
        for (int x = 0; x < g.nc; x++) {
            if (bget(pt.P, x)) typ[x] = 1; else if (bget(pt.K, x)) typ[x] = 2;
            if (typ[x]) rel.push_back(x);
        }
        for (int a = 0; a < 3; a++) {
            lc[a].assign(2 * g.r + 1, {});
            for (int x : rel) lc[a][g.L[a][x] + g.r].push_back(x);
        }
    }
    void initWeights() { wt.assign(g.nc, 0); for (int x : rel) wt[x] = typ[x] == 1 ? WP : WK; }
    inline int viol(int x, int t, int h) const {
        int ty = typ[x];
        int f = fx[bOf[t]][x];
        if (ty == 1) { if (f < 0 || t < f) return 0; if (t == f) return h == 2 ? 0 : wt[x]; return h >= 3 ? wt[x] : 0; }
        if (f >= 0 && t < f) return 0;
        return h == 2 ? wt[x] : 0;
    }
};

static long long modelCost(const Model& M, const Sched& S) {
    const Geo& g = M.g;
    long long c = 0;
    for (int t = M.t0; t <= M.T; t++) {
        int ln[6]; for (int e = 0; e < 6; e++) ln[e] = lineOf(g, e, S[t][e]);
        for (int x : M.rel) {
            int h = 0;
            for (int e = 0; e < 6; e++) if (ln[e] == g.L[AX[e]][x]) h++;
            c += M.viol(x, t, h);
        }
    }
    return c;
}

struct DPWork {
    vector<int> cur, nxt;
    vector<vector<int>> par;
    vector<vector<int>> d1, d2;
};

// optimise emitters X and X+3 given the others; returns model cost of the new schedule
static long long dpAxis(const Model& M, Sched& S, int X, mt19937* rng, int noise, DPWork& W) {
    const Geo& g = M.g; const int r = g.r;
    const int e1 = X, e2 = X + 3;
    const int nk = 2 * r + 1, NS1 = nk * 2, NS = NS1 * NS1, NL = 2 * r + 1;
    const int T = M.T, t0 = M.t0;
    W.d1.assign(T + 1, vector<int>(NL, 0));
    W.d2.assign(T + 1, vector<int>(NL, 0));
    long long base = 0;
    vector<int> fh(g.nc, 0);
    for (int t = t0; t <= T; t++) {
        int ln[6]; for (int e = 0; e < 6; e++) ln[e] = (e == e1 || e == e2) ? NONE : lineOf(g, e, S[t][e]);
        for (int x : M.rel) {
            int h = 0;
            for (int e = 0; e < 6; e++) if (ln[e] != NONE && ln[e] == g.L[AX[e]][x]) h++;
            fh[x] = h;
            base += M.viol(x, t, h);
        }
        for (int c = -(r - 1); c <= r - 1; c++) {
            int s1 = 0, s2 = 0;
            for (int x : M.lc[X][c + r]) {
                int v0 = M.viol(x, t, fh[x]);
                s1 += M.viol(x, t, fh[x] + 1) - v0;
                s2 += M.viol(x, t, fh[x] + 2) - v0;
            }
            if (noise && rng) { s1 += (int)((*rng)() % (noise + 1)); s2 += (int)((*rng)() % (noise + 1)); }
            W.d1[t][c + r] = s1; W.d2[t][c + r] = s2;
        }
    }
    const int INF = INT_MAX / 2;
    W.cur.assign(NS, INF); W.nxt.assign(NS, INF);
    if ((int)W.par.size() < T + 1) W.par.resize(T + 1);
    auto enc = [&](int k1, int a1, int k2, int a2) { return ((k1 + r) * 2 + a1) * NS1 + (k2 + r) * 2 + a2; };
    {
        const ES &s1 = S[t0 - 1][e1], &s2 = S[t0 - 1][e2];
        W.cur[enc(s1.k, s1.a, s2.k, s2.a)] = 0;
    }
    // per-state line cache
    vector<int> l1(NS), l2(NS);
    for (int s = 0; s < NS; s++) {
        int q1 = s / NS1, q2 = s % NS1;
        ES a; a.k = q1 / 2 - r; a.a = q1 & 1; ES b; b.k = q2 / 2 - r; b.a = q2 & 1;
        l1[s] = lineOf(g, e1, a); l2[s] = lineOf(g, e2, b);
    }
    for (int t = t0; t <= T; t++) {
        vector<int>& par = W.par[t];
        par.assign(NS, -1);
        std::fill(W.nxt.begin(), W.nxt.end(), INF);
        for (int s = 0; s < NS; s++) {
            int cv = W.cur[s];
            if (cv >= INF) continue;
            int q1 = s / NS1, q2 = s % NS1;
            int k1 = q1 / 2 - r, a1 = q1 & 1, k2 = q2 / 2 - r, a2 = q2 & 1;
            int o1[4][2], n1 = 0, o2[4][2], n2 = 0;
            o1[n1][0] = k1; o1[n1][1] = a1; n1++;
            if (k1 < r) { o1[n1][0] = k1 + 1; o1[n1][1] = a1; n1++; }
            if (k1 > -r) { o1[n1][0] = k1 - 1; o1[n1][1] = a1; n1++; }
            o1[n1][0] = k1; o1[n1][1] = a1 ^ 1; n1++;
            o2[n2][0] = k2; o2[n2][1] = a2; n2++;
            if (k2 < r) { o2[n2][0] = k2 + 1; o2[n2][1] = a2; n2++; }
            if (k2 > -r) { o2[n2][0] = k2 - 1; o2[n2][1] = a2; n2++; }
            o2[n2][0] = k2; o2[n2][1] = a2 ^ 1; n2++;
            for (int i = 0; i < n1; i++) {
                int base1 = ((o1[i][0] + r) * 2 + o1[i][1]) * NS1;
                for (int j = 0; j < n2; j++) {
                    int s2 = base1 + (o2[j][0] + r) * 2 + o2[j][1];
                    if (cv < W.nxt[s2]) { W.nxt[s2] = cv; par[s2] = s; }
                }
            }
        }
        const vector<int>& D1 = W.d1[t]; const vector<int>& D2 = W.d2[t];
        for (int s = 0; s < NS; s++) {
            if (W.nxt[s] >= INF) continue;
            int a = l1[s], b = l2[s], c;
            if (a != NONE && a == b) c = D2[a + r];
            else c = (a != NONE ? D1[a + r] : 0) + (b != NONE ? D1[b + r] : 0);
            W.nxt[s] += c;
        }
        swap(W.cur, W.nxt);
    }
    int bs = 0;
    for (int s = 1; s < NS; s++) if (W.cur[s] < W.cur[bs]) bs = s;
    int s = bs;
    for (int t = T; t >= t0; t--) {
        int q1 = s / NS1, q2 = s % NS1;
        S[t][e1].k = q1 / 2 - r; S[t][e1].a = q1 & 1;
        S[t][e2].k = q2 / 2 - r; S[t][e2].a = q2 & 1;
        s = W.par[t][s];
    }
    (void)base;
    return -1;
}

// ---------------------------------------------------------------------------------------------
// cleanup beam search over all 6 emitters
struct CNode { Conf cfg; Bits board; int parent; int key; };

static inline int edist(const Geo& g, const ES& s, int e, int line) {
    int d = abs(SG[e] * s.k - line) + (s.a ? 0 : 1);
    return d < 1 ? 1 : d;
}

static void cleanLB(const Geo& g, const Conf& c, const Bits& wP, const Bits& wK, int& hmax, int& hsum) {
    int ce[6], off[6];
    for (int e = 0; e < 6; e++) { ce[e] = SG[e] * c[e].k; off[e] = c[e].a ? 0 : 1; }
    hmax = 0; hsum = 0;
    auto dist = [&](int e, int line) { int d = abs(ce[e] - line) + off[e]; return d < 1 ? 1 : d; };
    bforeach(wP, [&](int x) {
        int d[6];
        for (int e = 0; e < 6; e++) d[e] = dist(e, g.L[AX[e]][x]);
        int a0 = min(d[0], d[3]), a1 = min(d[1], d[4]), a2 = min(d[2], d[5]);
        int mx = max(a0, max(a1, a2)), mn = min(a0, min(a1, a2));
        int sec = a0 + a1 + a2 - mx - mn;
        int dup = min(max(d[0], d[3]), min(max(d[1], d[4]), max(d[2], d[5])));
        int v = min(sec, dup);
        hmax = max(hmax, v); hsum += v;
    });
    bforeach(wK, [&](int x) {
        int d[6];
        for (int e = 0; e < 6; e++) d[e] = dist(e, g.L[AX[e]][x]);
        int a[3] = {min(d[0], d[3]), min(d[1], d[4]), min(d[2], d[5])};
        int v = max(a[0], max(a[1], a[2]));
        for (int A = 0; A < 3; A++) {
            int du = max(d[A], d[A + 3]);
            int other = min(a[(A + 1) % 3], a[(A + 2) % 3]);
            v = min(v, max(du, other));
        }
        hmax = max(hmax, v); hsum += v;
    });
}

// returns number of ticks appended (path), or -1
static int cleanup(const Geo& g, const Pat& pt, const Conf& c0, const Bits& b0, int maxDepth, int beamW,
                   vector<Conf>& path, Bits* fboard = nullptr) {
    const Bits Z = bz();
    if (bnone(andnot(pt.P, b0)) && bnone(pt.K & b0)) { path.clear(); return 0; }  // (should not happen)
    vector<vector<CNode>> lev;
    { CNode n; n.cfg = c0; n.board = b0; n.parent = -1; n.key = 0; lev.push_back({n}); }
    unordered_map<u64, int> seen;
    vector<CNode> cand;
    auto ckey = [](const Conf& c) { u64 h = 0; for (int e = 0; e < 6; e++) h = h * 64 + (u64)((c[e].k + 20) * 2 + c[e].a); return h; };
    const int nw = (g.nc + 63) / 64;
    const int slackW = 8;
    for (int d = 1; d <= maxDepth; d++) {
        cand.clear(); seen.clear();
        int rem = maxDepth - d;
        int minW = 1 << 20;
        const vector<CNode>& cur = lev.back();
        for (int ni = 0; ni < (int)cur.size(); ni++) {
            const CNode& nd = cur[ni];
            ES oc[6][4]; const Bits* om[6][4]; int no[6];
            for (int e = 0; e < 6; e++) {
                no[e] = 0;
                ES s = nd.cfg[e];
                ES opts[4]; int n = 0;
                opts[n++] = s;
                if (s.k < g.r) { ES q = s; q.k++; opts[n++] = q; }
                if (s.k > -g.r) { ES q = s; q.k--; opts[n++] = q; }
                { ES q = s; q.a ^= 1; opts[n++] = q; }
                for (int o = 0; o < n; o++) {
                    oc[e][no[e]] = opts[o];
                    int l = lineOf(g, e, opts[o]);
                    om[e][no[e]] = l == NONE ? &Z : &g.line[AX[e]][l + g.r];
                    no[e]++;
                }
            }
            int ix[6];
            Bits pa[7], pb[7], pc[7];
            pa[0] = Z; pb[0] = Z; pc[0] = Z;
            auto acc = [&](int lv, const Bits* m) {
                for (int w = 0; w < nw; w++) {
                    u64 x = m->w[w];
                    pc[lv + 1].w[w] = pc[lv].w[w] | (pb[lv].w[w] & x);
                    pb[lv + 1].w[w] = pb[lv].w[w] | (pa[lv].w[w] & x);
                    pa[lv + 1].w[w] = pa[lv].w[w] | x;
                }
            };
            for (ix[0] = 0; ix[0] < no[0]; ix[0]++) { acc(0, om[0][ix[0]]);
            for (ix[1] = 0; ix[1] < no[1]; ix[1]++) { acc(1, om[1][ix[1]]);
            for (ix[2] = 0; ix[2] < no[2]; ix[2]++) { acc(2, om[2][ix[2]]);
            for (ix[3] = 0; ix[3] < no[3]; ix[3]++) { acc(3, om[3][ix[3]]);
            for (ix[4] = 0; ix[4] < no[4]; ix[4]++) { acc(4, om[4][ix[4]]);
            for (ix[5] = 0; ix[5] < no[5]; ix[5]++) { acc(5, om[5][ix[5]]);
                Bits nb = Z, wP = Z, wK = Z;
                int Wc = 0;
                for (int w = 0; w < nw; w++) {
                    u64 cc = pc[6].w[w];
                    nb.w[w] = (nd.board.w[w] | (pb[6].w[w] & ~cc)) & ~cc;
                    wP.w[w] = pt.P.w[w] & ~nb.w[w];
                    wK.w[w] = pt.K.w[w] & nb.w[w];
                    Wc += __builtin_popcountll(wP.w[w]) + __builtin_popcountll(wK.w[w]);
                }
                if (Wc > minW + slackW) continue;
                minW = min(minW, Wc);
                Conf c;
                for (int e = 0; e < 6; e++) c[e] = oc[e][ix[e]];
                if (Wc == 0) {
                    path.clear(); path.push_back(c);
                    int pi = ni;
                    for (int dd = d - 1; dd >= 1; dd--) { path.push_back(lev[dd][pi].cfg); pi = lev[dd][pi].parent; }
                    reverse(path.begin(), path.end());
                    if (fboard) *fboard = andnot(nb, pt.P);
                    return d;
                }
                if (rem <= 0) continue;
                int hm, hs; cleanLB(g, c, wP, wK, hm, hs);
                if (hm > rem) continue;
                int key = (4 * hm + 2 * Wc) * 1024 + min(hs, 1023);
                u64 ck = ckey(c);
                auto it = seen.find(ck);
                if (it != seen.end()) {
                    CNode& s = cand[it->second];
                    if (s.key <= key) continue;
                    s.cfg = c; s.board = nb; s.parent = ni; s.key = key;
                } else {
                    seen[ck] = cand.size();
                    CNode s; s.cfg = c; s.board = nb; s.parent = ni; s.key = key;
                    cand.push_back(s);
                }
            }}}}}}
        }
        if (cand.empty()) return -1;
        if ((int)cand.size() > beamW) {
            nth_element(cand.begin(), cand.begin() + beamW, cand.end(), [](const CNode& a, const CNode& b) { return a.key < b.key; });
            cand.resize(beamW);
        }
        lev.push_back(cand);
    }
    return -1;
}

// ---------------------------------------------------------------------------------------------
struct Prm {
    int restarts = 6, rounds = 6, noise = 3;
    int trims = 2, delays = 1;
    int cleanW = 200, cleanD = 14;
    int WP = 10, WK = 10;
    int subCands = 16, maxClean = 6, maxCleanSub = 3, reweights = 3;
    double defW = 0.5;
};

// build a model for plan builds; tick ranges given by fires; builds[b] = list of fires; bEnd[b] = acceptance tick
static void setupModel(Model& M, int t0, const vector<vector<Fire>>& builds, const vector<int>& bEnd) {
    M.t0 = t0; M.T = bEnd.back();
    M.bOf.assign(M.T + 1, 0);
    int b = 0;
    for (int t = t0; t <= M.T; t++) { while (t > bEnd[b]) b++; M.bOf[t] = b; }
    M.fx.assign(builds.size(), vector<short>(M.g.nc, -1));
    for (size_t bb = 0; bb < builds.size(); bb++)
        for (auto& f : builds[bb])
            for (int x : M.lc[f.ax][f.line + M.g.r])
                M.fx[bb][x] = max<short>(M.fx[bb][x], (short)f.t);
}

struct Result { bool ok = false; Prog prog; Conf endCfg; Bits endBoard; int accDone = 0; };

// Exact evaluation of a planned schedule from t0.
//  ok=true if >= need acceptances happen (len = tick of the need-th acceptance).
//  else: failBuild = first build whose acceptance did not happen at its planned tick (boards at that tick
//  in failBoard, defects counted); defects = -1 if the deviation was an unexpected acceptance.
struct EvalRes { bool ok = false; int len = 0; int okBuilds = 0; int failBuild = -1; int defects = 1 << 30; Bits failBoard; };

static EvalRes evalSched(const Geo& g, const Pat& pt, const Sched& S, int t0, const Bits& board0, const vector<int>& bEnd, int need) {
    EvalRes R;
    vector<int> accT; Bits board = board0;
    int T = bEnd.back();
    const Bits Z = bz();
    size_t b = 0;
    for (int t = t0; t <= T; t++) {
        const Bits* m[6]; confMasks(g, S[t], m, Z);
        tickBoard(board, m);
        bool acc = bnone(andnot(pt.P, board)) && bnone(pt.K & board);
        if (acc) {
            board = andnot(board, pt.P); accT.push_back(t);
            if ((int)accT.size() >= need) { R.ok = true; R.len = t; R.okBuilds = need; return R; }
            if (t != bEnd[b]) { R.failBuild = (int)b; R.defects = -1; return R; }
            b++; R.okBuilds++;
            continue;
        }
        if (t == bEnd[b]) {
            R.failBuild = (int)b; R.failBoard = board;
            R.defects = bpop(andnot(pt.P, board)) + bpop(pt.K & board);
            return R;
        }
    }
    return R;
}

struct PlanEval { Sched S; EvalRes ev; long long mcost = 0; bool valid = false; };

static bool betterEval(const EvalRes& a, const EvalRes& b) {  // is a better than b
    if (a.ok != b.ok) return a.ok;
    if (a.ok) return a.len < b.len;
    if (a.okBuilds != b.okBuilds) return a.okBuilds > b.okBuilds;
    if ((a.defects < 0) != (b.defects < 0)) return a.defects >= 0;
    return a.defects < b.defects;
}

static void optimisePlan(const Geo& g, const Pat& pt, const Prm& prm, const Sched& prefix, int t0, const Bits& board0,
                         const vector<vector<Fire>>& builds, const vector<int>& bEnd, int need, int mainAx, PlanEval& best, mt19937& rng) {
    Model M(g, pt); M.WP = prm.WP; M.WK = prm.WK;
    setupModel(M, t0, builds, bEnd);
    DPWork W;
    for (int rs = 0; rs < prm.restarts; rs++) {
        M.initWeights();
        Sched S = prefix;
        S.resize(M.T + 1);
        for (int t = t0; t <= M.T; t++) S[t] = S[t - 1];
        int o[3] = {mainAx, (mainAx + 1) % 3, (mainAx + 2) % 3};
        if (rs & 1) swap(o[1], o[2]);
        bool noisy = rs >= 2;
        for (int rw = 0; rw <= prm.reweights; rw++) {
            long long prev = LLONG_MAX;
            for (int rd = 0; rd < prm.rounds; rd++) {
                bool nz = noisy && rw == 0 && rd < (prm.rounds + 1) / 2;
                for (int q = 0; q < 3; q++) dpAxis(M, S, o[q], nz ? &rng : nullptr, nz ? prm.noise : 0, W);
                long long c = modelCost(M, S);
                if (!nz && c >= prev) break;
                prev = c;
            }
            EvalRes ev = evalSched(g, pt, S, t0, board0, bEnd, need);
            if (VERB >= 3) fprintf(stderr, "    restart %d rw %d: model %lld ok %d len %d okB %d def %d\n", rs, rw, prev, ev.ok, ev.len, ev.okBuilds, ev.defects);
            if (!best.valid || betterEval(ev, best.ev)) { best.S = S; best.ev = ev; best.mcost = prev; best.valid = true; }
            if (ev.ok || ev.defects < 0) break;
            Bits bad = andnot(pt.P, ev.failBoard) | (pt.K & ev.failBoard);
            bforeach(bad, [&](int x) { M.wt[x] += (M.typ[x] == 1 ? M.WP : M.WK); });
        }
        if (best.ev.ok) break;
    }
}

// ---------------------------------------------------------------------------------------------
// plan generation
static void axisRange(const Geo& g, const Pat& pt, int A, int& lo, int& hi) {
    lo = 99; hi = -99;
    bforeach(pt.P, [&](int x) { lo = min(lo, g.L[A][x]); hi = max(hi, g.L[A][x]); });
}

// sweep plan: builds alternate direction over lines s..e of axis A, first fire at tick t1
static void sweepPlan(int A, int s, int e, int t1, int nb, vector<vector<Fire>>& builds, vector<int>& bEnd) {
    builds.clear(); bEnd.clear();
    int t = t1;
    for (int b = 0; b < nb; b++) {
        int from = (b % 2 == 0) ? s : e, to = (b % 2 == 0) ? e : s;
        int dir = to >= from ? 1 : -1;
        vector<Fire> fs;
        for (int c = from;; c += dir) { fs.push_back({A, c, t}); if (c == to) break; t++; }
        builds.push_back(fs); bEnd.push_back(t);
        t++;
    }
}

static int approachTicks(const Geo& g, const Conf& c, int A, int line) {
    int t = 1;
    for (int e : {A, A + 3}) {
        int d = abs(c[e].k - SG[e] * line) + (c[e].a ? 0 : 1);
        t = max(t, d);
    }
    return t;
}

static bool solveFrom(const Geo& g, const Pat& pt, const Prm& prm, const Sched& prefix, const Bits& board, int need, int limit,
                      Sched& out, mt19937& rng, int depth) {
    int t0 = (int)prefix.size();
    const Conf& c0 = prefix.back();
    struct Cand { int A, s, e, t1, T; vector<vector<Fire>> builds; vector<int> bEnd; PlanEval pe; double est; };
    vector<Cand> cands;
    for (int A = 0; A < 3; A++) {
        int lo, hi; axisRange(g, pt, A, lo, hi);
        int tr = depth == 0 ? prm.trims : min(prm.trims, 1);
        for (int ts = 0; ts <= tr; ts++)
            for (int te = 0; te <= tr; te++) {
                int a = lo + ts, b = hi - te;
                if (a > b) continue;
                for (int dir = 0; dir < 2; dir++) {
                    int s = dir ? b : a, e = dir ? a : b;
                    int ap = approachTicks(g, c0, A, s);
                    for (int dl = 0; dl <= prm.delays; dl++) {
                        Cand cd; cd.A = A; cd.s = s; cd.e = e;
                        cd.t1 = t0 - 1 + ap + dl;
                        int n = abs(e - s) + 1;
                        cd.T = cd.t1 + need * n - 1;
                        if (cd.T > limit) continue;
                        cands.push_back(cd);
                    }
                }
            }
    }
    sort(cands.begin(), cands.end(), [](const Cand& x, const Cand& y) { return x.T < y.T; });
    if (depth > 0 && (int)cands.size() > prm.subCands) cands.resize(prm.subCands);
    bool found = false;
    int bestLen = limit + 1;
    for (auto& cd : cands) {
        if (cd.T >= bestLen) { cd.est = 1e9; continue; }
        sweepPlan(cd.A, cd.s, cd.e, cd.t1, need, cd.builds, cd.bEnd);
        optimisePlan(g, pt, prm, prefix, t0, board, cd.builds, cd.bEnd, need, cd.A, cd.pe, rng);
        const EvalRes& ev = cd.pe.ev;
        if (ev.ok) {
            cd.est = 1e9;
            if (ev.len < bestLen) {
                bestLen = ev.len; out = cd.pe.S; out.resize(ev.len + 1); found = true;
                if (VERB) fprintf(stderr, "  [d%d t0=%d need=%d] A=%d %d..%d t1=%d -> clean %d\n", depth, t0, need, cd.A, cd.s, cd.e, cd.t1, ev.len);
            }
            continue;
        }
        if (ev.defects < 0) { cd.est = 1e9; continue; }
        int fbEnd = cd.bEnd[ev.failBuild];
        int remB = need - ev.failBuild - 1;
        cd.est = fbEnd + 1 + ev.defects * prm.defW + (remB ? cd.T - fbEnd : 0);
    }
    // cleanup for the most promising
    vector<int> ord;
    for (int i = 0; i < (int)cands.size(); i++) if (cands[i].est < 1e8) ord.push_back(i);
    sort(ord.begin(), ord.end(), [&](int a, int b) { return cands[a].est < cands[b].est; });
    int tried = 0;
    for (int i : ord) {
        if (tried >= (depth == 0 ? prm.maxClean : prm.maxCleanSub)) break;
        Cand& cd = cands[i];
        const EvalRes& ev = cd.pe.ev;
        int fbEnd = cd.bEnd[ev.failBuild];
        int remB = need - ev.failBuild - 1;
        int maxD = min(prm.cleanD, bestLen - 1 - fbEnd - remB);
        if (maxD < 1) continue;
        tried++;
        vector<Conf> path; Bits fboard;
        auto tc = chrono::steady_clock::now();
        int d = cleanup(g, pt, cd.pe.S[fbEnd], ev.failBoard, maxD, prm.cleanW, path, &fboard);
        if (VERB >= 2) fprintf(stderr, "  [d%d t0=%d need=%d] A=%d %d..%d t1=%d T=%d okB=%d def=%d est=%.1f cleanup %d (%.2fs)\n", depth, t0, need, cd.A, cd.s, cd.e,
                               cd.t1, cd.T, ev.okBuilds, ev.defects, cd.est, d, chrono::duration<double>(chrono::steady_clock::now() - tc).count());
        if (d < 0) continue;
        Sched np(cd.pe.S.begin(), cd.pe.S.begin() + fbEnd + 1);
        for (auto& c : path) np.push_back(c);
        int lenNow = (int)np.size() - 1;
        if (remB == 0) {
            if (lenNow < bestLen) { bestLen = lenNow; out = np; found = true;
                if (VERB) fprintf(stderr, "  [d%d t0=%d need=%d] A=%d %d..%d t1=%d -> %d (cleanup %d)\n", depth, t0, need, cd.A, cd.s, cd.e, cd.t1, lenNow, d); }
            continue;
        }
        Sched sub;
        if (solveFrom(g, pt, prm, np, fboard, remB, bestLen - 1, sub, rng, depth + 1)) {
            int L = (int)sub.size() - 1;
            if (L < bestLen) { bestLen = L; out = sub; found = true;
                if (VERB) fprintf(stderr, "  [d%d t0=%d need=%d] A=%d %d..%d -> %d (via recursion)\n", depth, t0, need, cd.A, cd.s, cd.e, bestLen); }
        }
    }
    return found;
}

static vector<int> parseList(const string& s) {
    vector<int> v; string cur;
    for (char ch : s + ",") { if (ch == ',') { if (!cur.empty()) v.push_back(stoi(cur)); cur.clear(); } else cur += ch; }
    return v;
}

static Prog readProg(const string& path) {
    auto v = parseOutput(path);
    if (v.empty()) return Prog();
    return v[0];
}

int main(int argc, char** argv) {
    if (argc < 4) { fprintf(stderr, "usage: lines <in> <best_out> <res_dir> [opts]\n"); return 1; }
    string inp = argv[1], bestp = argv[2], resdir = argv[3];
    Prm prm; vector<int> only; int threads = 2; int minP = 0;
    for (int i = 4; i + 1 < argc; i += 2) {
        string o = argv[i], v = argv[i + 1];
        if (o == "-only") only = parseList(v);
        else if (o == "-v") VERB = stoi(v);
        else if (o == "-t") threads = stoi(v);
        else if (o == "-rs") prm.restarts = stoi(v);
        else if (o == "-rd") prm.rounds = stoi(v);
        else if (o == "-noise") prm.noise = stoi(v);
        else if (o == "-tr") prm.trims = stoi(v);
        else if (o == "-dl") prm.delays = stoi(v);
        else if (o == "-cw") prm.cleanW = stoi(v);
        else if (o == "-cd") prm.cleanD = stoi(v);
        else if (o == "-WP") prm.WP = stoi(v);
        else if (o == "-WK") prm.WK = stoi(v);
        else if (o == "-minP") minP = stoi(v);
        else if (o == "-sc") prm.subCands = stoi(v);
        else if (o == "-rw") prm.reweights = stoi(v);
        else if (o == "-mc") prm.maxClean = stoi(v);
        else if (o == "-mcs") prm.maxCleanSub = stoi(v);
        else if (o == "-defW") prm.defW = atof(v.c_str());
    }
    auto pats = parseInput(inp);
    auto base = parseOutput(bestp);
    map<int, Geo> geos;
    for (auto& p : pats) if (!geos.count(p.r)) geos[p.r] = makeGeo(p.r);
    for (auto& p : pats) initPat(p, geos[p.r]);
    int n = pats.size();
    vector<int> todo;
    if (only.empty()) { for (int i = 0; i < n; i++) if (pats[i].np >= minP) todo.push_back(i); } else todo = only;
    sort(todo.begin(), todo.end(), [&](int a, int b) { return base[a].size() > base[b].size(); });
    int wins = 0, gain = 0;
#ifdef _OPENMP
    omp_set_num_threads(threads);
#endif
    #pragma omp parallel for schedule(dynamic, 1)
    for (int ii = 0; ii < (int)todo.size(); ii++) {
        int i = todo[ii];
        const Geo& g = geos[pats[i].r];
        string fn = resdir + "/" + to_string(i) + ".txt";
        int limit = base[i].size() - 1;
        Prog old = readProg(fn);
        if (!old.empty() && simulate(g, pats[i], old) >= pats[i].cnt) limit = min(limit, (int)old.size() - 1);
        string log;
        auto t0 = chrono::steady_clock::now();
        mt19937 rng(12345 + i);
        Sched pre(1, startConf()), out;
        bool ok = solveFrom(g, pats[i], prm, pre, bz(), pats[i].cnt, limit, out, rng, 0);
        Prog p;
        if (ok) p = toProg(out);
        bool win = ok && (int)p.size() <= limit && simulate(g, pats[i], p) >= pats[i].cnt;
        if (ok && !win && (int)p.size() <= limit) fprintf(stderr, "!! pattern %d: produced program fails exact check\n", i);
        if (win) { vector<Prog> v = {p}; writeOutput(fn, v); }
        #pragma omp critical
        {
            if (win) { wins++; gain += base[i].size() - p.size(); }
            fprintf(stderr, "pat %d r=%d cnt=%d P=%d K=%d best=%zu mine=%zu %s (%.1fs)\n", i, pats[i].r, pats[i].cnt, pats[i].np, pats[i].nk,
                    base[i].size(), ok ? p.size() : (size_t)0, win ? "WIN" : "", chrono::duration<double>(chrono::steady_clock::now() - t0).count());
            fflush(stderr);
        }
    }
    printf("wins %d gain %d\n", wins, gain);
    return 0;
}
