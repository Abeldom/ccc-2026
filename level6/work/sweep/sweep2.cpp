// sweep2: duplicate-pair line sweep with beam-searched cutters + axis-aware cleanup beam.
// usage: sweep2 <in> <best_out> <res_dir> [opts]
//   opts: -Bs n -Bc n -M n -Dc n -mc n -tr n -only i,j,.. -v n -min_gain n -order big
// For every pattern it tries to beat the program length in <best_out>; wins are written to
// <res_dir>/<idx>.txt (one program), only if shorter than an existing file there.
#include "common.h"

struct Prm {
    int Bs = 400, Bc = 150, M = 30, Dc = 16, maxCleanups = 8;
    int wP = 2, wK = 3;
    int la1 = 3, la2 = 2, la3 = 1;
    int trims = 2;
    int keyLB = 4, keyW = 2;
};
static int VERB = 0;

struct SNode {
    Cfg cfg;
    Bits board;
    int parent;
    int score;
    Ins ins;
};

struct BuildRes {
    bool ok = false;
    Prog prog;
    Cfg cfg;
    Bits board;
};

static inline int edist(const Geo& g, const Cfg& c, int e, int line) {
    // ticks needed until emitter e fires on 'line' (>=1)
    int d = abs(SG[e] * c.k[e] - line) + (c.on(e) ? 0 : 1);
    return d < 1 ? 1 : d;
}

static void dumpBoard(const Geo& g, const Pat& pt, const Bits& b) {
    vector<string> grid = pt.grid;
    for (int x = 0; x < g.nc; x++) {
        char& ch = grid[g.cy[x]][g.cx[x]];
        bool P = bget(pt.P, x), K = bget(pt.K, x), m = bget(b, x);
        ch = P ? (m ? 'M' : 'm') : K ? (m ? 'k' : 'o') : (m ? '#' : '.');
    }
    for (auto& r : grid) fprintf(stderr, "    %s\n", r.c_str());
}

// ---------------- cleanup -----------------
struct Cleaner {
    const Geo& g; const Pat& pt; const Prm& prm;
    Cleaner(const Geo& g_, const Pat& p_, const Prm& pr) : g(g_), pt(p_), prm(pr) {}

    inline void lb(const Cfg& c, const Bits& wP, const Bits& wK, int& hmax, int& hsum) const {
        int ce[6], off[6];
        for (int e = 0; e < 6; e++) { ce[e] = SG[e] * c.k[e]; off[e] = c.on(e) ? 0 : 1; }
        hmax = 0; hsum = 0;
        auto dist = [&](int e, int line) { int d = abs(ce[e] - line) + off[e]; return d < 1 ? 1 : d; };
        bforeach(wP, [&](int x) {
            int d[6];
            for (int e = 0; e < 6; e++) d[e] = dist(e, g.L[AX[e]][x]);
            int a0 = min(d[0], d[3]), a1 = min(d[1], d[4]), a2 = min(d[2], d[5]);
            // second smallest of a0,a1,a2
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

    struct CNode { Cfg cfg; Bits board; int parent; int key; Ins ins; };

    bool run(const vector<SNode>& init, int maxDepth, vector<Ins>& path, int& initIdx, Cfg& fcfg, Bits& fboard) const {
        const Bits Z = bz();
        vector<vector<CNode>> lev;
        {
            vector<CNode> l0;
            for (auto& s : init) { CNode c; c.cfg = s.cfg; c.board = s.board; c.parent = -1; c.key = 0; c.ins = s.ins; l0.push_back(c); }
            lev.push_back(l0);
        }
        unordered_map<u64, int> seen;
        vector<CNode> cand;
        for (int d = 1; d <= maxDepth; d++) {
            cand.clear(); seen.clear();
            const vector<CNode>& cur = lev.back();
            int rem = maxDepth - d;  // ticks left after this one
            for (int ni = 0; ni < (int)cur.size(); ni++) {
                const CNode& nd = cur[ni];
                Cfg oc[6][4]; char oi[6][4]; const Bits* om[6][4]; int no[6];
                for (int e = 0; e < 6; e++) {
                    no[e] = 0;
                    const char opts[4] = {'-', 'F', 'B', 'T'};
                    for (int o = 0; o < 4; o++) {
                        Cfg c = nd.cfg;
                        if (!applyIns(g, c, e, opts[o])) continue;
                        oc[e][no[e]] = c; oi[e][no[e]] = opts[o];
                        om[e][no[e]] = &emask(g, e, c.k[e], c.on(e), Z);
                        no[e]++;
                    }
                }
                int ix[6];
                for (ix[0] = 0; ix[0] < no[0]; ix[0]++)
                for (ix[1] = 0; ix[1] < no[1]; ix[1]++)
                for (ix[2] = 0; ix[2] < no[2]; ix[2]++)
                for (ix[3] = 0; ix[3] < no[3]; ix[3]++)
                for (ix[4] = 0; ix[4] < no[4]; ix[4]++)
                for (ix[5] = 0; ix[5] < no[5]; ix[5]++) {
                    const Bits* m[6];
                    Cfg c; c.act = 0;
                    for (int e = 0; e < 6; e++) {
                        m[e] = om[e][ix[e]];
                        c.k[e] = oc[e][ix[e]].k[e];
                        c.act |= oc[e][ix[e]].act & (1 << e);
                    }
                    Bits nb = nd.board;
                    tickBoard(nb, m);
                    Bits wP = andnot(pt.P, nb), wK = pt.K & nb;
                    Ins ins; for (int e = 0; e < 6; e++) ins[e] = oi[e][ix[e]];
                    if (bnone(wP) && bnone(wK)) {
                        path.clear(); path.push_back(ins);
                        int pi = ni;
                        for (int dd = d - 1; dd >= 1; dd--) { path.push_back(lev[dd][pi].ins); pi = lev[dd][pi].parent; }
                        reverse(path.begin(), path.end());
                        initIdx = pi;
                        fcfg = c; fboard = andnot(nb, pt.P);
                        return true;
                    }
                    if (rem <= 0) continue;
                    int hm, hs; lb(c, wP, wK, hm, hs);
                    if (hm > rem) continue;
                    int W = bpop(wP) + bpop(wK);
                    int key = (prm.keyLB * hm + prm.keyW * W) * 1024 + min(hs, 1023);
                    u64 ck = c.key();
                    auto it = seen.find(ck);
                    if (it != seen.end()) {
                        CNode& s = cand[it->second];
                        if (s.key <= key) continue;
                        s.cfg = c; s.board = nb; s.parent = ni; s.key = key; s.ins = ins;
                    } else {
                        seen[ck] = cand.size();
                        CNode s; s.cfg = c; s.board = nb; s.parent = ni; s.key = key; s.ins = ins;
                        cand.push_back(s);
                    }
                }
            }
            if (cand.empty()) return false;
            if ((int)cand.size() > prm.Bc) {
                nth_element(cand.begin(), cand.begin() + prm.Bc, cand.end(), [](const CNode& a, const CNode& b) { return a.key < b.key; });
                cand.resize(prm.Bc);
            }
            if (VERB >= 4) {
                int bi = 0;
                for (int i = 1; i < (int)cand.size(); i++) if (cand[i].key < cand[bi].key) bi = i;
                const CNode& b = cand[bi];
                Bits wP = andnot(pt.P, b.board), wK = pt.K & b.board;
                int hm, hs; lb(b.cfg, wP, wK, hm, hs);
                fprintf(stderr, "    clean d=%d n=%zu best: missP %d dirtyK %d hm %d hs %d ins %c%c%c%c%c%c cfg:", d, cand.size(), bpop(wP), bpop(wK), hm, hs, b.ins[0], b.ins[1], b.ins[2], b.ins[3], b.ins[4], b.ins[5]);
                for (int e = 0; e < 6; e++) fprintf(stderr, " %d%s", SG[e] * b.cfg.k[e], b.cfg.on(e) ? "*" : "");
                fprintf(stderr, "\n");
                if (VERB >= 5) dumpBoard(g, pt, b.board);
            }
            lev.push_back(cand);
        }
        return false;
    }
};

// ---------------- sweep -----------------
struct Plan { int A; vector<int> rows; int dir = 1; int pre = 0, post = 0, painter = 0, postOff = 1; };

struct SweepRun {
    Plan plan;
    bool done = false;
    BuildRes res;
    int Tend = 0;
    vector<vector<SNode>> lev;
    int est = 1 << 30;
};

// max matching size between cells (each with 4-bit mask of usable cutters) and 4 cutters
static inline int matchSize(const int* masks, int n) {
    // DP over set of used cutters
    unsigned reach = 1;  // bit s set if used-set s achievable
    for (int i = 0; i < n; i++) {
        unsigned nr = reach;
        for (int s = 0; s < 16; s++) if ((reach >> s) & 1) {
            int fr = masks[i] & ~s;
            while (fr) { int b = fr & -fr; nr |= 1u << (s | b); fr ^= b; }
        }
        reach = nr;
    }
    int best = 0;
    for (int s = 0; s < 16; s++) if ((reach >> s) & 1) best = max(best, __builtin_popcount(s));
    return best;
}

// instruction string moving emitter from (k0,on0) to (k1, active); if on0 and the path fires on lines
// outside [lo,hi] it switches off first.
static string emPath(const Geo& g, int e, int k0, bool on0, int k1, int lo, int hi) {
    int d = abs(k1 - k0);
    char mv = k1 > k0 ? 'F' : 'B';
    if (!on0) return string(d, mv) + "T";
    bool inrange = true;
    for (int j = 0; j < d; j++) {
        int kk = k0 + (k1 > k0 ? j : -j);
        int line = SG[e] * kk;
        if (abs(kk) < g.r && (line < lo || line > hi)) inrange = false;
    }
    if (inrange) return string(d, mv);
    return "T" + string(d, mv) + "T";
}

static void sweepRun(SweepRun& R, const Geo& g, const Pat& pt, const Cfg& c0, const Bits& b0, const Prm& prm, int bestLen) {
    const Bits Z = bz();
    const Plan& pl = R.plan;
    int A = pl.A;
    const vector<int>& rows = pl.rows;
    int em[2] = {A, A + 3};
    int nrows = rows.size();
    int lo = *min_element(rows.begin(), rows.end()), hi = *max_element(rows.begin(), rows.end());
    int dir = nrows > 1 ? (rows[1] > rows[0] ? 1 : -1) : pl.dir;
    // ---- schedule for the two dup emitters
    vector<string> sch(2);
    int painter = pl.painter;  // 0 or 1
    if (pl.pre > 0) {
        int preRow = rows[0] - dir;
        int e = em[painter];
        string s1 = emPath(g, e, c0.k[e], c0.on(e), SG[e] * preRow, lo - (dir > 0), hi + (dir < 0));
        string s = s1.empty() ? string(pl.pre, '-') : s1 + string(pl.pre - 1, '-');
        s += (SG[e] * dir > 0) ? 'F' : 'B';
        sch[painter] = s;
        int o = em[1 - painter];
        sch[1 - painter] = emPath(g, o, c0.k[o], c0.on(o), SG[o] * rows[0], lo, hi);
    } else {
        for (int i = 0; i < 2; i++) { int e = em[i]; sch[i] = emPath(g, e, c0.k[e], c0.on(e), SG[e] * rows[0], lo, hi); }
    }
    int T1 = max(1, (int)max(sch[0].size(), sch[1].size()));
    for (int i = 0; i < 2; i++) sch[i] = string(T1 - sch[i].size(), '-') + sch[i];
    int Tsw = T1 + nrows - 1;
    int Tend = Tsw + pl.post;
    R.Tend = Tend;
    if (Tend >= bestLen) return;
    vector<array<char, 2>> dup(Tend + 1);
    for (int t = 1; t <= T1; t++) for (int i = 0; i < 2; i++) dup[t][i] = sch[i][t - 1];
    for (int j = 1; j < nrows; j++) {
        int dr = rows[j] - rows[j - 1];
        for (int i = 0; i < 2; i++) dup[T1 + j][i] = (SG[em[i]] * dr > 0) ? 'F' : 'B';
    }
    for (int j = 1; j <= pl.post; j++) {
        int pe = em[painter], oe = em[1 - painter];
        dup[Tsw + j][painter] = j == 1 ? ((SG[pe] * dir > 0) ? 'F' : 'B') : '-';
        dup[Tsw + j][1 - painter] = (j == 1 && pl.postOff) ? 'T' : '-';
        (void)oe;
    }
    // ---- permanent region per tick
    vector<Bits> perm(Tend + 1);
    Bits base = bz();
    for (int x = 0; x < g.nc; x++) if (g.L[A][x] < lo || g.L[A][x] > hi) bset(base, x);
    for (int t = 0; t <= Tend; t++) {
        perm[t] = base;
        for (int j = 0; j < nrows && T1 + j <= t; j++) perm[t] = perm[t] | g.line[A][rows[j] + g.r];
    }
    // ---- per tick requirement lists: K cells that must be cut (sweep rows), P cells to paint
    vector<vector<int>> needK(Tend + 2), needP(Tend + 2);
    for (int j = 0; j < nrows; j++) bforeach(pt.K & g.line[A][rows[j] + g.r], [&](int x) { needK[T1 + j].push_back(x); });
    if (pl.pre > 0) {
        int preRow = rows[0] - dir;
        for (int t = T1 - pl.pre; t < T1; t++) if (t >= 1) bforeach(pt.P & g.line[A][preRow + g.r], [&](int x) { needP[t].push_back(x); });
    }
    if (pl.post > 0) {
        int postRow = rows.back() + dir;
        for (int t = Tsw + 1; t <= Tend; t++) bforeach(pt.P & g.line[A][postRow + g.r], [&](int x) { needP[t].push_back(x); });
    }
    int cut[4], ncut = 0;
    for (int e = 0; e < 6; e++) if (AX[e] != A) cut[ncut++] = e;
    vector<vector<SNode>>& lev = R.lev;
    SNode s0; s0.cfg = c0; s0.board = b0; s0.parent = -1; s0.score = 0;
    lev.push_back({s0});
    unordered_map<u64, int> seen;
    vector<SNode> cand;
    const int law[4] = {0, prm.la1, prm.la2, prm.la3};
    for (int t = 1; t <= Tend; t++) {
        cand.clear(); seen.clear();
        const vector<SNode>& cur = lev.back();
        for (int ni = 0; ni < (int)cur.size(); ni++) {
            const SNode& nd = cur[ni];
            Cfg cb = nd.cfg;
            applyIns(g, cb, em[0], dup[t][0]);
            applyIns(g, cb, em[1], dup[t][1]);
            const Bits* m[6];
            m[em[0]] = &emask(g, em[0], cb.k[em[0]], cb.on(em[0]), Z);
            m[em[1]] = &emask(g, em[1], cb.k[em[1]], cb.on(em[1]), Z);
            Cfg oc[4][4]; char oi[4][4]; const Bits* om[4][4]; int no[4];
            for (int q = 0; q < 4; q++) {
                int e = cut[q]; no[q] = 0;
                const char opts[4] = {'-', 'F', 'B', 'T'};
                for (int o = 0; o < 4; o++) {
                    Cfg c = cb;
                    if (!applyIns(g, c, e, opts[o])) continue;
                    oc[q][no[q]] = c; oi[q][no[q]] = opts[o];
                    om[q][no[q]] = &emask(g, e, c.k[e], c.on(e), Z);
                    no[q]++;
                }
            }
            int ix[4];
            for (ix[0] = 0; ix[0] < no[0]; ix[0]++)
            for (ix[1] = 0; ix[1] < no[1]; ix[1]++)
            for (ix[2] = 0; ix[2] < no[2]; ix[2]++)
            for (ix[3] = 0; ix[3] < no[3]; ix[3]++) {
                Cfg c = cb;
                Ins ins; ins[em[0]] = dup[t][0]; ins[em[1]] = dup[t][1];
                for (int q = 0; q < 4; q++) {
                    int e = cut[q];
                    m[e] = om[q][ix[q]];
                    c.k[e] = oc[q][ix[q]].k[e];
                    c.act = (c.act & ~(1 << e)) | (oc[q][ix[q]].act & (1 << e));
                    ins[e] = oi[q][ix[q]];
                }
                Bits nb = nd.board;
                tickBoard(nb, m);
                Bits wP = andnot(pt.P, nb), wK = pt.K & nb;
                if (bnone(wP) && bnone(wK)) {
                    Prog p; p.push_back(ins);
                    int pi = ni;
                    for (int tt = t - 1; tt >= 1; tt--) { p.push_back(lev[tt][pi].ins); pi = lev[tt][pi].parent; }
                    reverse(p.begin(), p.end());
                    R.res.ok = true; R.res.prog = p; R.res.cfg = c; R.res.board = andnot(nb, pt.P);
                    R.done = true; R.est = t;
                    return;
                }
                int score = prm.wP * bpop(wP & perm[t]) + prm.wK * bpop(wK & perm[t]);
                score *= 8;
                for (int dl = 1; dl <= 3 && t + dl <= Tend; dl++) {
                    const vector<int>& vk = needK[t + dl];
                    const vector<int>& vp = needP[t + dl];
                    if (vk.empty() && vp.empty()) continue;
                    int masks[64]; int n = 0;
                    for (int x : vk) {
                        int mk = 0;
                        for (int q = 0; q < 4; q++) {
                            int e = cut[q];
                            int dist = abs(SG[e] * c.k[e] - g.L[AX[e]][x]) + (c.on(e) ? 0 : 1);
                            if (dist <= dl) mk |= 1 << q;
                        }
                        if (n < 64) masks[n++] = mk;
                    }
                    int ms = matchSize(masks, n);
                    score += law[dl] * (n - ms);
                    if (!vp.empty()) {
                        int np = 0;
                        for (int x : vp) {
                            if (bget(nb, x)) continue;
                            int mk = 0;
                            for (int q = 0; q < 4; q++) {
                                int e = cut[q];
                                int dist = abs(SG[e] * c.k[e] - g.L[AX[e]][x]) + (c.on(e) ? 0 : 1);
                                if (dist <= dl) mk |= 1 << q;
                            }
                            if (n + np < 64) masks[n + np++] = mk;
                        }
                        int ms2 = matchSize(masks, n + np);
                        int want = min(np, 4 - ms);
                        if (want > 0) score += law[dl] * max(0, want - (ms2 - ms));
                    }
                }
                u64 ck = c.key();
                auto it = seen.find(ck);
                if (it != seen.end()) {
                    SNode& s = cand[it->second];
                    if (s.score <= score) continue;
                    s.cfg = c; s.board = nb; s.parent = ni; s.score = score; s.ins = ins;
                } else {
                    seen[ck] = cand.size();
                    SNode s; s.cfg = c; s.board = nb; s.parent = ni; s.score = score; s.ins = ins;
                    cand.push_back(s);
                }
            }
        }
        if ((int)cand.size() > prm.Bs) {
            nth_element(cand.begin(), cand.begin() + prm.Bs, cand.end(), [](const SNode& a, const SNode& b) { return a.score < b.score; });
            cand.resize(prm.Bs);
        }
        lev.push_back(cand);
    }
    Cleaner cl(g, pt, prm);
    int bs = 1 << 30;
    for (auto& nd : lev.back()) bs = min(bs, nd.score);
    for (auto& nd : lev.back()) {
        if (nd.score > bs + 16) continue;
        Bits wP = andnot(pt.P, nd.board), wK = pt.K & nd.board;
        int hm, hs; cl.lb(nd.cfg, wP, wK, hm, hs);
        int W = bpop(wP) + bpop(wK);
        R.est = min(R.est, Tend + max(hm, (W + 2) / 3));
    }
}

static BuildRes sweepFinish(SweepRun& R, const Geo& g, const Pat& pt, const Prm& prm, int bestLen) {
    if (R.done) return R.res;
    BuildRes res;
    auto& lev = R.lev;
    int Tend = R.Tend;
    if (lev.size() <= (size_t)Tend) return res;
    vector<int> order(lev.back().size());
    iota(order.begin(), order.end(), 0);
    sort(order.begin(), order.end(), [&](int a, int b) { return lev.back()[a].score < lev.back()[b].score; });
    if ((int)order.size() > prm.M) order.resize(prm.M);
    vector<SNode> fin; for (int i : order) fin.push_back(lev.back()[i]);
    if (VERB >= 3) {
        const SNode& b = lev.back()[order[0]];
        fprintf(stderr, "   best sweep end (A=%d score %d) missP %d dirtyK %d cfg:", R.plan.A, b.score, bpop(andnot(pt.P, b.board)), bpop(pt.K & b.board));
        for (int e = 0; e < 6; e++) fprintf(stderr, " %d%s", SG[e] * b.cfg.k[e], b.cfg.on(e) ? "*" : "");
        fprintf(stderr, "\n");
        dumpBoard(g, pt, b.board);
    }
    Cleaner cl(g, pt, prm);
    vector<Ins> path; int ii; Cfg fc; Bits fb;
    int maxD = min(prm.Dc, bestLen - 1 - Tend);
    if (maxD < 1) return res;
    if (!cl.run(fin, maxD, path, ii, fc, fb)) return res;
    Prog p;
    int pi = order[ii];
    if (VERB >= 2) { fprintf(stderr, "   sweep end state (score %d):\n", lev.back()[pi].score); dumpBoard(g, pt, lev.back()[pi].board); }
    for (int tt = Tend; tt >= 1; tt--) { p.push_back(lev[tt][pi].ins); pi = lev[tt][pi].parent; }
    reverse(p.begin(), p.end());
    for (auto& in : path) p.push_back(in);
    res.ok = true; res.prog = p; res.cfg = fc; res.board = fb;
    return res;
}

static vector<Plan> makePlans(const Geo& g, const Pat& pt, const Prm& prm) {
    vector<Plan> plans;
    for (int A = 0; A < 3; A++) {
        int lo = 99, hi = -99;
        bforeach(pt.P, [&](int x) { lo = min(lo, g.L[A][x]); hi = max(hi, g.L[A][x]); });
        auto npRow = [&](int c) { return bpop(pt.P & g.line[A][c + g.r]); };
        for (int dir = -1; dir <= 1; dir += 2)
            for (int ts = 0; ts <= prm.trims; ts++)
                for (int te = 0; te <= prm.trims; te++) {
                    int s = dir > 0 ? lo + ts : hi - ts;
                    int e = dir > 0 ? hi - te : lo + te;
                    if ((e - s) * dir < 0) continue;
                    int maxPre = ts == 0 ? 0 : min(2, (npRow(s - dir) + 3) / 4);
                    int maxPost = te == 0 ? 0 : min(2, (npRow(e + dir) + 3) / 4);
                    if (ts >= 2 && maxPre > 0 && npRow(s - 2 * dir) > 0) {}  // outer row left to cleanup
                    for (int pre = 0; pre <= maxPre; pre++)
                        for (int post = 0; post <= maxPost; post++) {
                            Plan p; p.A = A; p.dir = dir; p.pre = pre; p.post = post;
                            for (int c = s; c != e + dir; c += dir) p.rows.push_back(c);
                            plans.push_back(p);
                        }
                }
    }
    return plans;
}

static int CLEANONLY = 0;
static BuildRes bestBuild(const Geo& g, const Pat& pt, const Prm& prm, const Cfg& cfg, const Bits& board, int bestLen, string& log) {
    BuildRes best;
    if (CLEANONLY) {
        Cleaner cl(g, pt, prm);
        SNode s; s.cfg = cfg; s.board = board; s.parent = -1; s.score = 0;
        vector<Ins> path; int ii; Cfg fc; Bits fb;
        if (cl.run({s}, min(prm.Dc, bestLen - 1), path, ii, fc, fb)) { best.ok = true; best.prog = path; best.cfg = fc; best.board = fb; }
        return best;
    }
    auto plans = makePlans(g, pt, prm);
    vector<unique_ptr<SweepRun>> runs;
    for (auto& pl : plans) {
        runs.emplace_back(new SweepRun());
        runs.back()->plan = pl;
        sweepRun(*runs.back(), g, pt, cfg, board, prm, bestLen);
        if (runs.back()->done && (int)runs.back()->res.prog.size() < bestLen) { best = runs.back()->res; bestLen = best.prog.size(); }
    }
    vector<int> ord(runs.size()); iota(ord.begin(), ord.end(), 0);
    sort(ord.begin(), ord.end(), [&](int a, int b) { return runs[a]->est < runs[b]->est; });
    int tried = 0, bi = -1;
    for (int i : ord) {
        if (runs[i]->done || runs[i]->est >= bestLen) continue;
        if (tried >= prm.maxCleanups) break;
        tried++;
        auto t0 = chrono::steady_clock::now();
        BuildRes r = sweepFinish(*runs[i], g, pt, prm, bestLen);
        if (VERB) fprintf(stderr, "  A=%d rows %d..%d Tend %d est %d -> %s len %zu (%.2fs)\n", runs[i]->plan.A, runs[i]->plan.rows.front(), runs[i]->plan.rows.back(), runs[i]->Tend, runs[i]->est, r.ok ? "ok" : "fail", r.prog.size(), chrono::duration<double>(chrono::steady_clock::now() - t0).count());
        if (r.ok && (int)r.prog.size() < bestLen) { best = r; bestLen = r.prog.size(); bi = i; }
        runs[i]->lev.clear(); runs[i]->lev.shrink_to_fit();
    }
    if (best.ok) {
        if (bi >= 0) log += " [A" + to_string(runs[bi]->plan.A) + " " + to_string(runs[bi]->plan.rows.front()) + ".." + to_string(runs[bi]->plan.rows.back()) + " T" + to_string(runs[bi]->Tend) + "+" + to_string((int)best.prog.size() - runs[bi]->Tend) + "]";
        else log += " [swept " + to_string(best.prog.size()) + "]";
    }
    return best;
}

static Prog solvePattern(const Geo& g, const Pat& pt, const Prm& prm, int limit, string& log) {
    Prog total;
    Cfg cfg = startCfg();
    Bits board = bz();
    for (int b = 0; b < pt.cnt; b++) {
        int bestLen = limit - (int)total.size() + 1 - (pt.cnt - 1 - b);  // leave >=1 tick per remaining build
        if (bestLen <= 1) return Prog();
        BuildRes r = bestBuild(g, pt, prm, cfg, board, bestLen, log);
        if (!r.ok) return Prog();
        for (auto& in : r.prog) total.push_back(in);
        cfg = r.cfg; board = r.board;
    }
    return total;
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
    if (argc < 4) { fprintf(stderr, "usage: sweep2 <in> <best_out> <res_dir> [opts]\n"); return 1; }
    string inp = argv[1], bestp = argv[2], resdir = argv[3];
    Prm prm; vector<int> only; int threads = 4;
    for (int i = 4; i + 1 < argc; i += 2) {
        string o = argv[i], v = argv[i + 1];
        if (o == "-Bs") prm.Bs = stoi(v);
        else if (o == "-Bc") prm.Bc = stoi(v);
        else if (o == "-M") prm.M = stoi(v);
        else if (o == "-Dc") prm.Dc = stoi(v);
        else if (o == "-mc") prm.maxCleanups = stoi(v);
        else if (o == "-tr") prm.trims = stoi(v);
        else if (o == "-wP") prm.wP = stoi(v);
        else if (o == "-wK") prm.wK = stoi(v);
        else if (o == "-kLB") prm.keyLB = stoi(v);
        else if (o == "-kW") prm.keyW = stoi(v);
        else if (o == "-la") { auto l = parseList(v); prm.la1 = l[0]; prm.la2 = l[1]; prm.la3 = l[2]; }
        else if (o == "-only") only = parseList(v);
        else if (o == "-v") VERB = stoi(v);
        else if (o == "-t") threads = stoi(v);
        else if (o == "-cleanonly") CLEANONLY = stoi(v);
    }
    auto pats = parseInput(inp);
    auto base = parseOutput(bestp);
    map<int, Geo> geos;
    for (auto& p : pats) if (!geos.count(p.r)) geos[p.r] = makeGeo(p.r);
    for (auto& p : pats) initPat(p, geos[p.r]);
    int n = pats.size();
    vector<int> todo;
    if (only.empty()) for (int i = 0; i < n; i++) todo.push_back(i); else todo = only;
    // biggest first
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
        Prog p = solvePattern(g, pats[i], prm, limit, log);
        bool win = !p.empty() && (int)p.size() <= limit && simulate(g, pats[i], p) >= pats[i].cnt;
        if (win) {
            vector<Prog> v = {p};
            writeOutput(fn, v);
        }
        #pragma omp critical
        {
            if (win) { wins++; gain += base[i].size() - p.size(); }
            fprintf(stderr, "pat %d r=%d cnt=%d P=%d K=%d best=%zu mine=%zu %s%s (%.1fs)\n", i, pats[i].r, pats[i].cnt, pats[i].np, pats[i].nk, base[i].size(), p.size(), win ? "WIN" : "", log.c_str(), chrono::duration<double>(chrono::steady_clock::now() - t0).count());
            fflush(stderr);
        }
    }
    printf("wins %d gain %d\n", wins, gain);
    return 0;
}
