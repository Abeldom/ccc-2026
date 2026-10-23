// Sweep solver: duplicate-pair line sweep with beam-searched cutters + beam-searched cleanup.
#include "common.h"

struct Params {
    int Bs = 300;     // sweep beam width
    int Bc = 60;      // cleanup beam width
    int M = 40;       // nodes handed from sweep to cleanup
    int Dc = 30;      // max cleanup depth
    int wDirty = 2;   // weight of dirty K cell in sweep score
    int la1 = 3, la2 = 2, la3 = 1;
    int nPlansRange = 2;  // trims tried at each end
    int maxCleanups = 6;
};

static int VERB = 0, VERBB = 0;
struct BuildRes {
    bool ok = false;
    Prog prog;
    Cfg cfg;
    Bits board;
};

struct SNode {
    Cfg cfg;
    Bits board;
    int parent;
    int key;
    int score;
    Ins ins;
};

// ---------------- cleanup beam -----------------
// LB of ticks to fix wrong cells.
struct Cleaner {
    const Geo& g; const Pat& pt; const Params& prm;
    Cleaner(const Geo& g_, const Pat& p_, const Params& pr) : g(g_), pt(p_), prm(pr) {}

    // returns (h, sumLB) ; W computed outside
    inline void lb(const Cfg& c, const Bits& wrongP, const Bits& wrongK, int& hmax, int& hsum) const {
        int ce[6]; bool on[6];
        for (int e = 0; e < 6; e++) { ce[e] = SG[e] * c.k[e]; on[e] = c.on(e); }
        hmax = 0; hsum = 0;
        auto proc = [&](int x, int need) {
            int a[6];
            for (int e = 0; e < 6; e++) {
                int d = abs(ce[e] - g.L[AX[e]][x]) + (on[e] ? 0 : 1);
                a[e] = d < 1 ? 1 : d;
            }
            // need-th smallest (1-based)
            sort(a, a + 6);
            int v = a[need - 1];
            hmax = max(hmax, v); hsum += v;
        };
        bforeach(wrongP, [&](int x) { proc(x, 2); });
        bforeach(wrongK, [&](int x) { proc(x, 3); });
    }

    // run cleanup beam from initial nodes; returns true + path (instructions) + final cfg/board
    bool run(const vector<SNode>& init, int maxDepth, vector<Ins>& path, int& initIdx, Cfg& fcfg, Bits& fboard) const {
        const Bits Z = bz();
        vector<vector<SNode>> lev;
        lev.push_back(init);
        unordered_map<u64, int> seen;
        vector<SNode> cand;
        for (int d = 1; d <= maxDepth; d++) {
            cand.clear(); seen.clear();
            const vector<SNode>& cur = lev.back();
            for (int ni = 0; ni < (int)cur.size(); ni++) {
                const SNode& nd = cur[ni];
                // options per emitter
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
                    Cfg c;
                    c.act = 0;
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
                        // success: reconstruct
                        path.clear(); path.push_back(ins);
                        int pi = ni;
                        for (int dd = d - 1; dd >= 1; dd--) { path.push_back(lev[dd][pi].ins); pi = lev[dd][pi].parent; }
                        reverse(path.begin(), path.end());
                        initIdx = pi;
                        fcfg = c; fboard = andnot(nb, pt.P);
                        return true;
                    }
                    int hm, hs; lb(c, wP, wK, hm, hs);
                    int W = bpop(wP) + bpop(wK);
                    int key = (d + hm) * 1000000 + W * 1000 + min(hs, 999);
                    u64 ck = c.key();
                    auto it = seen.find(ck);
                    if (it != seen.end()) {
                        if (cand[it->second].key <= key) continue;
                        SNode& s = cand[it->second];
                        s.cfg = c; s.board = nb; s.parent = ni; s.key = key; s.ins = ins;
                    } else {
                        seen[ck] = cand.size();
                        SNode s; s.cfg = c; s.board = nb; s.parent = ni; s.key = key; s.score = W; s.ins = ins;
                        cand.push_back(s);
                    }
                }
            }
            if (cand.empty()) return false;
            if ((int)cand.size() > prm.Bc) {
                nth_element(cand.begin(), cand.begin() + prm.Bc, cand.end(), [](const SNode& a, const SNode& b) { return a.key < b.key; });
                cand.resize(prm.Bc);
            }
            lev.push_back(cand);
        }
        return false;
    }
};

// ---------------- sweep build -----------------
struct SweepRun {
    bool done = false;      // accepted during sweep
    BuildRes res;
    int Tend = 0;
    vector<vector<SNode>> lev;
    int est = 1 << 30;
};
static void dumpBoard(const Geo& g, const Pat& pt, const Bits& b) {
    vector<string> grid = pt.grid;
    for (int x = 0; x < g.nc; x++) {
        char& ch = grid[g.cy[x]][g.cx[x]];
        bool P = bget(pt.P, x), K = bget(pt.K, x), m = bget(b, x);
        ch = P ? (m ? 'M' : 'm') : K ? (m ? 'k' : 'o') : (m ? '#' : '.');
    }
    for (auto& r : grid) fprintf(stderr, "    %s\n", r.c_str());
}
static void sweepRun(SweepRun& R, const Geo& g, const Pat& pt, const Cfg& c0, const Bits& b0, int A, const vector<int>& rows, const Params& prm, int bestSoFar) {
    BuildRes& res = R.res;
    const Bits Z = bz();
    int em[2] = {A, A + 3};
    int nrows = rows.size();
    int kt[2], d[2], need[2], T1 = 1;
    for (int i = 0; i < 2; i++) {
        kt[i] = SG[em[i]] * rows[0];
        d[i] = abs(kt[i] - c0.k[em[i]]);
        need[i] = !c0.on(em[i]);
        T1 = max(T1, d[i] + need[i]);
    }
    int Tend = T1 + nrows - 1;
    R.Tend = Tend;
    if (Tend >= bestSoFar) return;
    vector<array<char, 2>> dup(Tend + 1);
    for (int t = 1; t <= T1; t++)
        for (int i = 0; i < 2; i++) {
            char ch = '-';
            if (t <= d[i]) ch = (kt[i] > c0.k[em[i]]) ? 'F' : 'B';
            else if (t == T1 && need[i]) ch = 'T';
            dup[t][i] = ch;
        }
    for (int j = 1; j < nrows; j++) {
        int dr = rows[j] - rows[j - 1];
        for (int i = 0; i < 2; i++) dup[T1 + j][i] = (SG[em[i]] * dr > 0) ? 'F' : 'B';
    }
    int lo = *min_element(rows.begin(), rows.end()), hi = *max_element(rows.begin(), rows.end());
    vector<Bits> dec(Tend + 1);
    Bits base = bz();
    for (int x = 0; x < g.nc; x++) if (g.L[A][x] < lo || g.L[A][x] > hi) bset(base, x);
    for (int t = 0; t <= Tend; t++) {
        dec[t] = base;
        for (int j = 0; j < nrows && T1 + j <= t; j++) dec[t] = dec[t] | g.line[A][rows[j] + g.r];
    }
    // K cells per sweep index
    vector<vector<int>> krow(nrows);
    for (int j = 0; j < nrows; j++) bforeach(pt.K & g.line[A][rows[j] + g.r], [&](int x) { krow[j].push_back(x); });
    int cut[4], nc = 0;
    for (int e = 0; e < 6; e++) if (AX[e] != A) cut[nc++] = e;

    vector<vector<SNode>>& lev = R.lev;
    SNode s0; s0.cfg = c0; s0.board = b0; s0.parent = -1; s0.key = 0; s0.score = 0;
    lev.push_back({s0});
    unordered_map<u64, int> seen;
    vector<SNode> cand;
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
                    // accepted at tick t
                    Prog p; p.push_back(ins);
                    int pi = ni;
                    for (int tt = t - 1; tt >= 1; tt--) { p.push_back(lev[tt][pi].ins); pi = lev[tt][pi].parent; }
                    reverse(p.begin(), p.end());
                    res.ok = true; res.prog = p; res.cfg = c; res.board = andnot(nb, pt.P);
                    R.done = true; R.est = t;
                    return;
                }
                int score = bpop(wP & dec[t]) + prm.wDirty * bpop(wK & dec[t]);
                // lookahead: unreachable K cells in next rows
                int la = 0;
                for (int dl = 1; dl <= 3; dl++) {
                    int j = t + dl - T1;
                    if (j < 0 || j >= nrows) continue;
                    int w = dl == 1 ? prm.la1 : dl == 2 ? prm.la2 : prm.la3;
                    if (!w) continue;
                    for (int x : krow[j]) {
                        bool ok = false;
                        for (int q = 0; q < 4 && !ok; q++) {
                            int e = cut[q];
                            int dist = abs(SG[e] * c.k[e] - g.L[AX[e]][x]) + (c.on(e) ? 0 : 1);
                            if (dist <= dl) ok = true;
                        }
                        if (!ok) la += w;
                    }
                }
                int key = score * 8 + la;
                u64 ck = c.key();
                auto it = seen.find(ck);
                if (it != seen.end()) {
                    SNode& s = cand[it->second];
                    if (s.key <= key) continue;
                    s.cfg = c; s.board = nb; s.parent = ni; s.key = key; s.score = score; s.ins = ins;
                } else {
                    seen[ck] = cand.size();
                    SNode s; s.cfg = c; s.board = nb; s.parent = ni; s.key = key; s.score = score; s.ins = ins;
                    cand.push_back(s);
                }
            }
        }
        if ((int)cand.size() > prm.Bs) {
            nth_element(cand.begin(), cand.begin() + prm.Bs, cand.end(), [](const SNode& a, const SNode& b) { return a.key < b.key; });
            cand.resize(prm.Bs);
        }
        lev.push_back(cand);
    }
    // estimate: Tend + LB of best few final nodes
    Cleaner cl(g, pt, prm);
    for (auto& nd : lev.back()) {
        if (nd.score > lev.back()[0].score + 6) continue;
        int hm, hs; cl.lb(nd.cfg, andnot(pt.P, nd.board), pt.K & nd.board, hm, hs);
        int W = bpop(andnot(pt.P, nd.board)) + bpop(pt.K & nd.board);
        R.est = min(R.est, Tend + max(hm, (W + 1) / 2));
    }
}

static BuildRes sweepFinish(SweepRun& R, const Geo& g, const Pat& pt, const Params& prm, int bestSoFar) {
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
    Cleaner cl(g, pt, prm);
    vector<Ins> path; int ii; Cfg fc; Bits fb;
    int maxD = min(prm.Dc, bestSoFar - Tend - 1);
    if (maxD < 1) return res;
    if (!cl.run(fin, maxD, path, ii, fc, fb)) return res;
    Prog p;
    int pi = order[ii];
    if (VERBB) { fprintf(stderr, "   sweep end state (score %d):\n", lev.back()[pi].score); dumpBoard(g, pt, lev.back()[pi].board); }
    for (int tt = Tend; tt >= 1; tt--) { p.push_back(lev[tt][pi].ins); pi = lev[tt][pi].parent; }
    reverse(p.begin(), p.end());
    for (auto& in : path) p.push_back(in);
    res.ok = true; res.prog = p; res.cfg = fc; res.board = fb;
    return res;
}

static vector<vector<int>> makePlans(const Geo& g, const Pat& pt, int A, const Params& prm) {
    int lo = 99, hi = -99;
    bforeach(pt.P, [&](int x) { lo = min(lo, g.L[A][x]); hi = max(hi, g.L[A][x]); });
    vector<vector<int>> plans;
    for (int dir = -1; dir <= 1; dir += 2)
        for (int ts = 0; ts < prm.nPlansRange; ts++)
            for (int te = 0; te < prm.nPlansRange; te++) {
                int s = dir > 0 ? lo + ts : hi - ts;
                int e = dir > 0 ? hi - te : lo + te;
                if ((e - s) * dir < 0) continue;
                vector<int> rows;
                for (int c = s; c != e + dir; c += dir) rows.push_back(c);
                plans.push_back(rows);
            }
    return plans;
}

static Prog solvePattern(const Geo& g, const Pat& pt, const Params& prm, int limit, string& log) {
    Prog total;
    Cfg cfg = startCfg();
    Bits board = bz();
    for (int b = 0; b < pt.cnt; b++) {
        BuildRes best; int bestLen = limit - (int)total.size() + 1;
        if (bestLen <= 0) return Prog();
        vector<unique_ptr<SweepRun>> runs;
        vector<pair<int, vector<int>>> info;
        for (int A = 0; A < 3; A++) {
            auto plans = makePlans(g, pt, A, prm);
            for (auto& rows : plans) {
                runs.emplace_back(new SweepRun());
                sweepRun(*runs.back(), g, pt, cfg, board, A, rows, prm, bestLen);
                info.push_back({A, rows});
                if (runs.back()->done && (int)runs.back()->res.prog.size() < bestLen) { best = runs.back()->res; bestLen = best.prog.size(); }
            }
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
            if (VERB) fprintf(stderr, "  build %d A=%d rows %d..%d Tend %d est %d -> %s len %zu (%.2fs)\n", b, info[i].first, info[i].second.front(), info[i].second.back(), runs[i]->Tend, runs[i]->est, r.ok ? "ok" : "fail", r.prog.size(), chrono::duration<double>(chrono::steady_clock::now() - t0).count());
            if (r.ok && (int)r.prog.size() < bestLen) { best = r; bestLen = r.prog.size(); bi = i; }
        }
        if (!best.ok) return Prog();
        if (bi >= 0) log += " [A" + to_string(info[bi].first) + " " + to_string(info[bi].second.front()) + ".." + to_string(info[bi].second.back()) + " T" + to_string(runs[bi]->Tend) + "+" + to_string((int)best.prog.size() - runs[bi]->Tend) + "]";
        else log += " [swept " + to_string(best.prog.size()) + "]";
        for (auto& in : best.prog) total.push_back(in);
        cfg = best.cfg; board = best.board;
    }
    return total;
}

int main(int argc, char** argv) {
    if (argc < 2) { fprintf(stderr, "usage: solver check <in> <out> | solve <in> <baseline> <out> [Bs Bc M]\n"); return 1; }
    string mode = argv[1];
    if (mode == "check") {
        auto pats = parseInput(argv[2]);
        auto progs = parseOutput(argv[3]);
        if (progs.size() != pats.size()) { printf("program count mismatch %zu vs %zu\n", progs.size(), pats.size()); return 1; }
        map<int, Geo> geos;
        long tot = 0; int bad = 0;
        for (size_t i = 0; i < pats.size(); i++) {
            int r = pats[i].r;
            if (!geos.count(r)) geos[r] = makeGeo(r);
            initPat(pats[i], geos[r]);
            int a = simulate(geos[r], pats[i], progs[i]);
            tot += progs[i].size();
            if (a < pats[i].cnt) { bad++; printf("case %zu: accepted %d < %d\n", i, a, pats[i].cnt); }
        }
        printf("score %ld invalid %d\n", tot, bad);
        return 0;
    }
    if (mode == "solve") {
        auto pats = parseInput(argv[2]);
        auto base = parseOutput(argv[3]);
        string outp = argv[4];
        Params prm;
        if (argc > 5) prm.Bs = atoi(argv[5]);
        if (argc > 6) prm.Bc = atoi(argv[6]);
        if (argc > 7) prm.M = atoi(argv[7]);
        int only = -1; if (argc > 8) only = atoi(argv[8]);
        if (argc > 9) VERB = atoi(argv[9]);
        if (argc > 10) VERBB = atoi(argv[10]);
        map<int, Geo> geos;
        for (auto& p : pats) if (!geos.count(p.r)) geos[p.r] = makeGeo(p.r);
        for (auto& p : pats) initPat(p, geos[p.r]);
        int n = pats.size();
        vector<Prog> out(n);
        vector<int> src(n, 0);
        #pragma omp parallel for schedule(dynamic, 1)
        for (int i = 0; i < n; i++) {
            if (only >= 0 && i != only) { out[i] = base[i]; continue; }
            const Geo& g = geos[pats[i].r];
            string log;
            Prog p = solvePattern(g, pats[i], prm, base[i].size() - 1, log);
            if (!p.empty() && simulate(g, pats[i], p) >= pats[i].cnt && p.size() < base[i].size()) { out[i] = p; src[i] = 1; }
            else out[i] = base[i];
            #pragma omp critical
            { fprintf(stderr, "pat %d r=%d cnt=%d P=%d K=%d base=%zu sweep=%zu %s%s\n", i, pats[i].r, pats[i].cnt, pats[i].np, pats[i].nk, base[i].size(), p.size(), src[i] ? "WIN" : "", log.c_str()); }
        }
        long tot = 0; for (auto& p : out) tot += p.size();
        writeOutput(outp, out);
        printf("total %ld\n", tot);
        return 0;
    }
    return 1;
}
