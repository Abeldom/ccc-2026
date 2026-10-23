// shots3: LNS window re-planning with a generalized backward beam, starting from existing programs.
// usage: shots3 <input> <start programs> <out> [-t sec] [-W0 w] [-Wmax w] [-win n] [-only i,j-k] [-threads n]
//                [-root n] [-seed s] [-cw closW] [-free 0/1] [-fixed 0/1] [-v]
#include "core3.h"

struct LParams { double tsec = 20; int W0 = 100, Wmax = 1600, maxWin = 8, rootN = 300; u64 seed = 1; double closW = 0.0;
    int doFree = 1, doFixed = 1, verbose = 0, maxAct = 1, sumMode = 0, cycle = 1, randH = 1, libdp = 0, libW = 400, libRounds = 2, libPer = 30; BParams hb; };

static atomic<long long> g_improve{0};

template<int NW>
struct LNS {
    typedef array<u64, NW> BS;
    const Pat& P; const Geo& G; Ctx<NW> C; LParams lp;
    chrono::steady_clock::time_point t0;
    LNS(const Pat& p, const Geo& g, const LParams& l) : P(p), G(g), lp(l) { C.init(p, g); }
    double el() const { return chrono::duration<double>(chrono::steady_clock::now() - t0).count(); }

    bool toProg(const vector<Cfg>& cf, vector<string>& prog) const {
        prog.clear();
        for (size_t t = 1; t < cf.size(); t++) { string s = C.stepInstr(cf[t - 1], cf[t]); if (s.empty()) return false; prog.push_back(s); }
        return true;
    }
    // verify & truncate; returns new length or -1
    int verify(vector<string>& prog) const {
        int t = simulateExact(G, P, prog, P.cnt);
        if (t < 0) return -1;
        prog.resize(t); return t;
    }
    bool appendTransit(vector<Cfg>& out, const Cfg& A, const Cfg& B) const {
        int liveE; int n = C.closure(A, B, 0, &liveE);
        vector<Cfg> tr; if (n > 0 && !C.transit(A, B, n, liveE, tr)) return false;
        for (auto& c : tr) out.push_back(c);
        return true;
    }
    // first tick > T that creates something on P∪K (start of the next build's work), or -1
    int firstWork(const vector<Cfg>& cf, int T) const {
        for (int t = T + 1; t < (int)cf.size(); t++) { BS S, D; C.effect(cf[t], S, D); for (int i = 0; i < NW; i++) if (S[i] & C.All[i]) return t; }
        return -1;
    }

    // randomized heuristic portfolio
    void randH(BParams& bp, u64 sd) const {
        if (!lp.randH) return;
        u64 h = mix64(sd * 0x9E37 + 12345);
        static const double AL[4] = {0, 0.5, 1.0, 0.5}; static const double GA[2] = {0.5, 0.7};
        bp.alpha = AL[h & 3]; bp.gamma = GA[(h >> 2) & 1];
    }
    // try improving prog; returns improved program (or the same)
    vector<string> run(vector<string> prog) {
        t0 = chrono::steady_clock::now();
        int curLen = verify(prog);
        if (curLen < 0) { fprintf(stderr, "start program invalid\n"); return prog; }
        int W = lp.W0;
        u64 seed = lp.seed;
        // cycle construction attempt (fresh) for count>1 / any count
        if (lp.cycle) tryCycle(prog, curLen, max(lp.W0, 400), seed);
        if (lp.libdp) libDP(prog, curLen, lp.libW, lp.libRounds, lp.libPer, seed + 99991);
        while (el() < lp.tsec) {
            bool imp = false;
            vector<Cfg> cf = C.parseCfgs(prog);
            vector<BS> F; vector<int> acc; C.forward(cf, F, acc);
            if ((int)acc.size() < P.cnt) { fprintf(stderr, "fwd mismatch\n"); break; }
            int nb = P.cnt;
            // builds in order of length (longest first)
            vector<int> order(nb); iota(order.begin(), order.end(), 0);
            auto bstart = [&](int b) { return b == 0 ? 0 : acc[b - 1]; };
            sort(order.begin(), order.end(), [&](int x, int y) { return acc[x] - bstart(x) > acc[y] - bstart(y); });
            for (int b : order) {
                if (imp || el() >= lp.tsec) break;
                int s = bstart(b), T = acc[b];
                // suffix-undecided sets within this build
                vector<BS> Ud(T + 1);
                Ud[T] = C.All;
                for (int t = T; t > s; t--) { BS S, D; C.effect(cf[t], S, D); for (int i = 0; i < NW; i++) Ud[t - 1][i] = Ud[t][i] & ~(S[i] | D[i]); }
                int w1 = (b + 1 < nb) ? firstWork(cf, T) : -1;
                if (b + 1 < nb && w1 < 0) continue;
                // ---- free-root suffix windows
                if (lp.doFree) {
                    for (int t1 = s; t1 < T && !imp && el() < lp.tsec; t1++) {
                        Spec<NW> sp; sp.X = cf[t1]; sp.F = F[t1]; sp.U0 = C.All; sp.freeRoot = true;
                        int oldCost = T - t1;
                        if (b + 1 < nb) { sp.nextMode = 1; sp.W1 = cf[w1]; oldCost += w1 - 1 - T; }
                        sp.bound = oldCost;
                        BParams bp = lp.hb; bp.W = W; bp.finW = max(20, W / 4); bp.rootN = lp.rootN; bp.seed = seed++; bp.closW = lp.closW; bp.maxAct = lp.maxAct; bp.sumMode = lp.sumMode; randH(bp, seed);
                        Beam<NW> bm(C); double tb = el(); auto res = bm.run(sp, bp);
                        if (lp.verbose > 1) fprintf(stderr, "    free b%d t1=%d bound=%d -> %d (%.2fs)\n", b, t1, oldCost, res.total, el() - tb);
                        if (res.total >= oldCost) continue;
                        vector<Cfg> ncf(cf.begin(), cf.begin() + t1 + 1);
                        if (!appendTransit(ncf, cf[t1], res.seq[0])) continue;
                        for (auto& c : res.seq) ncf.push_back(c);
                        if (b + 1 < nb) { if (!appendTransit(ncf, res.seq.back(), cf[w1])) continue; for (int t = w1; t < (int)cf.size(); t++) ncf.push_back(cf[t]); }
                        vector<string> np2; if (!toProg(ncf, np2)) { fprintf(stderr, "adjacency bug (free)\n"); continue; }
                        int L = verify(np2);
                        if (L > 0 && L < curLen) {
                            if (lp.verbose) fprintf(stderr, "  [free b%d t1=%d] %d -> %d (W=%d)\n", b, t1, curLen, L, W);
                            g_improve += curLen - L; curLen = L; prog = np2; imp = true;
                        }
                    }
                }
                // ---- fixed windows (t1, t2) inside the build
                if (lp.doFixed) {
                    for (int t2 = T - 1; t2 > s && !imp && el() < lp.tsec; t2--) {
                        for (int t1 = t2 - 1; t1 >= s && t1 >= t2 - lp.maxWin && !imp; t1--) {
                            int oldCost = t2 - t1;
                            // deletion check
                            {
                                bool adj = !C.stepInstr(cf[t1], cf[t2 + 1]).empty();
                                bool cons = true;
                                for (int i = 0; i < NW; i++) if ((Ud[t2][i] & C.Pm[i] & ~F[t1][i]) | (Ud[t2][i] & C.Km[i] & F[t1][i])) cons = false;
                                if (adj && cons) {
                                    vector<Cfg> ncf(cf.begin(), cf.begin() + t1 + 1);
                                    for (int t = t2 + 1; t < (int)cf.size(); t++) ncf.push_back(cf[t]);
                                    vector<string> np2; if (toProg(ncf, np2)) { int L = verify(np2);
                                        if (L > 0 && L < curLen) { if (lp.verbose) fprintf(stderr, "  [del b%d %d..%d] %d -> %d\n", b, t1, t2, curLen, L);
                                            g_improve += curLen - L; curLen = L; prog = np2; imp = true; break; } }
                                }
                            }
                            if (oldCost < 2) continue;
                            Spec<NW> sp; sp.X = cf[t1]; sp.F = F[t1]; sp.U0 = Ud[t2]; sp.freeRoot = false; sp.Cnext = cf[t2 + 1]; sp.bound = oldCost;
                            BParams bp = lp.hb; bp.W = W; bp.finW = 10; bp.rootN = 4096; bp.seed = seed++; bp.closW = lp.closW; bp.maxAct = lp.maxAct; bp.maxDepth = oldCost + 1; bp.sumMode = lp.sumMode; randH(bp, seed);
                            Beam<NW> bm(C); auto res = bm.run(sp, bp);
                            if (res.total >= oldCost) continue;
                            vector<Cfg> ncf(cf.begin(), cf.begin() + t1 + 1);
                            if (!appendTransit(ncf, cf[t1], res.seq[0])) continue;
                            for (auto& c : res.seq) ncf.push_back(c);
                            for (int t = t2 + 1; t < (int)cf.size(); t++) ncf.push_back(cf[t]);
                            vector<string> np2; if (!toProg(ncf, np2)) { fprintf(stderr, "adjacency bug (fixed)\n"); continue; }
                            int L = verify(np2);
                            if (L > 0 && L < curLen) {
                                if (lp.verbose) fprintf(stderr, "  [fix b%d %d..%d] %d -> %d (W=%d)\n", b, t1, t2, curLen, L, W);
                                g_improve += curLen - L; curLen = L; prog = np2; imp = true;
                            }
                        }
                    }
                }
            }
            if (!imp) { if (W >= lp.Wmax) break; W *= 2; }
        }
        return prog;
    }

    // library of builds + DP composition for count>1 (also count 1: best single build from I)
    void libDP(vector<string>& prog, int& curLen, int W, int rounds, int perRound, u64 seed) {
        typedef typename Beam<NW>::Res Res;
        vector<Res> lib;
        vector<Cfg> starts; starts.push_back(C.I);
        set<vector<int>> seenStart;
        auto key = [&](const Cfg& c) { vector<int> v; for (int e = 0; e < 6; e++) v.push_back(c.k[e]); v.push_back(c.a); return v; };
        seenStart.insert(key(C.I));
        for (int rd = 0; rd <= rounds && el() < lp.tsec; rd++) {
            vector<Cfg> cur = starts; starts.clear();
            for (auto& X : cur) {
                if (el() >= lp.tsec) break;
                Spec<NW> sp; sp.X = X; sp.F.fill(0); sp.U0 = C.All; sp.freeRoot = true; sp.libMode = true; sp.bound = 1000;
                BParams bp = lp.hb; bp.W = W; bp.finW = max(20, W / 4); bp.rootN = lp.rootN; bp.seed = seed++; bp.maxAct = lp.maxAct; randH(bp, seed);
                Beam<NW> bm(C); size_t before = lib.size();
                bm.run(sp, bp, &lib, perRound, 2);
                for (size_t i = before; i < lib.size(); i++) { Cfg y = lib[i].Y; if (seenStart.insert(key(y)).second) starts.push_back(y); }
            }
            // limit next-round starts: prefer Ys of the shortest builds
            if ((int)starts.size() > 6) starts.resize(6);
        }
        if (lib.empty()) return;
        int B = lib.size(), c = P.cnt;
        // DP over build sequences
        const int INF = 1 << 29;
        vector<vector<int>> f(c, vector<int>(B, INF)), par(c, vector<int>(B, -1));
        for (int e = 0; e < B; e++) f[0][e] = C.closure(C.I, lib[e].seq[0], lib[e].freeM) + (int)lib[e].seq.size();
        for (int i = 1; i < c; i++) for (int e = 0; e < B; e++) {
            int best = INF, bp = -1;
            for (int q = 0; q < B; q++) if (f[i - 1][q] < INF) { int v = f[i - 1][q] + C.closure(lib[q].Y, lib[e].seq[0], lib[e].freeM); if (v < best) { best = v; bp = q; } }
            f[i][e] = best + (int)lib[e].seq.size(); par[i][e] = bp; }
        int be = min_element(f[c - 1].begin(), f[c - 1].end()) - f[c - 1].begin();
        if (lp.verbose) fprintf(stderr, "  [libdp] lib=%d dp=%d cur=%d\n", B, f[c - 1][be], curLen);
        if (f[c - 1][be] >= curLen) return;
        vector<int> seqIdx(c); seqIdx[c - 1] = be; for (int i = c - 1; i > 0; i--) seqIdx[i - 1] = par[i][seqIdx[i]];
        vector<Cfg> ncf; ncf.push_back(C.I); Cfg cur = C.I;
        for (int i = 0; i < c; i++) {
            vector<Cfg> sq = lib[seqIdx[i]].seq; uint8_t fm = lib[seqIdx[i]].freeM;
            for (auto& x : sq) for (int e = 0; e < 6; e++) if ((fm >> e) & 1) { x.k[e] = cur.k[e]; x.a &= ~(1 << e); }
            if (!appendTransit(ncf, cur, sq[0])) return;
            for (auto& x : sq) ncf.push_back(x);
            cur = sq.back();
        }
        vector<string> np2; if (!toProg(ncf, np2)) { fprintf(stderr, "adjacency bug (libdp)\n"); return; }
        int L = verify(np2);
        if (lp.verbose) fprintf(stderr, "  [libdp] built %d (pred %d)\n", L, f[c - 1][be]);
        if (L > 0 && L < curLen) { g_improve += curLen - L; curLen = L; prog = np2; }
    }

    // fresh cycle construction (same build repeated); replaces prog if shorter
    void tryCycle(vector<string>& prog, int& curLen, int W, u64 seed) {
        Spec<NW> sp; sp.X = C.I; sp.F.fill(0); sp.U0 = C.All; sp.freeRoot = true; sp.nextMode = 2; sp.cnt = P.cnt; sp.bound = curLen;
        BParams bp = lp.hb; bp.W = W; bp.finW = max(20, W / 4); bp.rootN = lp.rootN; bp.seed = seed; bp.closW = lp.closW; bp.maxAct = lp.maxAct; bp.sumMode = lp.sumMode;
        Beam<NW> bm(C); auto res = bm.run(sp, bp);
        if (res.total >= curLen) return;
        vector<Cfg> ncf; ncf.push_back(C.I);
        Cfg cur = C.I;
        for (int b = 0; b < P.cnt; b++) {
            if (!appendTransit(ncf, cur, res.seq[0])) return;
            for (auto& c : res.seq) ncf.push_back(c);
            cur = res.seq.back();
        }
        vector<string> np2; if (!toProg(ncf, np2)) { fprintf(stderr, "adjacency bug (cycle)\n"); return; }
        int L = verify(np2);
        if (L > 0 && L < curLen) { if (lp.verbose) fprintf(stderr, "  [cycle] %d -> %d\n", curLen, L); g_improve += curLen - L; curLen = L; prog = np2; }
    }
};

template<int NW>
vector<string> runLNS(const Pat& P, const Geo& G, const vector<string>& prog, const LParams& lp) {
    LNS<NW> l(P, G, lp); return l.run(prog);
}

template<int NW>
int benchOne(const Pat& P, const Geo& G, const LParams& lp) {
    Ctx<NW> C; C.init(P, G);
    Spec<NW> sp; sp.X = C.I; sp.F.fill(0); sp.U0 = C.All; sp.freeRoot = true; sp.nextMode = 0; sp.bound = 100000;
    BParams bp = lp.hb; bp.W = lp.W0; bp.finW = max(20, lp.W0 / 4); bp.rootN = lp.rootN; bp.seed = lp.seed; bp.closW = lp.closW; bp.maxAct = lp.maxAct; bp.sumMode = lp.sumMode;
    Beam<NW> bm(C); auto res = bm.run(sp, bp);
    return res.total;
}

int main(int argc, char** argv) {
    if (argc < 4) { fprintf(stderr, "usage: shots3 in start out [opts]\n"); return 1; }
    string inp = argv[1], startf = argv[2], outf = argv[3];
    LParams lp; set<int> only; int threads = 4; double topFrac = 1.0; int bench = 0;
    for (int i = 4; i < argc; i++) {
        string o = argv[i]; auto nx = [&]() { return string(argv[++i]); };
        if (o == "-t") lp.tsec = stod(nx()); else if (o == "-W0") lp.W0 = stoi(nx()); else if (o == "-Wmax") lp.Wmax = stoi(nx());
        else if (o == "-win") lp.maxWin = stoi(nx()); else if (o == "-threads") threads = stoi(nx()); else if (o == "-root") lp.rootN = stoi(nx());
        else if (o == "-seed") lp.seed = stoull(nx()); else if (o == "-cw") lp.closW = stod(nx()); else if (o == "-free") lp.doFree = stoi(nx());
        else if (o == "-fixed") lp.doFixed = stoi(nx()); else if (o == "-v") lp.verbose = 1; else if (o == "-vv") lp.verbose = 2; else if (o == "-act") lp.maxAct = stoi(nx());
        else if (o == "-sum") lp.sumMode = stoi(nx()); else if (o == "-al") lp.hb.alpha = stod(nx()); else if (o == "-be") lp.hb.beta = stod(nx()); else if (o == "-ga") lp.hb.gamma = stod(nx()); else if (o == "-mu") lp.hb.mult = stod(nx()); else if (o == "-clw") lp.hb.cloW = stod(nx()); else if (o == "-bench") bench = 1; else if (o == "-rh") lp.randH = stoi(nx()); else if (o == "-lib") lp.libdp = stoi(nx()); else if (o == "-libW") lp.libW = stoi(nx()); else if (o == "-libR") lp.libRounds = stoi(nx()); else if (o == "-libP") lp.libPer = stoi(nx()); else if (o == "-cycle") lp.cycle = stoi(nx()); else if (o == "-top") topFrac = stod(nx());
        else if (o == "-only") { string s = nx(); stringstream ss(s); string t;
            while (getline(ss, t, ',')) { size_t d = t.find('-'); if (d != string::npos && d > 0) { int a = stoi(t.substr(0, d)), b = stoi(t.substr(d + 1)); for (int q = a; q <= b; q++) only.insert(q); } else only.insert(stoi(t)); } }
        else { fprintf(stderr, "unknown option %s\n", o.c_str()); return 1; }
    }
    auto pats = parseInput(inp);
    auto progs = readProgs(startf);
    int N = pats.size();
    if ((int)progs.size() != N) { fprintf(stderr, "program count mismatch %d vs %d\n", (int)progs.size(), N); return 1; }
    map<int, Geo> geos; for (auto& p : pats) if (!geos.count(p.r)) geos[p.r] = makeGeo(p.r);
    vector<int> jobs;
    for (int i = 0; i < N; i++) if (only.empty() || only.count(i)) jobs.push_back(i);
    sort(jobs.begin(), jobs.end(), [&](int a, int b) { return progs[a].size() > progs[b].size(); });
    if (topFrac < 1.0) jobs.resize(max(1, (int)(jobs.size() * topFrac)));
    if (bench) {
#ifdef _OPENMP
        omp_set_num_threads(threads);
#endif
        vector<int> tot(N, 0); auto tb = chrono::steady_clock::now();
        #pragma omp parallel for schedule(dynamic, 1)
        for (int j = 0; j < (int)jobs.size(); j++) {
            int i = jobs[j]; const Pat& P = pats[i]; const Geo& G = geos.at(P.r);
            int ncell = 3 * P.r * (P.r - 1) + 1; int nw = (ncell + 63) / 64; int v;
            switch (nw) { case 1: v = benchOne<1>(P, G, lp); break; case 2: v = benchOne<2>(P, G, lp); break; case 3: v = benchOne<3>(P, G, lp); break;
                case 4: v = benchOne<4>(P, G, lp); break; case 5: v = benchOne<5>(P, G, lp); break; case 6: v = benchOne<6>(P, G, lp); break;
                case 7: v = benchOne<7>(P, G, lp); break; default: v = benchOne<8>(P, G, lp); break; }
            tot[i] = v;
        }
        long long s = 0; for (int i : jobs) { s += tot[i]; printf("%d:%d ", i, tot[i]); }
        printf("\nBENCH total %lld time %.1f\n", s, chrono::duration<double>(chrono::steady_clock::now() - tb).count());
        return 0;
    }
    vector<vector<string>> out = progs;
    long long before = 0; for (auto& p : progs) before += p.size();
    auto tw = chrono::steady_clock::now();
#ifdef _OPENMP
    omp_set_num_threads(threads);
#endif
    #pragma omp parallel for schedule(dynamic, 1)
    for (int j = 0; j < (int)jobs.size(); j++) {
        int i = jobs[j];
        const Pat& P = pats[i]; const Geo& G = geos.at(P.r);
        int ncell = 3 * P.r * (P.r - 1) + 1; int nw = (ncell + 63) / 64;
        auto t0 = chrono::steady_clock::now();
        vector<string> res;
        switch (nw) {
            case 1: res = runLNS<1>(P, G, progs[i], lp); break;
            case 2: res = runLNS<2>(P, G, progs[i], lp); break;
            case 3: res = runLNS<3>(P, G, progs[i], lp); break;
            case 4: res = runLNS<4>(P, G, progs[i], lp); break;
            case 5: res = runLNS<5>(P, G, progs[i], lp); break;
            case 6: res = runLNS<6>(P, G, progs[i], lp); break;
            case 7: res = runLNS<7>(P, G, progs[i], lp); break;
            default: res = runLNS<8>(P, G, progs[i], lp); break;
        }
        double dt = chrono::duration<double>(chrono::steady_clock::now() - t0).count();
        #pragma omp critical
        {
            if (res.size() < out[i].size()) out[i] = res;
            fprintf(stderr, "case %d r=%d cnt=%d: %d -> %d (%.1fs) [total gain %lld]\n", i, P.r, P.cnt, (int)progs[i].size(), (int)out[i].size(), dt, (long long)g_improve);
            if (chrono::duration<double>(chrono::steady_clock::now() - tw).count() > 60) { writeProgs(outf, out); tw = chrono::steady_clock::now(); }
        }
    }
    writeProgs(outf, out);
    long long after = 0; for (auto& p : out) after += p.size();
    fprintf(stderr, "total %lld -> %lld\n", before, after);
}
