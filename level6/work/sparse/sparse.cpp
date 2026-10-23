// sparse solver driver.
// Usage: sparse <input> <best_out_file> <idxlist|all> <maxP> <mode> [options]
//   mode first : min first-build length from start (iterative deepening)
//   mode dp    : layered DP over acceptances with target best-1
// Writes results to stdout lines "RES idx len" followed by the program, and stats to stderr.
#include "core.h"
#include <chrono>

static double now() { return chrono::duration<double>(chrono::steady_clock::now().time_since_epoch()).count(); }

struct Opt {
    long long nodeLimit = 50000000;
    size_t maxEnds = 20000;
    int beam = 4000;
    int extraLen = 1;     // explore builds up to (shortest found from that node) + extraLen
    double timeLimit = 600;
};

static Cfg swapCols(const Cfg& c, int mask) {
    Cfg d = c;
    for (int a = 0; a < 3; a++) if (mask >> a & 1) swap(d[a], d[a + 3]);
    return d;
}
static int canonMask(const Cfg& c) {
    int m = 0; for (int a = 0; a < 3; a++) if (c[a] > c[a + 3]) m |= 1 << a; return m;
}

struct Solver {
    const PInfo& P;
    Opt opt;
    int mbl = 2;
    // build cache: canonical key -> per len results (relative to canonical real config)
    struct CacheEnt { map<int, vector<BuildRes>> byLen; map<int, bool> aborted; };
    unordered_map<u64, CacheEnt> cache;
    long long totalNodes = 0;
    double t0;
    bool timeUp = false;
    Solver(const PInfo& p, Opt o) : P(p), opt(o) { t0 = now(); }

    const vector<BuildRes>& builds(const Cfg& canon, int len, bool& ab) {
        u64 k = cfgKey(canon);
        CacheEnt& ce = cache[k];
        auto it = ce.byLen.find(len);
        if (it != ce.byLen.end()) { ab = ce.aborted[len]; return it->second; }
        BuildSearch bs(P);
        bs.nodeLimit = opt.nodeLimit; bs.maxEnds = opt.maxEnds;
        bs.run(canon, len);
        totalNodes += bs.nodes;
        vector<BuildRes> v; v.reserve(bs.ends.size());
        for (auto& kv : bs.ends) v.push_back(std::move(kv.second));
        ce.aborted[len] = bs.aborted; ab = bs.aborted;
        return ce.byLen[len] = std::move(v);
    }

    // min first build length via iterative deepening (returns -1 if > maxLen)
    int firstBuild(int maxLen, vector<Cfg>* path) {
        Cfg s = startCfgP(P);
        for (int len = 1; len <= maxLen; len++) {
            BuildSearch bs(P); bs.nodeLimit = opt.nodeLimit; bs.firstOnly = true;
            bs.run(s, len); totalNodes += bs.nodes;
            if (!bs.ends.empty()) { if (path) *path = bs.ends.begin()->second.path; return len; }
            if (bs.aborted) return -2 - len;
        }
        return -1;
    }

    struct DNode { Cfg cfg; int cost; int par; vector<Cfg> path; };
    // returns timeline (start + configs) of best program with total <= target, or empty
    vector<Cfg> dp(int target) {
        int c = P.cnt;
        vector<vector<DNode>> layers(c + 1);
        layers[0].push_back({startCfgP(P), 0, -1, {}});
        for (int a = 0; a < c; a++) {
            auto& Lc = layers[a];
            sort(Lc.begin(), Lc.end(), [](const DNode& x, const DNode& y) { return x.cost < y.cost; });
            if ((int)Lc.size() > opt.beam) Lc.resize(opt.beam);
            // reindex parents not needed: parents reference previous layer which is final
            unordered_map<u64, int> nxt;
            auto& Ln = layers[a + 1];
            for (int i = 0; i < (int)Lc.size(); i++) {
                if (now() - t0 > opt.timeLimit) { timeUp = true; break; }
                const DNode& X = Lc[i];
                int maxLen = target - X.cost - (c - a - 1) * mbl;
                int msk = canonMask(X.cfg);
                Cfg canon = swapCols(X.cfg, msk);
                int firstFound = -1;
                for (int len = 1; len <= maxLen; len++) {
                    if (firstFound >= 0 && len > firstFound + opt.extraLen) break;
                    bool ab;
                    const vector<BuildRes>& v = builds(canon, len, ab);
                    if (!v.empty() && firstFound < 0) firstFound = len;
                    for (const BuildRes& br : v) {
                        Cfg end = swapCols(br.end, msk);
                        u64 k = cfgKey(end);
                        int cost = X.cost + len;
                        auto it = nxt.find(k);
                        if (it != nxt.end() && Ln[it->second].cost <= cost) continue;
                        DNode nd; nd.cfg = end; nd.cost = cost; nd.par = i;
                        nd.path.resize(br.path.size());
                        for (size_t q = 0; q < br.path.size(); q++) nd.path[q] = swapCols(br.path[q], msk);
                        if (it != nxt.end()) Ln[it->second] = std::move(nd);
                        else { nxt[k] = Ln.size(); Ln.push_back(std::move(nd)); }
                    }
                }
            }
            if (timeUp) return {};
            fprintf(stderr, "  layer %d: %zu nodes (min cost %d) cache %zu nodes %lld t=%.1f\n", a + 1, Ln.size(),
                    Ln.empty() ? -1 : min_element(Ln.begin(), Ln.end(), [](const DNode& x, const DNode& y) { return x.cost < y.cost; })->cost,
                    cache.size(), totalNodes, now() - t0);
            if (Ln.empty()) return {};
        }
        // best final
        auto& Lf = layers[c];
        int bi = 0; for (int i = 1; i < (int)Lf.size(); i++) if (Lf[i].cost < Lf[bi].cost) bi = i;
        vector<vector<Cfg>> segs;
        int a = c, i = bi;
        while (a > 0) { segs.push_back(layers[a][i].path); i = layers[a][i].par; a--; }
        vector<Cfg> tl; tl.push_back(startCfgP(P));
        for (int q = segs.size() - 1; q >= 0; q--) for (auto& x : segs[q]) tl.push_back(x);
        return tl;
    }
};

int main(int argc, char** argv) {
    if (argc < 6) { fprintf(stderr, "usage\n"); return 1; }
    string inp = argv[1], bestf = argv[2], idxs = argv[3];
    int maxP = atoi(argv[4]); string mode = argv[5];
    Opt opt;
    int threads = 2;
    for (int i = 6; i + 1 < argc; i += 2) {
        string k = argv[i]; double v = atof(argv[i + 1]);
        if (k == "-nodes") opt.nodeLimit = (long long)v;
        else if (k == "-ends") opt.maxEnds = (size_t)v;
        else if (k == "-beam") opt.beam = (int)v;
        else if (k == "-extra") opt.extraLen = (int)v;
        else if (k == "-time") opt.timeLimit = v;
        else if (k == "-threads") threads = (int)v;
    }
    auto pats = parseInput(inp);
    auto best = parseOutput(bestf);
    vector<int> idx;
    if (idxs == "all") for (int i = 0; i < (int)pats.size(); i++) idx.push_back(i);
    else { stringstream ss(idxs); string t; while (getline(ss, t, ',')) idx.push_back(stoi(t)); }
    map<int, Geo> geos;
    for (auto& p : pats) if (!geos.count(p.r)) geos[p.r] = makeGeo(p.r);
    for (auto& p : pats) initPat(p, geos[p.r]);
    vector<int> work;
    for (int i : idx) if (pats[i].np <= maxP) work.push_back(i);
    vector<string> outs(work.size());
#ifdef _OPENMP
    omp_set_num_threads(threads);
#endif
#pragma omp parallel for schedule(dynamic, 1)
    for (int wi = 0; wi < (int)work.size(); wi++) {
        int pi = work[wi];
        const Pat& pt = pats[pi];
        const Geo& g = geos.at(pt.r);
        PInfo P = makePInfo(g, pt);
        if (P.N < 0) continue;
        Solver S(P, opt);
        S.mbl = hasStatic(P) ? 1 : 2;
        int bl = best[pi].size();
        double ts = now();
        string out;
        char buf[256];
        if (mode == "first") {
            vector<Cfg> path;
            int f = S.firstBuild(min(bl, 12), &path);
            snprintf(buf, sizeof buf, "FIRST %d r=%d P=%d K=%d cnt=%d best=%d static=%d first=%d LB=%d nodes=%lld t=%.1f\n", pi, pt.r, P.nP, P.nK,
                     pt.cnt, bl, S.mbl == 1, f, f > 0 ? f + (pt.cnt - 1) * S.mbl : -1, S.totalNodes, now() - ts);
            out = buf;
            if (f > 0 && pt.cnt == 1 && f < bl) {
                vector<Cfg> tl; tl.push_back(startCfgP(P)); for (auto& x : path) tl.push_back(x);
                Prog prog = reconstruct(P, tl);
                int a = simulate(g, pt, prog);
                if (a >= pt.cnt) {
                    snprintf(buf, sizeof buf, "RES %d %zu\n", pi, prog.size()); out += buf;
                    for (auto& in : prog) { snprintf(buf, sizeof buf, "%c %c %c %c %c %c\n", in[0], in[1], in[2], in[3], in[4], in[5]); out += buf; }
                } else { snprintf(buf, sizeof buf, "VERIFYFAIL %d acc=%d\n", pi, a); out += buf; }
            }
        } else if (mode == "univ") {
            Cfg U; for (int e = 0; e < 6; e++) U[e] = P.UNIV;
            int found = -1; size_t nends = 0; bool ab = false;
            for (int len = 1; len <= 4; len++) {
                BuildSearch bs(P); bs.nodeLimit = opt.nodeLimit; bs.maxEnds = opt.maxEnds;
                bs.run(U, len); S.totalNodes += bs.nodes;
                if (!bs.ends.empty()) { found = len; nends = bs.ends.size(); ab = bs.aborted; break; }
                if (bs.aborted) { ab = true; found = -len; break; }
            }
            snprintf(buf, sizeof buf, "UNIV %d r=%d P=%d K=%d cnt=%d best=%d minbuild=%d ends=%zu%s nodes=%lld t=%.1f\n", pi, pt.r, P.nP, P.nK,
                     pt.cnt, bl, found, nends, ab ? "+" : "", S.totalNodes, now() - ts);
            out = buf;
        } else if (mode == "dp") {
            vector<Cfg> tl = S.dp(bl - 1);
            if (!tl.empty()) {
                Prog prog = reconstruct(P, tl);
                int a = simulate(g, pt, prog);
                snprintf(buf, sizeof buf, "DP %d r=%d P=%d cnt=%d best=%d new=%zu acc=%d nodes=%lld t=%.1f\n", pi, pt.r, P.nP, pt.cnt, bl, prog.size(), a, S.totalNodes, now() - ts);
                out = buf;
                if (a >= pt.cnt && (int)prog.size() < bl) {
                    snprintf(buf, sizeof buf, "RES %d %zu\n", pi, prog.size()); out += buf;
                    for (auto& in : prog) { snprintf(buf, sizeof buf, "%c %c %c %c %c %c\n", in[0], in[1], in[2], in[3], in[4], in[5]); out += buf; }
                }
            } else {
                snprintf(buf, sizeof buf, "DP %d r=%d P=%d cnt=%d best=%d none%s nodes=%lld t=%.1f\n", pi, pt.r, P.nP, pt.cnt, bl, S.timeUp ? "(timeout)" : "", S.totalNodes, now() - ts);
                out = buf;
            }
        }
        fputs(out.c_str(), stderr);
        outs[wi] = out;
    }
    for (auto& s : outs) fputs(s.c_str(), stdout);
}
