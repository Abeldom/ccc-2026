// solve3: per pattern (count >= 2): descending search for cycles of general segment shapes (mean rate),
// warm-started from the best cycle found so far; then prefix SA (first build ending adjacent to a cycle
// build start) and tail SA (shorter free last build). Writes verified programs incrementally.
// Usage: solve3 <input> <best_lengths file> <hints file|none> <out sol> <iters per tick> <restarts>
//               <time per pattern s> <idx file|all> [seed sol files (prior results, "# idx" format)...]
#include "sa3.h"

struct Cycle {
    vector<int> lens;
    vector<Cfg> cf;
    double rate() const { return (double)cf.size() / lens.size(); }
    int start(int b) const { int t = 0; for (int i = 0; i < b; i++) t += lens[i]; return t; }
};

static double nowS() { return chrono::duration<double>(chrono::steady_clock::now().time_since_epoch()).count(); }

static vector<Cfg> simplePrefix(const Cfg& C1) {
    int p = 0;
    for (int e = 0; e < 6; e++) {
        int k = abs(C1[e].k);
        p = max(p, C1[e].a ? k : max(0, k - 1));
    }
    vector<Cfg> seq;
    Cfg cur = initCfg();
    for (int t = 0; t < p; t++) {
        for (int e = 0; e < 6; e++) {
            if (cur[e].k < C1[e].k) cur[e].k++;
            else if (cur[e].k > C1[e].k) cur[e].k--;
        }
        seq.push_back(cur);
    }
    return seq;
}

struct Opts { long long ipt; int restarts; double tlim; };

template <int W>
struct Solver {
    const Pattern& p;
    Opts o;
    SA<W> sa;
    Sim sim;
    double t0;
    string log;
    int ub;
    vector<Cycle> pool;
    int bestTotal = INT_MAX;
    vector<Cfg> bestSeq;

    Solver(const Pattern& pat, Opts op, u64 seed, int ub_) : p(pat), o(op), sa(pat, seed), sim(pat), ub(ub_) { t0 = nowS(); }
    double frac = 1.0;
    bool timeUp() const { return nowS() - t0 > o.tlim * frac; }

    bool consider(const vector<Cfg>& seq0) {
        Cfg prev = initCfg();
        for (auto& c : seq0) { if (!adjCfg(prev, c)) return false; prev = c; }
        auto acc = sim.run(seq0);
        if ((int)acc.size() < p.cnt) return false;
        int len = acc[p.cnt - 1];
        if (len >= bestTotal) return false;
        bestTotal = len;
        bestSeq.assign(seq0.begin(), seq0.begin() + len);
        return true;
    }
    // verify a cycle by simulation: run it 3 times after a canonical start... simplest: check each segment from
    // canonical with the exact cost model (sa) -> trust model; final programs are verified by Sim anyway.
    bool searchShape(const vector<int>& lens, const Cycle* warm, int restarts, Cycle& out) {
        Problem pr = makeCycle(lens);
        sa.setup(pr);
        long long iters = o.ipt * pr.T;
        for (int q = 0; q < restarts && !timeUp(); q++) {
            bool useWarm = warm && (q % 2 == 0);
            if (useWarm) {
                // repeat warm cycle segments cyclically, then fit each segment length by deleting/duplicating ticks
                vector<Cfg> v;
                int m0 = warm->lens.size();
                for (size_t i = 0; i < lens.size(); i++) {
                    int b = (i + q / 2) % m0;
                    int st = warm->start(b), L0 = warm->lens[b];
                    vector<Cfg> seg(warm->cf.begin() + st, warm->cf.begin() + st + L0);
                    while ((int)seg.size() > lens[i]) seg.erase(seg.begin() + sa.rng.rint(seg.size()));
                    while ((int)seg.size() < lens[i]) { int k = sa.rng.rint(seg.size()); seg.insert(seg.begin() + k, seg[k]); }
                    v.insert(v.end(), seg.begin(), seg.end());
                }
                sa.loadCfgs(v);
                if (!sa.repair()) { if (!sa.initRandom()) continue; }
                if (sa.anneal(iters, 0.8, 0.05) == 0) { out = {lens, sa.bestCfgs()}; return true; }
            } else {
                if (sa.run(iters, 2.0, 0.05) == 0) { out = {lens, sa.bestCfgs()}; return true; }
            }
        }
        return false;
    }
    // all shapes with m<=3 and rate in [lo, hi), lens >= Lmin, sorted by rate descending
    static vector<vector<int>> shapes(double lo, double hi, int Lmin, int maxM = 3) {
        vector<pair<double, vector<int>>> v;
        int Lmax = (int)ceil(hi) + 2;
        for (int a = Lmin; a <= Lmax; a++) {
            if (a >= lo - 1e-9 && a < hi - 1e-9) v.push_back({(double)a, {a}});
            if (maxM < 2) continue;
            for (int b = a; b <= Lmax; b++) {
                double rt = (a + b) / 2.0;
                if (b > a && rt >= lo - 1e-9 && rt < hi - 1e-9) v.push_back({rt, {a, b}});
                for (int c = b; c <= Lmax && maxM >= 3; c++) {
                    double r3 = (a + b + c) / 3.0;
                    if ((b > a || c > b) && r3 >= lo - 1e-9 && r3 < hi - 1e-9) {
                        v.push_back({r3, {a, b, c}});
                        if (a < b && b < c) v.push_back({r3, {a, c, b}});
                    }
                }
            }
        }
        // also uniform m=2 / m=3 of integer rates (different from m=1: more freedom)
        for (int a = Lmin; a <= Lmax; a++)
            if (a >= lo - 1e-9 && a < hi - 1e-9) {
                if (maxM >= 2) v.push_back({a - 1e-6, {a, a}});
                if (maxM >= 3) v.push_back({a - 2e-6, {a, a, a}});
            }
        sort(v.begin(), v.end(), [](auto& x, auto& y) { return x.first > y.first; });
        vector<vector<int>> res;
        for (auto& x : v) res.push_back(x.second);
        return res;
    }

    bool addCycle(const vector<int>& lens, const vector<Cfg>& cf, bool allowRepair) {
        Problem pr = makeCycle(lens);
        sa.setup(pr);
        sa.loadCfgs(cf);
        if (sa.total == 0 && sa.validAll()) { pool.push_back({lens, cf}); return true; }
        if (!allowRepair) return false;
        if (!sa.repair()) return false;
        if (sa.anneal(o.ipt * pr.T / 2, 0.5, 0.05) == 0) { pool.push_back({lens, sa.bestCfgs()}); return true; }
        return false;
    }
    vector<Cfg> seedFirst;
    void seedFromProgram(const vector<Cfg>& seq) {
        auto acc = sim.run(seq);
        int c = min((int)acc.size(), p.cnt);
        if (c < 2) return;
        seedFirst.assign(seq.begin(), seq.begin() + acc[0]);
        if (c >= p.cnt) consider(seq);
        int added = 0;
        for (int m = 1; m <= 3; m++)
            for (int j = 1; j + m <= c && added < 12; j++) {
                int st = acc[j - 1], en = acc[j + m - 1];
                vector<int> lens;
                for (int q = j; q < j + m; q++) lens.push_back(acc[q] - acc[q - 1]);
                vector<Cfg> cf(seq.begin() + st, seq.begin() + en);
                if (addCycle(lens, cf, true)) added++;
            }
        char buf[32];
        snprintf(buf, sizeof buf, " seedP%d", added);
        log += buf;
    }
    void cycleSearch(int Lmin, double Rlim) {
        // initial: smallest uniform L found descending from floor(Rlim)
        Cycle best;
        bool have = false;
        for (auto& c : pool) if (!have || c.rate() < best.rate()) { best = c; have = true; }
        if (!have) {
            // ascending search for a first cycle (m=1), fast restarts
            for (int L = max(Lmin, (int)floor(Rlim * 0.6)); L <= max(Lmin, (int)ceil(Rlim)) + 1 && !timeUp(); L++) {
                Cycle c;
                if (searchShape({L}, nullptr, 2, c)) { pool.push_back(c); best = c; have = true; break; }
            }
            if (!have) {
                // descend from above
                for (int L = (int)ceil(Rlim) + 2; L <= max(14, (int)ceil(Rlim) + 6) && !timeUp(); L++) {
                    Cycle c;
                    if (searchShape({L}, nullptr, 2, c)) { pool.push_back(c); best = c; have = true; break; }
                }
            }
            if (!have) return;
        }
        // descend (on failure retry the window with doubled restarts until the time is up)
        int rs = o.restarts;
        while (!timeUp()) {
            double rb = best.rate();
            if (rb <= Lmin + 1e-9) break;
            int maxM = rb <= 6.01 ? 3 : (rb <= 10.01 ? 2 : 1);
            auto sh = shapes(max((double)Lmin, rb - 1.0), rb, Lmin, maxM);
            if (sh.empty()) break;
            bool improved = false;
            for (auto& lens : sh) {
                if (timeUp()) break;
                Cycle c;
                if (searchShape(lens, &best, rs, c)) {
                    pool.push_back(c);
                    best = c;
                    improved = true;
                    char buf[64];
                    snprintf(buf, sizeof buf, " cyc%.2f", c.rate());
                    log += buf;
                    break;
                }
            }
            if (!improved) rs *= 2; else rs = o.restarts;
        }
        // extra diversity at the best rate
        for (int q = 0; q < 2 && !timeUp(); q++) {
            Cycle c;
            if (searchShape(best.lens, &best, 2, c)) pool.push_back(c);
        }
    }

    vector<Cfg> assemble(const vector<Cfg>& pre, const Cycle& c, int b0, int nb) {
        vector<Cfg> seq = pre;
        int m = c.lens.size();
        for (int j = 0; j < nb; j++) {
            int b = (b0 + j) % m;
            int st = c.start(b);
            for (int t = 0; t < c.lens[b]; t++) seq.push_back(c.cf[st + t]);
        }
        return seq;
    }
    // shortest build from config X (free end) of length < maxL; returns empty if none found
    vector<Cfg> tailBuild(const Cfg& X, int maxL, int minL) {
        vector<Cfg> bestv;
        for (int l = maxL - 1; l >= max(1, minL); l--) {
            if (timeUp()) break;
            Problem pr = makeChain({l}, X);
            sa.setup(pr);
            bool ok = false;
            for (int q = 0; q < o.restarts && !ok; q++)
                if (sa.run(o.ipt * l, 2.0, 0.05) == 0) { ok = true; bestv = sa.bestCfgs(); }
            if (!ok) break;
        }
        return bestv;
    }
    void prefixPhase(int Lmin) {
        if (pool.empty()) return;
        double rmin = 1e9;
        for (auto& c : pool) rmin = min(rmin, c.rate());
        // candidates: cycles with rate within slack (bigger slack for small counts)
        double slack = p.cnt <= 3 ? 2.0 : (p.cnt <= 6 ? 1.0 : 0.5);
        struct Cand { int est; int ci, b; vector<Cfg> pre; };
        vector<Cand> cands;
        for (int ci = 0; ci < (int)pool.size(); ci++) {
            const Cycle& c = pool[ci];
            if (c.rate() > rmin + slack + 1e-9) continue;
            for (int b = 0; b < (int)c.lens.size(); b++) {
                vector<Cfg> pre = simplePrefix(c.cf[c.start(b)]);
                for (int t = 0; t < c.lens[b]; t++) pre.push_back(c.cf[c.start(b) + t]);
                auto seq = assemble(pre, c, (b + 1) % c.lens.size(), p.cnt - 1);
                consider(seq);
                int est = (int)pre.size() + (int)llround((p.cnt - 1) * c.rate());
                cands.push_back({est, ci, b, pre});
            }
        }
        sort(cands.begin(), cands.end(), [](const Cand& a, const Cand& b) { return a.est < b.est; });
        if (cands.size() > 6) cands.resize(6);
        for (auto& cd : cands) {
            if (timeUp()) break;
            const Cycle& c = pool[cd.ci];
            int m = c.lens.size();
            // descending prefix length: first build ending adjacent to build cd.b's first config
            vector<Cfg> bestPre;  // prefix that ends right before build b (accepting at its end)
            int simpleT1 = cd.pre.size();
            // option A: prefix + cycle starting at b+1 (the simple one used build b itself as first build)
            // option B: SA prefix of length T1 ending adjacent to c.cf[start(b)], then cycle from b
            for (int T1 = simpleT1 - 1; T1 >= max(1, p.lb1); T1--) {
                if (timeUp()) break;
                int est = T1 + (int)llround((p.cnt - 1) * c.rate());
                if (est > bestTotal + 1) break;
                Problem pr = makeChain({T1}, initCfg());
                pr.fixEnd = true;
                pr.endCfg = c.cf[c.start(cd.b)];
                sa.setup(pr);
                bool ok = false;
                for (int q = 0; q < o.restarts && !ok; q++)
                    if (sa.run(o.ipt * T1, 2.0, 0.05) == 0) { ok = true; bestPre = sa.bestCfgs(); }
                if (!ok) break;
                consider(assemble(bestPre, c, cd.b, p.cnt - 1));
            }
            (void)m;
        }
        // bridge: seed program's first build, then one free build ending adjacent to a cycle build start
        if (!seedFirst.empty() && p.cnt >= 3) {
            Cfg X1 = seedFirst.back();
            for (auto& cd : cands) {
                if (timeUp()) break;
                const Cycle& c = pool[cd.ci];
                int m = c.lens.size();
                Cfg Cs = c.cf[c.start(cd.b)];
                if (adjCfg(X1, Cs)) consider(assemble(seedFirst, c, cd.b, p.cnt - 1));
                int maxl = c.lens[(cd.b + m - 1) % m] + 3;
                vector<Cfg> br;
                for (int l = maxl; l >= max(1, Lmin); l--) {
                    if (timeUp()) break;
                    Problem pr = makeChain({l}, X1);
                    pr.fixEnd = true;
                    pr.endCfg = Cs;
                    sa.setup(pr);
                    bool ok = false;
                    for (int q = 0; q < o.restarts && !ok; q++)
                        if (sa.run(o.ipt * l, 2.0, 0.05) == 0) { ok = true; br = sa.bestCfgs(); }
                    if (!ok) break;
                    vector<Cfg> pre = seedFirst;
                    pre.insert(pre.end(), br.begin(), br.end());
                    if (consider(assemble(pre, c, cd.b, p.cnt - 2))) log += " bridge";
                }
            }
        }
        // tail: shorten the last build of the current best program
        if (!bestSeq.empty() && p.cnt >= 2 && !timeUp()) {
            auto acc = sim.run(bestSeq);
            int s0 = acc[p.cnt - 2];  // ticks before the last build
            int lastLen = acc[p.cnt - 1] - s0;
            Cfg X = s0 > 0 ? bestSeq[s0 - 1] : initCfg();
            auto tb = tailBuild(X, lastLen, Lmin);
            if (!tb.empty()) {
                vector<Cfg> seq(bestSeq.begin(), bestSeq.begin() + s0);
                seq.insert(seq.end(), tb.begin(), tb.end());
                if (consider(seq)) log += " tail";
            }
        }
    }
};

static map<int, vector<vector<string>>> readSols(const string& path) {
    map<int, vector<vector<string>>> res;
    ifstream f(path);
    string line;
    int cur = -1;
    while (getline(f, line)) {
        if (line.empty()) continue;
        if (line[0] == '#') { cur = atoi(line.c_str() + 2); res[cur].push_back({}); continue; }
        if (cur >= 0) res[cur].back().push_back(line);
    }
    return res;
}

static vector<vector<vector<string>>> readOutFile(const string& path) {
    vector<vector<vector<string>>> progs;
    vector<vector<string>> cur;
    ifstream f(path);
    string line;
    while (getline(f, line)) {
        stringstream ss(line);
        vector<string> toks;
        string t;
        while (ss >> t) toks.push_back(t);
        if (toks.empty()) { if (!cur.empty()) { progs.push_back(cur); cur.clear(); } continue; }
        cur.push_back(toks);
    }
    if (!cur.empty()) progs.push_back(cur);
    return progs;
}
static bool progToCfgs(const vector<vector<string>>& prog, int r, vector<Cfg>& seq) {
    Cfg cur = initCfg();
    seq.clear();
    for (auto& line : prog) {
        if (line.size() != 6) return false;
        for (int e = 0; e < 6; e++) {
            char c = line[e][0];
            if (c == 'F') cur[e].k++;
            else if (c == 'B') cur[e].k--;
            else if (c == 'T') cur[e].a ^= 1;
        }
        for (int e = 0; e < 6; e++) if (abs(cur[e].k) > r) return false;
        seq.push_back(cur);
    }
    return true;
}
// exhL output lines: "<idx> ... L=<L> found=1 ... | k,a k,a ... ; (6 groups, sign +1 convention)"
static vector<Cycle> parseSeeds(const string& path, int want) {
    vector<Cycle> out;
    ifstream f(path);
    string line;
    while (getline(f, line)) {
        int idx;
        if (sscanf(line.c_str(), "%d", &idx) != 1 || idx != want) continue;
        auto bar = line.find('|');
        if (bar == string::npos) continue;
        auto pos = line.find("L=");
        if (pos == string::npos) continue;
        int L = atoi(line.c_str() + pos + 2);
        stringstream ss(line.substr(bar + 1));
        string tok;
        vector<vector<ES>> groups;
        vector<ES> cur;
        while (ss >> tok) {
            if (tok == ";") { groups.push_back(cur); cur.clear(); continue; }
            int k, a;
            if (sscanf(tok.c_str(), "%d,%d", &k, &a) == 2) cur.push_back(ES{(int8_t)k, (int8_t)a});
        }
        if (groups.size() != 6) continue;
        Cycle c{{L}, vector<Cfg>(L)};
        bool ok = true;
        for (int a = 0; a < 3; a++)
            for (int q = 0; q < 2; q++) {
                int e = a + 3 * q;
                auto& g = groups[2 * a + q];
                if ((int)g.size() != L) { ok = false; break; }
                for (int t = 0; t < L; t++) c.cf[t][e] = ES{(int8_t)(SGN[e] * g[t].k), g[t].a};
            }
        if (ok) out.push_back(c);
    }
    return out;
}

int main(int argc, char** argv) {
    if (argc < 9) { fprintf(stderr, "usage\n"); return 1; }
    auto pats = readInput(argv[1]);
    int n = pats.size();
    vector<int> ub(n, INT_MAX);
    {
        ifstream f(argv[2]);
        string line;
        while (getline(f, line)) {
            if (line.empty() || line[0] == '#') continue;
            stringstream ss(line);
            int i, r, P, c, b;
            if (ss >> i >> r >> P >> c >> b) if (i >= 0 && i < n) ub[i] = b;
        }
    }
    vector<int> Lmin(n, 2);
    if (string(argv[3]) != "none") {
        ifstream f(argv[3]);
        string line;
        while (getline(f, line)) {
            int idx;
            if (sscanf(line.c_str(), "%d", &idx) != 1 || idx < 0 || idx >= n) continue;
            auto gs = [&](const string& key) -> long long {
                auto pos = line.find(key + "=");
                if (pos == string::npos) return -2;
                return atoll(line.c_str() + pos + key.size() + 1);
            };
            long long s = gs("static"), t = gs("two");
            if (s >= 0) Lmin[idx] = 1;
            else if (t > 0) Lmin[idx] = 2;
            else if (t == 0) Lmin[idx] = 3;
        }
    }
    string outPath = argv[4];
    Opts o{atoll(argv[5]), atoi(argv[6]), atof(argv[7])};
    vector<int> idx;
    if (string(argv[8]) != "all") { ifstream f(argv[8]); int x; while (f >> x) idx.push_back(x); }
    else for (int i = 0; i < n; i++) idx.push_back(i);
    vector<vector<vector<string>>> seedProgs;
    if (argc > 9 && string(argv[9]) != "none") seedProgs = readOutFile(argv[9]);
    vector<string> exhFiles;
    for (int a = 10; a < argc; a++) exhFiles.push_back(argv[a]);
    FILE* fo = fopen(outPath.c_str(), "a");
    double T0 = nowS();
#pragma omp parallel for schedule(dynamic, 1)
    for (int ii = 0; ii < (int)idx.size(); ii++) {
        int i = idx[ii];
        const Pattern& p = pats[i];
        if (p.cnt < 2) continue;
        double t1 = nowS();
        int total = -1;
        string log;
        vector<string> prog;
        auto go = [&](auto* tag) {
            constexpr int W = std::remove_pointer_t<decltype(tag)>::value;
            Solver<W> s(p, o, 1234567ULL + i * 7919ULL, ub[i]);
            double Rlim = ub[i] == INT_MAX ? 10.0 : (double)(ub[i] - p.lb1) / (p.cnt - 1);
            s.frac = 0.55;
            for (auto& ef : exhFiles) for (auto& c : parseSeeds(ef, i)) s.addCycle(c.lens, c.cf, false);
            if ((int)seedProgs.size() == n) {
                vector<Cfg> sq;
                if (progToCfgs(seedProgs[i], p.r, sq)) s.seedFromProgram(sq);
            }
            s.cycleSearch(Lmin[i], Rlim);
            s.frac = 1.0;
            s.prefixPhase(Lmin[i]);
            log = s.log;
            double rmin = 0;
            if (!s.pool.empty()) { rmin = 1e9; for (auto& c : s.pool) rmin = min(rmin, c.rate()); }
            char buf[64];
            snprintf(buf, sizeof buf, " rate=%.2f pool=%zu", rmin, s.pool.size());
            log += buf;
            if (!s.bestSeq.empty()) { total = s.bestTotal; prog = toProgram(s.bestSeq); }
        };
        if (p.N <= 64) go((integral_constant<int, 1>*)nullptr);
        else if (p.N <= 128) go((integral_constant<int, 2>*)nullptr);
        else if (p.N <= 192) go((integral_constant<int, 3>*)nullptr);
        else if (p.N <= 256) go((integral_constant<int, 4>*)nullptr);
        else go((integral_constant<int, 6>*)nullptr);
#pragma omp critical
        {
            fprintf(stderr, "%d r=%d P=%d cnt=%d lb1=%d ub=%d -> %d%s %s (%.0fs, at %.0fs)\n", i, p.r, p.nP, p.cnt, p.lb1, ub[i], total,
                    log.c_str(), (total > 0 && total < ub[i]) ? "WIN" : "", nowS() - t1, nowS() - T0);
            if (total > 0) {
                fprintf(fo, "# %d total=%d\n", i, total);
                for (auto& l : prog) fprintf(fo, "%s\n", l.c_str());
                fflush(fo);
            }
        }
    }
    fclose(fo);
}
