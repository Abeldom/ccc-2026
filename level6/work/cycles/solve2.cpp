// Per pattern: (1) collect cycles of minimal rate L (seeds from exhaustive search files + cyclic SA with
// m = 1..3 builds per cycle), (2) for every cycle (and rotation) build the simple prefix and try to
// shorten the first build with a prefix SA that must end adjacent to the cycle's first config,
// (3) emit the shortest verified program.
// Usage: solve2 <input> <hints file> <out sol file> <iters> <restarts> <idx file or "all"> [seed files...]
#include "sa.h"

struct Hint { int stat = -1; long long two = 1; };
struct Cycle { int L, m; vector<Cfg> cf; };
struct Result { int total = INT_MAX; int T1 = 0, L = 0, m = 0; vector<string> prog; };

static vector<Cfg> simplePrefix(const Cfg& C1) {
    int p = 0;
    for (int e = 0; e < 6; e++) {
        int k = abs(C1[e].k);
        p = max(p, C1[e].a ? k : max(0, k - 1));
    }
    vector<Cfg> seq;
    Cfg cur;
    for (int e = 0; e < 6; e++) cur[e] = {0, 0};
    for (int t = 0; t < p; t++) {
        for (int e = 0; e < 6; e++) {
            if (cur[e].k < C1[e].k) cur[e].k++;
            else if (cur[e].k > C1[e].k) cur[e].k--;
        }
        seq.push_back(cur);
    }
    return seq;
}

template <int W>
Result solvePattern(const Pattern& p, const Hint& h, long long iters, int restarts, u64 seed, const vector<Cycle>& seeds,
                    string& log) {
    Result res;
    SA<W> sa(p, seed);
    Sim sim(p);
    auto tryProblem = [&](const Problem& pr, int rs, long long it) -> bool {
        sa.setup(pr);
        for (int q = 0; q < rs; q++)
            if (sa.run(it, 2.0, 0.05) == 0) return true;
        return false;
    };
    auto consider = [&](const vector<Cfg>& seq0, int T1, int L, int m) -> bool {
        auto acc = sim.run(seq0);
        if ((int)acc.size() < p.cnt) return false;
        int len = acc[p.cnt - 1];
        if (len >= res.total) return false;
        vector<Cfg> seq(seq0.begin(), seq0.begin() + len);
        res.total = len; res.T1 = T1; res.L = L; res.m = m;
        res.prog = toProgram(seq);
        return true;
    };
    // program: first build `pre` then cycle builds starting at build index b0
    auto assemble = [&](const vector<Cfg>& pre, const Cycle& c, int b0) {
        vector<Cfg> seq = pre;
        for (int j = 0; j + 1 < p.cnt; j++) {
            int b = (b0 + j) % c.m;
            for (int t = 0; t < c.L; t++) seq.push_back(c.cf[b * c.L + t]);
        }
        return seq;
    };
    // ---- first build only (cnt == 1) ----
    if (p.cnt == 1) {
        for (const Cycle& c : seeds) {  // simple: prefix + one build of any seed cycle
            for (int b = 0; b < c.m; b++) {
                vector<Cfg> seq = simplePrefix(c.cf[b * c.L]);
                for (int t = 0; t < c.L; t++) seq.push_back(c.cf[b * c.L + t]);
                consider(seq, seq.size(), 0, 1);
            }
        }
        for (int T1 = p.lb1; T1 < res.total; T1++) {
            Problem pr = makeJoint(T1, 0, 1);
            if (tryProblem(pr, restarts, iters)) { if (consider(sa.bestCfgs(), T1, 0, 1)) break; }
        }
        return res;
    }
    // ---- rate search ----
    int Lstart = h.stat >= 0 ? 1 : (h.two > 0 ? 2 : 3);
    vector<Cycle> cyc;
    int Lseed = INT_MAX;
    for (auto& c : seeds) Lseed = min(Lseed, c.L);
    int Lbest = -1;
    for (int L = Lstart; L <= 12 && Lbest < 0; L++) {
        if (L >= Lseed) { Lbest = L; break; }
        vector<int> ms = (L <= 2) ? vector<int>{1} : vector<int>{1, 2, 3};
        for (int m : ms) {
            Problem pr = makeJoint(0, L, m);
            sa.setup(pr);
            int got = 0;
            for (int q = 0; q < restarts && got < 3; q++)
                if (sa.run(iters, 2.0, 0.05) == 0) { cyc.push_back({L, m, sa.bestCfgs()}); got++; Lbest = L; }
        }
    }
    for (auto& c : seeds) if (Lbest < 0 || c.L <= Lbest) { cyc.push_back(c); Lbest = (Lbest < 0) ? c.L : min(Lbest, c.L); }
    if (Lbest < 0) return res;
    // extra cycles at Lbest for diversity
    {
        vector<int> ms = (Lbest <= 2) ? vector<int>{1} : vector<int>{1, 2, 3};
        for (int m : ms) {
            Problem pr = makeJoint(0, Lbest, m);
            sa.setup(pr);
            for (int q = 0; q < 2; q++) if (sa.run(iters, 2.0, 0.05) == 0) cyc.push_back({Lbest, m, sa.bestCfgs()});
        }
    }
    char buf[256];
    snprintf(buf, sizeof buf, " Lbest=%d cycles=%zu", Lbest, cyc.size());
    log += buf;
    // ---- prefixes ----
    struct Cand { int simpleT1; int ci, b; };
    vector<Cand> cands;
    for (int ci = 0; ci < (int)cyc.size(); ci++) {
        const Cycle& c = cyc[ci];
        if (c.L > Lbest) continue;
        for (int b = 0; b < c.m; b++) {
            const Cfg& C1 = c.cf[b * c.L];
            vector<Cfg> pre = simplePrefix(C1);
            int T1 = pre.size() + c.L;
            for (int t = 0; t < c.L; t++) pre.push_back(c.cf[b * c.L + t]);
            // after this first build the cycle continues with build b+1
            consider(assemble(pre, c, (b + 1) % c.m), T1, c.L, c.m);
            cands.push_back({T1, ci, b});
        }
    }
    sort(cands.begin(), cands.end(), [](const Cand& a, const Cand& b) { return a.simpleT1 < b.simpleT1; });
    if (cands.size() > 8) cands.resize(8);
    for (auto& cd : cands) {
        const Cycle& c = cyc[cd.ci];
        for (int T1 = p.lb1; T1 < cd.simpleT1; T1++) {
            if (T1 + (long long)(p.cnt - 1) * c.L >= res.total) break;
            Problem pr = makeJoint(T1, 0, 1);
            pr.fixEnd = true;
            pr.endCfg = c.cf[cd.b * c.L];
            if (tryProblem(pr, restarts, iters)) {
                vector<Cfg> pre = sa.bestCfgs();
                if (consider(assemble(pre, c, cd.b), T1, c.L, c.m)) break;
            }
        }
    }
    return res;
}

static vector<Cycle> parseSeeds(const string& path, int want, int r) {
    // lines from exhL: "<idx> ... L=<L> found=1 ... | k,a k,a ... ; (6 groups)"
    vector<Cycle> out;
    ifstream f(path);
    string line;
    while (getline(f, line)) {
        int idx;
        if (sscanf(line.c_str(), "%d", &idx) != 1 || idx != want) continue;
        auto bar = line.find('|');
        if (bar == string::npos) continue;
        auto pos = line.find("L=");
        int L = atoi(line.c_str() + pos + 2);
        string rest = line.substr(bar + 1);
        vector<vector<ES>> groups;
        stringstream ss(rest);
        string tok;
        vector<ES> cur;
        while (ss >> tok) {
            if (tok == ";") { groups.push_back(cur); cur.clear(); continue; }
            int k, a;
            sscanf(tok.c_str(), "%d,%d", &k, &a);
            cur.push_back(ES{(int8_t)k, (int8_t)a});
        }
        if (groups.size() != 6) continue;
        Cycle c{L, 1, vector<Cfg>(L)};
        for (int a = 0; a < 3; a++)
            for (int q = 0; q < 2; q++) {
                int e = a + 3 * q;
                auto& g = groups[2 * a + q];
                for (int t = 0; t < L; t++) c.cf[t][e] = ES{(int8_t)(SGN[e] * g[t].k), g[t].a};
            }
        (void)r;
        out.push_back(c);
    }
    return out;
}

int main(int argc, char** argv) {
    auto pats = readInput(argv[1]);
    int n = pats.size();
    vector<Hint> hints(n);
    {
        ifstream f(argv[2]);
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
            hints[idx].stat = (s == -2) ? -1 : s;
            hints[idx].two = (t == -2) ? 1 : t;
        }
    }
    string outPath = argv[3];
    long long iters = atoll(argv[4]);
    int restarts = atoi(argv[5]);
    vector<int> idx;
    if (string(argv[6]) != "all") { ifstream f(argv[6]); int x; while (f >> x) idx.push_back(x); }
    else for (int i = 0; i < n; i++) idx.push_back(i);
    vector<string> seedFiles;
    for (int a = 7; a < argc; a++) seedFiles.push_back(argv[a]);
    vector<Result> res(n);
    auto t0 = chrono::steady_clock::now();
#pragma omp parallel for schedule(dynamic, 1)
    for (int ii = 0; ii < (int)idx.size(); ii++) {
        int i = idx[ii];
        auto t1 = chrono::steady_clock::now();
        const Pattern& p = pats[i];
        vector<Cycle> seeds;
        for (auto& sf : seedFiles) { auto v = parseSeeds(sf, i, p.r); seeds.insert(seeds.end(), v.begin(), v.end()); }
        string log;
        Result r;
        u64 sd = 1000 + i;
        if (p.N <= 64) r = solvePattern<1>(p, hints[i], iters, restarts, sd, seeds, log);
        else if (p.N <= 128) r = solvePattern<2>(p, hints[i], iters, restarts, sd, seeds, log);
        else if (p.N <= 192) r = solvePattern<3>(p, hints[i], iters, restarts, sd, seeds, log);
        else if (p.N <= 256) r = solvePattern<4>(p, hints[i], iters, restarts, sd, seeds, log);
        else if (p.N <= 384) r = solvePattern<6>(p, hints[i], iters, restarts, sd, seeds, log);
        else r = solvePattern<8>(p, hints[i], iters, restarts, sd, seeds, log);
        res[i] = r;
        double dt = chrono::duration<double>(chrono::steady_clock::now() - t1).count();
#pragma omp critical
        fprintf(stderr, "%d r=%d P=%d cnt=%d lb1=%d seeds=%zu -> total=%d T1=%d L=%d m=%d%s (%.1fs)\n", i, p.r, p.nP, p.cnt, p.lb1,
                seeds.size(), r.total == INT_MAX ? -1 : r.total, r.T1, r.L, r.m, log.c_str(), dt);
    }
    FILE* f = fopen(outPath.c_str(), "w");
    for (int i : idx) {
        if (res[i].total == INT_MAX) continue;
        fprintf(f, "# %d T1=%d L=%d m=%d total=%d\n", i, res[i].T1, res[i].L, res[i].m, res[i].total);
        for (auto& s : res[i].prog) fprintf(f, "%s\n", s.c_str());
    }
    fclose(f);
    fprintf(stderr, "done %.1fs\n", chrono::duration<double>(chrono::steady_clock::now() - t0).count());
}
