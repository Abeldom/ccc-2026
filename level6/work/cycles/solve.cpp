// Per pattern: find the minimal cycle rate L (cyclic SA with m builds per cycle), then the shortest
// first build T1 joint with that cycle; emit program = first build + cycle repeated.
// Usage: solve <input> <two-hints file> <out programs file> <iters> <restarts> [idx file]
#include "sa.h"

struct Hint { int stat = -1; long long two = 0; };

struct Result {
    int total = INT_MAX;
    int T1 = 0, L = 0, m = 0;
    vector<string> prog;
    string info;
};

template <int W>
Result solvePattern(const Pattern& p, const Hint& h, long long iters, int restarts, u64 seed) {
    Result res;
    SA<W> sa(p, seed);
    Sim sim(p);
    auto tryProblem = [&](const Problem& pr, int rs, long long it) -> bool {
        sa.setup(pr);
        for (int q = 0; q < rs; q++)
            if (sa.run(it, 2.0, 0.05) == 0) return true;
        return false;
    };
    auto buildProg = [&](int T1, int L, int m, const vector<Cfg>& cf) -> vector<Cfg> {
        vector<Cfg> seq(cf.begin(), cf.begin() + T1);
        for (int j = 0; j + 1 < p.cnt; j++) {
            int b = j % m;
            for (int t = 0; t < L; t++) seq.push_back(cf[T1 + b * L + t]);
        }
        return seq;
    };
    auto consider = [&](int T1, int L, int m, const vector<Cfg>& cf) {
        vector<Cfg> seq = buildProg(T1, L, m, cf);
        auto acc = sim.run(seq);
        if ((int)acc.size() >= p.cnt && (int)seq.size() < res.total) {
            // trim after the cnt-th acceptance
            seq.resize(acc[p.cnt - 1]);
            res.total = seq.size();
            res.T1 = T1; res.L = L; res.m = m;
            res.prog = toProgram(seq);
            return true;
        }
        return false;
    };
    if (p.cnt == 1) {
        for (int T1 = p.lb1; T1 <= p.lb1 + 10 && T1 < res.total; T1++) {
            Problem pr = makeJoint(T1, 0, 1);
            if (tryProblem(pr, restarts, iters)) { if (consider(T1, 0, 1, sa.bestCfgs())) break; }
        }
        return res;
    }
    // ---- rate search ----
    int Lstart = h.stat >= 0 ? 1 : (h.two > 0 ? 2 : 3);
    int Lbest = -1;
    vector<pair<int, int>> found;  // (L, m)
    for (int L = Lstart; L <= 12 && Lbest < 0; L++) {
        vector<int> ms = (L <= 2) ? vector<int>{1} : vector<int>{1, 2, 3};
        for (int m : ms) {
            Problem pr = makeJoint(0, L, m);
            if (tryProblem(pr, restarts, iters)) { found.push_back({L, m}); Lbest = L; }
        }
    }
    if (Lbest < 0) return res;
    // ---- joint prefix + cycle ----
    for (auto [L, m] : found) {
        for (int T1 = p.lb1; T1 <= p.lb1 + 8; T1++) {
            if (T1 + (long long)(p.cnt - 1) * L >= res.total) break;
            Problem pr = makeJoint(T1, L, m);
            if (tryProblem(pr, restarts, iters)) { if (consider(T1, L, m, sa.bestCfgs())) break; }
        }
    }
    return res;
}

int main(int argc, char** argv) {
    auto pats = readInput(argv[1]);
    int n = pats.size();
    vector<Hint> hints(n);
    {
        ifstream f(argv[2]);
        string line;
        while (getline(f, line)) {
            int idx; if (sscanf(line.c_str(), "%d", &idx) != 1 || idx < 0 || idx >= n) continue;
            auto gs = [&](const string& key) -> long long {
                auto pos = line.find(key + "=");
                if (pos == string::npos) return -1;
                return atoll(line.c_str() + pos + key.size() + 1);
            };
            hints[idx].stat = gs("static");
            hints[idx].two = gs("two");
            if (hints[idx].two < 0) hints[idx].two = 1;  // unknown: try L=2
        }
    }
    string outPath = argv[3];
    long long iters = atoll(argv[4]);
    int restarts = atoi(argv[5]);
    vector<int> idx;
    if (argc > 6) { ifstream f(argv[6]); int x; while (f >> x) idx.push_back(x); }
    else for (int i = 0; i < n; i++) idx.push_back(i);
    vector<Result> res(n);
    auto t0 = chrono::steady_clock::now();
#pragma omp parallel for schedule(dynamic, 1)
    for (int ii = 0; ii < (int)idx.size(); ii++) {
        int i = idx[ii];
        auto t1 = chrono::steady_clock::now();
        const Pattern& p = pats[i];
        Result r;
        if (p.N <= 64) r = solvePattern<1>(p, hints[i], iters, restarts, 1000 + i);
        else if (p.N <= 128) r = solvePattern<2>(p, hints[i], iters, restarts, 1000 + i);
        else if (p.N <= 192) r = solvePattern<3>(p, hints[i], iters, restarts, 1000 + i);
        else if (p.N <= 256) r = solvePattern<4>(p, hints[i], iters, restarts, 1000 + i);
        else if (p.N <= 384) r = solvePattern<6>(p, hints[i], iters, restarts, 1000 + i);
        else r = solvePattern<8>(p, hints[i], iters, restarts, 1000 + i);
        res[i] = r;
        double dt = chrono::duration<double>(chrono::steady_clock::now() - t1).count();
#pragma omp critical
        fprintf(stderr, "%d r=%d P=%d cnt=%d lb1=%d -> total=%d T1=%d L=%d m=%d (%.1fs)\n", i, p.r, p.nP, p.cnt, p.lb1,
                r.total == INT_MAX ? -1 : r.total, r.T1, r.L, r.m, dt);
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
