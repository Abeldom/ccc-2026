// Test backward beam: cycle mode per pattern; prints rate and verified total with simple prefix.
// Usage: bbtest <input> <idx file|all> <B> <wD> <wK> [chain]
#include "bbeam.h"

static vector<Cfg> simplePrefix(const Cfg& C1) {
    int p = 0;
    for (int e = 0; e < 6; e++) { int k = abs(C1[e].k); p = max(p, C1[e].a ? k : max(0, k - 1)); }
    vector<Cfg> seq;
    Cfg cur = initCfg();
    for (int t = 0; t < p; t++) {
        for (int e = 0; e < 6; e++) { if (cur[e].k < C1[e].k) cur[e].k++; else if (cur[e].k > C1[e].k) cur[e].k--; }
        seq.push_back(cur);
    }
    return seq;
}

int main(int argc, char** argv) {
    auto pats = readInput(argv[1]);
    vector<int> idx;
    if (string(argv[2]) != "all") { ifstream f(argv[2]); int x; while (f >> x) idx.push_back(x); }
    else for (int i = 0; i < (int)pats.size(); i++) idx.push_back(i);
    int B = atoi(argv[3]);
    double wD = atof(argv[4]), wK = atof(argv[5]);
    bool chain = argc > 6;
#pragma omp parallel for schedule(dynamic, 1)
    for (int ii = 0; ii < (int)idx.size(); ii++) {
        int i = idx[ii];
        const Pattern& p = pats[i];
        auto t0 = chrono::steady_clock::now();
        int rate = -1, total = -1, acc = -1;
        auto go = [&](auto* tag) {
            constexpr int W = std::remove_pointer_t<decltype(tag)>::value;
            BBeam<W> bb(p);
            typename BBeam<W>::Opt o;
            o.B = B; o.wD = wD; o.wK = wK; o.cycle = !chain;
            for (int e = 0; e < 6; e++) o.X[e] = bb.G.id(ES{0, 0});
            auto res = bb.search(o);
            if (res.len == INT_MAX) return;
            rate = res.len;
            Sim sim(p);
            vector<Cfg> seq;
            if (!chain) {
                seq = simplePrefix(res.build[0]);
                for (int j = 0; j < p.cnt; j++) seq.insert(seq.end(), res.build.begin(), res.build.end());
            } else seq = res.build;
            auto a = sim.run(seq);
            acc = a.size();
            if ((int)a.size() >= (chain ? 1 : p.cnt)) total = a[(chain ? 1 : p.cnt) - 1];
        };
        if (p.N <= 64) go((integral_constant<int, 1>*)nullptr);
        else if (p.N <= 128) go((integral_constant<int, 2>*)nullptr);
        else go((integral_constant<int, 3>*)nullptr);
        double dt = chrono::duration<double>(chrono::steady_clock::now() - t0).count();
#pragma omp critical
        printf("%d r=%d P=%d K=%d cnt=%d rate=%d acc=%d total=%d (%.1fs)\n", i, p.r, p.nP, p.nK, p.cnt, rate, acc, total, dt);
        fflush(stdout);
    }
}
