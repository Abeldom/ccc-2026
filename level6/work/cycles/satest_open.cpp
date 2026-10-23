// Test: for each pattern find the minimal cycle length L (pure cyclic SA) in [Lmin, Lmax].
#include "sa.h"

template <int W>
int minCycle(const Pattern& p, int Lmin, int Lmax, long long iters, int restarts, u64 seed, vector<Cfg>& sol, int m) {
    SA<W> sa(p, seed);
    for (int L = Lmin; L <= Lmax; L++) {
        Problem pr = makeJoint(0, L, m); if (getenv("OPEN")) { pr.extra.clear(); }
        sa.setup(pr);
        for (int rs = 0; rs < restarts; rs++) {
            int c = sa.run(iters, 2.0, 0.05);
            if (c == 0) { sol = sa.bestCfgs(); return L; }
        }
    }
    return -1;
}

int main(int argc, char** argv) {
    auto pats = readInput(argv[1]);
    int Lmin = atoi(argv[2]), Lmax = atoi(argv[3]);
    long long iters = atoll(argv[4]);
    int restarts = atoi(argv[5]);
    int m = atoi(argv[6]);
    vector<int> idx; if (argc > 7) { ifstream f(argv[7]); int x; while (f >> x) idx.push_back(x); } else for (int i = 0; i < (int)pats.size(); i++) idx.push_back(i);
    int n = pats.size();
    vector<string> res(n);
    auto t0 = chrono::steady_clock::now();
#pragma omp parallel for schedule(dynamic, 1)
    for (int ii = 0; ii < (int)idx.size(); ii++) {
        int i = idx[ii];
        vector<Cfg> sol;
        int L;
        auto t1 = chrono::steady_clock::now();
        if (pats[i].N <= 64) L = minCycle<1>(pats[i], Lmin, Lmax, iters, restarts, i + 1, sol, m);
        else if (pats[i].N <= 128) L = minCycle<2>(pats[i], Lmin, Lmax, iters, restarts, i + 1, sol, m);
        else L = minCycle<4>(pats[i], Lmin, Lmax, iters, restarts, i + 1, sol, m);
        double dt = chrono::duration<double>(chrono::steady_clock::now() - t1).count();
        char buf[256];
        snprintf(buf, sizeof buf, "%d r=%d P=%d cnt=%d L=%d time=%.1f", i, pats[i].r, pats[i].nP, pats[i].cnt, L, dt);
        res[i] = buf;
    }
    for (int i : idx) puts(res[i].c_str());
    fprintf(stderr, "total %.1fs\n", chrono::duration<double>(chrono::steady_clock::now() - t0).count());
}
