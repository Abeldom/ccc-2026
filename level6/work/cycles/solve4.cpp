// solve4: chain improver. Takes an existing output file; for each pattern with count >= 2, splits the program
// into builds and greedily rebuilds builds 2..cnt (and optionally build 1): for build j, starting after the
// current end config X, search (SA, descending length, warm starts from the original build) for the shortest
// build from X. Writes verified programs that are shorter than the input.
// Usage: solve4 <input> <out file> <out sol> <iters per tick> <restarts> <time per pattern s> <idx file|all> [first]
#include "sa3.h"

static double nowS() { return chrono::duration<double>(chrono::steady_clock::now().time_since_epoch()).count(); }

static bool progToCfgs(const vector<vector<string>>& prog, int r, vector<Cfg>& seq) {
    Cfg cur = initCfg();
    seq.clear();
    for (auto& line : prog) {
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
static vector<vector<vector<string>>> readOut(const string& path) {
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

template <int W>
struct Chain {
    const Pattern& p;
    SA<W> sa;
    Sim sim;
    long long ipt;
    int restarts;
    double t0, tlim;
    Chain(const Pattern& pat, u64 seed, long long ipt_, int rs, double tl) : p(pat), sa(pat, seed), sim(pat), ipt(ipt_), restarts(rs), tlim(tl) { t0 = nowS(); }
    double deadline = 1e18;
    bool timeUp() const { return nowS() > deadline; }
    static int tdist(ES x, ES y) { if (x == y) return 0; if (x.k == y.k) return 1; return abs(x.k - y.k) + x.a + y.a; }
    // X -> ... -> C with inactive transit; returns the configs strictly between X and C
    static vector<Cfg> transit(const Cfg& X, const Cfg& C) {
        int T = 0;
        for (int e = 0; e < 6; e++) T = max(T, tdist(X[e], C[e]));
        vector<Cfg> v(max(0, T - 1));
        for (int e = 0; e < 6; e++) {
            ES x = X[e], y = C[e];
            vector<ES> path;
            if (x == y) { for (int t = 0; t < T; t++) path.push_back(x); }
            else if (x.k == y.k) { for (int t = 0; t < T - 1; t++) path.push_back(x); path.push_back(y); }
            else {
                ES cur = x;
                if (cur.a) { cur.a = 0; path.push_back(cur); }
                while (cur.k != y.k) { cur.k += (y.k > cur.k) ? 1 : -1; path.push_back(cur); }
                while ((int)path.size() < T - (y.a ? 1 : 0)) path.push_back(cur);
                if (y.a) path.push_back(y);
            }
            for (int t = 0; t + 1 < T; t++) v[t][e] = path[t];
        }
        return v;
    }
    // exact check: build from canonical accepts exactly at its last tick (not earlier)
    int buildAccept(const vector<Cfg>& b) {
        auto acc = sim.run(b);
        return acc.empty() ? -1 : acc[0];
    }
    // shortest build from X; `orig` (may be empty) is a known build (whose first config may not be adjacent to X)
    vector<Cfg> shortest(const Cfg& X, const vector<Cfg>& orig, int maxLen) {
        vector<Cfg> best;
        if (!orig.empty() && adjCfg(X, orig[0]) && buildAccept(orig) == (int)orig.size()) best = orig;
        if (best.empty() && !orig.empty()) {
            vector<Cfg> v = transit(X, orig[0]);
            v.insert(v.end(), orig.begin(), orig.end());
            int a = buildAccept(v);
            if (a > 0) { v.resize(a); best = v; }
        }
        int cur = best.empty() ? maxLen + 1 : (int)best.size();
        for (int l = cur - 1; l >= 1; l--) {
            if (timeUp()) break;
            Problem pr = makeChain({l}, X);
            sa.setup(pr);
            bool ok = false;
            for (int q = 0; q < restarts && !ok && !timeUp(); q++) {
                int c;
                if (q % 2 == 0 && !best.empty()) {
                    vector<Cfg> v = best;
                    while ((int)v.size() > l) v.erase(v.begin() + sa.rng.rint(v.size()));
                    sa.loadCfgs(v);
                    if (!sa.repair()) { if (!sa.initRandom()) continue; }
                    c = sa.anneal(ipt * l, 0.8, 0.05);
                } else c = sa.run(ipt * l, 2.0, 0.05);
                if (c == 0) {
                    auto v = sa.bestCfgs();
                    int a = buildAccept(v);
                    if (a > 0) { v.resize(a); best = v; ok = true; l = a; }
                }
            }
            if (!ok) break;
        }
        return best;
    }
};

int main(int argc, char** argv) {
    auto pats = readInput(argv[1]);
    int n = pats.size();
    auto progs = readOut(argv[2]);
    string outPath = argv[3];
    long long ipt = atoll(argv[4]);
    int restarts = atoi(argv[5]);
    double tlim = atof(argv[6]);
    vector<int> idx;
    if (string(argv[7]) != "all") { ifstream f(argv[7]); int x; while (f >> x) idx.push_back(x); }
    else for (int i = 0; i < n; i++) idx.push_back(i);
    bool doFirst = argc > 8;
    FILE* fo = fopen(outPath.c_str(), "a");
#pragma omp parallel for schedule(dynamic, 1)
    for (int ii = 0; ii < (int)idx.size(); ii++) {
        int i = idx[ii];
        const Pattern& p = pats[i];
        if (p.cnt < 2 && !doFirst) continue;
        vector<Cfg> seq;
        if (!progToCfgs(progs[i], p.r, seq)) continue;
        double t1 = nowS();
        vector<Cfg> out;
        string log;
        auto go = [&](auto* tag) {
            constexpr int W = std::remove_pointer_t<decltype(tag)>::value;
            Chain<W> ch(p, 98765ULL + i * 131ULL, ipt, restarts, tlim);
            auto acc = ch.sim.run(seq);
            if ((int)acc.size() < p.cnt) { log = "invalid input"; return; }
            vector<int> st(p.cnt), en(p.cnt);
            for (int j = 0; j < p.cnt; j++) { st[j] = j ? acc[j - 1] : 0; en[j] = acc[j]; }
            vector<Cfg> res;
            for (int j = 0; j < p.cnt; j++) {
                vector<Cfg> orig(seq.begin() + st[j], seq.begin() + en[j]);
                if (j == 0 && !doFirst) { res = orig; continue; }
                Cfg X = res.empty() ? initCfg() : res.back();
                int remB = p.cnt - j;
                ch.deadline = nowS() + max(0.0, (t1 + tlim - nowS())) / remB;
                // the remaining time is shared among remaining builds
                vector<Cfg> b = ch.shortest(X, orig, (int)orig.size() + 4);
                if (b.empty()) {
                    // fall back: transit to orig[0] by simple moves is not guaranteed safe; give up on this pattern
                    log += " fail";
                    res.clear();
                    break;
                }
                char buf[32];
                snprintf(buf, sizeof buf, " %d>%zu", (int)orig.size(), b.size());
                log += buf;
                res.insert(res.end(), b.begin(), b.end());
            }
            if (!res.empty()) {
                auto a2 = ch.sim.run(res);
                bool adjok = true;
                Cfg prev = initCfg();
                for (auto& c : res) { if (!adjCfg(prev, c)) adjok = false; prev = c; }
                if (adjok && (int)a2.size() >= p.cnt) { res.resize(a2[p.cnt - 1]); out = res; }
            }
        };
        if (p.N <= 64) go((integral_constant<int, 1>*)nullptr);
        else if (p.N <= 128) go((integral_constant<int, 2>*)nullptr);
        else if (p.N <= 192) go((integral_constant<int, 3>*)nullptr);
        else if (p.N <= 256) go((integral_constant<int, 4>*)nullptr);
        else go((integral_constant<int, 6>*)nullptr);
#pragma omp critical
        {
            fprintf(stderr, "%d r=%d P=%d cnt=%d len=%zu -> %zu%s%s (%.0fs)\n", i, p.r, p.nP, p.cnt, seq.size(), out.size(), log.c_str(),
                    (!out.empty() && out.size() < seq.size()) ? " WIN" : "", nowS() - t1);
            if (!out.empty() && out.size() < seq.size()) {
                auto prog = toProgram(out);
                fprintf(fo, "# %d total=%zu\n", i, prog.size());
                for (auto& l : prog) fprintf(fo, "%s\n", l.c_str());
                fflush(fo);
            }
        }
    }
    fclose(fo);
}
