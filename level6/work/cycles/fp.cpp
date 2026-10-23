// Polish programs with full-program SA: repeatedly delete one tick and repair with SA.
// Usage: fp <input> <hints file> <current out file> <out sol file> <iters> <attempts> <timeLimitSec> <idx file|all>
#include "fpsa.h"

static bool progToCfgs(const vector<vector<string>>& prog, int r, vector<Cfg>& seq) {
    Cfg cur;
    for (int e = 0; e < 6; e++) cur[e] = {0, 0};
    seq.clear();
    for (auto& line : prog) {
        for (int e = 0; e < 6; e++) {
            char c = line[e][0];
            if (c == 'F') cur[e].k++;
            else if (c == 'B') cur[e].k--;
            else if (c == 'T') cur[e].a ^= 1;
            if (abs(cur[e].k) > r) return false;
        }
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
vector<string> polish(const Pattern& p, const vector<Cfg>& seq0, int minLen, long long iters, int attempts, double tlim,
                      u64 seed, string& log) {
    Sim sim(p);
    auto acc = sim.run(seq0);
    if ((int)acc.size() < p.cnt) { log = "invalid start"; return {}; }
    vector<Cfg> seq(seq0.begin(), seq0.begin() + acc[p.cnt - 1]);
    vector<int> ends;
    for (int j = 0; j < p.cnt; j++) ends.push_back(acc[j] - 1);
    FPSA<W> fp(p, seed);
    fp.minLen = minLen;
    fp.load(seq, ends);
    if (fp.total != 0) { log = "cost model mismatch"; return {}; }
    vector<Cfg> best = seq;
    int startLen = seq.size();
    auto t0 = chrono::steady_clock::now();
    while (true) {
        if (chrono::duration<double>(chrono::steady_clock::now() - t0).count() > tlim) break;
        bool ok = false;
        vector<Cfg> cur = best;
        vector<int> curEnds = fp.segEnd;
        for (int a = 0; a < attempts && !ok; a++) {
            fp.load(cur, curEnds);
            // choose a deletable tick: weighted by segment length
            vector<int> cand;
            for (int t = 0; t < fp.T; t++) {
                int j = fp.segId[t];
                if (fp.segEnd[j] - fp.segStart(j) + 1 > minLen) cand.push_back(t);
            }
            if (cand.empty()) break;
            int d = cand[fp.rng.rint(cand.size())];
            if (!fp.deleteTick(d)) continue;
            if (fp.run(iters, 0.6, 0.05)) {
                vector<Cfg> c = fp.cfgs();
                auto ac = sim.run(c);
                if ((int)ac.size() >= p.cnt) {
                    c.resize(ac[p.cnt - 1]);
                    best = c;
                    vector<int> e2;
                    for (int j = 0; j < p.cnt; j++) e2.push_back(ac[j] - 1);
                    fp.load(best, e2);
                    ok = true;
                }
            }
        }
        if (!ok) break;
    }
    char buf[128];
    snprintf(buf, sizeof buf, "%d -> %zu", startLen, best.size());
    log = buf;
    return toProgram(best);
}

int main(int argc, char** argv) {
    auto pats = readInput(argv[1]);
    int n = pats.size();
    vector<int> minLen(n, 1);
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
            if (s >= 0) minLen[idx] = 1;
            else if (t > 0) minLen[idx] = 2;
            else if (t == 0) minLen[idx] = 3;
            else minLen[idx] = 2;
        }
    }
    auto progs = readOut(argv[3]);
    string outPath = argv[4];
    long long iters = atoll(argv[5]);
    int attempts = atoi(argv[6]);
    double tlim = atof(argv[7]);
    vector<int> idx;
    if (string(argv[8]) != "all") { ifstream f(argv[8]); int x; while (f >> x) idx.push_back(x); }
    else for (int i = 0; i < n; i++) idx.push_back(i);
    vector<vector<string>> res(n);
#pragma omp parallel for schedule(dynamic, 1)
    for (int ii = 0; ii < (int)idx.size(); ii++) {
        int i = idx[ii];
        const Pattern& p = pats[i];
        vector<Cfg> seq;
        string log = "skip";
        if (progToCfgs(progs[i], p.r, seq)) {
            u64 sd = 777 + i * 31;
            if (p.N <= 64) res[i] = polish<1>(p, seq, minLen[i], iters, attempts, tlim, sd, log);
            else if (p.N <= 128) res[i] = polish<2>(p, seq, minLen[i], iters, attempts, tlim, sd, log);
            else if (p.N <= 192) res[i] = polish<3>(p, seq, minLen[i], iters, attempts, tlim, sd, log);
            else if (p.N <= 256) res[i] = polish<4>(p, seq, minLen[i], iters, attempts, tlim, sd, log);
            else if (p.N <= 384) res[i] = polish<6>(p, seq, minLen[i], iters, attempts, tlim, sd, log);
            else res[i] = polish<8>(p, seq, minLen[i], iters, attempts, tlim, sd, log);
        }
#pragma omp critical
        fprintf(stderr, "%d r=%d P=%d cnt=%d: %s\n", i, p.r, p.nP, p.cnt, log.c_str());
    }
    FILE* f = fopen(outPath.c_str(), "w");
    for (int i : idx) {
        if (res[i].empty()) continue;
        fprintf(f, "# %d fp total=%zu\n", i, res[i].size());
        for (auto& s : res[i]) fprintf(f, "%s\n", s.c_str());
    }
    fclose(f);
}
