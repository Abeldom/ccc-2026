
template<int NW>
vector<string> solveT(const Pat& P, const Geo& G, const Params& prm, int modeMask, string& info) {
    int r = P.r;
    vector<pair<int,int>> Pc, Kc; set<pair<int,int>> Ps;
    for (int y = 0; y < (int)P.g.size(); y++) for (int x = 0; x < (int)P.g[y].size(); x++) if (P.g[y][x] == 'M') { Pc.push_back({y, x}); Ps.insert({y, x}); }
    int ny[6] = {0, 0, 1, 1, -1, -1}, nx[6] = {2, -2, 1, -1, 1, -1}; set<pair<int,int>> Ks;
    for (auto [y, x] : Pc) for (int d = 0; d < 6; d++) { pair<int,int> c = {y + ny[d], x + nx[d]}; if (G.interior.count(c) && !Ps.count(c)) Ks.insert(c); }
    Kc.assign(Ks.begin(), Ks.end());
    Solver<NW> S; S.init(r, Pc, Kc, G);
    vector<string> best; int bestLen = INT_MAX; string bestInfo;
    // mode 1: cycle
    if (modeMask & 2) {
        auto res = S.beam(S.I, 1, P.cnt, prm);
        if (res.total < INT_MAX) {
            vector<Cfg> seq = S.concretize(res, S.I);
            vector<string> prog; Cfg cur = S.I;
            S.emit(cur, seq, prog, cur);
            for (int b = 1; b < P.cnt; b++) S.emit(cur, seq, prog, cur);
            int t = simulate(G, P, prog, P.cnt);
            if (t > 0) { prog.resize(t); if ((int)prog.size() < bestLen) { bestLen = prog.size(); best = prog; bestInfo = "cycle L=" + to_string(seq.size()) + " pred=" + to_string(res.total); } }
            else bestInfo = "cyclefail";
        }
    }
    if ((modeMask & 1) && !(P.cnt == 1 && (modeMask & 2))) {
        vector<string> prog; Cfg cur = S.I; bool ok = true; int acc = 0;
        while (acc < P.cnt) {
            Params p2 = prm; p2.seed = prm.seed + acc * 1000 + 17;
            auto res = S.beam(cur, 0, 1, p2);
            if (res.total == INT_MAX) { ok = false; break; }
            vector<Cfg> seq = S.concretize(res, cur);
            size_t before = prog.size();
            Cfg endc; S.emit(cur, seq, prog, endc);
            int a2 = 0; int t = simulate(G, P, prog, acc + 1, &a2);
            if (t < 0) { ok = false; bestInfo += " chainfail"; break; }
            if (t < (int)before + 1) { ok = false; break; }
            prog.resize(t); acc++;
            cur = cfgAfter(prog, r, prog.size());
            if ((int)prog.size() >= bestLen) { ok = false; break; }
        }
        if (ok && (int)prog.size() < bestLen) { bestLen = prog.size(); best = prog; bestInfo = "chain"; }
    }
    info = bestInfo;
    return best;
}

int main(int argc, char** argv) {
    if (argc < 3) { fprintf(stderr, "usage: shots in out [beam finW rootN modeMask seed maxAct closW only]\n"); return 1; }
    Params prm; int modeMask = 3; int only = -1;
    if (argc > 3) prm.beam = atoi(argv[3]);
    if (argc > 4) prm.finW = atoi(argv[4]);
    if (argc > 5) prm.rootN = atoi(argv[5]);
    if (argc > 6) modeMask = atoi(argv[6]);
    if (argc > 7) prm.seed = atoll(argv[7]);
    if (argc > 8) prm.maxAct = atoi(argv[8]);
    if (argc > 9) prm.closW = atof(argv[9]);
    if (argc > 10) only = atoi(argv[10]);
    auto pats = parseInput(argv[1]);
    int N = pats.size();
    map<int, Geo> geos; for (auto& p : pats) if (!geos.count(p.r)) geos[p.r] = makeGeo(p.r);
    vector<vector<string>> out(N); vector<string> infos(N);
    vector<double> tm(N);
    #pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < N; i++) {
        if (only >= 0 && i != only) continue;
        auto t0 = chrono::steady_clock::now();
        const Pat& P = pats[i]; const Geo& G = geos.at(P.r);
        int ncell = 3 * P.r * (P.r - 1) + 1;
        int nw = (ncell + 63) / 64;
        vector<string> prog; string info;
        switch (nw) {
            case 1: prog = solveT<1>(P, G, prm, modeMask, info); break;
            case 2: prog = solveT<2>(P, G, prm, modeMask, info); break;
            case 3: prog = solveT<3>(P, G, prm, modeMask, info); break;
            case 4: prog = solveT<4>(P, G, prm, modeMask, info); break;
            case 5: prog = solveT<5>(P, G, prm, modeMask, info); break;
            case 6: prog = solveT<6>(P, G, prm, modeMask, info); break;
            case 7: prog = solveT<7>(P, G, prm, modeMask, info); break;
            default: prog = solveT<8>(P, G, prm, modeMask, info); break;
        }
        out[i] = prog; infos[i] = info;
        tm[i] = chrono::duration<double>(chrono::steady_clock::now() - t0).count();
        #pragma omp critical
        { fprintf(stderr, "case %d r=%d cnt=%d -> %d lines (%s) %.1fs\n", i, P.r, P.cnt, (int)prog.size(), info.c_str(), tm[i]); }
    }
    ofstream f(argv[2], ios::binary);
    long long tot = 0;
    for (int i = 0; i < N; i++) {
        if (out[i].empty()) out[i].push_back("- - - - - -");
        for (auto& s : out[i]) { f << s[0] << ' ' << s[1] << ' ' << s[2] << ' ' << s[3] << ' ' << s[4] << ' ' << s[5] << '\n'; }
        tot += out[i].size();
        if (i + 1 < N) f << '\n';
    }
    fprintf(stderr, "total %lld\n", tot);
}
