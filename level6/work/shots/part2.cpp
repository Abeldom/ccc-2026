    // closure ignoring emitters in freeMask (they only need to be off, always satisfiable directly)
    int closureF(const Cfg& X, const Cfg& C, uint8_t freeMask) const {
        bool dir = true; int m1 = 0, m2 = 0, i1 = -1;
        for (int e = 0; e < 6; e++) { if ((freeMask >> e) & 1) continue; int s = sidx(X, e), t = sidx(C, e);
            if (!directT[s * NS + t]) dir = false;
            int v = noffT[s * NS + t]; if (v > m1) { m2 = m1; m1 = v; i1 = e; } else if (v > m2) m2 = v; }
        if (dir) return 0;
        int best = m1;
        if (i1 >= 0) { int alt = max(nliveT[sidx(X, i1) * NS + sidx(C, i1)], m2); if (alt < best) best = alt; }
        return max(best, 1);
    }

    struct Node { Cfg c; uint8_t flo; int8_t yk[6]; BS U; int root; int par; int nP, nK; };
    struct EOpt { int8_t k; uint8_t a, flo, act, idx; };
    struct AxOpt { int8_t ki, kj; uint8_t ai, aj, fi, fj, act, ii, ij; BS m1, m2; };

    void eopts(const Node& nd, int e, vector<EOpt>& v) const {
        v.clear();
        if ((nd.flo >> e) & 1) {
            v.push_back({nd.c.k[e], 0, 1, 0, 0});
            for (int k = -(r - 1); k <= r - 1; k++) v.push_back({(int8_t)k, 1, 0, 1, (uint8_t)(k + r)});
        } else {
            int k = nd.c.k[e], a = (nd.c.a >> e) & 1;
            v.push_back({(int8_t)k, (uint8_t)a, 0, 0, 0});
            if (k + 1 <= r) v.push_back({(int8_t)(k + 1), (uint8_t)a, 0, 0, 1});
            if (k - 1 >= -r) v.push_back({(int8_t)(k - 1), (uint8_t)a, 0, 0, 2});
            v.push_back({(int8_t)k, (uint8_t)(a ^ 1), 0, 0, 3});
        }
    }
    // apply option index idx to emitter e of node -> (k,a,flo)
    EOpt eoptByIdx(const Node& nd, int e, int idx) const {
        if ((nd.flo >> e) & 1) { if (idx == 0) return {nd.c.k[e], 0, 1, 0, 0}; return {(int8_t)(idx - r), 1, 0, 1, (uint8_t)idx}; }
        int k = nd.c.k[e], a = (nd.c.a >> e) & 1;
        if (idx == 1) k++; else if (idx == 2) k--; else if (idx == 3) a ^= 1;
        return {(int8_t)k, (uint8_t)a, 0, 0, (uint8_t)idx};
    }

    Cfg yOf(const Node& nd, const vector<Cfg>& roots) const {
        Cfg y = roots[nd.root];
        for (int e = 0; e < 6; e++) if (!((y.a >> e) & 1) && !((nd.flo >> e) & 1)) y.k[e] = nd.yk[e];
        return y;
    }

    struct Result { int total = INT_MAX; vector<Node> path; vector<Cfg> roots; };

    Node child(const Node& nd, int combo, int pi) const {
        Node ch = nd; ch.par = pi;
        for (int e = 0; e < 6; e++) { int idx = (combo >> (5 * e)) & 31; EOpt o = eoptByIdx(nd, e, idx);
            ch.c.k[e] = o.k; if (o.a) ch.c.a |= 1 << e; else ch.c.a &= ~(1 << e);
            if (o.flo) ch.flo |= 1 << e; else ch.flo &= ~(1 << e);
            if (o.act) ch.yk[e] = o.k; }
        BS S, D; effect(ch.c, S, D);
        for (int i = 0; i < NW; i++) ch.U[i] = nd.U[i] & ~(S[i] | D[i]);
        BS t1, t2; for (int i = 0; i < NW; i++) { t1[i] = ch.U[i] & Pm[i]; t2[i] = ch.U[i] & Km[i]; }
        ch.nP = pc(t1); ch.nK = pc(t2);
        return ch;
    }
    u64 hashNode(const Node& n) const {
        u64 h = 0x12345;
        for (int a = 0; a < 3; a++) {
            int l1 = ((n.flo >> a) & 1) ? 999 : n.c.k[a] * SGN[a] * 2 + ((n.c.a >> a) & 1);
            int l2 = ((n.flo >> (a + 3)) & 1) ? 999 : n.c.k[a + 3] * SGN[a + 3] * 2 + ((n.c.a >> (a + 3)) & 1);
            if (l1 > l2) swap(l1, l2);
            h = mix64(h ^ (u64)(l1 + 100) ^ ((u64)(l2 + 100) << 16));
        }
        for (int i = 0; i < NW; i++) h = mix64(h ^ n.U[i]);
        return h;
    }

    // mode 0: chain from X: objective = closure(X,C1) + L
    // mode 1: cycle: objective = closure(I,C1) + cnt*L + (cnt-1)*closure(Y,C1)
    Result beam(const Cfg& X, int mode, int cnt, const Params& prm) {
        int wL = (mode == 0) ? 1 : cnt, wY = (mode == 0) ? 0 : cnt - 1;
        const Cfg& X0 = (mode == 0) ? X : I;
        Result res;
        BS All{}; for (int i = 0; i < NW; i++) All[i] = Pm[i] | Km[i];
        // ---- roots
        vector<pair<int,int>> opts; opts.push_back({99, 99});
        for (int a = -(r - 1); a <= r - 1; a++) { opts.push_back({a, 99}); opts.push_back({a, a}); for (int b = a + 1; b <= r - 1; b++) opts.push_back({a, b}); }
        int no = opts.size();
        vector<array<BS,2>> om(3 * no);
        for (int ax = 0; ax < 3; ax++) for (int o = 0; o < no; o++) {
            BS x1{}, x2{}; auto [l1, l2] = opts[o];
            BS a{}, b{}; a.fill(0); b.fill(0);
            if (l1 != 99) a = mask[ax * R2 + l1 * SGN[ax] + r];
            if (l2 != 99) b = mask[(ax + 3) * R2 + l2 * SGN[ax + 3] + r];
            for (int i = 0; i < NW; i++) { x1[i] = a[i] ^ b[i]; x2[i] = a[i] & b[i]; }
            om[ax * no + o] = {x1, x2};
        }
        typedef tuple<long long,int,int,int> RT;
        priority_queue<RT, vector<RT>, greater<RT>> rh;
        for (int oa = 0; oa < no; oa++) for (int ob = 0; ob < no; ob++) {
            const BS &x1 = om[oa][0], &x2 = om[oa][1], &y1 = om[no + ob][0], &y2 = om[no + ob][1];
            BS ab0, ab1, ab2, ab3;
            for (int i = 0; i < NW; i++) { u64 xa = x1[i] | x2[i], ya = y1[i] | y2[i];
                ab1[i] = (x1[i] & ~ya) | (y1[i] & ~xa); ab2[i] = (x2[i] & ~ya) | (y2[i] & ~xa) | (x1[i] & y1[i]);
                ab3[i] = (xa | ya) & ~ab1[i] & ~ab2[i]; ab0[i] = ~(xa | ya); }
            for (int oc = 0; oc < no; oc++) {
                const BS &z1 = om[2 * no + oc][0], &z2 = om[2 * no + oc][1];
                bool ok = true; int sc = 0;
                for (int i = 0; i < NW && ok; i++) { u64 za = z1[i] | z2[i];
                    u64 S = (ab2[i] & ~za) | (ab1[i] & z1[i]) | (ab0[i] & z2[i]);
                    u64 D = ab3[i] | (ab2[i] & za) | (ab1[i] & z2[i]);
                    if ((S & Km[i]) | (D & Pm[i])) ok = false; else sc += __builtin_popcountll(S & Pm[i]) * 64 + __builtin_popcountll(D & Km[i]); }
                if (!ok || sc < 64) continue;
                long long key = (long long)sc * 1000003 + (long long)(mix64(prm.seed * 7 + oa * 1000000 + ob * 1000 + oc) % 1000003);
                if ((int)rh.size() < prm.rootN) rh.push({key, oa, ob, oc});
                else if (key > get<0>(rh.top())) { rh.pop(); rh.push({key, oa, ob, oc}); }
            }
        }
        vector<Cfg>& roots = res.roots;
        vector<vector<Node>> layers(1);
        layers.reserve(prm.maxDepth + 4);
        while (!rh.empty()) { auto [key, oa, ob, oc] = rh.top(); rh.pop();
            Cfg c; c.a = 0; int oo[3] = {oa, ob, oc};
            for (int ax = 0; ax < 3; ax++) { auto [l1, l2] = opts[oo[ax]];
                c.k[ax] = (l1 != 99) ? l1 * SGN[ax] : 0; if (l1 != 99) c.a |= 1 << ax;
                c.k[ax + 3] = (l2 != 99) ? l2 * SGN[ax + 3] : 0; if (l2 != 99) c.a |= 1 << (ax + 3); }
            BS S, D; effect(c, S, D);
            Node nd; nd.c = c; nd.root = roots.size(); nd.par = -1; nd.flo = (~c.a) & 63;
            for (int e = 0; e < 6; e++) nd.yk[e] = c.k[e];
            for (int i = 0; i < NW; i++) nd.U[i] = All[i] & ~(S[i] | D[i]);
            BS t1, t2; for (int i = 0; i < NW; i++) { t1[i] = nd.U[i] & Pm[i]; t2[i] = nd.U[i] & Km[i]; }
            nd.nP = pc(t1); nd.nK = pc(t2);
            roots.push_back(c); layers[0].push_back(nd); }
        if (roots.empty()) return res;
        int bestDepth = -1, bestPar = -1, bestCombo = -1;
        auto clYof = [&](const Node& nd) { return closureF(yOf(nd, roots), nd.c, nd.flo); };
        auto evalFinish = [&](const Node& nd) { return closureF(X0, nd.c, nd.flo) + wY * ((wY > 0) ? clYof(nd) : 0); };
        for (int i = 0; i < (int)layers[0].size(); i++) if (layers[0][i].nP == 0) {
            int tot = wL * 1 + evalFinish(layers[0][i]);
            if (tot < res.total) { res.total = tot; bestDepth = 0; bestPar = i; bestCombo = -1; } }
        vector<AxOpt> axo[3]; vector<EOpt> e1, e2;
        int bestNPprev = np;
        for (int depth = 1; depth < prm.maxDepth; depth++) {
            if (wL * (depth + 1) >= res.total) break;
            vector<Node>& cur = layers[depth - 1];
            if (cur.empty()) break;
            int mnp = INT_MAX; for (auto& nd : cur) mnp = min(mnp, nd.nP);
            bestNPprev = mnp;
            double rate = max(1.0, (double)(np - bestNPprev) / depth);
            int rate16 = (int)(rate * 16 * prm.closW);
            typedef tuple<u64,int,int> CT;
            priority_queue<CT> h1, h2;
            size_t cap1 = prm.beam * 2, cap2 = prm.finW * 2;
            for (int pi = 0; pi < (int)cur.size(); pi++) {
                const Node& nd = cur[pi];
                BS UP, UK; for (int i = 0; i < NW; i++) { UP[i] = nd.U[i] & Pm[i]; UK[i] = nd.U[i] & Km[i]; }
                Cfg Ynd; if (mode == 1) Ynd = roots[nd.root];
                for (int a = 0; a < 3; a++) {
                    axo[a].clear(); eopts(nd, a, e1); eopts(nd, a + 3, e2);
                    for (auto& oi : e1) for (auto& oj : e2) {
                        if (oi.act + oj.act > prm.maxAct) continue;
                        AxOpt o; o.ki = oi.k; o.kj = oj.k; o.ai = oi.a; o.aj = oj.a; o.fi = oi.flo; o.fj = oj.flo; o.act = oi.act + oj.act; o.ii = oi.idx; o.ij = oj.idx;
                        BS A{}, B{}; A.fill(0); B.fill(0);
                        if (oi.a) A = mask[a * R2 + oi.k + r]; if (oj.a) B = mask[(a + 3) * R2 + oj.k + r];
                        for (int i = 0; i < NW; i++) { o.m1[i] = A[i] ^ B[i]; o.m2[i] = A[i] & B[i]; }
                        axo[a].push_back(o);
                    }
                }
                for (auto& X1 : axo[0]) {
                for (auto& Y1 : axo[1]) { if (X1.act + Y1.act > prm.maxAct) continue;
                    BS ab0, ab1, ab2, ab3;
                    for (int i = 0; i < NW; i++) { u64 xa = X1.m1[i] | X1.m2[i], ya = Y1.m1[i] | Y1.m2[i];
                        ab1[i] = (X1.m1[i] & ~ya) | (Y1.m1[i] & ~xa); ab2[i] = (X1.m2[i] & ~ya) | (Y1.m2[i] & ~xa) | (X1.m1[i] & Y1.m1[i]);
                        ab3[i] = (xa | ya) & ~ab1[i] & ~ab2[i]; ab0[i] = ~(xa | ya); }
                    for (auto& Z1 : axo[2]) { if (X1.act + Y1.act + Z1.act > prm.maxAct) continue;
                        bool ok = true; int nP = 0, nK = 0;
                        for (int i = 0; i < NW; i++) { u64 za = Z1.m1[i] | Z1.m2[i];
                            u64 S = (ab2[i] & ~za) | (ab1[i] & Z1.m1[i]) | (ab0[i] & Z1.m2[i]);
                            u64 D = ab3[i] | (ab2[i] & za) | (ab1[i] & Z1.m2[i]);
                            if ((S & UK[i]) | (D & UP[i])) { ok = false; break; }
                            nP += __builtin_popcountll(UP[i] & ~S); nK += __builtin_popcountll(UK[i] & ~D); }
                        if (!ok) continue;
                        int combo = X1.ii | (Y1.ii << 5) | (Z1.ii << 10) | (X1.ij << 15) | (Y1.ij << 20) | (Z1.ij << 25);
                        // child light-weight node for closure
                        Node cl; cl.c.k[0] = X1.ki; cl.c.k[3] = X1.kj; cl.c.k[1] = Y1.ki; cl.c.k[4] = Y1.kj; cl.c.k[2] = Z1.ki; cl.c.k[5] = Z1.kj;
                        cl.c.a = X1.ai | (Y1.ai << 1) | (Z1.ai << 2) | (X1.aj << 3) | (Y1.aj << 4) | (Z1.aj << 5);
                        cl.flo = X1.fi | (Y1.fi << 1) | (Z1.fi << 2) | (X1.fj << 3) | (Y1.fj << 4) | (Z1.fj << 5);
                        cl.root = nd.root;
                        for (int e = 0; e < 6; e++) cl.yk[e] = nd.yk[e];
                        if (X1.act | Y1.act | Z1.act) { for (int e = 0; e < 6; e++) if (((nd.flo >> e) & 1) && !((cl.flo >> e) & 1)) cl.yk[e] = cl.c.k[e]; }
                        u64 noise = mix64(prm.seed ^ ((u64)depth << 40) ^ ((u64)pi << 30) ^ (u64)combo);
                        if (nP > 0) {
                            u64 est = (u64)nP * 16;
                            if (rate16 > 0) {
                                int clo = (mode == 1 && cnt > 1) ? clYof(cl) : closureF(X0, cl.c, cl.flo);
                                est = max(est, (u64)rate16 * clo);
                            }
                            u64 key = (est << 40) | ((u64)min(nK, 65535) << 24) | (noise & 0xFFFFFF);
                            if (h1.size() < cap1) h1.push({key, pi, combo}); else if (key < get<0>(h1.top())) { h1.pop(); h1.push({key, pi, combo}); }
                        } else {
                            int rem = evalFinish(cl);
                            int tot = wL * (depth + 1) + rem;
                            if (tot < res.total) { res.total = tot; bestDepth = depth; bestPar = pi; bestCombo = combo; }
                            u64 key = ((u64)rem << 32) | (noise & 0xFFFFFFFF);
                            if (h2.size() < cap2) h2.push({key, pi, combo}); else if (key < get<0>(h2.top())) { h2.pop(); h2.push({key, pi, combo}); }
                        }
                    } } }
            }
            layers.emplace_back();
            vector<Node>& nxt = layers.back();
            unordered_set<u64> seen;
            auto take = [&](priority_queue<CT>& h, int lim) {
                vector<CT> v; while (!h.empty()) { v.push_back(h.top()); h.pop(); }
                sort(v.begin(), v.end());
                int got = 0;
                for (auto& [key, pi, combo] : v) { if (got >= lim) break;
                    Node ch = child(cur[pi], combo, pi);
                    u64 hh = hashNode(ch); if (!seen.insert(hh).second) continue;
                    nxt.push_back(ch); got++; } };
            take(h2, prm.finW); take(h1, prm.beam);
        }
        if (bestDepth < 0) return res;
        Node last;
        if (bestCombo < 0) last = layers[bestDepth][bestPar];
        else last = child(layers[bestDepth - 1][bestPar], bestCombo, bestPar);
        res.path.push_back(last);
        int p = last.par;
        for (int li = bestDepth - 1; li >= 0; li--) { const Node& nd = layers[li][p]; res.path.push_back(nd); p = nd.par; }
        return res;
    }

    // concrete config sequence (forward order) from a result path, given the build's start config X
    vector<Cfg> concretize(const Result& res, const Cfg& X) const {
        const Node& fin = res.path[0];
        vector<Cfg> seq;
        for (auto& nd : res.path) {
            Cfg c = nd.c;
            for (int e = 0; e < 6; e++) if ((nd.flo >> e) & 1) {
                c.k[e] = ((fin.flo >> e) & 1) ? X.k[e] : fin.yk[e];
                c.a &= ~(1 << e);
            }
            seq.push_back(c);
        }
        return seq;
    }

    // produce program lines for: transit X->seq[0], then seq
    void emit(const Cfg& X, const vector<Cfg>& seq, vector<string>& prog, Cfg& endCfg) const {
        int liveE; int n = closure(X, seq[0], &liveE);
        vector<Cfg> tr = (n > 0) ? transit(X, seq[0], n, liveE) : vector<Cfg>();
        Cfg cur = X;
        for (auto& c : tr) { prog.push_back(stepInstr(cur, c)); cur = c; }
        for (auto& c : seq) { prog.push_back(stepInstr(cur, c)); cur = c; }
        endCfg = cur;
    }
};
