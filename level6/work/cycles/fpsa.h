// Full-program SA: config sequence s[0..T-1] from the initial config, cnt segments with movable
// boundaries, each segment must build P from canonical and accept at its end tick.
#pragma once
#include "sa.h"

template <int W>
struct FPSA {
    const Pattern& p;
    int r, NL;
    EGraph G;
    BS<W> Pm, Km, ALL;
    vector<BS<W>> lm[3];
    int T = 0, nseg = 0, minLen = 1;
    vector<int> segEnd, segId;
    vector<array<int, 6>> s;
    vector<BS<W>> eq2, ge3;
    vector<int> segC;
    int total = 0;
    int initId;
    RNG rng;

    FPSA(const Pattern& pat, u64 seed) : p(pat), r(pat.r), NL(2 * pat.r - 1), G(pat.r), rng(seed) {
        Pm = Km = ALL = BS<W>::zero();
        for (int i = 0; i < p.N; i++) { ALL.set(i); if (i < p.nP) Pm.set(i); else Km.set(i); }
        for (int a = 0; a < 3; a++) {
            lm[a].assign(NL, BS<W>::zero());
            for (int i = 0; i < p.N; i++) lm[a][p.lineOf[i][a]].set(i);
        }
        initId = G.id(ES{0, 0});
    }
    inline int lineOf(int e, int sid) const {
        int k = sid / 2 - r, a = sid & 1;
        if (!a || abs(k) >= r) return -1;
        return SGN[e] * k + r - 1;
    }
    void tickMasks(int t) {
        BS<W> one[3], two[3];
        for (int a = 0; a < 3; a++) {
            int l1 = lineOf(a, s[t][a]), l2 = lineOf(a + 3, s[t][a + 3]);
            one[a] = two[a] = BS<W>::zero();
            if (l1 >= 0 && l1 == l2) two[a] = lm[a][l1];
            else {
                if (l1 >= 0) one[a] |= lm[a][l1];
                if (l2 >= 0) one[a] |= lm[a][l2];
            }
        }
        BS<W> aa = one[0] | two[0], bb = one[1] | two[1], cc = one[2] | two[2];
        eq2[t] = (two[0] & ~bb & ~cc) | (two[1] & ~aa & ~cc) | (two[2] & ~aa & ~bb) | (one[0] & one[1] & ~cc) |
                 (one[0] & one[2] & ~bb) | (one[1] & one[2] & ~aa);
        ge3[t] = (two[0] & (bb | cc)) | (two[1] & (aa | cc)) | (two[2] & (aa | bb)) | (one[0] & one[1] & one[2]);
    }
    int segStart(int j) const { return j == 0 ? 0 : segEnd[j - 1] + 1; }
    int segCost(int j) const {
        int a = segStart(j), b = segEnd[j];
        BS<W> und = ALL, bad = BS<W>::zero();
        for (int t = b; t >= a; t--) {
            BS<W> x = und & eq2[t];
            bad |= x & Km;
            und &= ~x;
            BS<W> y = und & ge3[t];
            bad |= y & Pm;
            und &= ~y;
        }
        bad |= und & Pm;
        return bad.pc();
    }
    void rebuildSegId() {
        segId.assign(T, 0);
        for (int j = 0; j < nseg; j++) for (int t = segStart(j); t <= segEnd[j]; t++) segId[t] = j;
    }
    void fullEval() {
        eq2.assign(T, BS<W>::zero());
        ge3.assign(T, BS<W>::zero());
        for (int t = 0; t < T; t++) tickMasks(t);
        rebuildSegId();
        segC.assign(nseg, 0);
        total = 0;
        for (int j = 0; j < nseg; j++) { segC[j] = segCost(j); total += segC[j]; }
    }
    void load(const vector<Cfg>& seq, const vector<int>& ends) {
        T = seq.size();
        s.assign(T, {});
        for (int t = 0; t < T; t++) for (int e = 0; e < 6; e++) s[t][e] = G.id(seq[t][e]);
        segEnd = ends;
        nseg = ends.size();
        fullEval();
    }
    vector<Cfg> cfgs() const {
        vector<Cfg> v(T);
        for (int t = 0; t < T; t++) for (int e = 0; e < 6; e++) v[t][e] = G.st(s[t][e]);
        return v;
    }
    bool okState(int e, int t, int sid, int lo, int hi) const {
        if (t == 0 && !G.adj[initId][sid]) return false;
        if (t - 1 >= 0 && (t - 1 < lo) && !G.adj[s[t - 1][e]][sid]) return false;
        if (t + 1 < T && (t + 1 > hi) && !G.adj[s[t + 1][e]][sid]) return false;
        return true;
    }
    bool resample(int e, int lo, int hi, vector<int>& out) {
        int len = hi - lo + 1, S = G.S;
        static thread_local vector<double> cnt;
        static thread_local vector<char> ok;
        cnt.assign((size_t)len * S, 0.0);
        ok.assign((size_t)len * S, 0);
        for (int i = 0; i < len; i++)
            for (int x = 0; x < S; x++) ok[i * S + x] = okState(e, lo + i, x, lo, hi);
        for (int x = 0; x < S; x++) cnt[(len - 1) * S + x] = ok[(len - 1) * S + x] ? 1.0 : 0.0;
        for (int i = len - 2; i >= 0; i--)
            for (int x = 0; x < S; x++) {
                if (!ok[i * S + x]) continue;
                double c = 0;
                for (int y : G.nb[x]) c += cnt[(i + 1) * S + y];
                cnt[i * S + x] = c;
            }
        double tot = 0;
        for (int x = 0; x < S; x++) tot += cnt[x];
        if (tot <= 0) return false;
        out.resize(len);
        double z = rng.rd() * tot;
        int cur = -1;
        for (int x = 0; x < S; x++) { if (cnt[x] <= 0) continue; z -= cnt[x]; cur = x; if (z < 0) break; }
        out[0] = cur;
        for (int i = 1; i < len; i++) {
            double c = 0;
            for (int y : G.nb[cur]) c += cnt[i * S + y];
            double z2 = rng.rd() * c;
            int nxt = -1;
            for (int y : G.nb[cur]) { if (cnt[i * S + y] <= 0) continue; z2 -= cnt[i * S + y]; nxt = y; if (z2 < 0) break; }
            cur = nxt;
            out[i] = cur;
        }
        return true;
    }
    bool adjOK() const {
        for (int e = 0; e < 6; e++) {
            if (!G.adj[initId][s[0][e]]) return false;
            for (int t = 0; t + 1 < T; t++) if (!G.adj[s[t][e]][s[t + 1][e]]) return false;
        }
        return true;
    }
    // delete tick d, shrinking its segment; repair adjacency. returns false if impossible
    bool deleteTick(int d) {
        int j = segId[d];
        if (segEnd[j] - segStart(j) + 1 <= minLen) return false;
        s.erase(s.begin() + d);
        T--;
        for (int q = j; q < nseg; q++) segEnd[q]--;
        vector<int> out;
        for (int e = 0; e < 6; e++) {
            bool bad = false;
            if (d == 0) bad = !G.adj[initId][s[0][e]];
            else if (d < T) bad = !G.adj[s[d - 1][e]][s[d][e]];
            if (!bad) continue;
            bool fixed = false;
            for (int w = 1; w <= 4 && !fixed; w++) {
                int lo = max(0, d - w), hi = min(T - 1, d + w - 1);
                if (resample(e, lo, hi, out)) {
                    for (int t = lo; t <= hi; t++) s[t][e] = out[t - lo];
                    fixed = true;
                }
            }
            if (!fixed) return false;
        }
        fullEval();
        return true;
    }
    // SA; returns true if cost 0 reached (state holds the solution)
    bool run(long long iters, double T0, double Tend) {
        if (total == 0) return true;
        vector<int> out;
        vector<array<int, 6>> saved;
        double lt0 = log(T0), lt1 = log(Tend);
        int bestTot = total;
        vector<array<int, 6>> bestS = s;
        vector<int> bestEnd = segEnd;
        for (long long it = 0; it < iters; it++) {
            double temp = exp(lt0 + (lt1 - lt0) * it / iters);
            int mv = rng.rint(10);
            if (mv == 0 && nseg > 1) {
                // boundary shift
                int j = rng.rint(nseg - 1);
                int d = rng.rint(2) ? 1 : -1;
                int ne = segEnd[j] + d;
                int l1 = ne - segStart(j) + 1, l2 = segEnd[j + 1] - ne;
                if (l1 < minLen || l2 < minLen) continue;
                int old = segEnd[j];
                segEnd[j] = ne;
                int c1 = segCost(j), c2 = segCost(j + 1);
                int nt = total - segC[j] - segC[j + 1] + c1 + c2;
                if (nt <= total || rng.rd() < exp((total - nt) / temp)) {
                    total = nt; segC[j] = c1; segC[j + 1] = c2;
                    if (d == 1) segId[ne] = j; else segId[old] = j + 1;
                } else segEnd[j] = old;
            } else {
                int ne = (rng.rint(4) == 0) ? 2 : 1;
                int len = 1 + rng.rint(min(T, 4));
                int lo = rng.rint(T - len + 1), hi = lo + len - 1;
                int es[2];
                es[0] = rng.rint(6);
                if (ne == 2) {
                    if (rng.rint(2)) es[1] = (es[0] + 3) % 6;
                    else { do es[1] = rng.rint(6); while (es[1] == es[0]); }
                }
                saved.assign(s.begin() + lo, s.begin() + hi + 1);
                bool okm = true;
                for (int q = 0; q < ne && okm; q++) {
                    if (!resample(es[q], lo, hi, out)) { okm = false; break; }
                    for (int t = lo; t <= hi; t++) s[t][es[q]] = out[t - lo];
                }
                if (!okm) { for (int t = lo; t <= hi; t++) s[t] = saved[t - lo]; continue; }
                for (int t = lo; t <= hi; t++) tickMasks(t);
                int j0 = segId[lo], j1 = segId[hi];
                int nt = total;
                static thread_local vector<int> nc;
                nc.resize(j1 - j0 + 1);
                for (int j = j0; j <= j1; j++) { nc[j - j0] = segCost(j); nt += nc[j - j0] - segC[j]; }
                if (nt <= total || rng.rd() < exp((total - nt) / temp)) {
                    total = nt;
                    for (int j = j0; j <= j1; j++) segC[j] = nc[j - j0];
                } else {
                    for (int t = lo; t <= hi; t++) { s[t] = saved[t - lo]; tickMasks(t); }
                }
            }
            if (total < bestTot) {
                bestTot = total; bestS = s; bestEnd = segEnd;
                if (total == 0) return true;
            }
        }
        s = bestS; segEnd = bestEnd; fullEval();
        return total == 0;
    }
};
