// Exhaustive search for L-cycles (cyclic schedules with one acceptance per period) via axis decomposition.
// Usage: exhL <input> <L> [idx list file or "all"] ; N <= 64 only.
#include "common.h"
typedef unsigned long long u64;

struct ESeq { vector<int> l; int d; vector<ES> raw; };
struct AOpt {
    vector<u64> one, two;
    int d;
    int ia, ib;           // emitter seq indices
    vector<int> code;     // per line: base-3 code of hits per tick
    unsigned touched = 0;
};

int main(int argc, char** argv) {
    auto pats = readInput(argv[1]);
    int L = atoi(argv[2]);
    vector<int> idx;
    if (argc > 3 && string(argv[3]) != "all") { ifstream f(argv[3]); int x; while (f >> x) idx.push_back(x); }
    else for (int i = 0; i < (int)pats.size(); i++) idx.push_back(i);
    int n = idx.size();
    vector<string> outLines(n);
    int P3 = 1; for (int t = 0; t < L; t++) P3 *= 3;
#pragma omp parallel for schedule(dynamic, 1)
    for (int ii = 0; ii < n; ii++) {
        int pi = idx[ii];
        const Pattern& p = pats[pi];
        int r = p.r, NL = 2 * r - 1;
        if (p.N > 64) { outLines[ii] = to_string(pi) + " skip"; continue; }
        auto t0 = chrono::steady_clock::now();
        u64 Pm = (p.nP == 64) ? ~0ULL : ((1ULL << p.nP) - 1);
        u64 ALL = (p.N == 64) ? ~0ULL : ((1ULL << p.N) - 1);
        u64 Km = ALL & ~Pm;
        u64 lm[3][24] = {};
        for (int i = 0; i < p.N; i++) for (int a = 0; a < 3; a++) lm[a][p.lineOf[i][a]] |= 1ULL << i;
        // emitter raw cyclic sequences for emitter with sign +1 (line = k)
        int S = 2 * (2 * r + 1);
        auto st = [&](int s) { return ES{(int8_t)(s / 2 - r), (int8_t)(s % 2)}; };
        map<vector<int>, int> seen;
        vector<ESeq> es;
        vector<int> cur(L);
        function<void(int)> rec = [&](int t) {
            if (t == L) {
                if (!adjES(st(cur[L - 1]), st(cur[0]))) return;
                vector<int> l(L);
                for (int q = 0; q < L; q++) { ES x = st(cur[q]); l[q] = (x.a && abs(x.k) < r) ? x.k + r - 1 : -1; }
                ES x0 = st(cur[0]);
                int d = x0.a ? abs(x0.k) : max(0, abs(x0.k) - 1);
                vector<ES> raw(L); for (int q = 0; q < L; q++) raw[q] = st(cur[q]);
                auto it = seen.find(l);
                if (it == seen.end()) { seen[l] = es.size(); es.push_back({l, d, raw}); }
                else if (d < es[it->second].d) { es[it->second].d = d; es[it->second].raw = raw; }
                return;
            }
            for (int s = 0; s < S; s++) {
                if (t > 0 && !adjES(st(cur[t - 1]), st(s))) continue;
                cur[t] = s;
                rec(t + 1);
            }
        };
        rec(0);
        vector<AOpt> ax[3];
        for (int a = 0; a < 3; a++) {
            map<vector<int>, int> sk;
            for (size_t i = 0; i < es.size(); i++)
                for (size_t j = i; j < es.size(); j++) {
                    vector<int> key(2 * L);
                    for (int t = 0; t < L; t++) { key[2 * t] = min(es[i].l[t], es[j].l[t]); key[2 * t + 1] = max(es[i].l[t], es[j].l[t]); }
                    int d = max(es[i].d, es[j].d);
                    auto it = sk.find(key);
                    if (it != sk.end()) {
                        AOpt& o = ax[a][it->second];
                        if (d < o.d) { o.d = d; o.ia = i; o.ib = j; }
                        continue;
                    }
                    AOpt o;
                    o.d = d; o.ia = i; o.ib = j;
                    o.one.assign(L, 0); o.two.assign(L, 0);
                    for (int t = 0; t < L; t++) {
                        int l = es[i].l[t], m = es[j].l[t];
                        if (l >= 0 && l == m) o.two[t] = lm[a][l];
                        else { if (l >= 0) o.one[t] |= lm[a][l]; if (m >= 0) o.one[t] |= lm[a][m]; }
                    }
                    o.code.assign(NL, 0);
                    for (int c = 0; c < NL; c++) {
                        int code = 0;
                        for (int t = L - 1; t >= 0; t--) code = code * 3 + (es[i].l[t] == c) + (es[j].l[t] == c);
                        o.code[c] = code;
                        if (code) o.touched |= 1u << c;
                    }
                    sk[key] = ax[a].size();
                    ax[a].push_back(o);
                }
        }
        long long cntValid = 0;
        int bestP = INT_MAX, bi = -1, bj = -1, bk = -1;
        vector<u64> s0(L), s1(L), s2(L), s3(L);
        vector<unsigned> allowedv(NL);
        vector<vector<u64>> allow(NL, vector<u64>((P3 + 63) / 64));
        vector<int> vv(L);
        for (size_t i = 0; i < ax[0].size(); i++) {
            const AOpt& H = ax[0][i];
            for (size_t j = 0; j < ax[1].size(); j++) {
                const AOpt& U = ax[1][j];
                int dHU = max(H.d, U.d);
                if (dHU >= bestP) continue;
                for (int t = 0; t < L; t++) {
                    u64 h1 = H.one[t], h2 = H.two[t], u1 = U.one[t], u2 = U.two[t];
                    u64 hh = h1 | h2, uu = u1 | u2;
                    s0[t] = ~(hh | uu);
                    s1[t] = (h1 & ~uu) | (u1 & ~hh);
                    s2[t] = (h2 & ~uu) | (u2 & ~hh) | (h1 & u1);
                    s3[t] = (h2 & uu) | (u2 & hh);
                }
                for (int c = 0; c < NL; c++) fill(allow[c].begin(), allow[c].end(), 0);
                for (int code = 0; code < P3; code++) {
                    int x = code;
                    for (int t = 0; t < L; t++) { vv[t] = x % 3; x /= 3; }
                    u64 und = ALL, bad = 0;
                    for (int t = L - 1; t >= 0; t--) {
                        int v = vv[t];
                        u64 eq = v == 0 ? s2[t] : v == 1 ? s1[t] : s0[t];
                        u64 ge = v == 0 ? s3[t] : v == 1 ? (s2[t] | s3[t]) : ~s0[t];
                        u64 xx = und & eq; bad |= xx & Km; und &= ~xx;
                        u64 yy = und & ge; bad |= yy & Pm; und &= ~yy;
                    }
                    bad |= und & Pm;
                    for (int c = 0; c < NL; c++) if (!(bad & lm[2][c])) allow[c][code >> 6] |= 1ULL << (code & 63);
                }
                unsigned req = 0; bool dead = false;
                for (int c = 0; c < NL; c++) {
                    bool anyb = false; for (u64 w : allow[c]) anyb |= (w != 0);
                    if (!anyb) { dead = true; break; }
                    if (!(allow[c][0] & 1)) req |= 1u << c;
                }
                if (dead || __builtin_popcount(req) > 4) continue;
                for (size_t k = 0; k < ax[2].size(); k++) {
                    const AOpt& V = ax[2][k];
                    if ((V.touched & req) != req) continue;
                    int d = max(dHU, V.d);
                    if (d >= bestP) continue;
                    bool ok = true;
                    unsigned tt = V.touched;
                    while (tt) { int c = __builtin_ctz(tt); tt &= tt - 1; int cd = V.code[c]; if (!(allow[c][cd >> 6] >> (cd & 63) & 1)) { ok = false; break; } }
                    if (!ok) continue;
                    cntValid++;
                    bestP = d; bi = i; bj = j; bk = k;
                }
            }
        }
        double dt = chrono::duration<double>(chrono::steady_clock::now() - t0).count();
        char buf[512];
        string real;
        if (bi >= 0) {
            // output raw realisation: per axis two emitter raw seqs (sign +1 convention, k = line offset)
            const AOpt* o[3] = {&ax[0][bi], &ax[1][bj], &ax[2][bk]};
            for (int a = 0; a < 3; a++) for (int q = 0; q < 2; q++) {
                const ESeq& E = es[q == 0 ? o[a]->ia : o[a]->ib];
                for (int t = 0; t < L; t++) { snprintf(buf, sizeof buf, " %d,%d", E.raw[t].k, E.raw[t].a); real += buf; }
                real += " ;";
            }
        }
        snprintf(buf, sizeof buf, "%d r=%d P=%d K=%d cnt=%d L=%d found=%d bestP=%d emitSeqs=%zu axOpts=%zu time=%.1f", pi, r, p.nP, p.nK, p.cnt, L,
                 bi >= 0, bestP == INT_MAX ? -1 : bestP, es.size(), ax[0].size(), dt);
        outLines[ii] = string(buf) + " |" + real;
        fprintf(stderr, "%s\n", buf);
    }
    for (auto& s : outLines) puts(s.c_str());
}
