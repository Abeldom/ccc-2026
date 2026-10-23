// Exhaustive search for static configs (L=1) and 2-builds (= 2-cycles) per pattern (N <= 64 cells).
// Output: one line per pattern: idx r |P| cnt lb1 staticCost(-1) twoCount bestP  + realization
#include "common.h"
typedef unsigned long long u64;

struct EOpt { int l1, l2; int d; };  // line index or -1, prefix need
struct AOpt {
    u64 one[2], two[2];
    int d;
    EOpt a, b;
    // per line codes (for V axis): v1*3+v2
    int code[16];
    unsigned touched;
};

int main(int argc, char** argv) {
    auto pats = readInput(argv[1]);
    int n = pats.size();
    vector<string> outLines(n);
#pragma omp parallel for schedule(dynamic, 1)
    for (int pi = 0; pi < n; pi++) {
        const Pattern& p = pats[pi];
        int r = p.r, NL = 2 * r - 1;
        if (p.N > 64) { outLines[pi] = to_string(pi) + " skip"; continue; }
        u64 Pm = (p.nP == 64) ? ~0ULL : ((1ULL << p.nP) - 1);
        u64 ALL = (p.N == 64) ? ~0ULL : ((1ULL << p.N) - 1);
        u64 Km = ALL & ~Pm;
        u64 lm[3][16] = {};
        for (int i = 0; i < p.N; i++)
            for (int a = 0; a < 3; a++) lm[a][p.lineOf[i][a]] |= 1ULL << i;
        auto absl = [&](int li) { return abs(li - (r - 1)); };
        // ---------- static ----------
        int bestStatic = INT_MAX;
        {
            vector<vector<int>> opts;
            opts.push_back({});
            for (int a = 0; a < NL; a++) {
                opts.push_back({a});
                opts.push_back({a, a});
                for (int b = a + 1; b < NL; b++) opts.push_back({a, b});
            }
            auto mk = [&](int ax, const vector<int>& o, u64& one, u64& two) {
                one = two = 0;
                if (o.size() == 1) one = lm[ax][o[0]];
                else if (o.size() == 2) {
                    if (o[0] == o[1]) two = lm[ax][o[0]];
                    else one = lm[ax][o[0]] | lm[ax][o[1]];
                }
            };
            auto mx = [&](const vector<int>& o) { int m = 0; for (int x : o) m = max(m, absl(x)); return m; };
            for (auto& A : opts) {
                u64 a1, a2; mk(0, A, a1, a2);
                for (auto& B : opts) {
                    u64 b1, b2; mk(1, B, b1, b2);
                    for (auto& C : opts) {
                        int m = max({mx(A), mx(B), mx(C)});
                        if (m >= bestStatic) continue;
                        u64 c1, c2; mk(2, C, c1, c2);
                        u64 aa = a1 | a2, bb = b1 | b2, cc = c1 | c2;
                        u64 eq2 = (a2 & ~bb & ~cc) | (b2 & ~aa & ~cc) | (c2 & ~aa & ~bb) | (a1 & b1 & ~cc) | (a1 & c1 & ~bb) | (b1 & c1 & ~aa);
                        if ((eq2 & Pm) == Pm && (eq2 & Km) == 0) bestStatic = m;
                    }
                }
            }
        }
        // ---------- 2-builds ----------
        vector<EOpt> eo;
        eo.push_back({-1, -1, 0});
        for (int c = 0; c < NL; c++) {
            int ac = absl(c);
            eo.push_back({c, c, ac});
            eo.push_back({c, -1, ac});
            eo.push_back({-1, c, max(0, ac - 1)});
            if (c + 1 < NL) eo.push_back({c, c + 1, ac});
            if (c - 1 >= 0) eo.push_back({c, c - 1, ac});
        }
        vector<AOpt> ax[3];
        for (int a = 0; a < 3; a++) {
            map<array<int, 4>, int> seen;
            for (size_t i = 0; i < eo.size(); i++)
                for (size_t j = i; j < eo.size(); j++) {
                    const EOpt &x = eo[i], &y = eo[j];
                    array<int, 4> key = {min(x.l1, y.l1), max(x.l1, y.l1), min(x.l2, y.l2), max(x.l2, y.l2)};
                    int d = max(x.d, y.d);
                    auto it = seen.find(key);
                    if (it != seen.end()) {
                        AOpt& o = ax[a][it->second];
                        if (d < o.d) { o.d = d; o.a = x; o.b = y; }
                        continue;
                    }
                    AOpt o{};
                    o.d = d; o.a = x; o.b = y;
                    int ls[2][2] = {{x.l1, y.l1}, {x.l2, y.l2}};
                    for (int t = 0; t < 2; t++) {
                        int l = ls[t][0], m = ls[t][1];
                        o.one[t] = o.two[t] = 0;
                        if (l >= 0 && m >= 0 && l == m) o.two[t] = lm[a][l];
                        else {
                            if (l >= 0) o.one[t] |= lm[a][l];
                            if (m >= 0) o.one[t] |= lm[a][m];
                        }
                    }
                    for (int c = 0; c < NL; c++) {
                        int v1 = (x.l1 == c) + (y.l1 == c), v2 = (x.l2 == c) + (y.l2 == c);
                        o.code[c] = v1 * 3 + v2;
                        if (o.code[c]) o.touched |= 1u << c;
                    }
                    seen[key] = ax[a].size();
                    ax[a].push_back(o);
                }
        }
        long long cntValid = 0;
        int bestP = INT_MAX;
        int bi = -1, bj = -1, bk = -1;
        for (size_t i = 0; i < ax[0].size(); i++) {
            const AOpt& H = ax[0][i];
            for (size_t j = 0; j < ax[1].size(); j++) {
                const AOpt& U = ax[1][j];
                int dHU = max(H.d, U.d);
                u64 s0[2], s1[2], s2[2], s3[2];
                for (int t = 0; t < 2; t++) {
                    u64 h1 = H.one[t], h2 = H.two[t], u1 = U.one[t], u2 = U.two[t];
                    u64 hh = h1 | h2, uu = u1 | u2;
                    s0[t] = ~(hh | uu);
                    s1[t] = (h1 & ~uu) | (u1 & ~hh);
                    s2[t] = (h2 & ~uu) | (u2 & ~hh) | (h1 & u1);
                    s3[t] = (h2 & uu) | (u2 & hh);
                }
                u64 eq2[2][3], ge3[2][3], le1[2][3];
                for (int t = 0; t < 2; t++) {
                    eq2[t][0] = s2[t]; ge3[t][0] = s3[t]; le1[t][0] = s0[t] | s1[t];
                    eq2[t][1] = s1[t]; ge3[t][1] = s2[t] | s3[t]; le1[t][1] = s0[t];
                    eq2[t][2] = s0[t]; ge3[t][2] = ~s0[t]; le1[t][2] = 0;
                }
                unsigned allowed[16];
                for (int c = 0; c < NL; c++) allowed[c] = 0;
                for (int v1 = 0; v1 < 3; v1++)
                    for (int v2 = 0; v2 < 3; v2++) {
                        u64 goodP = eq2[1][v2] | (eq2[0][v1] & le1[1][v2]);
                        u64 badK = eq2[1][v2] | (eq2[0][v1] & ~ge3[1][v2]);
                        u64 bad = (Pm & ~goodP) | (Km & badK);
                        for (int c = 0; c < NL; c++)
                            if (!(bad & lm[2][c])) allowed[c] |= 1u << (v1 * 3 + v2);
                    }
                unsigned req = 0;
                bool dead = false;
                for (int c = 0; c < NL; c++) {
                    if (!allowed[c]) { dead = true; break; }
                    if (!(allowed[c] & 1)) req |= 1u << c;
                }
                if (dead || __builtin_popcount(req) > 4) continue;
                for (size_t k = 0; k < ax[2].size(); k++) {
                    const AOpt& V = ax[2][k];
                    if ((V.touched & req) != req) continue;
                    bool ok = true;
                    unsigned t = V.touched;
                    while (t) {
                        int c = __builtin_ctz(t); t &= t - 1;
                        if (!(allowed[c] >> V.code[c] & 1)) { ok = false; break; }
                    }
                    if (!ok) continue;
                    cntValid++;
                    int d = max(dHU, V.d);
                    if (d < bestP) { bestP = d; bi = i; bj = j; bk = k; }
                }
            }
        }
        char buf[512];
        string real;
        if (bi >= 0) {
            const AOpt* o[3] = {&ax[0][bi], &ax[1][bj], &ax[2][bk]};
            for (int a = 0; a < 3; a++) {
                snprintf(buf, sizeof buf, " %d %d %d %d", o[a]->a.l1, o[a]->a.l2, o[a]->b.l1, o[a]->b.l2);
                real += buf;
            }
        }
        snprintf(buf, sizeof buf, "%d r=%d P=%d K=%d cnt=%d lb1=%d static=%d two=%lld bestP=%d", pi, r, p.nP, p.nK, p.cnt, p.lb1,
                 bestStatic == INT_MAX ? -1 : bestStatic, cntValid, bestP == INT_MAX ? -1 : bestP);
        outLines[pi] = string(buf) + " |" + real;
    }
    for (auto& s : outLines) puts(s.c_str());
}
