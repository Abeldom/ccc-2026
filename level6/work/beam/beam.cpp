// Level 6 tick-level beam search over the exact simulation.
// Usage: beam <input> <baseline output> <out file> [options]
//   options: -W w -K k -M m -A a -B b -al alpha -g gamma -de delta -s1 scale1 -only i,j,.. -threads t
//            -cfg "W,K,M,A,B,alpha,gamma,delta,s1;..." (several settings, best per pattern kept)
#include <bits/stdc++.h>
#include <omp.h>
using namespace std;
typedef unsigned long long u64;

static const int DY[6] = {0, 1, 1, 0, -1, -1}, DX[6] = {2, 1, -1, -2, -1, 1};
static const int AX[6] = {0, 1, 2, 0, 1, 2};       // H U V H U V
static const int SG[6] = {-1, 1, 1, 1, -1, -1};

struct Problem {
    int r, count;
    vector<string> g;
    int np = 0, nk = 0, n = 0, nw = 0;
    vector<int> H, U, V;                      // line coords per cell (P first, then K)
    vector<array<u64, 8>> em[6];              // emitter e at offset k (index k+r): mask
    array<u64, 8> Pm{}, Km{};
};

static void build(Problem& pb) {
    int r = pb.r, R = 2 * r + 1, C = 4 * r + 1;
    vector<vector<int>> type(R, vector<int>(C, 0));
    for (int i = 0; i < R; i++) {
        int d = abs(i - r), m = 2 * r - d;
        for (int j = 0; j <= m; j++) type[i][d + 2 * j] = (d == r || j == 0 || j == m) ? 2 : 1;
    }
    vector<pair<int, int>> per;
    for (int k = 0; k < r; k++) per.push_back({r - k, k});
    for (int k = 0; k < r; k++) per.push_back({0, r + 2 * k});
    for (int k = 0; k < r; k++) per.push_back({k, 3 * r + k});
    for (int k = 0; k < r; k++) per.push_back({r + k, 4 * r - k});
    for (int k = 0; k < r; k++) per.push_back({2 * r, 3 * r - 2 * k});
    for (int k = 0; k < r; k++) per.push_back({2 * r - k, r - k});
    vector<vector<int>> idx(R, vector<int>(C, -1));
    vector<pair<int, int>> cells;
    for (int y = 0; y < R; y++)
        for (int x = 0; x < C; x++)
            if (pb.g[y][x] == 'M') { idx[y][x] = cells.size(); cells.push_back({y, x}); }
    pb.np = cells.size();
    int nby[6] = {0, 0, 1, 1, -1, -1}, nbx[6] = {2, -2, 1, -1, 1, -1};
    for (int i = 0; i < pb.np; i++)
        for (int d = 0; d < 6; d++) {
            int y = cells[i].first + nby[d], x = cells[i].second + nbx[d];
            if (y < 0 || y >= R || x < 0 || x >= C || type[y][x] != 1 || idx[y][x] >= 0) continue;
            idx[y][x] = cells.size(); cells.push_back({y, x});
        }
    pb.n = cells.size(); pb.nk = pb.n - pb.np; pb.nw = (pb.n + 63) / 64;
    if (pb.nw > 8) { fprintf(stderr, "too many cells\n"); exit(1); }
    for (auto [y, x] : cells) { pb.H.push_back(y - r); pb.U.push_back((x - y - r) / 2); pb.V.push_back((x + y - 3 * r) / 2); }
    for (int i = 0; i < pb.n; i++) {
        if (i < pb.np) pb.Pm[i >> 6] |= 1ULL << (i & 63); else pb.Km[i >> 6] |= 1ULL << (i & 63);
    }
    for (int e = 0; e < 6; e++) {
        pb.em[e].assign(2 * r + 1, array<u64, 8>{});
        for (int k = -r; k <= r; k++) {
            int pos = ((e * r + k) % (6 * r) + 6 * r) % (6 * r);
            int cy = per[pos].first + DY[e], cx = per[pos].second + DX[e];
            int cnt = 0;
            while (cy >= 0 && cy < R && cx >= 0 && cx < C && type[cy][cx] != 0) {
                if (type[cy][cx] == 1) {
                    cnt++;
                    int line = AX[e] == 0 ? cy - r : AX[e] == 1 ? (cx - cy - r) / 2 : (cx + cy - 3 * r) / 2;
                    if (line != SG[e] * k) { fprintf(stderr, "geometry mismatch\n"); exit(1); }
                    int id = idx[cy][cx];
                    if (id >= 0) pb.em[e][k + r][id >> 6] |= 1ULL << (id & 63);
                }
                cy += DY[e]; cx += DX[e];
            }
            if (abs(k) >= r && cnt) { fprintf(stderr, "dead zone mismatch\n"); exit(1); }
        }
    }
}

static int g_verbose = 0;
struct Params {
    int W = 500, K = 24, M = 3000;
    double A = 1.0, B = 1.0, alpha = 0.6, gamma = 0.6, delta = 0.5, s1 = 0.34;
    double lamP = 0, lamK = 0, lamM = 0, lamKmin = 0, lamPmin = 0;
    int dupC = 0, dupK = 0, tcap = 12, actAway = 0;
    bool set(const string& k, double v) {
        if (k == "W") W = (int)v; else if (k == "K") K = (int)v; else if (k == "M") M = (int)v;
        else if (k == "A") A = v; else if (k == "B") B = v; else if (k == "al") alpha = v; else if (k == "g") gamma = v;
        else if (k == "de") delta = v; else if (k == "s1") s1 = v; else if (k == "lp") lamP = v; else if (k == "lk") lamK = v;
        else if (k == "lm") lamM = v; else if (k == "dc") dupC = (int)v; else if (k == "dk") dupK = (int)v; else if (k == "tc") tcap = (int)v; else if (k == "ad") actAway = (int)v;
        else if (k == "lkm") lamKmin = v; else if (k == "lpm") lamPmin = v;
        else return false;
        return true;
    }
};

static inline u64 mix(u64 x) {
    x += 0x9e3779b97f4a7c15ULL; x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL; return x ^ (x >> 31);
}

template <int NW>
struct BS { u64 b[NW]; uint8_t em[6]; int16_t acc; float score; };

struct Cand { float score; int parent; uint8_t c[3]; };

template <int NW>
struct Beam {
    const Problem& pb; Params pr;
    int r, L;                                    // L = number of lines per axis = 2r-1
    float wt[64];
    Beam(const Problem& p, const Params& q) : pb(p), pr(q) {
        r = pb.r; L = 2 * r - 1;
        for (int t = 0; t < 64; t++) wt[t] = pow(pr.gamma, max(t, 1));
    }
    static inline int kOf(uint8_t s, int r) { return (s >> 1) - r; }
    // distance of emitter (state s) to line l (index l in [0,L), line value l-(r-1)) of its axis
    inline void dists(int e, uint8_t s, int8_t* d) const {
        int k = kOf(s, r), a = s & 1;
        for (int l = 0; l < L; l++) { int tk = SG[e] * (l - (r - 1)); int dk_ = abs(k - tk); d[l] = dk_ + (a ? (dk_ ? pr.actAway : 0) : 1); }
    }
    float exactScore(const BS<NW>& s) const {
        int8_t d[6][32]; int mn[3][32], dp[3][32];
        for (int e = 0; e < 6; e++) dists(e, s.em[e], d[e]);
        for (int a = 0; a < 3; a++)
            for (int l = 0; l < L; l++) { mn[a][l] = min(d[a][l], d[a + 3][l]); dp[a][l] = max(d[a][l], d[a + 3][l]); }
        int off = r - 1;
        float pot = 0, pen = 0; int pP = 0, pK = 0;
        int dC = pr.dupC, dK = pr.dupK, tc = pr.tcap;
        int maxT = 0, minTK = 1000, minTP = 1000;
        for (int i = 0; i < pb.n; i++) {
            bool on = (s.b[i >> 6] >> (i & 63)) & 1;
            int h = pb.H[i] + off, u = pb.U[i] + off, v = pb.V[i] + off;
            if (i < pb.np) {
                if (on) { pP++; continue; }
                int a = mn[0][h], b = mn[1][u], c = mn[2][v];
                int med = max(min(a, b), min(max(a, b), c));
                int t = min({med, dp[0][h] + dC, dp[1][u] + dC, dp[2][v] + dC});
                pot += pr.alpha * wt[t] - pr.lamP * min(t, tc);
                maxT = max(maxT, t); minTP = min(minTP, t);
            } else {
                if (!on) continue;
                pK++;
                int t = min({max({mn[0][h], mn[1][u], mn[2][v]}), max(dp[0][h], mn[1][u]) + dK, max(dp[0][h], mn[2][v]) + dK,
                             max(dp[1][u], mn[0][h]) + dK, max(dp[1][u], mn[2][v]) + dK, max(dp[2][v], mn[0][h]) + dK,
                             max(dp[2][v], mn[1][u]) + dK});
                pen += pr.B + pr.lamK * min(t, tc) - pr.B * pr.delta * wt[t];
                maxT = max(maxT, t); minTK = min(minTK, t);
            }
        }
        if (minTK < 1000) pen += pr.lamKmin * min(minTK, tc);
        if (minTP < 1000) pen += pr.lamPmin * min(minTP, tc);
        return (float)(s.acc * (10.0 * pb.np + 1000.0) + pr.A * pP - pen + pot - pr.lamM * min(maxT, tc));
    }
    u64 hashOf(const BS<NW>& s) const {
        u64 h = s.acc;
        for (int e = 0; e < 6; e++) h = h * 64 + s.em[e];
        h = mix(h);
        for (int w = 0; w < NW; w++) h = mix(h ^ s.b[w]);
        return h;
    }

    // run; returns program (codes) or empty on failure
    vector<uint16_t> run(int maxDepth) {
        vector<BS<NW>> cur(1), nxt;
        memset(&cur[0], 0, sizeof(BS<NW>));
        for (int e = 0; e < 6; e++) cur[0].em[e] = (uint8_t)(r * 2);   // k=0 inactive
        vector<vector<pair<int, uint16_t>>> hist;
        struct Combo { uint8_t s1, s2; uint16_t code; u64 a1[NW], a2[NW]; int8_t mn[32], dp[32]; };
        vector<Combo> combos[3];
        vector<Cand> cands; cands.reserve((size_t)pr.W * pr.K);
        float Tab[3][16][16], TD[3][16];
        int off = r - 1;
        double big = 10.0 * pb.np + 1000.0;
        for (int depth = 1; depth <= maxDepth; depth++) {
            cands.clear();
            int found = -1; uint16_t foundCode = 0;
            for (int pi = 0; pi < (int)cur.size() && found < 0; pi++) {
                const BS<NW>& ps = cur[pi];
                // actions per emitter
                uint8_t ns[6][4], nc[6][4]; int na[6];
                for (int e = 0; e < 6; e++) {
                    uint8_t s = ps.em[e]; int k = kOf(s, r); na[e] = 0;
                    ns[e][na[e]] = s; nc[e][na[e]++] = 0;
                    if (k + 1 <= r) { ns[e][na[e]] = s + 2; nc[e][na[e]++] = 1; }
                    if (k - 1 >= -r) { ns[e][na[e]] = s - 2; nc[e][na[e]++] = 2; }
                    ns[e][na[e]] = s ^ 1; nc[e][na[e]++] = 3;
                }
                for (int a = 0; a < 3; a++) {
                    combos[a].clear();
                    int e1 = a, e2 = a + 3;
                    for (int i = 0; i < na[e1]; i++)
                        for (int j = 0; j < na[e2]; j++) {
                            Combo c; c.s1 = ns[e1][i]; c.s2 = ns[e2][j];
                            c.code = (uint16_t)((nc[e1][i] << (2 * e1)) | (nc[e2][j] << (2 * e2)));
                            const u64* m1 = (c.s1 & 1) ? pb.em[e1][c.s1 >> 1].data() : nullptr;
                            const u64* m2 = (c.s2 & 1) ? pb.em[e2][c.s2 >> 1].data() : nullptr;
                            for (int w = 0; w < NW; w++) {
                                u64 x = m1 ? m1[w] : 0, y = m2 ? m2[w] : 0;
                                c.a1[w] = x ^ y; c.a2[w] = x & y;
                            }
                            int8_t d1[32], d2[32];
                            dists(e1, c.s1, d1); dists(e2, c.s2, d2);
                            for (int l = 0; l < L; l++) { c.mn[l] = min(d1[l], d2[l]); c.dp[l] = max(d1[l], d2[l]); }
                            combos[a].push_back(c);
                        }
                }
                // potential tables over parent's missing P cells
                int n0 = combos[0].size(), n1 = combos[1].size(), n2 = combos[2].size();
                for (int a = 0; a < 3; a++) for (int i = 0; i < 16; i++) { TD[a][i] = 0; for (int j = 0; j < 16; j++) Tab[a][i][j] = 0; }
                {
                    static thread_local vector<int> mh, mu, mv;
                    mh.clear(); mu.clear(); mv.clear();
                    for (int i = 0; i < pb.np; i++)
                        if (!((ps.b[i >> 6] >> (i & 63)) & 1)) { mh.push_back(pb.H[i] + off); mu.push_back(pb.U[i] + off); mv.push_back(pb.V[i] + off); }
                    int nm = mh.size();
                    if (nm) {
                        for (int i = 0; i < n0; i++) for (int j = 0; j < n1; j++) {
                            float s = 0; const int8_t* A = combos[0][i].mn; const int8_t* Bm = combos[1][j].mn;
                            for (int q = 0; q < nm; q++) s += wt[max(A[mh[q]], Bm[mu[q]])];
                            Tab[0][i][j] = s;
                        }
                        for (int i = 0; i < n0; i++) for (int j = 0; j < n2; j++) {
                            float s = 0; const int8_t* A = combos[0][i].mn; const int8_t* Bm = combos[2][j].mn;
                            for (int q = 0; q < nm; q++) s += wt[max(A[mh[q]], Bm[mv[q]])];
                            Tab[1][i][j] = s;
                        }
                        for (int i = 0; i < n1; i++) for (int j = 0; j < n2; j++) {
                            float s = 0; const int8_t* A = combos[1][i].mn; const int8_t* Bm = combos[2][j].mn;
                            for (int q = 0; q < nm; q++) s += wt[max(A[mu[q]], Bm[mv[q]])];
                            Tab[2][i][j] = s;
                        }
                        for (int i = 0; i < n0; i++) { float s = 0; for (int q = 0; q < nm; q++) s += wt[combos[0][i].dp[mh[q]]]; TD[0][i] = s; }
                        for (int i = 0; i < n1; i++) { float s = 0; for (int q = 0; q < nm; q++) s += wt[combos[1][i].dp[mu[q]]]; TD[1][i] = s; }
                        for (int i = 0; i < n2; i++) { float s = 0; for (int q = 0; q < nm; q++) s += wt[combos[2][i].dp[mv[q]]]; TD[2][i] = s; }
                    }
                }
                // per-parent top-K
                int Kp = pr.K; Cand top[256]; int nt = 0; int minI = -1; float minS = -1e30f;
                float al = pr.alpha * pr.s1;
                for (int i = 0; i < n0 && found < 0; i++) {
                    const Combo& ch = combos[0][i];
                    for (int j = 0; j < n1 && found < 0; j++) {
                        const Combo& cu = combos[1][j];
                        u64 c0[NW], c1[NW], c2[NW], c3[NW];
                        for (int w = 0; w < NW; w++) {
                            u64 a1 = ch.a1[w], a2 = ch.a2[w], b1 = cu.a1[w], b2 = cu.a2[w];
                            u64 a0 = ~(a1 | a2), b0 = ~(b1 | b2);
                            c0[w] = a0 & b0; c1[w] = (a1 & b0) | (a0 & b1);
                            c2[w] = (a2 & b0) | (a1 & b1) | (a0 & b2); c3[w] = (a2 & (b1 | b2)) | (a1 & b2);
                        }
                        float baseHU = Tab[0][i][j] + TD[0][i] + TD[1][j];
                        for (int k = 0; k < n2; k++) {
                            const Combo& cv = combos[2][k];
                            u64 nb[NW]; int pP = 0, pK = 0;
                            for (int w = 0; w < NW; w++) {
                                u64 g1 = cv.a1[w], g2 = cv.a2[w], g0 = ~(g1 | g2);
                                u64 S = (c2[w] & g0) | (c1[w] & g1) | (c0[w] & g2);
                                u64 D = c3[w] | (c2[w] & (g1 | g2)) | (c1[w] & g2);
                                nb[w] = (ps.b[w] | S) & ~D;
                                pP += __builtin_popcountll(nb[w] & pb.Pm[w]);
                                pK += __builtin_popcountll(nb[w] & pb.Km[w]);
                            }
                            float sc;
                            int acc = ps.acc;
                            if (pP == pb.np && pK == 0) {
                                acc++;
                                if (acc >= pb.count) { found = pi; foundCode = ch.code | cu.code | cv.code; break; }
                                sc = (float)(acc * big);
                            } else {
                                sc = (float)(acc * big + pr.A * pP - pr.B * pK + al * (baseHU + Tab[1][i][k] + Tab[2][j][k] + TD[2][k]));
                            }
                            if (nt < Kp) {
                                top[nt].score = sc; top[nt].parent = pi; top[nt].c[0] = i; top[nt].c[1] = j; top[nt].c[2] = k; nt++;
                                if (nt == Kp) { minI = 0; for (int q = 1; q < nt; q++) if (top[q].score < top[minI].score) minI = q; minS = top[minI].score; }
                            } else if (sc > minS) {
                                top[minI].score = sc; top[minI].c[0] = i; top[minI].c[1] = j; top[minI].c[2] = k;
                                minI = 0; for (int q = 1; q < nt; q++) if (top[q].score < top[minI].score) minI = q; minS = top[minI].score;
                            }
                        }
                    }
                }
                if (found >= 0) break;
                // need combos for reconstruction later: store the chosen children fully now
                for (int q = 0; q < nt; q++) cands.push_back(top[q]);
                // stash combos info into the candidate by recomputing later from parent state (deterministic)
            }
            if (found >= 0) {
                hist.push_back({{found, foundCode}});
                vector<uint16_t> prog(hist.size());
                int idx = 0;
                for (int d = (int)hist.size() - 1; d >= 0; d--) { prog[d] = hist[d][idx].second; idx = hist[d][idx].first; }
                return prog;
            }
            if (cands.empty()) return {};
            // global top-M by stage-1 score
            size_t M = min((size_t)pr.M, cands.size());
            if (M < cands.size()) nth_element(cands.begin(), cands.begin() + M, cands.end(), [](const Cand& a, const Cand& b) { return a.score > b.score; });
            cands.resize(M);
            // rebuild children states and exact scores
            struct Child { BS<NW> s; uint16_t code; int parent; u64 h; };
            static thread_local vector<Child> ch;
            ch.clear(); ch.reserve(M);
            for (auto& c : cands) {
                const BS<NW>& ps = cur[c.parent];
                // regenerate actions for this parent (same order as above)
                uint8_t ns[6][4], nc[6][4]; int na[6];
                for (int e = 0; e < 6; e++) {
                    uint8_t s = ps.em[e]; int k = kOf(s, r); na[e] = 0;
                    ns[e][na[e]] = s; nc[e][na[e]++] = 0;
                    if (k + 1 <= r) { ns[e][na[e]] = s + 2; nc[e][na[e]++] = 1; }
                    if (k - 1 >= -r) { ns[e][na[e]] = s - 2; nc[e][na[e]++] = 2; }
                    ns[e][na[e]] = s ^ 1; nc[e][na[e]++] = 3;
                }
                Child x; x.parent = c.parent; x.code = 0;
                x.s.acc = ps.acc;
                u64 cnt1[NW] = {}, cnt2[NW] = {}, cnt3[NW] = {};   // saturating hit counter (bit-sliced)
                for (int a = 0; a < 3; a++) {
                    int e1 = a, e2 = a + 3, ci = c.c[a];
                    int i = ci / na[e2], j = ci % na[e2];
                    x.s.em[e1] = ns[e1][i]; x.s.em[e2] = ns[e2][j];
                    x.code |= (uint16_t)((nc[e1][i] << (2 * e1)) | (nc[e2][j] << (2 * e2)));
                    for (int e : {e1, e2}) {
                        uint8_t s = x.s.em[e];
                        if (!(s & 1)) continue;
                        const u64* m = pb.em[e][s >> 1].data();
                        for (int w = 0; w < NW; w++) {
                            u64 t = m[w];
                            cnt3[w] |= cnt2[w] & t; cnt2[w] |= cnt1[w] & t; cnt1[w] |= t;
                        }
                    }
                }
                int pP = 0, pK = 0;
                for (int w = 0; w < NW; w++) {
                    u64 S = cnt2[w] & ~cnt3[w], D = cnt3[w];
                    x.s.b[w] = (ps.b[w] | S) & ~D;
                    pP += __builtin_popcountll(x.s.b[w] & pb.Pm[w]);
                    pK += __builtin_popcountll(x.s.b[w] & pb.Km[w]);
                }
                if (pP == pb.np && pK == 0) { x.s.acc++; for (int w = 0; w < NW; w++) x.s.b[w] = 0; }
                x.s.score = exactScore(x.s);
                x.h = hashOf(x.s);
                ch.push_back(x);
            }
            sort(ch.begin(), ch.end(), [](const Child& a, const Child& b) { return a.s.score > b.s.score; });
            nxt.clear();
            vector<pair<int, uint16_t>> hl;
            static thread_local unordered_set<u64> seen;
            seen.clear();
            for (auto& x : ch) {
                if ((int)nxt.size() >= pr.W) break;
                if (!seen.insert(x.h).second) continue;
                nxt.push_back(x.s); hl.push_back({x.parent, x.code});
            }
            hist.push_back(hl);
            swap(cur, nxt);
            if (g_verbose && !cur.empty()) {
                const BS<NW>& b = cur[0]; int pP = 0, pK = 0;
                for (int w = 0; w < NW; w++) { pP += __builtin_popcountll(b.b[w] & pb.Pm[w]); pK += __builtin_popcountll(b.b[w] & pb.Km[w]); }
                fprintf(stderr, "d%d n=%d best=%.2f acc=%d pP=%d/%d pK=%d em:", depth, (int)cur.size(), b.score, b.acc, pP, pb.np, pK);
                for (int e = 0; e < 6; e++) fprintf(stderr, " %d%c", kOf(b.em[e], r), (b.em[e] & 1) ? '*' : ' ');
                fprintf(stderr, " K:");
                for (int i = pb.np; i < pb.n; i++) if ((b.b[i >> 6] >> (i & 63)) & 1) fprintf(stderr, " (h%d u%d v%d)", pb.H[i], pb.U[i], pb.V[i]);
                fprintf(stderr, "\n");
            }
        }
        return {};
    }
};

// combos index encoding: c[a] = i*na[e2]+j  -> fix the stage-1 storage to use that encoding
// (stage 1 stores i (combo index within axis list) which is exactly i*na[e2]+j by construction)

static vector<uint16_t> solveOne(const Problem& pb, const Params& pr, int maxDepth) {
    switch (pb.nw) {
        case 1: { Beam<1> b(pb, pr); return b.run(maxDepth); }
        case 2: { Beam<2> b(pb, pr); return b.run(maxDepth); }
        case 3: { Beam<3> b(pb, pr); return b.run(maxDepth); }
        case 4: { Beam<4> b(pb, pr); return b.run(maxDepth); }
        case 5: { Beam<5> b(pb, pr); return b.run(maxDepth); }
        case 6: { Beam<6> b(pb, pr); return b.run(maxDepth); }
        case 7: { Beam<7> b(pb, pr); return b.run(maxDepth); }
        default: { Beam<8> b(pb, pr); return b.run(maxDepth); }
    }
}

// exact re-simulation (independent of the beam) for safety
static int simulate(const Problem& pb, const vector<string>& prog) {
    int r = pb.r;
    int k[6] = {0}, a[6] = {0};
    vector<char> mat(pb.n, 0);
    int acc = 0;
    for (auto& line : prog) {
        for (int e = 0; e < 6; e++) {
            char c = line[2 * e];
            if (c == 'F') k[e]++; else if (c == 'B') k[e]--; else if (c == 'T') a[e] ^= 1;
        }
        vector<int> hit(pb.n, 0);
        for (int e = 0; e < 6; e++) {
            if (!a[e]) continue;
            int kk = k[e];
            // normalise to [-3r, 3r)
            int pos = ((e * r + kk) % (6 * r) + 6 * r) % (6 * r);
            int kn = pos - e * r; while (kn < -3 * r) kn += 6 * r; while (kn >= 3 * r) kn -= 6 * r;
            if (abs(kn) > r) continue;    // dead zone (verified: no interior cells for |k|>=r)
            const auto& m = pb.em[e][kn + r];
            for (int i = 0; i < pb.n; i++) if ((m[i >> 6] >> (i & 63)) & 1) hit[i]++;
        }
        for (int i = 0; i < pb.n; i++) { if (hit[i] == 2) mat[i] = 1; else if (hit[i] >= 3) mat[i] = 0; }
        bool ok = true;
        for (int i = 0; i < pb.n && ok; i++) ok = (i < pb.np) ? mat[i] : !mat[i];
        if (ok) { acc++; for (int i = 0; i < pb.np; i++) mat[i] = 0; }
    }
    return acc;
}

static vector<vector<string>> readProgs(const string& path) {
    ifstream in(path); string line; vector<vector<string>> res; vector<string> cur;
    while (getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        bool blank = line.find_first_not_of(" \t") == string::npos;
        if (blank) { if (!cur.empty()) { res.push_back(cur); cur.clear(); } continue; }
        cur.push_back(line);
    }
    if (!cur.empty()) res.push_back(cur);
    return res;
}

int main(int argc, char** argv) {
    if (argc < 4) { fprintf(stderr, "usage\n"); return 1; }
    string inp = argv[1], basef = argv[2], outf = argv[3];
    vector<Params> cfgs; Params def; set<int> only; int threads = 16; string prevf;
    for (int i = 4; i < argc; i++) {
        string o = argv[i];
        auto nx = [&]() { return string(argv[++i]); };
        if (o == "-threads") threads = stoi(nx());
        else if (o == "-v") g_verbose = 1;
        else if (o == "-prev") prevf = nx();
        else if (o == "-only") {
            string s = nx(); stringstream ss(s); string t;
            while (getline(ss, t, ',')) {
                size_t d = t.find('-');
                if (d != string::npos && d > 0) { int a = stoi(t.substr(0, d)), b = stoi(t.substr(d + 1)); for (int q = a; q <= b; q++) only.insert(q); }
                else only.insert(stoi(t));
            }
        }
        else if (o == "-d") {        // default params: "k=v,k=v"
            string s = nx(); stringstream s2(s); string u;
            while (getline(s2, u, ',')) { size_t e = u.find('='); if (e == string::npos || !def.set(u.substr(0, e), stod(u.substr(e + 1)))) { fprintf(stderr, "bad param %s\n", u.c_str()); return 1; } }
        }
        else if (o == "-cfg") {      // several settings "k=v,k=v;k=v,..." (each on top of defaults)
            string s = nx(); stringstream ss(s); string t;
            while (getline(ss, t, ';')) {
                Params p = def; stringstream s2(t); string u;
                while (getline(s2, u, ',')) { size_t e = u.find('='); if (e == string::npos || !p.set(u.substr(0, e), stod(u.substr(e + 1)))) { fprintf(stderr, "bad param %s\n", u.c_str()); return 1; } }
                cfgs.push_back(p);
            }
        }
        else { fprintf(stderr, "unknown option %s\n", o.c_str()); return 1; }
    }
    if (cfgs.empty()) cfgs.push_back(def);
    for (auto& p : cfgs) p.K = min(p.K, 256);
    ifstream in(inp); int n; in >> n;
    vector<Problem> pbs(n);
    for (int i = 0; i < n; i++) {
        int rows, cols; in >> pbs[i].r >> rows >> cols >> pbs[i].count;
        pbs[i].g.resize(rows); for (auto& s : pbs[i].g) in >> s;
        build(pbs[i]);
    }
    auto base = readProgs(basef);
    if ((int)base.size() != n) { fprintf(stderr, "baseline count mismatch\n"); return 1; }
    vector<vector<string>> best = base;
    if (!prevf.empty()) {
        auto pv = readProgs(prevf);
        if ((int)pv.size() == n) for (int i = 0; i < n; i++) if (pv[i].size() < best[i].size() && simulate(pbs[i], pv[i]) >= pbs[i].count) best[i] = pv[i];
    }
    vector<int> baseLen(n); for (int i = 0; i < n; i++) baseLen[i] = base[i].size();
    vector<pair<int, int>> jobs;
    for (int i = 0; i < n; i++) if (only.empty() || only.count(i)) for (int c = 0; c < (int)cfgs.size(); c++) jobs.push_back({i, c});
    // larger patterns first
    sort(jobs.begin(), jobs.end(), [&](auto& a, auto& b) { return pbs[a.first].n * pbs[a.first].count > pbs[b.first].n * pbs[b.first].count; });
    omp_set_num_threads(threads);
    vector<int> res(jobs.size(), -1);
    auto t0 = chrono::steady_clock::now();
#pragma omp parallel for schedule(dynamic, 1)
    for (int j = 0; j < (int)jobs.size(); j++) {
        int i = jobs[j].first; const Params& pr = cfgs[jobs[j].second];
        int maxDepth;
#pragma omp critical
        maxDepth = (int)best[i].size() - 1;
        auto codes = solveOne(pbs[i], pr, maxDepth);
        if (codes.empty()) { res[j] = -1; continue; }
        vector<string> prog;
        for (auto c : codes) {
            string s;
            for (int e = 0; e < 6; e++) { if (e) s += ' '; s += "-FBT"[(c >> (2 * e)) & 3]; }
            prog.push_back(s);
        }
        int acc = simulate(pbs[i], prog);
        res[j] = prog.size();
#pragma omp critical
        {
            if (acc < pbs[i].count) fprintf(stderr, "BUG: pattern %d cfg %d invalid program (acc %d)\n", i, jobs[j].second, acc);
            else if (prog.size() < best[i].size()) best[i] = prog;
        }
    }
    double el = chrono::duration<double>(chrono::steady_clock::now() - t0).count();
    long long tot = 0, totBase = 0;
    for (int i = 0; i < n; i++) { tot += best[i].size(); totBase += baseLen[i]; }
    // per-pattern report
    for (int i = 0; i < n; i++) {
        if (!only.empty() && !only.count(i)) continue;
        printf("p%d r=%d |P|=%d |K|=%d cnt=%d base=%d best=%d |", i, pbs[i].r, pbs[i].np, pbs[i].nk, pbs[i].count, baseLen[i], (int)best[i].size());
        for (int j = 0; j < (int)jobs.size(); j++) if (jobs[j].first == i) printf(" %d", res[j]);
        printf("\n");
    }
    long long sel = 0, selBase = 0; int fails = 0;
    for (int i = 0; i < n; i++) if (only.empty() || only.count(i)) {
        sel += best[i].size(); selBase += baseLen[i];
        bool any = false; for (int j = 0; j < (int)jobs.size(); j++) if (jobs[j].first == i && res[j] > 0) any = true;
        if (!any) fails++;
    }
    printf("selected %lld (baseline %lld, fails %d)  total %lld (baseline %lld) time %.1fs\n", sel, selBase, fails, tot, totBase, el);
    ofstream out(outf, ios::binary);
    for (int i = 0; i < n; i++) {
        for (auto& l : best[i]) out << l << "\n";
        if (i + 1 < n) out << "\n";
    }
    return 0;
}
