// Level 6: tick-level beam search over the exact simulation, used as a large-neighbourhood
// "tail re-solve" on existing programs.
// For each pattern, take its current best program P (length L), simulate it to tick t (for many t),
// run a beam search from that exact state, and accept a new program  P[0:t] + path + P[T_a:]  when it is
// shorter, where path ends with acceptance number a in a configuration from which the old suffix
// (after the old a-th acceptance at tick T_a) still completes the job (checked by exact simulation;
// after an acceptance P and K cells are empty, so only the emitter configuration matters).
// a = count means plain completion (no suffix).
//
// Usage: tail <input> <best file> <out file> [options]
//   -threads n   -only list   -cfg "k=v,..;k=v,.."   -d "k=v,.."   -passes p   -step s
//   -dcap D (max beam depth per run)   -tl seconds (time limit per pattern)   -tmin frac   -v
//   -cuts all|acc  (acc: only cut at t=0 and right after acceptances (+ few ticks))
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
    vector<array<u64, 8>> em[6];              // emitter e at kidx = k + 3r (k in [-3r,3r)): mask
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
    int N = 6 * r;
    for (int e = 0; e < 6; e++) {
        pb.em[e].assign(N, array<u64, 8>{});
        for (int kidx = 0; kidx < N; kidx++) {
            int k = kidx - 3 * r;
            int pos = ((e * r + k) % N + N) % N;
            int cy = per[pos].first + DY[e], cx = per[pos].second + DX[e];
            int cnt = 0;
            while (cy >= 0 && cy < R && cx >= 0 && cx < C && type[cy][cx] != 0) {
                if (type[cy][cx] == 1) {
                    cnt++;
                    int line = AX[e] == 0 ? cy - r : AX[e] == 1 ? (cx - cy - r) / 2 : (cx + cy - 3 * r) / 2;
                    if (line != SG[e] * k) { fprintf(stderr, "geometry mismatch\n"); exit(1); }
                    int id = idx[cy][cx];
                    if (id >= 0) pb.em[e][kidx][id >> 6] |= 1ULL << (id & 63);
                }
                cy += DY[e]; cx += DX[e];
            }
            if (abs(k) >= r && cnt) { fprintf(stderr, "dead zone mismatch\n"); exit(1); }
        }
    }
}

// ---------------- exact simulation on P u K cells ----------------
struct St { u64 b[8]; uint8_t em[6]; int acc; };
static St initSt(const Problem& pb) { St s; memset(&s, 0, sizeof s); for (int e = 0; e < 6; e++) s.em[e] = (uint8_t)(6 * pb.r); return s; }
static inline void stepSt(const Problem& pb, St& s, uint16_t code) {
    int N = 6 * pb.r;
    for (int e = 0; e < 6; e++) {
        int c = (code >> (2 * e)) & 3; int kidx = s.em[e] >> 1, a = s.em[e] & 1;
        if (c == 1) kidx = (kidx + 1) % N; else if (c == 2) kidx = (kidx + N - 1) % N; else if (c == 3) a ^= 1;
        s.em[e] = (uint8_t)(kidx * 2 + a);
    }
    u64 c1[8] = {}, c2[8] = {}, c3[8] = {};
    for (int e = 0; e < 6; e++) {
        if (!(s.em[e] & 1)) continue;
        const u64* m = pb.em[e][s.em[e] >> 1].data();
        for (int w = 0; w < pb.nw; w++) { u64 t = m[w]; c3[w] |= c2[w] & t; c2[w] |= c1[w] & t; c1[w] |= t; }
    }
    bool ok = true;
    for (int w = 0; w < pb.nw; w++) {
        s.b[w] = (s.b[w] | (c2[w] & ~c3[w])) & ~c3[w];
        if ((s.b[w] & pb.Pm[w]) != pb.Pm[w] || (s.b[w] & pb.Km[w])) ok = false;
    }
    if (ok) { s.acc++; for (int w = 0; w < pb.nw; w++) s.b[w] &= ~pb.Pm[w]; }
}
// independent string simulator (mirror of sim6.py) for final safety
static int simulateStr(const Problem& pb, const vector<uint16_t>& prog) {
    St s = initSt(pb);
    for (auto c : prog) stepSt(pb, s, c);
    return s.acc;
}

static int g_verbose = 0;
struct Params {
    int W = 500, K = 24, M = 3000;
    double A = 1.0, B = 1.0, alpha = 0.6, gamma = 0.6, delta = 0.5, s1 = 0.34;
    double lamP = 0, lamK = 0, lamM = 0, lamKmin = 0, lamPmin = 0;
    int dupC = 0, dupK = 0, tcap = 12, actAway = 0;
    double noise = 0; int seed = 1;
    bool set(const string& k, double v) {
        if (k == "W") W = (int)v; else if (k == "K") K = (int)v; else if (k == "M") M = (int)v;
        else if (k == "A") A = v; else if (k == "B") B = v; else if (k == "al") alpha = v; else if (k == "g") gamma = v;
        else if (k == "de") delta = v; else if (k == "s1") s1 = v; else if (k == "lp") lamP = v; else if (k == "lk") lamK = v;
        else if (k == "lm") lamM = v; else if (k == "dc") dupC = (int)v; else if (k == "dk") dupK = (int)v; else if (k == "tc") tcap = (int)v; else if (k == "ad") actAway = (int)v;
        else if (k == "lkm") lamKmin = v; else if (k == "lpm") lamPmin = v; else if (k == "nz") noise = v; else if (k == "seed") seed = (int)v;
        else return false;
        return true;
    }
};

static inline u64 mix(u64 x) {
    x += 0x9e3779b97f4a7c15ULL; x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL; return x ^ (x >> 31);
}

static inline u64 mix0(u64 x) {
    x += 0x9e3779b97f4a7c15ULL; x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL; return x ^ (x >> 31);
}
static inline u64 cfgKey(const uint8_t* em, int acc) {
    u64 h = acc;
    for (int e = 0; e < 6; e++) h = h * 256 + em[e];
    return mix0(h ^ 0x5555ULL);
}
// context of the old program for suffix reuse
struct Ctx {
    int t;                       // prefix length (start state = state after t lines)
    int L;                       // old program length (trimmed: acceptance #count at tick L)
    vector<int> accT;            // accT[a] = tick (1-based) of old acceptance a (a=1..count)
    vector<array<uint8_t, 6>> cfgAt;   // emitter config right after tick accT[a]
    const vector<uint16_t>* prog;
    int limit;                   // max allowed total length
    int dcap;
    const unordered_map<u64, int>* oldH = nullptr;   // hash of old state after tau lines -> tau
    const unordered_map<u64, vector<int>>* cfgIdx = nullptr;   // (em config, acc) key -> list of tau (descending)
    const vector<St>* ost = nullptr;                 // old states
    chrono::steady_clock::time_point deadline = chrono::steady_clock::time_point::max();
};
struct Found { int total = INT_MAX; vector<uint16_t> path; int a = -1; int resume = -1; };

template <int NW>
struct BS { u64 b[NW]; uint8_t em[6]; int16_t acc; float score; };

struct Cand { float score; int parent; uint8_t c[3]; };

template <int NW>
struct Beam {
    const Problem& pb; Params pr;
    int r, L, N;
    float wt[64];
    mutable u64 rng = 88172645463325252ULL;
    inline float urand() const { rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17; return (rng >> 11) * (1.0f / 9007199254740992.0f); }
    Beam(const Problem& p, const Params& q) : pb(p), pr(q) {
        rng = mix(0x1234567ULL + (u64)q.seed * 7919ULL);
        r = pb.r; L = 2 * r - 1; N = 6 * r;
        for (int t = 0; t < 64; t++) wt[t] = pow(pr.gamma, max(t, 1));
    }
    inline int kOf(uint8_t s) const { return (s >> 1) - 3 * r; }
    inline void dists(int e, uint8_t s, int8_t* d, int aw) const {
        int k = kOf(s), a = s & 1;
        for (int l = 0; l < L; l++) {
            int tk = SG[e] * (l - (r - 1)); int dk = abs(k - tk); dk = min(dk, N - dk);
            d[l] = dk + (a ? (dk ? aw : 0) : 1);
        }
    }
    inline void actions(uint8_t s, uint8_t* ns, uint8_t* nc, int& na) const {
        int kidx = s >> 1, k = kidx - 3 * r, a = s & 1; na = 0;
        ns[na] = s; nc[na++] = 0;
        { int ki = (kidx + 1) % N, kk = ki - 3 * r; if (abs(kk) <= r || abs(kk) < abs(k)) { ns[na] = (uint8_t)(ki * 2 + a); nc[na++] = 1; } }
        { int ki = (kidx + N - 1) % N, kk = ki - 3 * r; if (abs(kk) <= r || abs(kk) < abs(k)) { ns[na] = (uint8_t)(ki * 2 + a); nc[na++] = 2; } }
        ns[na] = s ^ 1; nc[na++] = 3;
    }
    // heuristic score + admissible lower bound (ticks to next acceptance)
    float exactScore(const BS<NW>& s, int& lb) const {
        int8_t d[6][32], d0[6][32]; int mn[3][32], dp[3][32], mn0[3][32], dp0[3][32];
        for (int e = 0; e < 6; e++) { dists(e, s.em[e], d[e], pr.actAway); dists(e, s.em[e], d0[e], 0); }
        for (int a = 0; a < 3; a++)
            for (int l = 0; l < L; l++) {
                mn[a][l] = min(d[a][l], d[a + 3][l]); dp[a][l] = max(d[a][l], d[a + 3][l]);
                mn0[a][l] = min(d0[a][l], d0[a + 3][l]); dp0[a][l] = max(d0[a][l], d0[a + 3][l]);
            }
        int off = r - 1;
        float pot = 0, pen = 0; int pP = 0;
        int dC = pr.dupC, dK = pr.dupK, tc = pr.tcap;
        int maxT = 0, minTK = 1000, minTP = 1000; lb = 0;
        for (int i = 0; i < pb.n; i++) {
            bool on = (s.b[i >> 6] >> (i & 63)) & 1;
            int h = pb.H[i] + off, u = pb.U[i] + off, v = pb.V[i] + off;
            if (i < pb.np) {
                if (on) { pP++; continue; }
                int a = mn[0][h], b = mn[1][u], c = mn[2][v];
                int med = max(min(a, b), min(max(a, b), c));
                int t = min({med, dp[0][h] + dC, dp[1][u] + dC, dp[2][v] + dC});
                pot += pr.alpha * wt[min(t, 63)] - pr.lamP * min(t, tc);
                maxT = max(maxT, t); minTP = min(minTP, t);
                a = mn0[0][h]; b = mn0[1][u]; c = mn0[2][v];
                med = max(min(a, b), min(max(a, b), c));
                int t0 = min({med, dp0[0][h], dp0[1][u], dp0[2][v]});
                lb = max(lb, t0);
            } else {
                if (!on) continue;
                int t = min({max({mn[0][h], mn[1][u], mn[2][v]}), max(dp[0][h], mn[1][u]) + dK, max(dp[0][h], mn[2][v]) + dK,
                             max(dp[1][u], mn[0][h]) + dK, max(dp[1][u], mn[2][v]) + dK, max(dp[2][v], mn[0][h]) + dK,
                             max(dp[2][v], mn[1][u]) + dK});
                pen += pr.B + pr.lamK * min(t, tc) - pr.B * pr.delta * wt[min(t, 63)];
                maxT = max(maxT, t); minTK = min(minTK, t);
                int t0 = min({max({mn0[0][h], mn0[1][u], mn0[2][v]}), max(dp0[0][h], mn0[1][u]), max(dp0[0][h], mn0[2][v]),
                              max(dp0[1][u], mn0[0][h]), max(dp0[1][u], mn0[2][v]), max(dp0[2][v], mn0[0][h]), max(dp0[2][v], mn0[1][u])});
                lb = max(lb, t0);
            }
        }
        lb = max(lb, 1);
        if (minTK < 1000) pen += pr.lamKmin * min(minTK, tc);
        if (minTP < 1000) pen += pr.lamPmin * min(minTP, tc);
        float nzv = pr.noise > 0 ? (float)(pr.noise * urand()) : 0.f;
        return (float)(s.acc * (10.0 * pb.np + 1000.0) + pr.A * pP - pen + pot - pr.lamM * min(maxT, tc)) + nzv;
    }
    u64 hashOf(const BS<NW>& s) const {
        u64 h = s.acc;
        for (int e = 0; e < 6; e++) h = h * 256 + s.em[e];
        h = mix(h);
        for (int w = 0; w < NW; w++) h = mix(h ^ s.b[w]);
        return h;
    }
    u64 hashSt(const St& s) const {
        u64 h = s.acc;
        for (int e = 0; e < 6; e++) h = h * 256 + s.em[e];
        h = mix(h);
        for (int w = 0; w < NW; w++) h = mix(h ^ s.b[w]);
        return h;
    }
    // does the old suffix (after old acceptance a) complete from emitter config em with acc=a ?
    bool suffixOK(const Ctx& cx, int a, const uint8_t* em) const {
        bool same = true;
        for (int e = 0; e < 6; e++) if (em[e] != cx.cfgAt[a][e]) same = false;
        if (same) return true;
        St s; memset(&s, 0, sizeof s); for (int e = 0; e < 6; e++) s.em[e] = em[e]; s.acc = a;
        const auto& P = *cx.prog;
        for (int i = cx.accT[a]; i < cx.L && s.acc < pb.count; i++) stepSt(pb, s, P[i]);
        return s.acc >= pb.count;
    }

    // run from start state; updates fd when a better total is found
    void run(const St& st0, Ctx& cx, Found& fd) {
        vector<BS<NW>> cur(1), nxt;
        memset(&cur[0], 0, sizeof(BS<NW>));
        for (int w = 0; w < NW; w++) cur[0].b[w] = st0.b[w];
        for (int e = 0; e < 6; e++) cur[0].em[e] = st0.em[e];
        cur[0].acc = st0.acc;
        if (st0.acc >= pb.count) return;
        vector<vector<pair<int, uint16_t>>> hist;
        struct Combo { uint8_t s1, s2; uint16_t code; u64 a1[NW], a2[NW]; int8_t mn[32], dp[32]; };
        vector<Combo> combos[3];
        vector<Cand> cands; cands.reserve((size_t)pr.W * pr.K);
        float Tab[3][16][16], TD[3][16];
        int off = r - 1;
        double big = 10.0 * pb.np + 1000.0;
        auto pathTo = [&](int idx, uint16_t lastCode) {
            vector<uint16_t> p(hist.size() + 1);
            p[hist.size()] = lastCode;
            for (int d = (int)hist.size() - 1; d >= 0; d--) { p[d] = hist[d][idx].second; idx = hist[d][idx].first; }
            return p;
        };
        for (int depth = 1; cx.t + depth <= cx.limit && depth <= cx.dcap; depth++) {
            if ((depth & 3) == 0 && chrono::steady_clock::now() > cx.deadline) return;
            cands.clear();
            int found = -1; uint16_t foundCode = 0;
            for (int pi = 0; pi < (int)cur.size() && found < 0; pi++) {
                const BS<NW>& ps = cur[pi];
                uint8_t ns[6][4], nc[6][4]; int na[6];
                for (int e = 0; e < 6; e++) actions(ps.em[e], ns[e], nc[e], na[e]);
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
                            dists(e1, c.s1, d1, pr.actAway); dists(e2, c.s2, d2, pr.actAway);
                            for (int l = 0; l < L; l++) { c.mn[l] = min(d1[l], d2[l]); c.dp[l] = max(d1[l], d2[l]); }
                            combos[a].push_back(c);
                        }
                }
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
                            int pP = 0, pK = 0;
                            for (int w = 0; w < NW; w++) {
                                u64 g1 = cv.a1[w], g2 = cv.a2[w], g0 = ~(g1 | g2);
                                u64 S = (c2[w] & g0) | (c1[w] & g1) | (c0[w] & g2);
                                u64 D = c3[w] | (c2[w] & (g1 | g2)) | (c1[w] & g2);
                                u64 nb = (ps.b[w] | S) & ~D;
                                pP += __builtin_popcountll(nb & pb.Pm[w]);
                                pK += __builtin_popcountll(nb & pb.Km[w]);
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
                for (int q = 0; q < nt; q++) cands.push_back(top[q]);
            }
            if (found >= 0) {
                int total = cx.t + depth;
                if (total < fd.total) { fd.total = total; fd.path = pathTo(found, foundCode); fd.a = pb.count; fd.resume = cx.L; cx.limit = total - 1; }
                return;
            }
            if (cands.empty()) return;
            size_t M = min((size_t)pr.M, cands.size());
            if (M < cands.size()) nth_element(cands.begin(), cands.begin() + M, cands.end(), [](const Cand& a, const Cand& b) { return a.score > b.score; });
            cands.resize(M);
            struct Child { BS<NW> s; uint16_t code; int parent; u64 h; };
            static thread_local vector<Child> chv;
            vector<Child>& ch = *reinterpret_cast<vector<Child>*>(&chv);
            ch.clear(); ch.reserve(M);
            for (auto& c : cands) {
                const BS<NW>& ps = cur[c.parent];
                uint8_t ns[6][4], nc[6][4]; int na[6];
                for (int e = 0; e < 6; e++) actions(ps.em[e], ns[e], nc[e], na[e]);
                Child x; x.parent = c.parent; x.code = 0;
                x.s.acc = ps.acc;
                u64 cnt1[NW] = {}, cnt2[NW] = {}, cnt3[NW] = {};
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
                if (pP == pb.np && pK == 0) {
                    x.s.acc++; for (int w = 0; w < NW; w++) x.s.b[w] &= ~pb.Pm[w];
                    int a = x.s.acc;
                    if (a < pb.count) {
                        int total = cx.t + depth + (cx.L - cx.accT[a]);
                        if (total < fd.total && total <= cx.limit && suffixOK(cx, a, x.s.em)) {
                            fd.total = total; fd.path = pathTo(c.parent, x.code); fd.a = a; fd.resume = cx.accT[a]; cx.limit = total - 1;
                        }
                    }
                }
                int lb;
                x.s.score = exactScore(x.s, lb);
                // admissible pruning: completion needs >= lb ticks + one per further acceptance
                if (cx.t + depth + lb + (pb.count - x.s.acc - 1) > cx.limit) continue;
                x.h = hashOf(x.s);
                if (cx.cfgIdx) {
                    // dominance splice: same emitter config and acc as old state at tau, P superset, K subset
                    auto it = cx.cfgIdx->find(cfgKey(x.s.em, x.s.acc));
                    if (it != cx.cfgIdx->end()) {
                        for (int tau : it->second) {
                            int total = cx.t + depth + (cx.L - tau);
                            if (total >= fd.total || total > cx.limit) break;   // list is descending in tau
                            const St& o = (*cx.ost)[tau];
                            bool same = o.acc == x.s.acc;
                            for (int e = 0; e < 6 && same; e++) same = o.em[e] == x.s.em[e];
                            bool dom = same;
                            for (int w = 0; w < NW && dom; w++) dom = !(o.b[w] & pb.Pm[w] & ~x.s.b[w]) && !(x.s.b[w] & pb.Km[w] & ~o.b[w]);
                            if (dom) { fd.total = total; fd.path = pathTo(c.parent, x.code); fd.a = -3; fd.resume = tau; cx.limit = total - 1; break; }
                        }
                    }
                }
                if (cx.oldH) {
                    auto it = cx.oldH->find(x.h);
                    if (it != cx.oldH->end()) {
                        int tau = it->second, total = cx.t + depth + (cx.L - tau);
                        if (total < fd.total && total <= cx.limit) {
                            fd.total = total; fd.path = pathTo(c.parent, x.code); fd.a = -2; fd.resume = tau; cx.limit = total - 1;
                        }
                    }
                }
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
            if (nxt.empty()) return;
            hist.push_back(hl);
            swap(cur, nxt);
        }
    }
};

static void runBeam(const Problem& pb, const Params& pr, const St& st, Ctx& cx, Found& fd) {
    switch (pb.nw) {
        case 1: { Beam<1> b(pb, pr); b.run(st, cx, fd); break; }
        case 2: { Beam<2> b(pb, pr); b.run(st, cx, fd); break; }
        case 3: { Beam<3> b(pb, pr); b.run(st, cx, fd); break; }
        case 4: { Beam<4> b(pb, pr); b.run(st, cx, fd); break; }
        case 5: { Beam<5> b(pb, pr); b.run(st, cx, fd); break; }
        case 6: { Beam<6> b(pb, pr); b.run(st, cx, fd); break; }
        case 7: { Beam<7> b(pb, pr); b.run(st, cx, fd); break; }
        default: { Beam<8> b(pb, pr); b.run(st, cx, fd); break; }
    }
}

static u64 hashState(const Problem& pb, const St& s) {
    switch (pb.nw) {
        case 1: { Beam<1> b(pb, Params()); return b.hashSt(s); }
        case 2: { Beam<2> b(pb, Params()); return b.hashSt(s); }
        case 3: { Beam<3> b(pb, Params()); return b.hashSt(s); }
        case 4: { Beam<4> b(pb, Params()); return b.hashSt(s); }
        case 5: { Beam<5> b(pb, Params()); return b.hashSt(s); }
        case 6: { Beam<6> b(pb, Params()); return b.hashSt(s); }
        case 7: { Beam<7> b(pb, Params()); return b.hashSt(s); }
        default: { Beam<8> b(pb, Params()); return b.hashSt(s); }
    }
}

static vector<vector<uint16_t>> readProgs(const string& path) {
    ifstream in(path); string line; vector<vector<uint16_t>> res; vector<uint16_t> cur;
    while (getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        bool blank = line.find_first_not_of(" \t") == string::npos;
        if (blank) { if (!cur.empty()) { res.push_back(cur); cur.clear(); } continue; }
        istringstream ss(line); string tok; uint16_t c = 0; int e = 0;
        while (ss >> tok && e < 6) { int v = tok[0] == 'F' ? 1 : tok[0] == 'B' ? 2 : tok[0] == 'T' ? 3 : 0; c |= v << (2 * e); e++; }
        cur.push_back(c);
    }
    if (!cur.empty()) res.push_back(cur);
    return res;
}
static string lineOf(uint16_t c) { string s; for (int e = 0; e < 6; e++) { if (e) s += ' '; s += "-FBT"[(c >> (2 * e)) & 3]; } return s; }

static void writeProgs(const string& path, const vector<vector<uint16_t>>& P) {
    string tmp = path + ".tmp";
    { ofstream out(tmp, ios::binary);
      for (size_t i = 0; i < P.size(); i++) { for (auto c : P[i]) out << lineOf(c) << "\n"; if (i + 1 < P.size()) out << "\n"; } }
    rename(tmp.c_str(), path.c_str()) == 0 || (remove(path.c_str()), rename(tmp.c_str(), path.c_str()));
}

int main(int argc, char** argv) {
    if (argc < 4) { fprintf(stderr, "usage: tail input best out [opts]\n"); return 1; }
    string inp = argv[1], basef = argv[2], outf = argv[3];
    vector<Params> cfgs; Params def; set<int> only; int threads = 4, passes = 2, step = 1, dcap = 1000; double tl = 1e9, tminf = 0, tmaxf = 1.0, tlL = 0;
    string cuts = "all"; bool useHash = true;
    for (int i = 4; i < argc; i++) {
        string o = argv[i];
        auto nx = [&]() { return string(argv[++i]); };
        if (o == "-threads") threads = stoi(nx());
        else if (o == "-v") g_verbose = 1;
        else if (o == "-passes") passes = stoi(nx());
        else if (o == "-step") step = stoi(nx());
        else if (o == "-dcap") dcap = stoi(nx());
        else if (o == "-tl") tl = stod(nx());
        else if (o == "-tmin") tminf = stod(nx());
        else if (o == "-tmax") tmaxf = stod(nx());
        else if (o == "-tlL") tlL = stod(nx());
        else if (o == "-cuts") cuts = nx();
        else if (o == "-nohash") useHash = false;
        else if (o == "-only") {
            string s = nx(); stringstream ss(s); string t;
            while (getline(ss, t, ',')) {
                size_t d = t.find('-');
                if (d != string::npos && d > 0) { int a = stoi(t.substr(0, d)), b = stoi(t.substr(d + 1)); for (int q = a; q <= b; q++) only.insert(q); }
                else only.insert(stoi(t));
            }
        }
        else if (o == "-d") {
            string s = nx(); stringstream s2(s); string u;
            while (getline(s2, u, ',')) { size_t e = u.find('='); if (e == string::npos || !def.set(u.substr(0, e), stod(u.substr(e + 1)))) { fprintf(stderr, "bad param %s\n", u.c_str()); return 1; } }
        }
        else if (o == "-cfg") {
            string s = nx(); stringstream ss(s); string t;
            while (getline(ss, t, ';')) {
                Params p = def; stringstream s2(t); string u;
                while (getline(s2, u, ',')) { if (u.empty()) continue; size_t e = u.find('='); if (e == string::npos || !p.set(u.substr(0, e), stod(u.substr(e + 1)))) { fprintf(stderr, "bad param %s\n", u.c_str()); return 1; } }
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
    auto best = readProgs(basef);
    if ((int)best.size() != n) { fprintf(stderr, "best count mismatch\n"); return 1; }
    vector<int> origLen(n);
    for (int i = 0; i < n; i++) {
        if (simulateStr(pbs[i], best[i]) < pbs[i].count) { fprintf(stderr, "input program %d invalid!\n", i); return 1; }
        // trim after final acceptance
        St s = initSt(pbs[i]); int T = 0;
        for (int k = 0; k < (int)best[i].size(); k++) { stepSt(pbs[i], s, best[i][k]); if (s.acc >= pbs[i].count) { T = k + 1; break; } }
        best[i].resize(T);
        origLen[i] = best[i].size();
    }
    vector<int> jobs;
    for (int i = 0; i < n; i++) if (only.empty() || only.count(i)) jobs.push_back(i);
    sort(jobs.begin(), jobs.end(), [&](int a, int b) { return pbs[a].n * (long)best[a].size() > pbs[b].n * (long)best[b].size(); });
    omp_set_num_threads(threads);
    auto t0 = chrono::steady_clock::now();
    auto lastWrite = t0;
    mutex mtx;
    int done = 0;
#pragma omp parallel for schedule(dynamic, 1)
    for (int jj = 0; jj < (int)jobs.size(); jj++) {
        int i = jobs[jj];
        const Problem& pb = pbs[i];
        vector<uint16_t> prog;
        { lock_guard<mutex> lk(mtx); prog = best[i]; }
        auto ts = chrono::steady_clock::now();
        double tlp = tlL > 0 ? tlL * prog.size() : tl;
        auto elapsed = [&]() { return chrono::duration<double>(chrono::steady_clock::now() - ts).count(); };
        for (int pass = 0; pass < passes && elapsed() < tlp; pass++) {
            bool imp = false;
            int Lp = prog.size();
            // states
            vector<St> st(Lp + 1); st[0] = initSt(pb);
            for (int k = 0; k < Lp; k++) { st[k + 1] = st[k]; stepSt(pb, st[k + 1], prog[k]); }
            auto rebuild = [&]() {
                Lp = prog.size(); st.assign(Lp + 1, St()); st[0] = initSt(pb);
                for (int k = 0; k < Lp; k++) { st[k + 1] = st[k]; stepSt(pb, st[k + 1], prog[k]); }
            };
            int tLow = (int)(tminf * Lp);
            int tHigh = min(Lp - 2, (int)(tmaxf * Lp));
            for (int t = tHigh; t >= tLow && elapsed() < tlp; t -= step) {
                if (t >= (int)prog.size() - 1) continue;
                if (cuts == "acc") {   // only t=0 or within 2 ticks after an acceptance
                    bool okc = (t == 0);
                    for (int q = 0; q <= 2 && !okc; q++) if (t - q >= 1 && st[t - q].acc > st[t - q - 1].acc) okc = true;
                    if (!okc) continue;
                }
                for (size_t ci = 0; ci < cfgs.size(); ci++) {
                    Ctx cx; cx.t = t; cx.L = prog.size(); cx.prog = &prog; cx.dcap = dcap; cx.limit = cx.L - 1;
                    cx.deadline = ts + chrono::milliseconds((long long)(min(tlp, 1e7) * 1000));
                    cx.accT.assign(pb.count + 1, 0); cx.cfgAt.assign(pb.count + 1, {});
                    for (int k = 1; k <= cx.L; k++) if (st[k].acc > st[k - 1].acc && st[k].acc <= pb.count) {
                        cx.accT[st[k].acc] = k; for (int e = 0; e < 6; e++) cx.cfgAt[st[k].acc][e] = st[k].em[e];
                    }
                    unordered_map<u64, int> oldH;
                    unordered_map<u64, vector<int>> cfgIdx;
                    if (useHash) {
                        for (int k = cx.L; k >= t + 2; k--) cfgIdx[cfgKey(st[k].em, st[k].acc)].push_back(k);
                        cx.cfgIdx = &cfgIdx; cx.ost = &st;
                    }
                    Found fd;
                    runBeam(pb, cfgs[ci], st[t], cx, fd);
                    if (fd.resume < 0) continue;
                    vector<uint16_t> np(prog.begin(), prog.begin() + t);
                    np.insert(np.end(), fd.path.begin(), fd.path.end());
                    np.insert(np.end(), prog.begin() + fd.resume, prog.end());
                    // trim + verify
                    St s = initSt(pb); int T = -1;
                    for (int k = 0; k < (int)np.size(); k++) { stepSt(pb, s, np[k]); if (s.acc >= pb.count) { T = k + 1; break; } }
                    if (T < 0) { fprintf(stderr, "BUG: p%d t=%d invalid new program (a=%d)\n", i, t, fd.a); continue; }
                    np.resize(T);
                    if (np.size() < prog.size()) {
                        if (g_verbose) fprintf(stderr, "p%d t=%d cfg%zu: %zu -> %zu (a=%d)\n", i, t, ci, prog.size(), np.size(), fd.a);
                        prog = np; imp = true; rebuild();
                        lock_guard<mutex> lk(mtx);
                        if (prog.size() < best[i].size()) best[i] = prog;
                        double sinceW = chrono::duration<double>(chrono::steady_clock::now() - lastWrite).count();
                        if (sinceW > 20) { writeProgs(outf, best); lastWrite = chrono::steady_clock::now(); }
                    }
                }
            }
            if (!imp) break;
        }
        {
            lock_guard<mutex> lk(mtx);
            done++;
            printf("p%d r=%d |P|=%d cnt=%d %d -> %d  (%.1fs) [%d/%zu]\n", i, pb.r, pb.np, pb.count, origLen[i], (int)best[i].size(), elapsed(), done, jobs.size());
            fflush(stdout);
        }
    }
    double el = chrono::duration<double>(chrono::steady_clock::now() - t0).count();
    long long tot = 0, tot0 = 0; int nimp = 0;
    for (int i = 0; i < n; i++) { tot += best[i].size(); tot0 += origLen[i]; if ((int)best[i].size() < origLen[i]) nimp++; }
    for (int i = 0; i < n; i++) if (simulateStr(pbs[i], best[i]) < pbs[i].count) fprintf(stderr, "FINAL BUG p%d\n", i);
    printf("total %lld -> %lld, improved %d patterns, time %.1fs\n", tot0, tot, nimp, el);
    writeProgs(outf, best);
    return 0;
}
