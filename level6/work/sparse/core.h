// sparse: exact-ish search for patterns with few matter cells.
// Works in "line coordinates": emitter e covers line c on axis AX[e] (c = SG[e]*k).
// Emitter state: active at line c (|c| <= r-1), or inactive with (p = toggle-off position, s = slack:
// number of idle ticks during which it could have moved).  Moves of inactive emitters are placed
// retroactively when the program is reconstructed.
#pragma once
#include "geo.h"
typedef unsigned __int128 M;
static inline int pc128(M x) { return __builtin_popcountll((u64)x) + __builtin_popcountll((u64)(x >> 64)); }

constexpr int MAXL = 23;

struct PInfo {
    int r = 0, cnt = 0, N = 0, nP = 0, nK = 0, NL = 0, SC = 0, NCODE = 0, UNIV = 0;
    M Pm = 0, Km = 0;
    M lineM[3][MAXL];
    M reach[3][MAXL][MAXL];  // [axis][center idx][d]: lines within distance d of center (d capped NL-1)
    vector<int> cellGeo;     // my cell -> geo cell index
    // code tables
    vector<int8_t> cAct, cC, cP, cS;       // decode
    vector<vector<uint16_t>> succ;         // successor codes
    int code(bool act, int c, int p, int s) const {
        if (act) return c + (r - 1);
        if (s >= (r - 1) + abs(p)) { p = 0; s = r - 1; }
        if (s > SC) s = SC;
        return NL + (p + r - 1) * (SC + 1) + s;
    }
    // mask of cells this emitter (axis a) can cover now-or-within d remaining ticks (see notes)
    inline M emask(int a, int cd, int d) const {
        if (cAct[cd]) return reach[a][cC[cd] + r - 1][min(d, NL - 1)];
        if (d < 1) return 0;
        int dd = cS[cd] + d - 1; if (dd > NL - 1) dd = NL - 1;
        return reach[a][cP[cd] + r - 1][dd];
    }
    inline M nowmask(int a, int cd) const { return cAct[cd] ? lineM[a][cC[cd] + r - 1] : (M)0; }
};

static PInfo makePInfo(const Geo& g, const Pat& pt) {
    PInfo P; P.r = g.r; P.cnt = pt.cnt; int r = g.r;
    P.NL = 2 * r - 1; P.SC = 2 * r - 2;
    vector<int> cells;
    bforeach(pt.P, [&](int c) { cells.push_back(c); });
    P.nP = cells.size();
    bforeach(pt.K, [&](int c) { cells.push_back(c); });
    P.N = cells.size(); P.nK = P.N - P.nP;
    if (P.N > 128) { P.N = -1; return P; }
    P.cellGeo = cells;
    for (int i = 0; i < P.N; i++) { if (i < P.nP) P.Pm |= (M)1 << i; else P.Km |= (M)1 << i; }
    for (int a = 0; a < 3; a++) for (int l = 0; l < MAXL; l++) P.lineM[a][l] = 0;
    for (int i = 0; i < P.N; i++) for (int a = 0; a < 3; a++) P.lineM[a][g.L[a][cells[i]] + r - 1] |= (M)1 << i;
    for (int a = 0; a < 3; a++) for (int c = 0; c < P.NL; c++) for (int d = 0; d < P.NL; d++) {
        M m = 0; for (int c2 = max(0, c - d); c2 <= min(P.NL - 1, c + d); c2++) m |= P.lineM[a][c2];
        P.reach[a][c][d] = m;
    }
    P.NCODE = P.NL + P.NL * (P.SC + 1) + 1;   // last code = UNIV (universal start: any state at tick 1)
    P.UNIV = P.NCODE - 1;
    P.cAct.assign(P.NCODE, 0); P.cC.assign(P.NCODE, 0); P.cP.assign(P.NCODE, 0); P.cS.assign(P.NCODE, 0);
    for (int c = -(r - 1); c <= r - 1; c++) { int cd = P.code(true, c, 0, 0); P.cAct[cd] = 1; P.cC[cd] = c; }
    for (int p = -(r - 1); p <= r - 1; p++) for (int s = 0; s <= P.SC; s++) {
        int cd = P.NL + (p + r - 1) * (P.SC + 1) + s; P.cAct[cd] = 0; P.cP[cd] = p; P.cS[cd] = s;
    }
    P.succ.assign(P.NCODE, {});
    P.cAct[P.UNIV] = 0; P.cP[P.UNIV] = 0; P.cS[P.UNIV] = P.SC;
    for (int c = -(r - 1); c <= r - 1; c++) P.succ[P.UNIV].push_back(P.code(true, c, 0, 0));
    P.succ[P.UNIV].push_back(P.code(false, 0, 0, P.SC));
    for (int cd = 0; cd < P.NCODE - 1; cd++) {
        auto& v = P.succ[cd];
        if (P.cAct[cd]) {
            int c = P.cC[cd];
            v.push_back(cd);
            if (c + 1 <= r - 1) v.push_back(P.code(true, c + 1, 0, 0));
            if (c - 1 >= -(r - 1)) v.push_back(P.code(true, c - 1, 0, 0));
            v.push_back(P.code(false, 0, c, 0));
        } else {
            int p = P.cP[cd], s = P.cS[cd];
            v.push_back(P.code(false, 0, p, s + 1));
            for (int c = max(-(r - 1), p - s); c <= min(r - 1, p + s); c++) v.push_back(P.code(true, c, 0, 0));
        }
    }
    return P;
}

typedef array<uint16_t, 6> Cfg;
static inline u64 cfgKey(const Cfg& c) {
    u64 k = 0;
    for (int a = 0; a < 3; a++) {
        uint16_t x = c[a], y = c[a + 3]; if (x > y) swap(x, y);
        k = (k << 10) | x; k = (k << 10) | y;
    }
    return k;
}
struct SKey { u64 c; M b; bool operator==(const SKey& o) const { return c == o.c && b == o.b; } };
struct SKeyH { size_t operator()(const SKey& k) const {
    u64 h = k.c * 0x9E3779B97F4A7C15ULL ^ (u64)k.b * 0xC2B2AE3D27D4EB4FULL ^ (u64)(k.b >> 64) * 0x165667B19E3779F9ULL;
    return h ^ (h >> 29); } };

// Tick: returns new board given new config.
static inline M tickBoard(const PInfo& P, const Cfg& c, M board) {
    M h1 = 0, h2 = 0, h3 = 0;
    for (int e = 0; e < 6; e++) {
        M m = P.nowmask(AX[e], c[e]);
        h3 |= h2 & m; h2 |= h1 & m; h1 |= m;
    }
    return (board | (h2 & ~h3)) & ~h3;
}

// ---------------- build search: from config X with clean board, all builds of exactly len ticks -----------
struct BuildRes { Cfg end; vector<Cfg> path; };   // path[0..len-1] = configs after each tick

struct BuildSearch {
    const PInfo& P;
    int len = 0;
    long long nodes = 0, nodeLimit = 2000000000LL;
    size_t maxEnds = 1000000;
    bool aborted = false;
    bool firstOnly = false;
    vector<unordered_set<SKey, SKeyH>> vis;
    unordered_map<u64, BuildRes> ends;
    vector<Cfg> stk;
    BuildSearch(const PInfo& p) : P(p) {}

    // per-level scratch
    struct Lvl {
        int nopt[6]; uint16_t opt[6][32]; M om[6][32]; M nm[6][32];
        M s1[7], s2[7], s3[7];
        M missP, presK, board; int R; Cfg cur, nw; int chosenIdx[6];
    };
    vector<Lvl> L;

    bool run(const Cfg& X, int len_) {
        len = len_; nodes = 0; aborted = false; ends.clear();
        vis.assign(len + 1, {});
        L.assign(len + 1, Lvl());
        stk.assign(len, Cfg());
        dfs(X, 0, len, 0);
        return !ends.empty();
    }
    // node prune check
    bool nodeOK(const Cfg& c, M board, int R) const {
        M missP = P.Pm & ~board, presK = P.Km & board;
        M a1 = 0, a2 = 0, a3 = 0;
        for (int e = 0; e < 6; e++) { M m = P.emask(AX[e], c[e], R); a3 |= a2 & m; a2 |= a1 & m; a1 |= m; }
        return !(missP & ~a2) && !(presK & ~a3);
    }
    void dfs(const Cfg& c, M board, int R, int depth) {
        if (aborted) return;
        if (++nodes > nodeLimit) { aborted = true; return; }
        Lvl& lv = L[depth];
        lv.cur = c; lv.board = board; lv.R = R;
        lv.missP = P.Pm & ~board; lv.presK = P.Km & board;
        for (int e = 0; e < 6; e++) {
            const auto& sv = P.succ[c[e]];
            lv.nopt[e] = sv.size();
            for (size_t k = 0; k < sv.size(); k++) {
                lv.opt[e][k] = sv[k];
                lv.om[e][k] = P.emask(AX[e], sv[k], R - 1);
                lv.nm[e][k] = P.nowmask(AX[e], sv[k]);
            }
        }
        lv.s1[6] = lv.s2[6] = lv.s3[6] = 0;
        for (int e = 5; e >= 0; e--) {
            M m = P.emask(AX[e], c[e], R);
            lv.s3[e] = lv.s3[e + 1] | (lv.s2[e + 1] & m);
            lv.s2[e] = lv.s2[e + 1] | (lv.s1[e + 1] & m);
            lv.s1[e] = lv.s1[e + 1] | m;
        }
        choose(depth, 0, 0, 0, 0, 0, 0, 0);
    }
    void choose(int depth, int j, M p1, M p2, M p3, M h1, M h2, M h3) {
        Lvl& lv = L[depth];
        if (aborted) return;
        if (j == 6) { leaf(depth, p2, p3, h2, h3); return; }
        for (int k = 0; k < lv.nopt[j]; k++) {
            // symmetry: identical partner states -> ordered choice
            if (j >= 3 && lv.cur[j] == lv.cur[j - 3] && k < lv.chosenIdx[j - 3]) continue;
            M m = lv.om[j][k];
            M q3 = p3 | (p2 & m), q2 = p2 | (p1 & m), q1 = p1 | m;
            M at2 = q2 | lv.s2[j + 1] | (q1 & lv.s1[j + 1]);
            if (lv.missP & ~at2) continue;
            M at3 = q3 | lv.s3[j + 1] | (q2 & lv.s1[j + 1]) | (q1 & lv.s2[j + 1]);
            if (lv.presK & ~at3) continue;
            M n = lv.nm[j][k];
            M g3 = h3 | (h2 & n), g2 = h2 | (h1 & n), g1 = h1 | n;
            if (lv.R == 1 && (P.Pm & g3)) continue;   // a P cell destroyed at the last tick
            lv.nw[j] = lv.opt[j][k]; lv.chosenIdx[j] = k;
            choose(depth, j + 1, q1, q2, q3, g1, g2, g3);
        }
    }
    void leaf(int depth, M q2, M q3, M h2, M h3) {
        Lvl& lv = L[depth];
        M nb = (lv.board | (h2 & ~h3)) & ~h3;
        stk[depth] = lv.nw;
        bool acc = ((nb & P.Pm) == P.Pm) && !(nb & P.Km);
        if (acc) {
            if (lv.R == 1) {
                u64 k = cfgKey(lv.nw);
                if (!ends.count(k)) {
                    BuildRes br; br.end = lv.nw; br.path.assign(stk.begin(), stk.begin() + len);
                    ends.emplace(k, std::move(br));
                    if (ends.size() >= maxEnds || firstOnly) aborted = true;
                }
            }
            return;
        }
        if (lv.R == 1) return;
        if ((P.Pm & ~nb & ~q2) || (P.Km & nb & ~q3)) return;
        SKey key{cfgKey(lv.nw), nb};
        if (!vis[lv.R - 1].insert(key).second) return;
        Cfg nw = lv.nw;
        dfs(nw, nb, lv.R - 1, depth + 1);
    }
};

// ---------------- reconstruction: config timeline -> instruction lines ----------------
// timeline[0] = start config (all inactive p=0 s=0), timeline[t] = config after tick t.
static Prog reconstruct(const PInfo& P, const vector<Cfg>& tl) {
    int T = tl.size() - 1;
    Prog prog(T, Ins{'-', '-', '-', '-', '-', '-'});
    for (int e = 0; e < 6; e++) {
        int pos = 0;          // true position (line coordinate)
        int offStart = 0;     // tick at which current inactive interval began (toggle-off tick, or 0)
        bool act = false;
        for (int t = 1; t <= T; t++) {
            int cd = tl[t][e];
            bool na = P.cAct[cd];
            if (act && na) {
                int c = P.cC[cd]; int d = c - pos;
                if (d != 0) { int dk = SG[e] * d; prog[t - 1][e] = dk > 0 ? 'F' : 'B'; }
                pos = c;
            } else if (act && !na) {
                prog[t - 1][e] = 'T'; act = false; offStart = t;
            } else if (!act && na) {
                int c = P.cC[cd]; int d = c - pos; int n = abs(d);
                // moves at ticks t-n .. t-1 (1-based), must be > offStart
                if (t - n <= offStart) { fprintf(stderr, "reconstruct: slack violation e=%d t=%d\n", e, t); }
                for (int q = 0; q < n; q++) { int dk = SG[e] * (d > 0 ? 1 : -1); prog[t - 1 - n + q][e] = dk > 0 ? 'F' : 'B'; }
                prog[t - 1][e] = 'T'; act = true; pos = c;
            }
        }
    }
    return prog;
}

static Cfg startCfgP(const PInfo& P) { Cfg c; for (int e = 0; e < 6; e++) c[e] = P.code(false, 0, 0, 0); return c; }

// static config existence (S ⊇ P, S ∩ K = ∅) -> min build length lower bound 1 else 2
static bool hasStatic(const PInfo& P) {
    int NL = P.NL;
    vector<pair<M, M>> opts[3];   // (one, two) masks per axis option
    for (int a = 0; a < 3; a++) {
        opts[a].push_back({0, 0});
        for (int x = 0; x < NL; x++) {
            opts[a].push_back({P.lineM[a][x], 0});
            opts[a].push_back({0, P.lineM[a][x]});
            for (int y = x + 1; y < NL; y++) opts[a].push_back({P.lineM[a][x] | P.lineM[a][y], 0});
        }
    }
    for (auto& A : opts[0]) for (auto& B : opts[1]) {
        M aa = A.first | A.second, bb = B.first | B.second;
        for (auto& C : opts[2]) {
            M cc = C.first | C.second;
            M eq2 = (A.second & ~bb & ~cc) | (B.second & ~aa & ~cc) | (C.second & ~aa & ~bb) |
                    (A.first & B.first & ~cc) | (A.first & C.first & ~bb) | (B.first & C.first & ~aa);
            if ((eq2 & P.Pm) == P.Pm && !(eq2 & P.Km)) return true;
        }
    }
    return false;
}
