// Probe: for each pattern, find a static configuration (<=2 lines per axis, duplicates allowed)
// whose exactly-2-hit set contains every pattern cell and no skirt cell. Report min max|offset|.
#include <bits/stdc++.h>
using namespace std;
int main(int argc, char** argv) {
    ifstream in(argv[1]);
    int n; in >> n;
    int found = 0; long long sumCostStatic = 0, sumCount = 0;
    for (int ci = 0; ci < n; ci++) {
        int r, rows, cols, cnt; in >> r >> rows >> cols >> cnt;
        vector<string> g(rows); for (auto& s : g) in >> s;
        // cells: (u,v) with |u|,|v|,|v-u| <= r-1
        set<pair<int,int>> P, K;
        for (int y = 0; y < rows; y++) for (int x = 0; x < cols; x++) if (g[y][x] == 'M') {
            int u = (x - y - r) / 2, v = (x + y - 3 * r) / 2; P.insert({u, v});
        }
        int du[6] = {1, -1, 0, 0, 1, -1}, dv[6] = {1, -1, 1, -1, 0, 0};
        for (auto [u, v] : P) for (int d = 0; d < 6; d++) {
            int a = u + du[d], b = v + dv[d];
            if (abs(a) <= r - 1 && abs(b) <= r - 1 && abs(b - a) <= r - 1 && !P.count({a, b})) K.insert({a, b});
        }
        vector<pair<int,int>> cells(P.begin(), P.end()); int np = cells.size();
        for (auto c : K) cells.push_back(c);
        // axis options: list of lines (0,1,2 entries; duplicates allowed)
        vector<vector<int>> opts; opts.push_back({});
        for (int a = -(r - 1); a <= r - 1; a++) { opts.push_back({a}); opts.push_back({a, a});
            for (int b = a + 1; b <= r - 1; b++) opts.push_back({a, b}); }
        auto mx = [](const vector<int>& o) { int m = 0; for (int x : o) m = max(m, abs(x)); return m; };
        int best = INT_MAX;
        for (auto& A : opts) { int mA = mx(A); if (mA >= best) continue;
            for (auto& B : opts) { int mB = max(mA, mx(B)); if (mB >= best) continue;
                for (auto& C : opts) { int m = max(mB, mx(C)); if (m >= best) continue;
                    bool ok = true;
                    for (int i = 0; i < (int)cells.size() && ok; i++) {
                        auto [u, v] = cells[i]; int h = v - u, hit = 0;
                        for (int x : A) hit += (x == h); for (int x : B) hit += (x == u); for (int x : C) hit += (x == v);
                        if (i < np) ok = (hit == 2); else ok = (hit != 2);
                    }
                    if (ok) best = m;
                } } }
        sumCount += cnt;
        if (best != INT_MAX) { found++; sumCostStatic += best + cnt; printf("case %d r=%d |P|=%d cnt=%d static maxoff=%d cost=%d\n", ci, r, np, cnt, best, best + cnt); }
    }
    printf("static-able: %d / %d\n", found, n);
}
