
// state of emitters after running prog from start
Cfg cfgAfter(const vector<string>& prog, int r, int upto) {
    Cfg c; for (int e = 0; e < 6; e++) c.k[e] = 0; c.a = 0;
    for (int t = 0; t < upto; t++) for (int e = 0; e < 6; e++) { char ch = prog[t][e];
        if (ch == 'F') c.k[e]++; else if (ch == 'B') c.k[e]--; else if (ch == 'T') c.a ^= 1 << e; }
    (void)r; return c;
}
