# Level 1: print a hexagonal board for each radius.
# Radius r -> 2r+1 rows. A row at distance d from the middle row is indented
# with d dashes on each side and has 2r+1-d tokens joined by '-'.
# The top/bottom rows (d == r) are all '*'; other rows are '*' + 'O'... + '*'.
import sys
def board(r):
    out = []
    for i in range(-r, r + 1):
        d = abs(i)
        n = 2 * r + 1 - d
        t = ['*'] * n if d == r else ['*'] + ['O'] * (n - 2) + ['*']
        out.append('-' * d + '-'.join(t) + '-' * d)
    return '\n'.join(out)
for name in sys.argv[1:]:
    nums = open('in_level-1_' + name + '.txt').read().split()
    n = int(nums[0])
    # each board is followed by a mandatory empty line
    res = ''.join(board(int(x)) + '\n\n' for x in nums[1:n + 1])
    open('out_level-1_' + name + '.txt', 'w', newline='\n').write(res)
