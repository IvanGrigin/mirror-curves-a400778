"""Exact checks for the closed form of OEIS A400778 (single-loop mirror curves).

Standard Python 3 only.  Run: python3 verify.py [b400778.txt]
"""
import sys
from fractions import Fraction
from itertools import product


def det(m):
    """Exact determinant by fraction-free (Bareiss) elimination."""
    m = [row[:] for row in m]
    n, sign, prev = len(m), 1, 1
    for c in range(n):
        p = next((r for r in range(c, n) if m[r][c]), None)
        if p is None:
            return 0
        if p != c:
            m[c], m[p], sign = m[p], m[c], -sign
        for r in range(c + 1, n):
            for j in range(c + 1, n):
                m[r][j] = (m[r][j] * m[c][c] - m[r][c] * m[c][j]) // prev
        prev = m[c][c]
    return sign * prev if n else 1


def face_det(n, k):
    """det of the (n-1)(k-1) face matrix M: 4 on the diagonal, -2 towards the
    vertical neighbours of even faces and the horizontal neighbours of odd faces."""
    faces = [(a, b) for a in range(n - 1) for b in range(k - 1)]
    idx = {f: i for i, f in enumerate(faces)}
    m = [[0] * len(faces) for _ in faces]
    for (a, b), i in idx.items():
        m[i][i] = 4
        step = [(0, 1), (0, -1)] if (a + b) % 2 == 0 else [(1, 0), (-1, 0)]
        for da, db in step:
            j = idx.get((a + da, b + db))
            if j is not None:
                m[i][j] = -2
    return det(m)


def closed_form(n, k):
    """2^(F + even - odd) * det(4I - A_d) on the odd faces, A_d = diagonal adjacency.

    This is the exact integer form of
    4^F * prod_{i<n/2, j<k/2} (1 - cos^2(i pi/n) cos^2(j pi/k)),  F = (n-1)(k-1)."""
    odd = [(a, b) for a in range(n - 1) for b in range(k - 1) if (a + b) % 2]
    even = (n - 1) * (k - 1) - len(odd)
    idx = {f: i for i, f in enumerate(odd)}
    m = [[0] * len(odd) for _ in odd]
    for (a, b), i in idx.items():
        m[i][i] = 4
        for da, db in product((1, -1), repeat=2):
            j = idx.get((a + da, b + db))
            if j is not None:
                m[i][j] = -1
    return 2 ** ((n - 1) * (k - 1) + even - len(odd)) * det(m)


def trace_loops(n, k, mirrors):
    """Number of closed paths of the mirror curve, straight from the definition.

    Squares are (r, c).  A ray state is (square, side it is heading to, turn),
    sides 0..3 = E, N, W, S; inside a square the ray goes from the midpoint of
    one side to the midpoint of a neighbouring side."""
    def blocked(r, c, side):
        dr, dc = ((0, 1), (-1, 0), (0, -1), (1, 0))[side]
        r2, c2 = r + dr, c + dc
        if not (0 <= r2 < n and 0 <= c2 < k):
            return True
        return frozenset(((r, c), (r2, c2))) in mirrors

    # Undirected segments: (square, i) joins side i and side i+1 of the square.
    parent = {}

    def find(x):
        while parent.setdefault(x, x) != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for r in range(n):
        for c in range(k):
            for side in range(4):
                a = (r, c, (side - 1) % 4)      # segment ending at `side`
                b = (r, c, side)                # segment starting at `side`
                if blocked(r, c, side):
                    parent[find(a)] = find(b)   # reflect: stay in the square
                else:                           # pass straight through
                    dr, dc = ((0, 1), (-1, 0), (0, -1), (1, 0))[side]
                    opp = (side + 2) % 4
                    # going straight on, the ray enters the next square through
                    # its side `opp` along the segment joining opp-1 and opp
                    parent[find(a)] = find((r + dr, c + dc, (opp - 1) % 4))
                    parent[find(b)] = find((r + dr, c + dc, opp))
    return len({find(x) for x in list(parent)})


def brute(n, k):
    edges = [frozenset(((r, c), (r, c + 1))) for r in range(n) for c in range(k - 1)]
    edges += [frozenset(((r, c), (r + 1, c))) for r in range(n - 1) for c in range(k)]
    total = 0
    for mask in range(1 << len(edges)):
        mirrors = {e for i, e in enumerate(edges) if mask >> i & 1}
        total += trace_loops(n, k, mirrors) == 1
    return total


def grid(n, k):
    """Vertices, edges (tail, head) and the matrices 4*Omega and D of the proof.

    Horizontal edges point east, vertical edges point north (towards row r-1)."""
    verts = [(r, c) for r in range(n) for c in range(k)]
    hor = [((r, c), (r, c + 1)) for r in range(n) for c in range(k - 1)]
    ver = [((r + 1, c), (r, c)) for r in range(n - 1) for c in range(k)]
    edges = hor + ver
    om = [[0] * len(edges) for _ in edges]
    for i, h in enumerate(hor):
        for j, v in enumerate(ver, len(hor)):
            for u in set(h) & set(v):
                om[i][j] = (-1) ** sum(u)       # 4 * Omega(h, v) = eps_u
                om[j][i] = -om[i][j]
    inc = [[(v == e[1]) - (v == e[0]) for e in edges] for v in verts[1:]]
    return edges, om, inc


def bordered(om, inc, keep, diag=0):
    """det [[4*Omega_S + diag*I, D_S^T], [-D_S, 0]] for the edge subset `keep`."""
    top = [[om[i][j] + diag * (i == j) for j in keep] + [row[i] for row in inc] for i in keep]
    bottom = [[-row[j] for j in keep] + [0] * len(inc) for row in inc]
    return det(top + bottom)


def check_lemma(n, k, subsets):
    """Single loop  <=>  the bordered determinant equals 4^(|S|-|V|+1), else it is 0."""
    edges, om, inc = grid(n, k)
    for keep in subsets:
        mirrors = {frozenset(edges[i]) for i in range(len(edges)) if i not in keep}
        single = trace_loops(n, k, mirrors) == 1
        want = 4 ** (len(keep) - len(inc)) if single else 0
        assert bordered(om, inc, keep) == want, (n, k, keep)


if __name__ == '__main__':
    import random
    random.seed(1)
    # 0. the determinant criterion of the proof, against the definition
    for n, k in [(2, 2), (2, 3), (2, 4), (3, 3)]:
        m = len(grid(n, k)[0])
        check_lemma(n, k, ([i for i in range(m) if mask >> i & 1] for mask in range(1 << m)))
        print('determinant criterion, all edge subsets', n, k)
    for n, k in [(3, 4), (4, 4), (3, 6), (5, 5)]:
        m = len(grid(n, k)[0])
        # bias towards few mirrors, where single loops actually occur
        check_lemma(n, k, ([i for i in range(m) if random.random() < p] for p in (0.6, 0.75, 0.9) for _ in range(150)))
        print('determinant criterion, 450 random edge subsets', n, k)
    for n, k in [(2, 2), (3, 3), (3, 4), (4, 4), (4, 5)]:
        edges, om, inc = grid(n, k)
        assert bordered(om, inc, list(range(len(edges))), diag=4) == 4 ** (len(edges) - len(inc)) * face_det(n, k)
    print('summed determinant equals the face determinant')
    # 1. brute force from the definition on small grids
    for n, k in [(1, 4), (2, 2), (2, 3), (2, 4), (3, 3), (2, 5), (3, 4)]:
        b, f = brute(n, k), face_det(n, k)
        assert b == f == closed_form(n, k), (n, k, b, f, closed_form(n, k))
        print('brute force', n, k, b)
    # 2. the published table (antidiagonals), if given
    if len(sys.argv) > 1:
        terms = [int(line.split()[1]) for line in open(sys.argv[1]) if line.strip() and line[0] != '#']
        pos = 0
        d = 1
        while pos < len(terms):
            for n in range(1, d + 1):
                k = d + 1 - n
                assert terms[pos] == face_det(n, k) == closed_form(n, k), (n, k)
                pos += 1
            d += 1
        print('all', pos, 'published terms agree (antidiagonals 1..%d)' % (d - 1))
    # 3. the product formula in floating point, and the statements in the OEIS entry
    from math import cos, pi, isclose
    for n in range(1, 10):
        for k in range(1, 10):
            p = 4.0 ** ((n - 1) * (k - 1))
            for i in range(1, (n + 1) // 2):
                for j in range(1, (k + 1) // 2):
                    p *= 1 - cos(i * pi / n) ** 2 * cos(j * pi / k) ** 2
            t = closed_form(n, k)
            assert isclose(p, t, rel_tol=1e-9), (n, k)
            assert (t == 4 ** ((n - 1) * (k - 1))) == (min(n, k) <= 2) and t <= 4 ** ((n - 1) * (k - 1))
    rec = {2: [4], 3: [16, -16], 4: [64, -512], 5: [256, -11776, 65536, -65536],
           6: [1024, -237568, 12582912, -150994944]}
    for n, coef in rec.items():
        row = [0] + [closed_form(n, k) for k in range(1, 15)]       # T(n,0) = 0
        start = 2 if n == 2 else len(coef)
        for k in range(start, 15):
            assert row[k] == sum(c * row[k - 1 - d] for d, c in enumerate(coef)), (n, k)
    print('product formula, inequality and the five row recurrences hold')
