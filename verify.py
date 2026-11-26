"""Exact checks for the mirror-curve theorems (OEIS A400778 and general cell regions).

Standard Python 3 only.  Run:  python3 verify.py [b400778.txt]
A region is a finite set of unit cells (r, c); r grows downwards.
"""
import random
import sys
from itertools import product
from math import cos, isclose, pi

DIRS = ((0, 1), (-1, 0), (0, -1), (1, 0))            # sides E, N, W, S of a cell


def det(m):
    """Exact integer determinant by fraction-free (Bareiss) elimination."""
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


def edges_of(cells):
    hor = [((r, c), (r, c + 1)) for r, c in sorted(cells) if (r, c + 1) in cells]     # pointing east
    ver = [((r + 1, c), (r, c)) for r, c in sorted(cells) if (r + 1, c) in cells]     # pointing north
    return hor, ver


# ---------------------------------------------------------------- the definition
def loops(cells, mirrors):
    """Number of closed paths of the mirror curve, by following the rays.

    Segment (cell, i) joins the midpoints of sides i and i+1 of the cell."""
    parent = {}

    def find(x):
        while parent.setdefault(x, x) != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for r, c in cells:
        for side, (dr, dc) in enumerate(DIRS):
            a, b = (r, c, (side - 1) % 4), (r, c, side)
            other = (r + dr, c + dc)
            if other not in cells or frozenset(((r, c), other)) in mirrors:
                parent[find(a)] = find(b)                       # reflect, stay in the cell
            else:                                               # go straight on
                opp = (side + 2) % 4
                parent[find(a)] = find(other + ((opp - 1) % 4,))
                parent[find(b)] = find(other + (opp,))
    return len({find(x) for x in list(parent)})


def brute(cells):
    hor, ver = edges_of(cells)
    edges = [frozenset(e) for e in hor + ver]
    return sum(loops(cells, {e for i, e in enumerate(edges) if mask >> i & 1}) == 1
               for mask in range(1 << len(edges)))


# ---------------------------------------------------------------- Theorem 1
def tours(cells):
    """Rooted cofactor of the out-Laplacian of the ray digraph (BEST theorem)."""
    arcs = []
    for r, c in cells:
        ports = [(2 * r + dr, 2 * c + dc) for dr, dc in DIRS]   # E, N, W, S: counterclockwise
        if (r + c) % 2 == 0:
            ports.reverse()
        arcs += zip(ports, ports[1:] + ports[:1])
    idx = {v: i for i, v in enumerate(sorted({v for a in arcs for v in a})[1:])}
    lap = [[0] * len(idx) for _ in idx]
    for a, b in arcs:
        if a in idx:
            lap[idx[a]][idx[a]] += 1
            if b in idx:
                lap[idx[a]][idx[b]] -= 1
    return det(lap)


def bordered(cells):
    """det [[I + Omega, D^T], [-D, 0]], computed as an integer determinant with 4I + 4 Omega."""
    hor, ver = edges_of(cells)
    edges = hor + ver
    m = len(edges)
    verts = sorted(cells)[1:]
    top = [[4 * (i == j) for j in range(m)] for i in range(m)]
    for i, h in enumerate(hor):
        for j, v in enumerate(ver, len(hor)):
            for u in set(h) & set(v):
                top[i][j] = (-1) ** sum(u)
                top[j][i] = -top[i][j]
    inc = [[(v == e[1]) - (v == e[0]) for e in edges] for v in verts]
    mat = [top[i] + [row[i] for row in inc] for i in range(m)]
    mat += [[-x for x in row] + [0] * len(verts) for row in inc]
    value, scale = det(mat), 4 ** (m - len(verts))
    assert value % scale == 0
    return value // scale


# ---------------------------------------------------------------- Theorems 2 and 3
def interior_points(cells):
    """Lattice points whose four surrounding cells all belong to the region."""
    return sorted({(r, c) for r, c in cells
                   if {(r - 1, c - 1), (r - 1, c), (r, c - 1)} <= cells})


def full(cells):
    """True if the cells outside the region are connected by king moves (no holes)."""
    rs, cs = [r for r, _ in cells], [c for _, c in cells]
    box = {(r, c) for r in range(min(rs) - 1, max(rs) + 2) for c in range(min(cs) - 1, max(cs) + 2)} - cells
    todo, seen = [min(box)], {min(box)}
    while todo:
        r, c = todo.pop()
        for dr, dc in product((-1, 0, 1), repeat=2):
            q = (r + dr, c + dc)
            if q in box and q not in seen:
                seen.add(q)
                todo.append(q)
    return seen == box


def stencil(cells):
    """det(4I - 2B): an even point is joined to its horizontal neighbours, an odd one to its vertical ones."""
    pts = interior_points(cells)
    idx = {p: i for i, p in enumerate(pts)}
    m = [[0] * len(pts) for _ in pts]
    for (r, c), i in idx.items():
        m[i][i] = 4
        for q in ([(r, c - 1), (r, c + 1)] if (r + c) % 2 == 0 else [(r - 1, c), (r + 1, c)]):
            if q in idx:
                m[i][idx[q]] = -2
    return det(m)


def forest_form(cells):
    """2^(|P| + even - odd) * det(4I - N) on the odd points, plus the digraph behind N.

    N[p][q] counts the two-step moves p -> q (one vertical step, then one horizontal step) inside P."""
    pts = set(interior_points(cells))
    odd = sorted(p for p in pts if sum(p) % 2)
    idx = {p: i for i, p in enumerate(odd)}
    m = [[4 * (i == j) for j in range(len(odd))] for i in range(len(odd))]
    arcs = {p: [] for p in odd}
    for (r, c), i in idx.items():
        for dr, dc in product((-1, 1), repeat=2):
            if (r + dr, c) in pts and (r + dr, c + dc) in pts:
                m[i][idx[(r + dr, c + dc)]] -= 1
                arcs[(r, c)].append((r + dr, c + dc))
    acyclic, left = True, dict(arcs)
    while left:                                               # peel off points without outgoing arcs
        sinks = [p for p, out in left.items() if not any(q in left for q in out)]
        if not sinks:
            acyclic = False
            break
        for p in sinks:
            del left[p]
    return 2 ** (len(pts) + (len(pts) - len(odd)) - len(odd)) * det(m), acyclic, len(pts)


def product_formula(n, k):
    p = 4.0 ** ((n - 1) * (k - 1))
    for i in range(1, (n + 1) // 2):
        for j in range(1, (k + 1) // 2):
            p *= 1 - cos(i * pi / n) ** 2 * cos(j * pi / k) ** 2
    return p


def rectangle(n, k):
    return {(r, c) for r in range(n) for c in range(k)}


def random_region(rng, height, width, size):
    cells = {(rng.randrange(height), rng.randrange(width))}
    while len(cells) < size:
        r, c = rng.choice(sorted(cells))
        dr, dc = rng.choice(DIRS)
        if 0 <= r + dr < height and 0 <= c + dc < width:
            cells.add((r + dr, c + dc))
    return cells


if __name__ == '__main__':
    rng = random.Random(2026)

    # A. Theorem 1 against the definition: every mirror placement of small regions is traced
    count = holes = 0
    while count < 150:
        cells = random_region(rng, 4, 5, rng.randrange(1, 13))
        if sum(map(len, edges_of(cells))) > 13:
            continue
        t = brute(cells)
        assert t == tours(cells) == bordered(cells), sorted(cells)
        count += 1
        holes += not full(cells)
    ring = rectangle(3, 3) - {(1, 1)}
    assert brute(ring) == tours(ring) == bordered(ring) == 8
    print(f'A. definition = Euler-tour cofactor = bordered determinant on {count} random regions '
          f'({holes} with holes) and the 8-cell ring')

    # B. Theorems 1 and 2 on larger regions (no enumeration)
    count = holes = 0
    for _ in range(400):
        cells = random_region(rng, 7, 7, rng.randrange(8, 40))
        t = tours(cells)
        assert t == bordered(cells)
        if full(cells):
            value, acyclic, npts = forest_form(cells)
            assert t == stencil(cells) == value
            assert t <= 4 ** npts and (t == 4 ** npts) == acyclic
            count += 1
        else:
            holes += 1
    print(f'B. cofactor = bordered determinant on 400 larger regions; point stencil, forest form and '
          f'the bound 4^|P| on the {count} without holes ({holes} have holes)')

    # C. Unicyclic regions: the count is the length of the cycle
    for h, w in [(3, 3), (3, 4), (4, 4), (3, 6), (5, 4)]:
        frame = rectangle(h, w) - {(r, c) for r in range(1, h - 1) for c in range(1, w - 1)}
        tail = frame | {(h, 0), (h + 1, 0)}
        assert tours(frame) == tours(tail) == 2 * (h + w) - 4
    print('C. frames with and without a tail: the count equals the cycle length')

    # D. Rectangles: stencil, product formula, bound, and the statements in the OEIS entry
    table = {}
    for n in range(1, 10):
        for k in range(1, 10):
            t = table[n, k] = stencil(rectangle(n, k))
            assert isclose(product_formula(n, k), t, rel_tol=1e-9)
            assert t <= 4 ** ((n - 1) * (k - 1)) and (t == 4 ** ((n - 1) * (k - 1))) == (min(n, k) <= 2)
            if n * k <= 30:
                assert t == tours(rectangle(n, k)) == bordered(rectangle(n, k))
    for n, k in [(1, 4), (2, 2), (2, 3), (2, 4), (3, 3), (2, 5), (3, 4)]:
        assert brute(rectangle(n, k)) == table[n, k]
    rec = {2: [4], 3: [16, -16], 4: [64, -512], 5: [256, -11776, 65536, -65536],
           6: [1024, -237568, 12582912, -150994944]}
    for n, coef in rec.items():
        row = [0] + [stencil(rectangle(n, k)) for k in range(1, 15)]
        for k in range(2 if n == 2 else len(coef), 15):
            assert row[k] == sum(c * row[k - 1 - d] for d, c in enumerate(coef)), (n, k)
    print('D. rectangles up to 9 x 9: product formula, bound, brute force on seven grids, five row recurrences')

    # E. the published table, if given
    if len(sys.argv) > 1:
        terms = [int(line.split()[1]) for line in open(sys.argv[1]) if line.strip() and line[0] != '#']
        pos, d = 0, 1
        while pos < len(terms):
            for n in range(1, d + 1):
                assert terms[pos] == stencil(rectangle(n, d + 1 - n)), (n, d + 1 - n)
                pos += 1
            d += 1
        print(f'E. all {pos} published terms agree (antidiagonals 1..{d - 1})')
