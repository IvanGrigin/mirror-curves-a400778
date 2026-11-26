# A closed form for single-loop mirror curves (OEIS A400778)

Put a dot in the centre of every square of an n-by-k rectangle. Each side shared by two squares may carry a two-sided mirror. A ray of light starts at the midpoint of a side at 45 degrees, reflects off the border and off every mirror, and passes through sides without a mirror. The rays form the *mirror curve*, a set of closed paths (Gerdes; Jablan, Radović, Sazdanović and Zeković). Let T(n,k) be the number of mirror placements whose mirror curve is a single closed path. This is OEIS A400778, the count of kolam, rangoli or sona designs on the grid that can be drawn in one stroke. [OEIS A400778](https://oeis.org/A400778)

**Theorem.** For all n, k ≥ 1,

$$
T(n,k)=4^{(n-1)(k-1)}\prod_{1\le i<n/2}\ \prod_{1\le j<k/2}\Bigl(1-\cos^2\frac{i\pi}{n}\,\cos^2\frac{j\pi}{k}\Bigr).
$$

The entry lists six conjectures, checked there for k ≤ 22 or n + k ≤ 17. All six follow.

**Corollary 1.** T(n,k) ≤ 4^((n−1)(k−1)), with equality exactly when n ≤ 2 or k ≤ 2. Every factor of the product lies strictly between 0 and 1, and the product is empty exactly in those cases.

**Corollary 2.** Fix n and write $s_i=\sin(i\pi/n)$. Then

$$
T(n,k)=4^{(k-1)[n\ \mathrm{even}]}\prod_{1\le i<n/2}u_i(k),\qquad
u_i(k)=\frac{(8+8s_i)^k-(8-8s_i)^k}{16\,s_i}.
$$

Each $u_i$ satisfies $u_i(k)=16u_i(k-1)-64\cos^2(i\pi/n)\,u_i(k-2)$ with $u_i(0)=0$, $u_i(1)=1$. So row n satisfies a linear recurrence of order $2^{\lfloor (n-1)/2\rfloor}$ whose characteristic roots are the products $4^{[n\ \mathrm{even}]}\prod_i(8\pm 8s_i)$. Explicitly:

| n | Result | Roots |
|---|---|---|
| 2 | T(2,k) = 4^(k−1) | 4 |
| 3 | T(3,k) = 16T(3,k−1) − 16T(3,k−2), which is 4^(k−1)·A001353(k) | 8 ± 4√3 |
| 4 | T(4,k) = 64T(4,k−1) − 512T(4,k−2) | 32 ± 16√2 |
| 5 | T(5,k) = 256T(5,k−1) − 11776T(5,k−2) + 65536T(5,k−3) − 65536T(5,k−4) | 64(1 ± sin 36°)(1 ± sin 72°) |
| 6 | T(6,k) = 1024T(6,k−1) − 237568T(6,k−2) + 12582912T(6,k−3) − 150994944T(6,k−4) | 256(1 ± 1/2)(1 ± √3/2) |

For n ≥ 3, with the convention T(n,0) = 0, each recurrence holds as soon as all its indices are non-negative, which is a slightly larger range than the entry conjectures. For n = 2 the recurrence T(2,k) = 4T(2,k−1) holds for k ≥ 2, starting from T(2,1) = 1.

*Derivation of Corollary 2 from the theorem.* For $c=\cos(i\pi/n)\ne0$, pairing j with k − j gives $\prod_{1\le j<k/2}(1-c^2\cos^2(j\pi/k))=\prod_{j=1}^{k-1}(1-c\cos(j\pi/k))=(c/2)^{k-1}U_{k-1}(1/c)$, where U is the Chebyshev polynomial of the second kind. Multiplying by $16^{k-1}$ gives $u_i(k)=(8c)^{k-1}U_{k-1}(1/c)$, and the three-term recurrence of U gives the recurrence and the closed form of $u_i$. For n = 5 the two sequences have $64\cos^2$ equal to $8(3\pm\sqrt5)$, with sum 48 and product 256. The product of two such sequences has characteristic polynomial $x^4-256x^3+(256\cdot48-2\cdot256)x^2-256\cdot256\,x+256^2$. For n = 6 the values are 48 and 16, and the extra factor 4^(k−1) rescales the roots by 4.

## Proof of the theorem

Let G be the grid graph: the dots are its vertices V, the sides shared by two squares are its edges E. Write F = (n−1)(k−1) for the number of unit faces of G. For a set of mirrors, let S be the set of edges *without* a mirror.

### Step 1. The mirror curve bounds a surface

Inside a square the ray runs between midpoints of adjacent sides, so the four possible segments form a diamond around the dot. At a side with a mirror, or on the border, the curve stays in its square. At a side without a mirror the two diamonds cross.

Take the diamond around each dot as a disc. For every edge in S join the two discs by a band with a half twist. The boundary of the resulting surface Σ_S follows each diamond and, at every band, crosses to the other side. That is the mirror curve. So the number of closed paths equals the number of boundary circles of Σ_S.

Colour the dots like a chessboard, with sign ε_v = ±1. Every edge joins dots of opposite colour. Turning over the disc of every odd dot removes all the half twists. Hence Σ_S is orientable. It is the ribbon graph with untwisted bands in which the edges at an even dot appear in the order east, north, west, south and at an odd dot in the reverse order.

### Step 2. A determinant that detects a single loop

Orient horizontal edges east and vertical edges north. Define a skew-symmetric matrix Ω indexed by E: if a horizontal edge h and a vertical edge e share the dot v, then Ω(h,e) = ε_v/4 and Ω(e,h) = −ε_v/4. All other entries are zero. Let D be the vertex-edge incidence matrix of G with the row of one dot removed, and let Ω_S, D_S denote the restrictions to S.

**Lemma.** The determinant of the block matrix

$$
N_S=\begin{pmatrix}\Omega_S & D_S^{\mathsf T}\\ -D_S & 0\end{pmatrix}
$$

is 1 if the mirror curve is a single closed path, and 0 otherwise.

*Proof.* If (V,S) is disconnected, the rows of D_S belonging to a component that misses the removed dot add up to zero, so the determinant is 0. The curve then has at least two paths. Assume (V,S) is connected and pick a spanning tree inside S. The tree columns of D_S form a matrix of determinant ±1. So there is an integer matrix X with D_S X = I, and together with a basis Z of the integer cycle lattice of (V,S) it forms a matrix (X | Z) of determinant ±1. In this basis the border becomes (I | 0), and block elimination leaves det N_S = det(Zᵀ Ω_S Z).

It remains to identify Zᵀ Ω_S Z. Two standard facts about a connected oriented ribbon graph are used. First, its first homology is the cycle lattice of the graph, and the intersection number of two cycles x, y is the sum over the vertices of ½ Σ (x_a y_b − x_b y_a), taken over pairs a < b of half-edges at the vertex in rotation order, where x_a is the flow of x leaving along a. Because a cycle has zero net flow at each vertex, the starting point of the order does not matter. Here every vertex has its edges among four directions, and averaging over the four starting points gives exactly the matrix Ω above, with the sign ε_v recording the reversed order at odd dots. Second, on a connected compact oriented surface with b ≥ 1 boundary circles, the classes of the boundary circles span the radical of the intersection form, and the form is unimodular on the quotient. So its determinant is 1 if b = 1 and 0 if b ≥ 2. ∎

The accompanying script tests the Lemma against direct ray tracing for every edge subset of the 2×2, 2×3, 2×4 and 3×3 grids and for random subsets of larger grids.

### Step 3. Summing over all placements

By the Lemma, T(n,k) is the sum of det N_S over all S ⊆ E. A determinant is linear in each diagonal entry, so adding 1 to the diagonal entries indexed by E sums the principal minors over all subsets S:

$$
T(n,k)=\det\begin{pmatrix}I+\Omega & D^{\mathsf T}\\ -D & 0\end{pmatrix}=\det\bigl(K^{\mathsf T}(I+\Omega)K\bigr).
$$

The second equality is the same change of basis as in the Lemma, now for the whole grid. Here K is the E-by-F matrix whose columns are the counterclockwise boundaries of the unit faces. They form a basis of the cycle lattice of a plane graph.

### Step 4. The face matrix

Two face boundaries share an edge, traversed in opposite directions, exactly when the faces are neighbours. So KᵀK = 4I − A_h − A_v, where A_h and A_v are the horizontal and vertical adjacency matrices of the (n−1)-by-(k−1) array of faces.

For Q = KᵀΩK, let f lie east of g, sharing a vertical edge m with lower endpoint v. Four pairs of perpendicular edges of f and g meet at an endpoint of m, and each contributes ε_v/4. So Q(f,g) = ε_v. The same count for vertical neighbours gives the opposite sign pattern. For faces touching only at a corner the two contributions cancel. Colouring the faces like a chessboard with signs σ_f, this says

Q(f,g) = σ_f for horizontal neighbours g, and Q(f,g) = −σ_f for vertical neighbours g.

Therefore M = Kᵀ(I+Ω)K has 4 on the diagonal, −2 from every even face to its vertical neighbours, −2 from every odd face to its horizontal neighbours, and zeros elsewhere. In other words M = 2(2I − B), where B sends an even face to its vertical neighbours and an odd face to its horizontal neighbours.

### Step 5. Evaluating the determinant

B exchanges the two colour classes. Eliminating the even faces gives

$$
\det(2I-B)=2^{\,e-o}\det(4I-A_d^{\rm odd}),
$$

where e and o count the even and odd faces. Indeed, a horizontal step followed by a vertical step leads from an odd face to a diagonal neighbour, each diagonal neighbour is reached exactly once, and so the product of the two blocks of B is the diagonal adjacency matrix restricted to the odd faces. Eliminating the odd faces instead gives $2^{\,o-e}\det(4I-A_d^{\rm even})$. Multiplying the two expressions,

$$
\det(2I-B)^2=\det(4I-A_d),
$$

with A_d the diagonal adjacency matrix of the whole array of faces. This matrix is the Kronecker product of the adjacency matrices of the paths with n−1 and k−1 vertices, so its eigenvalues are 4cos(iπ/n)cos(jπ/k) for 1 ≤ i ≤ n−1 and 1 ≤ j ≤ k−1. All of them are smaller than 4 in absolute value, so the determinants are positive and

$$
T(n,k)=2^{F}\det(2I-B)=4^{F}\prod_{i=1}^{n-1}\prod_{j=1}^{k-1}\Bigl(1-\cos\frac{i\pi}{n}\cos\frac{j\pi}{k}\Bigr)^{1/2}.
$$

Pairing i with n−i and j with k−j turns the product into the one in the theorem. ∎

**Remark.** Step 5 also gives an integer formula: T(n,k) = 2^(F+e−o) det(4I − A_d) on the odd faces. By the matrix-tree theorem this determinant counts the spanning trees of the graph formed by the odd faces with diagonal adjacency, with each face joined to an extra root by as many edges as it has missing diagonal neighbours.

## Verification

Run `python3 verify.py` with standard Python 3 (about ten seconds). To include the comparison with the published table, first download the b-file with `curl -O https://oeis.org/A400778/b400778.txt` and run `python3 verify.py b400778.txt`. The script

1. checks the Lemma against ray tracing, as described above;
2. counts single-loop placements by brute force from the definition for seven small grids;
3. compares both determinant formulas, in exact integer arithmetic, with all 136 terms of the OEIS b-file (antidiagonals 1 to 16, computed by Karan Sharma), when that file is given;
4. checks the product formula numerically for n, k ≤ 9, and checks Corollary 1 and the five recurrences for k ≤ 14.

## Scope

The six conjectures were open in the OEIS entry on 26 November 2026. A literature search for an earlier closed form found none, but it was not exhaustive. The two facts about ribbon graphs in Step 2 are quoted, not reproved here. An internal cross-check by a second AI system (OpenAI Codex) recomputed all 136 published terms and four larger grids, up to 10 by 11, by a different method: by the BEST theorem, T(n,k) is a cofactor of the Laplacian of the directed graph on side midpoints in which neighbouring diamonds turn in opposite directions. The note has not been refereed externally.

Prepared with AI assistance (Anthropic Claude); the cross-check mentioned above was run with OpenAI Codex. The sequence, its definition and the six conjectures are due to Karan Sharma.

## References

- P. Gerdes, [Mirror designs and mirror curves](https://www.mi.sanu.ac.rs/vismath/paulus/pg1.htm), Visual Mathematics.
- S. Jablan, L. Radović, R. Sazdanović and A. Zeković, [Mirror-curves and knot mosaics](https://arxiv.org/abs/1106.3784), arXiv:1106.3784.
- K. Sharma, [OEIS A400778](https://oeis.org/A400778), 2026.
