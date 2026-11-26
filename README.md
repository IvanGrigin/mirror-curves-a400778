# Single-loop mirror curves on regions of square cells (OEIS A400778)

A *region* R is a finite nonempty set of unit cells of the square grid. Put a dot in the centre of every cell. A side shared by two cells of R is *internal*; each internal side may carry a two-sided mirror. A ray of light starts at the midpoint of a side at 45 degrees. It reflects off every mirror and off the boundary of R, and passes straight through internal sides without a mirror. The rays form the *mirror curve*, a set of closed paths (Gerdes; Jablan, Radović, Sazdanović and Zeković). Let T(R) be the number of mirror placements for which the mirror curve is a single closed path. For the n-by-k rectangle this is T(n,k), OEIS A400778: the number of kolam, rangoli or sona designs on the grid that can be drawn in one stroke. [OEIS A400778](https://oeis.org/A400778)

The *cell graph* G has the cells of R as vertices and the internal sides as edges. Write v and m for their numbers. If G is disconnected, no placement gives a single path, so R is assumed *connected*, meaning that G is connected. Colour the cells like a chessboard: ε_u = +1 if row plus column of u is even, and −1 otherwise.

## Results

**Notation.** Orient the edges of G joining two cells of a row eastwards, and the others northwards. Let D be the vertex-edge incidence matrix of G with the row of one cell deleted. Let Ω be the skew-symmetric m-by-m matrix with Ω(h,e) = ε_u/4 and Ω(e,h) = −ε_u/4 whenever a horizontal edge h and a vertical edge e share the cell u, and zeros elsewhere.

**Theorem 1 (every connected region).** Let K be any matrix whose columns form an integral basis of the cycle lattice of G. Then

```math
T(R)=\det\begin{pmatrix}I+\Omega & D^{\mathsf T}\\ -D & 0\end{pmatrix}=\det\bigl(K^{\mathsf T}(I+\Omega)K\bigr).
```

**Corollary 1.** If G is a tree, T(R) = 1. If G has exactly one cycle, T(R) is the length of that cycle. If two connected regions are joined along a single internal side, the counts multiply.

An *interior point* of R is a grid point all four of whose surrounding cells belong to R. Let P be the set of interior points. Call a point even or odd according to the parity of its row plus its column, and write e and o for the numbers of even and odd points of P. Say that R *has no holes* if the cells outside R are connected by king moves.

**Theorem 2 (regions without holes).** Let B be the 0/1 matrix on P with B(p,q) = 1 when q is a horizontal neighbour of an even point p, or a vertical neighbour of an odd point p. Then

```math
T(R)=\det(4I-2B)=2^{\,|P|+e-o}\,\det(4I-N),
```

where N is the matrix on the odd points with N(p,q) the number of two-step moves from p to q inside P that take one vertical step and then one horizontal step. Consequently T(R) is 2^(|P|+e−o) times the number of spanning trees, directed towards the root, of the following digraph Γ: its vertices are the odd points and a root; each odd point p sends out four arcs, one for each choice of a vertical step followed by a horizontal step, ending at the point reached if the move stays inside P and at the root otherwise.

**Corollary 2.** For a region without holes, T(R) ≤ 4^|P|. Equality holds exactly when Γ has no directed cycle.

**Theorem 3 (rectangles).** For all n, k ≥ 1,

```math
T(n,k)=4^{(n-1)(k-1)}\prod_{1\le i\lt n/2}\ \prod_{1\le j\lt k/2}\Bigl(1-\cos^2\frac{i\pi}{n}\,\cos^2\frac{j\pi}{k}\Bigr).
```

**Corollary 3 (the six conjectures in the OEIS entry).** T(n,k) ≤ 4^((n−1)(k−1)), with equality exactly when n ≤ 2 or k ≤ 2. Moreover, with $s_i=\sin(i\pi/n)$,

```math
T(n,k)=4^{(k-1)[n\ \mathrm{even}]}\prod_{1\le i\lt n/2}u_i(k),\qquad
u_i(k)=\frac{(8+8s_i)^k-(8-8s_i)^k}{16\,s_i},
```

and $u_i(k)=16u_i(k-1)-64\cos^2(i\pi/n)\,u_i(k-2)$ with $u_i(0)=0$, $u_i(1)=1$. So row n satisfies a linear recurrence of order $2^{\lfloor (n-1)/2\rfloor}$:

| n | Result | Characteristic roots |
|---|---|---|
| 2 | T(2,k) = 4^(k−1) | 4 |
| 3 | T(3,k) = 16T(3,k−1) − 16T(3,k−2), which is 4^(k−1)·A001353(k) | 8 ± 4√3 |
| 4 | T(4,k) = 64T(4,k−1) − 512T(4,k−2) | 32 ± 16√2 |
| 5 | T(5,k) = 256T(5,k−1) − 11776T(5,k−2) + 65536T(5,k−3) − 65536T(5,k−4) | 64(1 ± sin 36°)(1 ± sin 72°) |
| 6 | T(6,k) = 1024T(6,k−1) − 237568T(6,k−2) + 12582912T(6,k−3) − 150994944T(6,k−4) | 256(1 ± 1/2)(1 ± √3/2) |

For n ≥ 3, with the convention T(n,0) = 0, each recurrence holds as soon as all its indices are non-negative. For n = 2 it holds for k ≥ 2, starting from T(2,1) = 1. The entry states these recurrences and the bound as conjectures, checked for k ≤ 22 or n + k ≤ 17.

## 1. Rays are Euler tours

Inside a cell a ray runs between the midpoints of two adjacent sides, so the four possible segments form a diamond around the dot. Orient the diamond of u clockwise if ε_u = +1 and counterclockwise otherwise. Let Δ be the directed graph whose vertices are the side midpoints of R and whose arcs are the oriented diamond segments. A boundary midpoint has one incoming and one outgoing arc. An internal midpoint has two of each, one from either cell.

Consider an internal midpoint between a clockwise cell on the left and a counterclockwise cell on the right. The arc entering it from the left cell heads down and to the right, and the arc leaving it into the right cell heads the same way. So "pass straight through" pairs the incoming arc of each cell with the outgoing arc of the other cell, and "reflect" pairs it with the outgoing arc of its own cell. These are the two bijections from the incoming arcs to the outgoing arcs. The other positions follow by symmetry. Hence a mirror placement is the same thing as a choice, at every vertex of Δ, of a bijection from incoming to outgoing arcs, and the closed paths of the mirror curve are the orbits. A placement gives a single path exactly when the orbit is one Euler tour of Δ, and distinct placements give distinct tours.

Δ is connected because G is, and every vertex has equal in- and outdegree. By the BEST theorem the number of Euler tours is the number of spanning trees directed towards a fixed root, times the product of (outdegree − 1)! over the vertices. All outdegrees are 1 or 2. By the matrix-tree theorem for digraphs,

**T(R) = det L_ρ,**

where L is the out-Laplacian of Δ (outdegrees on the diagonal, minus the number of arcs from x to y off the diagonal) and L_ρ deletes the row and column of a boundary midpoint ρ.

## 2. From midpoints to edges: proof of Theorem 1

**Lemma A.** Let A be a square integer matrix of order a, and C an integer a-by-c matrix. Suppose (N | X) is a square integer matrix of determinant ±1 with CᵀN = 0 and det(CᵀX) = ±1. Then

```math
\det\begin{pmatrix}A & C\\ -C^{\mathsf T} & 0\end{pmatrix}=\det(N^{\mathsf T}AN).
```

*Proof.* Multiply the block matrix by diag((N | X), I) on the right and by its transpose on the left. The determinant is unchanged. With G = XᵀC the result has block rows (NᵀAN, ∗, 0), (∗, ∗, G) and (0, −Gᵀ, 0). Since G is invertible, the last block row clears the upper ∗ and the last block column clears the left ∗. What remains is det(NᵀAN) times the determinant of the matrix with block rows (∗, G) and (−Gᵀ, 0), which is det(G)² = 1. ∎

*Ports.* Give every cell four ports, one per side, and let s(a) be the port after a on the oriented diamond. Let L₀ be the out-Laplacian of the disjoint union of the v oriented diamonds, so (L₀z)_a = z_a − z_{s(a)}. For an edge of G, call the port at its tail cell its tail port and the other one its head port. Let C be the matrix, ports by edges, with +1 at the tail port and −1 at the head port of each edge. Take ρ to be a boundary port. Deleting it from L₀ and C gives A and C′.

*First elimination.* Let N be the 0/1 matrix, ports by midpoints other than ρ, recording which midpoint each port lies on, and let X have one column per edge, the indicator of its tail port. Then (N | X) has determinant ±1, C′ᵀN = 0 and C′ᵀX = I. Gluing ports adds the corresponding rows and columns of a Laplacian, so NᵀAN = L_ρ. By Lemma A and Section 1,

```math
T(R)=\det\begin{pmatrix}A & C'\\ -C'^{\mathsf T} & 0\end{pmatrix}.
```

*Second elimination.* In every cell choose a reference port, taking ρ in its own cell. Replace the four potentials of a cell by a constant, the value at the reference port, and three deviations from it. For the cell of ρ the constant is absent. This is a change of basis of determinant ±1; apply it to rows and columns. Rows and columns of each diamond Laplacian sum to zero, so the constants drop out of L₀. On the deviations it becomes the direct sum H of the matrices H_u obtained from the diamond Laplacians by deleting the reference port. Listing the remaining ports 1, 2, 3 along the diamond, H_u is upper triangular with unit diagonal, and its inverse has all entries on and above the diagonal equal to 1. The matrix C′ splits into its part J on the deviations and the reduced incidence matrix D on the constants. With the variables ordered as deviations, constants, edges, the matrix has block rows (H, 0, J), (0, 0, D), (−Jᵀ, −Dᵀ, 0). Eliminating H, which has determinant 1, and reordering gives

```math
T(R)=\det\begin{pmatrix}S & D^{\mathsf T}\\ -D & 0\end{pmatrix},\qquad S=J^{\mathsf T}H^{-1}J .
```

*Third elimination.* A spanning tree of G shows that D has an integral right inverse X, and (K | X) has determinant ±1. Lemma A with C = Dᵀ and N = K gives T(R) = det(KᵀSK).

*The form S on cycles.* Let x and x′ be integer cycles of G. In a cell u they define port currents y and y′: the value at a port is +x_e at the tail port of e, −x_e at its head port, and 0 on the boundary. Conservation means that y and y′ sum to zero. By the shape of H_u⁻¹, the contribution of u to x′ᵀSx is the sum of y′_a y_b over 1 ≤ a ≤ b ≤ 3. For vectors of zero sum this equals y′ᵀz for any solution of z_a − z_{s(a)} = y_a, so it does not depend on the reference port. Averaging the analogous sums over the four possible reference ports and subtracting half of (Σy′)(Σy) = 0 gives

```math
\tfrac12\sum_a y'_a y_a+\tfrac14\sum_a\bigl(y'_a\,y_{s(a)}-y'_{s(a)}\,y_a\bigr).
```

Summed over the cells, the first term is x′ᵀx, because every edge has two ports. In the second term only pairs of perpendicular edges through a common cell survive, and checking the four positions of such a pair gives −ε_u/4 each time, for a horizontal x′-edge and a vertical x-edge. So KᵀSK = Kᵀ(I − Ω)K, the transpose of Kᵀ(I + Ω)K. A last application of Lemma A turns det(Kᵀ(I + Ω)K) into the bordered determinant. This proves Theorem 1. ∎

*Proof of Corollary 1.* For a tree K has no columns. For one cycle K is a single column x with entries ±1 on the cycle, and xᵀ(I + Ω)x = xᵀx because Ω is skew-symmetric. If two regions are joined along one side, that edge is a bridge of G. A basis of the cycle lattice is the union of bases for the two parts, no cell carries cycles of both parts, and Kᵀ(I + Ω)K is block diagonal. ∎

## 3. Regions without holes: proof of Theorem 2

Draw G with straight edges between cell centres. A bounded face containing an interior point is the unit square around it. If some bounded face contained no interior point, it would contain a grid point x next to a cell c outside R. No edge of G meets c or a grid point, so the centre of c lies in the same face, and so does every cell reached from c by king moves outside R. This contradicts the assumption that R has no holes. Hence the bounded faces of G are exactly the unit squares around the points of P, and their clockwise boundaries form an integral basis K of the cycle lattice.

Two of these boundaries share an edge, traversed in opposite directions, exactly when the points are neighbours. So KᵀK = 4I − A_h − A_v, with A_h and A_v the horizontal and vertical adjacency matrices of P.

For Q = KᵀΩK, let q be the east neighbour of p. Their faces share the vertical edge between the two cells that touch both points. Four pairs of perpendicular edges, one from each face, meet at the ends of that edge, and each contributes ε/4 with the same sign. So Q(p,q) = ±1. For vertical neighbours the same count gives the opposite sign pattern, and for points whose faces share only a cell the two contributions cancel. With a chessboard sign σ on P this reads Q(p,q) = σ_p for horizontal neighbours and −σ_p for vertical neighbours. Therefore Kᵀ(I + Ω)K has 4 on the diagonal, −2 from each point of one parity to its vertical neighbours and −2 from each point of the other parity to its horizontal neighbours. Exchanging the two parities transposes the matrix. This is det(4I − 2B).

B maps even points to odd points and back. Eliminating the even points from 2I − B gives det(2I − B) = 2^(e−o) det(4I − N), where N is the product of the two blocks: a vertical step from an odd point followed by a horizontal step. The matrix 4I − N is the out-Laplacian of Γ with the root deleted, and the matrix-tree theorem gives the statement about spanning trees. ∎

*Proof of Corollary 2.* A spanning tree directed towards the root selects one outgoing arc at each odd point, so there are at most 4^o of them, and 2^(|P|+e−o)·4^o = 4^|P|. Every selection is such a tree exactly when Γ has no directed cycle. ∎

## 4. Rectangles: proof of Theorem 3 and Corollary 3

For the n-by-k rectangle, P is an (n−1)-by-(k−1) array and F = (n−1)(k−1). A vertical step followed by a horizontal step reaches every diagonal neighbour exactly once, since the corner between them is again in P. So N is the diagonal adjacency matrix on the odd points. Eliminating the odd points instead gives the same with the even points and the factor 2^(o−e). Multiplying,

```math
\det(2I-B)^2=\det(4I-A_d),
```

with A_d the diagonal adjacency matrix of the whole array. It is the Kronecker product of the adjacency matrices of the paths with n−1 and k−1 vertices, so its eigenvalues are 4cos(iπ/n)cos(jπ/k) for 1 ≤ i ≤ n−1 and 1 ≤ j ≤ k−1, all smaller than 4 in absolute value. Hence

```math
T(n,k)=2^{F}\det(2I-B)=4^{F}\prod_{i=1}^{n-1}\prod_{j=1}^{k-1}\Bigl(1-\cos\frac{i\pi}{n}\cos\frac{j\pi}{k}\Bigr)^{1/2},
```

and pairing i with n−i and j with k−j gives the product in Theorem 3. ∎

*Proof of Corollary 3.* Every factor of the product lies strictly between 0 and 1, and the product is empty exactly when n ≤ 2 or k ≤ 2. For $c=\cos(i\pi/n)\ne0$, pairing j with k − j gives $\prod_{1\le j\lt k/2}(1-c^2\cos^2(j\pi/k))=\prod_{j=1}^{k-1}(1-c\cos(j\pi/k))=(c/2)^{k-1}U_{k-1}(1/c)$, with U the Chebyshev polynomial of the second kind. Multiplying by $16^{k-1}$ gives $u_i(k)=(8c)^{k-1}U_{k-1}(1/c)$, and the three-term recurrence of U gives the recurrence and the closed form of $u_i$. For n = 5 the two values of $64\cos^2$ are $8(3\pm\sqrt5)$, with sum 48 and product 256, and the product of the two sequences has characteristic polynomial $x^4-256x^3+(256\cdot48-2\cdot256)x^2-256\cdot256\,x+256^2$. For n = 6 the values are 48 and 16, and the factor 4^(k−1) rescales the roots by 4. ∎

## Remarks

*Earlier work.* The step from Euler tours to a determinant of the form "identity plus a skew-symmetric intersection matrix" is not new. For digraphs with all in- and outdegrees equal to 2 it is the formula of Macris and Pulé, further studied by Lauri. In the language of embedded graphs, single-loop placements are the spanning quasi-trees of an orientable ribbon graph, and Theorem 1 can be viewed as an instance of the matrix-quasi-tree theorem for such graphs (Merino, Moffatt and Noble); we have not compared the two statements in detail. What this note adds is the explicit evaluation for regions of square cells: the matrix on interior points in Theorem 2, its reading as a count of directed trees, and the product formula for rectangles with its consequences for A400778. We did not find these in the literature; the search was not exhaustive.

*Geometric meaning of Ω.* The discs bounded by the diamonds, joined by a half-twisted band across every side without a mirror, form an orientable surface whose boundary is the mirror curve. Kᵀ Ω K is the matrix of intersection numbers of cycles on that surface. The proof above does not use this.

## Verification

Run `python3 verify.py` with standard Python 3 (about fifteen seconds). To include the published table, download it with `curl -O https://oeis.org/A400778/b400778.txt` and run `python3 verify.py b400778.txt`. In exact integer arithmetic the script checks:

- on 150 random regions and an 8-cell ring, by tracing the rays for every mirror placement, that the count equals the cofactor of Section 1 and the bordered determinant of Theorem 1;
- on 400 larger random regions, that the cofactor equals the bordered determinant, and for those without holes that it equals det(4I − 2B) and the directed-tree form, with the bound and equality case of Corollary 2;
- Corollary 1 on rectangular frames with and without an attached tail;
- for rectangles up to 9 by 9, the product formula, the bound, seven brute-force counts and the five recurrences up to k = 14;
- all 136 terms of the OEIS b-file, antidiagonals 1 to 16, computed by Karan Sharma.

## Credits and scope

The sequence, its definition and the six conjectures are due to Karan Sharma. The note was prepared with AI assistance: the determinant on interior points, the product formula and this write-up with Anthropic Claude; the reduction through Euler tours in Sections 1 and 2 and the count for regions with one cycle with OpenAI Codex, each checked by the other system. The six conjectures were open in the OEIS entry on 26 November 2026. The note has not been refereed.

## References

- T. van Aardenne-Ehrenfest and N. G. de Bruijn, Circuits and trees in oriented linear graphs, Simon Stevin 28 (1951), 203–217.
- P. Gerdes, [Mirror designs and mirror curves](https://www.mi.sanu.ac.rs/vismath/paulus/pg1.htm), Visual Mathematics.
- S. Jablan, L. Radović, R. Sazdanović and A. Zeković, [Mirror-curves and knot mosaics](https://arxiv.org/abs/1106.3784), arXiv:1106.3784.
- J. Lauri, On a formula for the number of Euler trails for a class of digraphs, Discrete Mathematics (1997).
- N. Macris and J. V. Pulé, An alternative formula for the number of Euler trails for a class of digraphs, Discrete Mathematics 154 (1996), 301–305.
- C. Merino, I. Moffatt and S. Noble, [The critical group of a combinatorial map](https://arxiv.org/abs/2308.13342), Combinatorial Theory 5(3) (2025).
- K. Sharma, [OEIS A400778](https://oeis.org/A400778), 2026.
- W. T. Tutte, The dissection of equilateral triangles into equilateral triangles, Proceedings of the Cambridge Philosophical Society 44 (1948), 463–482.
