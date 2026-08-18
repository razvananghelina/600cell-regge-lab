# Protocol: all-resolution locality of the edgewise dual hierarchy

Date: 2026-08-12

This protocol is committed before evaluating the previously untested
resolution `q=3`, before constructing the finite relative-link certificate,
and before claiming an all-`q` theorem.  The values `(24,6,14)` were already
known at `q=1,2,4` in commit `9874e6b`; therefore this is **not** a blind
prediction of those values.  The evidential target is a proof covering every
positive integer `q`, not another finite match.

No spectrum, particle target, coupling, mass, `a1=5`, time scale or Planck
quantity may be consulted in this audit.

## Frozen carrier and notation

Let

\[
 B=\operatorname{sd}\partial\Delta^4,
 \qquad K_q=\operatorname{Esd}_q(B),\qquad q\in\mathbb N_{>0},
\]

where each chamber of `B` is ordered by face rank and is subdivided by the
same Edelsbrunner--Grayson edgewise rule already implemented in
`rank_edgewise_level`.  Shared faces are identified by their exact global
barycentric coordinates.

For a vertex `v` and an edge `e` of `K_q`, define

\[
 a_0(v)=\#\{\text{tetrahedra containing }v\},\qquad
 a_1(e)=\#\{\text{tetrahedra containing }e\},
\]

and

\[
 r_3(v)=\#\{\text{edges containing }v\}.
\]

The global maxima are denoted `a0(q)`, `a1(q)` and `r3(q)`.

## Complete hypotheses

1. `B` is exactly the order complex of all nonempty proper subsets of a
   five-element set.  A chamber is a maximal chain with ranks `1,2,3,4`.
2. `Esd_q` is the exact edgewise triangulation, restricts to the same
   triangulation on every common face, and introduces no cross-chamber
   simplex except through those common restrictions.
3. The vertex- and face-link classification of Jojić and Papaz,
   *Shelling of links and star clusters in edgewise subdivision of a
   simplex*, arXiv:2408.12756v3, applies to each ranked chamber.  In their
   notation the fixed simplex parameter is `k=4` and the resolution is `q`.
   Their Theorems 3.3/3.5 and 4.2 classify local links by partitions whose
   total size is four.
4. A local partition type is not yet a global type.  Promotion to an all-`q`
   statement additionally requires the relative-gluing lemma below; it may
   not be assumed from the cited classification.

## Candidate theorem and analytic vertex formula

For a base face represented by a rank chain

\[
 S=\{r_1<\cdots<r_s\}\subseteq\{1,2,3,4\},
\]

the number of chambers of `B` containing that face is

\[
 m_B(S)=r_1!\,(r_2-r_1)!\cdots(r_s-r_{s-1})!\,(5-r_s)!.
\]

The local edgewise vertex type in the relative interior of that face has
cyclic gap partition

\[
 \lambda(S)=
 (r_2-r_1,\ldots,r_s-r_{s-1},4-r_s+r_1).
\]

The cited product-of-chains description predicts

\[
 m_{\rm loc}(S)=\frac{4!}{\prod_i\lambda_i!}
\]

incident child tetrahedra per containing base chamber.  Hence the frozen
global prediction is

\[
 a_0(S)=m_B(S)m_{\rm loc}(S)
 =24\frac{r_1!(5-r_s)!}{(4-r_s+r_1)!}\in\{12,16,24\}.
\]

This formula must be derived from the actual coordinate convention rather
than accepted because it reproduces the already observed values.

## Relative-gluing lemma to prove

For degrees zero and one, the global link of a refined face is determined,
up to simplicial isomorphism, by:

1. its minimal base rank chain `S`;
2. the partition data of its local link inside one ranked chamber; and
3. the fixed coface incidence of `S` in `B`.

The determination must be independent of the numerical coordinate values
and of `q` once the partition type exists.  All local partition data have
total size four, so every type must occur by `q=4`.  Therefore, if the lemma
is proved, exact enumeration of the global types at `q=1,2,3,4` is exhaustive
for every `q>=1`.

A mere count of finitely many local types is insufficient: different
boundary identifications could produce different global links.  The
verifier must carry the boundary-stratum/rank-chain information through the
gluing.

## Frozen gates

The audit passes only if all of the following hold.

1. **Base flag formula.**  Exhaustive enumeration of every nonempty rank
   chain `S` proves the formula for `m_B(S)`.
2. **Local vertex formula.**  The product-of-chains/permutation certificate
   proves `m_loc(S)=4!/prod(lambda_i!)` in the repository's coordinate
   convention.
3. **Global vertex formula.**  Exact gluing proves
   `a0(v) in {12,16,24}` for every vertex and every `q`, hence `a0(q)=24`.
4. **Finite-type saturation.**  A machine-readable certificate enumerates
   every possible relative vertex and edge link type allowed by the
   partition classification and proves that it is represented among
   `q=1,2,3,4` with the same global gluing data.
5. **Edge bound.**  Every certified global edge link is a cycle of length
   four or six, hence `a1(q)=6`.
6. **Vertex-degree bound.**  Every global vertex link is a triangulated
   two-sphere.  From `F=a0(v)<=24`, `3F=2E` and `V-E+F=2`, prove
   `r3(v)=V=F/2+2<=14`, with equality represented; hence `r3(q)=14`.
7. **Untested control.**  Direct enumeration at `q=3` agrees exactly with
   the analytic/type certificate.  This is a falsification control, not the
   proof of the infinite claim.
8. **Implementation independence.**  The theorem certificate must not infer
   all-`q` coverage by extrapolating finite arrays or fitting a polynomial.

## Acceptance and kill boundaries

**Acceptance:** all eight gates pass.  The finite result from `9874e6b` may
then be upgraded to **DERIVED UNIFORM LOCALITY** under the complete hypotheses
above, with sharp all-resolution bounds

\[
 \boxed{a_0(q)=24,\quad a_1(q)=6,\quad r_3(q)=14}
 \qquad(q\geq1).
\]

**Kill/hold boundary:** if the relative-gluing lemma is false, if a legal
partition/gluing type is absent from `q<=4`, or if any certified type exceeds
the proposed bounds, the sharp theorem is rejected.  A weaker uniform bound
may be reported only as a separate theorem with its full hypotheses; it may
not be presented as proof of the sharp claim.

## Labels fixed in advance

- Exact finite enumerations: **DERIVED ON THE ENUMERATED CONTROLS**.
- All-`q` result after every gate: **DERIVED UNIFORM**, conditional on the
  cited edgewise-link theorems and the proved relative-gluing lemma.
- Finite agreement without the gluing lemma: **PATTERN**, not a theorem.
- Physical interpretation: **STRUCTURAL KINEMATICS ONLY**.
- Time, causality, inertia, mass, gauge dynamics and Planck units: **NOT
  CLAIMED**.
