# The dual constraint hierarchy is uniformly local on the full edgewise tower

Date: 2026-08-12

Protocol commit: `8d0c557`

Finite precursor commit: `9874e6b`

Registered verifier:
`reproducible/verify_whitney_dual_resolution_all_k.py`

Targeted result: **15/15 PASS**.  The full suite was not run, by explicit
user request.

## Verdict

Let

\[
 B=\operatorname{sd}\partial\Delta^4,
 \qquad K_q=\operatorname{Esd}_q(B),\qquad q\geq1,
\]

with the rank-compatible edgewise construction already certified in the
repository.  Then the sharp incidence maxima of the complete dual-cell
constraint resolution are

\[
 \boxed{a_0(q)=24,\qquad a_1(q)=6,\qquad r_3(q)=14}
 \qquad\text{for every integer }q\geq1.
\]

Here `a0` is the maximum number of tetrahedra containing a vertex, `a1` is
the maximum number containing an edge, and `r3` is the maximum number of
edges incident at a vertex.

> **DERIVED UNIFORM LOCALITY**, under the complete hypotheses below.  The
> earlier finite label `PATTERN TOWARD BOUNDED DUAL LOCALITY` is upgraded to
> an all-resolution theorem.

The proof is not a polynomial fit to `q=1,2,4`.  It reduces every possible
local link to a finite partition type and then proves exactly how that local
type glues through the fixed barycentric carrier.

## Complete hypotheses

1. `B` is the order complex of the nonempty proper subsets of a five-element
   set, so its chambers are maximal rank chains `1<2<3<4`.
2. `Esd_q` is the Edelsbrunner--Grayson edgewise subdivision used by the
   repository.  It is a genuine simplicial subdivision and agrees on common
   faces.
3. The fixed-dimension vertex- and face-link classification in Jojić and
   Papaz applies chamberwise.  Their simplex parameter is `4`, while the
   repository's resolution is their `q`.
4. Global cells are glued only through the exact shared-face restriction;
   this is checked parent by parent in the repository's coordinate
   convention.

Primary sources:

- Herbert Edelsbrunner and Daniel R. Grayson,
  [*Edgewise Subdivision of a Simplex*](https://doi.org/10.1007/s004540010063),
  Discrete & Computational Geometry 24 (2000), 707--719.
- Duško Jojić and Ognjen Papaz,
  [*Shelling of links and star clusters in edgewise subdivision of a
  simplex*](https://arxiv.org/abs/2408.12756), arXiv:2408.12756v3 (2025),
  especially the product-of-chains vertex theorem and the face-link join
  theorem.

## 1. Exact vertex incidence

A face of one chamber of `B` is encoded by a nonempty rank chain

\[
 S=\{r_1<\cdots<r_s\}\subseteq\{1,2,3,4\}.
\]

The number of barycentric chambers containing a fixed such face is the
number of permutations of five elements that contain the prescribed initial
sets:

\[
 m_B(S)=r_1!\,(r_2-r_1)!\cdots(r_s-r_{s-1})!\,(5-r_s)!.
\]

Every edgewise vertex in the relative interior of that face has cyclic-gap
partition

\[
 \lambda(S)=
 (r_2-r_1,\ldots,r_s-r_{s-1},4-r_s+r_1).
\]

The product-of-chains link theorem identifies the local incident child
tetrahedra with the shuffles of those blocks.  Therefore

\[
 m_{\rm loc}(S)=\frac{4!}{\prod_i\lambda_i!}.
\]

Top-dimensional children belonging to different base chambers have
disjoint interiors.  The global count is consequently the product

\[
 \begin{aligned}
 a_0(S)
 &=m_B(S)m_{\rm loc}(S)\\
 &=24\frac{r_1!(5-r_s)!}{(4-r_s+r_1)!}
 \in\{12,16,24\}.
 \end{aligned}
\]

The verifier exhausts all fifteen nonempty rank chains.  The only value
`16` comes from `S={2,3}`; the value `12` comes from old rank-two and
rank-three vertices; every other type has value `24`.

This proves `a0(q)<=24` independently of `q`.  Equality already occurs at
every resolution at an old rank-one or rank-four vertex, hence
`a0(q)=24`.

## 2. Exact edge incidence

For an edge in a subdivided tetrahedron, the face-link theorem writes its
local link as

\[
 K_{\sigma_1}*K_{\sigma_2},
\]

where `(lambda_1,lambda_2)` partitions four and each `sigma_i` partitions
`lambda_i`.  The number of local incident tetrahedra is

\[
 \prod_{i=1}^{2}
 \frac{\lambda_i!}{\prod_j\sigma_{i,j}!}.
\]

Exhausting the partitions of four gives exactly five relative isomorphism
types:

| number of vertices in the minimal base carrier | local tetrahedra |
|---:|---:|
| 2 | `1` |
| 3 | `2` or `3` |
| 4 | `4` or `6` |

The corresponding numbers of base chambers are:

- carrier size two: `4` or `6`;
- carrier size three: exactly `2`;
- carrier size four: exactly `1`.

Thus every possible global product is

\[
 4,\quad6,\quad2\cdot2=4,\quad2\cdot3=6,
 \quad1\cdot4=4,\quad1\cdot6=6.
\]

Because `K_q` is a subdivision of the closed combinatorial three-manifold
`B`, every edge link is a circle.  It is therefore exactly `C4` or `C6`.
Both occur at every resolution, so

\[
 a_1(q)=6.
\]

## 3. Vertex degree

Every vertex link is a triangulated two-sphere.  If `F=a0(v)` is its number
of triangular faces, then

\[
 3F=2E,\qquad V-E+F=2,
\]

and hence

\[
 r_3(v)=V=\frac{F}{2}+2.
\]

For `F in {12,16,24}` this gives degrees `{8,10,14}`.  Since `F=24`
occurs, the sharp result is

\[
 r_3(q)=14.
\]

## 4. Why four finite levels are exhaustive rather than an extrapolation

The partition data above always have total size four.  The link
classification proves that no new local combinatorial type can appear after
the resolution has enough coordinate levels to realize four parts.  The
`q=4` standard tetrahedron realizes exactly all five vertex types and all
five relative edge types.

The remaining possible loophole was global gluing.  The verifier retains,
for every refined vertex and edge:

- its exact sparse barycentric coordinates;
- its minimal base rank chain;
- every parent barycentric chamber; and
- the number of incident child tetrahedra contributed by each parent.

At `q=1,2,3,4`, the observed parent set is exactly the set of chambers
containing the base carrier.  Every parent contributes the same
partition-predicted local count.  This is the coordinate-level certificate
for the relative-gluing lemma.  For arbitrary `q`, the same conclusion
follows because the edgewise face restriction and the partition type, not
the numerical coordinate values, determine the local link.

Therefore `q<=4` is a finite set of universal type representatives supplied
by the classification theorem.  It is not being used as a numerical sample
from which the result is guessed.

## 5. Exact all-resolution census

The number of edgewise vertices in the relative interior of an
`(s-1)`-face is `binom(q-1,s-1)`.  Summing over the f-vector
`(30,150,240,120)` of `B`, and then using closed-manifold incidence and
Euler identities, gives

\[
 f(K_q)=
 (20q^3+10q,\;140q^3+10q,\;240q^3,\;120q^3).
\]

The complete vertex-incidence histogram is

\[
 \begin{array}{c|c}
 a_0(v)&\#v\\ \hline
 12&20\\
 16&30(q-1)\\
 24&20q^3-20q+10.
 \end{array}
\]

The complete edge-incidence histogram follows from
`E4+E6=f1` and `4 E4+6 E6=6 f3`:

\[
 \begin{array}{c|c}
 a_1(e)&\#e\\ \hline
 4&60q^3+30q\\
 6&80q^3-20q.
 \end{array}
\]

These are identities for every positive integer `q`, not fitted
polynomials.

## 6. Preregistered new control

Resolution `q=3` had not been evaluated before protocol commit `8d0c557`.
The direct exact enumeration gives

\[
 f(K_3)=(570,3810,6480,3240),
\]

with histograms

```text
tetrahedra per vertex: 12:20, 16:60, 24:490
tetrahedra per edge:    4:1710, 6:2100
edges per vertex:       8:20, 10:60, 14:490
```

It agrees exactly with the analytic certificate.  This is a falsification
control; it is not the proof of the infinite statement.

## 7. Attack on the framing and physical scope

The earlier wording risked treating “finitely many simplex link types” as
automatically proving global boundedness.  That implication is false in
general: gluing arbitrarily many parent cells around a common face can make
the global link arbitrarily large.  The present theorem works only because
the coarse carrier is fixed and its exact chamber multiplicities are the
factorial values above.  Both the local classification and the global
coface counts are load-bearing hypotheses.

The result closes a real mathematical concern:

> the complete canonical dual-cell resolution of the redundant broken-FEEC
> constraints has a uniformly bounded incidence stencil on the entire
> selected rank-edgewise tower.

It does **not** convert the constraints into physical gauge symmetry.  They
remain second class on the original phase space.  No symplectic auxiliary
sector, BRST charge, positive Hamiltonian, causal cone, clock rate, mass or
Planck scale is selected by this theorem.  The known minimal first-class
conversion still places a global Gram inverse in dressed observables.

Thus the correct label is **DERIVED UNIFORM KINEMATIC LOCALITY / OPEN
PHYSICAL DYNAMICS**.

## Status ledger

- **DERIVED:** the exact base chamber multiplicity formula for all fifteen
  rank carriers.
- **DERIVED FROM THE CITED LINK THEOREM:** exactly five local vertex types and
  five relative edge types exhaust every resolution.
- **DERIVED:** exact relative gluing multiplies local incidence by the fixed
  base coface count.
- **DERIVED UNIFORM:** `a0(q)=24`, `a1(q)=6`, `r3(q)=14` for all `q>=1`.
- **DERIVED:** the all-`q` f-vector and both complete incidence histograms.
- **DERIVED CONTROL:** the previously untested `q=3` census agrees exactly.
- **CLOSED:** growth of the complete dual-relation stencil on this selected
  tower.
- **OPEN:** a canonical local first-class/BRST dynamics without a dense Gram
  inverse or an arbitrary new scale.
- **NOT CLAIMED:** physical gauge symmetry, time, causality, inertia, mass,
  matter or Planck units.

## Reproduction

```bash
/home/razvan/science/.venv/bin/python -u \
  reproducible/verify_whitney_dual_resolution_all_k.py
```

Expected result: `15/15`.

## Subsequent Poisson--BRST result

The uniformly local dual hierarchy is used in
`whitney_reducible_poisson_brst_result.md`.  It supplies a canonical
reducible nilpotent Poisson--BRST differential with an exact physical
quotient and no spanning-tree choice.  The advance is kinematic: the local
symplectic-realisation and Hamiltonian-selection gates remain open.
