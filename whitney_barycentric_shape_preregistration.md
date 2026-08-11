# Preregistration: shape regularity of the Whitney barycentric tower

Date: 2026-08-11

## Question

The exact trace-stiffness second-refinement audit found improving degree
balance but accelerating worst-element Dirac norms and 24 element types.
Before interpreting any refinement factor physically, test the more basic
geometric hypothesis:

> Is iterated barycentric subdivision a uniformly shape-regular refinement
> tower for the regular tetrahedral carrier?

If not, stiffness growth may compensate element degeneration rather than
encode a physical renormalization group.

## Exact repeated-flag subroute

Let a parent tetrahedron have ordered vertices
((p_0,p_1,p_2,p_3)) and edge matrix

\[
 A=[p_1-p_0\;p_2-p_0\;p_3-p_0].
\]

Select at every barycentric step the same flag ordering `(0,1,2,3)`.  Its
child vertices are

\[
 p_0,
 \quad\frac{p_0+p_1}{2},
 \quad\frac{p_0+p_1+p_2}{3},
 \quad\frac{p_0+p_1+p_2+p_3}{4}.
\]

Therefore the child edge matrix must be

\[
 A'=AT,
 \qquad
 T=
 \begin{pmatrix}
 1/2&1/3&1/4\\
 0&1/3&1/4\\
 0&0&1/4
 \end{pmatrix}.
\]

The verifier will establish this identity exactly from barycentres before
using (T).

The frozen theorem candidate is:

1. (T) has exact eigenvalues (1/2,1/3,1/4);
2. after (n) repeated selections, (A_n=A_0T^n);
3. (kappa_2(T^n)\geq2^n), because the eigenvalue modulus ratio is (2^n);
4. for the regular reference tetrahedron, (kappa_2(A_0)=2), hence

\[
 \kappa_2(A_n)\geq2^{n-1},
 \qquad
 \kappa_2(G_n)\geq4^{n-1},
 \quad G_n=A_n^*A_n.
\]

If these statements hold, the tower is not uniformly shape regular.  One
nested subroute with unbounded aspect is sufficient; no average can rescue
the uniform condition.

## Exact type enumeration

Independently enumerate every exact Gram type generated from one regular
tetrahedron through four barycentric levels.  For each parent type and all 24
vertex orderings, compute

\[
 G'=T_\pi^*GT_\pi
\]

over the rationals, merge identical matrices exactly, and propagate exact
multiplicities.

Frozen levels: (0,1,2,3,4).  The total multiplicity must be (24^r) at
level (r).  No stochastic sample is allowed.

## Frozen shape diagnostics

For every exact type record:

1. affine condition number (sqrt{\kappa_2(G)});
2. diameter (h), the largest of all six edge lengths;
3. volume;
4. permutation-invariant tetrahedral mean-ratio quality

\[
 q=
 \frac{6(6\sqrt2)^{2/3}V^{2/3}}
 {\sum_{e}\ell_e^2},
\]

normalized so (q=1) for a regular tetrahedron;
5. numerical mass-orthonormal local Kähler--Dirac norm (a(G));
6. dimensionless shape factor (h\,a(G)).

The numerical local Whitney masses are reconstructed directly from (G) and
the defining barycentric moment integral.  They must reproduce the exact
symbolic level-zero and repeated first-child norms within (10^{-11}) before
the enumeration is accepted.

For each level record exact type count, total multiplicity, minimum (q),
maximum affine condition, and maximum (h,a).  No exponent is fitted.

Also record the repeated-flag chain through (n=8), including its exact Gram
matrix, condition number, mean ratio, and normalized Dirac factor.

## Decision protocol

- If the exact repeated-flag recurrence and spectral data hold, label
  **DERIVED NEGATIVE: the barycentric tower is not uniformly shape regular**.
- If the exact recurrence fails, discard the theorem and report the error.
- Numerical degradation of (q) or growth of (h,a) is supporting
  **DERIVED CONTROL**, not the proof itself.
- Do not infer a fitted divergence exponent from levels 0--4 or the chain.

## Scope

This result would kill iterated barycentric subdivision as a uniformly
controlled continuum carrier.  It would not kill:

- the finite base or first-refined Whitney geometry;
- exact trace stiffness on any fixed finite carrier;
- a shape-regular refinement rule;
- adaptive remeshing with a separately derived metric;
- a non-simplicial microscopic dynamics.

It also would not derive physical time, mass, (c), or a Planck scale.

## Status before execution

- **DERIVED INPUT:** exact barycentric child coordinates.
- **DERIVED INPUT:** observed type proliferation and local norm growth through
  level two.
- **OPEN:** exact repeated-flag theorem.
- **OPEN:** complete type/quality census through level four.
- **OPEN:** whether a shape-regular alternative is selected by the theory.
