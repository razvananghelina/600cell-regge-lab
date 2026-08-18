# Protocol: canonical rank-edgewise Whitney refinement

Date: 2026-08-11

This protocol is frozen before running the exact enumeration below.  It tests
a bounded construction; it does not search over spectra, penalty constants, or
physical targets.

## Question and hypotheses

Let `K` be a finite tetrahedral complex whose top cells are regular
tetrahedra, in particular the boundary complex of the 600-cell.  The candidate
tower is

\[
 K_n=\operatorname{Esd}_{2^n}(\operatorname{sd}K),\qquad n\geq 0.
\]

Here `sd` is applied **once**.  A chamber of `sd K` is a complete flag

\[
 \sigma_0\subset\sigma_1\subset\sigma_2\subset\sigma_3,
 \qquad \dim\sigma_i=i,
\]

and its four vertices are ordered by the derived ranks `(0,1,2,3)`.
`Esd_k` is the edgewise subdivision of Edelsbrunner and Grayson with that
order.  No vertex identifiers, longest-edge tie breakers, octahedral
diagonals, or optimized coefficients may enter the construction.

The relevant primary result is H. Edelsbrunner and D. R. Grayson,
“Edgewise Subdivision of a Simplex,” *Discrete & Computational Geometry* 24
(2000), 707–719, DOI `10.1007/s004540010063`.  It proves face compatibility,
composition `Esd_l Esd_k = Esd_lk`, and at most `d!/2` congruence classes in
dimension `d`, independently of `k`.

## Independent exact realization

The verifier will not import a mesh package.  It will enumerate the paper's
full color schemes directly.  For a tetrahedron and integer `k`, a full scheme
is a `k x 4` matrix with entries in `{0,1,2,3}`, nondecreasing in row-major
order, and four distinct columns.  A column represents the exact barycentric
lattice point obtained by counting its four colors and dividing by `k`.

All incidence and containment tests use rational arithmetic.  Floating point
is allowed only for reporting condition numbers and quality diagnostics after
the exact gates pass.

## Preregistered gates

1. **Color-scheme realization.**  For `k=1,2,3,4`, the enumeration must give
   exactly `k^3` top tetrahedra, `binomial(k+3,3)` vertices, disjoint interiors,
   and exact volume `1/k^3` per child.
2. **Uniform shapes.**  For the rank-ordered barycentric chamber of a regular
   tetrahedron, every child edge path must be an exact permutation of the
   three parent shape vectors divided by `k`.  Modulo all vertex relabellings,
   the normalized congruence signatures at `k=2,3,4` must agree and contain no
   more than three classes.  The general proof is then the color-transition
   argument, not extrapolation from three levels.
3. **Nesting.**  Every `Esd_4` tetrahedron must lie in exactly one `Esd_2`
   tetrahedron, with eight fine tetrahedra and the correct total volume in
   each coarse tetrahedron.
4. **Conformity.**  On two tetrahedra sharing a face, the subdivisions induced
   from the rank-ordered barycentric chambers must agree exactly on that face
   for `k=2,3`.
5. **Equivariance.**  For one regular tetrahedron, the union over all 24
   rank-ordered barycentric chambers followed by `Esd_k` must be invariant as
   an exact simplicial complex under all 24 vertex permutations for `k=2,3`.
   Since face dimension is preserved by every simplicial automorphism, this
   is also the local mechanism for 600-cell equivariance.
6. **Choice negative control.**  Applying ordered `Esd_2` directly to an
   unranked regular tetrahedron must produce three distinct subdivisions.  The
   rotational stabilizer `A4` must act without a fixed subdivision.  Thus the
   preliminary barycentric rank is essential rather than cosmetic.
7. **Symmetrization negative control.**  The fully tetrahedral-symmetric
   midpoint-plus-centroid 12-split must contain a repeatable central child
   whose exact singular-value ratio doubles on every repetition.  It is
   canonical but not uniformly shape regular.

## Acceptance and kill boundaries

**Acceptance for this numerical-geometry gate:** all seven gates pass.  This
establishes a canonical, conforming, nested, uniformly shape-regular tower for
the 600-cell carrier.  It selects the *refinement rule*, not an absolute length,
time step, stiffness coefficient, causal speed, or physical continuum limit.

**Kill for this construction:** any failure of equivariance, conformity,
nesting, or the general finite-shape argument rejects the rank-edgewise tower.

A failure of this construction is not a theorem against all possible
refinements.  Likewise, passing it is not evidence for a favorable spectrum;
the Whitney/stiffness calculation on the new carrier is a later, separately
preregistered experiment.

## Labels fixed in advance

- **DERIVED** if the exact gates establish the construction.
- **STRUCTURAL** for consequences inherited from the cited edgewise theorem.
- **PATTERN** only for optional finite-level numerical trends.
- **OPEN** for spectral convergence and all physical interpretation.

