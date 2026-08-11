# Preregistration: finite-stiffness Whitney conformity

Date: 2026-08-11

## Question and complete hypotheses

This test asks a narrower question than exact assembly:

> Can the exact Whitney Kähler--Dirac spectrum occur as the bounded
> low-energy limit of a finite, local, positive stiffness assigned to copy
> differences?

Fix the duplicated oriented element carrier, its positive block-local Whitney
mass matrix (M), its block-local weak Kähler--Dirac matrix (W=W^*), and the
assembly injection (J).  Retain **all** canonical neighbour rows (C): two
copies of the same simplex are compared exactly when their tetrahedra share a
triangle containing that simplex.  Every row has coefficients (+1,-1).
No spanning tree or independent-row basis is chosen.

The finite-stiffness generalized Hermitian pencil is

\[
 (W+\kappa C^*C)v=zMv,
 \qquad \kappa>0.
\]

Equivalently, in mass-orthonormal coordinates,

\[
 H_\kappa=A+\kappa L,
 \qquad
 A=M^{-1/2}WM^{-1/2},
 \qquad
 L=M^{-1/2}C^*CM^{-1/2}\geq0.
\]

The coefficient of every graph edge is one.  This makes (C^*C) invariant
under row ordering and row-sign conventions.  It does **not** select the
relative physical scale (kappa); that is one of the gates under test.

## Framing attack fixed before computation

This pencil is not an exact replacement for the Whitney operator at finite
(kappa): the existing universal pure-penalty theorem already proves that
the conforming subspace is not invariant.  The only admissible claim is an
effective low-energy limit.

Moreover, adding the positive degree-preserving block (kappa C^*C) breaks
the exact odd/chiral form of the Kähler--Dirac generator.  The construction is
a local Hermitian Hamiltonian with a stiffness sector, not itself a new
Kähler--Dirac operator.  Any physical interpretation must carry this caveat.

Finally, the algebraic normalization (+1,-1) does not supply physical units
for (kappa), nor a common refinement scaling across form degrees.  A good
finite-complex approximation alone will therefore be labelled STRUCTURAL,
not a derivation of time, mass, (c), or a Planck scale.

## Exact limit theorem to test

The neighbour certificate gives

\[
 \ker C=\operatorname{im}J.
\]

Let (P) be the orthogonal projector onto
(M^{1/2}\operatorname{im}J=\ker L), let (m=\dim\ker L), and define

\[
 a=\lVert A\rVert_2,
 \qquad
 g=\lambda_{m+1}(L)>0.
\]

The compression (PAP), restricted to (ker L), must have exactly the
assembled Whitney generalized spectrum

\[
 J^*WJ\,x=z\,J^*MJ\,x.
\]

For (kappa g>2a), Weyl separation gives exactly (m) eigenvalues below the
stiff sector.  Eliminating the orthogonal component gives the target-free
Schur bound

\[
 \operatorname{dist}
 \bigl(\sigma_{\rm bounded}(H_\kappa),\sigma(PAP)\bigr)
 \leq
 \frac{a^2}{\kappa g-2a}.
\]

The verifier will check the hypotheses and the corresponding numerical bound
on a completely independent small control.  This is a theorem-driven check;
no exponent will be fitted.

## Frozen small control

Use the boundary of a 4-simplex, with five tetrahedra and its complete
oriented cochain carrier:

- duplicated dimension: (75);
- assembled dimension: (30);
- all face-neighbour constraint rows: one for every nonempty subface of each
  shared triangle;
- expected row count: (70);
- expected rank: (45).

Local masses are recomputed by exact Whitney-form integration on the regular
reference tetrahedron.  The exact local weak operator is constructed from
the coboundary and those masses.  The canonical neighbour matrix is built
combinatorially, not imported from a prior numerical result.

The following dimensionless dyadic values are frozen before spectral
evaluation:

\[
 \kappa\in\{1,2,4,\ldots,2^{14}\}.
\]

For every value satisfying the independently computed separation condition,
record:

1. the 30 bounded eigenvalues selected by the separation ordering;
2. their maximum ordered error from the assembled spectrum;
3. the analytic Schur bound above;
4. the largest principal-angle sine from the conforming subspace;
5. the maximum constraint residual of the bounded eigenvectors;
6. the operator norm of (H_\kappa).

The full multiset is written to JSON.  Peaks, best points, and fitted slopes
will not be used.

## Locality and scale gates

The verifier must also record the exact sparsity support of (W+C^*C) on the
control.  Finite (kappa) passes the microscopic locality gate if every
nonzero couples either two degrees in the same tetrahedron or two copies
joined by a declared neighbour constraint.

At the same time, record the lower bound

\[
 \lVert H_\kappa\rVert_2
 \geq \kappa\lambda_{\max}(L)-a.
\]

If the bounded sector converges while this norm diverges, the correct verdict
is:

- **DERIVED:** a local finite-complex stiff realization of the Whitney
  low-energy limit exists;
- **DERIVED:** exact recovery requires the singular limit
  (kappa\rightarrow\infty);
- **STRUCTURAL NEGATIVE:** that limit has no uniformly bounded microscopic
  generator in the frozen normalization;
- **OPEN:** geometry selection of a finite stiffness and its refinement law.

## Acceptance and kill boundaries

The mathematical route advances only if all of the following hold:

1. (C^*C) is canonical and local under the stated neighbour rule;
2. its kernel is exactly the assembled subspace;
3. the compressed operator is exactly Whitney;
4. the separated finite-stiffness eigenvalues satisfy the preregistered
   target-free bound and approach the compressed spectrum.

The route is killed if any of these fail.  Even if all pass, the physical
tick does **not** advance unless a finite (kappa) and a refinement scaling
are selected independently by the geometry.  Choosing (kappa) to obtain a
desired hierarchy is forbidden fitting.

## Status before execution

- **DERIVED INPUT:** local Whitney masses and weak operator.
- **DERIVED INPUT:** canonical bounded-degree neighbour constraints.
- **DERIVED INPUT:** their kernel is exact conformity at the base and first
  barycentric levels.
- **OPEN:** finite-stiffness spectral convergence on the independent control.
- **OPEN:** selection and physical dimensions of (kappa).
- **OPEN:** uniform refinement scaling and causal continuum interpretation.
- **NOT CLAIMED:** Lorentzian time, inertia, mass, (c), (hbar), Newton's
  (G), or a Planck scale.
