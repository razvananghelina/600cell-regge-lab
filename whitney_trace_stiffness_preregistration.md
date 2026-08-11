# Preregistration: exact Whitney trace-jump stiffness

Date: 2026-08-11

## Motivation and framing attack

The complete first-refinement test rejected the unweighted occurrence-graph
energy (C^*C) as refinement-covariant.  Unit edge weights are canonical
combinatorially, but they ignore the metric size of a shared face and the
different scaling of Whitney form degrees.

This protocol tests one geometrically selected replacement, not a family of
weights:

> assign to every interelement jump its exact Whitney (L^2) trace norm on
> the shared triangle.

No diagonal mass lumping, inverse face size, tunable degree weight, or fitted
coefficient is admitted.  There is exactly one candidate bilinear form, up to
one overall scalar (kappa).

## Exact operator

Let (F) be a triangle shared by tetrahedra (T_-) and (T_+).  For form
degree (p=0,1,2), let (R_{F,p}) return the coefficient vector of the trace
jump

\[
 \operatorname{tr}_{F,T_-}u_-
 -\operatorname{tr}_{F,T_+}u_+,
\]

in the globally oriented Whitney (p)-simplex basis of (F).  The vector has
dimensions (3,3,1) for (p=0,1,2).

Let (M_{F,p}) be the exact Whitney (p)-form Gram matrix integrated with the
metric induced on (F).  Define

\[
 B_p=\sum_F R_{F,p}^*M_{F,p}R_{F,p},
 \qquad
 B=\operatorname{diag}(B_0,B_1,B_2,0).
\]

Every (M_{F,p}) is positive definite, so

\[
 \ker B=\ker R=\operatorname{im}J
\]

whenever the occurrence graphs are connected.  Thus the new form changes the
metric on constraint violations, not the exact physical kernel.

The finite-stiffness pencil remains

\[
 (W+\kappa B)v=zMv.
\]

It is Hermitian and element/face local.  Like the unweighted positive term,
it is degree-preserving and therefore breaks exact Kähler--Dirac oddness at
finite (kappa).

## Dimensional derivation

Under a uniform spatial dilation by (h) in dimension three, the element
Whitney mass of a (p)-form scales as

\[
 M_{T,p}\sim h^{3-2p}.
\]

On a two-dimensional face, its exact trace mass scales as

\[
 M_{F,p}\sim h^{2-2p}.
\]

Consequently the mass-orthonormal trace stiffness scales as

\[
 M_{T,p}^{-1/2}B_pM_{T,p}^{-1/2}\sim h^{-1},
\]

independently of (p).  The mass-orthonormal Kähler--Dirac operator also
scales as (h^{-1}).  Therefore the trace mass, with no extra power of (h),
is the unique direct (L^2)-trace candidate having the correct first-order
dimension for every form degree.

Adding the familiar (h^{-1}) interior-penalty factor would instead give a
second-order (h^{-2}) scale and is excluded from this first-order test.

## Frozen independent control

Use the boundary of a 4-simplex made from five congruent regular
tetrahedra, followed by its complete first barycentric subdivision.

At the base level:

- duplicated dimension: 75;
- assembled dimension: 30;
- tetrahedra: 5;
- triangles: 10.

At the refined level:

- duplicated dimension: 1,800;
- assembled f-vector: `(30,150,240,120)`;
- assembled dimension: 540;
- tetrahedra: 120;
- triangles: 240.

Every refined vertex is the barycentre of a base simplex.  Its coordinates
inside each parent tetrahedron are constructed exactly.  Face metrics seen
from the two parents must agree exactly; disagreement is a kill condition.

## Frozen outputs

For each level and (p=0,1,2), compute the complete generalized spectrum

\[
 B_pv=\lambda M_pv
\]

and record:

1. its exact nullity and expected conforming dimension;
2. smallest positive eigenvalue (g_{r,p});
3. largest eigenvalue (q_{r,p});
4. local Dirac norm (a_r);
5. degreewise scale (s_{r,p}=a_r/g_{r,p});
6. first-step ratios (R_p=s_{1,p}/s_{0,p});
7. ratio spread (max R_p/\min R_p).

Dense symmetric generalized diagonalization is used on the complete control,
not an iterative selected-mode search.  Eigen residuals must be below
(10^{-9}), and eigenvalues with magnitude below (10^{-10}) are classified
as zero only if the numerical nullity equals the exact assembly dimension.

## Decision protocol

- If every (R_p=1) within (10^{-8}), label exact trace stiffness
  **DERIVED NUMERICAL: first-step refinement-compatible on the control**.
- Otherwise label exact equality **DERIVED NUMERICAL NEGATIVE**, and print the
  full ratios without fitting a common exponent.
- The previously committed unweighted ratio spread is approximately 63.24.
  A smaller spread is reported only as **PATTERN: improved degree balance**,
  because this comparison is not blind and one control is insufficient.
- A larger or equal spread is a clean negative for the proposed geometric
  repair.
- No tolerance other than exact-equality and eigensolver certification gates
  determines success.

## Scope and kill boundary

This control can falsify the exact trace proposal cheaply.  It cannot prove a
complete 600-cell or continuum scaling law.

The route is killed if:

1. face metrics disagree across an identified shared face;
2. the kernel differs from conformity;
3. the form is nonlocal or non-Hermitian;
4. degree balance is no better than the unweighted candidate.

If the degree spread improves, the next required step is the complete
600-cell/base-refinement calculation.  It must be separately preregistered.

## Status before execution

- **DERIVED:** exact dimensional covariance of the face trace mass under
  uniform dilation.
- **STRUCTURAL:** using a trace-jump Hamiltonian as microscopic dynamics.
- **OPEN:** complete base/refined spectra on the independent control.
- **OPEN:** comparison with the unweighted degree mismatch.
- **OPEN:** one overall dimensionless stiffness.
- **OPEN:** complete 600-cell and repeated-refinement behavior.
- **NOT CLAIMED:** time, mass, inertia, (c), (hbar), Newton's (G), or a
  Planck scale.
