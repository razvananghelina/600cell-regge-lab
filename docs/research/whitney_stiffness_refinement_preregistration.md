# Preregistration: refinement scaling of canonical Whitney stiffness

Date: 2026-08-11

## Question

The finite-stiffness result proves a local low-energy Whitney limit on a
small control, but leaves the coefficient (kappa) unselected.  This test
asks the next necessary question:

> Does the complete 600-cell geometry determine a stable relative scaling of
> the canonical copy stiffness under its first barycentric refinement?

No measured mass, speed, time, cutoff, or phenomenological hierarchy is used.
The test compares only the base geometry with its derived first subdivision.

## Frozen operators and complete hypotheses

At each level (r=0,1), use:

- the complete closed 600-cell triangulation and its full first barycentric
  subdivision;
- the duplicated element carrier, ordered by tetrahedron and local simplex;
- the exact element Whitney metric (M_r);
- the exact element weak Kähler--Dirac matrix (W_r);
- every canonical face-neighbour copy-difference row (C_r), with
  coefficients (+1,-1).

The refined local metric is integrated on the flag orthoscheme

\[
 v\subset e\subset f\subset t
\]

inside the same regular reference tetrahedron.  Child tetrahedra are not
replaced by rescaled regular tetrahedra.

The tested mass-orthonormal matrices are

\[
 A_r=M_r^{-1/2}W_rM_r^{-1/2},
 \qquad
 L_r=M_r^{-1/2}C_r^*C_rM_r^{-1/2}.
\]

Because (M_r,W_r) are element block diagonal, the Dirac norm

\[
 a_r=\lVert A_r\rVert_2
\]

is computed exactly from one 15-dimensional element type at each level.

For form degree (p=0,1,2), define

\[
 g_{r,p}=\lambda_{min}^{+}(L_{r,p}),
 \qquad
 q_{r,p}=\lambda_{max}(L_{r,p}).
\]

The top degree has no copy constraints and is excluded from positive-gap
ratios.  The all-degree gap is

\[
 g_r=\min_{p=0,1,2}g_{r,p}.
\]

## Avoiding the large zero eigenspace

Directly asking an eigensolver for the smallest eigenvalues of (L_r) would
return tens of thousands of exact conforming zero modes.  That is not an
admissible gap calculation.

Instead use the equality of positive spectra

\[
 \sigma_+(L_{r,p})
 =
 \sigma_+\bigl(C_{r,p}M_{r,p}^{-1}C_{r,p}^*\bigr).
\]

For every global simplex, construct an orthonormal basis of the column space
of its connected occurrence-graph incidence matrix.  Their direct sum is an
orthonormal matrix (V_{r,p}) spanning the complete row image of
(C_{r,p}).  Compute the extremal eigenvalues of the strictly positive
operator

\[
 V_{r,p}^*C_{r,p}M_{r,p}^{-1}C_{r,p}^*V_{r,p}.
\]

This removes only the exactly known row-cycle nullspace.  It does not choose
an independent physical penalty or replace the canonical all-row matrix.

The method is calibrated before use on the already certified boundary of a
4-simplex.  It must reproduce the previously obtained all-degree values

\[
 g=7.5,
 \qquad q=45
\]

within (10^{-9}), with relative Ritz residual below (10^{-8}).

## Frozen refinement diagnostics

The target-free Schur estimate for relative spectral error
(epsilon>0) requires

\[
 \kappa
 \geq
 \frac{a_r}{g_r}\left(2+\frac1\epsilon\right).
\]

Therefore the geometry-dependent quantity relevant to refinement is

\[
 s_r=\frac{a_r}{g_r}.
\]

Record before interpretation:

1. (a_r);
2. every (g_{r,p}) and (q_{r,p});
3. the degree attaining (g_r);
4. (s_r=a_r/g_r);
5. the first-step ratio (R=s_1/s_0);
6. the corresponding sufficient (kappa) values for the frozen diagnostic
   tolerances (epsilon=10^{-1},10^{-2},10^{-3});
7. the degreewise ratios
   ((a_1/g_{1,p})/(a_0/g_{0,p})).

No exponent is fitted from two levels.  The three epsilon values are not
physical accuracy targets; they only display the exact same theorem at three
fixed resolutions.

## Decision labels frozen before computation

- If (R=1) within the numerical certificate tolerance, label a constant
  algebraic stiffness **PATTERN: first-step compatible**.  Do not claim a
  continuum theorem.
- If (R\neq1), label constant (kappa) **DERIVED NEGATIVE at the first
  refinement** for preserving the same relative Schur guarantee.
- If all three degreewise ratios agree within (10^{-8}), label their common
  rescaling **PATTERN**.  With only two levels it is not a law.
- If they disagree, label a single degree-balanced rescaling **DERIVED
  NEGATIVE at level one**; a worst-case common (kappa) may still control all
  degrees but over-stiffens the others.
- In every outcome, the absolute normalization of (kappa) remains **OPEN**.
  A first-step ratio cannot create a physical unit.

## Acceptance and kill boundaries

This route advances mathematically if:

1. the row-image method passes the independent control;
2. all complete base and refined extremal Ritz pairs meet the residual gate;
3. the positive gaps remain nonzero at both finite levels;
4. the first-step relative scaling is thereby determined without a fit.

It advances physically only if later work derives an absolute stiffness and
a repeated-refinement law compatible with a uniformly bounded causal
generator.  Neither is claimed here.

The route is killed at this stage if the canonical refined penalty loses its
positive quotient gap, or if the row-image calculation cannot distinguish
that gap reproducibly from the exact nullspace.

## Status before execution

- **DERIVED INPUT:** complete base and first-refined carriers.
- **DERIVED INPUT:** exact local child Whitney metrics.
- **DERIVED INPUT:** exact kernels of the canonical neighbour constraints.
- **DERIVED INPUT:** finite-complex Schur convergence theorem.
- **OPEN:** complete 600-cell stiffness gaps at both levels.
- **OPEN:** first-step ratio (R).
- **OPEN:** absolute stiffness, repeated-refinement law, time, and causal
  continuum dynamics.
- **NOT CLAIMED:** mass, inertia, (c), (hbar), Newton's (G), or a Planck
  scale.
