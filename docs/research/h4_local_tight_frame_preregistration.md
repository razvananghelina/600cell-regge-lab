# Preregistration: can the local H4 metric select a scalar tight frame?

Date: 2026-08-11

## Question frozen before computation

The connected three-bond walk of commits `9bf2fba`--`4ebec8c` is exactly
unitary and globally connected, but its one-period probability distribution
is strongly anisotropic.  Before introducing a new coin, test the smaller
necessary metric question:

> Do the four geometrically derived colour-neighbour steps of one barycentric
> H4 chamber admit **positive scalar weights** with zero drift and isotropic
> second moment?

No speed, mass, Standard-Model number or desired continuum exponent enters
this test.

## Frozen geometry

Reconstruct the 600-cell and its 14,400 complete flags exactly as in
`reproducible/verify_h4_three_bond_local_isotropy.py`.  Use chamber 0 only
after verifying that the four-step Gram matrix is chamber-independent under
the H4 action.

Let `x` be the normalized spherical centre of the reference flag chamber,
and let `s_i x` be its neighbour across rank colour `i`, for `i=0,1,2,3`.
Define the literal geodesic tangent steps

\[
d_i=\log_x(s_i x)\in T_xS^3,
\qquad \ell_i=\|d_i\|,
\]

and their unit directions

\[
u_i=d_i/\ell_i.
\]

The two variants are fixed in advance and both must be reported:

1. **literal-step variant:** `v_i=d_i`, retaining the unequal derived edge
   lengths;
2. **direction-only variant:** `v_i=u_i`, corresponding to an additional
   colour-dependent reparametrization to equal physical step length.

The direction-only variant is not allowed to rescue the literal variant
silently.  It represents a different physical hypothesis.

Thus the complete look-elsewhere count is `N=2`; no subset of colours and no
post-result rescaling is admitted.

## Frozen linear feasibility problem

For each of the two variants seek four weights `p_i` and a scalar `c` such
that

\[
p_i>0,\qquad \sum_i p_i=1,
\]

\[
\sum_i p_i v_i=0,
\]

and

\[
\sum_i p_i v_i v_i^T=cP_x,
\qquad P_x=I_4-xx^T.
\]

The last equality is coordinate-free isotropy on the three-dimensional
tangent space.  These equations are linear in `(p_0,p_1,p_2,p_3,c)`.
Compute and report:

- the numerical rank of the zero-drift system and whether its normalized
  solution is unique;
- all four drift-cancelling weights and their minimum;
- the three tangent eigenvalues, eigenvalue ratio and normalized traceless
  residual at those weights;
- the rank and relative least-squares residual of the complete simultaneous
  system;
- whether an unconstrained exact solution exists, before applying positivity;
- whether a positive simultaneous solution exists.

An inconsistent unconstrained system already proves that no positive
solution exists.  A positive result requires both the residual test and
strict positivity; merely minimizing anisotropy does not count.

## Numerical standard and calibration

Use two computations:

1. double precision from the repository geometry;
2. at least 80-decimal arithmetic after reconstructing each 600-cell
   coordinate from
   `{0, +/-1, +/-1/2, +/-phi/2, +/-1/(2phi)}`.

The conclusion must be stable between them.  Report condition numbers.  A
linear equality is accepted numerically only below relative residual
`1e-30` in the high-precision computation; a residual above `1e-12` is a
robust negative.  Any residual in between is **OPEN NUMERICAL**, not rounded
to either answer.

Before applying the estimator, calibrate it on the four regular-tetrahedron
directions

\[
(1,1,1),\ (1,-1,-1),\ (-1,1,-1),\ (-1,-1,1)
\]

divided by `sqrt(3)`.  It must recover unique weights `p_i=1/4`, zero drift,
and covariance `I_3/3`.

## Cheap extensions excluded in advance

Adding a stay-put probability cannot repair the tensor shape: it multiplies
both the drift and nonzero-step second moment by the same factor.  Verify this
algebraically, but do not count it as a third attempt.

Independent fitted step lengths, fitted colour multiplicities, signed
weights, subsets of colours and a posteriori choices among chambers are not
admissible.  A matrix-valued tetrad/coin is a genuinely larger construction
and remains outside this scalar test.

## Decision boundaries

- **DERIVED SCALAR METRIC TIGHT FRAME:** the literal geodesic steps have a
  unique positive simultaneous solution.
- **STRUCTURAL DIRECTION-ONLY TIGHT FRAME:** only the unit directions pass.
  This requires an extra colour-dependent clock/length prescription and is
  not yet physical dynamics.
- **DERIVED NUMERICAL SCALAR NO-GO AT FIRST SCALE:** neither complete system
  is compatible, with high-precision relative residual above `1e-12`.
- **OPEN NUMERICAL:** the residual lies in the frozen ambiguity interval.

## Hostile scope boundary

This is a necessary scalar diffusion/tight-frame gate, not a sufficient
quantum-walk continuum theorem.  In particular, a unitary 4-by-4 coin acting
on the maximally mixed internal state leaves that state maximally mixed, so
nonuniform `p_i` cannot automatically be interpreted as its four outgoing
probabilities.  Even a positive tight frame would still require a derived
unitary realization and a Dirac-dispersion calculation.

Conversely, a negative result closes only scalar reweighting of these four
fixed local directions.  It does not close position-dependent matrix coins,
larger ancillas or a genuine discrete tetrad construction.
