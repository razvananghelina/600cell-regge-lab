# Preregistration: one-period local isotropy of the three-bond H4 walk

Date: 2026-08-11

## Starting point

Commit `02e020c` proved exact strong support connectivity for all six orders
of the three Coxeter-bond translations.  Connectivity is necessary but says
nothing about whether the walk treats the three physical tangent directions
isotropically.

The single-bond orbit lengths (6,6,10) already disclose microscopic
combinatorial anisotropy.  No displacement moment has yet been computed.

## Geometry frozen before execution

1. Rebuild the regular 600-cell vertices on the unit three-sphere.
2. For every cell, use its normalized Euclidean vertex barycentre, exactly as
   in the repository's established first barycentric refinement.
3. For every complete flag (v\subset e\subset f\subset t), define its
   chamber representative as the normalized Euclidean mean of its four
   barycentric vertices.
4. Use the intrinsic rank-colour involutions and the three translations
   (T_{01},T_{12},T_{23}) frozen in commit `02e020c`.

No coordinate axes are fitted to an outcome.

## Coin and initial state

Use the paper's directional coin, with the standard spin-rotation convention

\[
C=e^{i\pi/3}R_{\sigma_z}(\pi/2)R_{\sigma_x}(\pi/2),
\qquad
R_{\sigma}(\theta)=e^{-i\theta\sigma/2},
\]

and

\[
\widehat C=I_2\otimes C.
\]

The global phase has no effect on the test but will be retained.  No coin
angle is varied.

Start at chamber index zero with the maximally mixed density on its four
active components:

\[
\rho_0=|k_0\rangle\langle k_0|\otimes I_4/4.
\]

This removes a fitted spin-polarization direction.  Check normalization by
evolving all four orthogonal input basis states and averaging probabilities.

## Complete schedule set

Retain the same preregistered (N=6) temporal orders.  The designated order
is

\[
(01)\to(12)\to(23).
\]

Apply one complete period, with the fixed coin after each translation.

## Observable fixed in advance

For every final chamber centre (y), compute the geodesic logarithm at the
initial centre (xin S^3):

\[
v_y=\log_x(y)
=\frac{\arccos(x\cdot y)}{\sqrt{1-(x\cdot y)^2}}
\left[y-(x\cdot y)x\right].
\]

Use the evolved chamber probability to compute:

\[
\mu=\mathbb E[v],
\qquad
\Sigma=\mathbb E[(v-\mu)(v-\mu)^T].
\]

Restrict (Sigma) to the three-dimensional tangent space (x^\perp) and
record all three eigenvalues, their ratio

\[
R=\lambda_{\max}/\lambda_{\min},
\]

and the normalized traceless residual

\[
A=\frac{\|\Sigma-(\operatorname{tr}\Sigma/3)I_T\|_F}
        {\operatorname{tr}\Sigma}.
\]

No peak, plateau, exponent or target value is fitted.

## Calibration

Before using the estimator, apply it to the four unit directions from the
centre of a regular tetrahedron with equal probabilities.  The control must
return zero mean, (R=1) and (A=0) to numerical precision.

## Decision and interpretation boundaries

- **ONE-PERIOD ISOTROPY:** (|\mu\|<10^{-12}), (|R-1|<10^{-10}) and
  (A<10^{-10}).
- **DERIVED NUMERICAL LOCAL ANISOTROPY:** any of those gates fails by more
  than the frozen tolerances.

Report the hit fraction among all six schedules.  A rare isotropic order
would be labelled PATTERN until independently selected; six hits out of six
would show that order is irrelevant.

A failure is not a continuum no-go.  It says only that the one-period walk
has drift or unequal tangent variance at the first barycentric scale.  A
multistep/refinement isotropization mechanism would remain OPEN.

This test does not compare with the Dirac dispersion, the Whitney operator,
particle masses, (c), or Planck units.
