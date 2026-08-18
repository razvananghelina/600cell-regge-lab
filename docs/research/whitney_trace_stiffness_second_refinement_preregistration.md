# Preregistration: second refinement of exact Whitney trace stiffness

Date: 2026-08-11

## Question

The complete 600-cell first step produced exact-trace stiffness ratios in the
narrow band `(3.25,4.38)`, but one refinement cannot distinguish a stable
flow from an accident.  This protocol performs a second barycentric
refinement on the independently calibrated closed boundary-of-4-simplex
control.

No factor near four is assumed or used as a target.

## Complete three-level control

Start with the boundary of a 4-simplex and apply barycentric subdivision
twice.  Expected exact carriers are:

| level | f-vector | tetrahedra | duplicated dimension |
|---:|---|---:|---:|
| 0 | `(5,10,10,5)` | 5 | 75 |
| 1 | `(30,150,240,120)` | 120 | 1,800 |
| 2 | `(540,3420,5760,2880)` | 2,880 | 43,200 |

Every new vertex is the barycentre of a simplex at the preceding level.
Coordinates are computed exactly in each parent tetrahedron chart.  All
shared-face metrics must agree exactly from both parents.

Unlike the first refinement of a regular tetrahedron, the second refinement
may contain several ordered element metric types.  Every type is enumerated
from its exact Gram matrix; replacing them by one regular or averaged element
is forbidden.

## Fixed operator

Use only the already derived exact trace-jump form

\[
 B_{r,p}=\sum_F R_{F,p}^*M_{F,p}R_{F,p},
 \qquad p=0,1,2,
\]

with the exact level-(r) duplicated Whitney element mass (M_{r,p}).

The candidate count remains one.  No degree weight, face exponent, lumping,
or overall stiffness is selected.

## Quotient spectra and frozen solvers

Remove the exact conforming kernel by the connected occurrence-graph
row-image construction already calibrated in
`verify_whitney_trace_stiffness_full.py`.

For degrees zero and one use generalized block LOBPCG with:

- block size 5;
- tolerance (10^{-9});
- maximum 2,000 iterations;
- seed `60020260811 + 100*level + degree`.

For degree two use the corrected full-row symmetric operator

\[
 H_2^{1/2}R_2M_2^{-1}R_2^*H_2^{1/2}
\]

and symmetric Lanczos with five eigenpairs at each edge, tolerance
(10^{-11}), and maximum 20,000 iterations.

Every reported spectral edge must have directly recomputed relative Ritz
residual below (10^{-7}).

## Mandatory calibration

Before level two, the implementation must reproduce all dense exact-trace
gaps and maxima at levels zero and one from
`whitney_trace_stiffness.json`, with:

- relative spectral error below (5\times10^{-7});
- relative Ritz residual below (10^{-7}).

Calibration failure invalidates the implementation and stops the run before
interpretation.

## Frozen outputs

For every level (r=0,1,2) and degree (p=0,1,2), record:

1. exact f-vector, row count, rank, redundancy, and kernel dimension;
2. number of exact element and face metric types;
3. five smallest and five largest positive quotient eigenvalues;
4. maximum Ritz residual;
5. local block-Dirac norm
   (a_r=\max_T\lVert A_{r,T}\rVert_2);
6. scale (s_{r,p}=a_r/g_{r,p}).

Then record the two step ratios

\[
 R^{01}_p=\frac{s_{1,p}}{s_{0,p}},
 \qquad
 R^{12}_p=\frac{s_{2,p}}{s_{1,p}},
\]

their degree spreads, and the componentwise flow ratios

\[
 Q_p=\frac{R^{12}_p}{R^{01}_p}.
\]

No average factor or exponent is fitted.

## Decision protocol

- If every (Q_p=1) within (5\times10^{-7}), label the step factor
  **DERIVED NUMERICAL: repeated once on the control**.  Do not call it an
  asymptotic law.
- Otherwise label exact step-factor repetition **DERIVED NUMERICAL NEGATIVE**.
- If the degree spread at step `1->2` is smaller than at step `0->1`, label
  **PATTERN: flow toward a common degree scaling**.
- If it is larger, label **PATTERN NEGATIVE: degree balance worsens**.
- If any level-two quotient gap vanishes, the finite-stiffness route reaches
  its kill boundary on this control.
- In every outcome, two steps remain insufficient for a continuum exponent
  or fixed-point theorem.

## Framing attacks fixed in advance

The local norm uses the worst exact element type.  This gives a sufficient
all-element Schur scale, not an average physical propagation speed.

The control topology is not the complete second-refined 600-cell.  A stable
factor here would justify, but not replace, that much larger calculation.

The overall (kappa), Kähler--Dirac oddness, the singular exact limit, and
Lorentzian time remain separate unsolved gates.

## Status before execution

- **DERIVED INPUT:** exact trace stiffness and uniform (h^{-1}) dimension.
- **DERIVED INPUT:** level-zero and level-one dense calibration spectra.
- **OPEN:** exact level-two metric-type census.
- **OPEN:** level-two quotient spectra and step ratios.
- **OPEN:** repeated-refinement law, absolute stiffness, chirality, and
  causal dynamics.
- **NOT CLAIMED:** mass, inertia, (c), (hbar), Newton's (G), or a Planck
  scale.
