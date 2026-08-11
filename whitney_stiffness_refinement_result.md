# Unweighted Whitney stiffness is not refinement-covariant

Date: 2026-08-11

Preregistration commit: `03e0abc`

Targeted verifier:
`reproducible/verify_whitney_stiffness_refinement.py`

Targeted result: **14/14 PASS**.  The verifier is registered.  The full suite
was not run by explicit user request.

## Headline

The complete 600-cell calculation reaches a sharp first-refinement negative:

> **DERIVED NUMERICAL NEGATIVE:** the canonical unweighted copy Laplacian
> (C^*C), with one common stiffness (kappa), is not invariant under the
> first barycentric refinement.  Preserving the same all-degree relative
> Schur guarantee requires increasing (kappa) by a factor
> (4.0816191349).

There is a stronger internal mismatch:

> **DERIVED NUMERICAL NEGATIVE:** the three constrained form degrees require
> incompatible first-step rescalings, approximately
> ((0.5403,3.8070,34.1729)).  No common rescaling keeps their separate
> stiffness balances fixed.

This does not kill every finite-stiffness route.  It kills the simplest
coefficient-free proposal: attach equal weight to every neighbour difference
and keep one algebraic (kappa) unchanged through subdivision.

## Complete calculation

Both carriers are complete, not sampled:

| level | f-vector | duplicated dimension | constraint quotient rank |
|---|---|---:|---:|
| base | `(120,720,1200,600)` | 9,000 | 6,360 |
| first barycentric | `(2640,17040,28800,14400)` | 216,000 | 153,120 |

The refined element metric is integrated exactly on every barycentric flag
orthoscheme inside the regular reference tetrahedron.  All 24 rank-ordered
child metrics agree exactly by tetrahedral symmetry.

The positive constraint spectra are computed without passing through the
large conforming zero eigenspace.  For every global simplex, an orthonormal
basis of the image of its occurrence-graph incidence matrix removes exactly
the row-cycle redundancy.  The remaining positive operator has dimensions
equal to the exact constraint quotient ranks above.

The method first reproduces the independently certified boundary-of-4-simplex
control:

\[
 g=7.5,
 \qquad
 q=45,
\]

with maximum Ritz residual below (1.5\times10^{-15}).

## Base and refined gaps

Let (a_r) be the mass-orthonormal local Kähler--Dirac norm and
(g_{r,p}) the smallest positive eigenvalue of the degree-(p) copy
penalty.

The complete results are:

| level | (a_r) | (g_{r,0}) | (g_{r,1}) | (g_{r,2}) | controlling degree |
|---|---:|---:|---:|---:|---:|
| base | 3.872983 | 1.433283 | 8.291796 | 12.000000 | 0 |
| refined | 11.936870 | 8.175537 | 6.712876 | 1.082291 | 2 |

Every quotient gap is strictly positive at both finite levels.  Therefore
the finite-complex stiff mechanism survives refinement algebraically.

But the bottleneck changes character completely: base conformity is limited
by vertex-copy modes, while refined conformity is limited by triangle-copy
modes.

The largest positive penalty eigenvalues are:

| degree | base | refined |
|---:|---:|---:|
| 0 | 39.270510 | 1080.000000 |
| 1 | 54.270510 | 266.102115 |
| 2 | 20.000000 | 14.248885 |

Thus refinement also broadens the stiff-sector scales very unevenly.

All complete extremal Ritz residuals are below
(3.1\times10^{-11}), versus the preregistered (10^{-8}) gate.  An
additional three-mode extremal audit, invoked with
`--multiplicity-audit`, reproduced the same edges and classifications.

## What the factor 4.0816 means

The target-free relative Schur bound requires

\[
 \kappa
 \geq
 \frac{a_r}{g_r}
 \left(2+\frac1\epsilon\right),
\]

where (g_r) is the smallest all-degree positive penalty gap.  The geometric
scale factors are

\[
 \frac{a_0}{g_0}=2.7021770885,
 \qquad
 \frac{a_1}{g_1}=11.0292577104.
\]

Their ratio is

\[
 R=4.0816191349.
\]

For the three frozen diagnostic accuracies, the sufficient stiffnesses are:

| relative bound (epsilon) | base (kappa) | refined (kappa) |
|---:|---:|---:|
| (10^{-1}) | 32.4261 | 132.3511 |
| (10^{-2}) | 275.6221 | 1124.9843 |
| (10^{-3}) | 2707.5814 | 11051.3162 |

These are sufficient theorem bounds, not fitted optimal values and not
physical accuracies.

## Why the form-degree mismatch matters

Applying the same comparison separately to the three constrained degrees
gives

\[
 R_0=0.5403316,
 \qquad
 R_1=3.8070171,
 \qquad
 R_2=34.1729051.
\]

A common worst-case rescaling can still control the full system, but it cannot
preserve the relative stiffness of all degrees.  It over-stiffens some
sectors and underlies a moving bottleneck in others.

This is not a small numerical discrepancy.  The largest and smallest ratios
differ by a factor greater than 63.

## Framing attack: this exposes a missing geometric weight

The result is conditional on the preregistered penalty

\[
 C^*C,
\]

which assigns unit energy to every neighbour difference.  That choice is
combinatorially canonical, but it is not the only defensible geometric
notion of stiffness.

Under refinement, faces shrink, Whitney (p)-form masses scale differently
with degree, and a continuum interface energy would normally carry
degree-dependent metric and face-size weights.  None of those weights is
present in the unweighted graph Laplacian.

Therefore the correct conclusion is not that all local stiffness is dead.
It is:

> **DERIVED DIAGNOSIS:** the unweighted occurrence-graph energy is not a
> refinement-covariant physical energy.  Any surviving route must derive its
> interface weights from the Whitney geometry before introducing an overall
> stiffness.

This also explains why the previous single-(kappa) model had unresolved
units.  The degree mismatch is quantitative evidence for that objection,
not merely dimensional rhetoric.

## What remains unselected

Even if a metric-weighted interface form supplies the relative degree and
refinement factors, it may still leave one overall positive dimensionless
constant.  Choosing that constant to reproduce a desired hierarchy would
remain fitting.

One refinement step also cannot establish an asymptotic power law.  No
scaling exponent is reported.

## Status ledger

- **DERIVED:** exact complete base/refined carriers and constraint ranks.
- **DERIVED:** exact positivity of every quotient penalty at both finite
  levels.
- **DERIVED NUMERICAL:** all six extremal positive spectra with Ritz residual
  below (3.1\times10^{-11}).
- **DERIVED NUMERICAL NEGATIVE:** constant unweighted (kappa) fails the
  same first-step relative guarantee by factor (4.0816191349).
- **DERIVED NUMERICAL NEGATIVE:** degreewise factors
  `(0.5403,3.8070,34.1729)` are incompatible.
- **DERIVED DIAGNOSIS:** unweighted copy energy lacks refinement covariance.
- **STRUCTURAL:** a common worst-case rescaling can preserve the sufficient
  all-degree bound.
- **OPEN:** metric/face-weighted canonical interface energy.
- **OPEN:** its absolute dimensionless coefficient.
- **OPEN:** repeated-refinement scaling and a uniformly bounded causal limit.
- **NOT CLAIMED:** physical time, mass, inertia, (c), (hbar), Newton's
  (G), or a Planck scale.

## Reproduction

Routine targeted verification:

```bash
/home/razvan/science/.venv/bin/python \
  reproducible/verify_whitney_stiffness_refinement.py
```

Expected result: `14/14`.

Optional slower three-mode edge audit:

```bash
/home/razvan/science/.venv/bin/python \
  reproducible/verify_whitney_stiffness_refinement.py \
  --multiplicity-audit
```
