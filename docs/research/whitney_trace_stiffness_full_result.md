# Complete exact trace stiffness nearly aligns the three refinement scales

Date: 2026-08-11

Original preregistration: `7ed6d49`  
Recorded preregistered solver failure: `9c04a94`  
Degree-two solver correction frozen before rerun: `62ca607`

Targeted verifier:
`reproducible/verify_whitney_trace_stiffness_full.py`

Corrected targeted result: **12/12 PASS**.  The verifier is registered.  The
full suite was not run by explicit user request.

## Headline

The exact Whitney trace-jump metric does not give an exact first-refinement
fixed point, but on the complete 600-cell it almost aligns the three form
degrees:

\[
 (R_0,R_1,R_2)
 =
 (4.37932,,3.24849,,3.80924).
\]

Their spread is

\[
 \frac{\max R_p}{\min R_p}=1.34811.
\]

The paired unweighted complete penalty had ratios

\[
 (0.540332,,3.80702,,34.1729)
\]

and spread

\[
 63.2443.
\]

Therefore:

> **PATTERN:** exact geometric trace weighting reduces the complete
> degree-refinement mismatch by a factor of about 46.9 in spread.

But:

> **DERIVED NUMERICAL NEGATIVE:** none of the three ratios is one and they are
> not equal to each other.  Exact first-step covariance remains false.

This is the strongest constructive evidence so far that the missing local
conformity energy should be metric and face based.  It is not yet a physical
renormalization law.

## Complete carriers and exact structure

The calculation uses every cell at both levels:

| level | f-vector | duplicated dimension | positive quotient rank |
|---|---|---:|---:|
| base | `(120,720,1200,600)` | 9,000 | 6,360 |
| first barycentric | `(2640,17040,28800,14400)` | 216,000 | 153,120 |

For each shared triangle (F), the penalty is the exact trace Gram form

\[
 B_p=\sum_F R_{F,p}^*M_{F,p}R_{F,p}.
\]

All occurrence graphs are connected, all shared-face metrics agree exactly
from both parents, and every row-image basis is orthonormal to within
(1.11\times10^{-15}).

The candidate count is one.  No degree coefficient, face-size exponent, or
spectral target was searched.

## Calibrated quotient solver

The large conforming nullspaces are removed exactly by passing to the image
of each occurrence-graph incidence matrix.  The positive trace spectrum is
then obtained from the reduced generalized pencil

\[
 V^*HRM^{-1}R^*HV,x
 =
 \lambda V^*HV,x.
\]

Before the complete calculation, the solver reproduced all 12 dense edges
from the base/refined boundary-of-4-simplex control:

- maximum relative eigenvalue error:
  (7.11\times10^{-15});
- maximum recomputed Ritz residual:
  (7.52\times10^{-9}).

The complete maximum recomputed residual is
(9.70\times10^{-9}), below the preregistered (10^{-7}) gate.

## Recorded failure and correction

The first preregistered run returned `10/11`.  Its refined degree-two largest
LOBPCG vector had residual (6.56\times10^{-6}), so the provisional ratios
were explicitly rejected and committed in `9c04a94`.

Degree two has no row-cycle redundancy: every triangle supplies one
independent jump.  Commit `62ca607` therefore froze a symmetric Lanczos solve
of

\[
 H_2^{1/2}R_2M_2^{-1}R_2^*H_2^{1/2}
\]

before recomputation, without changing the spectral definition or the
residual gate.

The corrected degree-two residual is (9.95\times10^{-12}), and its complete
base edge agrees with the already accepted calculation to
(1.22\times10^{-15}) relatively.

## Complete positive gaps

| degree | base gap | refined gap | base maximum | refined maximum |
|---:|---:|---:|---:|---:|
| 0 | 0.821319 | 0.578029 | 15.666547 | 62.925287 |
| 1 | 2.149657 | 2.039542 | 10.720025 | 47.837411 |
| 2 | 3.464102 | 2.802834 | 5.773503 | 22.312477 |

Every finite quotient remains positively gapped.  The lowest gap is controlled
by 0-form copy modes at both levels.

The local mass-orthonormal Dirac norms are

\[
 a_0=3.8729833462,
 \qquad
 a_1=11.9368695071.
\]

This gives degreewise scales (a/g):

| degree | base (a_0/g_{0,p}) | refined (a_1/g_{1,p}) | ratio |
|---:|---:|---:|---:|
| 0 | 4.715568 | 20.651001 | 4.379324 |
| 1 | 1.801675 | 5.852720 | 3.248489 |
| 2 | 1.118034 | 4.258856 | 3.809237 |

## What is genuinely new

The unweighted penalty changed its controlling degree and demanded wildly
different scalings.  Exact trace weighting keeps the same controlling degree
and compresses all three first-step ratios into a narrow order-four band.

A single worst-case factor

\[
 4.379324
\]

is therefore sufficient to preserve the same all-degree Schur guarantee at
the first refinement.  This factor is derived from the worst quotient gap,
not fitted to an average.

However, it is only a first-step relative factor.  It does not determine the
base value of (kappa), and two levels cannot show whether the factor
stabilizes, flows, or diverges.

## Framing attack

Calling the clustered ratios a renormalization-group law would be premature:

1. there is only one refinement step;
2. barycentric subdivision changes element shape as well as scale;
3. the three ratios still differ by about 35%;
4. no physical observable fixes the overall stiffness;
5. the positive term still breaks exact Kähler--Dirac oddness;
6. exact Whitney recovery still needs
   (kappa\rightarrow\infty), whose microscopic norm is unbounded.

Thus the result advances the mathematical architecture but not yet the
physical tick.

## Status ledger

- **DERIVED:** unique exact Whitney trace-jump energy.
- **DERIVED:** exact complete kernels, ranks, locality, and face agreement.
- **DERIVED NUMERICAL:** all complete base/refined positive spectral edges.
- **DERIVED NUMERICAL NEGATIVE:** exact first-step covariance fails.
- **PATTERN:** degreewise ratios cluster in `(3.25,4.38)`.
- **PATTERN:** spread improves from `63.2443` to `1.34811`.
- **DERIVED STRUCTURAL:** worst-case first-step factor `4.379324` preserves
  the same sufficient all-degree bound.
- **OPEN:** second and later refinements.
- **OPEN:** an asymptotic or fixed-point law.
- **OPEN:** absolute dimensionless stiffness.
- **OPEN:** chiral finite-stiffness dilation.
- **OPEN:** uniformly bounded causal dynamics.
- **NOT CLAIMED:** physical time, mass, inertia, (c), (hbar), Newton's
  (G), or a Planck scale.

## Reproduction

```bash
/home/razvan/science/.venv/bin/python -u \
  reproducible/verify_whitney_trace_stiffness_full.py
```

Expected corrected result: `12/12`.
