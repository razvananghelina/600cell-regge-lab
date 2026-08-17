# Preregistration: nonlinear boundary-coboundary test of the two schedules

Date: 2026-08-17

Prior-art commit: `5ccb29b`.

Status: frozen before evaluating any new perturbed Dirichlet solution,
on-shell action difference or mixed rectangle.

## 1. Frozen inputs and exclusions

Require the exact files and SHA-256 values

```text
nonlinear frozen cases/directions
2104c69ba6b21d3a3d92c7071d7f2702cb7d33f7f0e3ff17954f64c469f0c01d

audited canonical/action implementation
396c491fe51a9f5e04fa8402e2e5b16884fe23fc5057d8ded325e6064fbd3b9e

one-slab carrier/action source
ad93cdd08fabeeee56b009f23936696837c4362f88ae23f92a36d0395e61ffaf

preceding nonlinear result
a1e00071fa41f986dfaee84ea6e7689a14c50823f6c87d76889e6cb9346a7e3f.
```

Read from the frozen seed only:

- the accepted dynamic `base_old`, `base_x`, and `base_new` for each parity;
- the two calibrated base internal Jacobians;
- the unique physical-edge permutation;
- the four previously frozen zero-sum unit directions;
- `ETA=1e-4`.

Require that the odd direction is exactly the physical permutation of the
even direction.  Do not load a continuum spectrum, desired speed, chirality
sign, experimental value, full-carrier result or any fitted coefficient.

## 2. On-shell action

For `p=even,odd`, fixed logarithmic boundary values `(o,n)`, and internal
logarithms `y`, evaluate the existing action at

```text
q_old = exp(o),  x = exp(y),  q_new = exp(n).
```

Solve all 35 real equations

```text
Re gradient_internal = 0
```

with fixed-J Newton, once with the operational and once with the validation
`35 x 35` top-left blocks of the already certified canonical Jacobians.  Both
start from `log(base_x)`.  No adaptive Hessian, fitted seed, alternative root
method or failed-case amplitude reduction is allowed.

Use at most 20 iterations and backtracking factors
`1,1/2,...,2^-12`.  Accept a trial only if it remains on the certified
Lorentzian branch and satisfies

```text
||r_trial||_infinity <= (1-alpha/4) ||r||_infinity.
```

At convergence require

```text
residual infinity norm                 < 1e-55,
final fixed-J correction infinity norm < 1e-45,
maximum imaginary contamination        < 1e-70,
one negative Gram direction in all 2400 simplices,
minimum leading minor                  > 0,
minimum angle-argument modulus          > 1e-6.
```

After the final correction re-evaluate the action.  Define `W_p(o,n)` from
the pre-correction converged operational value; retain the operational versus
validation action difference and the action change under the final correction
as numerical calibration proxies.

## 3. Frozen boundary census

Let `d_i^p`, `i=1,...,4`, be the already frozen direction in parity `p`.
For every ordered pair `(i,j)`, signs `(sigma,tau)` and level

```text
i,j in {1,2,3,4},
sigma,tau in {-1,+1},
level in {1/2,1},
h = level*ETA,
```

use the four common physical Dirichlet boundary points

```text
(0,0):                 (o_base,             n_base),
(sigma h d_i,0):       (o_base+sigma h d_i, n_base),
(0,tau h d_j):         (o_base,             n_base+tau h d_j),
(sigma h d_i,tau h d_j):(o_base+sigma h d_i,n_base+tau h d_j).
```

Map these points from even to odd only through the already frozen direction
permutation.  Cache identical points, but evaluate every unique point for both
parities and both solver calibrations.

The census is exactly

```text
4 x 4 x 2 x 2 x 2 = 128 mixed rectangles,
161 unique boundary points,
644 internal solves including both parities and calibrations.
```

No direction, sign, level or failed solve may be removed after seeing output.

## 4. Difference and rectangle

At every point `z=(delta_o,delta_n)` define

```text
Delta(z) = W_odd(P z) - W_even(z).
```

For each frozen case define the operational rectangle

```text
R = Delta(sigma h d_i,tau h d_j)
  - Delta(sigma h d_i,0)
  - Delta(0,tau h d_j)
  + Delta(0,0).
```

For one boundary point define its empirical error proxy as the sum of

```text
|W_even,op-W_even,val|,
|W_odd,op-W_odd,val|,
the four absolute final-correction action changes,
1e-70.
```

Let `u` be the sum of those proxies over the four rectangle corners.  Classify
each of all 128 rectangles mechanically:

- `SEPARABLE_CONSISTENT` if `|R| <= 10*u`;
- `NONSEPARABLE` if `|R| > 100*u`;
- `OPEN` otherwise;
- `OPEN_SOLVE` if a required solve or branch control fails.

This is calibrated finite-precision evidence, not a formal interval proof.

## 5. Scaling diagnostic

For each fixed `(i,j,sigma,tau)`, if both levels are `NONSEPARABLE`, report

```text
order = log2(|R_full/R_half|).
```

Label it only diagnostically:

- `[1.5,2.5]`: `QUADRATIC_ACTION_COMPATIBLE`;
- `(2.5,3.5]`: `CUBIC_ACTION_COMPATIBLE`;
- otherwise: `OTHER_RESOLVED_ORDER`.

The label does not enter the main outcome.  In particular, a cubic-compatible
action difference is not by itself a physical interaction law.

## 6. Mechanical outcome

- any `NONSEPARABLE` rectangle:
  `BOUNDARY_COBOUNDARY_REFUTED_ON_FROZEN_RECTANGLES`;
- none nonseparable, but at least one `OPEN` or `OPEN_SOLVE`:
  `BOUNDARY_COBOUNDARY_OPEN`;
- all 128 `SEPARABLE_CONSISTENT` with every control passing:
  `BOUNDARY_COBOUNDARY_CONSISTENT_ON_FROZEN_RECTANGLES`;
- any provenance, census or implementation control failure:
  `BOUNDARY_COBOUNDARY_CONTROL_FAILED`.

The first outcome is a **DERIVED NEGATIVE** for endpoint separability of the
present two bare schedule actions on the fixed quotient.  The third is only a
**PATTERN** on the finite direction census; it is not a theorem of equivalence.

## 7. Scope and next decision

If endpoint separability is refuted, the nonlinear schedule dependence is not
removable by adding independent old/new boundary functions in the frozen
coordinates.  The next honest route is an improved/perfect action or a derived
nonlinear boundary-field redefinition; neither is thereby constructed.

If consistency survives, test the stronger homogeneous condition
`B_new=F`, `B_old=-F` before calling the schedules canonically equivalent.

This protocol does not cover the other 25 shape directions, arbitrary mixed
directions, the full 720-edge carrier, refinement, matter inhomogeneity or a
general canonical transformation mixing the two endpoints.  Only the new
targeted verifier will run; the full suite will not run.

## 8. Control-only implementation correction after the first launch

The first registered launch evaluated zero actions.  All 644 tasks stopped at
function entry because the audited angle function's standard-library symbol
`itertools.combinations` had not been imported into the new verifier's
namespace.  Independently, the direction provenance check demanded `1e-65`
normalization residuals although the frozen JSON stores the entries to roughly
60 decimal places; the observed residuals were `1.3e-60` to `2.2e-60`.

Before any Dirichlet solution or rectangle existed, the implementation was
corrected by:

1. importing `combinations`;
2. setting the direction serialization control to `1e-58`, above the analytic
   `O(30*10^-60)` decimal-storage floor.

No point, direction, sign, amplitude, solver criterion, uncertainty rule,
classification threshold or outcome rule changed.  The failed control artifact
is not scientific output and is not used as an input.
