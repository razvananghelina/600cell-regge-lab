# Result: calibrated complete-action solver continuation

Date: 2026-08-16

Machine-result commit: `94675be`

Machine-result SHA-256:
`14bac60c3c561db74290105d4049f15bd44ff056c597a72882aa16a4ee6d7719`

Prior-art gate: `6b7f9e4`

Base preregistration: `4b6b10c`

Failed calibration result: `87c289f`

Calibrated continuation preregistration: `cf27934`

Continuation implementation: `5346d3a`

Multiprocessing-only repair: `5d058bb`

Frozen upstream numerical boundary: `64a13f6`

Targeted solver verifier:
`reproducible/verify_gravity_600cell_dust_action_solver_repair_continuation.py`

Post-result auditor:
`reproducible/verify_gravity_600cell_dust_action_solver_repair_continuation_result.py`

Machine-readable result:
`reproducible/gravity_600cell_dust_action_solver_repair_continuation.json`

## Executive verdict

The calibrated 100-decimal continuation completed all 80 preregistered
states.  The targeted solver ended with `11/11` implementation checks and the
post-result artifact auditor ended with `14/14` checks.

The global outcome is

```text
ACTION_REDUCED_SCAN_NUMERICALLY_UNRESOLVED
```

The complete result is:

```text
80 frozen grid states
63 independently validated transverse roots
63/63 reduced scalars resolved nonzero and positive
 9 states with no robust action descent
 8 states at the frozen 12-iteration limit
 0 validated scalar zeros
 0 forced bisection midpoints
```

Therefore:

- **DERIVED COMPUTATIONAL:** 63 frozen states solve all 34 transverse
  complete-action equations.
- **DERIVED COMPUTATIONAL NEGATIVE:** none of those 63 states solves the full
  35-component internal stationary equation, because its remaining scalar is
  nonzero by at least `1.83e14` estimated errors.
- **PATTERN-informed:** all 63 resolved scalars have the same positive sign.
  The sign had been inspected before the calibrated preregistration, so its
  provenance cannot be upgraded.
- **OPEN NUMERICALLY:** 17 states have no validated transverse root under the
  frozen solver budget.  Solver failure is not root nonexistence.
- **OPEN:** no continuous-interval no-root theorem, second slab, full carrier,
  physical tick, causal speed, mass scale or Planck unit follows.

This is useful, but it is not a discovered tick.  It says that the tested
prescribed-boundary one-slab continuations do not yield a full stationary
internal configuration at any of the 63 validated grid points.

## Status ledger

| Item | Status | Evidence |
|---|---|---|
| Published complete-action control | **DERIVED CONTROL** | retained `14/14` upstream certificates |
| Even derivative calibration | **DERIVED CONTROL** | every amended and inherited gate passes |
| Odd derivative calibration | **DERIVED CONTROL** | every amended and inherited gate passes |
| Frozen target carrier | **DERIVED COUNT** | `2 x 4 x 2 x 5 = 80` unique states |
| Complete target evaluation | **DERIVED COMPUTATIONAL** | all `16/16` signed cases completed |
| Validated transverse roots | **DERIVED COMPUTATIONAL** | `63/80` |
| No robust descent | **OPEN NUMERICALLY** | `9/80` |
| Iteration limit | **OPEN NUMERICALLY** | `8/80` |
| Validated scalar zeros | **DERIVED NEGATIVE COUNT** | `0/63` |
| Validated nonzero scalars | **DERIVED COMPUTATIONAL** | `63/63` |
| Common positive sign | **PATTERN-informed** | `63/63`; sign inspected upstream |
| Forced sign-change brackets | **DERIVED NEGATIVE COUNT** | `0` |
| Signed-case hits | **DERIVED NEGATIVE COUNT** | `0/16` |
| Direction/parity-pair hits | **DERIVED NEGATIVE COUNT** | `0/8` |
| Phase-contrast hits | **DERIVED NEGATIVE COUNT** | `0/4` |
| Fully resolved signed cases | **DERIVED COMPUTATIONAL LOCAL** | `6/16` have no grid hit |
| Numerically open signed cases | **OPEN NUMERICALLY** | `10/16` |
| Continuous interval | **OPEN** | five-point grid is not an interval proof |
| Second slab / momentum matching | **NOT TESTED** | outside frozen carrier |
| Full 600-cell carrier | **NOT TESTED** | order-24 invariant carrier only |
| Physical time, inertia, masses, `c`, Planck units | **NOT DERIVED** | no multi-step dynamics or dimensional scale |
| External novelty | **OPEN** | focused search, not a dedicated review |

No full-suite run was performed.  This obeys the user's instruction to run
only the verifier relevant to the active calculation.

## 1. Exact equation and why the scalar matters

For each frozen state the verifier solves

```text
F(t,z,b) = Q^T E_action(t,z,b) = 0,       z in R^34,
```

where the 35-vector is

```text
E_action,i = (1/24) d S_total / d log(u_i)
```

for the complete Lorentzian Regge-plus-dust action.  After a transverse root
is independently validated, it computes

```text
g(t,z,b) = w(t)^T E_action(t,z,b).
```

The post-result auditor checked the possible framing failure directly.  The
upstream boundary-Legendre construction distinguishes:

```text
35 internal variables       -> internal equations of motion
30 final-boundary variables -> post-momenta
30 old-boundary derivatives -> pre-momenta
```

Thus `g` is not an endpoint momentum that was incorrectly required to vanish.
It is the last component of the internal equation.

The auditor also reconstructed `[Q w(t)]` on every frozen `t`.  It has rank 35
throughout, with minimum singular value

```text
0.999998258139.
```

Hence imposing both `F=0` and `g=0` is a well-conditioned decomposition of
the full internal condition `E_action=0`.  Conversely, a validated `F=0` with
nonzero `g` is genuinely not a stationary internal slab on this carrier.

## 2. Numerical strength of the 63 validated results

Every validation gate passes for all 63 states:

- disjoint operational and validation action rows agree;
- every branch audit passes;
- all imaginary contaminations remain below `1.05e-78`, versus the frozen
  `1e-70` threshold;
- every raw and preconditioned transverse zero gate passes;
- the maximum transverse-norm/error ratio is `5.33e-4`;
- the maximum preconditioned-norm/error ratio is `0.966`.

The reduced scalar range is

```text
minimum g = 8.421488290121002e-17   (odd, d1, sign +, t=+0.10)
maximum g = 1.037882148866696e-16   (even, d4, sign +, t=-0.10)
```

Its estimated errors lie between `4.16e-31` and `4.59e-31`, giving

```text
min |g|/error = 1.8334e14
max |g|/error = 2.4972e14.
```

The statement `g != 0` is therefore not marginal numerical thresholding at
the 63 validated states.

## 3. The 17 unresolved states

The unresolved cases are not silently counted as negative results.

Nine stop because no damping in

```text
1, 1/2, ..., 1/1024
```

produces a robust reduction of the complete-action merit.  Their final
merit/error ratios range from `10.38` to `19.95`.

Eight accept twelve iterations but remain nonzero.  Their final merit/error
ratios range from `388` to `31596`.  They are numerically farther from the
zero band than the no-descent states, but the frozen iteration limit still
forbids declaring root absence.

The outcome distribution is:

| Parity | Validated | No descent | Iteration limit |
|---|---:|---:|---:|
| even | 35 | 5 | 0 |
| odd | 28 | 4 | 8 |
| total | 63 | 9 | 8 |

By direction the validated/unresolved counts are `17/3`, `16/4`, `15/5`,
and `15/5`.  The central grid value `t=0` is numerically hardest: `10/16`
validate there, compared with `14/16` at both endpoints.

## 4. Complete signed-case table

| Parity | Direction | Sign | Validated | No descent | Limit | Case outcome |
|---|---:|---:|---:|---:|---:|---|
| even | 1 | - | 5 | 0 | 0 | no grid stationary point |
| even | 1 | + | 5 | 0 | 0 | no grid stationary point |
| even | 2 | - | 5 | 0 | 0 | no grid stationary point |
| even | 2 | + | 3 | 2 | 0 | numerically unresolved |
| even | 3 | - | 5 | 0 | 0 | no grid stationary point |
| even | 3 | + | 5 | 0 | 0 | no grid stationary point |
| even | 4 | - | 4 | 1 | 0 | numerically unresolved |
| even | 4 | + | 3 | 2 | 0 | numerically unresolved |
| odd | 1 | - | 4 | 1 | 0 | numerically unresolved |
| odd | 1 | + | 3 | 0 | 2 | numerically unresolved |
| odd | 2 | - | 5 | 0 | 0 | no grid stationary point |
| odd | 2 | + | 3 | 1 | 1 | numerically unresolved |
| odd | 3 | - | 4 | 0 | 1 | numerically unresolved |
| odd | 3 | + | 1 | 2 | 2 | numerically unresolved |
| odd | 4 | - | 4 | 0 | 1 | numerically unresolved |
| odd | 4 | + | 4 | 0 | 1 | numerically unresolved |

`No grid stationary point` means only the five frozen `t` values in that
signed case.  It is not a continuous no-root statement.

## 5. Framing attack: what this experiment actually measures

The internal-versus-boundary audit rescues the mathematical stationary
condition, but a second framing limitation remains.

The 80 tests are boundary-value problems: an old 600-cell boundary and a
deformed final 600-cell boundary are prescribed, and the 35 bulk variables of
the intervening slab are varied.  A generic prescribed pair of boundaries
need not admit a stationary Regge filling.  Therefore failure to find one is
not by itself failure of discrete dynamics.

The repository's published symmetric control is important here.  It already
certifies that the undeformed De Felice--Fabri dust sandwich solves all 35
complete one-slab internal equations in both schedule parities.  So:

- **DERIVED:** the action and internal stationary equations can possess a
  solution at the symmetric control;
- **DERIVED COMPUTATIONAL LOCAL:** the 63 validated deformed boundary states
  are not solutions;
- **STRUCTURAL:** the common positive `g` can be read as a one-sided internal
  action slope for those prescribed boundary data;
- **OPEN:** whether it is a physical obstruction, a pseudo-constraint fixing
  admissible boundary data, or a restriction artifact of the order-24
  carrier.

Calling the positive `g` a "tick", "time arrow" or expansion law would be
incorrect.  It is a residual action slope, not an elapsed time or a solved
next state.

## 6. Post-result prior-art update

The focused post-result search weakens any novelty claim.

De Felice and Fabri already report that 600-cell Regge evolutions reach a
stopping point while the spatial volume is nonzero, interpreted as a
causality-breaking point.  Their generalized calculation allows more free
variables but reports the same qualitative issue:

- A. De Felice and E. Fabri, [*The Friedmann universe of dust by Regge
  Calculus: study of its ending point*](https://arxiv.org/abs/gr-qc/0009093).
- A. De Felice and E. Fabri, [*Singularities of the closed RW metric in Regge
  Calculus: a generalized evolution of the
  600-cell*](https://arxiv.org/abs/gr-qc/0106077).

Bahr and Dittrich show that curvature can break exact discrete gauge
symmetries and replace constraints by pseudo-constraints that depend on
next-step lapse/shift-like data.  This is a plausible structural category for
an internal residual, but it does not identify our `g` without a dedicated
canonical derivation:

- B. Bahr and B. Dittrich, [*(Broken) Gauge Symmetries and Constraints in
  Regge Calculus*](https://arxiv.org/abs/0905.1670).

Dittrich and Hoehn formulate canonical simplicial evolution in which data
introduced by one move may be fixed by constraints from later moves.  This
supports testing slab composition rather than interpreting one prescribed
boundary problem as a complete evolution law:

- B. Dittrich and P. A. Hoehn, [*Canonical simplicial
  gravity*](https://arxiv.org/abs/1108.1974).

Jercher and Steinhaus find in a different Lorentzian Regge cosmology that
causal regularity and inequalities between geometric and matter boundary data
control whether a height variable has a solution.  Their result reinforces
the need to audit causal cell type and admissible boundary data:

- A. F. Jercher and S. Steinhaus, [*Cosmology in Lorentzian Regge calculus:
  causality violations, massless scalar field and discrete
  dynamics*](https://arxiv.org/abs/2312.11639).

None of these papers supplies the repository's 63-state scalar calculation.
Conversely, this calculation does not establish a new theorem about Regge
causality or pseudo-constraints.  External novelty remains **OPEN**.

## 7. Next logical test

Do not spend the next mission merely increasing the iteration limit on the 17
states.  That could complete this finite scan, but even a `80/80` positive
sign would still not produce evolution.

The higher-value next test is a **two-slab homogeneous control** using the same
complete action:

1. glue two certified 600-cell dust slabs along a shared 600-cell slice;
2. fix only the two outer boundary slices;
3. vary every symmetry-reduced bulk variable and the shared-slice scale;
4. require the shared-boundary post-momentum of the first slab to equal the
   pre-momentum of the second;
5. first reproduce the time-symmetric published control;
6. only then perturb the outer boundaries and ask whether an expanding and a
   contracting branch are selected;
7. compare the resulting discrete scale sequence with the Friedmann dust
   control without fitting a time step.

This is the first test that can distinguish:

- a genuine discrete evolution law;
- arbitrary incompatible boundary data;
- a causal/pseudo-constraint obstruction;
- an artifact of the one-slab order-24 restriction.

It still will not derive `c`, Planck time or particle masses.  Those require a
multi-step causal propagation law and at least one dimensional scale after
the dynamics itself survives.

