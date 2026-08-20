# Protocol: all-schedule refined H4 internal-Jacobian census

Date: 2026-08-20

Prior-art gate commit: `4ea4430`.

This protocol is frozen before constructing any ten-by-ten matrix.

## 1. Frozen inputs

Use and hash exactly:

```text
reproducible/verify_gravity_600cell_refined_h4_stationary_fill.py
  89aab727792e20a81e7577e0425f8fa4b1e84e2a7ae66caa9e79a4aebf3581e7
reproducible/gravity_600cell_refined_h4_stationary_fill.json
  283be37bc7530a3cc4fce9e279272359f107f09fb7b1b0eaff141059bfb4e018
docs/gravity/gravity_600cell_refined_h4_stationary_fill_result.md
  899143cfcf75ce08a8dd6daca776cd0e38e02e38593f4a171d3d3f851dae7d91
```

Load definitions only from the verifier's AST; do not execute its top-level
body.  Reconstruct its exact rank geometry and parse all 24 certified abstract
simplex/triangle incidence records from the frozen artifact.  No rounded
spatial coordinate is an action input.

Require the upstream outcome
`REFINED_H4_INDUCED_FILL_OFF_SHELL`, `12/12`, 24 schedules, 240 internal
entries, 96 certified vertical nonzeros, zero certified cross nonzeros and one
residual-vector class.

## 2. Matrix definition

Keep all twelve boundary coordinates fixed.  Let

```text
y=(log(d_01/d_01^0),...,log(d_23/d_23^0),
   log(rho_0/rho_0^0),...,log(rho_3/rho_3^0)).
```

For schedule `sigma`, let `G_sigma(y)` be the ten **total orbit** log
derivatives of the complete gravity-plus-dust action.  Define

```text
H_sigma = dG_sigma/dy at y=0.
```

This total-orbit convention is required because `H_sigma` is then the true
internal log-coordinate action Hessian and must be symmetric.  Per-edge row
normalization may be reported for conditioning, but must not replace the
matrix in symmetry, inertia or nullity claims.

## 3. High-precision derivative and error envelope

At 100 decimal digits, compute centred gradient differences at

```text
h0=1e-10, h1=5e-11, h2=2.5e-11.
```

For every column form

```text
D(h)=(G(+h e_j)-G(-h e_j))/(2h),
H_A=(4D(h1)-D(h0))/3,
H_B=(4D(h2)-D(h1))/3.
```

Use the symmetrized real part of `H_B` only after recording the raw complex
matrix.  Define

```text
matrix_scale=max(1,max_abs(H_B)),
entry_error=100*max_abs(H_A-H_B)+1e-50*matrix_scale,
spectral_error=10*entry_error.
```

Mandatory controls:

- raw imaginary entries and raw antisymmetry are each no larger than
  `entry_error`;
- repeat `H_B` at 140 decimal digits and require its difference from the
  100-digit value to lie within `entry_error`;
- on the lexicographically first schedule and its reverse, independently
  reproduce `v^T H v` from Richardson-extrapolated centred second differences
  of the action for three frozen directions: `cross_01`, `rho_0`, and the
  normalized induced common-lapse tangent.  Maximum relative disagreement
  must be below `1e-28`.

## 4. Rank, inertia and schedule comparison

Use high-precision symmetric eigendecomposition.  An eigenvalue is certified
nonzero only if its absolute value exceeds `spectral_error`; otherwise it is
unresolved/zero-compatible.  Report rank, inertia and all ten eigenvalues for
every schedule.  Do not use a hand-selected relative SVD cutoff.

Two matrices are the same class only when their maximum entrywise difference
is at most the larger of their two entry envelopes.  Report the number of
matrix classes, complete class membership and all 12 time-reversal-pair
differences.

The exact induced-lapse tangent at the inherited fill is

```text
t_cross(rs)=-rho0/d_rs,
t_rho(r)=1.
```

Normalize it in Euclidean log-coordinate norm.  Report `||H t||`, its Rayleigh
quotient and its overlap with every zero-compatible eigenvector.  Call it a
null/gauge candidate only when `||H t|| <= 10*spectral_error` and the largest
overlap is at least `1-1e-20` for every schedule.

For conditioning only, compute the un-applied linear Newton proposal from the
certified eigendecomposition and report its norm and linear residual.  Do not
update a coordinate in this mission.

## 5. Frozen outcomes

Use the first applicable outcome:

1. `REFINED_H4_INTERNAL_JACOBIAN_CONTROL_FAILED` if any provenance,
   reconstruction, derivative, symmetry, precision, directional or
   time-reversal control fails.
2. `REFINED_H4_INTERNAL_JACOBIAN_MIXED_RANK` if schedules have different
   certified ranks.
3. `REFINED_H4_INTERNAL_JACOBIAN_FULL_RANK_SINGLE_CLASS` if every rank is ten
   and there is one matrix class.
4. `REFINED_H4_INTERNAL_JACOBIAN_FULL_RANK_MULTIPLE_CLASSES` if every rank is
   ten and more than one class remains.
5. `REFINED_H4_INTERNAL_JACOBIAN_INDUCED_LAPSE_NULL` if every rank is nine and
   the sole zero-compatible direction passes the induced-lapse test.
6. `REFINED_H4_INTERNAL_JACOBIAN_RANK_DEFICIENT_OTHER` for any other common
   deficient rank or unresolved null geometry.

Full rank licenses an ungauged square ten-equation local solve.  Outcome 5
licenses a separately preregistered dust-clock gauge.  Multiple classes
require solving every class; they do not by themselves establish different
on-shell dynamics.

## 6. Exclusions and deliverables

Forbidden: applying Newton, finding a root, eliminating an internal variable,
forming an effective boundary Hessian, diagonalizing a boundary/mode operator,
selecting a schedule, or comparing with `c`, `G`, Planck or continuum targets.

Deliver a registered verifier, deterministic JSON, result note, identical
targeted rerun and static registry audit.  Do not run the full suite.
