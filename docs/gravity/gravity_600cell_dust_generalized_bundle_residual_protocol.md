# Protocol: residual-certified finite-family bundle comparison

Date: 2026-08-18

Prior-art/framing gate commit: `7e815de`.

## 1. Frozen question

For the already accepted three-slab 600-cell Regge--dust trajectory, determine
whether the rank-15 generalized negative-mode projectors at the old and
shifted centered times are separated in every one of the four disclosed
finite-difference schedules, after certifying each schedule's source ball and
eigensolver residual at high precision.

The result is explicitly a statement about the complete frozen four-schedule
family.  It is not an analytic enclosure of the exact continuum Hessian and
cannot prove exact equality of two subspaces.

No reduced recurrence, product, characteristic root, wave speed, mass or
particle interpretation may be computed in this mission.

## 2. Frozen inputs and arithmetic

Use exactly the inputs and geometry sources already frozen by
`verify_gravity_600cell_dust_generalized_bundle_direct_precision.py`.  Require
that verifier to reproduce `13/13` and the committed diagnostic artifact hash

```text
8ded406366dbf291da02dfbf995c4e37036cc6ce745d9240d14905664ba6042a.
```

Replay the direct slab construction at `100` decimal digits with Flint ball
arithmetic at `80` decimal digits.  Preserve both parities, sectors `4,5`, and
all four already disclosed derivative schedules.  Convert each Flint entry to
an `mpmath` midpoint by its decimal ball midpoint, not through binary64.  Carry
the Frobenius norm of the Flint component radii as the schedule-local source
error.

Use an arithmetic floor of

```text
1e-70 * max(1, relevant operator norm)
```

for each high-precision construction.  Do not add the earlier binary64
`1000*eps` floor and do not collapse the four derivative schedules into one
unstructured family-norm ball.

## 3. Deterministic carriers

### 3.1 Conformal incidence image

Construct the `30 x 120` compressed incidence matrix directly from the exact
integer incidence relation and the high-precision one-dimensional sector
basis.  Scan its columns lexicographically and accept a column after modified
Gram--Schmidt exactly when its residual norm exceeds

```text
1e-50 * max(1, ||compressed incidence||).
```

Stop at five accepted columns.  Reject the calculation unless the rank is
exactly `5`, the orthonormality residual is below `1e-60`, and the compressed
incidence residual outside their span is below `1e-60` relative to scale.
This construction uses geometry only and is shared by the old and shifted
times.

### 3.2 Action-selected shape nullspace

For every schedule-local midpoint `M`, form `R=U* M`.  Select five pivot
columns by the same lexicographic modified Gram--Schmidt rule.  With pivot
block `R_P` and free block `R_F`, construct the nullspace graph

```text
W_raw[pivots,:] = -R_P^{-1} R_F,
W_raw[free,:]   = I_25.
```

Orthonormalize it using the positive eigendecomposition of
`W_raw* W_raw`.  Reject unless there are five pivots, `W` has 25 columns,
`||U* M W||` is below `1e-60` relative to scale, and `||W*W-I||<1e-60`.

For the source ball, let `epsilon_M` and `epsilon_V` be the Flint radius norms
plus the arithmetic floor.  Bound the shape rotation with the already used
rank-gap formula

```text
eta_S = 2 epsilon_M / (sigma_min(U* M) - 2 epsilon_M),
```

rejecting the cell unless the denominator is positive.

## 4. Hermitian-definite pencil and residual certificate

Restrict and Hermitian-symmetrize

```text
A = -W* V W,       B = -W* M W.
```

Propagate the source balls by

```text
epsilon_B = epsilon_M + 2 eta_S (||M|| + epsilon_M) + floor,
epsilon_A = epsilon_V + 2 eta_S (||V|| + epsilon_V) + floor.
```

Require `B>0` after subtracting `epsilon_B`.  Let `B=L L*` be its lower
Cholesky factor and form

```text
H = L^{-1} A L^{-*}.
```

Diagonalize `H` at high precision, ordered increasingly, and require exactly
15 negative and 10 positive eigenvalues and a positive gap

```text
g = lambda_16 - lambda_15.
```

With `Q=[Q_- Q_+]`, compute the actual midpoint off-diagonal residual

```text
R_H = Q_+* H Q_-.
```

Its a posteriori Davis--Kahan contribution is

```text
eta_res = sin(0.5 atan(2 ||R_H|| / g)).
```

The schedule-local Flint source-ball contribution retains the conservative
whole-pencil bound, now without derivative-family variation:

```text
epsilon_pencil = epsilon_A / b_lower
                + ||A|| epsilon_B / (b_min b_lower) + floor,

eta_source = 2 epsilon_pencil / (g - 2 epsilon_pencil).
```

Reject unless `g>2 epsilon_pencil`.  Lift the lower cluster back to the
30-dimensional Euclidean edge carrier and form its basis-independent
orthogonal projector

```text
Z = W L^{-*} Q_-,
P = Z (Z*Z)^{-1} Z*.
```

The complete per-schedule projector error is

```text
eta_P = 2 eta_S
      + sqrt(b_max / b_lower) (eta_source + eta_res)
      + epsilon_B / b_lower
      + floor.
```

Record every additive contribution, the residual, gap, source radii,
projector Hermiticity/idempotence residual and `eta_P`.

## 5. Controls against the committed calculation

For each of the 32 high-precision projectors, compare its binary conversion
with the corresponding projector reconstructed by the committed direct
verifier.  It must fall inside the latter's already committed broad projector
error.  For each of the 16 matched old/shifted cells, the new distance must
also lie inside the old distance plus its committed comparison error.

Failure of either control is `RESIDUAL_BUNDLE_CONTROL_FAILED`; no scientific
classification survives it.

## 6. Complete comparisons and fixed labels

At each time/parity/sector compare all `C(4,2)=6` schedule pairs, for exactly

```text
2 times * 2 parities * 2 sectors * 6 = 48
```

within-time comparisons.  Compare every old schedule with every shifted
schedule, for exactly

```text
2 parities * 2 sectors * 4 * 4 = 64
```

cross-time comparisons.

For a projector distance `d` and complete error
`e=eta_left+eta_right+floor`, assign only:

```text
d <= 10 e       ZERO_CONSISTENT,
d >  100 e      ROTATION_RESOLVED,
otherwise       OPEN.
```

No matched schedule may receive preferential weight.

## 7. Frozen outcome hierarchy

Apply the first matching branch:

1. any provenance, rank, positivity, gap, residual, overlap, count or finite
   arithmetic control fails:
   `RESIDUAL_BUNDLE_CONTROL_FAILED`;
2. any of the 48 within-time comparisons is not `ZERO_CONSISTENT`:
   `RESIDUAL_FINITE_FAMILY_SCHEME_DEPENDENT`;
3. all 64 old/shifted comparisons are `ROTATION_RESOLVED`:
   `RESIDUAL_FINITE_FAMILY_ROTATION_RESOLVED`;
4. all 64 are `ZERO_CONSISTENT`:
   `RESIDUAL_FINITE_FAMILY_ZERO_CONSISTENT`;
5. otherwise:
   `RESIDUAL_FINITE_FAMILY_ROTATION_OPEN`.

The branch names do not assert analytic continuum rotation or a physical
connection.  Even outcome 3 remains **DERIVED COMPUTATIONAL, conditional on
the frozen derivative family**.

## 8. Required artifact and verifier checks

Write one deterministic JSON artifact containing:

- all input hashes and prior-art/protocol commits;
- 32 carrier/pencil/projector records;
- 48 within-time records and label counts;
- 64 old/shifted records and label counts;
- complete controls, outcome and status ledger.

Register the verifier in `reproducible/run_all.py` before its first scientific
execution.  Execute it twice and require byte-identical artifacts.  Run no
unrelated verifier and do not run the full suite.
