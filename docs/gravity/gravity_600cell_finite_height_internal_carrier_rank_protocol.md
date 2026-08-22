# Preregistration: finite-height internal carrier rank

Date: 2026-08-22.

Prior-art/framing commit: `3da6ec6`.

Status: **FROZEN BEFORE THE FIRST EVALUATION OF AN INTERNAL-CARRIER RANK,
SINGULAR VALUE, MINOR OR NULL VECTOR AT THE FINITE-HEIGHT BACKGROUND.**

## 1. Frozen inputs

The verifier must reject any hash mismatch.

```text
docs/gravity/gravity_600cell_finite_height_internal_carrier_rank_prior_art.md
  f3d04db084a63944e6963687747dcbe510d910d12f8b65b612ae50f6a1d89696

reproducible/gravity_600cell_finite_height_carrier_quadratic.json
  0ec142bfc68d04498992a6cdba7437933560b860244573d187cb6e018ece78f9

reproducible/gravity_600cell_finite_height_carrier_quadratic_adversarial.json
  54915cf364c36af6bbc8e1dbd36433079269d293453478bfdf589e547d462ad6

reproducible/gravity_600cell_finite_height_carrier_quadratic_adversarial_matrices.npy
  8a3ea0c3b8ee720d8ffdf07e7486aefdd0247ca1cfdbeb99f443091376f31729

reproducible/verify_gravity_600cell_finite_height_carrier_quadratic_adversarial.py
  8d37012f556ce5be0bb863ad12d4572d197c90a5b96974912e81a98c1956a8f8

reproducible/verify_gravity_600cell_dust_full_boundary_tangent.py
  c8662bb0835865aac6696fc3f474ed668fed3fe393b9c32a59e709a984c35571

reproducible/verify_gravity_global_regge_orbits.py
  ad93cdd08fabeeee56b009f23936696837c4362f88ae23f92a36d0395e61ffaf

reproducible/gravity_600cell_full_scale_strut_symbolic_gap_resolution.json
  ea2c52f0cd227516734defc509330e528b140f71bfd0f50e87036f3fa9832179
```

Require the two accepted quadratic artifacts to retain their `22/22` primary
and `18/18` adversarial outcomes, the exact carrier to retain rank 240, and
the primary internal gradients to remain below `1e-25` for both parities.

The finite-height background is read from the frozen primary artifact and
independently checked against the exact homogeneous formulas before use.

## 2. Exact matrix and row split

For each parity reconstruct the complete `2,400`-simplex slab, its `95`
free regular `2T` edge orbits and the high-precision identity-row Hessian
kernels.  Let

```text
T_p = H_p[active,active] G_p       shape 1560 x 240,
R_p = H_p[internal,active] G_p     shape  840 x 240.
```

`R_p` is not symmetrized.  Its rows are equations, not a bilinear form.
Classify two nested maps:

```text
R_diag : the 720 oriented cross-diagonal rows,
R_full : all 720 diagonal plus 120 pole rows.
```

Identify row kinds from physical edge labels.  Do not assume the pole orbits
occur first or last within the 35 internal orbit types.

Compose each response directly from the orbit kernel and sparse carrier:

```text
H[(r,a),(c,b)] = K[r,c,a^-1*b],
T[(r,a),j]     = sum_(c,b) H[(r,a),(c,b)] G[(c,b),j].
```

Never materialize a complete `2280 x 2280` Hessian.  Keep all accumulation in
`mpmath` at 180 decimal digits.  The primary quadratic matrices and any
continuum data remain unread until every rank classification is fixed in
memory.

## 3. Frozen derivative hierarchy and reality controls

Use the already accepted derivative steps

```text
h0 = 1e-25,
h1 = 5e-26,
h2 = 2.5e-26,
h3 = 1.25e-26.
```

For either `R_diag` or `R_full`, construct

```text
M01=(4 M1-M0)/3,
M12=(4 M2-M1)/3,
M23=(4 M3-M2)/3.
```

Every base and displaced simplex must retain Lorentzian inertia `(3,1)`, the
accepted leading-minor and angle-argument margins, and the local derivative
tables must retain their accepted step hierarchy.  Require the final physical
response kernels to have imaginary Frobenius residue below `1e-140` relative
to `max(1,||M12||_F)`.  Raw Lorentzian boost-angle imaginary components are
diagnostic, not contamination.

## 4. Seven target-free minimal sectors

Rebuild the seven deterministic minimal sectors of the regular `2T` action at
180 digits.  Require

```text
irrep dimensions       [1,1,1,2,2,2,3],
isotypic dimensions    [1,1,1,4,4,4,9],
sum d^2                24,
```

and every orthogonality, central-splitter, right-action leakage and conjugate-
pair residual below `1e-140`.

Compress `R_diag` and `R_full` independently.  In a sector of dimension `d`
their shapes are

```text
R_diag,s : (30d) x (10d),
R_full,s : (35d) x (10d).
```

The ten source copies are ordered as five vertex-scale and five strut
regular orbits using the pole identity, not by a guessed vertex ordering.

## 5. Frozen scaling and singular-value classifier

For each parity, sector and row scope, derive two positive source scalars from
`M12`:

```text
s_scale = 1/max(1,||M12[:,scale]||_F),
s_strut = 1/max(1,||M12[:,strut]||_F).
```

Apply the same block-diagonal source scaling to `M01`, `M12` and `M23`.
It is invertible and symmetry-preserving, so it cannot change a nullity.

Set

```text
N      = max(1,||M12||_F),
e_step = max(||M01-M12||_F,||M12-M23||_F)/N,
e_total = e_step+1e-135.
```

At 180 digits compute every singular value of all three scaled matrices.
Using the common normalized threshold, classify each as

```text
ZERO     sigma/N <= 10 e_total,
NONZERO  sigma/N >  100 e_total,
OPEN     otherwise.
```

A block rank is resolved only if no level has an `OPEN` singular value and
all three Richardson levels give the same nullity.  Record every singular
value and label.  No floating binary64 singular value can override this
classification.

## 6. Direct minor certificate

For every resolved block with rank `r>0`, convert only `M12` to binary64 and
select an `r x r` submatrix by this frozen deterministic procedure:

1. pivoted QR of `M12` selects `r` columns;
2. pivoted QR of the selected-column transpose selects `r` rows;
3. ties inherit SciPy's deterministic index order.

The columns should be all `10d` columns for a full-column-rank block, but the
general rule also covers a positive nullity.  Evaluate the same selected
minor directly in 180-digit arithmetic at `M01`, `M12` and `M23`.  With

```text
e_det = max pairwise determinant difference + 1e-120,
```

require

```text
min determinant modulus > 100 e_det.
```

Together, the singular upper bound and direct-minor lower bound certify the
reported numerical rank.  A failed minor is `NUMERICALLY_OPEN`, not evidence
for an extra null vector.  The row and column indices are part of the artifact.

## 7. Global reconstruction and kernel comparison

After all sector ranks are fixed in memory, retain the global `M12` matrices
in the common physical data order

```text
sigma_0,...,sigma_119,c_0,...,c_119.
```

The global nullity must equal

```text
sum_s d_s * nullity_s.
```

If the complete nullity is positive, compute its right-kernel projector from
the global binary64 matrix using exactly the sector-certified nullity.  The
projector uncertainty is the maximum within-parity difference among the
`M01`, `M12` and `M23` projectors plus

```text
500*eps_machine*max(matrix dimensions).
```

Even and odd kernels are `AGREE` only if their projector two-norm difference
is at most ten times that uncertainty, `DEPENDENT` only if it exceeds one
hundred times it, and otherwise `OPEN`.  A zero kernel needs no projector
comparison.

## 8. Accepted-quadratic closure control

Only after the rank census is frozen in memory, load the accepted adversarial
quadratic matrices.  Reconstruct

```text
Q'_p = (G_p^* T_p + T_p^* G_p)/2
```

from the response kernel and require relative Frobenius agreement with the
accepted `R12` quadratic form below `1e-10` for both parities.  This is a
group-product/carrier closure control, not a rank target.

## 9. Hostile controls

1. Under the actual `e_total`, a zero matrix must have full nullity and an
   embedded identity must have nullity zero.
2. Add `+1/10` to the source-scale coefficient on all 24 translates of the
   lexicographically first oriented-diagonal orbit.  At least one sector and
   the global response must change by more than `1e-12` relatively.
3. Add a conceptual `+1/10` to the matching 24 diagonal entries of the
   internal Hessian.  Its induced response must change by more than `1e-12`
   relatively.
4. Reversing `a^-1*b` to `b*a^-1` is recorded as a convention diagnostic but
   is not an acceptance gate.
5. A deterministic positive upper-triangular source change must preserve all
   resolved nullities.

No failed hostile control may be removed or replaced after execution.

## 10. Outcome hierarchy

### `FINITE_HEIGHT_INTERNAL_CARRIER_CONTROL_FAILED`

Use if provenance, background, branch, group, carrier, reality, closure or a
hostile control fails.

### `FINITE_HEIGHT_INTERNAL_CARRIER_NUMERICALLY_OPEN`

Use if any actual sector has an open singular value, inconsistent Richardson
nullity, failed direct minor, or an unresolved positive-kernel parity
comparison.

### `FINITE_HEIGHT_INTERNAL_CARRIER_FULL_COLUMN_RANK_PRIMARY`

Use only if every complete sector in both parities has nullity zero with both
the singular and direct-minor certificates.  The diagonal-only nullities are
reported but do not weaken this complete kill result.

### `FINITE_HEIGHT_INTERNAL_CARRIER_KERNEL_SELECTED_PRIMARY`

Use only if both complete maps have the same positive global nullity and
their physical-data kernel projectors `AGREE`.

### `FINITE_HEIGHT_INTERNAL_CARRIER_SCHEDULE_DEPENDENT_PRIMARY`

Use only if complete nullities differ, or equal positive-nullity projectors
are classified `DEPENDENT`.

Every material `PRIMARY` outcome requires a separately preregistered,
mechanically different replication before consolidation.  No primary result
alone is promoted to a programme-level theorem.

## 11. Interpretation firewall

Full column rank is a bounded negative for this exact geometry-selected
carrier at this exact finite-height slab.  It does not refute the full Regge
boundary space, a nonlocal/higher-rank carrier, refinement or general
relativity.

A selected kernel supplies candidates only.  It does not establish gauge
invariance, a graviton, a wave equation, stability, a continuum limit, a
physical tick, `c`, `G`, Planck units, particle masses or Standard-Model
physics.

Run only the new verifier and static registry checks.  Do not run the full
suite.
