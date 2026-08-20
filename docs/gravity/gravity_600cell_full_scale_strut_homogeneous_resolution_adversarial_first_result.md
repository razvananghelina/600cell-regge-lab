# First adversarial result: line confirmed, hostile absolute thresholds failed

Date: 2026-08-20  
Status: **DERIVED FAILURE of the preregistered control; replication not yet accepted**

## Frozen outcome

Source commit `230061e` produced

```text
HOMOGENEOUS_ADVERSARIAL_CONTROL_FAILED
6/7 tests passed
```

Artifact SHA-256:

```text
6ac12d07c7912addcd152f0de3ca019fa765fa32bc479a11ee7007ea8676e531
```

Only the targeted verifier ran.  The full suite did not run.

## What passed

- all frozen provenance controls;
- the independent full-action rebuild without artifact mutation;
- all 50 single-column-deleted D/K Gram determinants, giving the same rank
  lower bounds in both parities;
- deterministic normal-equation null extraction, independent of the primary
  symbolic generator and stored SVD candidate;
- D residuals `6.99e-83` and `6.40e-83`;
- K residuals `6.09e-83` and `5.58e-83`;
- D/K projector agreement at approximately `1.20e-85`;
- even/odd D projector agreement at `6.23e-82`;
- agreement with the primary physical ratio at `1.50e-76`.

Thus the direct matrices contain the same unique numerical line to much higher
accuracy than the preregistered `1e-30` comparison thresholds.

## What failed

The bundled hostile-control Boolean required absolute normalized residuals

```text
missing lambda > 1e-10,
wrong sign     > 1e-3.
```

Both parities instead gave

```text
missing lambda = 4.655748101552639...e-13,
wrong sign     = 2.988229695853413...e-7.
```

The controls do destroy the line: relative to the accepted D residual, their
separations are about `6.7e69` and `4.3e75`.  The preregistered absolute bounds
were nevertheless false because the Frobenius normalization is dominated by the
near-pole matrix scale.  This scaling issue should have been derived before
choosing thresholds.

The formal outcome therefore remains `CONTROL_FAILED`.  The 6/7 run cannot be
renamed a replication after seeing the data.

## Permitted next test

Any repair must be frozen before source modification and must use a previously
unevaluated precision level.  A defensible scale-free hostile criterion compares
each corrupted residual to the simultaneously reconstructed correct residual,
rather than to an arbitrary absolute number.  Because the P160 ratios are now
disclosed, a repaired run must use a fresh P200G reconstruction and retain both
absolute values in its artifact.

Until such a run succeeds, the ledger is:

- **DERIVED PRIMARY:** exact homogeneous weak-pole line and uniqueness;
- **PATTERN/STRONG COMPUTATIONAL CONFIRMATION:** the direct P160 matrices recover
  the same line;
- **OPEN:** formal adversarial replication;
- **OPEN:** the omitted pole equation and physical interpretation.

