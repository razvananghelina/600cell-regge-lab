# Repair protocol: scale-free hostile controls at fresh P200G

Date: 2026-08-20  
Status: **preregistered after disclosed P160 control failure; not pristine blind evidence**

## Why a repair is permitted

The first adversarial artifact and result note in commit `5d43620` are immutable
inputs to this repair.  The P160 line, all 50 rank certificates and all line
comparison tests passed, but the formal outcome was
`HOMOGENEOUS_ADVERSARIAL_CONTROL_FAILED` because two arbitrary absolute hostile
thresholds were incompatible with the Frobenius scale of the near-pole matrix.

This repair does not relabel that run.  It evaluates a previously unrun P200G
matrix reconstruction and writes a new artifact.  Because the P160 values are
known, this is a disclosed post-failure validation, not a blind preregistration.

## Frozen fresh computation

Use:

```text
mpmath dps = 200
Arb dps    = 180
conversion digits = 185
```

Retain the accepted P200G derivative variants
`1e-40,1e-30,3e-40,3e-30`.  Re-execute the full-action builder, rebuild the two
homogeneous D/K matrices and extract the line with the same fixed-last-component
normal-equation method.  Retain every original rank, residual, group-spread,
D/K-projector, parity-projector and primary-ratio criterion.

Pin the failed P160 artifact and result note by SHA-256.  Write the repaired
result to

```text
gravity_600cell_full_scale_strut_homogeneous_resolution_adversarial_p200g.json
```

without overwriting the failure artifact.

## Scale-free hostile criteria

For each parity, let `r` be the normalized residual of the directly extracted D
line, and let `r_missing`, `r_sign` be the residuals of the missing-lambda and
wrong-sign directions.  Require

```text
r < 1e-30,
r_missing > 1e-20,
r_sign    > 1e-20,
r_missing/r > 1e40,
r_sign/r    > 1e40.
```

The absolute floor only rejects numerical zero; the evidential criterion is the
dimensionless separation ratio.  The `1e40` ratio is intentionally more than 29
orders below the disclosed P160 separations, so this test asks for a qualitative
destruction of the line rather than reproducing the observed digits.

No threshold may change after the P200G run.

## Outcome

All original tests plus the repaired controls pass:

```text
HOMOGENEOUS_WEAK_POLE_LINE_REPLICATED_AFTER_CONTROL_REPAIR
```

Otherwise retain the existing preregistered failure hierarchy with a `_P200G`
diagnostic suffix in the artifact.  A positive result is labelled **DERIVED
COMPUTATIONAL, adversarially replicated after a disclosed control repair**.  It
still says nothing about the omitted pole equation or physical interpretation.

Only the targeted verifier and a static registry check may run.  No full suite.

