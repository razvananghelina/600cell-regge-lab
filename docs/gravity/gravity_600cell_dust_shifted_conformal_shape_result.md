# Result: the conformal/shape split survives the shifted slab

Date: 2026-08-18

## Headline

**DERIVED COMPUTATIONAL.**  On the next independently accepted nonstationary
dust-Regge middle slice, the action-selected conformal and shape carriers are
again reducing for both centered recurrence coefficients.  All

```text
2 schedules * 7 sectors * 4 variants
* 2 operators * 2 carriers = 224
```

required residuals are `ZERO_CONSISTENT`, and all `112/112` even/odd schedule
comparisons are `SCHEDULE_ROBUST`.  The certified outcome is

```text
SHIFTED_CONFORMAL_SHAPE_DYNAMICS_DECOUPLED.
```

This removes the first failure branch of the persistence mission.  It does
not yet say that the old `30` negative-stiffness directions survive.

## Provenance

| stage | commit |
|---|---|
| prior-art and framing gate | `33da8dd` |
| target-disclosed protocol | `62810b4` |
| verifier registered before first execution | `6b9eef8` |
| byte-identical certified artifact | `a850fd9` |

Verifier:

```text
reproducible/verify_gravity_600cell_dust_shifted_conformal_shape.py
SHA-256 95b5c282c18cbdcf1fb9c87e5c8e62605062063f8c553f8b6f81f90c178c2303
```

Artifact:

```text
reproducible/gravity_600cell_dust_shifted_conformal_shape.json
SHA-256 acf9029114fee28800ad4cf6bff131fe59c5a8a22f7c8f7fa334b5a454842ec2
```

Two targeted executions returned byte-identical artifacts and `11/11 PASS`.
The full suite was not run.

## Complete hypothesis boundary

The result assumes the fixed labelled 600-cell edge carrier, the first three
accepted fixed-mass homothetic slabs, the literal edge identification,
both staircase schedules, all seven minimal symmetry sectors, all four frozen
derivative variants, the canonical unsigned vertex--edge incidence and the
committed midpoint/radius enclosures.  No constraint quotient, dust
perturbation, physical clock or continuum refinement is included.

## What passed

- exact carrier geometry and all `24` group actions;
- conformal rank `120` and shape rank `600`;
- all `56` action-relative shape carrier reconstructions;
- all `224` conformal/shape invariance residuals;
- all `112` schedule comparisons.

The high-precision symmetry-basis residual remains below `1.55e-98`.

## Interpretation

- **DERIVED COMPUTATIONAL:** conformal/shape closure is present at two
  consecutive centered recurrences, not only at the first one.
- **STRUCTURAL:** this is evidence for a persistent action-relative split on
  the fixed finite carrier.
- **OPEN:** shifted stiffness inertia, persistence or rotation of the two old
  negative subspaces, a canonical reduced two-step product, continuum wave
  physics and a limiting speed.
