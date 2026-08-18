# Third-tick stationary roots versus canonical target

Date: 2026-08-16

## Provenance and process disclosure

- target-independent roots committed first: `3401137`;
- comparison specification with disclosure: `3f665fc`;
- registered verifier: `9bad4b0`;
- targeted verifier:
  `reproducible/verify_gravity_600cell_dust_third_tick_stationary_target.py`;
- artifact:
  `reproducible/gravity_600cell_dust_third_tick_stationary_target.json`;
- artifact SHA-256:
  `4d1f81dafcab9d3aa40ff08fdaaad90b80235809dd32becd790bdee1704ab6cf`.

Only this targeted verifier was run.  It returned **5/5**.  The full suite was
not run.

**PROCESS DISCLOSURE:** after the clean target-independent root commit, a
read-only shell diagnostic printed the scalar target residuals before the
formal comparison specification was committed.  Therefore the comparison
arithmetic is reproducible but not cleanly preregistered.  The root-enumeration
firewall remains intact and no root or bound was changed.

## Mechanical verdict

```text
STATIONARY_THIRD_TICK_NO_HIT
hits = 0/2
```

| root | structural label | residual/component | norm | bound | hit |
|---|---|---:|---:|---:|---:|
| 0 | contracting | `+1.6161137465e-9` | `8.8518195448e-9` | `3.6513653962e-21` | no |
| 1 | time reversal | `-9.0808989801e-3` | `4.9738132138e-2` | `3.6513653962e-21` | no |

Both roots were compared on all 30 components and both schedule parities.
The contracting mismatch is uniform and small in scale but exceeds the
certified norm bound by about `2.42e12`; it is not an approximate pass.

## Interpretation

- **DERIVED COMPUTATIONAL:** holding the inherited second-tick lapse exactly
  fixed does not give a third canonical seam.
- **DERIVED:** the time-reversal root has the wrong momentum sign and is not a
  forward solution.
- **STRUCTURAL:** the contracting root was selected geometrically before the
  target and has a resolved nonzero `(G,P)` determinant.
- **OPEN:** whether a cleanly preregistered local `(C,R)` correction produces
  the third tick.

The sole legitimate next seed is committed root 0.  The already committed
Jacobian and mismatch predict, before any new action evaluation,

```text
delta_C approximately +9.9819435e-11,
delta_R approximately -1.7796186e-5,
u3/u1 approximately 3.0000107,
v3/v1 approximately 4.9999775.
```

These are diagnostics for a later clean protocol, not accepted results.
