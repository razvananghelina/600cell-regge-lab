# Result: shifted negative-shape persistence is precision-open

Date: 2026-08-18

## Headline

The action-selected conformal/shape decomposition survives the next slab, but
the old `30`-position negative-stiffness subsystem is neither certified nor
refuted there.

In every one of the `16` disclosed sector/schedule/variant cells, the shifted
Hermitian midpoint still has exactly

```text
15 negative eigenvalues + 10 positive eigenvalues,
```

but the `15` negative signs no longer cross the frozen `100 x error` boundary.
Their certified labels are

```text
15 OPEN + 10 POSITIVE_RESOLVED.
```

The final preregistered rank-gate outcome is therefore

```text
SHIFTED_NEGATIVE_RANK_OR_SECTOR_OPEN,
```

not `CHANGED` and not `COMMON_CARRIER_CERTIFIED`.

## Provenance ledger

| stage | commit |
|---|---|
| primary-literature and framing gate | `33da8dd` |
| persistence protocol | `62810b4` |
| shifted split verifier before execution | `6b9eef8` |
| shifted split artifact | `a850fd9` |
| blind shifted stiffness verifier before execution | `f7bb477` |
| disclosed zero-consistent semantics repair | `926e47d` |
| blind shifted stiffness artifact, before target comparison | `5b474c2` |
| disclosed open-vs-changed outcome repair | `c9ccf0f` |
| target-disclosed rank verifier before execution | `163e0a6` |
| deterministic rank-gate artifact | `730fa87` |

The two semantic repairs are recorded in the protocol.  Neither changed a
midpoint, error ball, carrier or sign count.

## Registered verifiers and artifacts

| object | SHA-256 |
|---|---|
| `verify_gravity_600cell_dust_shifted_conformal_shape.py` | `95b5c282c18cbdcf1fb9c87e5c8e62605062063f8c553f8b6f81f90c178c2303` |
| `gravity_600cell_dust_shifted_conformal_shape.json` | `acf9029114fee28800ad4cf6bff131fe59c5a8a22f7c8f7fa334b5a454842ec2` |
| `verify_gravity_600cell_dust_shifted_shape_stiffness.py` | `031d0dd1cab45d0093015fcab7ce7b56e098a5742895eed71f5a531aee31c2a6` |
| `gravity_600cell_dust_shifted_shape_stiffness.json` | `14fe5bc91e3ae4712c6ea19b8120785e2facd364e1ceb194009123fa353a4315` |
| `verify_gravity_600cell_dust_shifted_negative_persistence.py` | `ca9c347c42dab8b207197783ba95c1d1dfba68f2861286b9509003ff9980130e` |
| `gravity_600cell_dust_shifted_negative_persistence.json` | `3f241a21dc9c97b187d1a1b7e30dd9580f5b4605873d2efe705216d84f971b71` |

Each scientific verifier was executed twice with a byte-identical artifact.
Only targeted verifiers were run; the full suite was not run.

## Stage A: the split really persists

The shifted conformal/shape verifier reports:

```text
224/224 dynamic residuals      ZERO_CONSISTENT,
112/112 schedule comparisons  SCHEDULE_ROBUST,
11/11 checks                   PASS.
```

Thus the shape carrier remains a canonical reducing space for both shifted
`Gamma` and `Omega`.  The failure below is not leakage into the conformal
sector.

## Stage B: blind shifted inertia

The complete shifted census reports:

| full-multiplicity Hermitian label | count |
|---|---:|
| positive resolved | 2,000 |
| negative resolved | 0 |
| zero-consistent | 2,336 |
| open | 464 |

All `56/56` kinetic forms remain positive-definite-resolved, all normalized
eigenvalues remain real-consistent, all `56/56` action-compatibility
residuals are zero-consistent and all `56/56` schedule comparisons are
robust.  The blind outcome is

```text
SHIFTED_SHAPE_STIFFNESS_SIGN_OPEN.
```

This is an uncertainty result, not evidence of positive stiffness.

## Stage C: what happened to the old two sectors

The old cells have the exact resolved inertia `15 negative + 10 positive`.
The shifted cells have `15 open + 10 positive`, so their certified negative
rank interval is

```text
[0,15]
```

in all `16/16` comparisons.  Rank `15` remains allowed everywhere, but is
proved nowhere at the shifted tick.

The post-blind midpoint diagnostic is striking but not decisive:

```text
shifted midpoint negative count, all cells       15
maximum old/new ordered-spectrum distance        9.5424e-7
shifted/old complete error ratio                  4.55632649
shifted negative/positive cluster gap             1.75995e-2
```

**POST-BLIND STRUCTURAL PATTERN:** the central spectra barely move and the
`15`-dimensional low cluster remains sharply separated from the upper `10`.
But choosing midpoint signs would violate the error protocol, so this pattern
did not select the outcome.

## Why no projector or two-step product appears

The protocol required a completely resolved shifted `15+10` sign split before
calling the shifted fiber a negative-stiffness carrier.  Since that gate is
open, the verifier deliberately performed

```text
projector comparisons       0
reduced products            0.
```

This prevents a visually plausible midpoint cluster from being promoted into
a canonical propagating subsystem.  It also means that rotation versus common
carrier remains **OPEN**.

## Status ledger

- **DERIVED COMPUTATIONAL:** the conformal/shape split persists across two
  consecutive centered recurrences.
- **DERIVED COMPUTATIONAL:** the shifted kinetic form is definite and the
  shifted stiffness census is complete.
- **OPEN:** the old rank-15 negative carrier in each of sectors `4` and `5` is
  not resolved at the shifted tick.
- **NOT REFUTED:** rank `15` lies inside every shifted certified rank interval.
- **POST-BLIND STRUCTURAL PATTERN:** the two midpoint clusters and spectra are
  nearly unchanged.
- **OPEN:** common carrier, bundle rotation, nonautonomous reduced product,
  physical instability, graviton interpretation, continuum limit and `c`.

## Next load-bearing calculation

Audit the shifted restricted-form error budget without looking at a desired
sign: separate source-radius, carrier-lift and arithmetic contributions, then
preregister a higher-precision reconstruction only if an already available
action derivative can reduce the dominant term.  An arbitrary smaller error
bar or use of midpoint signs is forbidden.
