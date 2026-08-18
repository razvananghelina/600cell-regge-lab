# Result: the rank-30 negative-shape subsystem persists

Date: 2026-08-18

## Headline

After removing the selected binary64 tangent-serialization interface and
applying the disclosed canonical Hessian-symmetry control, the two old
rank-`15` negative-stiffness sectors persist at the shifted tick.

Every one of the `16 = 2 parities * 2 sectors * 4 derivative schedules`
cells is completely resolved as

```text
15 NEGATIVE_RESOLVED + 10 POSITIVE_RESOLVED.
```

The preregistered outcome is

```text
SHIFTED_DIRECT_NEGATIVE_RANK_PERSISTS.
```

This certifies a persistent rank-`30` negative subsystem on the finite
action-selected shape carrier.  It does not yet certify that the two
rank-`15` fibers are the same subspaces or are connected by an
action-selected transport.

## Complete provenance ledger

| stage | commit |
|---|---|
| primary-literature and precision gate | `52f337b` |
| serialization-attribution artifact | `67e99fb` |
| original direct-precision protocol | `4db96b0` |
| registered verifier before first execution | `4e731d9` |
| first `CONTROL_FAILED` artifact | `cab30d3` |
| target-independent identity diagnostics | `d2069e0` |
| disclosed protocol amendment | `2ba64c7` |
| canonical Hessian-symmetry implementation | `fb49173` |
| final deterministic artifact | `99d2766` |

The final verifier source has SHA-256

```text
1b54cd25899037fc66c2b58e01ef3bac267c6ebf2c6917d2a05ac4ac0feed1c5
```

and its JSON artifact has SHA-256

```text
86b53f228d6cfa7326a677d881463f1b849e76bc6c9ac2b0e8aa6fd427042944.
```

The corrected verifier was executed twice with a byte-identical artifact and
reported `14/14` checks passed both times.  Only targeted verifiers were run;
the full suite was not run.

## Why the first run was rejected

The first direct run already produced `15+10` in all `16` cells, but failed a
preregistered control, so its signs were not accepted.  All recovery
identities passed; only the three principal-function adjoint identities
failed.  Their residuals scaled exactly as the square of the complex-step
size:

```text
h near 1e-20: residual about 1e-33 ... 1e-31,
h near 1e-15: residual about 1e-23 ... 1e-21.
```

The diagnosis exposed a real methodological distinction: Flint propagated
balls rigorously after receiving the finite-step Hessian, but the numerical
finite-step Hessian was not itself an interval enclosure of the exact
derivative.  The former broad binary balls had hidden its tiny antisymmetric
truncation defect.

The repair was fixed before the corrected execution:

1. record the raw antisymmetric defect and the complete four-schedule family
   variation;
2. reject any defect outside that family variation;
3. use the unique Hermitian projection `(H+H*)/2`, selected by the exact
   symmetry of an action Hessian;
4. add the ordered stiffness-eigenvalue variation across the four schedules
   to the final sign error.

No coefficient was fitted and no derivative step, precision or sign threshold
was changed.

## Complete numerical controls

The corrected execution reports:

```text
4/4 direct slab branch reconstructions                 PASS
32/32 raw Hessian defects inside schedule variation   PASS
32/32 boundary-twist determinants exclude zero        PASS
32/32 principal identity families                     PASS
16/16 M,V balls overlap the old broad controls        PASS
16/16 V-radius reduction factors > 100                PASS
16/16 action-selected shape carriers                  PASS
16/16 restricted kinetic forms positive definite     PASS
16/16 action-compatibility residuals zero-consistent  PASS
```

The final `V` radius is reduced relative to the archived route by a factor of
about

```text
1.973094e9.
```

The raw Hessian symmetry correction is at most about `2.51e-5` of the full
four-schedule variation.  At binary output resolution, the ordered stiffness
eigenvalues agree across schedules; the added observed eigenvalue variation
is zero.  The complete restricted error is nevertheless retained at about
`9.72e-8`, and the least-separated stiffness sign remains about `2014.7`
error units from zero, well beyond the frozen `100`-unit resolution boundary.

The aggregate minimal-sector census over the `16` cells is

```text
240 NEGATIVE_RESOLVED + 160 POSITIVE_RESOLVED.
```

Because each target sector has irrep multiplicity one here, this is the same
rank-`15` result per cell, not a multiplicity-weighting artefact.

## Status ledger

- **DERIVED COMPUTATIONAL, conditional on the frozen derivative family:**
  sectors `4,5` each retain negative stiffness rank `15` at the shifted tick,
  for both parities and all four schedules.
- **DERIVED COMPUTATIONAL:** the previous OPEN verdict was caused by the
  binary tangent-serialization interface, whose removal improves the final
  coefficient radius by about nine orders of magnitude.
- **STRUCTURAL:** the unique Hermitian projection is fixed by action-Hessian
  symmetry and is much smaller than the schedule-family variation.
- **OPEN:** a formal analytic/automatic ball derivative independent of the
  finite-step family.
- **OPEN:** equality or canonical transport of the old and shifted negative
  fibers, projector rotation, reduced propagator and longer-time persistence.
- **OPEN:** physical instability, graviton or wave interpretation, inertial
  mass, continuum/refinement behavior and an effective limiting speed.

Here “inertia” still means the signature of the quadratic action—the number
of positive and negative curvature directions—not particle inertia.

## What follows if this survives

The next load-bearing gate is the fiber, not another rank count.  Reconstruct
the old and shifted rank-`15` spectral projectors with certified gaps and ask
whether the action itself selects a connection between them.  An arbitrary
Procrustes, polar or basis alignment is forbidden.  If no action-selected
connection exists, the result remains a persistent count rather than a
propagating physical subsystem.
