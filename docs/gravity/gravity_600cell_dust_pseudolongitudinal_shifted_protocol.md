# Protocol: shifted temporal persistence of the pseudo-longitudinal defect

Date: 2026-08-19

This protocol is committed before the shifted pseudo-longitudinal residuals
are calculated or inspected.  The current-tick result and the desired
falsifier are disclosed.  This is a target-disclosed temporal control, not a
blind discovery calculation.

## Frozen question and scope

On the accepted shifted centered recurrence of the same fixed 600-cell
carrier, test whether the tangential vertex-displacement image remains
non-invariant under the generalized action stiffness in each of the two
rank-15 negative sectors.

The result may establish persistence across two consecutive accepted
backgrounds.  It may not be called a curvature limit, a refinement limit,
continuum recovery, propagation or a lifetime.

## Frozen provenance

| input | SHA-256 |
|---|---|
| prior-art gate | `740eefaee14ea3ff634f8cff237041cecd675c4ceaf3d5be6ccdb9a3778a57ef` |
| shifted centered JSON | `265bd863de2365f19f7679373155fdaa23fb0bb3e75c221cfd9d9ec5b6ac2a47` |
| shifted centered NPZ | `c000f4fcae67e6c0648046878c2bd1ffd0616c38510ccf788c67cf99832397b8` |
| shifted direct-rank artifact | `86b53f228d6cfa7326a677d881463f1b849e76bc6c9ac2b0e8aa6fd427042944` |
| shifted direct-rank verifier | `1b54cd25899037fc66c2b58e01ef3bac267c6ebf2c6917d2a05ac4ac0feed1c5` |
| current scale-invariant artifact | `aee9088c1b0bb1cd4ebec12707014ed187f5c6cbe7ad84cc6cd74db705b8d20c` |
| current scale-invariant verifier | `2525ef28bace9c9e2a21ca715c25b0524d11353d3a3f5d649649bebaffe87ba6` |
| rigidity/York verifier | `deba8d9f9bca4a5848134943ec77544e5487d44a59c44234f632b6f2aeb51382` |
| conformal verifier | `d77dc8853826d9aecc4395fc4aae405d0505bbd644ec3a3229f640b2e980bcb4` |
| boundary-tangent verifier | `c8662bb0835865aac6696fc3f474ed668fed3fe393b9c32a59e709a984c35571` |
| binary-orbit verifier | `ad93cdd08fabeeee56b009f23936696837c4362f88ae23f92a36d0395e61ffaf` |
| 600-cell implementation | `ea5bce4b6c52e0834539ca4b1df9c6a67a3a5ed4da32f4e0298a493fc5315c7f` |
| consolidated current result | `1196c7155eb1057c8c6159e63ade56846e70395d8321702c4ef839b480bbbd83` |
| shifted direct result | `904732e25f2dc49c40557cf3da4daa4da172894f2990433522599435d63fbba0` |

The verifier must reject any provenance mismatch.  It must also require
`SHIFTED_CENTERED_CERTIFIED`, the shifted direct census
`15 negative + 10 positive` in all 16 cells and the preserved current
scale-invariant outcome `SCALE_INVARIANT_DIRECT_REFUTATION_CONFIRMED`.

## Mechanically frozen construction

1. Import the literal binary-orbit verifier and retain all 43 controls.
2. Rebuild normalized 600-cell vertices, adjacency and both parity edge
   orders independently from the shifted archive.
3. Rebuild the full rigidity matrix `R`, its tangential image `D=R T` and
   cell-incidence carrier `C`.  Require ranks `470`, `354`, `470` for
   `R`, `D`, `[C D]`, hence intersection dimension four.
4. Reconstruct all seven binary symmetry sectors at 100 decimal digits and
   use only sectors 4 and 5.
5. For each parity, sector and derivative variant, project the shifted
   archived `M_midpoint,V_midpoint` into the same action-selected
   conformal/shape carrier.  No shifted residual may be loaded from another
   artifact.
6. Use column-pivoted QR, not an eigenvector alignment, to construct the
   15-dimensional tangential image `L` inside the 25-dimensional shape
   carrier.
7. Compute the full shifted census before loading the current-tick residual
   values for comparison.

The 16 cells are exactly

```text
2 parities x sectors {4,5} x
{operational_primary, operational_shadow,
 validation_primary, validation_shadow}.
```

## Observables and resolution rule

With

```text
B = -W* M W,
A = -W* V W,
```

compute

```text
rho_span = ||(1-P_BL) A L||_2 / ||A L||_2,
rho_comm = ||(1-P_L) B^-1 A L||_2 / ||B^-1 A L||_2.
```

Both denominators must exceed `1e-12`.  The exact norm inequalities connecting
the unnormalized span and inverse residuals must hold within the same frozen
roundoff floor used by the current scale-invariant audit.

For each parity/sector family, the error for each dimensionless observable is

```text
maximum variation from operational_primary across four variants
+ maximum conditioning-derived relative floor across those variants.
```

Classification is frozen as

```text
value <= 10 error     ZERO_CONSISTENT
value > 100 error     NONZERO_RESOLVED
otherwise             OPEN.
```

The relative augmented rank of `[B L, A L]` uses threshold
`100 * max(observable errors) * largest singular value`.  Only the robust
comparison with 15 is interpreted.

## Post-census comparison

Only after all shifted cells have been classified may the verifier load the
current scale-invariant artifact.  For matching cells and variants it will
report, without fitting or a stability threshold,

```text
shifted rho_span / current rho_span,
shifted rho_comm / current rho_comm,
absolute shifted-current differences.
```

These ratios are descriptive.  Two time points cannot establish a trend.

## Frozen outcome hierarchy

1. `SHIFTED_PSEUDOLONGITUDINAL_CONTROL_FAILED` if any provenance, geometry,
   carrier, denominator, conditioning or upstream-outcome control fails.
2. `SHIFTED_PSEUDOLONGITUDINAL_DEFECT_PERSISTS` if all 16 shifted cells have
   both observables `NONZERO_RESOLVED` and augmented rank strictly above 15.
3. `SHIFTED_PSEUDOLONGITUDINAL_DEFECT_CLOSES` if all 16 shifted cells have
   both observables `ZERO_CONSISTENT` and augmented rank exactly 15.
4. `SHIFTED_PSEUDOLONGITUDINAL_OPEN` otherwise.

## Interpretation firewall

- `PERSISTS` is **DERIVED COMPUTATIONAL / STRUCTURAL** only for these two
  finite backgrounds and this derivative family.
- `CLOSES` would refute temporal persistence but would not prove continuum
  gauge symmetry.
- Mixed cells are `OPEN`, not majority-voted.
- The computation cannot decide whether the subsystem is physical,
  unstable, a pseudo-constraint or a coarse-lattice artefact.
- A material result requires a subsequent mechanically independent audit.
- No full-suite run is authorized; only this verifier and static registry
  guards may be run.

