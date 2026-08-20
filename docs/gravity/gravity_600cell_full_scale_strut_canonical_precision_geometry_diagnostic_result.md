# Geometry-control diagnostic result

Date: 2026-08-20  
Status: **DERIVED diagnostic; resolver conclusion remains OPEN**

## Frozen execution

This execution used the diagnostic preregistered in commit `e6701d8` and the
instrumented verifier committed in `3dbe10c`. It changed no matrix, precision,
threshold, candidate, classification, or outcome rule.

- verifier SHA-256:
  `ed0a56bb8bf901b22cb4373bb3cc5e7446e2538fa0b4a1abd0ba26dc9ff32efc`;
- diagnostic artifact SHA-256:
  `97c7d46c7851b0bce6f8eef82ee196389f4ca4f45a22e2f19d0a822cd643e42a`;
- execution: `12/15`, outcome
  `FULL_SCALE_STRUT_CANONICAL_PRECISION_CONTROL_FAILED`.

## Exact localization

Every serialized diagnostic is internally consistent with the original aggregate
Boolean. At P100, all seven geometry conjuncts pass in both parities. At P160,
exactly one conjunct fails in both parities:
`maximum_imaginary_below_floor`.

| level/parity | maximum imaginary residue | frozen threshold | other conjuncts |
|---|---:|---:|---|
| P100/even | `3.8006165646679586e-79` | `1e-70` | all pass |
| P100/odd | `4.6280636534726750e-79` | `1e-70` | all pass |
| P160/even | `3.5680070678277676e-119` | `1e-130` | all pass |
| P160/odd | `7.4154130875215017e-119` | `1e-130` | all pass |

The other conjuncts are the complete 120-column coverage, 4440-entry carrier,
irrep dimensions `[1,1,1,2,2,2,3]`, branch-entry validation, and the two exact
negative-direction histograms. All pass at both levels and parities.

## Interpretation

- **DERIVED:** the failure is not a generic geometry-control failure. It is
  localized to the imaginary-residue threshold.
- **STRUCTURAL:** the observed scale is consistent with roundoff amplified by
  the smallest central-difference step: approximately `10^-100 / 10^-20 =
  10^-80` at P100 and `10^-160 / 10^-40 = 10^-120` at P160. This scale argument
  does not by itself certify a replacement threshold.
- **OPEN:** the primary intersection classification remains unaccepted until a
  separately preregistered, precision-aware arithmetic control passes.

No physical claim follows from this diagnostic.
