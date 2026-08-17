# Frozen nonlinear boundary-covariance cases

Date: 2026-08-17

This is the Stage-A record committed before evaluating any perturbed nonlinear
canonical equation or comparing any nonlinear output.

## Provenance

- prior-art gate: `526a202`
- protocol: `05f76c3`
- first implementation, before its interrupted run: `fc45d7e`
- precision-initialization correction: `8981046`
- corrected implementation, before the successful run: `0048f06`

Artifact:
`reproducible/gravity_600cell_dust_nonlinear_boundary_covariance_seeds.json`

Artifact SHA-256:

```text
2104c69ba6b21d3a3d92c7071d7f2702cb7d33f7f0e3ff17954f64c469f0c01d
```

The corrected targeted run passed `10/10`.  The full suite was not run.

## Calibration

Both schedules reproduce

```text
s_min(J)       = 4.2445618107e-9
epsilon_J      = 6.8022e-23 even
                 6.8024e-23 odd
response error = 4.0286e-16.
```

All 1,522 base/derivative action evaluations remain on the Lorentzian branch,
all Hessian entries pass operational/validation calibration, both complete
Hessians pass reciprocity, and both canonical matrices have resolved rank 65.

The first attempt was interrupted before the odd Hessian completed because its
precision context was initialized in the wrong order.  It wrote no artifact,
derived no case and evaluated no nonlinear perturbation.  The correction and
observed failed calibration are recorded separately rather than erased.

## Frozen cases

The common momentum unit is

```text
p_star = 0.00090810444890653157621006281318...
```

Exactly 32 paired inputs are committed:

```text
4 older Helmert directions
x 2 sectors: position, momentum
x 2 signs
x 2 levels: half, full.
```

The full-level amplitudes derived from the fixed unknown-response norm
`ETA=1e-4` are:

| direction | position amplitude | dimensionless momentum amplitude |
|---:|---:|---:|
| 1 | `9.71701466696e-12` | `2.02979199813e-4` |
| 2 | `9.71701467062e-12` | `2.02979200121e-4` |
| 3 | `9.71701467427e-12` | `2.02979193803e-4` |
| 4 | `9.71701467792e-12` | `2.02979180725e-4` |

The artifact records

```text
nonlinear_perturbed_action_evaluations = 0
nonlinear_outputs_compared             = false.
```

## Status before Stage B

- **DERIVED:** all 32 inputs and both operational/validation response seeds are
  fixed by the calibrated linear system.
- **STRUCTURAL:** the large position/momentum amplitude difference follows from
  the chosen dimensionless momentum unit and the canonical response; it was not
  tuned against a nonlinear result.
- **OPEN:** nonlinear solvability and even/odd output covariance for every case.
