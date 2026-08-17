# Preregistered protocol: identify the homogeneous curvature-kernel line

Date: 2026-08-17

Prior-art gate commit: `9177531`.

This is a target-disclosed confirmatory protocol.  It was committed before a
right-singular kernel vector was compared with any tangent, uniform, weak or
geometric-lapse candidate.

## Frozen inputs

| input | SHA-256 |
|---|---|
| `gravity_600cell_dust_internal_curvature_response.json` | `95b6edd8e21ad20a0db97a7c8e7027db7da6547b2b994ad1eb595cf2307f29dc` |
| `verify_gravity_600cell_dust_internal_curvature_response.py` | `276982879fae5f8fa735f27a6fa30bfe965dc3e41c169d8a229a61c23511ae66` |
| `gravity_600cell_dust_full_boundary_tangent.json` | `4da8bcd2890a54bc9d3b60c6195df2933ea56194d942ab0285b51599ba287bd5` |
| `gravity_600cell_dust_full_boundary_tangent.npz` | `816c605da2a655442bbadce7a23965f0822f99e7bdc1d0a4a27af548de85446b` |
| `verify_gravity_600cell_dust_full_boundary_tangent.py` | `c8662bb0835865aac6696fc3f474ed668fed3fe393b9c32a59e709a984c35571` |
| `verify_gravity_600cell_dust_hyperbolic_lapse_alignment.py` | `e461296a965c9b80fb89fae5660ce642858f3d3dfa0b24ccdecc2aced53c7047` |
| `verify_gravity_600cell_dust_full_anisotropic_legendre_rank.py` | `834b97c85e386def853b6308e65e831c52d62d7cbcc4b23118602120d6c676e5` |
| `verify_gravity_global_regge_orbits.py` | `ad93cdd08fabeeee56b009f23936696837c4362f88ae23f92a36d0395e61ffaf` |
| `gravity_600cell_dust_homothetic_canonical_lapse.json` | `4b1c59c0518eec11b88b140cdecdf558d762c0d70b4826a758f67544e14ac5b9` |

Both schedule parities and the same four 100-decimal/80-decimal-ball
derivative variants are mandatory.  Only the trivial `2T` minimal sector is
needed; no nontrivial sector or full suite will run.

## Reconstruction and controls

For each schedule parity:

1. Reconstruct the full Hessian, the canonical response `Y`, the
   incidence-derived internal deficit Jacobian and

   ```text
   F : C^60 -> C^160
   ```

   exactly as in the passing curvature verifier.
2. Require the trivial sector by constant-vector overlap above `1/2`, all
   branch/reality/equivariance controls, all four Flint determinants excluding
   zero and reproduction of the committed `rank 59 / nullity 1` ledger.
3. Define `K` from the last right singular vector of `F`.  Require one
   singular value under the fixed zero gate, the other 59 over the fixed
   resolved-nonzero gate, and a smallest/next-smallest separation above
   `1e4` in every variant.
4. Compare the four kernel lines with one another.  The operational/shadow
   and operational/validation line distances enter every later uncertainty.

## Near-`-1` tangent selection

For each tangent variant, order all 60 eigenvalues by `abs(lambda + 1)` and
select the first two.  This is the disclosed post-result target.  Require:

- exactly one selected eigenvalue has modulus below one and one above one;
- their imaginary parts are inside the calibrated eigensolver error;
- their reciprocal-conjugate defect is inside the calibrated error;
- the ratio between the third and second `abs(lambda+1)` distances exceeds
  10 in every variant.

The individual direct eigenvector lines are named `near_minus_contracting`
and `near_minus_expanding`.  Their two-plane is reconstructed independently
as an ordered complex-Schur subspace selected by the geometric-mean boundary
between the second and third `abs(lambda+1)` distances.  Direct-versus-Schur
plane distance is a control, not a hit.

Tangent invariance of `K` is the line distance between `K` and `T K`.  If
identified, the Rayleigh multiplier and residual are reported.  An eigenline
candidate cannot receive the final `IDENTIFIED_EIGENLINE` signature unless
this invariance test is also identified in both schedules.

## Frozen candidate census

For the source kernel line `K`, compare separately with:

```text
1  near_minus_contracting line,
2  near_minus_expanding line,
3  their near-minus two-plane,
4  uniform position line (1_30,0_30),
5  uniform momentum line (0_30,1_30),
6  their uniform phase two-plane,
7  pure-position C^30 subspace,
8  pure-momentum C^30 subspace.
```

Transport `K` by `Y` and compare separately with:

```text
9  the five-dimensional canonical weak Schur lift,
10 the five-dimensional independently derived geometric lapse subspace.
```

This gives exactly 10 comparisons per parity and 20 total.  Also report the
one cross-schedule comparison between the two source kernel lines in the
literal common boundary-phase ordering.  No best-of-candidates combination
is permitted.

## Distance calibration and labels

All comparisons use the sine of the largest principal angle from the line to
the candidate subspace.  For each comparison form an absolute uncertainty
from all applicable terms:

1. the three frozen derivative-variant distance changes;
2. kernel-line perturbation bounded by the full-response perturbation divided
   by its singular gap;
3. tangent-line/plane perturbation from tangent-ball radius, eigensystem
   condition and selected spectral separation;
4. Flint response/lift radius propagated through `Y`;
5. direct-eigenvector versus Schur-plane discrepancy where applicable;
6. a dimension- and condition-scaled binary64 QR/SVD/eigensolver floor;
7. `1e-70`.

Assign:

```text
IDENTIFIED       distance <= 10 epsilon,
SEPARATED        distance > 100 epsilon,
NUMERICALLY_OPEN otherwise.
```

Because these distances lie in `[0,1]`, an accumulated angular uncertainty
of `1e-2` or larger is by definition insufficient for identification and is
forced to `NUMERICALLY_OPEN`.  This prevents a very loose perturbation bound
from turning every possible line into a false hit.

No raw angle, coordinate correlation or visual similarity overrides these
labels.

## Outcome hierarchy

1. `HOMOGENEOUS_CURVATURE_KERNEL_CONTROL_FAILED` for any provenance,
   reconstruction, nullity, selection or calibration failure.
2. `HOMOGENEOUS_CURVATURE_KERNEL_SCHEDULE_DEPENDENT` if the two kernel lines
   are resolved separated in their literal common ordering.
3. `HOMOGENEOUS_CURVATURE_KERNEL_EIGENLINE_IDENTIFIED` if the same individual
   near-`-1` branch is identified in both schedules and `T K = K` is
   identified in both.
4. `HOMOGENEOUS_CURVATURE_KERNEL_LINE_IDENTIFIED` if the same exact uniform
   position or momentum line is identified in both schedules.
5. `HOMOGENEOUS_CURVATURE_KERNEL_SUBSPACE_LOCALIZED` if no line is identified
   but at least one of the preregistered higher-dimensional candidates is
   identified consistently in both schedules.
6. `HOMOGENEOUS_CURVATURE_KERNEL_UNIDENTIFIED_OR_OPEN` otherwise.

Membership in the geometric lapse or canonical weak subspace is reported
independently even when an earlier outcome applies.  It does not by itself
justify the word “gauge.”

## Explicit exclusions

- no fitted basis change between schedules;
- no candidate invented from the observed kernel coordinates;
- no nonlinear continuation in this mission;
- no claim that a deficit-kernel line solves all canonical constraints;
- no continuum dispersion, speed or Planck-scale comparison;
- no full-suite run.
