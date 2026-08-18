# Preregistered blind protocol: centered Jacobi coefficient census

Date: 2026-08-18

Prior-art gate commit: `0b69556`.

No centered coefficient, determinant, inertia, generalized stiffness
eigenvalue or desired spatial spectrum has been inspected before this commit.

## Frozen inputs

| input | SHA-256 |
|---|---|
| `gravity_600cell_dust_three_slice_jacobi.json` | `514e01937d621e82c240ea5cad621fb2bc699d09c4940b9be46fa1498152d90c` |
| `gravity_600cell_dust_three_slice_jacobi.npz` | `63d95e79c11b25cada660f9a2422654eb92180263dad64e1cbf0ecc30b67d7f8` |
| `verify_gravity_600cell_dust_three_slice_jacobi.py` | `a875751eebb202dbb0b92780c3f48e7e275470360442f3ed1154e310cc36a884` |

Require outcome `THREE_SLICE_JACOBI_CERTIFIED`, `560` stored arrays, the
recorded archive SHA, `spatial_target_loaded=false`, and both schedules with
the seven dimensions `3,2,2,2,1,1,1`.

Re-enclose every binary midpoint using its stored Flint radius plus one
half-ULP in each real and imaginary component.  Arithmetic uses 80-decimal
Flint complex balls.  All four derivative variants are mandatory.

## 1. Unique centered decomposition

For every schedule, sector and variant construct

```text
M = (K_- + K_+)/2,
N = (K_+ - K_-)/2,
V =  K_- + K_0 + K_+.
```

The following residuals must contain zero entrywise:

```text
K_- - (M-N),
K_+ - (M+N),
K_0 - (V-2M).
```

Record midpoint Frobenius norms and complete radius envelopes.

## 2. Centered mass regularity

Compute all

```text
2 schedules * 7 sectors * 4 variants = 56
```

Flint determinant balls of `M`.  A determinant ball, not a binary singular
threshold, decides invertibility.  Only if every determinant excludes zero
may the normalized operators be constructed:

```text
Gamma = M^-1 N,
Omega = M^-1 V.
```

Require the normalized recurrence residuals to contain zero entrywise:

```text
M^-1 K_- - (I-Gamma),
M^-1 K_+ - (I+Gamma),
M^-1 K_0 - (-2I+Omega).
```

The output NPZ contains midpoint and radius arrays for exactly five matrices
`M,N,V,Gamma,Omega`:

```text
2 schedules * 7 sectors * 4 variants * 5 matrices * 2 fields = 560 arrays.
```

## 3. Hermitian-part inertia

For the operational-primary midpoint of each `M`, form

```text
H_M = (M+M*)/2.
```

Its error `epsilon_H` is the sum of:

1. maximum operator-norm midpoint variation over all four variants;
2. maximum Flint-radius Frobenius norm over those variants;
3. `10 eps_machine max(1,||H_M||_2)`.

Classify each eigenvalue of `H_M` as

```text
POSITIVE_RESOLVED   lambda >  100 epsilon_H,
NEGATIVE_RESOLVED   lambda < -100 epsilon_H,
ZERO_CONSISTENT    |lambda| <  10 epsilon_H,
OPEN               otherwise.
```

Restore full-carrier counts with representation weight `d`.  Label a sector
`POSITIVE_DEFINITE` or `NEGATIVE_DEFINITE` only when every eigenvalue is
resolved with that sign, `INDEFINITE` when both resolved signs occur, and
`INERTIA_OPEN` otherwise.

Also report the adjoint-defect ratios

```text
||X-X*||_F / max(1,||X||_F),  X=M,N,V.
```

These are **STRUCTURAL DIAGNOSTICS**.  An indefinite Hermitian part warns
against a naive no-ghost reading, but cannot be promoted to a coordinate-free
ghost theorem because the time-fibre identification is not a derived
superspace parallel transport.

## 4. Blind normalized-operator census

For `Gamma` and `Omega`, report:

- Frobenius norm, singular spectrum and condition number;
- eigenvalues, eigenvector condition number and variant/ball error;
- for every `Omega` eigenvalue, `REAL_CONSISTENT`, `RESOLVED_COMPLEX` or
  `COMPLEX_OPEN`, using respectively `<10`, `>100`, or the intermediate
  multiple of the Bauer--Fike eigenvalue error;
- full-carrier weighted reality counts.

The matrix error is defined by the same variant, radius and binary floor sum
as above.  The eigenvalue error is that matrix error times the operational
eigenvector condition number, plus the binary eigenvalue floor.  If the
eigenvector matrix is singular or nonfinite, every eigenvalue is open.

No target number of real modes is specified.  The outcome of the construction
does not depend on the reality count.

## 5. Schedule comparison

For each of the seven sectors compare the ordered singular spectra of
`Gamma` and `Omega`, giving fourteen primary comparisons.  The comparison
error is the sum of both schedules' matrix errors plus
`10 eps_machine` times the largest singular value.

Use

```text
SCHEDULE_ROBUST       distance <= 10 epsilon,
SCHEDULE_DEPENDENT    distance > 100 epsilon,
SCHEDULE_OPEN         otherwise.
```

Also report the optimally matched `Omega` eigenvalue distance with the sum of
both Bauer--Fike errors.  It is secondary because the matrices may be
nonnormal; it may veto only if it is resolved dependent under its own error.

## Frozen outcome hierarchy

1. `CENTERED_JACOBI_CONTROL_FAILED` for any provenance, carrier, source
   archive or output-cardinality failure.
2. `CENTERED_JACOBI_MASS_SINGULAR` if any of the 56 `M` determinant balls
   contains zero.
3. `CENTERED_JACOBI_IDENTITY_FAILED` if any algebraic or normalized
   recurrence ball identity fails.
4. `CENTERED_JACOBI_SCHEDULE_DEPENDENT` for a resolved primary comparison or
   a resolved secondary eigenvalue comparison.
5. `CENTERED_JACOBI_SCHEDULE_OPEN` if no comparison is dependent but a
   primary comparison is open.
6. `CENTERED_JACOBI_CERTIFIED` only if all controls, mass determinants,
   identities and all fourteen primary schedule comparisons pass robustly.

Inertia and generalized-stiffness reality remain reported scientific results
under every structurally certified outcome; they do not move the hierarchy.

## Explicit exclusions

- no Whitney/Hodge, graph, Kähler--Dirac, Regge `curl^T curl`, Lichnerowicz or
  continuum harmonic spectrum;
- no scalar/vector/tensor or graviton labels;
- no dropping of `N` and no fitted field redefinition;
- no frequency, dispersion, proper-time normalization, `c` or Planck scale;
- no refinement or full-suite run.

