# Preregistered blind protocol: three-slice Jacobi operator

Date: 2026-08-18

Prior-art gate commit: `aca7971`.

No boundary twist determinant, reconstructed principal-function block,
three-slice coefficient, recurrence spectrum or spatial comparison has been
inspected before this protocol commit.

## Frozen inputs

| input | SHA-256 |
|---|---|
| `gravity_600cell_dust_full_boundary_tangent.json` | `4da8bcd2890a54bc9d3b60c6195df2933ea56194d942ab0285b51599ba287bd5` |
| `gravity_600cell_dust_full_boundary_tangent.npz` | `816c605da2a655442bbadce7a23965f0822f99e7bdc1d0a4a27af548de85446b` |
| `verify_gravity_600cell_dust_full_boundary_tangent.py` | `c8662bb0835865aac6696fc3f474ed668fed3fe393b9c32a59e709a984c35571` |
| `gravity_600cell_dust_two_step_full_tangent.json` | `f7fbf18535cc00dacec9a9ffa95f97f2d1847ac83073f27d39fcdb7968b0bafc` |
| `gravity_600cell_dust_two_step_full_tangent.npz` | `ce78ebf415584b1cdcf1d2cb07687135b624ad4939e0a4e54650653f7b384e6d` |
| `verify_gravity_600cell_dust_two_step_full_tangent.py` | `c1a3fb09146188c1932ab81629ab69817f2a2f19108fdf8d9e89d78b6de8f717` |

Require the upstream outcomes

```text
FULL_BOUNDARY_TANGENT_BLIND_CENSUS_CERTIFIED,
TWO_STEP_FULL_TANGENT_COCYCLE_CERTIFIED,
```

the array counts `224` and `448`, and their recorded archive hashes.  Only
stored target-free tangent balls may be used; no spectrum from another
operator may enter.

Arithmetic uses 80-decimal Flint complex balls.  Every binary64 midpoint is
re-enclosed with its stored Flint radius plus one half-ULP independently in
its real and imaginary components.

## Carrier and mandatory variants

For each of the two schedule parities and seven minimal sectors use the
position dimensions

```text
n = 90, 60, 60, 60, 30, 30, 30
```

with representation weights

```text
d = 3, 2, 2, 2, 1, 1, 1.
```

All four derivative variants are mandatory:

```text
operational_primary, operational_shadow,
validation_primary,  validation_shadow.
```

Partition each committed tangent ball as

```text
T_i = [ A_i  B_i ]
      [ C_i  D_i ].
```

## 1. Boundary twist gate

For every schedule, sector, variant and slab, compute the Flint determinant
ball of `B_i`.  There are

```text
2 schedules * 7 sectors * 4 variants * 2 slabs = 112
```

twist determinants.  Every determinant must exclude zero before an inverse
or recurrence is accepted.

Also report binary midpoint singular extrema and condition numbers, with
uncertainty from all variants, stored ball radii and the serialization floor.
The determinant ball, not a binary SVD threshold, decides regularity.

Full pre-Legendre regularity is not allowed to substitute for this gate.

## 2. Reconstruct the quadratic principal functions

For each regular tangent block construct

```text
S_i,01 = -B_i^-1,
S_i,00 =  B_i^-1 A_i,
S_i,10 =  C_i - D_i B_i^-1 A_i,
S_i,11 =  D_i B_i^-1.
```

The following Flint-ball residuals must contain zero entrywise:

```text
S_i,00 - S_i,00*,
S_i,11 - S_i,11*,
S_i,10 - S_i,01*,
A_i - B_i S_i,00,
C_i - (S_i,10 + S_i,11 A_i),
D_i - S_i,11 B_i.
```

The star is conjugate transpose in a complex minimal sector.  Record
midpoint Frobenius norms and radius envelopes rather than printing only a
Boolean.

## 3. Three-slice Hessian and normalized recurrence

Construct, without rescaling individual terms,

```text
K_- = S_1,10,
K_0 = S_1,11 + S_2,00,
K_+ = S_2,01,

P = -K_+^-1 K_0,
Q = -K_+^-1 K_-.
```

The natural unnormalized equation is

```text
K_- delta q_0 + K_0 delta q_1 + K_+ delta q_2 = 0.
```

The `P,Q` form is only its solved coordinate recurrence.  No matrix norm or
eigenvalue of `P,Q` is a physical frequency at this gate.

## 4. Independent product equivalence

The two-step archive contains a separately committed rigorous product ball

```text
C_21 = [ A_21  B_21 ]
       [ C_21  D_21 ].
```

For every schedule, sector and variant require the two linearized seam
identities to contain zero entrywise:

```text
K_- + K_0 A_1 + K_+ A_21 = 0,
      K_0 B_1 + K_+ B_21 = 0.
```

Also require the equivalent solved identities

```text
P A_1 + Q = A_21,
P B_1     = B_21.
```

These controls falsify an incorrect momentum sign, phase ordering, boundary
identification or block convention.  A good spectrum cannot override a
failed variational identity.

## 5. Target-free census and schedule comparison

For each schedule and sector, report for the operational primary midpoint:

- singular extrema and condition numbers of `B_1`, `B_2`;
- Frobenius norms of `K_-`, `K_0`, `K_+`, `P`, `Q`;
- singular spectra of the two horizontal operators
  `[K_- K_0 K_+]` and `[P Q]`;
- the background-asymmetry diagnostic
  `||K_+ - K_-*||_F / max(||K_+||_F,||K_-||_F)`;
- all ball residual norms and their radius envelopes.

The asymmetry is a blind diagnostic.  No small/large target and no physical
label is attached to it.

For each of the fourteen schedule comparisons (`7` sectors times the two
horizontal operators), compare the ordered singular spectra.  The
uncertainty is the sum of:

1. the maximum operator-norm displacement among all four derivative
   variants;
2. the Frobenius norm of the complete Flint radius matrix;
3. `10 eps_machine` times the largest compared singular value.

Use the frozen labels

```text
SCHEDULE_ROBUST       distance <= 10 epsilon,
SCHEDULE_DEPENDENT    distance > 100 epsilon,
SCHEDULE_OPEN         otherwise.
```

The output NPZ stores midpoint and radius arrays for exactly five matrices
`K_-,K_0,K_+,P,Q`:

```text
2 schedules * 7 sectors * 4 variants * 5 matrices * 2 fields = 560 arrays.
```

## Frozen outcome hierarchy

1. `THREE_SLICE_JACOBI_CONTROL_FAILED` for any provenance, carrier, archive
   or upstream-outcome failure.
2. `THREE_SLICE_JACOBI_TWIST_SINGULAR` if any of the 112 determinant balls
   contains zero.
3. `THREE_SLICE_JACOBI_VARIATIONAL_IDENTITY_FAILED` if a regular construction
   violates any adjoint, tangent-recovery or product-equivalence ball
   identity.
4. `THREE_SLICE_JACOBI_SCHEDULE_DEPENDENT` if any primary schedule comparison
   is resolved dependent.
5. `THREE_SLICE_JACOBI_SCHEDULE_OPEN` if none is dependent but at least one is
   open.
6. `THREE_SLICE_JACOBI_CERTIFIED` only if all controls, twists and
   variational identities pass and all fourteen schedule comparisons are
   robust.

The outcome is independent of whether any coefficient happens to resemble a
finite difference or a spatial Laplacian.

## Explicit exclusions

- no graph, Hodge, Kähler--Dirac or continuum `S^3` spectrum;
- no scalar/vector/tensor labels or desired degeneracies;
- no fitted temporal/spatial factorization;
- no frequency, dispersion, limiting speed or Planck unit;
- no refinement, third slab or nonlinear anisotropic solve;
- no full-suite run.

