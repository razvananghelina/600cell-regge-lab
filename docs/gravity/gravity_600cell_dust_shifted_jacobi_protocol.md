# Preregistered blind protocol: shifted three-slice Jacobi operator

Date: 2026-08-18

Prior-art and framing gate: `920ce5d`.

The third-slab tangent artifact was committed in `d666315` before this
protocol.  No shifted boundary-twist determinant, principal-function block,
Jacobi coefficient, spectrum or comparison with an old mode carrier has been
inspected.

## Frozen inputs

| input | SHA-256 |
|---|---|
| `gravity_600cell_dust_two_step_full_tangent.json` | `f7fbf18535cc00dacec9a9ffa95f97f2d1847ac83073f27d39fcdb7968b0bafc` |
| `gravity_600cell_dust_two_step_full_tangent.npz` | `ce78ebf415584b1cdcf1d2cb07687135b624ad4939e0a4e54650653f7b384e6d` |
| `verify_gravity_600cell_dust_two_step_full_tangent.py` | `c1a3fb09146188c1932ab81629ab69817f2a2f19108fdf8d9e89d78b6de8f717` |
| `gravity_600cell_dust_later_slab_tangent.json` | `58a95d90d569b25a3aa396346f5198472c8aed706846b4057182b04c9f7480c4` |
| `gravity_600cell_dust_later_slab_tangent.npz` | `77b4dd54a5dcba9d1aa12870b361c9d7d7dde11ccaaa558361b9c1dc24768196` |
| `verify_gravity_600cell_dust_later_slab_tangent.py` | `3740622494b88bafc55839e7a53d9c4a335d4783503d236b0f2e9637429d9eea` |

Require

```text
TWO_STEP_FULL_TANGENT_COCYCLE_CERTIFIED  16/16,
LATER_SLAB_TANGENT_CERTIFIED             16/16,
```

their `448`-array counts and their recorded NPZ hashes.  Only the stored blind
tangent balls may be used.  In particular, files containing `centered`,
`conformal`, `shape`, `negative`, `root_count` or `kinetic_metric` results are
forbidden inputs.

Arithmetic uses 80-decimal Flint complex balls.  Re-enclose every binary64
midpoint with its stored radius plus one half-ULP independently in the real
and imaginary components.

## Carrier and mandatory variants

For both schedule parities and all seven minimal sectors use

```text
irrep dimensions      3,2,2,2,1,1,1,
position dimensions  90,60,60,60,30,30,30,
full position count  sum 30 d^2 = 720.
```

All four frozen derivative variants are mandatory.

Partition the committed balls as

```text
T_2 = [A_2 B_2; C_2 D_2],
T_3 = [A_3 B_3; C_3 D_3].
```

## Boundary-twist and principal-function gates

For every schedule, sector, variant and slab compute the Flint determinant of
`B_i`.  All

```text
2 schedules * 7 sectors * 4 variants * 2 slabs = 112
```

determinant balls must exclude zero.

For each regular tangent reconstruct

```text
S_i,01 = -B_i^-1,
S_i,00 =  B_i^-1 A_i,
S_i,10 =  C_i - D_i B_i^-1 A_i,
S_i,11 =  D_i B_i^-1.
```

Require entrywise Flint containment of zero for all adjoint identities and
for reconstruction of `A_i,C_i,D_i`.  Record full midpoint and radius norms;
a Boolean alone is insufficient.

## Shifted Jacobi operator

Construct without termwise normalization

```text
K_-^(2) = S_2,10,
K_0^(2) = S_2,11 + S_3,00,
K_+^(2) = S_3,01,

P_2 = -(K_+^(2))^-1 K_0^(2),
Q_2 = -(K_+^(2))^-1 K_-^(2).
```

The natural equation is

```text
K_-^(2) delta q_1 + K_0^(2) delta q_2
                         + K_+^(2) delta q_3 = 0.
```

No rescaling or comparison with an earlier coefficient is admitted.

## Independent product equivalence

Use the separately committed rigorous product `C_32=T_3T_2`.  In every cell
require the implicit and solved identities

```text
K_- + K_0 A_2 + K_+ A_32 = 0,
      K_0 B_2 + K_+ B_32 = 0,

P_2 A_2 + Q_2 = A_32,
P_2 B_2       = B_32
```

to contain zero entrywise.  These are the mandatory controls on momentum
sign, phase ordering and slice identification.

## Target-free census and archive

Record only:

- singular extrema and condition numbers of `B_2,B_3`;
- Frobenius norms of `K_-^(2),K_0^(2),K_+^(2),P_2,Q_2`;
- singular spectra of `[K_- K_0 K_+]` and `[P_2 Q_2]`;
- the background-asymmetry diagnostic
  `||K_+ - K_-^*||/max(||K_+||,||K_-||)`;
- all complete ball residual norms and radius envelopes.

For the fourteen schedule comparisons use the inherited singular-spectrum
uncertainty and labels `SCHEDULE_ROBUST`, `SCHEDULE_OPEN` and
`SCHEDULE_DEPENDENT` with the same `10/100` bands.

The deterministic archive stores exactly

```text
2 schedules * 7 sectors * 4 variants
* 5 matrices * 2 fields = 560 arrays.
```

## Frozen outcome hierarchy

1. `SHIFTED_JACOBI_CONTROL_FAILED` for any provenance, carrier or archive
   failure.
2. `SHIFTED_JACOBI_TWIST_SINGULAR` if any determinant ball contains zero.
3. `SHIFTED_JACOBI_VARIATIONAL_IDENTITY_FAILED` if any principal-function or
   product-equivalence identity fails.
4. `SHIFTED_JACOBI_SCHEDULE_DEPENDENT` if a primary comparison is resolved
   dependent.
5. `SHIFTED_JACOBI_SCHEDULE_OPEN` if none is dependent but at least one is
   open.
6. `SHIFTED_JACOBI_CERTIFIED` only if every structural control passes and all
   fourteen comparisons are robust.

## Explicit exclusions

- no centered `M,N,V` spectrum at this stage;
- no old/new subspace comparison;
- no conformal, shape or negative-mode carrier;
- no physical mode, wave, speed or growth interpretation;
- no full-suite run.

