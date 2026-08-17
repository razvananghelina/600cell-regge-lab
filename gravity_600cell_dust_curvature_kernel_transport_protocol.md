# Preregistered protocol: transport of the curvature-kernel line

Date: 2026-08-17

Prior-art gate commit: `ed51853`.

The target equation `T_1 K_1 = K_2` is disclosed before constructing the
second-slab curvature operator.  No second-slab kernel coordinate or transport
distance has been inspected before this commit.

## Frozen inputs

| input | SHA-256 |
|---|---|
| `gravity_600cell_dust_homogeneous_two_by_two.json` | `d0017d4cfdf3a8833cf19bfcd287b21ac91a7f631c803d5d67114fdf64b77622` |
| `verify_gravity_600cell_dust_homogeneous_two_by_two.py` | `b97793c99ad2a24d5fd744f6a2e029b8fb51a40c632598c3860aeea602f6c816` |
| `gravity_600cell_dust_homothetic_canonical_lapse.json` | `4b1c59c0518eec11b88b140cdecdf558d762c0d70b4826a758f67544e14ac5b9` |
| `gravity_600cell_dust_second_tick_local_correction.json` | `936984bc84a714140ce16917ee559b346b3c0d4a5ba92d8fb723398a120f8e70` |
| `verify_gravity_600cell_dust_second_tick_local_correction.py` | `cef59fa0bc3a1c8fa3be0193234371b7dda303a0ec72683ddcdd88bcb40f3725` |
| `verify_gravity_global_regge_orbits.py` | `ad93cdd08fabeeee56b009f23936696837c4362f88ae23f92a36d0395e61ffaf` |
| `verify_gravity_600cell_dust_full_anisotropic_legendre_rank.py` | `834b97c85e386def853b6308e65e831c52d62d7cbcc4b23118602120d6c676e5` |
| `verify_gravity_600cell_dust_full_boundary_tangent.py` | `c8662bb0835865aac6696fc3f474ed668fed3fe393b9c32a59e709a984c35571` |
| `verify_gravity_600cell_dust_internal_curvature_response.py` | `276982879fae5f8fa735f27a6fa30bfe965dc3e41c169d8a229a61c23511ae66` |

Both independently derived schedule parities and all four previously frozen
derivative variants are mandatory.  Arithmetic remains at 100 decimal digits
with 80-decimal Flint balls.  Only the trivial binary-tetrahedral sector is
needed; the full suite is excluded.

## Slab reconstruction

For each parity take the committed states, without averaging parities:

```text
slab 1:
  q_old      = exp(2 a_0) L_0^2,       a_0 = 0,
  q_new      = exp(2 a_1) L_0^2,
  q_diagonal = exp(a_0+a_1) L_0^2-rho_1,
  q_pole     = -rho_1;

slab 2:
  q_old      = exp(2 a_1) L_0^2,
  q_new      = exp(2 a_2) L_0^2,
  q_diagonal = exp(a_1+a_2) L_0^2-rho_2,
  q_pole     = -rho_2.
```

The total dust mass is the same committed `MASS` in both Hessians.  Replacing
`rho_2` by a mass rescaled from `L_2` is forbidden.  Require every base and
displaced simplex to remain on the certified one-negative-direction branch.

Independently verify the complete thirty-component seam between the committed
first-slab post-momentum and second-slab pre-momentum under the literal orbit
map.  Failure of this background control prevents any transport conclusion.

## Operators and fixed carrier

Use the same preregistered basis

```text
U = [(1_30,0_30)/sqrt(30), (0_30,1_30)/sqrt(30)].
```

For each slab reconstruct its `F_n : C^60 -> C^160` exactly as in the frozen
parent verifier.  For slab 1 also reconstruct the action-generated `T_1`.
Require independently for each `F_n`:

```text
full rank/nullity       59 / 1,
rank(F_n U)              1,
unique kernel line K_n   contained in U.
```

The kernel vector is the smaller high-precision eigendirection of
`(F_n U)*(F_n U)`; no 60-dimensional eigensystem and no target-dependent
rotation is allowed.

For every derivative variant form

```text
v_1       = U k_1,
v_2       = U k_2,
w_1       = T_1 v_1,
d_line    = sin angle(w_1,v_2),
r_curve   = ||F_2 w_1||_2 / ||w_1||_2.
```

The line distance is the primary test.  The fresh second-slab curvature
response is an independent diagnostic.

## Frozen uncertainty and labels

Every uncertainty is the sum of:

1. operational-primary versus operational-shadow change;
2. validation-primary versus validation-shadow change;
3. operational-primary versus validation-primary change;
4. propagated Flint radii for `F_1`, `F_2` and `T_1`;
5. kernel angular errors obtained by dividing the operator error by the
   resolved `F_n U` singular gap;
6. a dimension-scaled 100-decimal serialization floor;
7. `1e-70`.

For the bounded line distance:

```text
TRANSPORT_IDENTIFIED   d_line <= 10 epsilon and epsilon < 1e-2,
TRANSPORT_SEPARATED    d_line > 100 epsilon,
NUMERICALLY_OPEN       otherwise.
```

Classify `r_curve` by the analogous `ZERO/NONZERO/NUMERICALLY_OPEN` rule.  It
cannot override the primary line label, but disagreement is a control failure
because `K_2` is required to be the unique kernel of `F_2`.

Compare the two schedule `K_1`, `K_2` and transported lines using the literal
orbit-set permutation.  Because that permutation fixes `U`, their coordinates
in the `(q,p)` basis may also be compared directly.  Any resolved parity
difference outranks a transport hit.

## Frozen outcome hierarchy

1. `CURVATURE_KERNEL_TRANSPORT_CONTROL_FAILED` for any provenance, seam,
   branch, determinant, rank, uniqueness or diagnostic-consistency failure.
2. `CURVATURE_KERNEL_TRANSPORT_SCHEDULE_DEPENDENT` for a resolved parity
   difference in `K_1`, `K_2`, `T_1K_1` or the transport label.
3. `CURVATURE_KERNEL_TRANSPORT_IDENTIFIED` only when both schedules identify
   `T_1K_1=K_2` and the fresh `F_2T_1K_1` response is zero.
4. `CURVATURE_KERNEL_TRANSPORT_REFUTED` when both schedules resolve the line
   separation and the fresh curvature response is nonzero.
5. `CURVATURE_KERNEL_TRANSPORT_NUMERICALLY_OPEN` otherwise.

## Interpretation firewall

A positive result establishes only a transported, background-dependent
curvature-blind tangent line.  It does not establish a first-class constraint,
gauge symmetry, lapse, clock or physical time direction.  A negative result
closes only this candidate distribution; it does not prove absence of every
possible pseudo-constraint.

No later tick, nonlinear continuation, refinement, continuum dispersion or
full-suite run is included.
