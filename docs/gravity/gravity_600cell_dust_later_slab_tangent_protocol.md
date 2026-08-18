# Preregistered blind protocol: third-slab tangent and two-update product

Date: 2026-08-18

Prior-art gate commit: `920ce5d`.

This protocol is blind to the earlier centered Jacobi spectrum, the
conformal/shape split and all negative-mode artifacts.  No `T_3` matrix,
product spectrum or later Jacobi coefficient has been evaluated before this
commit.

## Frozen inputs

| input | SHA-256 |
|---|---|
| `gravity_600cell_dust_two_step_full_tangent.json` | `f7fbf18535cc00dacec9a9ffa95f97f2d1847ac83073f27d39fcdb7968b0bafc` |
| `gravity_600cell_dust_two_step_full_tangent.npz` | `ce78ebf415584b1cdcf1d2cb07687135b624ad4939e0a4e54650653f7b384e6d` |
| `verify_gravity_600cell_dust_two_step_full_tangent.py` | `c1a3fb09146188c1932ab81629ab69817f2a2f19108fdf8d9e89d78b6de8f717` |
| `gravity_600cell_dust_second_tick_local_correction.json` | `936984bc84a714140ce16917ee559b346b3c0d4a5ba92d8fb723398a120f8e70` |
| `gravity_600cell_dust_third_tick_local_correction.json` | `ebf2f1a11b9a4e9c76fb1ce33066c0782429cf6500770df7bbe4d92de4a050c0` |
| `verify_gravity_600cell_dust_third_tick_local_correction.py` | `d0a1bf9bd3beb9e7ed3c805a12dcd63c7c83db977a3a1d2cd127d9b9c492a7a9` |
| `verify_gravity_600cell_dust_full_boundary_tangent.py` | `c8662bb0835865aac6696fc3f474ed668fed3fe393b9c32a59e709a984c35571` |
| `verify_gravity_global_regge_orbits.py` | `ad93cdd08fabeeee56b009f23936696837c4362f88ae23f92a36d0395e61ffaf` |
| `verify_gravity_600cell_dust_full_anisotropic_legendre_rank.py` | `834b97c85e386def853b6308e65e831c52d62d7cbcc4b23118602120d6c676e5` |

Require the upstream outcomes and pass counts

```text
TWO_STEP_FULL_TANGENT_COCYCLE_CERTIFIED  16/16,
SECOND_HOMOTHETIC_TICK_ACCEPTED           6/6,
THIRD_HOMOTHETIC_TICK_ACCEPTED            6/6.
```

Both later tick artifacts must declare fixed mass and
`mass_recomputed_from_later_scale=false`.  The existing two-map archive must
contain exactly `448` arrays and reproduce its recorded SHA-256.

Arithmetic remains at `100` decimal digits for midpoint construction and
`80` decimal digits for Flint balls.  Every stored binary64 midpoint from the
`T_2` archive is re-enclosed with the stored radius plus one half-ULP in each
real and imaginary component.

## Exact third-slab geometry and seam

For each schedule parity use its own committed values without averaging:

```text
a_old = a_2,
a_new = a_3,
r     = r_3,

q_old      = exp(2 a_2) L_0^2,
q_new      = exp(2 a_3) L_0^2,
q_diagonal = exp(a_2+a_3) L_0^2-rho_3,
q_pole     = -rho_3.
```

The dust Hessian uses the original conserved total mass.  Re-derive the
literal old-to-new boundary map from the slab-3 geometry.  Before any tangent
is accepted, require the complete thirty-component equality

```text
p_pre,3[i] = p_post,2[mapping[i]]
```

inside the larger of the committed seam bounds.  Also require all base and
displaced Lorentzian branch, reality, carrier and local derivative controls.

## Complete tangent and product

Derive the same seven deterministic minimal sectors from the literal group
action.  Their irrep dimensions must be

```text
3,2,2,2,1,1,1,
```

and the restored full phase dimension must be

```text
sum 60 d^2 = 1440.
```

For every parity, sector and all four frozen derivative variants:

1. project the complete slab-3 action Hessian;
2. require the full pre-Legendre determinant ball to exclude zero;
3. solve the canonical response and form `T_3`;
4. require `T_3^* Omega T_3-Omega` to contain zero entrywise;
5. re-enclose the matching committed `T_2` ball;
6. form the rigorous product `C_32=T_3T_2`;
7. require `C_32^* Omega C_32-Omega` to contain zero entrywise.

No earlier negative carrier, fitted projector or selected Schur block may
enter.

## Blind diagnostics and schedule comparison

For `T_3` and `C_32` separately record the same target-free diagnostics as
upstream:

- singular values and determinant modulus;
- symplectic and reciprocal-singular defects;
- eigenvalues, reciprocal-conjugate defect and conditioning;
- unit-consistent, open and resolved-off-unit counts with representation
  multiplicity restored;
- minimum and maximum resolved moduli.

These branch counts are diagnostics only.  Their values do not affect the
outcome and are not called Lyapunov exponents.

For each of seven sectors and both maps compare the two schedule parities by
their ordered singular spectra.  Use the inherited calibrated uncertainty and
labels

```text
SCHEDULE_ROBUST       distance <= 10 epsilon,
SCHEDULE_DEPENDENT    distance > 100 epsilon,
SCHEDULE_OPEN         otherwise.
```

Optimally matched eigenvalue distance is secondary and uses its full
conditioning envelope.

The deterministic output archive stores both maps:

```text
2 schedules * 7 sectors * 4 variants
* 2 maps * 4 fields = 448 arrays.
```

## Frozen outcome hierarchy

1. `LATER_SLAB_TANGENT_CONTROL_FAILED` for any provenance, seam, branch,
   carrier, sector, archive, determinant or canonicality failure.
2. `LATER_SLAB_TANGENT_SCHEDULE_DEPENDENT` if any primary singular-spectrum
   comparison is resolved dependent, or a resolved eigenvalue dependence is
   not contradicted by conditioning.
3. `LATER_SLAB_TANGENT_SCHEDULE_OPEN` if none is dependent but at least one
   primary singular comparison is open.
4. `LATER_SLAB_TANGENT_CERTIFIED` only if all fourteen primary comparisons
   are robust and every structural control passes.

The outcome is independent of how many eigenvalues are inside or outside the
unit circle.

## Explicit exclusions

- no shifted Jacobi construction in this first artifact;
- no comparison with the old conformal, shape or negative carriers;
- no physical mode labels, frequency, energy norm or limiting speed;
- no Lyapunov claim from a two-update product;
- no tick four, refinement or nonlinear anisotropic solve;
- no full-suite run.

