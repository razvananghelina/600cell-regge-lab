# Adversarial protocol: dense second finite-height tangent and two-step map

Date: 2026-08-22

Primary artifact commit: `2a80690`.

Status: **FROZEN BEFORE IMPLEMENTING OR EXECUTING THE DENSE SECOND-SLAB
HESSIAN, FULL PRE-LEGENDRE SOLVE OR DENSE TWO-STEP PRODUCT.**

Only the new targeted verifier may run.  The full suite is forbidden.

## 1. Frozen provenance

```text
docs/gravity/gravity_600cell_second_full_boundary_tangent_prior_art.md
  d3740e0b08b2f3ec6adf2c69c762e5e5dc0cdd87a571d6d27bc62e78518e70be

docs/gravity/gravity_600cell_second_full_boundary_tangent_protocol.md
  c6615e5011f0a07e5ddaccd00c63d7ed9419e4058451afc3e2243423098c7024

docs/gravity/gravity_600cell_second_full_boundary_tangent_protocol_correction.md
  cd8b99a64e3c88c91781188db0e08e5f8bf2cfc1d450e305949f7350288cda80

docs/gravity/gravity_600cell_second_full_boundary_tangent_history_provenance_correction.md
  fc321949e07de6ec551743abf5274c337c9ed0adb8628142fa1e83f218f1164b

docs/gravity/gravity_600cell_second_full_boundary_tangent_first_run_failure.md
  d9cdc909aeddb8a01b0051420c2dead72007905d3e8420fe259b843d1aefd0ef

docs/gravity/gravity_600cell_second_full_boundary_tangent_second_run_failure.md
  58eda42f4babb4ea3c47d3ade1d5e5d327595ce1edb55127c4e488229bf804a6

reproducible/verify_gravity_600cell_finite_height_second_full_boundary_tangent.py
  8c29c66eb4ec253229685cdbe56eb0371fb00860f0f2e15804fca0c7c64ec536

reproducible/gravity_600cell_finite_height_second_full_boundary_tangent.json
  f97db13031fd366b74e7d327abf61d8d23c24ee2889e500a2cfc747ed7dd6990

reproducible/gravity_600cell_finite_height_second_full_boundary_tangent.npz
  c80a1af93deddc9526c362a8a76eec5dfc4b8360440a1a69f73ca61a419aac9a

reproducible/verify_gravity_600cell_finite_height_full_boundary_tangent_adversarial.py
  6e4fdbe6822ac024d2bad9aa22fe73f85aa9805a018146482ace47cfbefc43b6

reproducible/gravity_600cell_finite_height_full_boundary_tangent_adversarial.json
  ee9491b2ae5fdf3f2a9d0d78c0e837c8c2692797d87ccd8e1757efeadd8060e7

reproducible/verify_gravity_600cell_finite_height_full_boundary_tangent.py
  c4e60d6ef87131d87a93b64d5381d16d8de8d3990340efd5405ec983f64db94d

reproducible/gravity_600cell_finite_height_full_boundary_tangent.npz
  0c34f179821f9d0b74de4906051bbcb7149b4e79881410ea662241adc0aa19bf

reproducible/gravity_600cell_finite_height_composition_adversarial.json
  d50e87f736e51585596aa1d7778238febaf7422840d668499878d8bd917f99e9

reproducible/verify_gravity_600cell_finite_height_composition_adversarial.py
  8395e921ab1c1f518abb567a114f1eb8bfdf2068be031bff55c8d2f0cff56c2b

reproducible/verify_gravity_600cell_dust_full_anisotropic_legendre_rank.py
  834b97c85e386def853b6308e65e831c52d62d7cbcc4b23118602120d6c676e5

reproducible/verify_gravity_global_regge_orbits.py
  ad93cdd08fabeeee56b009f23936696837c4362f88ae23f92a36d0395e61ffaf
```

Require the primary outcome to remain
`TWO_STEP_FULL_BOUNDARY_TANGENT_SCHEDULE_ROBUST_PRIMARY` with `31/31`, 672
numeric arrays and the frozen archive hash.  Require the old dense first-slab
control to retain its `22/22` adversarially replicated outcome.  Hashes may be
checked before the dense decision; no primary tangent entry, sector result or
schedule comparison may be read until section 8.

## 2. Mechanically different route

The dense decision must not use:

- the primary 95-by-95 orbit-convolution kernels;
- a binary-tetrahedral sector projection;
- a minimal-block determinant;
- a primary rank, scale-lift or schedule label;
- a primary tangent entry or singular/eigenvalue target.

Instead assemble the complete real 2280-edge Hessian directly from all 2400
four-simplices and 6240 triangles, form each complete 1560-by-1560
pre-Legendre matrix, and solve all 1440 right-hand sides in physical edge
coordinates.  The decisive labels are assigned in this real-space basis.

## 3. Independent history and three slab assemblies

At 120 decimal digits reconstruct `q1,h1,r1,m1,pi1,q2,h2,r2,m2,pi2` from the
equations and brackets in the primary protocol.  Compare only after solving
with the frozen 180-digit adversarial branch-B record below `1e-65`.

For each staircase parity assemble three independent actions:

1. first physical slab: `old=1`, `new=r1^2`, `rho=h1^2`, mass `m0`;
2. second normalized slab: `old=1`, `new=r2^2`, `rho=h2^2`, mass `m1`;
3. second direct physical slab: every squared length in item 2 multiplied by
   `c=r1^2`, mass `m0`.

Use the dense route's independently established centered logarithmic steps

```text
1e-18, 5e-19, 2.5e-19, 1.25e-19
```

and the three adjacent Richardson levels.  Require the same Lorentzian
inertia and angle-argument controls.  The physical second evaluator must rerun
every local derivative and complete assembly; it may not multiply a stored
normalized matrix.

Process stages sequentially and release no-longer-needed 2280-edge Hessians so
the replication does not depend on excessive simultaneous memory.

## 4. Raw reciprocity and scale lift

Before the unique licensed symmetrization `H -> (H+H^T)/2`, require every raw
Hessian to satisfy the dense Richardson-plus-roundoff reciprocity envelope.

At all four raw and all three Richardson levels compare

```text
H2_physical = r1^2 H2_normalized.
```

The uncertainty is the sum of both adjacent-level variations, the direct
assembly roundoff envelopes and a conditioning-aware binary64 term.  Use the
same frozen classifier:

```text
AGREES   distance <= 10 uncertainty
REFUTED  distance > 100 uncertainty
OPEN     otherwise.
```

Any refuted level gives `SCALE_LIFT_REFUTED`; no refuted level and at least
one open level gives `SCALE_LIFT_OPEN`; otherwise the dense lift is confirmed.

## 5. Full rank and canonical maps

For every stage, parity and Richardson level form the full pre-Legendre matrix
and compute both `gesvd` and `gesdd` singular values.  A stage is `REGULAR`
only if every normalized smallest singular value exceeds 100 times the sum of
Richardson, driver, assembly and binary-roundoff uncertainties.

For every regular block solve all 1440 columns, retain the exact physical
final-to-old edge shift, and test the three real symplectic block identities.
Require the directly solved physical second tangent to agree with

```text
D_c T2_normalized D_c^-1,  D_c=diag(I,c I), c=r1^2,
```

under the 10/100 classifier.  Require the identity lift and the `r1` momentum
lift to be refuted in at least one complete real map.  Omitting `K_NO` must
also fail the symplectic gate.

## 6. Dense schedule and two-step products

Before opening a primary tangent entry:

1. compare the even and odd direct physical second tangents in their common
   lexicographic 720-edge basis;
2. form all four dense products
   `T2_physical[p2] T1_physical[p1]` at all three Richardson levels;
3. check every product's symplectic identities;
4. compare all six pairs of the four products;
5. require a common `1e-3` synthetic corruption to be detected;
6. require identity and `r1` hostile scale lifts to be refuted in at least one
   product.

Use the same 10/100 schedule classifier and propagate both step variations,
assembly envelopes and binary roundoff.  No tangent eigenvalue or singular-
value spectrum is computed.

## 7. First-slab control

The dense route recomputes the first physical tangent rather than importing a
primary map.  Its three complete array hashes must reproduce the already
accepted dense first-slab hashes.  This is a known-answer implementation
control, not a target for the second slab.

## 8. Delayed primary reconciliation

Freeze in memory all dense rank, canonicality, scale-lift, second-schedule and
two-step-schedule labels.  Only then open the two primary NPZ archives.

Construct the deterministic minimal bases solely for closure.  Project each
dense `H12` map from the lexicographic physical edge basis into the matching
minimal basis, then compare entries with:

- all 14 primary second physical sector maps;
- all 28 primary two-step product sector maps.

Use primary ball radii, dense adjacent-level variations and binary roundoff in
the 10/100 classifier.  These projections may confirm or refute the primary
result but may not alter a dense label.  No spectrum is licensed.

## 9. Outcome hierarchy

Assign exactly one:

1. `SECOND_FULL_TANGENT_DENSE_CONTROL_FAILED`;
2. `SECOND_FULL_TANGENT_DENSE_SCALE_LIFT_REFUTED`;
3. `SECOND_FULL_TANGENT_DENSE_SCALE_LIFT_OPEN`;
4. `SECOND_FULL_TANGENT_DENSE_RANK_OPEN`;
5. `SECOND_FULL_TANGENT_DENSE_CANONICALITY_FAILED`;
6. `SECOND_FULL_TANGENT_DENSE_SCHEDULE_DEPENDENT`;
7. `SECOND_FULL_TANGENT_DENSE_OPEN`;
8. `SECOND_FULL_TANGENT_PRIMARY_REFUTED`;
9. `TWO_STEP_FULL_BOUNDARY_TANGENT_SCHEDULE_ROBUST_ADVERSARIALLY_REPLICATED`.

Outcome 9 requires every dense gate and every delayed primary comparison to
agree.  It remains a two-step linearized canonical response, not a physical
mode spectrum, graviton, stability theorem, wave equation, limiting speed,
`G`, Planck scale or particle result.
