# Preregistered protocol: internal curvature response of the 119 strong modes

Date: 2026-08-17

Prior-art gate commit: `564b3a9`.

This protocol was committed before assembling or evaluating any deficit-angle
response on an expanding or contracting tangent subspace.  The geometric
triangle census in the prior-art note was allowed because it did not load a
tangent eigenvector.

## Frozen inputs

| input | SHA-256 |
|---|---|
| `gravity_600cell_dust_homothetic_canonical_lapse.json` | `4b1c59c0518eec11b88b140cdecdf558d762c0d70b4826a758f67544e14ac5b9` |
| `gravity_600cell_dust_full_boundary_tangent.json` | `4da8bcd2890a54bc9d3b60c6195df2933ea56194d942ab0285b51599ba287bd5` |
| `gravity_600cell_dust_full_boundary_tangent.npz` | `816c605da2a655442bbadce7a23965f0822f99e7bdc1d0a4a27af548de85446b` |
| `verify_gravity_600cell_dust_full_boundary_tangent.py` | `c8662bb0835865aac6696fc3f474ed668fed3fe393b9c32a59e709a984c35571` |
| `verify_gravity_600cell_dust_full_anisotropic_legendre_rank.py` | `834b97c85e386def853b6308e65e831c52d62d7cbcc4b23118602120d6c676e5` |
| `verify_gravity_global_regge_orbits.py` | `ad93cdd08fabeeee56b009f23936696837c4362f88ae23f92a36d0395e61ffaf` |

The verifier must reject any provenance mismatch.  It uses the required
interpreter `/home/razvan/science/.venv/bin/python`, 100-decimal `mpmath`,
80-decimal Flint complex balls and the same four frozen finite-difference
steps as the full tangent audit:

```text
operational primary    1e-20
operational shadow     1e-15
validation primary     3e-20
validation shadow      3e-15
```

## Construction fixed before the result

### 1. Triangle carrier and real-branch deficit coordinates

For each schedule parity, reconstruct all 6,240 triangle hinges and their
free `2T` orbits directly from the slab incidence data.  Require:

- 260 orbits of size 24;
- 100 boundary-spacelike orbits;
- 100 internal-spacelike orbits;
- 60 internal-timelike orbits;
- no null triangle;
- the causal and boundary label constant on every orbit.

Only the 160 internal orbit types enter the response.  For each internal
triangle use

```text
kappa_h = -i (2*pi + sum theta_s,h)  when area-square(h) > 0,
kappa_h =     2*pi + sum theta_s,h   when area-square(h) < 0.
```

The row phase is held fixed during differentiation.  It must make the base
deficits and their four Jacobians real to the already-used `1e-70` arithmetic
gate.  Since the phase is invertible, zero/nonzero and rank are unchanged by
it.

### 2. Incidence-derived curvature Jacobian

For each local 4-simplex pattern and derivative variant, reuse the audited
Lorentzian angle derivative

```text
partial theta_(s,h) / partial log(q_e).
```

Sum it over the literal 4-simplices incident on each internal triangle.  No
area weight, dual-volume weight, fitted coefficient or boundary angle is
allowed.  This gives

```text
D_kappa : C^2280 -> R^3840.
```

As an independent equivariance control, assemble the sparse direct derivative
at all 3,840 internal triangles and compare every populated row/column entry
with the representative double-orbit kernel using
`g_row^-1 * g_column`.  The maximum residual must be below `1e-70` for all
four variants.

Project the representative kernel into each of the seven deterministic
minimal `2T` sectors.  In a sector of irrep dimension `d`, the block has shape

```text
(160 d) x (95 d).
```

### 3. Canonical slab response

Independently reconstruct the committed canonical solve

```text
Y : boundary phase C^(60d) -> (internal,new) edges C^(65d)
```

from the full Hessian in each sector and variant, using Flint balls.  Stack
the literal old-boundary position projection above it:

```text
Z = ([I,0], Y)^T : C^(60d) -> slab edges C^(95d).
```

The internal-curvature response on boundary phase space is

```text
F = D_kappa Z : C^(60d) -> C^(160d).
```

All pre-Legendre Flint determinant balls must exclude zero.  No tangent
eigenvector enters this construction.

### 4. Strong subspaces frozen from the blind tangent census

For each tangent matrix and each derivative variant:

- identify the trivial sector mechanically by constant-vector overlap above
  `1/2`;
- use `k = 5d` in every nontrivial sector and `k = 4` in the trivial sector;
- define `E_plus` as the ordered complex-Schur subspace of the `k` largest
  eigenvalue moduli;
- define `E_minus` as the ordered complex-Schur subspace of the `k` smallest
  eigenvalue moduli;
- independently form the corresponding direct-eigenvector spans and compare
  them with the Schur spans;
- require selected count `k` and a selected/unselected modulus ratio above 2
  for every strong selection.

With regular-representation multiplicity `d`, the frozen count is

```text
sum_nontrivial d*(5d) + 1*4 = 119
```

for each branch.  The unresolved fifth homogeneous reciprocal pair is
reported separately and is forbidden from entering the strong verdict.

### 5. Rank questions and exact attempt count

There are exactly:

```text
14 full phase-space maps = 2 parities x 7 sectors,
28 strong restrictions  = 2 parities x 7 sectors x 2 branches.
```

For every variant compute singular values of

```text
F,
R_plus  = F E_plus,
R_minus = F E_minus,
```

with orthonormal Schur bases.  For the strong restrictions also compute the
singular spectra using the independently orthonormalized direct-eigenvector
spans.

For each singular spectrum define a conservative absolute uncertainty from:

1. operational-primary versus operational-shadow change;
2. validation-primary versus validation-shadow change;
3. operational-primary versus validation-primary change;
4. Schur-versus-direct spectrum change for a strong restriction;
5. the Flint response radius propagated by `||D_kappa||`;
6. the tangent-ball subspace bound using eigenvector condition and selected
   spectral separation;
7. a dimension-scaled binary64 SVD floor;
8. the fixed `1e-70` arithmetic floor.

The same uncertainty applies to every singular value in that one spectrum.
Classify a singular value mechanically as:

```text
ZERO             sigma <= 10 epsilon,
RESOLVED_NONZERO sigma > 100 epsilon,
NUMERICALLY_OPEN otherwise.
```

Thus:

- a restriction is `ZERO` only if every singular value is `ZERO`;
- it is `INJECTIVE` only if every singular value is `RESOLVED_NONZERO`;
- otherwise its resolved rank, zero count and open count are reported without
  reinterpretation.

Operator norms and response gains are printed, but because no physical
curvature norm has been derived, their magnitudes and plus/minus ordering are
**PATTERN**, never an acceptance gate.

## Preregistered outcome hierarchy

Apply the first matching outcome:

1. `INTERNAL_CURVATURE_RESPONSE_CONTROL_FAILED` if provenance, carrier,
   branch, reality, equivariance, Flint determinant or strong-selection
   controls fail.
2. `FULL_BOUNDARY_PHASE_CURVATURE_INJECTIVE` if all 14 full maps and all 28
   strong restrictions are injective.
3. `STRONG_TANGENT_CURVATURE_INJECTIVE` if all 28 strong restrictions are
   injective but at least one full map is not.
4. `STRONG_TANGENT_CURVATURE_ZERO` if all 28 strong restrictions are zero.
5. `STRONG_TANGENT_CURVATURE_PARTIAL_OR_OPEN` for every mixed, deficient-rank
   or numerically open result.

Outcomes 2 or 3 are a **DERIVED COMPUTATIONAL NEGATIVE** for the exact claim
that the 119 strong directions are curvature-preserving lapse/gauge modes.
They are not a proof that the modes are physical gravitons.  Outcome 4
supports, but does not prove, a pseudo-gauge interpretation on the curved
background.  Outcome 5 remains honest and carries the complete rank ledger.

## Explicit exclusions

- no comparison with a continuum harmonic spectrum;
- no fit to a desired speed, growth rate or Planck scale;
- no interpretation of a one-step eigenvalue as a Lyapunov exponent;
- no boundary intrinsic-curvature observable in this mission;
- no Ricci/Weyl, constraint-satisfying or polarization classification;
- no full-suite run: only the newly registered targeted verifier will run.
