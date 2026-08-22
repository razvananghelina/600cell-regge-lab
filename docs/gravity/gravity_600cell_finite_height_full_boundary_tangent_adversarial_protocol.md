# Adversarial preregistration: dense real-space finite-height canonical map

Date: 2026-08-22.

Primary result commit: `88833b0`.

Status: **FROZEN BEFORE THE FIRST DENSE REAL-SPACE PRE-LEGENDRE MATRIX,
SINGULAR VALUE, SOLVE, SYMPLECTIC DEFECT OR PARITY COMPARISON.**

The disclosed primary result is

```text
21/21 PASS,
both parities 7/7 REGULAR and 7/7 CANONICAL,
FINITE_HEIGHT_FULL_BOUNDARY_TANGENT_SCHEDULE_ROBUST_PRIMARY.
```

This protocol is allowed to try to refute that result.  No primary rank,
condition number, tangent entry, singular value or parity distance may be
used to set an adversarial step, tolerance or outcome gate.

## 1. Frozen inputs

```text
docs/gravity/gravity_600cell_finite_height_full_boundary_tangent_prior_art.md
  6fe3e10daf97fd60849a837e56716ced594e19c77117ecc14f862822edc10074

docs/gravity/gravity_600cell_finite_height_full_boundary_tangent_protocol.md
  373cfd80a6e41993157e240313874de47317436839bcdccb7d5ae79b78855235

docs/gravity/gravity_600cell_finite_height_full_boundary_tangent_first_failure.md
  45533840dadfa37f64c7688cd5a09de335ead20a63264ac45af428308e85fcdc

reproducible/verify_gravity_600cell_finite_height_full_boundary_tangent.py
  c4e60d6ef87131d87a93b64d5381d16d8de8d3990340efd5405ec983f64db94d

reproducible/gravity_600cell_finite_height_full_boundary_tangent.json
  266638aeaa825b327b63a84eda36a499456dc4b4f9a86f964cee5f79d6d6e930

reproducible/gravity_600cell_finite_height_full_boundary_tangent.npz
  0c34f179821f9d0b74de4906051bbcb7149b4e79881410ea662241adc0aa19bf

reproducible/gravity_600cell_finite_height_carrier_quadratic_adversarial.json
  54915cf364c36af6bbc8e1dbd36433079269d293453478bfdf589e547d462ad6

reproducible/verify_gravity_600cell_dust_full_anisotropic_legendre_rank.py
  834b97c85e386def853b6308e65e831c52d62d7cbcc4b23118602120d6c676e5

reproducible/verify_gravity_global_regge_orbits.py
  ad93cdd08fabeeee56b009f23936696837c4362f88ae23f92a36d0395e61ffaf
```

Require the primary artifact to retain its disclosed outcome and hashes,
and the independently accepted finite-height quadratic artifact to retain
`18/18`.  The primary numerical archive remains unopened until the dense
rank, canonicality and schedule labels have all been frozen in memory.

The adversarial verifier must be registered exactly once with no registry
duplicates.  Run only it and the directly imported 43 geometry controls.

## 2. Mechanically different decisive construction

The primary route assembled only one identity-group Hessian row per regular
`2T` orbit at 180 digits, reconstructed a group convolution, projected it to
seven complex minimal sectors and solved 42 Flint ball systems of sizes at
most `195`.

The adversarial route is forbidden to call any of those decisive functions:

```text
high_precision_sector_bases,
assemble_full_representative_kernels,
project_full_kernel,
canonical_tangent_ball.
```

Instead, for each parity it must:

1. assemble every local contribution directly into four dense real
   `2280 x 2280` Hessians in physical edge coordinates;
2. form three dense Richardson Hessians;
3. extract and classify the full real `1560 x 1560` pre-Legendre matrix;
4. solve the full real system with `1440` right-hand sides;
5. construct the full `1440 x 1440` boundary tangent;
6. test the real symplectic block identities without representation theory;
7. compare the two complete physical-coordinate maps directly.

This shares the audited local Regge area/angle primitives and exact carrier,
but not the primary decisive assembly, group reduction, arithmetic type,
linear-system size or canonicality calculation.

## 3. Independent precision and dense Hessian hierarchy

Use 120 decimal digits for local Lorentzian geometry and the independently
frozen centered logarithmic angle steps

```text
h0 = 1e-18,
h1 = 5e-19,
h2 = 2.5e-19,
h3 = 1.25e-19.
```

Convert local contributions to binary64 only when inserting them into the
dense ambient matrices.  Let these be `H0,H1,H2,H3` and form

```text
H01=(4*H1-H0)/3,
H12=(4*H2-H1)/3,
H23=(4*H3-H2)/3.
```

Require every base and displaced simplex to retain the accepted Lorentzian
branch, positive leading-minor and angle-argument margins, and the local
entrywise step hierarchy.  Require the final physical dense matrices to be
real under the existing assembly control.

For each level set

```text
N_H=max(1,||H12||_F),
e_step,H=max(||H01-H12||_F,||H12-H23||_F)/N_H,
e_round,H=2*(bound_even+bound_odd)/N_H+100*eps_binary64,
e_H=e_step,H+e_round,H.
```

Here `bound_parity` is the directly assembled forward-summation
antisymmetry bound.  Require

```text
max_level ||H-H^T||_F/N_H <= 10*e_H.
```

Only after this raw reciprocity gate passes, replace each Hessian by
`(H+H^T)/2` for the dense solve.  This is the unique orthogonal projection
onto the exact Hessian symmetry and is frozen before seeing any rank or map.

## 4. Full real pre-Legendre rank

In physical edge order `O=720`, `X=840`, `N=720`, form at every level

```text
J = [[ K_XX,  K_XN],
     [-K_OX, -K_ON]].
```

Set

```text
N_J=max(1,||J12||_F),
e_step,J=max(||J01-J12||_F,||J12-J23||_F)/N_J,
e_svd=max(max_abs(svd_gesvd-svd_gesdd)/N_J,
          10*eps_binary64*max(1,||J12||_2)/N_J),
e_J=e_step,J+e_svd+e_H.
```

Require all three full real SVDs to satisfy

```text
sigma_min(J)/N_J > 100*e_J.
```

Record the complete singular spectra and condition estimates.  `REGULAR`
requires this full-space gap; any other case is `NUMERICALLY_OPEN`, not a
proof of singularity.

## 5. Direct full tangent and real symplectic identities

Only after a parity is `REGULAR`, solve for all `1440` columns:

```text
R = [[-K_XO,0],
     [ K_OO,I]],
Y = J^-1 R,

T_raw = [[Y_N],
         [[K_NO,0]+K_NX*Y_X+K_NN*Y_N]].
```

Derive the final-to-old permutation from all physical edge labels.  Apply it
to both output configuration and output momentum, then reorder both input
and output by the common lexicographic list of 720 old physical edges.  No
orbit or representation data is used.

For every level partition `T=[[A,B],[C,D]]` and evaluate the full real
identities

```text
A^T C-C^T A = 0,
B^T D-D^T B = 0,
A^T D-C^T B = I.
```

Combine their Frobenius norms into `r_sym`.  Define `e_sym` from both
adjacent tangent differences propagated as

```text
2*max(1,||T12||_2)*delta_T+delta_T^2,
```

plus the dense Hessian roundoff contribution, the two adjacent defect
differences and `100*eps_binary64*max(1,||T12||_2)^2`.  Require
`r_sym<=10*e_sym`.

As hostile controls on the `K12` even solve:

1. omit the direct `K_NO` term and require the symplectic defect to exceed
   `100*e_sym`;
2. cyclically shift only the output boundary edge labels and require the map
   to change above `100` times its direct comparison uncertainty;
3. reverse the pre-momentum identity sign in a separately solved exact
   scalar quadratic control and require disagreement with the known map.

No failed hostile control may be replaced.

## 6. Direct schedule verdict before opening the primary archive

For each Richardson level define

```text
N_T=max(1,||T12_even||_F,||T12_odd||_F),
d_level=||T_level,even-T_level,odd||_F/N_T.
```

Define

```text
e_T=max_parity(
      ||T01-T12||_F+||T12-T23||_F
    )/N_T
    +100*eps_binary64*max_parity(max(1,||T12||_2))/N_T
    +e_H.
```

Assign from all three distances:

```text
SCHEDULE_ROBUST     if max(d_level)<=10*e_T,
SCHEDULE_DEPENDENT  if min(d_level)>100*e_T,
SCHEDULE_OPEN       otherwise.
```

The rank, canonicality and schedule labels must be fixed in memory before
reading a primary tangent entry or singular value.

## 7. Post-classification primary reconciliation

Only now open the primary NPZ archive.  For each parity, compute the singular
values of every stored `K12` minimal tangent and repeat that list by the
preregistered irrep dimension `d`.  Their sorted union has length `1440` and
is the representation-theoretic prediction for the singular spectrum of the
full real map.

Compare it with the direct full-map singular spectrum.  The uncertainty is
the sum of:

- the maximum adjacent dense tangent spectral variation;
- the maximum primary Flint-radius Frobenius norm;
- `100*eps_binary64*max(1,cond(J12))*max(1,||T12||_2)`;

all divided by `max(1,||T12||_2)`, plus `1e-15`.  Assign

```text
PRIMARY_AGREES  if normalized distance <=10*uncertainty,
PRIMARY_REFUTED if normalized distance >100*uncertainty,
OPEN            otherwise.
```

This is a unitary-invariant comparison.  No basis alignment or post-result
matching is allowed beyond sorting singular values.

## 8. Outcome hierarchy

1. `FINITE_HEIGHT_FULL_BOUNDARY_TANGENT_ADVERSARIAL_CONTROL_FAILED` for any
   provenance, carrier, branch, reciprocity, physical permutation, scalar or
   hostile-control failure;
2. `FINITE_HEIGHT_FULL_BOUNDARY_TANGENT_ADVERSARIAL_RANK_OPEN` if either full
   real `J` is unresolved;
3. `FINITE_HEIGHT_FULL_BOUNDARY_TANGENT_ADVERSARIAL_CANONICALITY_FAILED` if
   both are regular but either full real map fails the symplectic gate;
4. `FINITE_HEIGHT_FULL_BOUNDARY_TANGENT_PRIMARY_REFUTED` only for a resolved
   direct schedule dependence or a `PRIMARY_REFUTED` singular-spectrum
   comparison, with all controls passing;
5. `FINITE_HEIGHT_FULL_BOUNDARY_TANGENT_ADVERSARIAL_OPEN` for either
   classifier gap;
6. `FINITE_HEIGHT_FULL_BOUNDARY_TANGENT_SCHEDULE_ROBUST_ADVERSARIALLY_REPLICATED`
   only if both full maps are regular and canonical, the direct schedule is
   robust, and both primary singular spectra agree.

No outcome derives a physical perturbation spectrum.  A replicated map is
still one-step local canonical regularity with frozen matter, not a graviton,
wave equation, tick, limiting speed, `G`, Planck scale or particle mass.

Only the new adversarial verifier is run.  The full suite remains forbidden.
