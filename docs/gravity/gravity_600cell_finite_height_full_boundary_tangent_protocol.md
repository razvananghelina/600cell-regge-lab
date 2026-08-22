# Preregistration: finite-height full-boundary canonical tangent

Date: 2026-08-22.

Prior-art gate commit: `c3fc22d`.

Status: **FROZEN BEFORE THE FIRST FINITE-HEIGHT PRE-LEGENDRE RANK,
DETERMINANT, CONDITION NUMBER, FULL TANGENT OR SCHEDULE-MAP COMPARISON.**

The previous finite-height quadratic verifier assembled a complete Hessian
transiently, but it neither formed nor inspected the `1560 x 1560`
pre-Legendre matrix defined below and did not construct a full boundary
map.  No result from that uncomputed matrix is known at freeze time.

## 1. Frozen inputs and source provenance

The verifier must reject any mismatch in these SHA-256 values:

```text
docs/gravity/gravity_600cell_finite_height_full_boundary_tangent_prior_art.md
  6fe3e10daf97fd60849a837e56716ced594e19c77117ecc14f862822edc10074

reproducible/gravity_600cell_finite_height_carrier_quadratic.json
  0ec142bfc68d04498992a6cdba7437933560b860244573d187cb6e018ece78f9

reproducible/gravity_600cell_finite_height_carrier_quadratic_adversarial.json
  54915cf364c36af6bbc8e1dbd36433079269d293453478bfdf589e547d462ad6

reproducible/gravity_600cell_finite_height_internal_carrier_rank.json
  513fdea33f6b868efa6d6f2b2526bade7ce615ea949f955588916a8d0baee0c8

reproducible/gravity_600cell_finite_height_internal_kernel_canonical_reconciliation.json
  81ec0379247023451e82ab42f5beb026ee2d1b083aa5e2553e42b894554266f6

reproducible/gravity_600cell_dust_full_boundary_tangent.json
  4da8bcd2890a54bc9d3b60c6195df2933ea56194d942ab0285b51599ba287bd5

reproducible/gravity_600cell_dust_full_boundary_tangent.npz
  816c605da2a655442bbadce7a23965f0822f99e7bdc1d0a4a27af548de85446b

reproducible/verify_gravity_600cell_dust_full_boundary_tangent.py
  c8662bb0835865aac6696fc3f474ed668fed3fe393b9c32a59e709a984c35571

reproducible/verify_gravity_600cell_dust_full_anisotropic_legendre_rank.py
  834b97c85e386def853b6308e65e831c52d62d7cbcc4b23118602120d6c676e5

reproducible/verify_gravity_global_regge_orbits.py
  ad93cdd08fabeeee56b009f23936696837c4362f88ae23f92a36d0395e61ffaf
```

Require the finite-height primary/adversarial quadratic outcomes to remain
accepted, the restricted internal-rank result and its exact canonical
reconciliation to remain accepted, and the old-background full tangent to
remain `19/19 FULL_BOUNDARY_TANGENT_BLIND_CENSUS_CERTIFIED`.  The last item
is a formula/provenance control only; no old numerical rank or spectrum is a
target for the new background.

The source files are implementation libraries.  Their bottom-level
scientific calculations and old background constants must not execute.

The new verifier must be registered exactly once in `run_all.py`, and the
registry must have no duplicate names.  Only the new targeted verifier and
its directly imported 43 geometry controls may run.

## 2. Independent finite-height background reconstruction

At 180 decimal digits define

```text
z(q)       = (q^2+2)/(2(q^2+3)),
epsilon(q) = 2*pi-5*acos(z(q)),
mu(q)      = 180*epsilon(q)/(pi*sqrt(q^2+4)),
p(q)       = 180*q*epsilon(q)/sqrt(q^2+4)
             -600*sqrt(3)*asinh(q/sqrt(8(q^2+3))).
```

Set `v=3/2`, solve by deterministic bisection on the frozen bracket `(9,10)`

```text
E(q)=4*pi*(mu(q)-mu(v))+q*(p(q)-p(v))=0,
```

and reconstruct

```text
M      = mu(v),
h      = (p(q)-p(v))/(2*pi*mu(q)),
lambda = 1+h*q,
rho    = h^2.
```

Require residual below `1e-140`, bracket width below `1e-150`, agreement
with the committed state below `1e-70`, and

```text
h>0, lambda>0, rho>0, lambda-rho>0.
```

No other root or finite-height history may replace this background.

## 3. Complete carrier and common boundary coordinates

For both frozen staircase parities reconstruct exactly

```text
old O = 720, internal X = 840, new N = 720,
four-simplices = 2400, triangles = 6240,
free 2T edge-orbit types = 30+35+30 = 95.
```

Use logarithmic signed-squared-edge variables and the base values

```text
old       1,
diagonal  lambda-rho,
pole     -rho,
new       lambda^2.
```

Derive the old-to-final boundary map solely from physical labels.  For every
old orbit type and all 24 group elements require the unique relation

```text
{u,v} -> {u+120,v+120}.
```

Require this to cover all 720 edges.  Reorder both parity outputs into the
same lexicographic physical old-edge order before a direct schedule
comparison.  No a posteriori orbit permutation, irrep permutation, phase or
basis rotation is allowed.

Require the two parities to have literally the same order-24 action table.
Construct the seven deterministic minimal right-regular sectors before any
matrix result, with

```text
irrep dimensions       [1,1,1,2,2,2,3],
isotypic dimensions    [1,1,1,4,4,4,9],
sum d^2 = 24,
sum 60*d^2 = 1440.
```

All basis, central-splitter, right-action leakage and conjugate-pair
residuals must be below `1e-140`.  The even and odd basis matrices themselves
must agree below `1e-140`, not only their dimensions or central eigenvalues.

## 4. Hessian kernel and Richardson hierarchy

Use the complete Lorentzian Regge--dust Hessian with the already accepted
branch and boundary term.  At 180 digits evaluate centered logarithmic angle
derivatives at

```text
h0 = 1e-25,
h1 = 5e-26,
h2 = 2.5e-26,
h3 = 1.25e-26.
```

Assemble only the 95-by-95 orbit-convolution kernel, retaining every `oo`,
`ox`, `on`, `xx`, `xn`, `no`, `nx` and `nn` contribution.  Form three
Richardson levels entrywise:

```text
K01=(4*K1-K0)/3,
K12=(4*K2-K1)/3,
K23=(4*K3-K2)/3.
```

Require all base and displaced simplices to retain Lorentzian inertia
`(3,1)`, positive leading-minor and angle-argument margins, the frozen
entrywise step hierarchy, and physical kernel imaginary residue below
`1e-140`.

Project each complete kernel to every deterministic minimal sector.  For a
sector of dimension `d`, the complete Hessian block has size `95d`, while

```text
O=30d, X=35d, N=30d.
```

Require projected Hermitian reciprocity inside ten times the directly
measured Richardson variation plus `1e-135`.  Do not symmetrize a block to
make this gate pass.

## 5. Pre-Legendre regularity classifier

For every parity, sector and Richardson level form

```text
J = [[ K_XX,  K_XN],
     [-K_OX, -K_ON]],

shape(J)=(65d,65d).
```

From `J12`, define a positive diagonal right scaling

```text
D_jj = 1/max(1, ||J12[:,j]||_2).
```

Apply the same invertible `D` to all three levels.  Set

```text
N_J = max(1,||J12*D||_F),
e_step = max(||(J01-J12)*D||_F,
             ||(J12-J23)*D||_F)/N_J,
e_svd = max(max_abs(svd_gesvd-svd_gesdd)/N_J,
            10*eps_binary64*max(1,||J12*D||_2)/N_J),
e_J = e_step+e_svd+1e-135.
```

Here both LAPACK drivers are run on `J12*D`; the driver discrepancy is an
observed numerical floor, not a theorem about the true Hessian.  Binary64
SVD supplies a backward-stable diagnostic gap, not the sole rank
certificate.  Require at all three Richardson levels

```text
sigma_min(J*D)/N_J > 100*e_J.
```

Independently convert each unscaled 180-digit `J` to a 140-decimal Flint
complex ball matrix.  Every determinant ball must exclude zero.  Record
determinant balls, all singular values, condition estimates and the complete
error ledger.

Only the conjunction of the SVD gap and all determinant exclusions receives
`REGULAR`.  Any other case is `NUMERICALLY_OPEN`; this protocol does not call
a determinant ball containing zero a proof of singularity.

## 6. Forced canonical differential

Construct a tangent only after its sector is `REGULAR`.  For each level set

```text
R = [[-K_XO, 0],
     [ K_OO, I]],

Y = J^-1 R,
delta x = Y_X (delta o,delta p_pre),
delta n = Y_N (delta o,delta p_pre),

delta p_post = [K_NO,0]
               +K_NX Y_X+K_NN Y_N.
```

Apply the exact physical final-to-old boundary relabelling to both
`delta n` and `delta p_post`.  Solve in Flint balls; store midpoint and
entrywise radii for every `60d x 60d` minimal tangent.

The complex minimal basis obeys `W^T conjugate(W)=I`.  Hence the correct
within-sector real canonical identity is

```text
T^* Omega T = Omega,
Omega=[[0,I_30d],[-I_30d,0]].
```

For each sector define the symplectic uncertainty from the sum of:

- both adjacent Richardson tangent differences;
- both adjacent Richardson defect differences;
- the maximum Flint-radius Frobenius norm;
- `100*eps_binary64*max(1,||T12||_2)^2`;
- `1e-135`.

Require the operational `K12` defect norm to be below ten times this
calibrated uncertainty.  Also require `abs(log(abs(det(T12))))` to lie below
ten times its identically constructed cross-level and conditioning
uncertainty.  These are formula controls, not physical results.

## 7. Exact known control and hostile failures

Before using a Regge block, test the same implementation on the exact scalar
quadratic action whose Hessian entries are

```text
K_xx=2, K_xo=K_ox=3, K_xn=K_nx=5,
K_oo=7, K_on=K_no=11, K_nn=13.
```

The frozen exact result is

```text
J = [[2,5],[-3,-11]], det(J)=-7,
Y = [[2/7,5/7],[-5/7,-2/7]],
T = [[-5/7,-2/7],[22/7,-1/7]],
T^T Omega T = Omega.
```

Omitting the direct `K_NO` term gives

```text
T_bad=[[-5/7,-2/7],[-55/7,-1/7]]
```

and a nonzero symplectic defect with off-diagonal magnitude `22/7`; this
hostile control must fail the canonicality gate.

Further require:

1. reversing the pre-momentum sign changes the scalar tangent and fails the
   exact expected matrix;
2. a synthetic `+1e-3` change to the first real coordinate of one copy of a
   schedule tangent is classified `SCHEDULE_DEPENDENT` under the actual
   frozen uncertainty;
3. replacing the physical old/final edge map by a cyclic orbit shift changes
   the actual tangent above the measured uncertainty; it is a corruption
   diagnostic and cannot replace the physical map.

If an actual uncertainty is so large that a hostile control is not detected,
the outcome is `CONTROL_FAILED`; no threshold may be changed after seeing
the result.

## 8. Direct schedule comparison

The common physical boundary labels and literally equal minimal bases make
a direct block comparison meaningful.  For every sector use the `K12`
midpoints after physical relabelling and set

```text
N_T=max(1,||T_even||_F,||T_odd||_F),
d_T=||T_even-T_odd||_F/N_T.
```

The uncertainty `e_T` is the sum, for both parities, of the two adjacent
Richardson tangent variations, maximum Flint-radius Frobenius norm and
`100*eps_binary64*max(1,||T||_2)`, divided by `N_T`, plus `1e-135`.

Assign without any post-result matching:

```text
SCHEDULE_ROBUST     if d_T <= 10*e_T,
SCHEDULE_DEPENDENT  if d_T > 100*e_T,
SCHEDULE_OPEN       otherwise.
```

The global label is `DEPENDENT` if any sector is dependent, otherwise `OPEN`
if any is open, otherwise `ROBUST`.  Direct robustness is stronger than
matching eigenvalue multisets.  No eigenvalues are computed in this mission.

## 9. Frozen outcome hierarchy

Assign exactly one outcome:

1. `FINITE_HEIGHT_FULL_BOUNDARY_TANGENT_CONTROL_FAILED` for any provenance,
   registry, background, carrier, branch, basis, reciprocity, scalar-control
   or hostile-control failure;
2. `FINITE_HEIGHT_FULL_BOUNDARY_TANGENT_RANK_OPEN` if any actual sector is
   not certified regular;
3. `FINITE_HEIGHT_FULL_BOUNDARY_TANGENT_CANONICALITY_FAILED` if every sector
   is regular but any canonical control fails;
4. `FINITE_HEIGHT_FULL_BOUNDARY_TANGENT_SCHEDULE_DEPENDENT_PRIMARY` if every
   canonical map passes but the direct schedule label is dependent;
5. `FINITE_HEIGHT_FULL_BOUNDARY_TANGENT_SCHEDULE_OPEN` for the comparison
   gap;
6. `FINITE_HEIGHT_FULL_BOUNDARY_TANGENT_SCHEDULE_ROBUST_PRIMARY` if all seven
   direct comparisons are robust.

A passing verifier means the frozen hierarchy was executed honestly; a
rank-open, canonicality-failed or schedule-dependent result is a valid
negative.  Any material `PRIMARY` outcome requires a separately
preregistered mechanically different replication before consolidation.

## 10. Claim boundary

Even outcome 6 proves only that this finite-height slab has a locally unique,
action-generated, schedule-robust differential for arbitrary boundary
configuration and pre-momentum perturbations with matter frozen.

It does not identify physical gauge-invariant modes, count two graviton
polarizations, prove stability, define an invariant one-step spectrum,
construct a second anisotropic tick, show convergence under refinement,
derive a wave equation, a limiting speed, `G`, Planck units or particle
masses.  External novelty remains **OPEN**.

Only the registered mission-specific verifier is run.  The full suite is
forbidden unless the user later asks for it explicitly.
