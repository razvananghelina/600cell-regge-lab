# Preregistration: blind full-boundary canonical tangent census

Date: 2026-08-17

Prior-art commit: `5dfc2a7`.

Status: frozen before constructing any full-boundary tangent block, singular
value, eigenvalue or schedule comparison.  The already committed
orbit-quotient tangent is a control input, not a target.

## 1. Frozen inputs and exclusions

Require these exact SHA-256 values:

```text
homothetic dynamic tick:
4b1c59c0518eec11b88b140cdecdf558d762c0d70b4826a758f67544e14ac5b9

two-slab geometric identification:
a5a22d219b71e49c154c1ef80ed9da93b1aef0b93cd2d6ed22f041b71f62db77

committed orbit-quotient tangent:
1ed8d63b4c8a6a4530570a2894820962c7c3c7852747a1112cdf1b242253dbb5

full canonical-rank artifact:
7dc33fcebe8e2cb62be9bba5dfd1fca06fa176a06afe3717d2e9e866f67a7226

full vertex-lapse Schur artifact:
4a441ce6b328ffcbb1b673e1c932d411c6a8a00434107bc010e44537190a9349

audited full-rank source:
834b97c85e386def853b6308e65e831c52d62d7cbcc4b23118602120d6c676e5

audited one-slab geometry source:
ad93cdd08fabeeee56b009f23936696837c4362f88ae23f92a36d0395e61ffaf
```

The accepted tick must pass `7/7`; the orbit tangent must carry one blind
map per parity; the full canonical operator must be resolved rank
`1560/1560`; and the independent pole-Schur audit must be regular rank
`120/120` in both parities.

No continuous harmonic, desired degeneracy, dispersion relation, value of
`c`, Planck scale or particle datum may be parsed.  Only the mission-specific
verifier is run.

## 2. Exact boundary identification

Reconstruct the same order-24 stabilizer and all free edge orbits directly
from the carrier.  For every old orbit type `i` and every group element `g`,
require exactly

```text
shift( old_edge[i,g], +120 ) = final_edge[P(i),g],
```

where `P` is the committed old-to-final orbit map.  Thus all 720 edges are
checked and the group coordinate is fixed, not fitted.  Require the even and
odd stabilizer action sets to be literally equal.

The output relabelling in each irreducible block is therefore

```text
P_d = P_30 tensor I_d,
```

applied to both final coordinates and post-momenta.

## 3. Complete Hessian kernels and the 2T reduction

Use 100 decimal digits and the four already calibrated local angle derivative
steps

```text
operational primary = 1e-20,
operational shadow  = 1e-15,
validation primary  = 3e-20,
validation shadow   = 3e-15.
```

Re-evaluate all twenty local simplex patterns and assemble, for every pair of
the 95 orbit types, the full group-convolution kernel

```text
K_ab(g) = K[(a,e),(b,g)].
```

This includes all `oo`, `ox`, `on`, `xx`, `xn`, `no`, `nx` and `nn`
blocks, not only the pre-Legendre submatrix.  Require the Lorentzian branch,
reality and step controls used by the full-rank audit.

Construct the seven deterministic minimal right-regular bases `W_r`, of
dimensions

```text
d = 1,1,1,2,2,2,3.
```

For each sector form the complete `95d x 95d` block

```text
K_r = (I_95 tensor W_r)^* K (I_95 tensor W_r)
```

from the convolution kernel.  Require basis residuals below `1e-70`, the
dimension identity

```text
sum_r d_r * (60 d_r) = 1440,
```

calibrated Hermitian reciprocity of every projected Hessian, and reproduction
of every committed pre-Legendre minimal singular spectrum within normalized
error `2e-10`.

## 4. Tangent formula in every minimal block

With block coordinates `O=30d`, `X=35d`, `N=30d`, set

```text
J_r = [[ K_XX,  K_XN],
       [-K_OX, -K_ON]],

R_r = [[-K_XO, 0],
       [ K_OO, I]].
```

Convert every 100-digit matrix entry to an 80-decimal complex Flint ball.
For every derivative variant require that `det(J_r)` excludes zero and solve

```text
Y_r = J_r^-1 R_r.
```

Before relabelling, form

```text
delta n      = (Y_r)_N,
delta p_post = [K_NO,0] + K_NX (Y_r)_X + K_NN (Y_r)_N.
```

Finally apply `diag(P_d,P_d)`.  Store the midpoint and the maximum entrywise
ball radii of all four `60d x 60d` maps.

The unique sector with constant overlap one must reproduce the committed
real `60 x 60` orbit-quotient tangent within normalized Frobenius error
`2e-8`.  This is an independent formula/order/sign control.

## 5. Correct complex-sector symplectic control

The chosen minimal bases can be complex.  Their conjugates satisfy

```text
W_r^T conjugate(W_r) = I_d,
```

and the real canonical form pairs a sector with its conjugate copy.  Hence
the correct minimal-block identity is

```text
T_r^* Omega_r T_r = Omega_r,
Omega_r = [[0,I_(30d)],[-I_(30d),0]].
```

It is not `T_r^T Omega_r T_r=Omega_r` for a genuinely complex block.

For every sector and derivative variant compute the midpoint defect.  Define

```text
epsilon_T = ||Top-Top_shadow||2
          + ||Tval-Tval_shadow||2
          + ||Top-Tval||2
          + maximum Flint-radius Frobenius norm + 1e-70,
```

and define `epsilon_sym` by the same three differences among the four
symplectic-defect matrices, plus their maximum Flint-radius Frobenius norm
and `1e-70`.  Require

```text
||Top^* Omega Top - Omega||2 <= 10 epsilon_sym.
```

Also require reciprocal singular-value products to be consistent within ten
times their identically constructed variation proxy, and require
`abs(abs(det(Top))-1)` to satisfy the corresponding calibrated gate.

These are controls implied by the generating-function construction, not
new physical effects.

## 6. Blind spectra and numerical labels

For every operational minimal block record:

- all singular values and all complex eigenvalues;
- spectral radius, minimum eigenvalue modulus and `abs(det T)`;
- eigenvector condition number and maximum eigenpair residual;
- reciprocal singular-value products;
- the optimal matching distance from the spectrum to
  `{1/conjugate(lambda)}`.

Match spectra using the Hungarian algorithm.  Define the empirical
eigenvalue uncertainty by the largest optimal-matching distance from the
operational spectrum to each shadow/validation spectrum, plus

```text
10 eps_machine ||Top||2 max(1,cond(V)).
```

This is a conservative sensitivity label, not a pseudospectral theorem.
Classify each eigenvalue by `abs(abs(lambda)-1)` as

```text
UNIT_CONSISTENT  below 10 epsilon_eig,
RESOLVED_OFF_UNIT above 100 epsilon_eig,
OPEN             otherwise.
```

Because a real multi-step background and a physical norm have not been
constructed, `RESOLVED_OFF_UNIT` means only one-step spectral amplification
in these canonical coordinates.

For singular values use ordered componentwise variation plus the Flint
radius and `10 eps_machine ||Top||2`.  No continuum comparison follows in
this artifact.

## 7. Schedule comparison

The two schedule parities have literally the same stabilizer actions and
therefore the same deterministic seven sector bases.  Compare corresponding
sector spectra by optimal matching, and corresponding ordered singular
spectra componentwise.  Use the sum of the two within-parity uncertainties.

Assign per sector and globally:

```text
SCHEDULE_ROBUST     if both distances <= 10 times uncertainty,
SCHEDULE_DEPENDENT  if either distance > 100 times uncertainty,
SCHEDULE_OPEN       otherwise.
```

No permutation among equal-dimensional irreducibles is allowed after seeing
the spectra.

## 8. Mechanical outcome hierarchy

Assign exactly one combined outcome:

1. `FULL_BOUNDARY_TANGENT_CONTROL_FAILED` if provenance, carrier, branch,
   basis, kernel, boundary or reproduction controls fail;
2. `FULL_BOUNDARY_TANGENT_RANK_OPEN` if a Flint determinant contains zero
   despite the upstream regularity certificates;
3. `FULL_BOUNDARY_TANGENT_CANONICALITY_FAILED` if any calibrated
   symplectic, reciprocal-singular or determinant-modulus gate fails;
4. `FULL_BOUNDARY_TANGENT_BLIND_CENSUS_CERTIFIED` otherwise, with schedule
   robustness and unit-circle counts reported separately.

The verifier passes when it reconstructs the frozen object and assigns this
hierarchy honestly.  A negative or schedule-dependent spectrum is a valid
passing scientific result.

## 9. Claim boundary

Even outcome 4 does not identify gauge-invariant tensor modes, prove a
graviton dispersion law, yield a physical instability, construct a second
tick, select proper time or derive `c`.  Those require, in order, curvature
observables/gauge quotient, a second dynamically solved slab, and refinement.

Only the new targeted verifier and its direct imported geometry controls are
to be run.  The full suite is deliberately excluded.
