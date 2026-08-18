# Preregistered protocol: direct high-precision shifted stiffness

Date: 2026-08-18

Prior-art gate commit: `52f337b`.

Precision-attribution artifact commit: `67e99fb`.

Status: **TARGET-DISCLOSED, PREREGISTERED BEFORE THE DIRECT RECONSTRUCTION IS
EXECUTED.**

The old rank-`15` result, the shifted midpoint `15+10` pattern and the shifted
OPEN labels are disclosed.  The purpose of this calculation is not a blind
discovery.  It is a controlled attempt to remove the binary serialization
interface selected by the preceding `16/16` attribution audit.

## Complete hypothesis list

The claim being tested is conditional on all of the following:

1. the three accepted fixed-mass homothetic dust slabs and their frozen
   states;
2. the literal 600-cell one-slab Regge action and the two orientation
   schedules;
3. the seven deterministic binary-icosahedral sector bases;
4. the four already frozen complex-step derivative schedules;
5. the action-selected conformal/shape restriction used in the committed
   stiffness census;
6. sectors `4,5`, both parities and all four derivative variants;
7. rigorous 80-decimal Flint ball arithmetic through the canonical tangent,
   principal-function and centered-coefficient reconstruction;
8. one binary64 conversion only after the final centered matrices have been
   formed, with the mandatory half-ULP re-enclosure included in the final
   error.

No statement outside this finite carrier follows from the result.

## Frozen inputs

| input | SHA-256 |
|---|---|
| `gravity_600cell_dust_homothetic_canonical_lapse.json` | `4b1c59c0518eec11b88b140cdecdf558d762c0d70b4826a758f67544e14ac5b9` |
| `gravity_600cell_dust_second_tick_local_correction.json` | `936984bc84a714140ce16917ee559b346b3c0d4a5ba92d8fb723398a120f8e70` |
| `gravity_600cell_dust_third_tick_local_correction.json` | `ebf2f1a11b9a4e9c76fb1ce33066c0782429cf6500770df7bbe4d92de4a050c0` |
| `verify_gravity_global_regge_orbits.py` | `ad93cdd08fabeeee56b009f23936696837c4362f88ae23f92a36d0395e61ffaf` |
| `verify_gravity_600cell_dust_full_anisotropic_legendre_rank.py` | `834b97c85e386def853b6308e65e831c52d62d7cbcc4b23118602120d6c676e5` |
| `verify_gravity_600cell_dust_full_boundary_tangent.py` | `c8662bb0835865aac6696fc3f474ed668fed3fe393b9c32a59e709a984c35571` |
| `verify_gravity_600cell_dust_conformal_supermetric.py` | `d77dc8853826d9aecc4395fc4aae405d0505bbd644ec3a3229f640b2e980bcb4` |
| `commons/cell600.py` | `ea5bce4b6c52e0834539ca4b1df9c6a67a3a5ed4da32f4e0298a493fc5315c7f` |
| `gravity_600cell_dust_shifted_centered.json` | `265bd863de2365f19f7679373155fdaa23fb0bb3e75c221cfd9d9ec5b6ac2a47` |
| `gravity_600cell_dust_shifted_centered.npz` | `c000f4fcae67e6c0648046878c2bd1ffd0616c38510ccf788c67cf99832397b8` |
| `gravity_600cell_dust_shifted_precision_audit.json` | `409e428ca05b4b6c6e380d7af6d84fc3834afc1fdb95065d6f0a05618e1d2cee` |

The old centered archive is a control only.  Its stiffness eigenvalues and
sign labels are not loaded.

## Direct reconstruction

For each parity, construct the second slab at `(a1,a2,r2)` and the third slab
at `(a2,a3,r3)` directly from the accepted states.  Differentiate the local
Regge angle patterns at 100 decimal digits using the four frozen schedules,
assemble the full representative Hessian and project it into sectors `4,5`.

For each slab and cell, form the canonical tangent as a Flint ball directly
from that projected high-precision Hessian.  Without converting that tangent
to binary64, split

```text
T_j = [A_j B_j; C_j D_j]
```

and reconstruct the Hamilton principal-function blocks

```text
S_j,00 = B_j^-1 A_j,
S_j,01 = -B_j^-1,
S_j,11 = D_j B_j^-1,
S_j,10 = C_j - S_j,11 A_j.
```

Require both `B_j` determinant balls to exclude zero and all adjoint/recovery
identities to contain zero entrywise.  Then form, still in Flint arithmetic,

```text
Kminus = S_2,10,
Kzero  = S_2,11 + S_3,00,
Kplus  = S_3,01,
M      = (Kminus + Kplus)/2,
V      = Kminus + Kzero + Kplus.
```

Only `M,V` are converted to binary64 for the final linear-algebra census.
Their component radii must include both their stored Flint radii and half an
ULP of each real and imaginary binary64 midpoint.

## Independent controls

For all `16 = 2*2*4` cells:

- reproduce the exact branch and geometry controls of the source tangent
  construction;
- require the direct `M,V` balls to overlap entrywise with the corresponding
  previously committed broad balls;
- require the new complete `V`-radius Frobenius norm to be smaller than the
  old one by a resolved factor greater than `100`;
- rebuild the rank-`5` conformal image and the rank-`25` action-selected shape
  complement, including its subspace-error bound;
- require the restricted kinetic form `B=-M_S` to be positive-definite
  resolved and the action compatibility residual to be zero-consistent.

The factor `100` is a numerical-method control, not a sign threshold chosen
from a desired eigenvalue.

## Frozen sign and rank gate

On the restricted Hermitian stiffness `A=-V_S`, label every eigenvalue with
the existing thresholds:

```text
POSITIVE_RESOLVED   lambda >  100 epsilon,
NEGATIVE_RESOLVED   lambda < -100 epsilon,
ZERO_CONSISTENT    |lambda| < 10 epsilon,
OPEN                otherwise.
```

Use the first applicable complete outcome:

1. `SHIFTED_DIRECT_PRECISION_CONTROL_FAILED` if any geometric, branch,
   determinant, principal-identity, overlap, carrier, kinetic,
   compatibility or precision-reduction control fails;
2. `SHIFTED_DIRECT_NEGATIVE_RANK_CHANGED` if any cell is completely resolved
   but does not have exactly `15` negative and `10` positive directions;
3. `SHIFTED_DIRECT_NEGATIVE_RANK_OPEN` if no cell refutes the rank but at
   least one cell contains `OPEN` or `ZERO_CONSISTENT` directions;
4. `SHIFTED_DIRECT_NEGATIVE_RANK_PERSISTS` only if all `16/16` cells are
   completely resolved as exactly `15` negative plus `10` positive.

Outcome 4 certifies only persistence of the rank and selected sectors.  It
authorizes a separate projector/connection comparison; it is not itself a
common-fiber, wave, graviton, mass or physical-instability result.

## Prohibited adaptations

- no new derivative step or precision chosen after seeing a sign;
- no deletion of half-ULP or stored ball radii;
- no midpoint-only sign;
- no rotation of fibers by a fitted Procrustes or polar map;
- no use of the old stiffness eigenvalues to compute the new result;
- no full-suite run.

## Disclosed amendment after the first execution

The first execution, recorded at commit `cab30d3`, returned

```text
SHIFTED_DIRECT_PRECISION_CONTROL_FAILED.
```

All `16/16` stiffness cells were individually resolved as `15+10`, but those
signs are ignored because the preregistered principal-identity control failed.
Target-independent diagnostics were added at commit `d2069e0`.  They locate
the failure entirely in the three adjoint identities; all recovery identities
contain zero entrywise.  The adjoint residuals scale as `h^2`:

```text
h = 1e-20 or 3e-20:  Frobenius residual about 1e-33 ... 1e-31,
h = 1e-15 or 3e-15:  Frobenius residual about 1e-23 ... 1e-21.
```

Thus the original wording conflated rigorous Flint propagation *after* the
finite-difference derivative with a rigorous enclosure of the exact action
Hessian.  The latter was not constructed.  Broad binary serialization balls
had previously hidden this finite-step antisymmetric defect.

Before a corrected execution, freeze the following target-independent repair:

1. for every parity/slab/sector, construct all four raw projected action
   Hessians;
2. record `delta_sym = ||H-H*||_F/2` for every raw Hessian and the complete
   family variation `delta_family = max_v ||H_v-H_primary||_F`;
3. require every `delta_sym <= delta_family + 1e-70`; otherwise the control
   fails;
4. use the unique Hermitian projection `(H+H*)/2`, fixed by the exact symmetry
   of an action Hessian, before constructing the tangent;
5. require all principal adjoint and recovery identities to contain zero
   entrywise after that projection;
6. for each parity/sector, add the maximum ordered stiffness-eigenvalue
   variation across all four derivative schedules to every cell's restricted
   error before applying the frozen `10/100` sign thresholds.

This projection has no fitted coefficient and cannot select a desired sign.
The enlarged schedule envelope is more conservative than the first run.

Even if the corrected run passes, the result is **DERIVED COMPUTATIONAL
conditional on the frozen derivative family**, not a formal interval proof of
the exact analytic derivative.  A hyper-dual/automatic or analytic
ball-derivative implementation remains the stronger independent confirmation.
