# Target-disclosed protocol: direct-precision action-York residual

Date: 2026-08-19

Prior-art gate commit: `1b076f2`.

The approximate target values `4.07e-4` and `2.96e-5` are already known.
This protocol is therefore a falsification/confirmation test with its complete
construction and outcome hierarchy frozen before the direct Hessians are
evaluated.  It is not a blind discovery test.

## 1. Frozen provenance

Require byte-exact hashes:

| input | SHA-256 |
|---|---|
| `gravity_600cell_dust_action_york_direct_precision_prior_art.md` | `6cff580ed92893dada0099ad3330afbadc4d7a4d67b34984747a699e20a3b057` |
| `gravity_600cell_dust_action_york_negative_result.md` | `338b2fca275c8a7af0f866a680512b88b1be5698a028e0588580b47c9e463c87` |
| `gravity_600cell_dust_action_york_negative.json` | `fd0763af779cb02d96f7e1d7a8856b117dd4bf2c9413f01de6246c597743df27` |
| `verify_gravity_600cell_dust_action_york_negative.py` | `370bae86c27e82f9dda4592e8db1774786a1d2c1919ed96e3fceb6e372e6be7b` |
| `gravity_600cell_dust_centered_jacobi.npz` | `1077fb562abd4b16a9b5d664d5b7669e2ace0344022aa12bc071fcc4fd4691ef` |
| `verify_gravity_600cell_dust_full_boundary_tangent.py` | `c8662bb0835865aac6696fc3f474ed668fed3fe393b9c32a59e709a984c35571` |
| `verify_gravity_600cell_dust_two_step_full_tangent.py` | `c1a3fb09146188c1932ab81629ab69817f2a2f19108fdf8d9e89d78b6de8f717` |
| `verify_gravity_global_regge_orbits.py` | `ad93cdd08fabeeee56b009f23936696837c4362f88ae23f92a36d0395e61ffaf` |
| `verify_gravity_600cell_dust_full_anisotropic_legendre_rank.py` | `834b97c85e386def853b6308e65e831c52d62d7cbcc4b23118602120d6c676e5` |
| `gravity_600cell_dust_homothetic_canonical_lapse.json` | `4b1c59c0518eec11b88b140cdecdf558d762c0d70b4826a758f67544e14ac5b9` |
| `gravity_600cell_dust_second_tick_local_correction.json` | `936984bc84a714140ce16917ee559b346b3c0d4a5ba92d8fb723398a120f8e70` |
| `commons/cell600.py` | `ea5bce4b6c52e0834539ca4b1df9c6a67a3a5ed4da32f4e0298a493fc5315c7f` |

The preceding result must remain exactly
`NEGATIVE_SHAPE_LONGITUDINAL_IDENTITY_OPEN`, with `13/15` checks and sixteen
`EQUALITY_CONSISTENT` projector labels.

## 2. Direct local reconstruction

Use 120-decimal `mpmath` arithmetic.  Reconstruct both schedule parities and
both physical slabs from their local 4-simplex data and the frozen state
strings.  Do not load a midpoint or radius of `M`, `V`, `Gamma` or `Omega`
from the centered archive.

Use four symmetric logarithmic angle-derivative steps:

```text
operational primary = 1e-25
operational shadow  = 1e-18
validation primary  = 3e-25
validation shadow   = 3e-18.
```

The two scales differ by seven orders, and the validation scale differs by a
factor three.  Require the Lorentzian branch, reality and derivative
convergence controls on all twenty local patterns for each slab and schedule.

For only the two already disclosed one-dimensional sectors `4,5`, project the
complete `95 x 95` Hessian and solve the `65 x 65` implicit Legendre system
directly in `mpmath`.  No Flint matrix ball and no binary tangent archive may
intervene.  Construct the two `60 x 60` tangent maps and their principal
function blocks, then

```text
K_- = S_1,10,
K_0 = S_1,11 + S_2,00,
K_+ = S_2,01,
M   = (K_- + K_+)/2,
V   = K_- + K_0 + K_+,
A   = -Herm(V),
B   = -Herm(M).
```

Reconstruct all four derivative variants independently.  Only after all
direct matrices exist, open the hash-frozen centered archive for the formula
control.  The direct primary `M,V` matrices must reproduce its midpoints to
operator distance at most `1e-9`; this is only a formula/order control and is
not the residual classifier.

## 3. Exact golden-ratio geometric carrier

Construct the 120 unit quaternions from the exact algebraic coordinate forms

```text
(+-1,0,0,0),
(+-1/2,+-1/2,+-1/2,+-1/2),
even permutations of (0,+-1/2,+-phi/2,+-1/(2phi)).
```

Sort them by the existing rounded keys solely to recover the frozen vertex
labels.  Require exact agreement of the resulting 720-edge graph and both
thirty-orbit edge orders with the slab carrier.  Form the fractional edge
rigidity differential analytically from these coordinates and project it with
the 120-decimal binary-tetrahedral sector vectors.  The geometry construction
must not read a negative eigenvector.

For every schedule, sector and derivative variant reconstruct

```text
K = im C,
S_H = ker(C* Herm(M)),
L_H = im(P_S D),
T_H = ker(L_H* B),
```

and require the already disclosed dimensions `5+25`, then `15+10`, with all
rank gaps at least `1e-6` and `B` positive definite by at least `1e-3`.

## 4. Two independent invariance residuals

Compute for all sixteen cells:

```text
r_cross = ||L* A T||_2.
```

Independently avoid `T` and generalized eigenvectors.  With `L` Euclidean
orthonormal, set

```text
X = (L* B L)^-1 (L* A L),
r_image = ||A L - B L X||_2.
```

The exact generalized-invariance condition implies both residuals vanish.
Record their ratio after the deterministic conditioning factors are removed;
do not require equality of the two norms.

For each scalar diagnostic use the frozen empirical convergence error

```text
epsilon = maximum displacement over the four derivative variants
        + 1000 eps_binary * 30 * kappa_carrier
          * max(1, ||A||_2, ||B||_2, ||D||_2).
```

Here `kappa_carrier` is the maximum condition number of the direct-sum,
longitudinal-rank and restricted-`B` solves actually used in the cell.  This
is an empirical high-precision convergence certificate, not a formal interval
proof.  Classify zero below `10 epsilon`, nonzero above `100 epsilon`, and the
intermediate band OPEN.

Also record the state residuals and junction bounds from the frozen tick
artifacts.  If a state/junction control fails, no identity verdict is allowed.

## 5. Sign and negative control

Require, under the same convergence construction:

```text
max eig(L* A L) < 0,
min eig(T* A T) > 0.
```

Repeat the deterministic one-vector longitudinal/transverse rotation from the
preceding protocol.  Its distance from the negative generalized projector
must be nonzero resolved.  The direct calculation is allowed to use
generalized eigenvectors only for this sign/control section, never to define
`L_H`.

## 6. Outcome hierarchy

Return exactly one:

1. `DIRECT_ACTION_YORK_CONTROL_FAILED` if provenance, geometry, state,
   branch, derivative, archive-reproduction, rank or kinetic controls fail.
2. `DIRECT_LONGITUDINAL_IDENTITY_NUMERICALLY_REFUTED` if all controls pass
   and either `r_cross` or `r_image` is nonzero resolved in every one of the
   sixteen cells.
3. `DIRECT_LONGITUDINAL_IDENTITY_NUMERICALLY_RESOLVED` if all sixteen cells
   have both residuals zero-consistent, the longitudinal/transverse signs are
   resolved, the projectors agree and the rotated controls separate.
4. `DIRECT_LONGITUDINAL_IDENTITY_OPEN` otherwise.

A refutation means that the negative carrier is only pseudo-longitudinal on
this fixed curved slab.  It is a DERIVED COMPUTATIONAL / STRUCTURAL negative,
not a symbolic theorem and not a refutation of continuum diffeomorphism
symmetry.  A positive outcome has the same limited scope.  No continuum,
polarization, speed, mass or particle target is loaded.

Only this verifier, documentation layout and static registry guards may run.
The full suite remains deliberately excluded.
