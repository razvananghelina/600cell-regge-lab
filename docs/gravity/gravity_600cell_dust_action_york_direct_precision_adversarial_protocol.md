# Adversarial protocol: independent pseudo-longitudinal residual

Date: 2026-08-19

The primary target-disclosed calculation has already returned
`DIRECT_LONGITUDINAL_IDENTITY_NUMERICALLY_REFUTED`.  This protocol is frozen
after that result and exists to attack it mechanically.  It may not call or
AST-extract any function from the primary direct-precision verifier.

## Frozen inputs

| input | SHA-256 |
|---|---|
| `verify_gravity_600cell_dust_action_york_direct_precision.py` | `73d852d58b21a9a15306a565d5cf4fb998b159fadb82830739ab0996ac07270e` |
| `gravity_600cell_dust_action_york_direct_precision.json` | `d57351e852ab40eb7809397c84e5f57ff58e5ae0bd31f9dcaf87efdc84be76b5` |
| `gravity_600cell_dust_action_york_direct_precision_protocol.md` | `e655d6025e790ff2beb653c5e9f4c2f38233606c3607f9219ae222bafdfed36e` |
| `gravity_600cell_dust_centered_jacobi.json` | `fe0c2d231c2b7eaa8a96cc051de8b3a9b034e384589ab6411db81562af0d9b56` |
| `gravity_600cell_dust_centered_jacobi.npz` | `1077fb562abd4b16a9b5d664d5b7669e2ace0344022aa12bc071fcc4fd4691ef` |
| `verify_gravity_600cell_dust_rigidity_york.py` | `deba8d9f9bca4a5848134943ec77544e5487d44a59c44234f632b6f2aeb51382` |
| `verify_gravity_600cell_dust_conformal_supermetric.py` | `d77dc8853826d9aecc4395fc4aae405d0505bbd644ec3a3229f640b2e980bcb4` |
| `verify_gravity_600cell_dust_full_boundary_tangent.py` | `c8662bb0835865aac6696fc3f474ed668fed3fe393b9c32a59e709a984c35571` |
| `verify_gravity_global_regge_orbits.py` | `ad93cdd08fabeeee56b009f23936696837c4362f88ae23f92a36d0395e61ffaf` |
| `commons/cell600.py` | `ea5bce4b6c52e0834539ca4b1df9c6a67a3a5ed4da32f4e0298a493fc5315c7f` |

Require the primary artifact to pass `17/17`, carry all sixteen
`NONZERO_RESOLVED` labels for both direct invariance residuals and have the
stated refutation outcome.

## Independent construction

1. Load `M,V` only from the previously committed centered archive.  Do not
   reconstruct a local Regge Hessian, an implicit Legendre inverse or a
   principal-function block.
2. Rebuild the 600-cell using the repository's rounded binary carrier and
   normalize every vertex.  Construct `R`, the vertex-tangent projector and
   `D=R P_tan` using the independently audited rigidity helper.
3. Reconstruct the two schedule edge orders and seven symmetry sectors using
   the conformal/rigidity pipeline, not the primary exact-golden-ratio code.
4. In sectors `4,5`, build the action shape space as before, but obtain the
   longitudinal basis by column-pivoted QR of `P_S D`, not by the primary SVD.

The dimensions must again be `5+25` and `15+10`, with the global controls
`rank R=470`, `rank D=354` and `dim(im C intersection im D)=4`.

## Independent residuals

Do not reuse the primary formulas `L*AT` or
`A L-B L(L*BL)^-1(L*AL)` as the deciding observables.

With `L` the QR-selected orthonormal longitudinal basis, form an orthonormal
basis `Q_BL` of `im(BL)` by another pivoted QR and compute

```text
r_span = ||(I-Q_BL Q_BL*) A L||_2.
```

Independently form the generalized operator `G=B^-1 A` and compute

```text
r_comm = ||(I-L L*) G L||_2.
```

Exact generalized invariance requires both quantities to vanish.  Also form
the ordered singular values of `[B L, A L]`; record how many exceed the
rank-15 threshold, but do not use a desired rank to select `L`.

For each schedule/sector family define the empirical error of a scalar as the
maximum displacement across the four frozen derivative variants plus

```text
1000 eps_binary * 30 * kappa
* max(1, ||A||_2, ||B||_2, ||D||_2),
```

where `kappa` is the largest condition number of the QR diagonal block,
`B`, and `B L`.  Apply the frozen `10/100` bands.  Record the ratio between
the adversarial primary residuals and their counterparts inferred from the
primary artifact; agreement is a control, not an outcome criterion.

## Outcome hierarchy

1. `ADVERSARIAL_PSEUDOLONGITUDINAL_CONTROL_FAILED` for any provenance,
   upstream, global-rank, carrier, archive or conditioning failure.
2. `ADVERSARIAL_DIRECT_REFUTATION_CONFIRMED` if both new residuals are
   `NONZERO_RESOLVED` in all sixteen cells and every augmented span has rank
   greater than fifteen.
3. `ADVERSARIAL_DIRECT_REFUTATION_REFUTED` if either new residual is
   `ZERO_CONSISTENT` in any well-conditioned cell.
4. `ADVERSARIAL_DIRECT_REFUTATION_OPEN` otherwise.

Confirmation remains DERIVED COMPUTATIONAL / STRUCTURAL, not a symbolic or
interval theorem.  It supports pseudo-longitudinal rather than exact-gauge
language only for this fixed curved slab.  No full suite is run.
