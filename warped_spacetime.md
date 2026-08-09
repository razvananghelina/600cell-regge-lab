# Warped inflation tower as a spacetime direction

Date: 2026-07-24

## Audit provenance correction

The current executable fixes its criteria before constructing spectra, but
there is no version-control or external timestamp record proving
pre-result registration of the numerical thresholds, warp roster, or
boundary choices.  “Frozen” therefore means immutable during a rerun, not a
verified preregistration.  The scoped negative remains reproducible for the
displayed finite constructions.

## Verdict

**DERIVED negative under the frozen finite spectral protocol.** The physically
warped operator was tested, rather than the earlier direct product. Neither
the golden warp nor the McKay warp produces a qualifying four-dimensional
plateau in either estimator at any of `N=8,16,24`. The result is unchanged
when both finite-end boundary conditions are changed from Dirichlet to
Neumann. The registered verdict is

`NO_N_STABLE_4D_PLATEAU_WARPED_NEGATIVE`.

This closes the proposed finite 600-cell tower-spacetime route for the stated
Jacobi operator and frozen plateau test. It is not a theorem excluding every
continuum limit, altered radial operator, or different discretization.

## Construction frozen before evaluation

Let `ell=log(phi)` be the proper level spacing. For each of the 52 verified
distinct spatial levels `lambda`, with its exact full-complex multiplicity,
the warped block is

`T_lambda = L_radial + lambda diag(r^(-2n)), n=0,...,N-1`.

The roster is:

| label | `r` | role |
|---|---:|---|
| w1 golden | `phi` | **DERIVED** quasicrystal-inflation warp; primary |
| w2 McKay | `2` | **DERIVED** McKay/Perron--Frobenius warp |
| w3 no warp | `1` | consistency control |

`L_radial` is the uniform second difference divided by `ell^2`, with
off-diagonal `-1/ell^2`. The primary Dirichlet choice uses diagonal
`2/ell^2` at every retained site (zero ghost values). The Neumann sensitivity
choice uses the path quadratic form, with endpoint diagonal `1/ell^2` and
interior diagonal `2/ell^2`. Both conditions are imposed at both finite
ends. They were frozen before spectral evaluation.

The plateau criteria are imported unchanged from
`verify_holographic_dimension.py`: minimum width `0.50` decade, target
tolerance `0.35`, counting log-RMSE at most `0.08`, counting local-dimension
standard deviation at most `0.35`, heat standard deviation at most `0.35`,
and heat derivative at most `1.00`. The counting scan is accelerated with
prefix sums; the verifier checks that it reproduces the original frozen scan
on the 600-cell spectrum, including the selected endpoints and statistics.
No scientific gate is changed.

## Exact mode decomposition

In a basis diagonalizing `D3^2`, the warped spatial term is diagonal in the
spatial-mode label and the radial operator acts as the identity on spatial
cochains. Therefore the full matrix is the direct sum of the `T_lambda`
blocks. This is an exact block-structure statement, not an approximation.

The verifier independently assembles the full warped matrix for the
30-cochain boundary of the 5-cell at `N=4`. Its directly diagonalized
spectrum and the union of its Jacobi-block spectra agree with maximum sorted
eigenvalue error `0.000e+00`. The production calculation uses only the 52
invariant levels and their multiplicities; it never rediagonalizes the
2640-dimensional 600-cell complex. Every dataset closes at exactly `2640 N`
states.

## Frozen results

Every entry below is the result for the target `d=4`. `NONE/NONE` means that
neither counting nor heat has a qualifying interval.

| boundary | N | w1 golden | w2 McKay | w3 no warp |
|---|---:|---|---|---|
| Dirichlet | 8 | `NONE/NONE` | `NONE/NONE` | `NONE/NONE` |
| Dirichlet | 16 | `NONE/NONE` | `NONE/NONE` | `NONE/NONE` |
| Dirichlet | 24 | `NONE/NONE` | `NONE/NONE` | `NONE/NONE` |
| Neumann | 8 | `NONE/NONE` | `NONE/NONE` | `NONE/NONE` |
| Neumann | 16 | `NONE/NONE` | `NONE/NONE` | `NONE/NONE` |
| Neumann | 24 | `NONE/NONE` | `NONE/NONE` | `NONE/NONE` |

The single-floor control remains three-dimensional in exactly the calibrated
sense available for this finite complex: its counting plateau has width
`0.795` decade, `d_N=3.0688`, log-RMSE `0.0659`, and local standard deviation
`0.2303`; it has no qualifying 4D interval. Full counting and heat curves for
every dataset are emitted with:

```text
PYTHONPATH=/tmp/science-python-deps \
python3 reproducible/verify_warped_spacetime.py --full-curves
```

## Geometry fingerprints

The verifier implements the exact flat diagonal kernel

`K_R4(t)=(4 pi t)^(-2)`

and the exact Plancherel representation of the hyperbolic diagonal kernel

`K_H4(t)=exp[-9 t/(4 R^2)]/(8 pi^2 R^4)`
` times integral_0^infinity q(q^2+1/4)tanh(pi q)exp[-t q^2/R^2] dq`,

with the registered curvature radius `R=ell=log(phi)`. It checks the flat UV
normalization and the curvature-gap bottom `9/(4R^2)` (dimensionless value
`9/4`). At `t=R^2`, the implemented per-volume values are
`K_R4=0.118095978309` and `K_H4=0.0152177336345`.

No 4D window appeared, so shape fitting was not activated. There is therefore
no `R4`, `H4`, Euclidean-AdS4, or holographic-bulk fingerprint to report.
Consequently the AdS-versus-observed-dS cosmology tension is not promoted to
a claim here. If a later construction does produce an H4 fingerprint, that
tension must be stated explicitly.

## Boundary and index checks

**STRUCTURAL boundary statement:** on floor `n`, the spatial block is
`r^(-2n) D3^2`. Multiplication by the conformal factor `r^(2n)` recovers
exactly the same 52 `D3^2` levels and every original multiplicity. Thus the
large-scale end retains the 600-cell spectral data precisely as a conformal
boundary spectrum. Without that conformal rescaling the spatial operator
tends to zero; it is not literally an unscaled copy at infinity.

**INDEX OPEN:** the old vertex-operator index `-4` is not the Fredholm,
Witten, or boundary index of this finite positive `D4^2` construction. The
warped computation supplies no new chiral map whose index equals `-4`.
Accordingly the original `d_ST=4` index remains unexplained and is not forced
into the failed warped route.

## Status synthesis

- **DERIVED:** exact Jacobi mode decomposition, multiplicity closure, all
  spectra and frozen plateau negatives, endpoint-condition insensitivity,
  and conformally rescaled recovery of the full 600-cell spectrum.
- **STRUCTURAL:** interpreting the recovered large-scale spectral data as a
  conformal `S3` boundary.
- **DOWNGRADED / rejected:** a dynamical fourth dimension from either
  registered warp at these truncations; any H4/Euclidean-AdS4 identification.
- **OPEN:** a separately preregistered continuum-limit theorem or a different
  derived radial discretization, and reconciliation of the old index `-4`.

No unregistered numerical constants were searched and no PDF build was
attempted.
