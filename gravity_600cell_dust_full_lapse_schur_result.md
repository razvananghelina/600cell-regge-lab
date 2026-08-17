# Result: the complete 120-dimensional vertex-lapse Schur sector is regular

Date: 2026-08-17  
Targeted verifier: `18/18`.  No full-suite run was requested or performed.  
Result artifact SHA-256:
`4a441ce6b328ffcbb1b673e1c932d411c6a8a00434107bc010e44537190a9349`.

## Provenance

- prior-art gate: `58f14e1`;
- preregistered protocol: `15a0699`;
- registered implementation before evaluation: `ead771f`;
- preserved first `16/18` control-failure artifact: `c494aa7`;
- explicit Lorentzian-reality clarification: `6e64120`;
- corrected implementation before rerun: `7ce3995`;
- passing machine result: `df9e4bb`.

The first failure was not deleted.  It incorrectly required raw Lorentzian
boost angles to be real.  The correction applies the already intended reality
gate to the assembled canonical kernel, whose imaginary contamination is
below `4.63e-79`.  No Schur entry, determinant, rank threshold, or subspace
threshold changed.

## Complete object

At the accepted non-static dust tick, for each of the two schedule parities,
the full canonical pre-Legendre Jacobian is

```text
J : (840 internal + 720 new-boundary variations)
      -> (840 internal equations + 720 old momenta).
```

The 120 timelike pole edges are exactly five free regular `2T` orbits.  They
were selected by edge geometry, not by singular vectors.  In each irreducible
sector of dimension `d`, the precision calculation partitioned

```text
J_d = [ A_d  B_d ],      A_d : 60d x 60d,
      [ C_d  D_d ]       S_d : 5d x 5d,

S_d = D_d-C_d A_d^(-1)B_d.
```

All local geometry and angle derivatives were rebuilt at 100 decimal digits.
All `A_d` solves and Schur matrices were then evaluated with complex Flint
balls at 80 decimal digits.  No binary64 inverse enters the Schur result.

## Controls

**DERIVED COMPUTATIONAL.**

- The high-precision `2T` basis has orthonormality and invariance residuals
  between `1.0e-100` and `1.6e-98`.
- The 65-row representative kernels reconstruct every stored complete-block
  singular multiset within `1.70e-13` normalized error in the even parity and
  `8.99e-14` in the odd parity.
- All four independent angle-derivative step choices produce 1,500
  representative kernel entries and retain the same Lorentzian branch.
- The 120 independently defined vertex-lapse columns have rank 120, and their
  sum reproduces the frozen collective lapse direction with printed error
  exactly zero.
- Every one of the seven strong blocks is resolved and invertible.  Its weakest
  nonzero margin above the preregistered `100 epsilon` boundary is at least
  `1.26e10`.

## Rank result

**DERIVED COMPUTATIONAL — `FULL_LAPSE_SCHUR_REGULAR`.**

For both schedule parities,

```text
resolved pole-Schur rank = 120 / 120,
error-consistent nullity = 0,
numerically open count   = 0.
```

Every Flint determinant ball for every strong and Schur block, under all four
derivative estimates, excludes zero.  The weakest Schur singular value clears
the calibrated nonzero boundary by at least `4.08e12` in the even parity and
`5.97e12` in the odd parity.  The old marginal factor `1.092` was therefore a
binary64 full-SVD limitation, not evidence for a true null direction.

The common Schur scale is

```text
|s_lapse| = 4.24456181727093e-9.
```

It is tiny compared with the next full-J singular value near `43`, but it is
not numerically compatible with zero.

## Stronger scalar pattern

**PATTERN, observed after the rank calculation.**  Every `5d x 5d` Schur
midpoint is numerically

```text
S_d = -4.244561817270933e-9 I_(5d).
```

Across all seven irreps and both parities, the fitted scalar spread is at most
`8.28e-25`.  Within a block, the non-scalar spectral norm is at most about
`8.28e-25` in binary64 inspection; the high-precision derivative uncertainty
is of order `1e-23` or smaller.

This is much stronger than `2T` symmetry alone, because the commutant on five
regular copies is not one-dimensional.  It was noticed after seeing the
result, so it is not promoted to an exact identity or a theorem.  It now needs
an independent derivation from the vertex-local incidence and the
Schlaefli/dust terms.

## Are these exactly geometric lapse displacements?

**OPEN under the preregistered thresholds.**

The canonical Schur-lift subspaces and the frozen geometric vertex-lapse
subspaces have projector distances between approximately `1.14e-5` and
`1.30e-5`.  The protocol classified `distance < 1e-8` as identified and
`distance > 1e-4` as separated; therefore all these comparisons remain
`NUMERICALLY_OPEN`.

The canonical lifts and the weakest `5d` full-J right-singular subspaces are
much closer, at binary64 distances `0` to `2.98e-8`, but most again sit on or
above the conservative identification boundary.  This supports the Schur
localization but does not establish an exact generator identity.

## Prior-art reconciliation

The result is structurally consistent with curvature lifting exact Regge
vertex-displacement constraints into background-dependent
[pseudo-constraints](https://arxiv.org/abs/0905.1670), and with the canonical
pre/post-Legendre framework of [Dittrich and
Hoehn](https://arxiv.org/abs/0912.1817).  The post-result search found no
primary source reporting this exact scalar-like 120-dimensional Schur sector
on a dust 600-cell.  Absence from this search does not prove novelty;
external novelty remains **OPEN**.

## Physical meaning and next gate

**DERIVED:** the complete one-tick canonical relation is locally invertible
also in all 120 lapse-like directions.  The 120 weak modes are not exact gauge
nulls at this curved dust background.

**STRUCTURAL / plausible:** their tiny common stiffness is the signature of
discretization-lifted pseudo-constraints.

**NOT DERIVED:** a physical dust clock, a gauge quotient, two graviton
polarizations, stable propagation, a continuum dispersion relation, a
limiting speed, an absolute tick, or Planck units.

The precision objection to the full tangent map is now removed.  Before a
large inversion, the most economical next work is to derive or refute the
post-observed scalar identity.  If it is structural, the 120 dangerous
directions can be eliminated analytically; then the complete `1440 x 1440`
boundary tangent can be built sector by sector and tested for symplecticity,
stability, and only afterwards spatial dispersion.
