# Result: corrected strut carrier is not the frozen canonical/dynamic subspace

Date: 2026-08-19

## 1. Status ledger

| stage | commit | result |
|---|---:|---|
| prior-art gate | `f419238` | exact object and literature boundary frozen |
| target-blind carrier protocol | `5e46f63` | no Hessian/Schur/tangent target loaded |
| protocol correction | `67f4e0d` | removed the unjustified uniform-singular-vector requirement |
| carrier verifier registration | `fe5116d` | registered before execution |
| first symbolic-control failure | `2ceb0df` | preserved; no scientific artifact emitted |
| predicate-only correction | `10f54b3` | structural SymPy equality replaced by simplified difference |
| target-blind carrier artifact | `dab941b` | 13/13, `CORRECTED_STRUT_CARRIER_FROZEN` |
| alignment protocol | `d019ba6` | 42 comparisons and thresholds preregistered |
| alignment verifier | `d765a9f` | registered before target comparison |
| upstream-exit failure | `45dec07` | first wrapper execution preserved |
| upstream-exit audit fix | `d092f5d` | only intercepted a successful upstream `sys.exit(0)` |
| primary alignment artifact | `7ef7a7b` | 15/15, all 42 comparisons `SEPARATED` |
| adversarial protocol | `19ea7d3` | post-result target disclosed |
| adversarial verifier | `2c72165` | polar/projector route registered before execution |
| adversarial artifact | `f2a35f6` | 14/14, separation corroborated |

Frozen artifact hashes:

```text
target-blind carrier
e8035fb9c35ad693d1dd2adbda79485b6dd8d42bdf40a95b70a92466e47027d7

primary alignment
5652b1371563ff11919be130af15f5b48850e2cc65a50ec35e5de85fdb587f90

adversarial alignment
3b0fd6da76195279f1beac540c326c61eff5e3172a63bb89baf69502254c5b1f
```

## 2. Target-blind geometric result

For a lateral homothetic Lorentzian trapezoid and `lambda != 1`, direct
differentiation gives

```text
delta d_(u->v)^2 = (-delta s_u + lambda delta s_v)/(lambda-1).
```

In the logarithmic action variables, with

```text
q_diag = lambda L_minus^2-rho,
kappa  = rho/((lambda-1)q_diag),
```

the geometry-selected column data are

```text
delta log q_(u,v+120) = kappa(c_u-lambda c_v),
delta log rho_v       = c_v,
delta log q_new       = 0.
```

**DERIVED.** The resulting `1560 x 120` carrier has literal pole identity,
rank 120, exact schedule-stabilizer equivariance, exact incidence support and
sums to the old collective lapse column.  Both rational local response blocks
reproduce the same formula.  One row with source/target roles reversed is
rejected.

At the accepted tick its singular values give:

```text
full carrier gain       1.00000000011 .. 17.22881603935
uniform gain            1.00000000055
119-complement gain     6.64539855977 .. 17.22881603935
```

The nonzero uniform/complement coupling is about `4.98e-4`.  Thus the
uniform vector is not a singular vector; the preregistered correction was
necessary.

## 3. Primary comparison

The carrier was committed before any dynamic target was loaded.  Only then
were the following fixed comparisons made in every one of two staircase
parities and seven minimal `2T` sectors:

```text
corrected carrier versus canonical pole-Schur lift,
corrected carrier versus transported largest-modulus branch,
corrected carrier versus transported smallest-modulus branch.
```

This is exactly `2 x 7 x 3 = 42` comparisons and one carrier candidate per
parity; no mixtures or retained best cases were allowed.

All controls passed 15/15 and all 42 labels were `SEPARATED`.  The projector
distances were

```text
canonical  0.997798521 .. 0.998315953
plus       0.997798415 .. 0.998315885
minus      0.997794964 .. 0.998312729
```

corresponding to maximum principal angles of about `86.19` to `86.67`
degrees.  The largest calibrated comparison error was below `6.5e-8`.

The frozen gap gate did not pass globally.  Six sectors have minimum modulus
gaps at least `8.5098`; the last one-dimensional sector has gap only
`1.006134 < 2`.  The preregistered outcome is therefore

```text
CORRECTED_STRUT_EXTREME_SELECTION_OPEN
```

rather than a retroactively stronger refutation label.

## 4. Mechanically different adversarial audit

The adversarial verifier did not call the primary verifier or its decisive
QR/SVD/Schur routines.  It instead used:

1. Hermitian eigendecomposition of each Gram matrix to form the polar basis;
2. the largest absolute eigenvalue of the explicit Hermitian projector
   difference;
3. direct tangent eigenvectors sorted by modulus instead of ordered Schur
   vectors.

It also tested deterministic invertible column mixing, an orthogonal
complement corruption, consistent complex conjugation, both staircase
parities and the frozen source/target role reversal.

All 14 controls passed.  The independent distances remained

```text
0.997794964 .. 0.998315953
```

and differed from the primary values by at most `1.353e-13`.  The maximum
polar orthonormality residual was `7.409e-13`, the maximum direct-eigenvector
residual was `2.477e-19`, and the deliberate role reversal changed a sector
projector by `9.812e-2` and a target distance by `6.993e-7`.

Accepted adversarial outcome:

```text
CORRECTED_STRUT_SEPARATION_ADVERSARIALLY_CORROBORATED
```

## 5. Exact claim boundary

- **DERIVED COMPUTATIONAL, adversarially corroborated:** the complete
  corrected carrier is not equal to the canonical pole-Schur lift in any of
  the 14 parity-sector cases.
- **DERIVED COMPUTATIONAL for the frozen candidates:** it is also not equal
  to either fixed-count transported extreme candidate in any case.
- **OPEN:** a globally well-separated physical extreme branch, because the
  homogeneous fifth-pair gap fails the preregistered `>2` gate.
- **OPEN:** the dimension of the intersections.  A projector distance near
  one is the sine of the *largest* principal angle.  It proves non-equality
  and at least one almost-orthogonal direction, not zero intersection, not
  near-transversality of every direction, and not that all 119 relative modes
  differ.
- **OPEN:** gauge, pseudo-constraint, curvature response, propagating tensor
  modes and physical instability.

The `119+1` dimension match selected the hypothesis but did not survive the
operator comparison.  It is not evidence that the 119 relative strut modes
are gravitons or a lapse sector.

## 6. Post-result primary-literature search

The post-result search used the more precise terms `corrected strut carrier`,
`canonical response lift`, `projector separation`, `direct extreme tangent
branch`, `Regge tent move`, `pseudo constraint` and `600-cell evolution`.

- Dittrich and Hoehn,
  [arXiv:0912.1817](https://arxiv.org/abs/0912.1817), distinguish exact
  flat-background vertex-displacement constraints from curved-background
  pseudo-constraints.
- Dittrich and Hoehn,
  [arXiv:1108.1974](https://arxiv.org/abs/1108.1974), derive canonical
  simplicial evolution from Hamilton's principal function and pre/post
  constraints.
- Hoehn,
  [arXiv:1411.5672](https://arxiv.org/abs/1411.5672), separates lattice
  gravitons from lapse/shift data in flat linearized Regge calculus.
- De Felice and Fabri,
  [arXiv:gr-qc/0009093](https://arxiv.org/abs/gr-qc/0009093) and
  [arXiv:gr-qc/0106077](https://arxiv.org/abs/gr-qc/0106077), evolve a
  dust-filled 600-cell and study its causal stopping point, but do not supply
  this carrier or this sectorwise canonical comparison.
- Jercher and Steinhaus,
  [arXiv:2312.11639](https://arxiv.org/abs/2312.11639), show in Lorentzian
  frustum cosmology that heights can be dynamical and that matter and causal
  structure affect the available evolution branches.

No located primary source computes the object tested here.  Search absence is
not a proof of novelty; external novelty remains **OPEN**.

## 7. Consequence and next falsifiable question

This is an honest negative, not a gravity discovery.  Kinematic admissibility
of 120 arbitrary strut data does not make them canonical lapse freedom.  The
action selects a different graph in the full response space.

The immediate target-disclosed calculation is now the kernel of

```text
G_corrected - C_canonical
```

in all 14 sectors.  Because both matrices have the same literal pole identity,
this kernel is exactly the coefficient space whose corrected geometric
variation also satisfies the canonical strong equations.  Unlike the largest
principal angle, its rank decides the intersection dimension.  Only after
that census should the full 240-dimensional scale-plus-strut carrier be
pulled through the action Hessian.

No full-suite run was performed.  The scoped static guard after registration
reported 350 distinct registered verifiers, 352 verifier files including two
documented exclusions, zero duplicates, zero unregistered files and zero
missing files.

