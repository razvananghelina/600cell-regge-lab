# Result: all 119 strong tangent directions change internal Regge curvature

Date: 2026-08-17

## Provenance

```text
prior-art gate                                      564b3a9
preregistered protocol                              2629f5e
registered implementation                           cb31ee9
missing Flint scalar import correction              14ce3f7
full-spectrum conditional control correction        bd043b6
passing result artifact                             a3d9e7d
```

The two implementation corrections were runtime-only harness fixes.  The
first run stopped before the first Flint solve because `acb` was not imported.
The second stopped before producing the first rank record because a
Schur-versus-direct control was incorrectly requested for the unreduced full
map.  Neither correction changed an operator, derivative, selected subspace,
threshold or outcome rule.

Targeted verifier:
`reproducible/verify_gravity_600cell_dust_internal_curvature_response.py`.

Artifact:

```text
reproducible/gravity_600cell_dust_internal_curvature_response.json
SHA-256 95b6edd8e21ad20a0db97a7c8e7027db7da6547b2b994ad1eb595cf2307f29dc
```

Only this targeted verifier and its direct 43-control geometry import were
run.  The full suite was not run.  Final result:

```text
14/14 PASS
STRONG_TANGENT_CURVATURE_INJECTIVE
```

## What was measured

For each of the two frozen schedule parities, the verifier rebuilt the 6,240
triangle hinges and the incidence-derived Jacobian of the 3,840 **internal**
Regge deficits with respect to all 2,280 logarithmic squared-edge variables.
Boundary exterior angles were deliberately excluded.

The causal row convention

```text
kappa_h = -i epsilon_h  for spacelike hinges,
kappa_h =    epsilon_h  for timelike hinges
```

made the underlying operator real to `1.58e-76`.  This phase is invertible,
so it does not affect kernels or ranks.  Every sparse entry reconstructed
from all 3,840 rows agreed with the independently assembled free-orbit kernel
to `1.58e-76`.

The canonical Legendre response then mapped the 1,440-dimensional boundary
phase tangent into slab-edge variations, and the deficit Jacobian mapped
those variations to internal curvature.  The calculation was block-reduced
only after the full literal incidence check.

## Complete preregistered ledger

There were exactly:

```text
14 full maps          = 2 parities x 7 minimal sectors,
28 strong restrictions = 2 parities x 7 sectors x 2 branches.
```

The labels were:

| object | INJECTIVE | PARTIAL/OPEN | ZERO |
|---|---:|---:|---:|
| full boundary-phase response | 12/14 | 2/14 | 0/14 |
| expanding strong restriction | 14/14 | 0/14 | 0/14 |
| contracting strong restriction | 14/14 | 0/14 | 0/14 |

After restoring regular-representation multiplicities, each strong branch
has rank exactly

```text
119 / 119.
```

The smallest singular values over all sectors and both schedules were:

```text
expanding restriction       0.2870106 ... 0.5021315,
contracting restriction     0.0099664 ... 0.0136858.
```

Their most conservative calibrated uncertainties were at most
`1.94e-9` and `2.33e-9`, respectively.  Even the weakest contracting margin
is about `4.55e6` uncertainty units above zero.  Thus this is not a marginal
rank decision.

## The one-dimensional exception in the full phase space

The 12 nontrivial-sector full maps are injective.  Each trivial-sector full
map has the same calibrated ledger:

```text
rank 59, resolved zero count 1, open count 0, out of 60 columns.
```

Its smallest singular value is about `6.04e-14`, below the fixed zero gate
from an uncertainty of `2.55e-9`.  With multiplicities restored, the complete
boundary-phase map therefore has

```text
rank 1439 / 1440,
calibrated nullity 1.
```

This null direction is homogeneous because it lies in the trivial `2T`
sector.  Its identity is not established by this mission.

There is now an obvious post-result numerical coincidence: the blind tangent
census also left exactly one homogeneous reciprocal pair outside the 119
strong pairs, with moduli approximately

```text
0.9939037271 and 1.0061336654.
```

Whether the one-dimensional curvature kernel is one member, a combination,
or neither member of that two-dimensional tangent plane was **not** compared
in the preregistered test.  The equality of counts is therefore only a
**PATTERN / OPEN** clue and must be tested in a new target-disclosed protocol.

## Scientific verdict

**DERIVED COMPUTATIONAL NEGATIVE:** under the complete hypothesis list in the
protocol, the 119 strong expanding directions and the 119 reciprocal
contracting directions are not curvature-preserving lapse/gauge directions.
The restriction of the internal deficit-angle response is injective on every
one of the 28 preregistered subspaces.

This closes the easiest explanation of the large tangent eigenvalues.  Their
near-alignment with the weak Schur lift is an amplification mechanism, but it
does not place them in the kernel of Regge curvature.

It does **not** prove that the 119 modes are physical gravitational waves:

- the background is curved, where exact discrete diffeomorphism symmetry is
  generically broken;
- a constraint-violating direction can also change deficits;
- deficit response alone does not separate Ricci, Weyl, scalar, vector and
  tensor content;
- the Euclidean row norm used for singular values is structural, not a
  derived physical curvature norm;
- one-step expansion is not a multi-tick Lyapunov exponent.

Thus the physically honest statement is narrower: these modes carry a
resolved, injective internal Regge-curvature response and cannot be discarded
as exact curvature-preserving gauge artifacts.

## Post-result prior-art reconciliation

The technical terms exposed by the calculation do not produce a published
identification of this exact one-dimensional kernel or the 119-mode response.
The closest primary results remain:

- Hoehn's
  [canonical linearized Regge analysis](https://arxiv.org/abs/1411.5672),
  which identifies lattice gravitons as curvature degrees of freedom on flat
  backgrounds;
- Bahr and Dittrich's
  [broken-symmetry analysis](https://arxiv.org/abs/0905.1670), which explains
  why curvature removes exact vertex-displacement gauges and produces
  pseudo-constraints;
- Dittrich and Hoehn's
  [covariant-to-canonical construction](https://arxiv.org/abs/0912.1817),
  which locates the same issue in the discrete Legendre evolution;
- Christiansen's
  [linearization of Regge calculus](https://arxiv.org/abs/1106.4266), which
  relates a three-dimensional Regge quadratic form to the continuum
  `curl^T curl` complex but does not cover this four-dimensional Lorentzian
  dust slab.

External novelty remains **OPEN**.  A literature search is not a novelty
proof.

## Next falsifiable gate

The next calculation should not jump to a continuum dispersion relation.  It
should identify the calibrated one-dimensional homogeneous kernel by comparing
it, in a preregistered and target-disclosed order, with:

1. the two individual near-unit homogeneous tangent directions;
2. their invariant two-plane;
3. the canonical weak homogeneous lift subspace;
4. the independently derived geometric lapse subspace.

If none matches, the kernel is a new algebraic combination and its nonlinear
integrability must be tested.  If one matches, we have isolated the only
curvature-preserving boundary-phase direction before asking whether the other
119 curvature-carrying modes satisfy the linearized constraints.
