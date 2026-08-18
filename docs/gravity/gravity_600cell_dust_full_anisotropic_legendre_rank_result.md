# Full anisotropic Legendre-rank result for the 600-cell dust tick

Date: 2026-08-17  
Status: target verifier `18/18`; no full-suite run was requested or performed.

## Question and complete hypotheses

At each of the two already-derived schedule parities, take the accepted
homothetic dust solution with the fixed De Felice--Fabri mass normalization,
the complete 2,400-simplex Lorentzian Regge slab, the uniform 120-pole dust
action, and logarithmic signed-squared edge variables.  Do not restrict
perturbations to schedule-invariant edge orbits.  The carrier is

```text
720 old-boundary edges + 840 internal edges + 720 new-boundary edges
= 2280 variables.
```

The canonical pre-Legendre Jacobian has 1,560 rows and columns: 840 internal
equations plus 720 old momenta, differentiated against 840 internal plus 720
new-boundary variables.  The question is whether this complete Jacobian is
regular, has a resolved nullspace, or remains numerically undecidable under
the preregistered precision hierarchy.

This question does **not** assume a continuum graviton, a gauge quotient, two
polarizations, a dispersion relation, a limiting speed, an absolute tick, or
a Planck scale.

## Provenance ledger

- `f266cc2`: prior-art map, before calculation.
- `6b61e70`: target-free protocol, before calculation.
- `ce13a0e`: frozen full-action directional controls, before calculation.
- `ed0cd61`: registered implementation, pushed before its first execution.
- `59b9a69`: preserved first result, `16/18`, with both reciprocity controls
  failing because the code omitted the binary64 assembly-roundoff term that
  the protocol already required.
- `5bd2cd1`: post-failure specification of that omitted term using a standard
  forward-summation envelope; no observed singular value enters the bound.
- `4183bea`: implementation of the fixed envelope, pushed before rerunning.
- `c874ef9`: corrected result artifact, `18/18`.

No continuum target or desired wave speed was parsed by the verifier.

## Controls

**DERIVED COMPUTATIONAL.**  In each parity:

- all 840 individual internal equations have maximum residual
  `4.26e-31`;
- all 1,440 individual boundary momenta reproduce the accepted orbit result
  within `3.80e-53`;
- the maximum within-orbit gradient spread is below `4.13e-96`;
- the independently assembled binary64 gradient differs by `3.47e-18`;
- six preregistered full-action directional differences reproduce `H w` with
  maximum relative error `3.96e-16`;
- restricting the complete Hessian to orbit-constant vectors reproduces the
  previously committed canonical singular spectrum within `6.93e-14` and the
  `60 x 60` tangent map within `1.79e-10` relative error.

The four high-precision local derivative estimates become identical after
binary64 conversion.  The raw Hessian antisymmetry is `6.85e-13`.  The
independently frozen forward-rounding bound is `9.05e-11`, based on at most 73
summands per entry, so reciprocity passes without symmetrizing the matrix.

The schedule stabilizer is the binary tetrahedral group `2T`.  Its seven
isotypic components have regular-representation dimensions

```text
1, 1, 1, 4, 4, 4, 9,
```

and the deterministic minimal canonical blocks have sizes

```text
65, 65, 65, 130, 130, 130, 195.
```

Their spectra repeat with the required irrep multiplicities.  Maximum
off-block leakage is `5.01e-12`; maximum repetition error is `2.55e-11`.

## Result

**DERIVED COMPUTATIONAL — `FULL_CANONICAL_LEGENDRE_REGULAR`.**

For both schedule parities the complete canonical rank is

```text
rank = 1560 / 1560,
error-consistent nullity = 0,
numerically open count = 0.
```

Therefore the fixed discrete action defines a local implicit canonical update
for every anisotropic boundary direction at this background.  This closes the
specific underdetermination question that could not be answered by the old
order-24-invariant quotient.

This is not automatically good news about continuum gravity.  On a curved
Regge background, full rank can mean that discretization has lifted continuum
diffeomorphism gauge directions.  It establishes local invertibility of this
finite canonical map, not the number of physical continuum degrees of
freedom.

## The weak 120-dimensional sector

**DERIVED COMPUTATIONAL.**  Every block contains a sharply separated weak
cluster:

- the cluster lies between `4.2443e-9` and `4.2450e-9`;
- the next singular value is at least `43.0`;
- the resulting condition numbers are about `6.7e11` to `7.0e11`;
- the cluster has 120 full-carrier directions in each parity.

The minimal blocks contain respectively `5d` weak directions for irrep
dimension `d`.  After restoring representation multiplicity, the count is

```text
5 * sum_irreps d^2 = 5 * 24 = 120.
```

**STRUCTURAL / PATTERN.**  Thus the weak carrier has the representation count
of five regular representations of `2T`.  This matches the earlier appearance
of five collective/relative lapse directions, but equality of counts is not an
identification of generators.

**OPEN.**  The weakest accepted direction clears the frozen nonzero boundary
by only a factor `1.092` in the odd parity.  The preregistered verdict is full
rank, but directly inverting a matrix with condition number near `7e11` in
binary64 would be scientifically weak.  A high-precision Schur-complement
calculation confined to this target-independently selected 120-dimensional
cluster is required before constructing the full tangent map.

## Prior-art reconciliation

**KNOWN.**  Bahr and Dittrich explain that curvature in Regge calculus can
break exact discrete gauge symmetries and replace constraints by
[pseudo-constraints](https://arxiv.org/abs/0905.1670).  Dittrich, Freidel and
Speziale analyze Hessian zero modes and their relation to gauge in
[linearized Regge dynamics](https://arxiv.org/abs/0707.4513).  De Felice and
Fabri evolved a highly symmetric dust 600-cell and studied its stopping point
and causal singularity in [the symmetric model](https://arxiv.org/abs/gr-qc/0009093).

**OPEN NOVELTY.**  The pre- and post-result searches found no primary source
that computes this exact complete `2280`-edge dust-slab Jacobian, resolves it
under the binary tetrahedral stabilizer, or reports its 120-dimensional weak
cluster.  Absence from this search is not proof of novelty.

## Physical status and next gate

**DERIVED:** the theory now has a locally invertible, fully anisotropic
one-tick canonical relation on the fixed discrete background.

**NOT DERIVED:** gravitational waves, Einstein's two polarizations, stable
propagation, a continuum limit, Lorentz invariance, `c`, `t_P`, `m_P`, or
particle masses.

The next defensible calculation is:

1. certify the 120-dimensional weak Schur sector at higher precision and test
   explicit vertex-displacement/lapse generators against it;
2. only if that survives, solve the seven `2T` blocks for the complete
   `1440 x 1440` boundary tangent map;
3. test its symplecticity and stability before comparing its anisotropic modes
   with any spatial Laplacian or continuum dispersion law.

The positive result opens that calculation.  It does not yet open the physics
gate.
