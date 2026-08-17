# Direct Regge acceleration under projected 600-cell refinement

Date: 2026-08-17

Prior-art commit: `81db1ec`  
Frozen protocol: `23157e2`  
Stage-A coefficient commit: `98a1e1d`  
Stage-B verifier-registration commit: `33ee717`

Frozen Stage-A input:

```text
reproducible/gravity_600cell_projected_refinement_acceleration_blind.json
SHA-256 640bc0dd3d6f1ae727f8113bf29514878874effffd14f539f5a43e3c3b18d069
```

Stage-B comparison artifact:

```text
reproducible/gravity_600cell_projected_refinement_acceleration_comparison.json
SHA-256 eccefbc6d386abcaaf968c5d7d2371e4d8a20746d86ce04251dd02584cfa8513
```

Targeted verifier:

```text
reproducible/verify_gravity_600cell_projected_refinement_acceleration_comparison.py
```

It passes `5/5` and returns

```text
PROJECTED_REGGE_ACCELERATION_ROBUST_TWO_LEVEL_IMPROVEMENT.
```

No full-suite result is claimed here.  In accordance with the current run
policy, only the mission-specific verifiers were run.

## Exact comparison

The coefficient `a` is defined by

```text
log(R_1/R_0)=a*(tau/R_0)^2+O((tau/R_0)^4).
```

The time-symmetric closed-dust Friedmann half-step gives `a_FLRW=-1/2`.
The table reports `a/a_FLRW`, its absolute distance from one, and the two
successive error-reduction factors.

```text
variant                    K0 ratio/error       K1 ratio/error       K2 ratio/error       reductions
legacy_float_shortest      1.07897949/.07897949 1.01913429/.01913429 1.00474446/.00474446 4.12764, 4.03297
first_tie_rank_0           1.07897949/.07897949 1.01913442/.01913442 1.00474389/.00474389 4.12761, 4.03348
first_tie_rank_1           1.07897949/.07897949 1.01913440/.01913440 1.00474391/.00474391 4.12762, 4.03347
first_tie_rank_2           1.07897949/.07897949 1.01913434/.01913434 1.00474402/.00474402 4.12763, 4.03336
```

**DERIVED NUMERICAL:** every preregistered refinement tower strictly improves
at both levels.  The coarse `7.8979%` acceleration excess becomes about
`1.9134%` after one refinement and `0.4744%` after two.

**PATTERN:** errors decrease by factors near four.  This is the finite
two-level signature expected of second-order convergence.  Two refinements
do not prove an asymptotic order or an infinite-mesh limit.

## Regulator ambiguity

The exact regular parent leaves three equal central-octahedron diagonals.
Four rules fixed before target comparison resolve this ambiguity differently.
Their coefficient spreads and their ratios to the smallest target distance
are

```text
level 1: spread 6.478125e-8, spread/min-distance 6.771221e-6
level 2: spread 2.828360e-7, spread/min-distance 1.192421e-4
```

**DERIVED NUMERICAL:** this disclosed regulator ambiguity is far smaller than
the refinement error and cannot explain the improvement.  The four towers
are not an exhaustive census of the `3^600` local diagonal assignments, so
full regulator independence remains **OPEN**.

## Prior-art reconciliation after disclosure

The broad outcome is not novel.  Regge solutions are expected, and in known
tests observed, to be second-order accurate; Brewin and Gentle discuss this
distinction between equation residuals and solution convergence in
[*On the convergence of Regge calculus to general relativity*](https://arxiv.org/abs/gr-qc/0006017).
Tsuda and Fujiwara already showed convergence toward closed FLRW for a
pseudo-regular, angularly averaged geodesic 4-dome sequence built from
projected 600-cell subdivisions in
[*Oscillating 4-Polytopal Universe in Regge Calculus*](https://doi.org/10.1093/ptep/ptab079).

**KNOWN:** convergence toward Friedmann under a suitable Regge refinement is
part of the established programme.

**DERIVED NUMERICAL:** the narrower statement established here is that the
complete, direct, non-averaged irregular cellular action on these two
projected-red carriers follows the same improvement pattern for the
homogeneous conserved-dust acceleration coefficient.

**OPEN:** a primary search found no published calculation identical to this
direct coefficient audit.  Failure to locate one is not evidence of novelty.

## Physical verdict and boundary

This result materially strengthens one interpretation of the earlier exact
fixed-600-cell tick law: its `7.90%` excess is a finite-resolution Regge error
that rapidly moves toward the closed-dust Friedmann value, rather than a new
600-cell cosmological acceleration law.

It does **not** establish a new theory of gravity.  Only one homogeneous scale
degree of freedom was varied, with one global lapse and one conserved total
dust mass selected separately on every carrier.  Symmetry has removed local
gravitational degrees of freedom.

The following remain **OPEN**:

- an infinite-refinement convergence theorem;
- local lapse and edge equations on the irregular carriers;
- a locally conserved refined dust discretization;
- anisotropic stability and propagating tensor modes;
- a generalized fluctuation equation `K h = omega^2 M h` with gauge removed;
- a limiting speed, a dynamically selected tick, Planck units, and particle
  masses.

The next discriminating calculation is therefore not another homogeneous
tick.  It is the quadratic Regge action for inhomogeneous boundary-edge
perturbations, followed by the physical generalized spectrum after removing
constraint and gauge directions.  Only a stable low-mode dispersion relation
under refinement could support an effective propagation speed.
