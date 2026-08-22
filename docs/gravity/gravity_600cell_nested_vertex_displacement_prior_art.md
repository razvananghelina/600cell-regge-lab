# Prior-art gate: nested tangential vertex-displacement carrier

Date: 2026-08-22

Status: completed before constructing or evaluating a new prolongation
matrix, edge differential or Hessian response.

## Exact question and complete hypotheses

Let

```text
K0 = P(sd K_600),
K1 = P(Esd_2(sd K_600)),
P(x) = x / ||x||.
```

The repository has already certified that `K0` and `K1` are closed,
`H4`-equivariant triangulations of the same declared round three-sphere.  The
vertices of `K1` are exactly the old vertices of `K0` together with one
radially projected midpoint for every edge of `K0`.

The present gate asks only whether this construction selects a canonical
linear prolongation from tangential vertex displacements on `K0` to
tangential vertex displacements on `K1`.  It does not ask whether a Regge
Hessian has a small eigenvalue.

For coarse unit vertices `x_i in R^4`, define the refined background points

```text
y_(i,i) = x_i,
y_(i,j) = P(x_i + x_j),  i != j,
```

for the unordered keys already selected by the edgewise carrier.  For
coarse tangents `u_i` satisfying `x_i . u_i = 0`, the candidate prolongation
is the derivative of this same geometric construction:

```text
U_(i,i) = u_i,
U_(i,j) = (I - y_(i,j) y_(i,j)^T) (u_i + u_j) / ||x_i + x_j||.
```

No temporal staircase, action coefficient, singular-value threshold,
continuum target or desired mode enters this definition.

## Repository facts that must not be recomputed or overstated

- **DERIVED COMPUTATIONAL:** `K0` has f-vector
  `(2640,17040,28800,14400)` and `K1` has f-vector
  `(19680,134880,230400,115200)`.
- **DERIVED COMPUTATIONAL:** the `19680` fine vertices are exactly `2640`
  old vertices plus `17040` projected edge midpoints.
- **DERIVED NEGATIVE:** spatial geometry and time orientation do not select
  one of the 24 temporal staircase schedules.  Their internal temporal
  diagonals are distinct.
- **DERIVED EXACT / STRUCTURAL:** the cotangent pullback from six refined
  homogeneous momentum orbits to one coarse momentum has rank one and an
  inverse fiber of dimension five.  Geometry and symplectic pairing alone
  therefore do not select a reverse momentum lift.
- **DERIVED NEGATIVE:** at fixed radius the accepted coarse and refined bare
  stationary covectors differ by about `1.58096%`; mass normalization repairs
  only the declared homogeneous scalar comparison.
- **DERIVED COMPUTATIONAL:** on the single stationary `K0` product, the
  internal kernel is only the product-duration line for all 24 schedules.
  This is not a two-resolution restoration result.

These facts leave a configuration-space feasibility question open.  They do
not license a phase-space or constraint-convergence claim.

## Known framework from primary literature

- Edgewise subdivision is a standard canonical subdivision of a ranked
  simplex, not a new physical principle: Edelsbrunner and Grayson,
  [*Edgewise Subdivision of a Simplex*](https://doi.org/10.1007/s004540010063).
- Embedding maps between discretizations are necessary data in a
  cylindrically consistent dynamics and should ultimately be selected or
  corrected by the dynamics: Dittrich,
  [arXiv:1205.6127](https://arxiv.org/abs/1205.6127).
- Exact vertex-displacement constraints occur in flat linearized Regge
  calculus, while curved finite Regge backgrounds generically replace them
  by pseudo-constraints: Hoehn,
  [arXiv:1411.5672](https://arxiv.org/abs/1411.5672), and Bahr and Dittrich,
  [arXiv:0905.1670](https://arxiv.org/abs/0905.1670).
- Improved or perfect actions can restore discrete gauge symmetry in special
  settings; a bare subdivision does not do this by itself: Bahr and Dittrich,
  [arXiv:0907.4323](https://arxiv.org/abs/0907.4323).

The literature establishes the framework, not the project-specific
prolongation certificate below.  Search absence is not evidence of external
novelty; novelty remains **OPEN**.

## Framing attack

### What this map can establish

Because every old vertex is literally retained in `K1`, the displayed
vertex map is injective if it is well defined: its old-vertex rows reproduce
every coarse tangent exactly.  It would therefore identify one common
geometric carrier across two spatial resolutions without matching
eigenvectors after seeing their eigenvalues.

The derivative is also forced by the already selected radial midpoint rule.
If it is tangent and `O(4)`-equivariant, no parent-simplex, vertex label,
inner product, temporal schedule or Schur coefficient remains to be chosen.

### What this map cannot establish

1. It transports spatial tangential displacements only.  It does not provide
   the normal/lapse part of a four-dimensional vertex-displacement gauge
   generator.
2. It is a configuration prolongation, not an inverse cotangent lift.  The
   already proved five-parameter momentum ambiguity remains.
3. It does not make the bare actions cylindrically consistent and does not
   prove that their Hessians converge.
4. It does not supply an on-shell `K1` finite-height background with matched
   conserved matter.  That is a separate prerequisite for a dynamical
   comparison.
5. Restricting a future Hessian to this carrier can falsify restoration of
   these tangential shift-like directions.  A softening result would not by
   itself recover the full Hamiltonian and diffeomorphism constraint set.

Accordingly, a positive feasibility result is infrastructure, not gravity.

## Target-free feasibility theorem to test

Before opening any action singular value, certify all of the following by two
mechanically distinct constructions:

1. the fine-key census is exactly `2640` singleton/old keys plus `17040`
   two-endpoint edge keys;
2. every two-endpoint denominator `||x_i+x_j||` is nonzero and the analytic
   derivative agrees with a centered finite-difference derivative of the
   normalized midpoint map;
3. the derivative sends every coarse tangent to a fine tangent;
4. the old-vertex restriction composed with prolongation is the identity,
   hence the prolongation has full column rank on the `3*2640` tangential
   carrier without a numerical rank threshold;
5. ambient `O(4)` covariance holds for independent fixed rotations and
   reflections;
6. the construction is independent of parent chamber and of all 24 temporal
   staircase schedules;
7. the induced fine-edge variation agrees with direct differentiation of
   fine squared chord lengths;
8. a deliberately corrupted midpoint weight and a deliberately omitted old
   vertex are both detected.

If any of items 1--7 fails, the declared carrier is not canonically matchable
and the present refinement route closes.  If all pass, the next admissible
physics calculation is not a spectrum: it is first the construction of a
matched on-shell `K1` background and matter normalization for the same
finite-height action prescription.

## Status before calculation

- **KNOWN:** the subdivision, radial derivative and discrete-gravity
  embedding-map framework.
- **DERIVED:** the two finite carriers and their vertex nesting.
- **OPEN:** the complete mechanically replicated prolongation certificate.
- **OPEN:** matched `K0/K1` on-shell finite-height backgrounds.
- **OPEN:** constraint restoration, physical gravitons, a wave equation,
  `c`, `G`, Planck units and particle masses.
