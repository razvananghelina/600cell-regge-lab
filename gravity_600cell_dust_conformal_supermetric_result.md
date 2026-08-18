# Result: the canonical conformal image carries the maximal minority kinetic sign

Date: 2026-08-18

## Provenance

```text
prior-art gate                         96dd1ff
blind protocol                         4d23b25
geometry-only orbit-order correction   298035f
registered implementation              b86014b
first certified result artifact        a202e72
```

The correction is part of the evidential record.  Before any centered
eigenspace was read, the geometry control showed that the two staircase
schedules order the same 720 literal edges differently.  The first protocol
had incorrectly demanded identical row sequences.  The corrected protocol
requires and verifies the exact row permutation and pairs each centered
matrix with its own incidence ordering.  No scientific target was inspected
to make that correction.

Targeted verifier:
`reproducible/verify_gravity_600cell_dust_conformal_supermetric.py`.

Artifact:

```text
reproducible/gravity_600cell_dust_conformal_supermetric.json
SHA-256 b38d55f9f575ddffd34edeaa5e835d9e10919e6d96a0c284d73c31a072675025
```

Two complete targeted executions returned

```text
11/11 PASS
CONFORMAL_MAXIMAL_MINORITY_CERTIFIED
```

The second took `36.95 s` and reproduced the JSON byte for byte with the
same SHA-256 hash.  The full suite was not run.

## Canonical map and exact carrier controls

On the literal regular 600-cell boundary, define

```text
C : R^120 -> R^720,
(C sigma)_uv = sigma_u + sigma_v,
```

where the edge coordinate is the variation of the logarithmic squared edge
magnitude.  This is the tangent of the standard perpendicular-bisector
discrete conformal scaling, up to an irrelevant global factor.

The verifier independently reconstructed and certified:

- 120 vertices, 720 edges and degree 12 at every vertex;
- exactly two unit entries in every edge row of `C`;
- `C^T C = 12 I + A` for the literal adjacency matrix;
- a connected graph with 1,200 triangular faces;
- exact equivariance under all 24 binary-tetrahedral actions;
- five free vertex orbits and thirty free edge orbits;
- an exact row permutation between the two schedule-specific edge orders.

Because `C sigma=0` gives opposite endpoint values on every edge, any odd
cycle forces its values to vanish; connectedness then forces all 120 values
to vanish.  Thus

```text
rank C = 120
```

is an exact graph-theoretic result, not a numerical rank guess.  The binary
SVD control agrees.

In every minimal irrep sector of dimension `d`, the projected conformal
image has exactly dimension `5d`.  All fourteen sector tests resolved this
rank.  The smallest nonzero sector singular value is at least about `2.8795`,
whereas the largest nominal zero is below `1.26e-15`.

## Primary result

For each schedule, sector and derivative variant, let

```text
H = (M+M*)/2,
G = U_C* H U_C,
```

where `U_C` is an orthonormal basis of the canonical conformal image in that
minimal sector.  There are

```text
2 schedules * 7 sectors * 4 variants = 56
```

independent audits of the single geometry-selected map.

Every complete sector matrix reproduced the blind inertia

```text
H : (5d positive, 25d negative, 0 zero, 0 open),
```

and every conformal restriction gave

```text
G : (5d positive, 0 negative, 0 zero, 0 open).
```

Restoring representation multiplicity gives, per schedule,

```text
full edge carrier       120 positive + 600 negative,
conformal vertex image  120 positive,
shape complement          0 positive + 600 negative.
```

The restricted eigenvalues lie between about `45.5355` and `164.861` in the
fixed action convention.  The weakest positive eigenvalue is more than
`2.18e7` times its complete preregistered error envelope.  The restricted
condition numbers range from `1.0000` to `3.3324`.  All 28 schedule
comparisons are `SCHEDULE_ROBUST`; the largest distance is below
`6.97e-7` of its comparison error.

**DERIVED COMPUTATIONAL:** on this fixed carrier and background, the
canonical vertex-conformal image is a maximal subspace carrying the minority
inertia sign of the centered action-derived kinetic bilinear form.  This is
strictly stronger than the previously observed `120:600 = 1:5` count:
geometry now identifies which 120-dimensional family carries the exceptional
sign.

Equivalently, because the restriction is nondegenerate, there is a canonical
bilinear-form orthogonal direct sum

```text
R^720 = im C  direct-sum  ker(C^T H),
```

with `H` positive on the first factor and negative on the 600-dimensional
second factor.  This conclusion follows algebraically from the certified
restriction and full inertia; it does not require choosing Euclidean
eigenvectors.

Changing the sign of the whole action swaps the words positive and negative
but leaves the invariant statement unchanged: the one-per-vertex conformal
family carries the minority sign, while the five-per-vertex shape family
carries the opposite sign.

## Why the stronger Euclidean claim fails without harming the result

All fourteen operational diagnostics returned

```text
SPECTRAL_SEPARATED
INVARIANCE_SEPARATED.
```

The distance from the conformal image to the Euclidean positive spectral
subspace ranges from `0.11489` to `0.28092`, corresponding to maximum
principal angles of about `6.60` to `16.32` degrees.  The leakage

```text
||(I-P_C) H U_C||_2
```

ranges from about `26.28` to `54.89` and is resolved nonzero.

This is a **DERIVED STRUCTURAL NEGATIVE:** the conformal image is neither an
invariant subspace of the matrix `H` nor its Euclidean positive eigenspace.
It does not refute the primary result.  An indefinite bilinear form has many
maximal definite subspaces, while matrix eigenspaces depend on an auxiliary
positive norm and change under non-orthogonal coordinates.  The physically
relevant finite DeWitt-type statement is definiteness of the restricted
bilinear form, not Euclidean spectral equality.

## Prior-art comparison after the result

- [Glickenstein](https://arxiv.org/abs/0906.1560) makes the vertex-conformal
  map **KNOWN** for piecewise-flat two- and three-manifolds.
- [Champion--Glickenstein--Young](https://arxiv.org/abs/1007.0048) study the
  Einstein--Hilbert--Regge functional within discrete conformal classes on a
  double tetrahedron.  That is an important control, but it is not the
  present two-slab centered kinetic bilinear form or the 600-cell.
- [Hartle--Miller--Williams](https://arxiv.org/abs/gr-qc/9609028) show that
  the Lund--Regge supermetric can have triangulation-dependent degeneracy,
  signature change and additional physical timelike directions.  The
  present result is therefore not a generic consequence of saying "Regge".
- [Barrett et al.](https://arxiv.org/abs/gr-qc/9411008),
  [De Felice--Fabri](https://arxiv.org/abs/gr-qc/0009093) and their
  [generalized 600-cell evolution](https://arxiv.org/abs/gr-qc/0106077)
  study 600-cell cosmological evolution, but the located material does not
  compute this canonical conformal restriction.

No located primary source reports the present complete statement.  External
novelty remains **OPEN**, because search failure is not a novelty proof.

## Physical status and limit

- **DERIVED:** the fixed Regge dust action has a canonical one-conformal plus
  five-shape kinetic-sign decomposition per vertex, after symmetry
  restoration.
- **DERIVED:** the result is schedule robust and not produced by fitting,
  target rotation or a numerical tolerance.
- **STRUCTURAL:** this is the finite analogue of the continuum DeWitt local
  `1:5` trace/traceless signature.
- **OPEN:** `H` is not yet proven to be the Lund--Regge supermetric; it is the
  Hermitian part of a centered coefficient obtained after a literal
  identification of adjacent time fibres.
- **OPEN:** the 600 shape directions are not 600 gravitons.  No Hamiltonian
  and diffeomorphism constraint quotient has yet reduced them to two tensor
  polarizations per continuum point.
- **OPEN:** refinement, continuum convergence, physical proper time,
  dispersion, limiting speed and Planck units.

This is genuine progress toward gravitational kinematics, not yet a theory
of observed gravity.  The next load-bearing question is whether the
action-derived constraints or pseudo-constraints select a canonical gauge
and physical tensor quotient inside the 600-dimensional shape factor.
