# Adversarial protocol: projective-coordinate audit of prism rigidity

Date: 2026-08-19

Primary protocol commit: `360d2a5`.
Primary verifier commit: `8825c76`.
Primary artifact commit: `f5d7c9b`.

Primary artifact SHA-256:

```text
ce9eb1917dd647c6dd8155a0f9646a72dc7734c0310f763ec31e070403230db8
```

## Claims under attack

The primary exact-minor calculation reported:

1. the six quadrilateral-planarity conditions refute the naive six-flex
   graph count for unequal-scale and generic projective prisms;
2. the equal-scale cell nevertheless has three genuine non-isometric modes;
3. these modes translate the top tetrahedron tangentially and admit finite
   same-length/different-volume realizations;
4. six diagonal lengths restore rigidity but do not select a staircase.

The audit must not reuse the 24 determinant-minor planarity equations.

## Independent projective parameterization

The labelled face lattice of `Delta_3 x I` is preserved by projective
transformations.  Embed every affine vertex as homogeneous `(x,1)` in
projective four-space.  For an infinitesimal projective transformation

```text
Y -> (I + epsilon A)Y,
```

the induced affine velocity is computed directly after dividing by the last
homogeneous coordinate.  Fix one scalar entry of `A` to remove the projective
multiple, leaving 24 independent projective directions.

This gives a 24-parameter tangent chart of planar-faced prism realizations
without writing a single face-planarity minor.  Differentiate the 16 natural
squared lengths with respect to those parameters.  The quotient shape count
is obtained by subtracting the ten Euclidean/Lorentzian isometry directions:

```text
projective tangent dimension 24
minus length-map rank
minus isometry dimension 10.
```

The audit will use both `diag(1,1,1,1)` and `diag(1,1,1,-1)` on the frozen
scale family `q={1,9/10,11/10,2}`.

## Finite audit

Rebuild the two preregistered finite pairs without importing the primary
constraint matrices.  Check:

1. equality of the complete 16-length vectors;
2. inequality of at least one cross-diagonal length;
3. inequality of the oriented four-volume determinant;
4. preservation of parallel top/bottom tetrahedra and parallelogram side
   faces by direct vector identities.

The cross-diagonal comparison is important: it shows concretely which metric
information the natural edges omit.

## Controls

- **positive determination control:** every `q!=1` member must have
  projective length-map rank 14, leaving only ten isometries;
- **negative determination control:** `q=1` must have rank 11, leaving three
  additional modes;
- **singular-stratum control:** the rank must return to 14 on both rational
  sides `q=9/10` and `q=11/10`;
- **signature attack:** all conclusions must agree in Euclidean and
  Lorentzian signature.

## Frozen verdict

If the projective-coordinate ranks disagree with the primary quotient-flex
counts, the result becomes **OPEN**.  If they agree and the finite pairs
survive, the mixed result is adversarially corroborated:

- generic unequal-scale local polytopal geometry is infinitesimally
  determined by the 16 lengths;
- equal-scale geometry is not;
- the three missing data have the kinematics of a discrete shift, but calling
  them gauge remains **OPEN** until a full-carrier constraint analysis.

The audit is registered.  No full suite is run.

