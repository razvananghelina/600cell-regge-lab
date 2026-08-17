# Preregistration: exact overlay face-poset and order-complex census

Date: 2026-08-17

Prior-art commit: `4a85d25`.

Status: frozen before solving any four-wall intersection.

## 1. Frozen inputs and exclusions

Read only

```text
reproducible/gravity_600cell_universal_staircase_overlay.json
SHA-256 0dd03eed878f599463a44160484c74ddeaa0511fc70c8b2e77bc05a2f36dd3dc
```

and require its certified 148 sign words.  Use the already certified spatial
600-cell f-vector

```text
(120,720,1200,600).
```

No gravity action, metric, edge length, desired face count, continuum target
or fitted coordinate may be read.

## 2. Exact affine representation

Eliminate `lambda_3` by

```text
lambda_3 = 1-lambda_0-lambda_1-lambda_2
```

and work in rational coordinates

```text
(lambda_0,lambda_1,lambda_2,t).
```

Represent the 20 labelled affine forms exactly:

- six nonnegative prism-boundary forms
  `lambda_0,...,lambda_3,t,1-t`;
- 14 internal forms `h_A=t-sum_(i in A)lambda_i`.

Require all 20 labelled hyperplanes to be distinct and require the known
`S4 x C2` action to permute them with the prescribed sign reversal on
internal complement pairs.

## 3. Exact arrangement vertices

Exhaust all `binomial(20,4)=4845` four-hyperplane subsets.  For each subset:

1. solve the affine equations over the rationals if their coefficient matrix
   has rank four;
2. retain the solution iff all six prism boundary forms are nonnegative;
3. deduplicate identical rational points.

For every retained point, record all vanishing labelled forms and all 14
internal signs.  Require the active coefficient vectors to have rank four.
Store every exact coordinate and active set in the artifact.

## 4. Chamber polytopes and complete faces

For each of the 148 frozen chamber sign words, select every arrangement vertex
whose nonzero internal signs agree with the word.  Require its convex hull to
have affine dimension four.

For each of the 20 labelled walls, collect the compatible chamber vertices on
that wall.  It is a chamber facet iff their affine hull has dimension three.
Generate all nonempty intersections of chamber-facet vertex sets.  In a
bounded polytope every face is an intersection of facets and is determined by
its vertices.  Deduplicate the union over all 148 chambers by exact global
vertex sets.

For each face record:

- exact affine dimension;
- sorted global vertex set;
- containing chamber count;
- minimal spatial support, the set of indices `i` for which `lambda_i` is not
  identically zero on the face;
- membership in each of the six prism boundary facets.

Require every face to be contained in a chamber, every chamber-generated face
to occur globally, and inclusion of vertex sets to define a graded poset.

## 5. Convex-cell and topology controls

For every positive-dimensional face `F`, take all strict global subfaces with
vertex sets contained in `F`.  Require their dimension-wise alternating count
to equal the Euler characteristic of the boundary sphere:

```text
sum_(G proper face of F) (-1)^dim(G) = 1+(-1)^(dim(F)-1).
```

Require:

- exactly 148 four-faces;
- local polyhedral-complex Euler characteristic 1;
- local boundary-subcomplex Euler characteristic 0;
- the two time-end subcomplexes are tetrahedral face lattices with
  f-vector `(4,6,4,1)`;
- every codimension-one interior face lies in two four-cells and every
  boundary codimension-one face lies in one.

## 6. Functorial symmetry

Act on rational vertices by all 24 spatial permutations and by `t -> 1-t`.
Require all 48 transformations to preserve the vertex set, the face set,
dimensions, inclusions and boundary membership.  Compute exact face-orbit
counts by dimension, but do not preregister their values.

## 7. Local order complex

The barycentric subdivision has one vertex for every nonempty polyhedral face
and one `k`-simplex for every strict chain of `k+1` faces.  Count all chains by
exact dynamic programming over face inclusion.  Report the complete local
order-complex f-vector and require dimension four and Euler characteristic 1.

Repeat the chain count on the local boundary and on each time end; require
boundary Euler characteristic 0 and the barycentric tetrahedron vector

```text
(15,50,60,24)
```

on each time end.

## 8. Exact global assembly over the 600-cell

Classify each local face and each chain by the spatial support of its maximal
face.  Require `S4` symmetry to make the counts identical for supports of the
same cardinality `s=1,...,4`.  If `a_(s,d)` is the number per one fixed support
of size `s`, assemble

```text
global_f[d] = sum_(s=1)^4 f_K[s-1] * a_(s,d),
f_K=(120,720,1200,600).
```

Apply the same formula to chains to obtain the global barycentric f-vector.
Require both global complexes to have dimension four and Euler characteristic
0, as required for `S^3 x I`.

For each time boundary, require the already certified barycentric 600-cell
boundary vector

```text
(2640,17040,28800,14400)
```

and Euler characteristic 0.

## 9. Mechanical outcome

- Every exact affine, face-poset, convex-cell, topology, symmetry, chain and
  global-assembly control passes:
  `UNIVERSAL_OVERLAY_FACE_POSET_CERTIFIED`.
- Any control fails:
  `UNIVERSAL_OVERLAY_FACE_POSET_CONTROL_FAILED`.

All local/global f-vectors and orbit counts are blind outputs; none is frozen
as a desired result.  A passing outcome is **DERIVED COMBINATORIAL**, not a
Lorentzian or dynamical result.  Run only this registered verifier, never the
full suite.

