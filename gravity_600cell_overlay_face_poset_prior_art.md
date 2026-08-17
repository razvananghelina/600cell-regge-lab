# Prior-art gate: complete face poset of the universal staircase overlay

Date: 2026-08-17

Status: written before enumerating any arrangement vertex or proper face.

## 1. Exact object and carrier

Use the certified arrangement inside the tetrahedral prism

```text
P = Delta^3 x I
```

cut by the 14 internal walls

```text
h_A = t-sum_(i in A) lambda_i = 0,
empty != A != {0,1,2,3},
```

together with the six boundary facets

```text
lambda_0=...=lambda_3=0, t=0, t=1.
```

The object to be computed is the complete finite face poset of the resulting
closed polyhedral complex, followed by its order complex (the barycentric
subdivision).  Locally this is a triangulation of the four-ball
`Delta^3 x I`.  Applying the same `S4 x C2`-functorial construction to every
tetrahedron of the 600-cell boundary gives a proposed triangulation of
`S^3 x I`.

The frozen source is the 148-chamber artifact from commit `94c42ad`; no
gravity action, metric or transfer weight is part of this mission.

## 2. Primary prior art

- Hyperplane arrangements and their cell decompositions are standard
  combinatorial objects.  Kastner and Panizzut describe an implemented exact
  cell-decomposition algorithm in
  [*Hyperplane arrangements in polymake*](https://arxiv.org/abs/2003.13548).
  The local construction here is an affine arrangement restricted to a
  bounded prism.
- Zaslavsky-type face enumeration and its extensions are reviewed in
  Randriamaro,
  [*Face Counting for Topological Hyperplane
  Arrangements*](https://arxiv.org/abs/2003.02241).  General face-counting
  theory does not by itself print the face poset of this restricted,
  nongeneric arrangement.
- The barycentric subdivision of a regular cell complex is its order complex:
  vertices are nonempty cells and simplices are strict chains.  Brenti and
  Welker study its face-vector transformation in
  [*f-Vectors of Barycentric Subdivisions*](https://arxiv.org/abs/math/0606356).
- The 14-wall chamber set was recognized post-result as the 148 nonconstant
  positive four-variable threshold functions.  Work on positive threshold
  functions, including their extremal Boolean-lattice points, is represented
  by Lozin et al.,
  [*Specifying a positive threshold function via extremal
  points*](https://arxiv.org/abs/1706.01747).  Chamber counts do not determine
  the lower-dimensional restricted face poset required here.
- A common refined carrier remains distinct from an improved/perfect gravity
  action; see Bahr and Dittrich,
  [*Improved and Perfect Actions in Discrete
  Gravity*](https://arxiv.org/abs/0907.4323).

No primary source located in this gate supplies the requested exact local and
600-cell-global face vectors for this arrangement.  That absence is not proof
of novelty; external novelty remains **OPEN**.

## 3. Why the full face poset is load-bearing

The previous mission proved that coarse incidence cannot select additive fine
weights: `rank(R)=15` and 133 fine directions are invisible.  A defensible
continuation must therefore define a new action directly on a derived fine
carrier.  For ordinary simplicial Regge calculus this requires, before any
metric calculation:

1. all vertices and incidences of the refined carrier;
2. a canonical triangulation of every polyhedral chamber;
3. exact agreement of those triangulations on shared faces;
4. a tractable global simplex count.

The order complex supplies items 2 and 3 functorially once item 1 is exact.
The face census decides item 4 without already committing to a physical
action.

## 4. KNOWN / CONTROL / OPEN

### KNOWN

- The restricted arrangement has exactly 148 full-dimensional chambers.
- It is invariant under `S4 x C2` and restricts compatibly to all four spatial
  prism faces.
- The face poset of a bounded polyhedral complex determines a functorial
  barycentric subdivision.
- The 600-cell boundary has spatial f-vector `(120,720,1200,600)`.

### CONTROL

- Enumerate arrangement vertices from exact rational intersections of four
  independent walls among the 20 internal/boundary hyperplanes.
- Recover each chamber polytope, its facets and all their intersections using
  exact affine ranks; deduplicate globally by vertex sets.
- Require every face interval and every chamber closure to be consistent.
- Check local ball Euler characteristic 1 and boundary Euler characteristic
  0.
- Classify faces by their minimal spatial support and use that classification
  to assemble exact global counts from `(120,720,1200,600)`.
- Count all strict face chains exactly and check that the two time boundaries
  reproduce the already certified barycentric 600-cell boundary vector
  `(2640,17040,28800,14400)`.

### OPEN

- The local arrangement and barycentric-subdivision f-vectors.
- The global `S^3 x I` refined f-vector and Euler characteristic.
- The number of new four-simplices and whether a direct Regge calculation is
  computationally practical.
- A nondegenerate Lorentzian realization of every new four-simplex.
- Regge hinge incidences, deficit angles, dust action and canonical evolution.

## 5. Framing attack

A successful census is still combinatorics, not dynamics.  Barycentric
subdivision is canonical relative to the certified polyhedral overlay, but it
is not guaranteed to preserve desirable Lorentzian simplex shapes or to give
a perfect action.  Conversely, a very large simplex count is a legitimate
practical negative: it would show that the universal-overlay route is formally
clean but computationally ill-conditioned for the next gravity calculation.

