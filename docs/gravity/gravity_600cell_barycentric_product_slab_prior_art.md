# Prior-art gate: canonical barycentric product slab

Date: 2026-08-17

Status: completed before constructing the product face poset, enumerating its
chains, or testing staircase containment.

## 1. Exact object

Let `K` be the certified simplicial boundary complex of the 600-cell, with

```text
f(K)=(120,720,1200,600),   |K| homeomorphic to S^3.
```

Let `I` be the regular CW interval with two zero-cells and one one-cell.  Form
the regular product CW complex

```text
X=K x I.
```

Its nonempty cells are exactly pairs `(sigma,tau)`, where `sigma` is a
nonempty simplex of `K` and `tau` is one of the three cells of `I`.  Order
them by componentwise closure inclusion.  Define

```text
B=OrderComplex(FacePoset(X)),
```

the barycentric subdivision of `X`: a vertex for each product cell and a
simplex for each strict chain.

This carrier uses only cell incidence.  It uses no five-colouring, vertex
ordering, even/odd schedule, action value, desired spectrum or physical
constant.

## 2. Why it is being tested

The two bare staircase slabs are triangulations selected by vertex orders and
give different nonlinear canonical maps.  The order complex `B` is a
schedule-free triangulation of the unsplit cylinder.  Its boundary vertices
are the cells of `K`, exactly the 2640-element carrier already used by the
repository's oriented cochain/Kähler--Dirac construction.

The load-bearing question is stronger than canonicity:

> Does `B` geometrically refine every staircase triangulation of a local
> tetrahedral prism `Delta^3 x I`?

If yes, it is a direct common refinement on which the two schedule actions can
in principle be compared or coarse-grained.  If no, it remains a canonical
third discretization but does not by itself solve the schedule ambiguity.

## 3. Mathematical prior art

For a finite regular CW complex, the order complex of its nonempty face poset
is its barycentric subdivision and has the same realization.  This is standard
CW-poset topology; see the use of Björner's theorem in
[Clark--Tchernev, *Regular CW-complexes and poset resolutions of monomial
ideals*](https://arxiv.org/abs/1310.2315), and the modern formulation in
[*The discrete flow category: structure and
computation*](https://doi.org/10.1007/s41468-024-00194-5).

Thus `B` triangulates `S^3 x I` and is functorial under every automorphism of
`K` and under reflection of `I`.  This topology/canonicity is **KNOWN**, not a
new physics result.

For the gravity interpretation:

- Bahr and Dittrich, [*Improved and Perfect Actions in Discrete
  Gravity*](https://arxiv.org/abs/0907.4323), define an improved action by
  solving a finer triangulation at fixed coarse boundary data.  Merely choosing
  another triangulation or averaging two actions is not that construction.
- Dittrich and Steinhaus, [*Path integral measure and triangulation
  independence in discrete gravity*](https://arxiv.org/abs/1110.6866), show
  that four-dimensional Regge Hamilton--Jacobi data are generally
  triangulation dependent.
- Bahr, Dittrich and Hellmann, [*Discretization independence implies
  non-locality in 4D discrete quantum gravity*](https://arxiv.org/abs/1404.5288),
  show why a genuinely discretization-independent 4D description is expected
  to acquire nonlocal structure.

No located primary source constructs this exact 600-cell product carrier or
tests it against the repository's staircase triangulations.  External novelty
is **OPEN**.

## 4. Exact local containment test

Realize one product polytope `Delta^3 x I` with its eight vertices in rational
coordinates.  Its product cells are `A x T`, for every nonempty face `A` of
the tetrahedron and each interval cell `T`.  Place the barycentric vertex of a
product cell at the product of its ordinary face barycentres.

Enumerate every maximal face-poset chain, hence every barycentric 4-simplex.
For each of all 24 total orders of the four tetrahedron vertices, construct the
standard four-simplex staircase triangulation

```text
{(v_0,0),...,(v_k,0),(v_k,1),...,(v_3,1)}, k=0,...,3.
```

Using exact rational barycentric coordinates, a barycentric simplex refines
that staircase iff all five of its vertices lie in one staircase simplex.
No floating tolerance is needed.

Because the product barycentric complex is invariant under `S4`, containment
for one vertex order is equivalent to containment for all 24.  The exhaustive
enumeration will nevertheless report every order.

## 5. Existing controls

- **DERIVED:** `K` has the exact 600-cell f-vector and simplicial dimension 3.
- **DERIVED:** its oriented cochain carrier has 2640 cells and Betti numbers of
  `S^3`.
- **DERIVED NEGATIVE:** even and odd staircase schedules agree at first order
  but differ quadratically in nonlinear canonical evolution.
- **DERIVED NEGATIVE:** an anchored endpoint counterterm does not remove the
  configuration difference on any of 16 pure-momentum rays.
- **OPEN:** a schedule-free common refinement and any effective action on it.

## 6. Framing attack

Even a successful common-refinement result supplies only a carrier.  It does
not select edge lengths, a Lorentzian branch, dust transport, an action,
coarse-graining constraints, a measure, a physical lapse or a continuum
limit.  Calling it a perfect action would be false.

Conversely, failure of the refinement property kills only this cheap direct
bridge.  `B` remains a canonical triangulation of the cylinder, and an overlay
or further derived subdivision might still refine both staircase complexes.
Adding arbitrary vertices until containment works would not be canonical and
is forbidden.

## 7. Status before enumeration

- **KNOWN:** topology and functoriality of the product face-poset order
  complex.
- **CONTROL:** exact 600-cell incidence and local rational product geometry.
- **OPEN:** global f-vector and exact local staircase-refinement census.
- **NOT TESTED:** Lorentzian realizability, Regge action, dynamics and
  coarse-graining.
