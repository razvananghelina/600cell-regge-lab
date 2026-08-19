# Prior-art gate: complete variable-face closure on the 600-cell dual graph

Date: 2026-08-19

## Exact object and hypotheses

The carrier is the complete regular 600-cell boundary complex:

```text
120 vertices, 720 edges, 1200 triangular faces, 600 tetrahedra,
connected four-regular tetrahedron-dual graph.
```

Each tetrahedron labels one homothetic Lorentzian frustum and carries the
already derived six-dimensional upper-edge-plus-strut infinitesimal kernel.
Each shared spatial triangle labels one lateral triangular-frustum face.  The
two-cell theorem derives, rather than fits, a one-dimensional infinitesimal
transition line on that face.

The new operator is the complete linearized face-matching matrix with

```text
6 variables per tetrahedral frustum,
1 derived transition coefficient per shared face,
exact matching of the three shared upper vertices.
```

The background transition between adjacent local frames must be the full
four-dimensional Lorentz isometry of the shared **lateral** frustum face.  A
spatial 600-cell face reflection embedded trivially in time is not sufficient
when `lambda != 1`; using it would silently freeze or misidentify the
expanding connection.

The question is the exact global kernel dimension at static and expanding
homothetic representatives.

## Body--hinge prior art

### KNOWN

An infinitesimal motion of a body--hinge framework assigns a screw coordinate
to every rigid body and permits the difference across a hinge to lie on one
hinge line.  A single hinge leaves one relative degree of freedom.

Tay and Whiteley characterized generic body--hinge rigidity by spanning-tree
packings.  In dimension three, where a rigid body has six screw coordinates,
each hinge is equivalent to five independent body--bar constraints; generic
rigidity requires `5G` to contain six edge-disjoint spanning trees.  Primary
sources and a later proof of the molecular version include:

- [Whiteley, *The union of matroids and the rigidity of frameworks*,
  1988](https://doi.org/10.1137/0401025);
- [Tay, *Linking (n-2)-dimensional panels in n-space II*,
  1989](https://doi.org/10.1007/BF01788696);
- [Katoh--Tanigawa, *A Proof of the Molecular Conjecture*,
  2009](https://arxiv.org/abs/0902.0236).

The 600-cell tetrahedron-dual graph is four-regular.  If it is four-edge
connected, every partition into `k` parts has at least `2k` crossing edges;
therefore

```text
5 E_cross >= 10 k >= 6(k-1).
```

The Nash-Williams/Tutte condition then predicts generic body--hinge rigidity.
This is a combinatorial control, not a verdict for the special realization.

### Why the theorem does not decide this mission

The present six-dimensional local space is not the full ten-dimensional
Poincare screw space of a generic four-dimensional rigid body.  It is the
strut-preserving subspace of a Lorentzian frustum.  Its hinge lines are fixed
by regular 600-cell incidence and homothetic geometry, hence highly symmetric
and nongeneric.  Moreover, adjacent local spaces are related by curved
Lorentz transports rather than one global Euclidean frame.

Generic body--hinge rigidity can therefore be used only as a rank control and
a warning that positive modes would be special.  It cannot falsify such
modes.

## Regge/connection prior art

Discrete connection matrices live on codimension-one faces and their ordered
products give hinge curvature in connection formulations of Regge calculus:
[Khatsymovsky 2015](https://arxiv.org/abs/1509.04974).  The connection is
metric-dependent in a discrete Levi-Civita construction:
[Khatsymovsky 2019](https://arxiv.org/abs/1906.11805).  Shape matching and
Levi-Civita conditions constrain the larger connection phase space:
[Dittrich--Ryan 2008](https://arxiv.org/abs/0807.2806) and
[Anza--Speziale 2014](https://arxiv.org/abs/1409.0836).

These sources justify including the transition variation.  They do not give
the complete matrix or its rank on this carrier.

## CONTROL

1. Reconstruct the exact golden-field 600-cell and require
   `(120,720,1200,600)`, two tetrahedra per face and connected dual graph.
2. For every directed adjacent pair, construct the full Lorentzian lateral-
   face transition directly from the six shared bottom/top vertices.  Require
   exact metric preservation, inverse pairing and mapping of all six face
   vertices.
3. Every isolated face block must reproduce the accepted local theorem:

   ```text
   fixed transition compatible dimension    6,
   variable transition compatible dimension 7,
   relative quotient dimension              1.
   ```

4. On a dual spanning tree, the reduced body--hinge matrix must have nullity

   ```text
   6 + (600-1) = 605.
   ```

5. In an untwisted generic-control realization on the same dual graph, the
   appropriate body--hinge matrix must retain exactly the six global trivial
   screw motions when the Tay--Whiteley condition is satisfied.
6. Reverse every face orientation and use an odd relabelling of the canonical
   tetrahedron; global kernel dimensions must be unchanged.

## OPEN

- The exact complete variable-connection kernel dimension.
- Whether static and expanding strata agree.
- Whether any positive kernel decomposes into scale, lapse, anisotropic or
  gauge directions.
- Whether the actual symmetric carrier attains the generic body--hinge rank.
- Finite reconstruction, an action, dynamics and continuum convergence.

## Proposed difference

The proposed result is the exact complete compatibility rank for the derived
six-dimensional frustum bodies and one-dimensional face-transition lines on
the regular 600-cell.  No located source states that specialized rank.

External novelty remains **OPEN**; search absence is not proof.

## Decision boundary

- zero global kernel: variable face connections are nevertheless killed by
  complete closure, rescuing infinitesimal metric rigidity by a stronger and
  correctly formulated theorem;
- positive global kernel: record every dimension before interpretation and
  classify it in a separate target-disclosed test;
- failed transition, tree, generic or convention control: the construction
  remains **OPEN** and no physical Hessian is authorized.

Only preregistered targeted verifiers and static registry guards may be run.
