# Result: a cellular tetrahedral frustum is missing exactly six local shapes

Date: 2026-08-19

## Headline

The schedule-free cellular tetrahedral-frustum graph does not determine an
anisotropic flat four-dimensional block.  Its 16 squared lengths leave
exactly six non-isometric infinitesimal flexes.  Adding one diagonal on each
of the six quadrilateral faces removes all six flexes, but the existing
geometry selects none of the possible diagonal sets.

The correct classification is:

```text
DERIVED EXACT / STRUCTURAL NEGATIVE:
the present length-only cellular data do not authorize a refined
anisotropic Hessian.
```

The homogeneous cellular action remains valid inside its homothetic ansatz.

## Provenance ledger

| stage | commit | outcome |
|---|---|---|
| prior-art gate | `55277dc` | six-shape prediction disclosed |
| primary exact protocol | `adc7243` | three rational frusta, 24 completions |
| registered primary verifier | `9f80044` | no rank evaluated yet |
| primary artifact | `fa65ac7` | `11/11`, six shapes underdetermined |
| adversarial protocol | `b39cd71` | irregular tetrahedron, polynomial Jacobian |
| registered adversarial verifier | `976e4b7` | no adversarial rank evaluated yet |
| adversarial artifact | `26b7d98` | `9/9`, exact corroboration |

Artifact hashes:

```text
primary
c55f98313121018ff5ca1fc834260e8f2f075248a21fd7b99a356d89b2d18255

adversarial
7763287a12075a911134b24e5f23c3c682198923bda1ab8f75ac1d9541540fc1
```

## Primary exact calculation

For eight vertices in four dimensions, the local configuration space modulo
Lorentz isometries has dimension

```text
8*4 - 10 = 22.
```

The tetrahedral-prism graph contains only

```text
6 lower edges + 6 upper edges + 4 struts = 16 lengths.
```

At three exact rational homothetic representatives, including static and
expanding cases with timelike struts, the Minkowski rigidity matrix gave

```text
cellular rank                         16
full-coordinate kernel               16
independent Lorentz isometries        10
non-isometric flexes                   6.
```

Fixing the lower tetrahedron before calculation gave the same result in a
different form:

```text
top coordinates                       16
top-edge plus strut constraints        10
exact Jacobian rank                    10
remaining flexes                        6.
```

The Euclidean and Lorentz signatures produced identical exact ranks, as
expected from multiplication by a nonsingular metric matrix.

For each of the 24 staircase orders, the six selected cross diagonals acted
with rank six on the flex kernel and raised the complete rank to 22.

## Mechanically independent audit

The adversarial verifier did not import the primary rigidity construction.
It fixed a different, irregular equal-radius rational tetrahedron, wrote the
ten squared-length polynomials explicitly in 16 symbolic top coordinates,
differentiated them with SymPy and evaluated three different rational
representatives.

It returned

```text
base polynomial Jacobian rank/nullity    10/6  in 3/3 cases
staircase determinants nonzero            72/72
cross action on flex kernel rank six       72/72.
```

The preregistered unpredicted census considered all `2^6=64` independent
choices of one diagonal per quadrilateral face.  Every one of all 64 choices
gave full fixed-bottom rank 16 at all three representatives.

Thus local rigidity does not select the 24 staircase completions.  It admits
64 local metric completions; global conforming product topology reduces them
to the already certified 24 total orders, and spatial `H4` symmetry plus time
orientation selects none of those 24.

## Relation to prior art

Symmetry-reduced polytopal cosmology legitimately varies common edge and
strut lengths, but does not supply a general anisotropic uniqueness theorem
for the 16-edge frustum graph.  Non-simplicial blocks may be triangulated,
which supplies the missing lengths, while four-dimensional Regge results are
not generically triangulation independent:

- [Tsuda--Fujiwara, higher-dimensional polytopal universe](https://arxiv.org/abs/2109.01075);
- [Tsuda--Fujiwara, oscillating 4-polytopal universe](https://arxiv.org/abs/2011.04120);
- [Dittrich--Steinhaus, triangulation independence](https://arxiv.org/abs/1110.6866).

First-order Regge calculus permits additional angle/connection variables,
but they are new geometric data with compatibility equations, not values
derived from the insufficient 16 lengths:

- [Barrett, First order Regge calculus](https://arxiv.org/abs/hep-th/9404124).

## What is closed

- **DERIVED EXACT NEGATIVE:** the current 16 cellular graph lengths do not
  determine local anisotropic hinge geometry.
- **DERIVED EXACT:** exactly six additional independent local data suffice at
  the tested representatives.
- **DERIVED NEGATIVE:** neither local rigidity, global conformity, spatial
  `H4` nor time orientation selects one completion.
- **NOT AUTHORIZED:** a refined length-only anisotropic Hessian on one chosen
  schedule.

This reaches the kill boundary for the classical schedule-free,
length-only refinement route.  Computing a large Hessian before resolving
the six shapes would manufacture its answer.

## What remains open

The number six is structurally suggestive: it equals the dimension of the
Lorentz group and the number of independent components of a three-dimensional
extrinsic-curvature tensor.  After the lower and upper intrinsic tetrahedral
metrics are fixed, the unresolved data describe relative placement of the
two tetrahedra.  However:

- **PATTERN:** the dimension matches a relative Lorentz frame or discrete
  connection;
- **OPEN:** a canonical identification of the six flexes with a Lorentz
  holonomy or extrinsic curvature;
- **OPEN:** gluing/closure and first-order equations for such variables;
- **STRUCTURAL unless separately derived:** adding area-angle variables,
  choosing a staircase order or uniformly averaging the 24 actions;
- **OPEN:** every physical spectrum, continuum limit and propagation speed.

The most economical next route is therefore not another length Hessian.  It
is a local first-order test: determine whether a relative `SO(3,1)` frame
maps bijectively onto the six exact flex directions and whether face-gluing
selects or constrains it without fitted coefficients.

