# Result: two-frustum gluing removes every relative length-flex mode

Date: 2026-08-19

## Headline

Two homothetic tetrahedral frusta sharing a face do not retain an independent
relative mode inside the six local length-plus-strut flexes.  The compatible
pair space is exactly the six-dimensional diagonal:

```text
(local motion on left) = (local motion on right).
```

This is not caused by an inability to see a face holonomy.  Before the
strut constraints, the full Poincare algebra has exactly the expected
one-dimensional pointwise stabilizer of the shared spacelike triangle.  The
strut-preserving local kernel has zero intersection with that direction.

The correct verdict is:

```text
DERIVED EXACT LOCAL NEGATIVE:
the six cellular flexes do not themselves supply an independent face
connection or holonomy variable.
```

## Provenance ledger

| stage | commit | outcome |
|---|---|---|
| prior-art gate | `4107248` | fixed-frame scope separated from connection theory |
| primary protocol | `e5bf53e` | one-dimensional positive control disclosed |
| primary registration | `ba95728` | no gluing matrix evaluated yet |
| primary artifact | `c3adf97` | `9/9`, diagonal only |
| adversarial protocol | `32adc05` | irregular direct five-vertex union |
| adversarial registration | `f760642` | no union Jacobian evaluated yet |
| adversarial artifact | `708bbf8` | `11/11`, diagonal only corroborated |

Artifact hashes:

```text
primary
0e09c3f8f38c8158deff5b81bc6fe4d5d6dd685a24cce83e015fb95e3f26a70e

adversarial
0f8e70ef89b7fd5a8995349d40c77f6d3f637f2d9ce137ce2c9ff07b2fed2542
```

## Primary exact calculation

The carrier was an exact regular triangular bipyramid: two reflected regular
tetrahedra sharing a face.  For each of

```text
(lambda,tau)=(1,5),(2,5),(3,11),
```

both local six-dimensional Poincare kernels agreed exactly with the direct
ten-length flex kernels and with the analytic static/expanding forms.

A general Poincare Killing field evaluated on the three shared upper
vertices gave

```text
evaluation matrix rank                      9,
pointwise face stabilizer dimension         1.
```

Consequently, two unrestricted Poincare bodies glued on that face have an
eleven-dimensional compatible space: ten common motions plus one relative
face stabilizer.

After restricting both bodies to their accepted six-dimensional
strut-preserving kernels, every representative instead gave

```text
constrained gluing rank                      6,
compatible pair dimension                   6,
relative-difference rank                    0,
relative face stabilizer dimension          0,
compatible space                            exact diagonal.
```

The result is identical on the static and expanding strata.

## Mechanically independent audit

The adversarial verifier used an irregular rational triangular bipyramid and
did not construct two local kernel bases.  It represented the five physical
upper vertices only once, wrote directly

```text
9 distinct union top-edge polynomials
+ 5 strut polynomials
= 14 squared-length constraints
```

in twenty upper coordinates, and differentiated the union polynomial map.
At

```text
(lambda,tau)=(1,7),(2,7),(3,13)
```

it obtained exact rank/nullity

```text
14/6
```

in all three cases.  Each direct six-dimensional kernel equalled exactly the
common Poincare displacement image on all five vertices.

As a mechanically different positive control, the audit used sixteen
unconstrained entries of `A`, four translations and ten separately imposed
Lorentz equations.  Adding pointwise vanishing on the face produced a
rank-19 system in twenty variables, hence the same one-dimensional full-
Poincare stabilizer.  Metric-sign reversal and a different exact rational
boost preserved all decisions.

## Why the face mode disappears

A full Poincare Killing field can fix a spacelike triangle pointwise by a
boost in its two-dimensional normal plane, with an accompanying translation.
But an expanding homothetic frustum requires

```text
b(A)=tau/(lambda-1) A n,
```

while the face stabilizer requires the translation that places the fixed
plane on the actual triangular face.  For nonzero lapse these conditions are
incompatible except at `A=0`.

On the static stratum the local kernel has `A n=0` and spatial `b`.  A
three-dimensional Euclidean Killing field vanishing at three non-collinear
points is also zero.  Thus the intersection vanishes on both strata for
different but exact reasons.

## Relation to prior art

This result is consistent with the distinction in connection formulations:
independent Lorentz/affine transition matrices are additional variables, and
their equations recover Regge geometry only after metricity, torsion or
shape-matching constraints.  Relevant primary formulations include
[Dittrich--Ryan](https://arxiv.org/abs/0807.2806),
[Anza--Speziale](https://arxiv.org/abs/1409.0836),
[Asante et al.](https://arxiv.org/abs/1908.05970), and the affine-connection
Regge formulation of
[Khatsymovsky](https://arxiv.org/abs/1509.04974).

The post-result search found no source identifying the present six
length-flexes with an independent face connection.  Search absence is not a
novelty proof; external novelty remains **OPEN**.

## Important correction to the earlier framing

The diagonal-only result kills the hidden-connection interpretation, but it
also shows that the six local shape freedoms are not independent from cell
to cell.  Therefore a local cell-rigidity failure alone is not yet a proof
that the complete 600-cell slab is globally flexible.

Shared-face consistency propagates one local motion across adjacent cells.
On a connected dual graph, all local motions must be parallel transports of
one seed.  Curvature holonomy around spatial edges may then remove that last
six-dimensional seed entirely.

Thus the earlier statement

```text
local six-shape failure immediately kills every global length-only route
```

was too broad.  What remains killed is a **local schedule-free action whose
hinge data are supposed to come from one cell in isolation**.  A global
closure reconstruction is newly authorized and may either rescue or finally
close the length-only anisotropic geometry.

## Status ledger

| Claim | Status |
|---|---|
| Full Poincare algebra has a face stabilizer | **DERIVED EXACT, dimension 1** |
| That stabilizer survives the strut-preserving six-flex kernel | **REFUTED DERIVED EXACT** |
| Two compatible local six-flexes can differ | **REFUTED DERIVED EXACT** |
| Existing six flexes are an independent face connection | **REFUTED LOCAL, adversarially corroborated** |
| Independent first-order face holonomies are forbidden | **NOT CLAIMED** |
| Global face consistency can reduce local flexibility | **DERIVED EXACT for two cells** |
| Full 600-cell global length geometry is rigid | **OPEN** |
| Curvature holonomy kills the final common seed | **OPEN** |
| Global schedule-free anisotropic action | **OPEN pending closure** |
| Action, dynamics, propagation and continuum limit | **NOT TESTED** |

## Next discriminating step

Build the exact Levi-Civita face transports of the regular 600-cell spatial
triangulation.  Propagate one local six-flex seed around closed dual loops
encircling spatial edges and compute the common fixed subspace of the
resulting affine holonomies.

- fixed subspace dimension zero: global closure removes all local shapes and
  the global length-only geometry is infinitesimally rigid;
- positive fixed dimension: those modes remain genuine global
  underdetermination;
- dependence on a chosen development or schedule: the route remains
  structural and no anisotropic action is authorized.

The known regular-tetrahedron deficit supplies a cheap exact control:

```text
theta = arccos(1/3),
delta = 2 pi - 5 theta,
cos(delta)=241/243,
sin(delta)=22 sqrt(2)/243.
```

These values must be rederived, not assumed, in the global protocol.
