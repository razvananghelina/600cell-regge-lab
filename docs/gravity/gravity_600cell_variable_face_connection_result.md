# Result: a derived variable face connection restores exactly one relative mode

Date: 2026-08-19

## Headline

The fixed-frame two-frustum theorem is correct but is not the complete metric
gluing problem.  If the face transition is allowed to vary by the unique
Poincare generator fixing the shared lower triangle pointwise, the compatible
space increases exactly from six to seven dimensions:

```text
fixed face transition:      6 common + 0 relative modes;
derived variable transition: 6 common + 1 relative mode.
```

The primary exact verifier returned `11/11`.  A mechanically independent
52-variable polynomial Jacobian on two irregular carriers returned `9/9`.

The accepted verdict is:

```text
DERIVED EXACT, ADVERSARIALLY CORROBORATED:
one connection-coupled relative mode survives across every tested
nondegenerate homothetic-frustum face.

DERIVED CORRECTION:
zero fixed space under the frozen background holonomy does not imply zero
kernel for metric gluing when the face connection is allowed to vary.
```

This correction preserves every exact frozen-connection matrix previously
computed.  It refutes the broader rigidity interpretation, not the arithmetic.

## Complete local theorem

Let `V` be four-dimensional Minkowski space, `n` a unit timelike vector and
`P` a nondegenerate spacelike affine triangle in `n^perp`.  Let

```text
Q = lambda P + tau n,  tau != 0,
```

be its parallel homothetic upper triangle.  Embed `P` and `Q` as the shared
lower and upper triangles of two nondegenerate homothetic tetrahedral frusta.
Restrict each frustum to its exact six-dimensional upper-edge-plus-strut
kernel `K`.

Let `S_minus` and `S_plus` be the Poincare Lie-algebra lines fixing `P` and
`Q` pointwise.  Then

```text
dim(S_minus) = dim(S_plus) = 1,
dim span(S_minus,S_plus) = 2,
dim(K intersect span(S_minus,S_plus)) = 1.
```

Consequently, matching the upper triangles modulo the derived transition
variation in `S_minus` produces a seven-dimensional compatible space.  Its
six-dimensional diagonal is the common local motion, and the quotient is
exactly one-dimensional.  Both the relative local motion and the connection
variation have rank one on that quotient.

### Why the intersection is one-dimensional

Choose the boost generator `B` in the two-dimensional normal plane of `P`.
The pointwise stabilizers of the two parallel planes have the same linear
part `B` and translations differing by `B d`, where `d` translates `P` to
`Q`.

For `lambda != 1`, the local kernel is the graph

```text
b(A) = tau/(lambda-1) A n.
```

Restricted to `A` proportional to `B`, the two stabilizer translations span
exactly the required `b(B)`, giving one intersection line.  At `lambda=1`,
the boost itself is excluded by `A n=0`, but the difference of the two
parallel-plane stabilizers is a pure spatial-normal translation.  That
translation lies in the static kernel, again giving exactly one line.

The same dimension on the static and expanding strata therefore has two
different exact mechanisms; it is not inferred from continuity through the
singular Lorentz chart.

## Provenance ledger

| stage | commit | result |
|---|---|---|
| prior-art gate | `e00791d` | frozen versus variable connection separated |
| primary protocol | `2260b72` | `6 -> 7` prediction frozen before ranks |
| primary registration | `0115377` | exact Poincare intersection verifier frozen |
| primary artifact | `dab8198` | `11/11`, one connection-coupled mode |
| adversarial protocol | `3867f97` | direct 52-variable polynomial system frozen |
| adversarial registration | `e61f32a` | no Poincare chart or analytic kernel imported |
| adversarial artifact | `40a1285` | `9/9`, two irregular carriers corroborate |

Artifact hashes:

```text
primary
001212016553d006862e68edc4f780f37ca1476110b6e0aed3e987f52a43b5e3

adversarial
c8c8c58711e5bf4e49c110e84518ddf643b75cc4377d05fb5f577003b8395466
```

## Primary exact calculation

For the regular reflected pair and

```text
(lambda,tau)=(1,5),(2,5),(3,11),
```

the fixed-frame matrix reproduced the old compatible dimension six and the
exact diagonal.  Both lower and upper triangle stabilizers had dimension one,
their span had dimension two, and the lower stabilizer acted nontrivially on
the upper triangle.

On every representative the variable-transition matrix gave

```text
rank                                      6,
compatible dimension                     7,
relative local-parameter rank            1,
connection-parameter rank                1,
K/inter-stabilizer intersection          1,
relative dimension modulo diagonal       1.
```

An exact rational Lorentz boost intertwined the complete compatible
subspaces, not merely their dimensions.  Reversing the sign convention of
the Minkowski metric preserved every subspace.

## Mechanically independent audit

The adversarial verifier used two irregular rational reflected bipyramids
and three representatives on each.  It introduced directly

```text
16 left upper coordinates
+ 16 right upper coordinates
+ 16 unconstrained entries of A
+ 4 translations
= 52 variables.
```

It wrote the upper-edge, strut, Lorentz, lower-face-fixing and nonlinear
transition equations as polynomials and only then evaluated their Jacobian.
It did not construct a Poincare basis or import the six-dimensional kernel.

All six systems gave

```text
each separate local cell               rank/nullity 10/6,
lower-face connection stabilizer       rank/nullity 19/1,
stabilizer also fixing upper face      rank/nullity 20/0,
connection frozen in complete system   rank/nullity 46/6,
connection variable in complete system rank/nullity 45/7,
connection projection on kernel        rank 1,
shared displacement difference         rank 1.
```

All ranks were unchanged by `eta -> -eta`.

## Correction to the frozen-holonomy result

The complete 720-loop audit remains a valid theorem about a seed parallel
under a fixed background connection:

```text
the frozen connection has no nonzero parallel section in the accepted local
six-flex family.
```

What no longer follows is:

```text
therefore the complete metric-gluing kernel is zero.
```

A metric deformation can vary the face transition.  The newly derived
one-dimensional stabilizer variation supplies exactly the term omitted from
the frozen propagation equation.  Thus the earlier global result must be
read conditionally, with “fixed-frame” or “frozen-connection” stated in every
hypothesis.

This is an example of the distinction emphasized in connection and
shape-matching formulations of Regge geometry: the connection is tied to the
metric/gluing equations and cannot be held fixed without justification.

## Relation to prior art

Discrete connection matrices on codimension-one faces and curvature as their
hinge products are standard in
[Khatsymovsky](https://arxiv.org/abs/1509.04974).  The metric-dependent
discrete Levi-Civita construction is explicit in
[Khatsymovsky 2019](https://arxiv.org/abs/1906.11805).  The larger
connection/area-angle phase space before gluing constraints is discussed by
[Dittrich--Ryan](https://arxiv.org/abs/0807.2806), while Lorentzian
Levi-Civita and shape-matching constraints on general cellular decompositions
are treated by
[Anza--Speziale](https://arxiv.org/abs/1409.0836).

The post-result search also found modern Lorentzian frustum cosmology with a
dynamical height variable
([Jercher--Steinhaus](https://arxiv.org/abs/2312.11639)), but no primary
source giving the exact `6 -> 7` theorem above for a tetrahedral homothetic
frustum face.

Search absence is not novelty proof.  External novelty remains **OPEN**.

## Status ledger

| Claim | Status |
|---|---|
| Fixed-frame compatible pair is exactly diagonal dimension six | **DERIVED EXACT** |
| Lower shared triangle has a one-dimensional Poincare stabilizer | **DERIVED EXACT** |
| That derived transition variation acts on the upper triangle | **DERIVED EXACT** |
| Variable-transition compatible dimension is seven | **DERIVED EXACT, ADVERSARIALLY CORROBORATED** |
| Relative quotient and connection projections each have rank one | **DERIVED EXACT, ADVERSARIALLY CORROBORATED** |
| Frozen 720-loop holonomy matrices are incorrect | **REFUTED; THEY REMAIN EXACT** |
| Frozen holonomy alone proves metric global rigidity | **REFUTED DERIVED EXACT** |
| One independent fitted scalar per face is authorized | **REFUTED FRAMING; the line is derived, its coefficient is constrained globally** |
| Complete variable-connection 600-cell kernel | **OPEN** |
| The mode is already extrinsic curvature or a spin connection | **OPEN** |
| Schedule-free action and anisotropic Hessian | **NOT YET AUTHORIZED** |
| Dynamics, waves, effective `c` or continuum GR | **NOT TESTED** |

## Next discriminating step

Assemble the complete body-hinge-like compatibility matrix on the 600-cell
dual graph.  Associate a six-dimensional local kernel variable to every one
of the 600 tetrahedral frusta and the uniquely derived one-dimensional
transition line to every one of the 1200 shared faces.  The face equations
must be obtained by exact vertex evaluation and the actual face transports;
no independent coefficient may be fitted.

The alternatives are:

- global nullity zero: loops constrain all local face modes, rescuing metric
  infinitesimal rigidity by a stronger variable-connection theorem;
- positive nullity: classify the surviving modes before constructing an
  action;
- convention-dependent nullity: the connection/gluing definition remains
  structural and the route is not ready for dynamics.

Only after this complete variable-connection closure gate may an implicit
reconstruction map or anisotropic Hessian be attempted.
