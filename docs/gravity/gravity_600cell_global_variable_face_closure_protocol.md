# Protocol: complete variable-face closure on the regular 600-cell

Date: 2026-08-19

This protocol is committed before constructing any four-dimensional lateral-
face transition or evaluating any global rank.  No target kernel dimension is
assumed.

## Frozen provenance

| input | SHA-256 |
|---|---|
| global variable-face prior-art gate | `a76a28e8247e2fd1d0ea3536e6a345ba4a091ce931bc0e9c570f3234286b5014` |
| consolidated local variable-face theorem | `2db55cb87ec1c01d537cdbc11010bc9ea740762c598108e4c2de0f3acca72cc8` |
| primary local verifier | `69a5d7479a5df427cead76f82db31fe62a9190c28c967f699c846881634fb0f6` |
| primary local artifact | `001212016553d006862e68edc4f780f37ca1476110b6e0aed3e987f52a43b5e3` |
| adversarial local verifier | `9a3c6985eb4d833ef4ecc21f9e964577d102d5d2e7beac6a4163c4225faa5984` |
| adversarial local artifact | `c8c8c58711e5bf4e49c110e84518ddf643b75cc4377d05fb5f577003b8395466` |
| complete frozen-holonomy verifier | `54fa9775a2f14d708359167d3f8b81e03d985f24594b453f16028d9981d9be0d` |
| complete frozen-holonomy artifact | `f224fe123c882ccda97d4ca6ec67c9fd810d58ed8377c5afb457a1dec69f4b87` |

Both local artifacts must preserve their one-mode outcomes.  The old global
artifact must preserve its exact frozen-loop census while being treated only
as a conditional control.

## Exact carrier and canonical local cell

Rebuild the 120 golden-field vertices and derive all simplices by exact clique
incidence.  Require

```text
f=(120,720,1200,600),
two tetrahedra per face,
five tetrahedra per spatial edge,
connected four-regular dual graph,
dual edge-connectivity four.
```

For every sorted tetrahedron assign its four global vertices, in sorted
order, to

```text
p0=( 1, 1, 1,0), p1=( 1,-1,-1,0),
p2=(-1, 1,-1,0), p3=(-1,-1, 1,0).
```

Use `eta=diag(1,1,1,-1)`, `n=(0,0,0,1)` and

```text
q_i=lambda p_i+tau n
```

at the exact representatives

```text
(lambda,tau)=(1,5),(2,5),(3,11).
```

## Full lateral-face Lorentz transition

For a directed adjacent pair `T <- T'`, order the three shared global
vertices identically in both cells.  Four affinely independent points of the
shared lateral face are

```text
p_a, p_b, p_c, q_a.
```

In the source frame derive the Lorentz normal of their affine span.  Reflect
the source lower apex across that three-plane using the exact Minkowski
reflection formula.  The unique affine map sending

```text
target:  p'_a,p'_b,p'_c,q'_a,target lower apex
source:  p_a, p_b, p_c, q_a, reflected source lower apex
```

defines the directed transition.

For all 2400 directed transitions require exactly:

1. homogeneous affine form and `L^T eta L=eta`;
2. mapping of all three shared lower and all three shared upper vertices;
3. mapping of the target upper apex to the reflection of the source upper
   apex;
4. inverse equality for the reverse directed face;
5. nontrivial time--space mixing on at least one expanding face, so the test
   cannot collapse back to the old spatial transition by accident.

## Local face block and elimination

Let `K` be the exact six-dimensional local upper-edge-plus-strut Poincare
kernel.  In the source frame evaluate `K` on the three shared upper vertices.
Evaluate the target `K` in its own frame and transport its displacement
vectors by the transition linear part.

Let `S_f` be the exact one-dimensional Poincare stabilizer of the source
shared lower triangle.  The raw face block is

```text
B_f = [E_source K,
       -L_f E_target K,
       -E_source S_f],
```

with shape `12 x 13`.  Every face block must have rank six, kernel dimension
seven, rank-one transition projection and a six-dimensional common-motion
subspace after transport.

Eliminate the transition coefficient exactly.  Project the seven-dimensional
raw kernel to its first twelve cell coordinates and take its algebraic
annihilator.  This produces a rank-five reduced block

```text
C_f [z_source,z_target]=0,
```

equivalent to the existence of the unique face coefficient.  No row may be
selected by numerical pivot tolerance.

## Complete matrices

Assemble the reduced physical matrix

```text
C : Q^(6*600) -> Q^(5*1200),
shape 6000 x 3600.
```

Also assemble:

1. a breadth-first dual spanning-tree restriction with 599 faces;
2. an untwisted deterministic generic body--hinge control on the same dual
   graph.  For face number `r`, use hinge line

   ```text
   (1,r,r^2,r^3,r^4,r^5)
   ```

   and its five consecutive annihilators.

The disclosed control ranks/nullities are

```text
physical spanning tree: rank 2995, nullity 605;
generic full control:    rank 3594, nullity 6.
```

The generic control is evidence only that the graph/rank engine sees the
Tay--Whiteley prediction.  It is not evidence about the physical special
realization.

## Exact modular certificates

Every rational matrix entry must be reduced modulo each of the disclosed
primes

```text
p1=1000003,
p2=1000033,
```

after verifying that no denominator vanishes.  Rank is computed by exact
finite-field elimination; no floating rank participates in a verdict.

A modular rank of 3600 is an exact certificate of rational full column rank,
because it exhibits a nonzero minor over the rationals.  A modular rank
deficit alone is **not** proof of a rational kernel: bad-prime rank loss is
possible.  Any positive-kernel outcome additionally requires an explicit
rational null vector checked against the original matrix.

Record both modular ranks for every representative, the tree control and the
generic control.

## Convention attacks

1. Reverse every undirected face orientation and rebuild all blocks.
2. Apply the odd canonical relabelling `(0 1)` to every tetrahedron and
   rebuild transitions and blocks from the relabelled coordinates.
3. Repeat the static calculation with `eta -> -eta`.

All complete modular ranks must remain unchanged.  The actual reduced row
spaces need not be literally equal across gauges.

## Outcome hierarchy

1. `GLOBAL_VARIABLE_FACE_CONTROL_FAILED` if provenance, incidence,
   four-dimensional transition, local block, tree, generic, denominator or
   convention controls fail.
2. `GLOBAL_VARIABLE_FACE_KERNEL_ZERO` if controls pass and both primes give
   rank 3600 at all three representatives.
3. `GLOBAL_VARIABLE_FACE_POSITIVE_KERNEL` only if controls pass, both primes
   have the same deficit and at least one explicit nonzero rational kernel
   vector is verified for every positive representative.
4. `GLOBAL_VARIABLE_FACE_MODULAR_DEFICIT_OPEN` if controlled modular deficits
   occur without a rational kernel certificate.
5. `GLOBAL_VARIABLE_FACE_OPEN` otherwise.

## Interpretation firewall

Full rank would establish a candidate infinitesimal rigidity theorem for the
complete **variable-transition** cellular system.  It would require a
mechanically independent replication before acceptance and still would not
prove finite uniqueness.

A positive kernel would be recorded before labeling it scale, lapse, gauge or
physics.  Classification is a separate preregistered mission.

No action, Hessian, evolution equation, wave speed or full suite is
authorized by this protocol.
