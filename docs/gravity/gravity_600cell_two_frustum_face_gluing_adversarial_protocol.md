# Adversarial protocol: direct five-vertex audit of two-frustum gluing

Date: 2026-08-19

This protocol is committed after the primary `TWO_FRUSTUM_DIAGONAL_ONLY`
artifact and before the independent union Jacobian is evaluated.  It does
not import the primary gluing construction or form a pair of local
Poincare-kernel bases.

## Frozen inputs

| input | SHA-256 |
|---|---|
| primary gluing protocol | `7d6f6028b6585bc472ee25aca455194d2ac13ed61fc31f7e0f339f4a9bf697f8` |
| primary gluing verifier | `52636ae59bd4e4568df175e32b7c3aeae4fbfbc3d475d255131b6db671c41ae7` |
| primary gluing artifact | `0e09c3f8f38c8158deff5b81bc6fe4d5d6dd685a24cce83e015fb95e3f26a70e` |
| consolidated local theorem | `436fb57037e491b6bdb8fee9ad8b10ab8da1621fd9ecda73e1fcac3fa616fa29` |

The frozen artifact must retain `9/9`, the one-dimensional unrestricted
face stabilizer, and `TWO_FRUSTUM_DIAGONAL_ONLY` at all three representatives.

## Irregular reflected carrier

Use the irregular shared face and left apex

```text
p0=(5,0,0,0), p1=(0,5,0,0), p2=(0,0,5,0),
p3=(3,4,0,0).
```

Reflect `p3` exactly across the plane of `(p0,p1,p2)` to obtain

```text
p4=(5/3,8/3,-4/3,0).
```

The two tetrahedra are `(0,1,2,3)` and `(0,1,2,4)`.  At

```text
(lambda,tau)=(1,7),(2,7),(3,13),
q_i=lambda p_i+tau n.
```

all five struts must be timelike.  No equal-strut assumption is made; the
reflected apex has a different radius from the four original equinorm
vertices.

## Direct union-polynomial construction

Use twenty symbolic upper coordinates, one four-vector for each of the five
physical upper vertices.  Construct exactly fourteen squared-length
polynomials:

```text
9 distinct top edges in the two-tetrahedron union,
5 corresponding struts.
```

Differentiate them directly to form one `14 x 20` polynomial Jacobian.  The
shared face vertices occur only once; gluing is built into the carrier rather
than imposed by comparing two local parameter spaces.

Require at every representative:

```text
rank=14, nullity=6.
```

Construct the analytic common Poincare displacement space on all five upper
vertices using

```text
lambda!=1: b(A)=tau/(lambda-1) A n;
lambda=1:  A n=0, <b,n>=0.
```

The direct twenty-coordinate polynomial kernel must equal this image exactly.
Dimension equality alone is insufficient.

## Mechanically independent positive face control

Use all sixteen entries of an unconstrained matrix `A` plus four translations
`b`.  Impose separately:

```text
the ten independent equations A^T eta+eta A=0,
A q_i+b=0 for the three shared face vertices.
```

The resulting redundant exact system must have rank nineteen and nullity one.
This reconstructs the full-Poincare pointwise-face stabilizer without the
primary six-generator evaluation matrix.

## Convention attacks

1. Replace `eta` by `-eta`; the fourteen-row direct kernel must be identical.
2. Apply the rational boost
   `boost_z(cosh=5/4,sinh=3/4)` to all five bottom and upper vertices; the
   boosted direct kernel must equal the block-Lorentz transport of the
   original kernel.
3. The full-Poincare pointwise-face stabilizer must remain one-dimensional
   after the boost.

## Outcome hierarchy

1. `ADVERSARIAL_TWO_FRUSTUM_CONTROL_FAILED` if provenance, reflection,
   Lorentzian struts, polynomial completeness, full-Poincare positive control
   or convention attacks fail.
2. `ADVERSARIAL_TWO_FRUSTUM_DIAGONAL_ONLY` if every direct union kernel is
   exactly the six-dimensional common Poincare image.
3. `ADVERSARIAL_TWO_FRUSTUM_DISAGREEMENT` if controls pass but any extra or
   missing direct union direction appears.
4. `ADVERSARIAL_TWO_FRUSTUM_OPEN` otherwise.

## Interpretation firewall

Corroboration closes only a hidden-connection reading of the current
length-plus-strut flexes.  It supports, but does not itself select, moving to
a genuine first-order carrier with independent face Lorentz/Poincare
transition data and closure/shape-matching constraints.  No global carrier,
action or dynamics is tested here.

Only this verifier and static registry guards may be run.
