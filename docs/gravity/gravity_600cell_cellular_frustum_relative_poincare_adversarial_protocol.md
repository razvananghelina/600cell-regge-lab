# Adversarial protocol: irregular-frustum Poincare stratification

Date: 2026-08-19

This protocol is committed after the invariant correction result but before
the adversarial matrices are evaluated.  It does not import either previous
Poincare implementation and does not reuse the regular tetrahedron.

## Frozen inputs

| input | SHA-256 |
|---|---|
| prior-art gate | `a8811a441fecd137b37085e4018fea7abb3f365750dfc426d16f1b46e5282e7c` |
| correction protocol | `ad569b4c4ecbfb4b6d3db7d0225dd26f0426e7e973ddb9c82c21daff5b404b3c` |
| correction verifier | `e85e2df690234e19e0343183499c6f8465bc149bc66a87c5246dfe0bda4c1d61` |
| correction artifact | `f571869be3341b74b2341c2bf776e99b21174f9f0fb0c5d02e42585c2f3ebaa2` |
| independent rigidity verifier supplying the irregular control geometry | `ecc5e0cb5f8913325f00137245f33299c7607b395d219161bf7e0e806068c18a` |

The correction artifact must retain `13/13` and
`STATIC_STABILIZER_AND_EXPANDING_LORENTZ_CHART_CORROBORATED`.

## Mechanically different carrier

Use the irregular equal-radius rational tetrahedron

```text
p0=(5,0,0,0), p1=(0,5,0,0),
p2=(0,0,5,0), p3=(3,4,0,0)
```

and the three independently frozen representatives

```text
(lambda,tau)=(1,7),(2,7),(3,13),
q_i=lambda p_i+tau e_t.
```

Rebuild the six top-edge and four strut squared-length polynomials in 16
symbolic top coordinates and differentiate them exactly.  Their direct
ten-row Jacobian is the flex-space definition for this audit.

## Redundant affine-Lorentz reconstruction

Do not insert a six-generator Lorentz basis into the displacement matrix.
Instead use twenty unknowns:

```text
all 16 entries of A and all 4 entries of b,
delta q_i=A q_i+b.
```

Impose the ten independent linear equations in

```text
A^T eta + eta A = 0
```

together with the four differentiated strut equations.  The resulting
`14 x 20` exact system must have nullity six.  Map that kernel to the 16
vertex displacements and compare its image directly with the polynomial
ten-length kernel.  Equality of dimension alone is insufficient.

This construction has different parameter redundancy, matrix shapes and
kernel route from the primary calculation.

## Frozen classification tests

For each representative compute:

- the dimension of the image of the `A` part in the 16-dimensional matrix
  coordinate space;
- the dimension of the pure-translation kernel obtained by setting `A=0`;
- the rank and determinant of the four-strut translation block.

Require:

```text
lambda=1:   dim image(A)=3, pure translations=3, rank(T)=1;
lambda!=1:  dim image(A)=6, pure translations=0, rank(T)=4.
```

At the static representative, the `A` image must equal exactly

```text
{A : A^T eta+eta A=0 and A e_t=0}.
```

At expanding representatives it must equal all of `so(3,1)`.

As a convention attack, apply the different rational boost

```text
L = boost_y(cosh=13/12,sinh=5/12)
```

and require direct equality of the boosted static `A` image with the
stabilizer of `L e_t`.  Repeat all decision ranks with metric `-eta`; they
must be unchanged.  The coordinate rotation/boost split is deliberately not
used.

## Symbolic determinant control

For general `lambda,tau`, factor the irregular tetrahedron's translation-
block determinant.  It must be a nonzero constant times

```text
tau (lambda-1)^3.
```

This tests whether the regular result was caused by tetrahedral symmetry.

## Outcome hierarchy

1. `ADVERSARIAL_POINCARE_STRATIFICATION_CONTROL_FAILED` if provenance,
   irregular-tetrahedron nondegeneracy, polynomial construction, metric-sign
   or boost controls fail.
2. `ADVERSARIAL_POINCARE_STRATIFICATION_CORROBORATED` if the direct
   polynomial kernel equals the redundant affine-Lorentz kernel and every
   static/expanding classification above passes.
3. `ADVERSARIAL_POINCARE_STRATIFICATION_DISAGREEMENT` if controls pass but
   either kernel equality or classification differs.
4. `ADVERSARIAL_POINCARE_STRATIFICATION_OPEN` otherwise.

## Interpretation firewall

Corroboration would establish an exact local kinematic theorem for this
homothetic frustum family, not a gravitational connection.  It would also
establish that the static singular stratum is generic within the tested
family rather than an artifact of the regular tetrahedron.  Face gluing,
closure, torsion, simplicity, shape matching, action and dynamics remain
**OPEN**.

Only this verifier and static registry guards may be run.
