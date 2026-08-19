# Adversarial protocol: direct polynomial variable-face gluing

Date: 2026-08-19

This protocol is committed after the primary variable-transition result and
before evaluating any independent polynomial Jacobian.  It must not import
the primary Poincare basis, analytic six-kernel, stabilizer intersection or
compatible-space construction.

## Frozen inputs

| input | SHA-256 |
|---|---|
| prior-art gate | `2ed809fedad24fa15977b39e4dd6fec386e9080c123208d54fd089554ce44d2d` |
| primary protocol | `f6b91206a857cda6ebfe5cb9988110de5f12a9c1ca51bcbdb733a8429682ca6a` |
| primary verifier | `69a5d7479a5df427cead76f82db31fe62a9190c28c967f699c846881634fb0f6` |
| primary artifact | `001212016553d006862e68edc4f780f37ca1476110b6e0aed3e987f52a43b5e3` |
| fixed-frame result | `b5bb18c75ea1359d33b9985ad5816c21f437960c06f8c4eae793a3505509add3` |
| irregular fixed-frame artifact | `0f8e70ef89b7fd5a8995349d40c77f6d3f637f2d9ce137ce2c9ff07b2fed2542` |

The primary artifact must retain `11/11` and
`ONE_CONNECTION_COUPLED_RELATIVE_MODE`.  The old irregular artifact must
retain `11/11` and `ADVERSARIAL_TWO_FRUSTUM_DIAGONAL_ONLY`.

## Two independent irregular carriers

Use two rational triangular bipyramids whose shared face lies in `z=0` and
whose apices are exact reflections across that plane:

```text
carrier A:
p0=( 0, 0, 0,0)  p1=( 2, 0, 0,0)  p2=(1,3,0,0)
p3=( 1,-1, 4,0)  p4=( 1,-1,-4,0)

carrier B:
p0=(-1, 1, 0,0)  p1=( 3, 1, 0,0)  p2=(0,4,0,0)
p3=( 2,-2, 7,0)  p4=( 2,-2,-7,0).
```

Require nonzero tetrahedral affine volumes and nonuniform edge-length
multisets.  Evaluate

```text
(lambda,tau)=(1,7),(2,7),(4,13),
q_i=lambda p_i+tau n.
```

No floating arithmetic or rank tolerance is permitted.

## Direct polynomial system

Introduce 52 independent variables:

```text
16 left upper-vertex coordinates,
16 right upper-vertex coordinates,
16 unconstrained entries of A,
4 entries of b.
```

Construct the equations as polynomials before differentiating:

1. six upper-edge and four strut squared lengths for the left cell;
2. the same ten polynomials for the right cell;
3. the ten independent entries of `A^T eta + eta A`;
4. twelve equations `A p_i+b=0` on the shared lower triangle;
5. twelve exact transition equations

   ```text
   x_i(left) = (I+A) x_i(right)+b
   ```

   on the shared upper triangle.

Evaluate the full polynomial Jacobian at the homothetic point and `A=b=0`.
The decisive expected rank/nullity is

```text
45/7.
```

This method has no ten-parameter Poincare coordinate chart and no imported
six-dimensional kernel.

## Controls and projections

For each cell separately, the direct ten-polynomial Jacobian must have
rank/nullity `10/6`.

For the connection variables alone, the Lorentz equations plus lower-face
fixing must have rank/nullity `19/1`.  Adding pointwise fixation of the upper
shared triangle must raise this to `20/0`.  This is the positive/negative
control that the one connection variation is real and acts nontrivially.

Freeze all twenty connection variables at zero in the complete system.  The
result must have rank/nullity `46/6`, reproducing the old fixed-frame count in
redundant variables.

On the seven-dimensional variable-system kernel require:

```text
rank of the 20 connection coordinates              = 1,
rank of the 12 shared-upper displacement difference = 1.
```

Repeat every rank after `eta -> -eta`.  Because the entire calculation is
exact, a sign disagreement is a control failure rather than a tolerance
issue.

## Outcome hierarchy

1. `ADVERSARIAL_VARIABLE_FACE_CONTROL_FAILED` if provenance, irregularity,
   local ranks, stabilizer controls, frozen count or metric-sign control
   fails.
2. `ADVERSARIAL_ONE_CONNECTION_MODE` if all six carrier/representative
   systems have `45/7` and both decisive projections have rank one.
3. `ADVERSARIAL_CONNECTION_FORCED_ZERO` if all controls pass but every
   variable system remains at nullity six with zero connection projection.
4. `ADVERSARIAL_VARIABLE_FACE_UNDERDETERMINED` if any controlled system has
   nullity greater than seven or either decisive projection rank exceeds one.
5. `ADVERSARIAL_VARIABLE_FACE_OPEN` otherwise.

## Interpretation firewall

Corroboration refutes the inference

```text
zero frozen-holonomy fixed space => zero metric-gluing kernel.
```

It does not decide the complete 600-cell system with 1200 variable face
connections.  That global incidence system requires a new preregistration;
neither independent face coefficients nor arbitrary holonomies may be fitted.

No action, Hessian, dynamics or full suite is authorized here.
