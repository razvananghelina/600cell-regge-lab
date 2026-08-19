# Protocol: exact variable-transition gluing of two homothetic frusta

Date: 2026-08-19

This protocol is committed before evaluating any variable-transition rank.
It asks whether the fixed-frame diagonal-only theorem survives when the
derived pointwise stabilizer of the shared lower triangle is allowed to vary
the face connection.

## Frozen provenance

| input | SHA-256 |
|---|---|
| variable-connection prior-art gate | `2ed809fedad24fa15977b39e4dd6fec386e9080c123208d54fd089554ce44d2d` |
| local relative-Poincare theorem | `436fb57037e491b6bdb8fee9ad8b10ab8da1621fd9ecda73e1fcac3fa616fa29` |
| two-frustum fixed-frame result | `b5bb18c75ea1359d33b9985ad5816c21f437960c06f8c4eae793a3505509add3` |
| two-frustum fixed-frame artifact | `0e09c3f8f38c8158deff5b81bc6fe4d5d6dd685a24cce83e015fb95e3f26a70e` |
| frozen-connection global result | `72c8b2c0ffbd9d13aef8f14404270cac29896c876ed6f015a4dc7a41a89b6535` |
| complete fixed-connection artifact | `f224fe123c882ccda97d4ca6ec67c9fd810d58ed8377c5afb457a1dec69f4b87` |

The fixed-frame artifact must retain `9/9` and `TWO_FRUSTUM_DIAGONAL_ONLY`.
The complete global artifact must retain `11/11` and
`ADVERSARIAL_GLOBAL_FLEX_SEED_KILLED`.  Those controls preserve the old
conditional results; they do not decide the present extension.

## Exact carrier

Use `eta=diag(1,1,1,-1)`, `n=(0,0,0,1)` and the same reflected regular
triangular bipyramid as the fixed-frame audit:

```text
p0=( 1, 1, 1,0)   p1=( 1,-1,-1,0)
p2=(-1, 1,-1,0)   p3=(-1,-1, 1,0)
p4=(5/3,5/3,-5/3,0).
```

The cells `(0,1,2,3)` and `(0,1,2,4)` share `(0,1,2)`.  Evaluate the exact
representatives

```text
(lambda,tau)=(1,5),(2,5),(3,11),
q_i=lambda p_i+tau n.
```

All ranks and intersections use exact SymPy rational arithmetic.

## Local kernels and fixed-frame control

Independently reconstruct for each cell the direct Jacobian of its six upper
edges and four struts.  Require rank ten, nullity six, and equality with the
evaluation of the analytic local Poincare kernel `K`:

```text
lambda != 1: b(A)=tau/(lambda-1) A n;
lambda = 1:  A n=0, <b,n>=0.
```

Let `F_plus` evaluate a Poincare generator on the three shared upper
vertices.  The frozen matrix

```text
G_fixed = [F_plus K, -F_plus K]
```

must have rank six, kernel dimension six, and exact diagonal compatible
space.  This reproduces rather than assumes the earlier outcome.

## Derived connection variation

Let `F_minus` evaluate a Poincare generator on the three shared lower
vertices and define

```text
S_minus = ker(F_minus),
S_plus  = ker(F_plus).
```

Require both to be one-dimensional pointwise triangle stabilizers with
nonzero Lorentz part.  The two shared triangles are parallel but distinct,
so require

```text
dim span(S_minus,S_plus)=2,
rank(F_plus S_minus)=1.
```

The second condition proves that the lower-face connection variation acts
nontrivially on the upper triangle; it is not an unobservable coefficient.

Allow exactly one scalar `xi` multiplying `S_minus`.  The variable-transition
matching matrix is

```text
G_var = [F_plus K, -F_plus K, -F_plus S_minus].
```

For its exact kernel record:

1. total compatible dimension;
2. rank of the relative local parameter `y_left-y_right`;
3. rank of the connection coefficient `xi`;
4. dimension of

   ```text
   K intersect span(S_minus,S_plus);
   ```

5. whether quotienting by the six-dimensional common diagonal leaves zero,
   one or more relative directions.

## Frozen prediction and outcome hierarchy

The disclosed analytic prediction is

```text
dim(K intersect span(S_minus,S_plus)) = 1,
rank(G_var)                             = 6,
compatible dimension                   = 7,
relative local-parameter rank          = 1,
connection-coefficient rank            = 1,
relative dimension modulo diagonal     = 1.
```

This prediction follows from the two parallel triangle stabilizers differing
by a translation and is frozen before matrix evaluation.

Assign exactly one outcome:

1. `VARIABLE_FACE_CONNECTION_CONTROL_FAILED` if provenance, carrier,
   local-kernel, fixed-frame or triangle-stabilizer controls fail.
2. `VARIABLE_FACE_CONNECTION_FORCED_ZERO` if controls pass and the variable
   compatible space remains exactly the six-dimensional diagonal with
   `xi=0`.
3. `ONE_CONNECTION_COUPLED_RELATIVE_MODE` if all six disclosed numerical
   predictions hold on all three strata.
4. `VARIABLE_FACE_CONNECTION_UNDERDETERMINED` if more than one relative
   direction survives.
5. `VARIABLE_FACE_CONNECTION_OPEN` otherwise.

## Covariance and falsification controls

Apply one exact rational Lorentz boost to all points, local kernels and
stabilizers.  Rebuild every evaluation matrix in the boosted coordinates and
require identical invariant dimensions and exact intertwining of the
compatible spaces.

Repeat with the metric convention `eta -> -eta`; all squared-length
Jacobians may change sign, but every subspace and rank decision must remain
unchanged.

## Interpretation firewall

`ONE_CONNECTION_COUPLED_RELATIVE_MODE` would refute only the broad inference
from frozen-connection holonomy to metric rigidity.  It would not refute any
of the exact 720 frozen-loop matrices.  It would instead show that global
closure must include one derived transition variation per shared face and
the linearized compatibility equations that relate them.

Conversely, `VARIABLE_FACE_CONNECTION_FORCED_ZERO` would justify the local
fixed-frame step but would still not prove finite nonlinear uniqueness.

No global Hessian, dynamics, propagation speed or full suite is authorized
by this protocol.
