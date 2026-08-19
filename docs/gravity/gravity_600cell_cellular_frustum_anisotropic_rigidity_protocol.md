# Protocol: exact local rigidity of the cellular tetrahedral frustum

Date: 2026-08-19

This target-disclosed protocol is committed before evaluating any rigidity
matrix.  It tests local metric determinacy only; no Regge action, dust term,
physical Hessian or spectrum is evaluated.

## Frozen provenance

| input | SHA-256 |
|---|---|
| prior-art gate | `92c88042e8233a542b9f21e96a99bc0d09cf13cff89a8e243354f97984baaaab` |
| balanced-slab result | `6cd15954c73129fad4ac5905bdbe4440e9ef2a748b8a41c712662a34da3599bc` |
| balanced-slab verifier | `f59b8fc89106b42077eca281ff3d956a5a5d6fb4be70b73465133035b1ce0f57` |
| balanced-slab artifact | `0a9e9e796cd671c82f2e428bfa21ba63ccb07fe76867e4553979c3c54b22a0d5` |
| tournament adversarial verifier | `794797d15f37887eb09c9aa168db73286705c47e1421c38e58a71d607099f2a8` |
| tournament adversarial artifact | `dd1043a8cb712adb4f0717f95024b9ce62132501198938bb997e7ab3dad8bf65` |
| homogeneous cellular verifier | `e88111adaeb333abf80b68e06e23d7840ef14399238ada9d0f3cd722d7934e50` |
| homogeneous cellular artifact | `640bc0dd3d6f1ae727f8113bf29514878874effffd14f539f5a43e3c3b18d069` |

The upstream balanced result must retain 24 distinct spatially invariant
orders and classification `STRUCTURAL`.  The homogeneous cellular artifact
must retain its derived coefficient outcome; it is a control, not evidence
for anisotropic rigidity.

## Exact representative family

Use the Minkowski metric

```text
eta = diag(1,1,1,-1)
```

and the centered regular tetrahedron

```text
p0=( 1, 1, 1,0)   p1=( 1,-1,-1,0)
p2=(-1, 1,-1,0)   p3=(-1,-1, 1,0).
```

For the three preregistered rational controls

```text
(lambda,tau) = (1,5), (2,5), (3,11),
```

set

```text
qi=(lambda*pi_spatial,tau).
```

All four struts are timelike, both boundary tetrahedra are nondegenerate and
the eight vertices affinely span four dimensions.  The middle point is the
expanding homothetic representative; the other two test that the rank is not
a one-point accident.

## Exact rigidity matrices

For an edge `(i,j)`, the row for its squared Minkowski length is the exact
Jacobian

```text
d ell_ij^2 / d x_i =  2 eta (x_i-x_j),
d ell_ij^2 / d x_j = -2 eta (x_i-x_j).
```

All ranks and nullspaces must be computed by SymPy over the rationals.  No
floating threshold is permitted.

The cellular graph is exactly

```text
all 6 bottom edges,
all 6 top edges,
the 4 corresponding struts (pi,qi).
```

The verifier must establish for every representative:

1. `n=8`, `m=16`, affine dimension four;
2. cellular rigidity rank 16 and full-coordinate kernel dimension 16;
3. the four translations plus six Lorentz generators are independent,
   annihilated by the rigidity matrix and span dimension 10;
4. the quotient infinitesimal-flex dimension is therefore `16-10=6`;
5. after fixing all bottom coordinates, the Jacobian of six top-edge plus
   four strut constraints has rank 10 on 16 top-coordinate variables, again
   leaving six flexes.

As a signature control, replace `eta` by the Euclidean identity and require
the same ranks.  This is expected because right multiplication by a
nonsingular block metric cannot change Jacobian rank.

## The 24 staircase completions

For each of the 24 total orders of `{0,1,2,3}`, add one cross diagonal for
each unordered pair `{i,j}`:

```text
pi -> qj  if i precedes j,
pj -> qi  otherwise.
```

Require:

1. six distinct added diagonals and 24 distinct diagonal sets;
2. full 22-edge rigidity rank 22 for every order and every representative;
3. with the bottom fixed, the stacked 16-by-16 Jacobian of top edges, struts
   and those six diagonals has rank 16;
4. the six diagonal differentials restricted to the six-dimensional fixed-
   bottom flex kernel have rank six;
5. time reversal maps every diagonal set to the set associated with the
   reversed order, giving twelve two-element orbits and no fixed order.

The third condition directly tests that the six chosen diagonals remove the
six missing local shapes rather than merely increasing a row count.

## Outcome hierarchy

1. `CELLULAR_FRUSTUM_RIGIDITY_CONTROL_FAILED` if provenance, topology,
   signature, affine-span, isometry or exact-arithmetic controls fail.
2. `CELLULAR_FRUSTUM_SIX_SHAPES_UNDERDETERMINED` if the cellular graph leaves
   exactly six non-isometric flexes and all 24 six-diagonal completions remove
   them at all three representatives.
3. `CELLULAR_FRUSTUM_UNDERDETERMINATION_WORSE` if the cellular graph leaves
   more than six flexes or any staircase completion remains flexible.
4. `CELLULAR_FRUSTUM_RIGIDITY_OPEN` otherwise.

## Interpretation firewall

The expected negative closes the claim that the existing 16 cellular lengths
define a unique anisotropic classical block.  It does not invalidate the
homogeneous action, where the homothetic ansatz supplies the missing shape.
It also does not select a triangulation, authorize an equal-weight ensemble
or forbid a future first-order area/angle theory.  Each of those adds a new
principle and must be preregistered separately.

Only this verifier and static registry guards may be run.

