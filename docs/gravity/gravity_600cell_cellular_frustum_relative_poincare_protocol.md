# Protocol: exact Poincare decomposition of the six cellular-frustum flexes

Date: 2026-08-19

This protocol is committed before evaluating the Poincare constraint matrix.
It tests a local kinematic identification only.  It contains no Regge action,
dust, fitted coefficient, global 600-cell gluing or physical spectrum.

## Frozen provenance

| input | SHA-256 |
|---|---|
| relative-Poincare prior-art gate | `a8811a441fecd137b37085e4018fea7abb3f365750dfc426d16f1b46e5282e7c` |
| cellular-frustum consolidated result | `7221fd948c6e21aa59ea2738d4ef6a224f13b674689a051837366fcfa76f203b` |
| primary rigidity verifier | `2f766503296aa43f6192d2cce6ce44faac3b7fb57ba131ba0fbf393a2da80f60` |
| primary rigidity artifact | `c55f98313121018ff5ca1fc834260e8f2f075248a21fd7b99a356d89b2d18255` |
| adversarial rigidity verifier | `ecc5e0cb5f8913325f00137245f33299c7607b395d219161bf7e0e806068c18a` |
| adversarial rigidity artifact | `7763287a12075a911134b24e5f23c3c682198923bda1ab8f75ac1d9541540fc1` |

The two upstream artifacts must retain their exact six-flex outcomes.  They
are controls; their six-dimensional count is not evidence for the new
interpretation.

## Exact carrier and conventions

Use

```text
eta = diag(1,1,1,-1)

p0=( 1, 1, 1,0)   p1=( 1,-1,-1,0)
p2=(-1, 1,-1,0)   p3=(-1,-1, 1,0)

q_i = lambda p_i + tau e_t
```

at the already frozen rational representatives

```text
(lambda,tau) = (1,5), (2,5), (3,11).
```

The six Lorentz generators are ordered as three spatial rotations followed
by three boosts:

```text
J01, J02, J12, K03, K13, K23,
```

where every generator `A` obeys `A^T eta + eta A = 0`.  Append four common
translations in coordinate order.  These ten columns define a matrix `U`
from Poincare parameters to the stacked 16 top-vertex displacements:

```text
U(A,b)_i = A q_i + b.
```

All matrices, ranks, nullspaces, determinants and factorisations must use
exact SymPy rational/polynomial arithmetic.  No numerical rank tolerance is
permitted.

## Control: every top-edge flex is a restricted Poincare motion

Let `E` be the exact `6 x 16` Jacobian of the six top squared lengths, and
let `S` be the exact `4 x 16` Jacobian of the four strut squared lengths.
For every representative require:

1. `rank(E)=6`, hence `dim ker(E)=10`;
2. `rank(U)=10`;
3. `E U=0` exactly;
4. `col(U)=ker(E)`, checked by both containment and dimension;
5. `rank([E;S])=10`, with a six-dimensional fixed-bottom flex kernel.

Failure here invalidates every interpretation below.

## Decisive Poincare constraint map

Construct, rather than postulate,

```text
C = S U,
C(A,b)_i = 2 <q_i-p_i,A q_i+b>_eta.
```

Let `C_L` be its first six Lorentz columns and `T` its last four translation
columns.  Let `K` be an exact basis matrix for `ker(C)`.  Require:

1. `rank(C)=4`, `dim ker(C)=6`;
2. `U K` spans exactly `ker([E;S])`;
3. the rank of the Lorentz projection is the rank of the first six rows of
   `K`;
4. the pure-translation flex dimension is `4-rank(T)`.

The test is basis-independent because only images, kernels and exact ranks
are classified.

## Frozen static-versus-expanding predictions

No corresponding matrix has been evaluated while writing this protocol.
The disclosed geometry predicts a stratified answer.

At `lambda=1`:

```text
rank(T)                              1
pure-translation flex dimension      3
rank(ker(C) -> so(3,1))              3
Lorentz image                        spatial rotations only.
```

Thus the static six-flex space should be rotations plus common spatial
translations of the upper tetrahedron relative to the lower one, not a copy
of all `so(3,1)`.

At each `lambda!=1` representative:

```text
rank(T)                              4
pure-translation flex dimension      0
rank(ker(C) -> so(3,1))              6.
```

Then `T` is invertible and

```text
b(A) = -T^(-1) C_L A
```

gives a unique local graph over `so(3,1)`.  The symbolic determinant of `T`
must also be computed for general `lambda,tau` and factorized; its zero set
must explain, rather than numerically hide, the static degeneration.

## Covariance and anti-convention controls

Repeat the exact classification after each of the following transformations:

1. shift both tetrahedra by the rational origin vector `(2,-1,3,4)` and
   verify the parameter law `A'=A`, `b'=b-A r` exactly;
2. apply the rational Lorentz boost with `(cosh,sinh)=(5/4,3/4)` in the
   first-space/time plane and verify conjugation `A'=L A L^(-1)`, `b'=L b`;
3. simultaneously relabel the four paired vertices by `(0 2 3 1)`.

The physical displacement space and all decision ranks must be unchanged.
The compensating translation coefficients themselves are not invariants and
must not be promoted to observables.

As a falsifying control, the static case must retain its nonzero pure-
translation kernel.  A verifier reporting six Lorentz directions there has
mixed a dimension count with an isomorphism or silently fixed translations.

## Outcome hierarchy

1. `RELATIVE_POINCARE_CONTROL_FAILED` if provenance, exact arithmetic,
   rigid-motion completeness, covariance or upstream controls fail.
2. `UNIFORM_RELATIVE_LORENTZ_CHART` if the Lorentz projection is an
   isomorphism at all three representatives.
3. `STRATIFIED_RELATIVE_LORENTZ_CHART` if both expanding representatives
   are exact Lorentz graphs but the static representative has the frozen
   three-dimensional pure-translation sector and rank-three rotational
   image.
4. `RELATIVE_LORENTZ_INTERPRETATION_REFUTED` if the controls pass but either
   expanding representative is not a Lorentz graph.
5. `RELATIVE_LORENTZ_INTERPRETATION_OPEN` otherwise.

## Interpretation firewall

Even `UNIFORM` would establish only a local coordinate chart on the flex
space.  `STRATIFIED` would refute a uniform identification through the
static slice while retaining a chart on the expanding stratum.  Neither
outcome derives:

- a Lorentz connection or spin connection on the 600-cell;
- extrinsic curvature, a symplectic form or canonical momentum;
- simplicity, closure, torsion-free or shape-matching equations;
- a first-order action, evolution, propagation speed or continuum limit.

Only after this local result survives an independent construction may a
separate face-gluing/closure protocol be written.

Only the new verifier and static registry guards may be run.
