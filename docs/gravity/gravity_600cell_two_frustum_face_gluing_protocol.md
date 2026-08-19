# Protocol: exact fixed-frame gluing of two homothetic frusta

Date: 2026-08-19

This target-disclosed protocol is committed before evaluating any two-body
gluing matrix.  It asks only whether the accepted six local flexes already
contain an independent shared-face mode.

## Frozen provenance

| input | SHA-256 |
|---|---|
| two-frustum prior-art gate | `d38994ada998df4736858b1a242802097134383d870135bd542178c5343ff63b` |
| consolidated local theorem | `436fb57037e491b6bdb8fee9ad8b10ab8da1621fd9ecda73e1fcac3fa616fa29` |
| invariant correction verifier | `e85e2df690234e19e0343183499c6f8465bc149bc66a87c5246dfe0bda4c1d61` |
| invariant correction artifact | `f571869be3341b74b2341c2bf776e99b21174f9f0fb0c5d02e42585c2f3ebaa2` |
| adversarial local verifier | `8f9bd9882f2efa9f7fdc415ad3c6ca13927283ddb548a8d779c8955ad8ff7e21` |
| adversarial local artifact | `b750943349dc60ee42d08c0ba61d9a7a0838e3f9346ac1d4397000519b4d6395` |

Both accepted local artifacts must retain `13/13` and their exact
static-versus-expanding stratification.

## Exact triangular-bipyramid carrier

Use `eta=diag(1,1,1,-1)`, `n=(0,0,0,1)` and

```text
p0=( 1, 1, 1,0)   p1=( 1,-1,-1,0)
p2=(-1, 1,-1,0)   p3=(-1,-1, 1,0)
p4=(5/3,5/3,-5/3,0).
```

The left tetrahedron is `(0,1,2,3)` and the right tetrahedron is
`(0,1,2,4)`.  They are congruent regular tetrahedra on opposite sides of the
shared face `(0,1,2)`; `p4` is the exact reflection of `p3` across that face.

Set

```text
q_i=lambda p_i+tau n
```

at

```text
(lambda,tau)=(1,5),(2,5),(3,11).
```

All geometry and ranks must use exact SymPy rational arithmetic.

## Local kernels

For each tetrahedron independently reconstruct:

1. its ten-dimensional Poincare displacement matrix `U_T`;
2. its four differentiated strut constraints `C_T`;
3. the exact six-dimensional parameter kernel `K_T=ker(C_T)`;
4. equality with the direct kernel of its six top-edge plus four strut
   squared-length Jacobian.

In the common developed frame, require `K_L=K_R` as subspaces of the same
ten Poincare parameters.  Also compare them with the analytic forms:

```text
lambda!=1: b(A)=tau/(lambda-1) A n;
lambda=1:  A n=0 and <b,n>=0.
```

Failure is a control failure, not evidence about gluing.

## Positive full-Poincare face control

Restrict a general Poincare Killing field `(A,b)` to the three shared upper
vertices.  Its `12 x 10` evaluation matrix `F` must have

```text
rank(F)=9, dim ker(F)=1.
```

The one-dimensional kernel must have Lorentz rank one and must fix the three
face points exactly.  For two unrestricted Poincare bodies, the equality
matrix `[F,-F]` must therefore have kernel dimension eleven: ten common
motions plus one relative pointwise-face stabilizer.

This is the positive control that the verifier can see a face-holonomy-like
relative direction before strut constraints are imposed.

## Decisive constrained gluing matrix

Let the columns of `K_L` and `K_R` parameterize the two local six-flex
spaces.  Impose equality of the three shared upper-vertex displacements:

```text
G = [F K_L, -F K_R].
```

Compute `ker(G)` exactly and map it back to the two ten-parameter Poincare
vectors.  Record:

1. total compatible-pair dimension;
2. rank of the difference `(A_L,b_L)-(A_R,b_R)` on that kernel;
3. dimension of `ker(F K_T)`, the relative face stabilizer inside one local
   strut-preserving kernel;
4. equality of the compatible-pair space with the diagonal
   `{(z,z): z in K_L=K_R}`.

The disclosed prediction is

```text
rank(G)=6,
compatible-pair dimension=6,
relative difference rank=0,
dim ker(F K_T)=0.
```

## Outcome hierarchy

1. `TWO_FRUSTUM_FACE_GLUING_CONTROL_FAILED` if provenance, congruence,
   local-kernel reconstruction or the one-dimensional full-Poincare face
   stabilizer fails.
2. `TWO_FRUSTUM_HIDDEN_FACE_MODE` if the constrained compatible space has
   dimension seven and its unique relative direction is the full-Poincare
   pointwise face stabilizer.
3. `TWO_FRUSTUM_DIAGONAL_ONLY` if the compatible space has dimension six,
   every pair is diagonal and the constrained relative stabilizer is zero.
4. `TWO_FRUSTUM_FACE_UNDERDETERMINED` if more than one relative direction
   survives.
5. `TWO_FRUSTUM_FACE_GLUING_OPEN` otherwise.

## Interpretation firewall

`TWO_FRUSTUM_DIAGONAL_ONLY` reaches a local kill boundary only for the claim
that the six already-derived length flexes themselves provide an independent
face connection.  It does not forbid a first-order theory with a new frame
transition or holonomy variable.  Conversely, one surviving direction would
be only a candidate face mode; it would still require an independent audit,
closure around edges and an action.

No global 600-cell propagation is authorized before adversarial replication
of a material two-frustum result.  Only the new verifier and static registry
guards may be run.
