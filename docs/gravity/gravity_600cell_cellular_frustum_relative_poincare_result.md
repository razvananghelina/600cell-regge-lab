# Result: the six frustum flexes form a stratified relative-Poincare kernel

Date: 2026-08-19

## Headline

The six missing infinitesimal shapes of a homothetic tetrahedral spacetime
frustum are not merely a dimension match.  They are exactly the relative
Poincare motions of the upper tetrahedron that preserve the four struts.
However, their relation to the Lorentz algebra is stratified:

```text
lambda != 1:
    the six-flex kernel is a graph over all so(3,1);

lambda = 1:
    the kernel is three observer-normal rotations
    plus three relative spatial translations.
```

Therefore:

```text
DERIVED EXACT LOCAL:
the six cellular flexes have a constrained relative-Poincare meaning.

DERIVED EXACT NEGATIVE:
they are not a uniform Lorentz-frame variable through the static slice.

OPEN:
whether global gluing turns the local variables into a Levi-Civita/spin
connection or extrinsic curvature.
```

## Provenance ledger

| stage | commit | outcome |
|---|---|---|
| prior-art gate | `403059c` | no automatic connection interpretation |
| primary protocol | `f2dc6d3` | static/expanding split disclosed |
| primary registration | `7a7e411` | no matrix evaluated yet |
| preserved first artifact | `57d4a7f` | `12/13`, covariance control failed |
| correction protocol | `5bd3b5b` | invariant observer-stabilizer criterion |
| correction registration | `0de0a74` | failed result kept frozen |
| correction artifact | `618455a` | `13/13`, classifier error isolated |
| adversarial protocol | `6250871` | irregular tetrahedron, redundant variables |
| adversarial registration | `a0f6274` | no adversarial matrix evaluated yet |
| adversarial artifact | `a2524c3` | `13/13`, stratification corroborated |

Artifact hashes:

```text
preserved failed primary
3ac5cce9db2b2f828e0ced2114f301f761dd9371847b712ad47119709396cf7d

invariant correction
f571869be3341b74b2341c2bf776e99b21174f9f0fb0c5d02e42585c2f3ebaa2

independent adversarial audit
b750943349dc60ee42d08c0ba61d9a7a0838e3f9346ac1d4397000519b4d6395
```

## Exact local theorem

Let `V` be Minkowski four-space, let `n` be a unit timelike vector, and let
`p_0,...,p_3` be affinely independent points in `n^perp`.  Define

```text
q_i = lambda p_i + tau n,    tau != 0.
```

Fix the six intrinsic upper-tetrahedron lengths and the four corresponding
strut lengths.  Every infinitesimal upper displacement preserving its six
intrinsic lengths has a unique representation

```text
delta q_i = A q_i + b,
A in so(3,1), b in V.
```

The four struts impose

```text
C(A,b)_i = 2 <q_i-p_i,A q_i+b> = 0.
```

The solution space has dimension six.  More precisely:

1. if `lambda != 1`, its projection to `so(3,1)` is an isomorphism;
2. if `lambda = 1`, it is

```text
{A : A n=0} direct-sum {b : <b,n>=0},
```

with dimensions `3+3`.

### Proof sketch

The differences `q_i-q_0` span a nondegenerate spacelike three-plane.  The
linearized six-distance equations say that the induced map on this plane has
zero metric-symmetric part.  Its tangential skew part has dimension three
and its normal part has dimension three; these extend uniquely to
`A in so(3,1)`.  The displacement of one vertex then fixes `b`, proving the
unique Poincare representation.

The translation block of `C` has rows proportional to

```text
((lambda-1) p_i + tau n)^T eta.
```

Its determinant is a nonzero multiple of

```text
tau (lambda-1)^3 Vol_aff(p_0,p_1,p_2,p_3).
```

Thus it is invertible for `lambda!=1`, and every `A` has a unique
compensating `b(A)`.

At `lambda=1`, subtracting two strut equations gives

```text
<n,A(p_i-p_j)> = -<A n,p_i-p_j> = 0.
```

The differences span `n^perp`, while `A n` already lies in `n^perp`, so
`A n=0`.  The remaining common equation is `<n,b>=0`.  This yields exactly
the two three-dimensional summands above.

## Exact computations

For the centered regular tetrahedron, the symbolic translation determinant
was

```text
det(T) = -256 tau (lambda-1)^3.
```

The primary ranks were

| `(lambda,tau)` | `rank(T)` | pure translations | Lorentz image |
|---|---:|---:|---:|
| `(1,5)` | 1 | 3 | 3 |
| `(2,5)` | 4 | 0 | 6 |
| `(3,11)` | 4 | 0 | 6 |

The first verifier incorrectly demanded that the separate coordinate ranks
called rotations and boosts remain unchanged under an observer boost.  It
returned `12/13`.  The result was not overwritten.  The correction instead
verified exact Poincare intertwiners and the covariant stabilizer

```text
so(3)_n = {A : A n=0}.
```

Under the rational boost, the static coordinate split changed from
`[3,0]` to `[3,2]`, while the total three-dimensional subalgebra transformed
exactly by conjugation and fixed the boosted normal.  This identifies the
original failure as a noninvariant classifier.

The independent audit used an irregular rational tetrahedron, sixteen
unconstrained matrix entries for `A`, four translations and ten separately
imposed Lorentz equations.  It compared the resulting `20`-parameter kernel
directly with the polynomial squared-length kernel.  At three representatives
the two displacement spaces agreed exactly.  Its symbolic determinant was

```text
det(T) = 800 tau (lambda-1)^3.
```

A different rational boost and the convention change `eta -> -eta` left all
invariant decisions unchanged.

## What the result means

The earlier six-shape count now has a precise kinematic explanation.  A
tetrahedron with fixed intrinsic metric is a rigid body in four-dimensional
Minkowski space.  Its relative placement has ten Poincare parameters, and
the four struts impose four independent conditions, leaving six.

For a changing scale, the four strut vectors span spacetime.  They fix the
translation required by any infinitesimal Lorentz motion, so the remaining
kernel can be charted by `so(3,1)`.  At equal scale, all four strut vectors
are parallel; the chart degenerates and three relative spatial translations
replace the three boost directions.

This degeneration is a failure of the Lorentz chart, not a divergence of the
full six-dimensional kernel.  Calling the six variables “extrinsic
curvature” before gluing would still be an invention.

## Post-result literature check

The technical terms learned from the calculation led to three relevant
existing frameworks:

- Poincare 2-group/BFCG formulations use Lorentz connection and tetrad-like
  translational data together, rather than treating a six-dimensional count
  as a Lorentz connection:
  [Mikovic--Vojinovic, 2011](https://arxiv.org/abs/1110.4694) and
  [Asante--Dittrich--Girelli--Riello--Tsimiklis,
  2019](https://arxiv.org/abs/1908.05970).
- Connection phase spaces are larger than length-Regge phase space until
  gluing/metricity constraints are imposed:
  [Dittrich--Ryan, 2008](https://arxiv.org/abs/0807.2806).
- Lorentzian secondary simplicity constraints can recover shape matching and
  discrete extrinsic geometry only on shell:
  [Anza--Speziale, 2014](https://arxiv.org/abs/1409.0836).

No located primary source states the exact four-strut stratification theorem
above.  That absence is not a novelty proof.  External novelty is **OPEN**;
the algebra itself is elementary enough that independent rediscovery is
likely.

## Status ledger

| Claim | Status |
|---|---|
| Six cellular flexes are only a dimension coincidence | **REFUTED** |
| Flexes equal the constrained relative-Poincare kernel | **DERIVED EXACT LOCAL** |
| Expanding stratum is a graph over `so(3,1)` | **DERIVED EXACT LOCAL, adversarially corroborated** |
| Static stratum is rotations plus spatial translations | **DERIVED EXACT LOCAL, adversarially corroborated** |
| Uniform Lorentz-frame variable through `lambda=1` | **REFUTED DERIVED EXACT** |
| Six flexes are extrinsic curvature | **OPEN** |
| Six flexes define a spin/Levi-Civita connection | **OPEN** |
| Face gluing and closure select the missing shapes | **OPEN** |
| First-order action and equivalence with the homogeneous Regge tick | **OPEN** |
| Propagation, continuum limit, `c`, `G` or Planck units | **NOT TESTED** |

## Next discriminating step

Use two adjacent tetrahedra sharing a triangular face.  Transport both local
Poincare kernels into a common face frame and impose equality of the three
shared upper-vertex displacements.  The preregistered alternatives must be:

1. gluing leaves exactly the expected one-dimensional pointwise stabilizer
   of the shared triangle, giving a candidate face holonomy/dihedral mode;
2. gluing removes all relative freedom, so no independent connection is
   present;
3. more freedom survives, so closure/shape matching is still insufficient.

Only after the two-tetrahedron test may the construction be propagated to
all 600 tetrahedra of the spatial carrier.
