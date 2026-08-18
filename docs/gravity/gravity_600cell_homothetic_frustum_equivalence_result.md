# Homothetic staircase schedules are one Lorentzian 4-frustum

Date: 2026-08-17

## 1. Provenance

- prior-art gate: `e4428bf`;
- frozen protocol: `affb02e`;
- registered implementation before first evaluation: `0edd49d`;
- targeted verifier:
  `reproducible/verify_gravity_600cell_homothetic_frustum_equivalence.py`;
- result artifact:
  `reproducible/gravity_600cell_homothetic_frustum_equivalence.json`;
- artifact SHA-256:
  `7e7c23efaf24a2c99a68f3b302b9ef575e0f777ef46f73ccaea9f99e1ecd58dc`.

Only the targeted verifier was run.  It returns **9/9**.  The full suite was
not run.

## 2. Mechanical verdict

```text
HOMOTHETIC_SCHEDULES_ONE_LORENTZIAN_FRUSTUM
```

**DERIVED.**  For every

```text
R_minus>0, R_plus>0, rho>0,
```

all 24 staircase schedules are nonoverlapping triangulations of the same flat
tetrahedral Lorentzian 4-frustum.  Their pullback metric components differ at
the same parameter point, but the metrics are exactly isometric through PL
changes of coordinates that fix both spatial time boundaries pointwise.

This refutes the physical no-go interpretation of the immediately preceding
`gcd=r-1` result.  It does not refute its arithmetic.

## 3. One physical affine four-plane

Let the regular tetrahedron directions satisfy

```text
u_i.u_i=1,
u_i.u_j=c=phi/2  (i!=j),
D=1+3c.
```

For

```text
p_i^-=(R_minus u_i,0),
p_i^+=(R_plus u_i,T),
T^2=rho+(R_plus-R_minus)^2,
```

define

```text
n=(u_0+u_1+u_2+u_3)/D.
```

The exact Gram calculation gives

```text
n.u_i=1,
n.n=4/D.
```

All eight vertices therefore satisfy the single affine equation

```text
n.x_space-[(R_plus-R_minus)/T] x_time=R_minus.
```

An independent homogeneous-coordinate calculation gives rank five and
nullity one, hence an affine hull of dimension four and no second hidden
constraint.

The Minkowski normal has square

```text
4/D-(R_plus-R_minus)^2/[rho+(R_plus-R_minus)^2]
 = [4/D-1]+rho/[rho+(R_plus-R_minus)^2] > 0,

4/D-1=3(1-c)/D>0.
```

The normal is spacelike, so the common four-plane has Lorentzian signature.

## 4. Independent intrinsic signature

Use coordinates `(y_0,y_1,y_2,z)` on the affine plane, with

```text
sum_i y_i=R_minus+(R_plus-R_minus)z,
x_time=T z.
```

The exact induced `4 x 4` Gram matrix has spatial leading block

```text
(1-c)(I_3+J_3),
det=4(1-c)^3>0,
```

and its final Schur complement is

```text
-rho-(3/4)(1-c)(R_plus-R_minus)^2 < 0.
```

Thus the signature is **exactly `(3,1)`** for all admitted parameters.  No
sampled eigenvalue or numerical threshold was used.

## 5. The common polytope

In these coordinates the physical cell is

```text
Q={0<=z<=1,
   y_i>=0,
   sum_i y_i=R_minus+(R_plus-R_minus)z}.
```

Exhausting intersections of its six facet hyperplanes gives exactly:

```text
8 vertices,
2 tetrahedral time facets with 4 vertices each,
4 lateral triangular-frustum facets with 6 vertices each.
```

Its face lattice is that of `Delta^3 x I`, and its exact coordinate
four-volume is

```text
Vol(Q)
 = (R_minus^3+R_minus^2 R_plus
    +R_minus R_plus^2+R_plus^3)/24
 = (R_minus+R_plus)(R_minus^2+R_plus^2)/24.
```

## 6. Every staircase is a geometric triangulation

All 96 labelled `(order,split)` simplices were checked.  The physical-to-
parameter oriented determinant ratio for split `k` is exactly

```text
R_minus^k R_plus^(3-k)>0.
```

Therefore no simplex degenerates or reverses orientation.  For each of the 24
schedules:

```text
internal tetrahedral facets : 3, each used twice
boundary tetrahedra         : 14
  bottom                    : 1
  top                       : 1
  each of 4 lateral facets  : 3
sum of 4 simplex volumes    : Vol(Q).
```

Coherent orientation, boundary and volume still do not alone exclude a fold.
The verifier closes that gap with the explicit projective homeomorphism from
the standard prism:

```text
d=(1-tau)/R_minus+tau/R_plus,
y_i=lambda_i/d,
z=(tau/R_plus)/d,
```

whose exact inverse is

```text
lambda_i=y_i/[R_minus+(R_plus-R_minus)z],
tau=z R_plus/[R_minus(1-z)+R_plus z].
```

Both denominators are positive throughout the admitted domains.  Both
compositions, all eight vertices and all six facets were verified exactly.
Projective maps preserve convex simplices and their intersections, so every
standard staircase maps to a nonoverlapping triangulation of the same `Q`.

## 7. The correct reconciliation of the metrics

There are 32 distinct local affine map types.  For each pair, with Jacobians
`A` and `B`, the verifier forms

```text
H=B^-1 A
```

and checks

```text
H^T (B^T G_Q B) H = A^T G_Q A.
```

All `32 x 32` identities pass exactly.  The bottom and top restrictions are
independent of schedule:

```text
bottom: lambda -> R_minus sum_i lambda_i u_i,
top:    lambda -> R_plus  sum_i lambda_i u_i.
```

Hence the schedule transition is an isometry and fixes both physical Cauchy
boundaries pointwise.

The previous audit instead demanded

```text
A^T G_Q A = B^T G_Q B
```

without the coordinate Jacobian `H`.  Its exact obstruction

```text
gcd=r-1
```

says that the identity map of the chosen parameter prism is an isometry only
in the static case.  That is a coordinate statement, not a physical
distinction between the frusta.

## 8. Edge-length reconciliation

The common frustum independently reproduces every inherited homothetic length:

```text
bottom spatial edge^2 = L_minus^2,
top spatial edge^2    = L_plus^2,
same-vertex strut^2   = -rho,
cross diagonal^2      = L_minus L_plus-rho,

R_minus=phi L_minus,
R_plus =phi L_plus.
```

Thus this is not a different geometry invented to rescue the result.  It is
the geometry already encoded by the repository's Regge simplices.

## 9. Consequences for the programme

### Homogeneous sector

**DERIVED GEOMETRIC POSITIVE.**  Schedule ambiguity is gauge on the regular
homothetic sector.  The four already accepted homogeneous dust ticks are not
invalidated by a metric-choice ambiguity between the two independently
derived parities.  Their observed parity agreement now has an exact geometric
explanation.

This does not yet prove equality of every action and canonical boundary term;
that is the next explicit comparison.

### Universal overlay

The 148-chamber overlay is unnecessary for identifying the homogeneous
physical cell.  It is a common refinement only after choosing a common
parameter identification, and that identification is not physical data.

The simpler product face poset returns in a different role: because `Q` has
face lattice `Delta^3 x I`, its order complex is a canonical barycentric
triangulation of the physical frustum.  It need not refine every old
staircase to be a legitimate schedule-free carrier.  Whether its simplicial
Regge action equals a direct cellular-frustum action remains **OPEN**.

### Anisotropic sector

No conclusion extends automatically to arbitrary independent boundary edge
data.  Homothety is what placed all eight vertices in one flat four-plane.
Perturbations may violate flat-frustum embeddability or shape matching, making
the schedule choice physical again.  Anisotropic stability and gravitational-
wave modes remain **OPEN**.

## 10. Prior-art correction

The physical cell itself is **KNOWN**, not a new discovery.  Tsuda and
Fujiwara explicitly use the tetrahedral frustum as the fundamental block of a
regular 4-polytopal universe, including the 600-cell case:
<https://arxiv.org/abs/2011.04120>.  De Felice and Fabri supply the inherited
five-dimensional Lorentzian embedding:
<https://arxiv.org/abs/gr-qc/0009093>.  Modern Lorentzian frustum cosmology and
its causal gates are developed by Jercher and Steinhaus:
<https://arxiv.org/abs/2312.11639>.

The exact local projective reconciliation of this repository's 24 schedules
and the correction of its `r-1` interpretation were not located explicitly.
External novelty of that narrow statement remains **OPEN**.

## 11. Next falsification test

Do not enumerate another carrier first.  Compare the **homogeneous action**:

1. derive the direct cellular tetrahedral-frustum Regge action with its
   trapezoidal and spatial triangular hinges;
2. substitute the same `(L_minus,L_plus,rho,M)` into both committed staircase
   actions;
3. compare the complete action, internal variation and pre/post boundary
   momenta exactly or at certified arbitrary precision;
4. require equality for every positive Lorentzian point in a stated symbolic
   domain, not only at the four accepted ticks.

If the actions agree, homogeneous dynamics is schedule-independent and the
frustum is the correct canonical carrier.  If they differ, the common metric
is real but the simplicial action contains a subdivision artifact that must be
removed by a cellular/perfect action before the ticks can be called canonical.
