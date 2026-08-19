# Protocol: unequal-scale extension of the 600-cell prism shift

Date: 2026-08-19

Prior-art/framing gate commit: `e455921`.

This protocol is frozen before evaluating the exact branch ideal, tangent
cone or complete-carrier consequence.

## 1. Local rational realization

Use the centered edge-`sqrt(8)` regular tetrahedron

```text
b0=( 1, 1, 1),  b1=( 1,-1,-1),
b2=(-1, 1,-1),  b3=(-1,-1, 1).
```

All four radii equal `3`.  Let

```text
a=q-1,
s=(s1,s2,s3),
B_i=b_i,
T_i=q*b_i+s+N*n,
n.n=-1,  n.b_i=n.s=0.
```

For each corresponding strut define the exact Lorentzian square

```text
ell_i^2 = |a*b_i+s|^2-N^2.
```

For every one of the six lateral quadrilaterals, verify exact affine
planarity from the four declared vertices.  This guards against proving a
statement about a graph framework that is not a polytopal prism/frustum.

## 2. Common-strut ideal

Form the three independent differences

```text
c_i=ell_i^2-ell_0^2,  i=1,2,3.
```

Let `E` have rows `(b_i-b_0)^T`.  Require, as exact polynomial identities,

```text
det(E)=-16,
c=2*a*E*s,
E^-1*c/2=a*s.
```

Consequently the common-strut ideal is, after an invertible rational row
operation,

```text
I=<a*s1,a*s2,a*s3>
 = <a> intersection <s1,s2,s3>.
```

Verify both ideal containments by explicit polynomial reductions.  The real
solution set must therefore be exactly

```text
q=1  union  s=0.
```

No numerical root finder or tolerance is permitted in this decisive gate.

## 3. Singular anchor and tangent cone

At `(a,s)=(0,0)`, compute the exact `3 x 4` constraint Jacobian.  It must be
zero, so the Zariski tangent space has dimension four.

For a tangent direction `(alpha,beta)` compute the second directional
derivative:

```text
d^2 c/dt^2 at t=0 = 4*alpha*E*beta.
```

Require:

- the pure-scale direction and all three pure-shift directions lie in the
  quadratic tangent cone;
- every frozen mixed direction `(alpha=1,beta=e_i)` is rejected;
- invertibility of `E` proves the cone is exactly the union
  `alpha=0` or `beta=0`.

Thus the large linear tangent at the static point may not be promoted to a
smooth four-dimensional configuration chart.

## 4. Global 600-cell consequence

Reconstruct the literal 600-cell and its tetrahedra.  Verify exact counts

```text
(V,E,F,T)=(120,720,1200,600)
```

and graph connectedness.  On a homogeneous unequal-scale slab, apply the
local common-strut result to all 600 cells.  Each local tangential shift must
vanish.  Under the already derived potential gluing map this gives

```text
d0 phi=0.
```

Compute the exact incidence rank independently over two primes and require
rank `119`; hence `phi` is constant and the quotient shift carrier has
dimension zero for `q!=1`.

The equal-scale positive control must retain rank `119` modulo constants.

## 5. Nonuniform-strut recovery and controls

For arbitrary strut differences `c` and `a!=0`, reconstruct

```text
s=E^-1*c/(2*a).
```

Use the exact frozen controls

```text
a in {-1/2, 1/3, 1},
s in {(1,2,3), (-2,1,4), (3,-1,-2)}.
```

Generate `c` from the geometry and require exact recovery.  Check

```text
det(2*a*E)=(2*a)^3 det(E),
```

so the inverse necessarily has a pole at equal scale.

Negative geometric control: replace only

```text
b3=(-1,-1,1) by b3'=(-1,-1,2),
```

while keeping homothety about the declared origin.  At `a=1`, solve the
common-strut equations and require a nonzero `s` that equalizes all four
struts.  This demonstrates that equal radii/concentric regularity, rather
than the implementation, forces zero shift.

## 6. Verdicts

Return

```text
DYNAMIC_SHIFT_EXTENSION_OBSTRUCTED
```

only if all exact local, ideal, tangent-cone, global and negative-control
gates pass.

Return `UNEQUAL_SCALE_SHIFT_SURVIVES` if a nonzero common-strut shift exists
for the regular control at `q!=1`.

Return `BRANCH_GEOMETRY_OPEN` on any failure of planarity, exact ideal
reduction, incidence rank or negative-control sensitivity.

## 7. Interpretation boundary

A passing obstruction means:

- **DERIVED:** the equal-scale potential sector is a singular static branch;
- **DERIVED:** it cannot be copied as an independent field onto successive
  homogeneous unequal-scale ticks with one common strut/lapse;
- **STRUCTURAL:** nonuniform strut data can encode a derived shift away from
  equal scale;
- **OPEN:** the status of that derived variable after the complete internal
  Regge equations and constraints are solved.

It does not refute the already derived static Hessian.  It refutes only its
premature interpretation as the spatial half of a freely propagated scalar
field.  No dispersion, `c`, mass or Planck target is used.

Only this mission verifier and static guards may run; the full suite is
excluded.
