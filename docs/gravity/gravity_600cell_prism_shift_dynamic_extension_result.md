# Result: the equal-scale prism shift is a singular static branch

Date: 2026-08-19

## Verdict

**DERIVED EXACT KINEMATIC OBSTRUCTION.**  For a regular tetrahedral frustum
with homogeneous scale ratio `q`, one common timelike strut square and a
tangential translation `s` of the top tetrahedron, the complete common-strut
condition is exactly

```text
(q-1)*s=0.
```

Consequently:

```text
q=1:       the three local shift components survive;
q!=1:      the local shift is exactly zero.
```

On the connected 600-cell, the equal-scale family has 119 vertex-potential
directions modulo constants, while a homogeneous unequal-scale common-strut
slab has zero such directions.

The positive Regge Hessian

```text
H_phi = [2*pi-5*acos(1/3)]/(L*sqrt(rho)) * Delta_0
```

remains valid on its declared equal-scale branch.  What fails is the proposed
interpretation of `phi` as an independent field that can simply be copied
onto successive unequal-scale ticks.

## Provenance

| stage | commit |
|---|---|
| prior-art and framing gate | `e455921` |
| frozen protocol | `b251d5b` |
| registered verifier | `2031d85` |
| preserved first `12/13` artifact | `fcb3a14` |
| disclosed comparator correction | `7380dcf` |
| exact entrywise repair | `b6cffb4` |
| passing `13/13` artifact | `25a5070` |

The final artifact is

```text
reproducible/gravity_600cell_prism_shift_dynamic_extension.json
SHA-256 32d5269b27756a4c6fec4603855db643106e571007d3f3dd1a0a6c69d33a0095.
```

The first artifact remains in history with SHA-256

```text
134f8e68335f3acdd40eb909d4dcb4fae4361329c4540a65182fe061affb499a.
```

Its only failed gate compared mutable and immutable SymPy matrices by
structural container equality.  The printed expressions and the independent
transformed identity were already exact.  The preregistered correction
replaced that comparator by three exact simplified entrywise equalities; no
tolerance or verdict was changed.

## Exact local theorem

Use the centered rational regular tetrahedron

```text
b0=( 1, 1, 1),  b1=( 1,-1,-1),
b2=(-1, 1,-1),  b3=(-1,-1, 1).
```

All four squared radii are `3`.  Put `a=q-1` and realize the corresponding
vertices as

```text
B_i=b_i,
T_i=q*b_i+s+N*n,
n.n=-1.
```

Every lateral quadrilateral is exactly planar: its three affine difference
vectors satisfy

```text
T_j-B_i = (T_i-B_i)+q*(B_j-B_i).
```

Thus the calculation concerns a genuine frustum rather than an incomplete
edge framework.

The four strut squares are

```text
ell_i^2=|a*b_i+s|^2-N^2.
```

Let `E` have rows `(b_i-b_0)^T` for `i=1,2,3`.  The exact audit obtains

```text
det(E)=-16,
c_i=ell_i^2-ell_0^2,
c=2*a*E*s,
E^-1*c/2=a*s.
```

Therefore equal struts give `a*s=0`.  Exact Groebner elimination also gives

```text
<a*s1,a*s2,a*s3>
  = <a> intersection <s1,s2,s3>.
```

The admissible algebraic set is precisely the union of the equal-scale branch
and the zero-shift branch.

## Why the linearized picture is misleading

At `(a,s)=(0,0)` the constraint Jacobian is the zero `3 x 4` matrix.  Its
linear tangent space therefore contains one scale direction and three shift
directions.

For a tangent `(alpha,beta)`, however, the exact second directional
derivative is

```text
d^2 c/dt^2 = 4*alpha*E*beta.
```

Since `E` is invertible, the quadratic tangent cone is only

```text
alpha=0  union  beta=0.
```

Every pure scale or pure shift direction integrates along its own branch,
but a generic scale-plus-shift direction does not.  This is the precise
reason a static Hessian cannot automatically be completed by inventing a
mixed temporal block at the singular anchor.

## Complete-carrier consequence

The verifier independently reconstructs

```text
(V,E,F,T)=(120,720,1200,600),
```

five tetrahedra per edge and a connected vertex graph.  The oriented
vertex-edge incidence has rank `119` over both `F_101` and `F_1000003`.

At equal scale, gradients of arbitrary vertex potentials supply the 119
shape-matched shifts.  At unequal homogeneous scale, every local shift
vanishes, so every potential difference vanishes.  Connectivity then forces
the potential to be constant.  Hence the quotient carrier collapses from
dimension `119` to `0`.

## What nonuniform struts do

If the four strut squares are allowed to differ, their three independent
differences satisfy

```text
c=2*(q-1)*E*s,
s=E^-1*c/[2*(q-1)].
```

All three frozen rational controls recover `s` exactly.  The map determinant
is

```text
det(2*(q-1)*E)=-128*(q-1)^3,
```

so the inverse necessarily becomes singular at equal scale.

Thus away from equal scale a tangential shift can reappear only as a variable
encoded by inhomogeneous strut/lapse data; it is not an independent invisible
cell modulus.

The geometric negative control changes one tetrahedron radius while retaining
homothety about the declared origin.  It finds the nonzero compensating shift

```text
s=(3/10,3/10,-3/10)
```

and makes all four strut squares equal.  The code can therefore detect a
surviving shift when the equal-radius hypothesis is removed.

## Physical reading

- **DERIVED:** the 119-dimensional `phi` carrier is confined to the
  equal-scale branch of the common-lapse homogeneous geometry.
- **DERIVED:** the static graph-Laplacian Hessian is not by itself the spatial
  half of a scalar wave equation.
- **DERIVED:** at unequal scale, nonuniform strut differences determine the
  corresponding local shifts with a `1/(q-1)` singular map.
- **STRUCTURAL:** the variables look more like lapse/shift constraint data
  than propagating tensor modes.
- **OPEN:** whether eliminating the complete internal Regge equations leaves
  any physical longitudinal response on the boundary.
- **OPEN:** tensor propagation, a limiting speed, continuum convergence,
  inertia and masses.

This is not a failure of the theory.  It prevents a false shortcut from a
static spatial Hessian to a d'Alembertian.

## Post-result literature reconciliation

The learned-term search again found the standard neighbouring constructions:

- homogeneous common-strut Collins--Williams models and their local/global
  variation issue in
  [Liu--Williams](https://arxiv.org/abs/1501.07614);
- regular/homothetic polytopal frusta in
  [Tsuda--Fujiwara](https://arxiv.org/abs/2011.04120);
- background-dependent lapse/shift and curvature modes in canonical
  linearized Regge calculus in
  [Hoehn](https://arxiv.org/abs/1411.5672).

No located primary source states the exact branch ideal above.  Search
absence is not a novelty proof, so external novelty remains **OPEN**.

## Corrected next step

Do not construct a free recurrence for one copy of `phi` per homogeneous
tick.  The next admissible test is to embed the exact nonuniform-strut map
into the already certified complete `840`-internal-edge slab carrier:

1. use vertex strut-square differences to reconstruct the 119 longitudinal
   local shifts at the first accepted unequal-scale background;
2. form this canonical 119-dimensional subspace before reading the action
   Hessian;
3. restrict the complete internal equations and their boundary coupling to
   it;
4. determine whether it is eliminated as an auxiliary constraint sector or
   induces a nonzero Schur-complement response on boundary perturbations;
5. only then compare any surviving boundary operator with `Delta_0`.

This test can genuinely distinguish a constraint/lapse sector from a
propagating one.  A desired wave speed or spectrum must not be loaded.
