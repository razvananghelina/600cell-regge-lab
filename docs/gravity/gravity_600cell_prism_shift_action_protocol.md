# Protocol: cellular Regge Hessian on the prism-shift family

Date: 2026-08-19

Prior-art gate commit: `29fac99`.

This protocol is frozen before evaluating a nonconstant shift potential,
assembling the Hessian or comparing it with the graph Laplacian.

## 1. Frozen carrier and normalization

Reconstruct the regular 600-cell boundary from `commons/cell600.py`.  Rescale
its Euclidean chord coordinates so every spatial edge has length

```text
L=1.
```

Use

```text
rho=1
```

for direct action controls.  The analytic Hessian must retain symbolic
positive `L,rho`.

For each tetrahedron sorted as `(v0,v1,v2,v3)`, use affine prism coordinates
`(lambda1,lambda2,lambda3,z)` and metric

```text
H=[G a; a^T -rho],
a_i=phi(vi)-phi(v0).
```

The inward facet conormals are frozen as

```text
bottom=(0,0,0,+1),       top=(0,0,0,-1),
n0=(-1,-1,-1,0),         n1=(1,0,0,0),
n2=(0,1,0,0),            n3=(0,0,1,0).
```

For two facets use

```text
cos(theta)=-n_i^T H^-1 n_j
            /sqrt[(n_i^T H^-1 n_i)(n_j^T H^-1 n_j)].
```

On a negative real square-root argument use the already certified
Lorentzian branch

```text
sqrt(-x-i0)=-i sqrt(x).
```

At `phi=0`, every lateral angle must be `acos(1/3)` and every bottom/top
angle `pi/2`.

## 2. Complete action

For every spatial edge `e=(u,v)`, let

```text
x_e=phi(v)-phi(u),
A_e=i*sqrt(rho+x_e^2)
```

in the direct `L=rho=1` control.  Sum the six appropriate lateral-facet
angles from each of the five tetrahedra incident on `e`:

```text
epsilon_e=2*pi-sum_(sigma contains e) theta_(sigma,e).
```

For each of the 1,200 spatial triangles, sum its two incident bottom angles
and its two incident top angles:

```text
psi_bottom=pi-sum theta_bottom,
psi_top   =pi-sum theta_top.
```

With the common equilateral triangle area `sqrt(3)/4`, evaluate

```text
S_grav=-i [sum_edges A_e epsilon_e
           +sum_faces A_triangle(psi_bottom+psi_top)].
```

The 120-strut dust action is constant in `phi` and is recorded but omitted
from derivatives.

## 3. Mandatory branch and boundary controls

Before interpreting a Hessian:

1. recover the 600-cell f-vector and five tetrahedra per spatial edge;
2. certify signature `(3,1)` for every cell at every direct control point;
3. recover the static lateral and boundary angle anchors;
4. verify cell by cell that

```text
theta_top=pi-theta_bottom
```

on corresponding facets to absolute tolerance `2e-11`;
5. require the complete top-plus-bottom boundary action below `2e-9` in
   absolute value at every control point;
6. require imaginary contamination of the total gravitational action below
   `2e-9`;
7. reproduce the exact static value

```text
720*[2*pi-5*acos(1/3)]
```

to relative tolerance `2e-12`.

The direct nonconstant controls are the five frozen integer directions:

```text
one vertex delta,
graph distance from that vertex,
squared graph distance,
first ambient coordinate rank,
(17*i^2+3*i+5) mod 101,
```

each centered by its mean and normalized to Euclidean norm one.  Evaluate
them at amplitudes `0`, `+1e-3` and `-1e-3`.

## 4. Schläfli-gradient control

At the base potential

```text
phi_base=1e-3*(normalized squared graph distance),
```

compute the complete lateral deficit angles.  The preregistered reduced
gradient is

```text
g_v=sum_(oriented e incident v)
    sign(v,e)*epsilon_e*x_e/sqrt(1+x_e^2).
```

Compare its directional contractions with centered differences of the full
area-angle action along the other four normalized directions, using steps

```text
1e-5 and 5e-6.
```

The maximum relative discrepancy must be below `3e-7`, with an absolute
floor `3e-9`.  Failure blocks the Hessian interpretation.

## 5. Exact Hessian prediction

Independently differentiate symbolically

```text
sqrt(rho*L^2+x^2)
```

twice at `x=0`.  It must give `1/(L*sqrt(rho))`.

Using only the source edge incidence, assemble

```text
H_pred=kappa*d0^T*d0,
kappa=[2*pi-5*acos(1/3)]/(L*sqrt(rho)).
```

At `L=rho=1`, compare directional second differences of the complete action
for all five frozen directions at

```text
h=1e-3, 5e-4, 2.5e-4.
```

Use the two finest values for the second-order Richardson extrapolate.  The
maximum relative discrepancy from `d^T H_pred d` must be below `2e-7`, with
absolute floor `2e-9`.

## 6. Exact spectral census

Verify entrywise that `d0^T*d0=12I-A_600`.  Compute the symmetric spectrum
and compare all nine clusters, without relabelling, with the already frozen
exact values and multiplicities from the prior-art note.  After removing the
constant vector require:

```text
rank=119,
minimum eigenvalue=kappa*(12-6*varphi)>0.
```

No fitted coefficient, diagonal shift or spectral reordering is allowed.

## 7. Verdicts frozen before execution

Return

```text
SHIFT_HESSIAN_IS_GRAPH_LAPLACIAN
```

only if every carrier, branch, boundary, Schläfli, direct-Hessian and exact-
spectrum gate passes.

Return `SHIFT_DIRECTIONS_EXACTLY_NULL` only if all complete-action controls
pass but the quotient Hessian is zero.

Return `PREDICTION_REFUTED_NONZERO_OTHER_OPERATOR` if controls pass and a
different nonzero Hessian is resolved.  Return `ACTION_BRANCH_OPEN` if any
angle, boundary or Schläfli control fails.

## 8. Classification boundary

A graph-Laplacian result is **DERIVED IN THIS RESTRICTED SECTOR**.  It means
the cellular Regge action locally selects constant shift potential.  It does
not make this scalar-potential sector a graviton, a matter field or the full
ADM shift, and it supplies no propagation law or limiting speed.

Only the registered mission verifier and static guards may run.  The full
suite is excluded.

