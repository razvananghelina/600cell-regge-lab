# Prior-art gate: variable face connection in homothetic-frustum gluing

Date: 2026-08-19

## Exact object and hypotheses

This mission audits one assumption in the fixed-frame two-frustum and global
holonomy results before any global reconstruction Jacobian is built.

Use two nondegenerate homothetic tetrahedral Lorentzian frusta sharing one
lateral triangular-frustum face.  Their lower tetrahedra are fixed, their
upper vertices are

```text
q_i = lambda p_i + tau n,  tau != 0,
```

and each cell is restricted to the exact six-dimensional kernel preserving
its six upper-edge and four strut squared lengths to first order.

The earlier fixed-frame gluing imposed equality of the three shared upper
vertex displacements after the background lower-slice development was fixed.
The present object instead permits the face transition itself to vary by an
infinitesimal Poincare isometry `X` which fixes the three shared lower
vertices pointwise.  Thus

```text
X p_i = 0                                      (shared lower triangle),
delta q_i(left) = delta q_i(right) + X q_i     (shared upper triangle).
```

No arbitrary Schur or connection coefficient is permitted.  The admissible
space for `X` is exactly the pointwise Poincare stabilizer of the lower
triangle, derived from the same carrier.  It must have dimension one.

The question is whether allowing this derived variation changes the earlier
diagonal-only compatible-pair space.

## Why this gate is necessary

The complete 600-cell holonomy audit correctly computed the common fixed
space of the local flex seed under a **frozen** background connection.  A
metric deformation, however, can also vary the Levi-Civita face transition.
Holding that variation at zero is an additional hypothesis unless the gluing
equations themselves force it.

A zero parallel section for a frozen connection is therefore not by itself a
proof that the full metric reconstruction has zero infinitesimal kernel.  The
two statements coincide only if the variable-transition test below reduces
again to the fixed-frame diagonal.

## Primary literature

### KNOWN

1. In connection formulations of Regge calculus, matrices on codimension-one
   simplices are discrete connection variables; the curvature around a hinge
   is their ordered product.  Eliminating the connection by its equations
   recovers length Regge calculus:
   [Khatsymovsky, *Affine connection form of Regge calculus*,
   2015](https://arxiv.org/abs/1509.04974).

2. Edge lengths and piecewise-affine coordinates determine the discrete
   metric-compatible connection; consequently the connection changes when
   the metric data change:
   [Khatsymovsky, *On the discrete Christoffel symbols*,
   2019](https://arxiv.org/abs/1906.11805).

3. Connection/area-angle phase spaces are larger than length-Regge phase
   space before gluing and metricity constraints are imposed:
   [Dittrich--Ryan, *Phase space descriptions for simplicial 4d
   geometries*, 2008](https://arxiv.org/abs/0807.2806).

4. In Lorentzian and general cellular settings, secondary simplicity and
   Levi-Civita conditions are tied to shape matching; a face connection
   cannot in general be frozen independently of those conditions:
   [Anza--Speziale, *A note on the secondary simplicity constraints in loop
   quantum gravity*, 2014](https://arxiv.org/abs/1409.0836).

### CONTROL

- Setting `X=0` must reproduce the already certified fixed-frame result:
  compatible-pair dimension six and zero relative difference.
- Before the local strut restriction, the pointwise stabilizer of either
  nondegenerate shared triangle must be exactly one-dimensional.
- A simultaneous Poincare change of both cell frames is gauge and must not be
  counted as a relative physical mode.

### OPEN

- Whether the derived one-dimensional lower-face stabilizer has a nonzero
  action on the shared upper triangle inside the difference of the two local
  six-flex kernels.
- Whether the variable-transition compatible-pair dimension is six, seven or
  larger.
- If a relative mode appears, whether global edge-star closure constrains it
  or leaves a genuine metric deformation.
- Whether the connection is uniquely Levi-Civita for the nonsimplicial
  cellular frustum data without additional shape variables.

## Proposed difference from prior art

The standard frameworks establish why the connection must be included or
eliminated consistently.  The proposed calculation is the exact rank and
intersection theorem for this specific homothetic tetrahedral-frustum carrier
and its already certified six-dimensional local kernel.

No located primary source states that rank.  Search absence is not a novelty
proof; external novelty remains **OPEN**.

## Decision boundary

- If the variable-transition compatible space is still exactly diagonal,
  the frozen-connection propagation is locally justified and the global
  reconstruction regularity test remains the next step.
- If a derived relative mode appears, the broad global-rigidity
  interpretation is refuted.  The fixed-connection holonomy calculation
  remains valid only conditionally, and global closure must be rebuilt with
  connection variations included.
- If more than the derived face-stabilizer freedom appears, the cellular
  gluing data are further underdetermined and the schedule-free metric route
  remains unclosed.

Only a preregistered targeted verifier and static registry guards may be run.
