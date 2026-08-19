# Prior-art gate: the 119-dimensional static kernel as vertex gradients

Date: 2026-08-19

## Exact object

The complete variable-face closure census found, before interpretation,

```text
lambda=1: modular rank 3481 in 3600 columns, deficit 119;
lambda=2,3: modular rank 3600, deficit 0.
```

Both disclosed primes gave the same ranks and every construction control
passed.  A modular deficit is not yet a rational kernel theorem.

At `lambda=1`, the local six-dimensional kernel is

```text
three spatial rotations + three spatial translations.
```

The one relative face mode derived in the two-cell theorem is a translation
normal to the shared spatial triangle.  This suggests the target-blind map

```text
continuous piecewise-linear vertex scalar
    -> its constant spatial gradient on each tetrahedron
    -> the translation sector of the local six-kernel.
```

The domain has 120 vertex values.  On a connected triangulation the only
piecewise-linear scalar with zero gradient on every tetrahedron is a global
constant, so the gradient image has dimension `120-1=119`.

The mission is to prove or refute

```text
ker(C_static) = grad(P1)  and dim = 119
```

over the rationals, without reading a physical interpretation into the
scalar.

## Primary literature

### KNOWN

Continuous piecewise-linear finite elements on a simplicial mesh are
determined by one value per vertex.  Their gradients are piecewise constant;
continuity across a face forces the tangential gradient components to agree,
so the gradient jump is normal to that face.

This is the first map in the finite-element de Rham complex.  Relevant
primary references include:

- [Arnold--Falk--Winther, *Finite element exterior calculus: from Hodge
  theory to numerical stability*, 2009](https://arxiv.org/abs/0906.4325);
- [Arnold--Falk--Winther, *Finite element exterior calculus, homological
  techniques, and applications*, 2006](https://doi.org/10.1017/S0962492906210018);
- [Nedelec, *Mixed finite elements in R3*,
  1980](https://doi.org/10.1007/BF01396415).

The generic body--hinge theorem explains why normal-jump compatibility is a
rigidity constraint, but it does not identify this special symmetric kernel.

### CONTROL

1. Reconstruct the exact static reduced matrix `C_static` from the frozen
   protocol and reproduce both modular ranks `3481`.
2. Build the local gradient map directly from the four canonical tetrahedron
   vertices; do not fit a null vector to `C_static`.
3. Assemble the global `3600 x 120` map `G`, with zero rotation coordinates
   and the local gradient in the translation coordinates.
4. Verify every face identity `C_f G_f=0` exactly over the rationals.
5. Verify `G 1=0` and exact rank 119.  A nonzero minor modulo each disclosed
   prime plus the explicit constant kernel supplies the rational rank proof.
6. An intentionally discontinuous one-vertex value assignment evaluated with
   the wrong neighboring vertex correspondence must fail at least one face;
   this guards against a verifier which annihilates every cellwise gradient.

### OPEN

- Whether `grad(P1)` exhausts the rational static kernel.
- Whether the scalar is lapse, time-reparametrization gauge, a bending mode
  or merely a kinematic finite-element potential.
- Why nonzero homothetic expansion lifts the entire gradient space.
- Whether an action assigns zero, positive or negative quadratic weight to
  these modes.

## Exact proof boundary

If

```text
C_static G = 0,
rank_Q(G)=119,
rank_Q(C_static)>=3481,
```

then rank-nullity forces

```text
ker_Q(C_static)=im_Q(G),
dim ker_Q(C_static)=119.
```

The modular rank supplies the lower bound on `rank_Q(C_static)`; the explicit
rational image supplies the matching upper bound.  No rational row reduction
of the full matrix is required.

## Interpretation firewall

Equality with a standard discrete gradient is a mathematical identification,
not a physical gauge theorem.  The name “lapse” is forbidden until an action
or canonical constraint annihilates the same directions.

Only a preregistered targeted verifier and static registry guards may be run.
