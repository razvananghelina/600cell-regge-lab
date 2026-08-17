# Preregistration: barycentric product slab and staircase-refinement census

Date: 2026-08-17

Prior-art commit: `dec110d`.

Status: frozen before constructing any product cell, chain or containment
matrix.

## 1. Frozen source and exclusions

Rebuild the 600-cell boundary only through the audited functions `parity` and
`boundary_600_cell` from

```text
reproducible/verify_dimension_reconciliation.py
SHA-256 a819ae9d472317d456cf7d67f588b31586b7aed2400a2540b8352f4661b39d45.
```

Require the exact f-vector `(120,720,1200,600)` and simplicial incidences.
Do not import a gravity action, schedule output, desired continuum result,
physical constant or fitted coordinate.

## 2. Global product face poset

Represent each cell of `X=K x I` as

```text
(sigma,T),
sigma a nonempty simplex of K,
T in {bottom vertex, top vertex, interval edge}.
```

Order cells by componentwise face inclusion.  Require

```text
number of product cells = 3*(120+720+1200+600)=7920.
```

Enumerate strict chains by dynamic programming over cell dimension.  The
number of chains of length `r+1` is the `r`-simplex count of the order complex
`B`.  Report the complete `f(B)=(f0,...,f4)`.

Before enumeration, two analytic controls are frozen:

```text
f0 = 7920,
f4 = 600 * 8 * 4! = 115200.
```

The second identity uses that each four-dimensional simple polytope
`Delta^3 x I` has eight vertices and `4!` maximal flags at each vertex.
Require dimension four and Euler characteristic

```text
f0-f1+f2-f3+f4=0,
```

as expected for `S^3 x I`.

## 3. Boundary controls

Restrict to cells whose interval component is one of the two endpoints.
Require two disjoint induced boundary components.  Each must have

```text
2640 vertices,
600*4! = 14400 top-dimensional tetrahedra,
Euler characteristic 0.
```

Interval reflection must swap the two components, preserve every vertical
product cell setwise and preserve every chain incidence.  This is a
combinatorial involution; no metric embedding is used.

## 4. Exact local enumeration

For one `Delta^3 x I`, enumerate all

```text
15 nonempty tetrahedron faces x 3 interval cells = 45 product cells
```

and every maximal chain.  Require exactly

```text
8*4! = 192
```

barycentric four-simplices.

Embed the eight prism vertices in rational coordinates
`(lambda_1,lambda_2,lambda_3,t)`.  Put each product-cell barycentre at the
product of ordinary face barycentres.

For every one of the 24 total vertex orders construct the four standard
staircase four-simplices.  Using exact rational affine coordinates, assign a
barycentric four-simplex to a staircase simplex only if all five barycentric
vertices lie in it.  Record for each order:

```text
contained_count,
mixed_count = 192-contained_count,
multiplicity distribution of containing staircase simplices.
```

No epsilon or geometric tolerance is permitted.

## 5. Functorial canonicity

Check mechanically that permuting the four tetrahedron vertices maps product
cells to product cells and maximal chains to maximal chains.  Together with
the interval reflection this realizes the local `S4 x C2` action.  Globally,
the same componentwise construction proves functoriality under every
automorphism of `K`; an enumeration of all 14400 `H4` elements is unnecessary
and would not strengthen the poset proof.

## 6. Mechanical outcome

- all 24 orders have `contained_count=192` and `mixed_count=0`:
  `BARYCENTRIC_PRODUCT_IS_UNIVERSAL_STAIRCASE_REFINEMENT`;
- at least one order has a mixed simplex, with all controls passing:
  `BARYCENTRIC_PRODUCT_NOT_A_STAIRCASE_COMMON_REFINEMENT`;
- any source, incidence, count, rational-containment or symmetry control
  failure:
  `BARYCENTRIC_PRODUCT_CONTROL_FAILED`.

The first outcome is a **DERIVED COMBINATORIAL** permission to attempt a
common-refinement action.  The second is a **DERIVED NEGATIVE** only for this
direct common-refinement claim; the carrier remains a canonical triangulation
of the cylinder.

## 7. Scope

No outcome constructs an improved/perfect action.  Lorentzian edge lengths,
nondegenerate simplex signatures, dust transport, the collective constraint,
coarse-to-fine maps, integration measure and canonical evolution are all
**OPEN**.  Only the targeted combinatorial verifier will run; the full suite
will not run.
