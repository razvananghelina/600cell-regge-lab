# Protocol: ordered balanced slab on the projected rank-edgewise carrier

Date: 2026-08-19

Prior-art commit: `dabc098`.

## Scope

This protocol freezes a purely combinatorial selection test before any local
Regge or dust equation is evaluated.  It uses only:

```text
K0 = P(sd K_600),
K1 = P(Esd_2(sd K_600)),
r(v)=dim(face represented by a K0 vertex),
c(v,w)=r(v)+r(w) mod 4,
past < future.
```

The radial coordinates are irrelevant to the product topology.  No continuum
target, action value, lapse, Hessian or spectrum may be imported.

## Independent reconstruction

The verifier will rebuild the 600-cell from `commons/cell600.py`, find its
600 tetrahedra from adjacency, construct the barycentric chamber complex and
then form the eight-child direct rank split

```text
v0 m01 m02 m03        v1 m01 m12 m13
v2 m02 m12 m23        v3 m03 m13 m23
m01 m02 m03 m13       m01 m02 m12 m13
m02 m03 m13 m23       m02 m12 m13 m23.
```

This direct split is intentionally different from the colour-scheme
enumeration used by the original carrier verifier.

## Frozen tests

### A. Spatial controls

For both `K0` and `K1`:

1. reproduce the certified f-vector;
2. certify that the tetrahedron dual graph is connected;
3. certify that every triangle has incidence two.

### B. Colouring

1. On `K0`, use the literal face rank.
2. On `K1`, use the endpoint-rank residue `c(v,w)` above.
3. Verify that every spatial edge has differently coloured endpoints and
   every tetrahedron contains `{0,1,2,3}`.
4. Starting from every one of the `4!` assignments on one seed tetrahedron,
   propagate across shared triangles.  Record the exact number of complete
   labelled four-colourings.  Any branch requiring a choice after the seed is
   rejected.
5. Reconstruct all 120 left, 120 right and conjugation actions and verify that
   the declared rank-residue colouring is preserved.

### C. Product topology

For each of the 24 linear orders of the four colour classes, split every
tetrahedral prism into four pentachora by the standard staircase rule.

For the rank order `0<1<2<3` on each carrier:

1. top pentachora are distinct;
2. every codimension-one face has incidence one or two;
3. the incidence-one faces are exactly the two copies of the input spatial
   tetrahedra;
4. the bottom and top boundary counts agree and there are no side-boundary
   faces because the spatial carrier is closed;
5. the expected top count is `4*f3(Ki)`.

For the complete set of 24 colour orders, record:

1. the number of distinct labelled slab complexes;
2. the number preserved by the spatial `H4` action;
3. the permutation induced by time reversal;
4. fixed points and orbit sizes under time reversal.

The look-elsewhere number `N_order` is part of the result and may not be
discarded after execution.

## Selection boundary

There are two different standards, frozen in advance:

- **existence:** a proper rank-derived colouring and a conforming staircase
  slab exist;
- **canonical selection:** spatial rank geometry plus declared time
  orientation leave exactly one legitimate ordered slab.

Passing existence but finding `N_order>1` means the temporal carrier is only
**STRUCTURAL**.  In that case no local lapse or Hessian will be run on one
chosen schedule in this mission, because doing so would silently fit a point
inside an unresolved finite family.

A reverse-time partner is not identified with the original after time
orientation has been declared.  It is reported as the oppositely oriented
carrier.  Spatial `H4` equivalence is allowed and will be tested exactly.

## Artifact and reproducibility

The registered verifier will be

```text
reproducible/verify_gravity_600cell_projected_rank_edgewise_balanced_slab.py
```

and will write

```text
reproducible/gravity_600cell_projected_rank_edgewise_balanced_slab.json.
```

Only this verifier and the static registry/documentation guards are run.
The full suite is explicitly excluded by the user's instruction.

