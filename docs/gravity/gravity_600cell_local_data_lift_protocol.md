# Protocol: exact universal local canonical-data lift

Date: 2026-08-19

Freeze this target-disclosed protocol before constructing the first 48-column
affine constraint or comparing any support set with a vertex star.

## Frozen inputs

- local-lift prior-art gate SHA-256
  `91a98568926afa6c556b143050a081354bcce28bc02da673d8863aa1aadc8aa7`;
- exact global rational-lift source SHA-256
  `65a097cd11dea830fd16bad988cd6d1b88ce4b84e0700b7dfaf7477a2c198ecb`;
- its first artifact SHA-256
  `1b6ac46a0ea4889f476cc71d51ac464c27caa6d4b6a9b2f2d74ff93da77b123f`;
- exact global result note SHA-256
  `c69e367fed93498705a30058134000bd77be4845589e7942a90b059e53aa3ecc`;
- complete admissibility source SHA-256
  `4d3595fbf418fc0876dba5a1129bdbcbd49d43a68ef9e6fd5fba2f0cb6e6873e`;
- refuted old local-lift source and artifact SHA-256
  `3adb80448e19fd99f0b8ec205497f11325b0a6f7a72c9a2785f9b65778707750`
  and
  `4065950aaac4180ec1cdd0b82f7a8bc403b2969c50d26cf14cc28592085cb2c5`.

## Complete hypotheses and disclosed target

For a tetrahedral cell with its deterministic local vertex order
`(v0,v1,v2,v3)`, posit one 6-by-8 rational block `X`:

```text
cell flex = X (sigma_v0,...,sigma_v3, strut_v0,...,strut_v3)^T.
```

Within each legitimate construction, the same `X` must be used for all 600
cells.  No cell-dependent coefficient is allowed.  The target is:

1. the affine global face constraints determine all 48 entries uniquely;
2. all 48 entries are nonzero;
3. substitution of the induced 3600-by-240 local-star lift into every one of
   the 6000 original face rows gives exact zero;
4. hence every cell-flex row has exactly the eight data at its cell's four
   vertices and every datum has exactly its vertex's 20-cell star as support.

Test both rational representatives `(2,5),(3,11)`, both exact local
right-inverse graphs, reversed face orientation, odd canonical relabelling,
and reversed metric sign.

## Mechanically independent construction

Do not call or copy the 3600-variable global rational elimination.

For every original face row and every global data coordinate, substitute the
universal local block symbolically.  This produces exact affine constraints
on only 48 unknowns.  Perform a new deterministic sparse rational elimination
on those 48 unknowns and the constant column.

Require rank 48 and no inconsistent zero-unknown row.  Reverse-substitute the
unique block, build the entire star-supported lift from that block, and check
all original rows directly.  Canonically serialize and print all 48 rational
entries in the artifact.

The existing exact 3600-variable result may be used only as frozen provenance
and as a uniqueness implication after this independent residual succeeds.

## Support and incidence controls

For the constructed lift, require exactly, not just by count:

```text
support(cell flex row) = {sigma_v, strut_v : v in that cell};
support(data column v) = {cells containing v}.
```

Check that each cell set has four vertices, each vertex star has 20 cells,
every local row has eight nonzero coefficients, and every data column has 120
nonzero flex coefficients across those 20 cells.

## Corrupted-image attack

Use the same frozen one-row corruption of the unsigned incidence map: on the
lexicographically first edge `{u0,v0}`, replace `(u0:1,v0:1)` by `(u0:1)`.
Apply the same 48-unknown universal-block ansatz.  It must leave an exact
affine inconsistency for every legitimate construction.  Record the first
obstruction.  A corrupted positive is a control failure.

## Reconciliation with the old local block

Independently reconstruct the old radial-scale plus normal-strut physical
displacement block.  For each construction:

1. verify that it induces the same ten local squared-length data;
2. solve its six old flex coefficients exactly from the local kernel;
3. require the new and old 6-by-8 blocks to differ;
4. verify that their physical displacement difference lies exactly in the
   six-dimensional local Poincare kernel;
5. for each representative, require baseline and alternate-right-inverse
   versions of the *new physical displacement response* to agree exactly,
   even though their flex-coordinate blocks differ.

This identifies the repair as a local Poincare gauge correction selected by
global face gluing.  It does not undo the frozen negative residual of the old
block.

## Outcome hierarchy

- `LOCAL_DATA_LIFT_CONTROL_FAILED`: provenance, geometry, affine rank,
  residual, exact support, corrupted-image, or old-block reconciliation fails.
- `LOCAL_DATA_LIFT_DISAGREEMENT_OPEN`: legitimate constructions disagree on
  existence of a unique universal block.
- `UNIVERSAL_LOCAL_DATA_LIFT_REFUTED`: all controlled constructions agree
  that no universal 6-by-8 block satisfies the complete equations.
- `LOCAL_STAR_DATA_LIFT_DERIVED`: every construction has a unique exact
  all-nonzero 6-by-8 block, exact star support, zero global residual, rejected
  corruption, and the stated old-block reconciliation.

The positive branch remains a first-order kinematic result.  The action
Hessian, symplectic form, constraint/gauge split, tick, tensor propagation,
`c`, `G`, and Planck scales remain **OPEN**.

## Reproducibility discipline

Register and commit the verifier before its first execution.  Run only this
targeted verifier and freeze its first artifact before beginning the boundary
action/Hessian calculation.

