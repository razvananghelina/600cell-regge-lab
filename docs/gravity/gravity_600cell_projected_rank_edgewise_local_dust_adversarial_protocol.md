# Adversarial protocol: consistent-matrix reconstruction of local P1 dust

Date: 2026-08-19

Primary artifact commit: `de0c220`  
Primary artifact SHA-256:
`53463e5271301ae41eb26564875d26991ddea8024a9e09ae3c302d428ad39779`

## Independent construction

The primary verifier computed each tetrahedral volume with a Gram determinant
and directly added `V_t/4` to its four vertices.

The audit will instead:

1. reconstruct the 600 source tetrahedra with NetworkX maximal cliques;
2. rebuild the projected carriers with the direct rank split;
3. compute every tetrahedral volume from the Cayley--Menger determinant;
4. assemble the standard consistent local `P1` mass matrix

   ```text
   M_t = V_t/20 * [[2,1,1,1],
                   [1,2,1,1],
                   [1,1,2,1],
                   [1,1,1,2]];
   ```

5. derive the lumped nodal weights only as global row sums of the assembled
   sparse matrix.

The row sums should recover `V_t/4` algebraically, but the carrier, volume
formula and assembly route differ from the primary calculation.

## Frozen gates

- A symbolic local row sum must equal `V/4` in every row.
- Cayley--Menger and Gram volumes must agree on a regular tetrahedron to
  relative error `<2e-14`.
- The reconstructed primary Gram-weight byte digest must equal the frozen
  digest stored in the primary artifact at both levels.
- Cayley--Menger consistent-matrix row sums and Gram direct weights must agree
  pointwise to relative error `<2e-10`.
- Their sums, minima and maxima must reproduce the frozen primary scalars to
  relative error `<2e-10`.
- All assembled consistent matrices must be symmetric, all local eigenvalues
  positive and every row sum positive.

## Positive and negative controls

Positive control: on the regular 600-cell, the row-sum weights must be uniform
under the same numerical tolerance.

Negative control: replace the local weights on each irregular projected
carrier by the globally uniform rule `V_total/N_vertices`.  This rule is
positive, conserves total mass and is symmetry invariant, but it is not exact
on the global `P1` hat function at a vertex whose true assembled weight differs
maximally from the uniform value.  Require the maximum relative discrepancy
to exceed `0.1` on both irregular carriers.

This negative control demonstrates why positivity, total conservation and
symmetry alone do not select the P1 rule.

## Outcome

Write

```text
reproducible/gravity_600cell_projected_rank_edgewise_local_dust_adversarial.json
```

and assign exactly one outcome:

1. `ADVERSARIAL_P1_LOCAL_DUST_DISAGREEMENT` if any gate fails;
2. `ADVERSARIAL_P1_LOCAL_DUST_CORROBORATED` otherwise.

Corroboration remains conditional on the **STRUCTURAL** continuous nodal P1
dust ansatz.  It still does not construct local lapse dynamics.
