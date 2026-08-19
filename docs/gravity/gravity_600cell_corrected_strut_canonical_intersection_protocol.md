# Protocol: corrected-strut/canonical-graph intersection census

Date: 2026-08-19

Status: no singular value, rank or nullity of an actual
`G_corrected-C_canonical` block has been evaluated.

## 1. Frozen provenance

Require exact SHA-256 values:

```text
prior-art gate
a56c5042be8596876ca1b2b45085069bcb1da7d1d859ac12d98f5e7e8420fc4f

target-blind corrected carrier
e8035fb9c35ad693d1dd2adbda79485b6dd8d42bdf40a95b70a92466e47027d7

primary non-equality artifact
5652b1371563ff11919be130af15f5b48850e2cc65a50ec35e5de85fdb587f90

adversarial non-equality artifact
3b0fd6da76195279f1beac540c326c61eff5e3172a63bb89baf69502254c5b1f

frozen response source
e461296a965c9b80fb89fae5660ce642858f3d3dfa0b24ccdecc2aced53c7047

frozen tangent archive
816c605da2a655442bbadce7a23965f0822f99e7bdc1d0a4a27af548de85446b
```

The non-equality artifacts disclose the target but do not contain or imply an
intersection nullity.  This is a target-disclosed completion, not a blind
prediction.

## 2. Exact graph reduction

For both parities, reconstruct the 65 orbit types in the response domain and
the five pole positions from edge geometry.  In every minimal `2T` sector of
dimension `d`, require both projected matrices

```text
G, C : C^(5d) -> C^(65d)
```

to have rank `5d`, and their pole row blocks to equal the same literal
identity below `1e-13`.  Only after this control may the graph intersection be
computed as

```text
ker(G-C).
```

If the pole identities disagree, stop with a control failure; do not report a
nullity from the difference.

## 3. Frozen numerical variants and zero calibration

Reconstruct all four frozen derivative variants.  For each sector form
`D_variant = G-C_variant`.  Compute all singular values with the divide-and-
conquer SVD and define a matrix uncertainty bound as the sum of:

1. the maximum Frobenius radius of the four Flint lift balls;
2. the maximum spectral norm `||D_variant-D_operational_primary||_2`;
3. the committed corrected-carrier binary discrepancy times `||G||_2`;
4. `50 eps_machine` times the largest matrix norm and condition number used.

Call this `epsilon_matrix`.  Classify each operational-primary singular value:

```text
ZERO       sigma <= 10 epsilon_matrix,
NONZERO    sigma > 100 epsilon_matrix,
OPEN       otherwise.
```

The sector nullity is resolved only if no singular value is `OPEN` and all
four variants give the same zero count under their common bound.  Record every
singular value, the bound, gap ratio `sigma_min/epsilon_matrix`, rank and a
basis of any resolved kernel through its projector (not arbitrary basis
coefficients).

No expected nullity is preregistered.

## 4. Fixed controls

1. **Positive intersection control:** replace `C` by `G`; require nullity
   exactly `5d`.
2. **Negative intersection control:** replace `C` by `G+E`, where `E` embeds
   the `5d x 5d` identity in the first `5d` non-pole rows; require nullity
   zero.
3. Apply the same frozen nonsingular upper-triangular coefficient change `T`
   to both graphs; `ker((G-C)T)` must have the same nullity.
4. Check the independent image formula
   `dim(im G intersect im C) = rank(G)+rank(C)-rank([G,C])` with one common
   calibrated threshold.
5. Reverse the source/target roles on the first committed staircase diagonal;
   require the difference matrix and at least one singular value to change
   above the carrier roundoff.
6. Complex conjugation and even/odd parities are recorded separately.  They
   may agree, but are not combined to choose a verdict.

## 5. Complete ledger and outcome

The ledger is exactly 14 actual sectors plus the fixed controls; there is no
search over carrier mixtures, graph conventions or retained small singular
values.

- `CORRECTED_STRUT_CANONICAL_INTERSECTION_RESOLVED` if all 14 nullities are
  calibrated and every control passes;
- `CORRECTED_STRUT_CANONICAL_INTERSECTION_NUMERICALLY_OPEN` if a singular
  value falls in the open band or the four variants disagree without a hard
  control failure;
- `CORRECTED_STRUT_CANONICAL_INTERSECTION_CONTROL_FAILED` otherwise.

The verifier passes for either a resolved or honest numerical-OPEN outcome.
Any material resolved claim requires a separate mechanically different
adversarial verifier before consolidation.

## 6. Claim boundary

A zero intersection would refute pure-strut canonical freedom on this fixed
slab.  A nonzero intersection would select only those common directions; it
would not identify them as gauge or physical.  Neither result decides the
full scale-plus-strut carrier, curvature propagation, a graviton dispersion,
an absolute tick, `c`, `G`, Planck units or particle masses.

Run only this targeted verifier and the static registry guard.

