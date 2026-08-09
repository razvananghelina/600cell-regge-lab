# Corrected target-independent `(2,3,5)` incidence census

Date: 2026-08-09

## Status and provenance

This file corrects the enumeration in commit `36bd6825b6b64edc6c95c1ac67a21a80693318fc`.
That commit is preserved as audit history.  Its reported `N=3` is false.
Commit `d2b02eb5fc1a42428bc06ff67b0cef4c6d2da295` performed a comparison against
that incomplete list and is therefore also superseded.

The corrected enumeration code contains no external comparison module.  Its
full output is frozen in
`reproducible/orbifold_incidence_preregistered.json` and regenerates exactly.

This correction cannot restore genuine blinding: the desired comparison was
already known before it was committed, and local Git history is mutable.  It
does, however, freeze the complete corrected list before the corrected
comparison verifier and conclusion are added.

## Why the old twisted sum vanished

For a representative `d` of a double coset `H d K`, the old code summed over
all `(h,k)` in `H x K` while testing whether `h d k` remained in `H d K`.
That membership test is tautologically true.  The coefficient therefore
factorized as two full character sums and vanished whenever either character
was nontrivial.

That is not a twisted Mackey kernel.  The orbit is `(H x K)/L_d`, where
`L_d=H intersect d K d^-1`.  A covariant line-valued kernel exists when the
two characters agree on `L_d`; it is then one-dimensional and nonzero.  In
all three cross-cell cases every `L_d` is the binary center `C2`, so existence
is exactly equality of central parity.

**DERIVED NEGATIVE:** the old claim that every nontrivial twisted incidence
operator vanishes.

## Exact geometric certificate

The independent construction uses the exact quaternion action over
`Q(phi)`.  The twelve vertex rays have normalized mutual dot products

`-1, -sqrt(5)/5, +sqrt(5)/5`.

Edges are the thirty maximal-dot-product pairs; faces are the twenty graph
triangles.  No rotation eigensolver or geometric tolerance selects a cell.
The resulting double-coset data are:

| pair | double cosets | group-set sizes | incidence counts by orbit |
|---|---:|---|---|
| vertex--edge | 6 | `20` each | one `60`, five `0` |
| vertex--face | 4 | `30` each | one `60`, three `0` |
| edge--face | 10 | `12` each | one `60`, nine `0` |

The integer index assigned to a double coset is not invariant.  With exact
oriented quaternion axes and the old right-coset ordering, the incidence
labels are `(0,1,8)`.  Replacing only the `C6` axis by its antipode reproduces
the eigenvector convention's `(0,3,6)`.  In both cases the geometric statement
is identical: exactly one orbit contains all sixty incidences and no orbit is
mixed.  **DERIVED.**

## Hom-space and canonicity audit

For the twenty induced line modules, the exact Gram matrix has diagonal Hom
dimensions

`3 (x8), 4 (x2), 7 (x4), 8 (x2), 15 (x2), 16 (x2)`

and ordered off-diagonal histogram

`0:200, 2:32, 3:8, 4:60, 6:48, 7:4, 10:24, 14:2, 15:2`.

There are no one-dimensional full Hom spaces.  Nevertheless, support on the
single actual incidence orbit selects one Mackey summand, hence one line in
the full Hom space.  Changing a base-fiber trivialization only rescales the
whole map and does not change its kernel or cokernel.  **DERIVED canonicity,
conditional on the stated support definition.**

A map using a linear combination of incidence and non-incidence orbits, or
freely chosen Schur-channel coefficients, is excluded as **STRUCTURAL**, not
canonical.

## Corrected complete census

Matching central parity leaves exactly

- `C4 -> C10`: `4*10/2 = 20` map/adjoint pairs;
- `C6 -> C10`: `6*10/2 = 30` map/adjoint pairs;
- `C6 -> C4`: `6*4/2 = 12` map/adjoint pairs.

Thus there are **62** primitive incidence-map adjoint pairs.  All exact
kernel and cokernel characters are recorded in the JSON file.  Their channel
ranks attain the representation-theoretic upper bound modulo both `601` and
`1801`; this sandwiches the characteristic-zero ranks exactly.

Of the sixty parity-compatible normalized pairs `F -> E -> V`, exactly one
has zero composition in `Z[z]/Phi_60(z)`:

`F0 -> E2 -> V0`.

It is the oriented cellular complex, with

`H2=rho1, H1=0, H0=rho1`.

Counting one map with its adjoint as one object and this complex with its
adjoint as one object gives `N=63`.  Their full virtual-index multiset has 28
distinct characters and is frozen in the JSON file.  Virtual indices are
recorded only as a required census: they are endpoint differences and contain
no operator-dependent information.

## Target-independent ledger

- **DERIVED:** exact f-vector, incidence-orbit purity, and double-coset counts.
- **DERIVED:** Gram/Hom census; relation-lattice rank `20-9=11`.
- **DERIVED NEGATIVE:** `N=3` and “all nontrivial twists vanish”.
- **DERIVED:** 62 nonzero primitive incidence-map adjoint pairs.
- **DERIVED:** one short complex, `F0 -> E2 -> V0`, with zero middle cohomology.
- **STRUCTURAL:** incidence-coset support as the definition of canonicity.
- **OPEN IN THIS FILE:** every external module comparison.

Verifier: `reproducible/verify_incidence_operator_enumeration.py`.
