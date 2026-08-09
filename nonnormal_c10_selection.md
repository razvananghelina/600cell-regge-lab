# Non-normal phase contexts select `2I`

## Conditional theorem

Let `G` be a finite subgroup of `SU(2)`.  If `G` contains a non-normal cyclic
subgroup `C10`, then `G` is the binary icosahedral group `2I`.

### Proof by the finite-`SU(2)` classification

Finite subgroups of `SU(2)` are cyclic, binary dihedral, or the exceptional
groups `2T,2O,2I`.

- In a cyclic group every subgroup is normal.
- In `Dic_n`, every element outside the normal rotation subgroup `C_(2n)`
  has order four.  Hence every `C10` lies inside `C_(2n)`.  A cyclic group has
  a unique subgroup of each order, so that `C10` is characteristic in the
  rotations and normal in `Dic_n`.
- `|2T|=24` and `|2O|=48`, so neither contains an element of order 10.
- Exact enumeration in `2I` gives six `C10` subgroups in one conjugacy orbit;
  therefore each is non-normal.

This proves the stated classification theorem.  **DERIVED conditional on the
standard finite-`SU(2)` classification.**

## Relation to the non-binary principle

The Fibonacci phase lift produces `C10`, but `Dic_5/D7` proves that `C10`
alone does not select `2I`.  The additional condition can be formulated as:

> **Non-normal context principle.**  The elementary phase context is not a
> globally invariant subsystem; all of its conjugate contexts belong to the
> same reality.

Under this principle the single normal `C10` of `Dic_5` is excluded, while
the six conjugate `C10` contexts of `2I` survive.  A chosen context has index
12, producing the homogeneous space `2I/C10` underlying the twelve
icosahedral vertex directions.

The theorem is mathematical.  The non-normal context principle is a new
**STRUCTURAL AXIOM**, motivated by “reality is non-binary”; it is not derived
from S01.  It should not be advertised as a theorem of physics.

## Strong relational closure

The six contexts satisfy stronger exact properties:

- their total intersection is precisely the binary center `{+1,-1}`;
- every pair of distinct contexts also intersects precisely in that center;
- any two distinct contexts generate the full group `2I`;
- the union of all six contexts contains only 50 elements, so the remaining
  70 arise from mixed relational products rather than mere collection.

Thus one context is incomplete, while the relation between any two distinct
contexts closes to the whole 120-element reality.  **DERIVED.**  Reading this
as a formal realization of “reality is relational rather than binary” is
**STRUCTURAL interpretation**.

## The `Dic_5` counterexample becomes the local normalizer

For a chosen `C10=H`, its normalizer inside `2I` has order 20 and element-order
histogram

`{1:1, 2:1, 4:10, 5:4, 10:4}`,

the exact signature of `Dic_5`.  Hence

`N_2I(H) isomorphic to Dic_5`.

The earlier global counterexample is therefore precisely the symmetry that
preserves one local phase axis.  Passing from the local normalizer to all
conjugate contexts gives

`|2I/N_2I(H)|=6` unoriented axes,

while retaining phase orientation gives

`|2I/H|=12` oriented directions.

This turns the negative control into a local-to-global theorem.  **DERIVED.**

## Context graph and bootstrap closure

Join two contexts when together they generate `2I`.  Every pair does, so the
context graph is the complete graph `K6`.  Its adjacency spectrum is

`{5,-1,-1,-1,-1,-1}`,

and its Laplacian spectrum is

`{0,6,6,6,6,6}`.

Consequently its degree and spectral gap are exactly `(5,6)=(a1,b1)`.  This
is an exact **DERIVED internal closure** after `C10` and `2I` are already in
place.  It is not a non-circular derivation of the original `a1=5`, because
`C10` was reached using the fifth-root phase.  Calling it a foundational
derivation would be false; calling it a bootstrap consistency loop is fair.

## Status ledger

- **DERIVED:** exact `2I` group, 24 order-ten elements, six conjugate `C10`
  subgroups, non-normality and index 12.
- **DERIVED:** explicit `Dic_5` countercontrol with a unique normal `C10`.
- **DERIVED:** `N_2I(C10)=Dic_5`, central pairwise intersections, and
  two-context generation of all `2I`.
- **DERIVED:** context graph `K6` and its exact `(degree,gap)=(5,6)`.
- **DERIVED CONDITIONAL:** uniqueness of `2I` among finite `SU(2)` subgroups
  once a non-normal `C10` is required.
- **STRUCTURAL:** elevating non-normality to the foundational relational
  principle.
- **OPEN:** deriving that principle from something weaker and connecting its
  twelve cosets to a physical observable rather than only icosahedral
  geometry.
- **NOT CLAIMED:** Standard-Model selection, a matter spectral triple or a
  dynamical gauge field.

Exact verifier:
`reproducible/verify_nonnormal_c10_selection.py`.
