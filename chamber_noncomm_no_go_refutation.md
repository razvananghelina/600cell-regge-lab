# Refutation of the proposed chamber commutativity theorem

Date: 2026-08-09

## Proposed statement

Let the Hilbert space be the 120 oriented icosahedral chambers, with the
derived three-regular chamber adjacency `D`, orientation grading `gamma`, and
geometric reflection `J`.  The proposed theorem asserted that every
faithfully and unitally represented finite star algebra satisfying order
zero, first order, `[gamma,A]=0`, and nondegenerate intersection form must be
commutative.

## Verdict

The statement is **REFUTED**.

There is an exact counterexample with

`A = M2(C) direct-sum C direct-sum C direct-sum C`.

It uses the fixed `H,D,gamma,J`, satisfies every listed premise, has nonzero
one-forms and has strict integral Poincare duality.  It fails orientability
and connectedness, neither of which was included in the proposed statement.

The counterexample is a nonlocal **STRUCTURAL** algebra on the chamber
carrier.  It is not selected by chamber locality or by the founding group.
That does not weaken it as a counterexample to a theorem quantifying over
all faithful representations.

## The error in the proposed reduction

From `[gamma,A]=0` one obtains two representations

`rho_+: A -> M60(C)` and `rho_-: A -> M60(C)`.

It does not follow that either restriction is faithful, nor that they are
equivalent copies of the total faithful representation.  Only

`ker(rho_+) intersect ker(rho_-)=0`

is forced by total faithfulness.

In the counterexample, the `M2` summand is absent from one sheet and one
scalar summand is absent from the other.  Their kernels are complementary.
Consequently the suggested reduction to two mutually commuting faithful
copies of the same `M_k` does not apply.  The bound `k^2<=60` is therefore not
a proof strategy for the stated general theorem.  Even in the simple-factor
case, representation multiplicities and divisibility, not only this
inequality, must be tracked.

The existing contraction search also cannot support the unrestricted
theorem: its variables encode set partitions, so every algebra in that search
is commutative by construction.  Finding only commutative algebras in a
commutative search space is not evidence against arbitrary star algebras.

## A valid local no-go for the exact C5 witness

There is a narrower theorem which the data do prove.  For the 30-block exact
`C5` witness, the common-refinement dimensions

`c_ij=dim(P_i J(P_j)H_+)`

have census

`70 x 1, 10 x 2, 10 x 3`

over the 90 nonzero cells.  The graph joining `i` and `j` whenever `c_ij=1`
is connected on all 30 nodes.

If the scalar block on node `i` were amplified to `M_(n_i)`, order zero would
require every nonzero cell dimension to be a multiple of `n_i*n_j`.  Along a
unit cell this forces `n_i=n_j=1`; connectivity forces all 30 factors to be
scalar.

Thus this particular partition witness admits no noncommutative matrix
amplification which retains the same central supports.  **DERIVED scoped
no-go.**  It says nothing about nonlocal central projectors or different
Krajewski multiplicities, which is exactly where the counterexample lives.

## Exact multiplicity certificate

Index the four simple summands by dimensions

`n=(2,1,1,1)`.

On the positive sheet take the nonzero bimodule multiplicities

`mu_01=1, mu_12=3, mu_23=1, mu_11=54`.

The weighted sheet dimension is

`1*(2*1) + 3*(1*1) + 1*(1*1) + 54*(1*1) = 60`.

The negative sheet contains the transposed cells, as required by `J`.  For
minimal `K0` projections the intersection form is

```text
Cap = mu-mu^T
    = [[ 0,  1,  0,  0],
       [-1,  0,  3,  0],
       [ 0, -3,  0,  1],
       [ 0,  0, -1,  0]].
```

It is antisymmetric, has rank four, Pfaffian one and determinant one.
**DERIVED exact certificate.**

## First-order Dirac block

Order the first six positive-sheet dimensions as

`(01 of dimension 2) + (23 of dimension 1) + (12 of dimension 3)`.

On them use

```text
T_path = [[0,I3],
          [I3,0]].
```

Each nonzero block changes only one Krajewski index, so first order holds.
The remaining 54 dimensions are the scalar-scalar `(1,1)` cell.  Any
symmetric `54 x 54` operator is first-order legal there.

The verifier checks order zero and first order on a full seven-element
linear basis of `M2(C)+C^3`, not only on central projectors.

## Why this is the exact fixed chamber `D`

Identify the negative chamber sheet with the positive one using `J` and set

`S=(D J)|H_+`.

For the exact integer chamber matrices:

- `S` is symmetric and invertible;
- `nullity(S+I)=6` exactly over `Q`;
- `nullity(S-I)=0`.

Thus `S` splits orthogonally as

`S=(-I6) direct-sum R54`.

The displayed `T_path` is explicitly unitarily congruent to `-I6`.  The
`R54` block is placed on the scalar-scalar cell, where it is unrestricted.
Therefore the abstract counterexample is transported to the exact chamber
`D` by a unitary of the form

`diag(U,conjugate(U))`.

Such a unitary preserves the standard sheet-swap antiunitary `J` and the
grading.  Hence this is an exact existence proof on the fixed carrier, not a
change of spectrum or enlargement of `H`.  **DERIVED.**

## Failed additional gates

The 54-dimensional diagonal `(1,1)` cell occurs with both grading signs.
Every represented Hochschild zero-cycle has the same value on a paired
positive/negative coordinate in that cell, whereas `gamma` has opposite
values.  Metric-dimension-zero orientability therefore fails.

The represented element `(0,1,0,1)` is non-scalar and commutes with `D`, so
connectedness also fails.

These are **DERIVED NEGATIVES**.  They do not rescue the proposed theorem,
because neither was one of its premises.

## Correct boundary

- **REFUTED:** order zero + first order + grading compatibility + faithful
  total representation + nondegenerate or unimodular intersection form
  forces commutativity.
- **OPEN:** the same statement after adding metric-dimension-zero
  orientability and connectedness.
- **OPEN:** whether sheetwise faithfulness or equivalence follows from some
  independently derived physical principle.
- **STRUCTURAL:** the counterexample's nonlocal block allocation.
- **NOT CLAIMED:** a manifold-like matter triple, a canonical algebra, or a
  Standard Model sector.

Exact verifier:
`reproducible/verify_chamber_noncomm_no_go_refutation.py`.
