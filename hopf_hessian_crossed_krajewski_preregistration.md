# Preregistration: blind Krajewski support census for `M5+M5+M5+M15`

Date: 2026-08-11

Parent crossed-product verdict: `d886acc`.

## Purpose and evidential limit

The five-point crossed product canonically fixes the algebra

```text
B=M5_0+M5_+ +M5_- +M15_3,
```

but not a physical bimodule.  This protocol asks only whether abstract
central Krajewski supports exist and how many there are.  Existence is not
selection.  A large survivor count is a negative evidential result even if
some supports pass all combinatorial gates.

No Hessian matrix, `Box` configuration or selector character may be used in
STEP 1.

## Complete central-level hypotheses

Let `mu_ij` be the positive-grading multiplicity of the irreducible
`B-B^op` cell `(i,j)`, with node sizes

```text
n=(5,5,5,15).
```

For the minimal-support census restrict to binary multiplicities.  This is
exhaustive for minimal Hilbert dimension because duplicating an occupied
cell changes neither central support, intersection rank nor first-order
legality and only increases dimension.

Impose:

1. KO6 reality supplies the negative sheet with multiplicity `mu^T`;
2. metric-dimension-zero orientability forbids diagonal cells and forbids
   simultaneous occupation of `(i,j)` and `(j,i)`;
3. the intersection form `Q=mu-mu^T` has rank four;
4. every node is represented faithfully;
5. first-order odd blocks may connect a positive cell `(i,j)` to the reverse
   of another positive cell `(k,l)` only when `i=l` or `j=k`;
6. the graph seen by commutators with the centre is connected: for every
   directed path `i->j->l`, a legal odd block may link central characters
   `i` and `l`, and these outer-node links must connect all four nodes.

Gate 6 is a necessary central connectedness condition, not a proof of full
matrix-algebra connectedness or existence of numerical Dirac coefficients.

## Exact finite census

There are six unordered node pairs.  Each is absent, oriented one way, or
oriented the other way:

```text
N_total=3^6=729.
```

For every support record:

- occupied directed edges;
- total KO6 Hilbert dimension
  `2 sum_(i,j) mu_ij n_i n_j`;
- `rank(Q)` and `det(Q)`;
- central first-order link graph and its number of components;
- the number of legal odd cell-pair blocks;
- invariance under the exact conjugation permutation swapping `M5_+` and
  `M5_-` while fixing `M5_0,M15_3`;
- whether that conjugation preserves or reverses the grading.

Write the complete survivor multiset and the exact number of distinct
supports.  Identify the minimum Hilbert dimension and count every minimizer.

## Decision boundary

- If there are no survivors, close the entire crossed-product Krajewski route
  at the necessary central level.
- If there are many survivors or no symmetry-selected unique orbit, label the
  result **STRUCTURAL EXISTENCE / SELECTION NEGATIVE**.  Do not choose one for
  its later Hessian behavior.
- Only if a small, independently symmetry-selected set survives may a
  separate committed STEP 2 compare its exact map space with the Hessian.

## Scope

Passing this census would not establish order zero beyond the standard
bimodule construction, a real antiunitary on matrix factors, full first
order, nonzero one-forms, full connectedness, orientability at the matrix
level, or Poincare duality over the relevant real form.  Those gates remain
separate.  Conversely, failing a necessary gate here is decisive for this
algebra type.
