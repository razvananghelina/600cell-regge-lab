# Preregistration: spectrally ordered four-node KO6 carrier

Date: 2026-08-11

## Disclosed candidate

This protocol follows, and is not blind to, the exact node separator committed
in `b74b491`.  That result gives the four real Wedderburn nodes

```text
sizes:               6, 6, 12, 12
joint (u_edge,v_ref): (2,5), (2,-5), (phi-1,0), (-phi,0).
```

The candidate is the off-diagonal part of the full enveloping bimodule:

```text
H_off = direct_sum_(i != j) C^(n_i) tensor C^(n_j)*.
```

The algebra acts on the left tensor factor, the opposite algebra on the
right, and `J` exchanges `(i,j)` with `(j,i)` followed by coefficient
conjugation.  A spectral ordering of the joint labels orients every unordered
node pair; `gamma=+1` on the chosen direction and `-1` on its transpose.

This differs from the already failed standard sheet double of the enveloping
module.  There `J` exchanged two copies with the same central cell profile,
which made a metric-dimension-zero orientation cycle impossible.  Here `J`
transposes the central cells themselves.

## Canonicity attack frozen before calculation

The pair `(u_edge,v_ref)` does not by itself specify whether `u_edge` or
`v_ref` has lexicographic priority.  Nor does it privilege ascending over
descending order.  The verifier must therefore enumerate, rather than hide,
all eight signed lexicographic readings:

```text
priority in {u then v, v then u},
direction of u in {ascending, descending},
direction of v in {ascending, descending}.
```

Global reversal of both directions is allowed to give the grading-reversed
carrier, but it will not be silently quotiented out in the raw count.

No individual reading may be called geometrically selected.  A conclusion is
robust only if it holds for all eight.

## Frozen exact gates

For every one of the eight readings, compute exactly:

1. the six oriented off-diagonal central cells;
2. the Hilbert-space dimension
   `2*sum_(i<j)n_i*n_j`;
3. faithfulness of the left and opposite representations;
4. `J^2=+1` and `J gamma=-gamma J`;
5. order zero for the complete matrix algebra, proved at cell level and
   checked on the matrix-unit action rule;
6. the explicit metric-dimension-zero Hochschild cycle made from central
   block units and its equality to `gamma`;
7. the minimal-projector intersection matrix `Q`, its Pfaffian and
   determinant;
8. the complete set of first-order-compatible odd cell-to-cell block
   positions and whether its induced central-link graph is connected.

The intersection convention is the one already used by the repository's
blind central Krajewski census: one copy of `(i,j)` contributes `+1` to
`Q_ij`, its `J` partner contributes the antisymmetric transpose, and
`Q=mu-mu^T`.

## Acceptance and kill boundaries

- If all eight readings give a faithful, orientable order-zero KO6 carrier
  with nondegenerate `Q`, record **STRUCTURAL ROBUST CARRIER EXISTENCE**.
  This means the node-order ambiguity does not obstruct Poincare duality.
- If any gate depends on choosing a favorable ordering, report the exact hit
  fraction and label it **PATTERN**, not selection.
- If every reading is degenerate or non-orientable, close this continuation.

Even full acceptance does **not** produce a finite spectral triple.  No
Dirac matrix is selected here.  Nonempty legal first-order blocks and a
connected possible-link graph are necessary possibilities, not evidence that
a canonical nonzero, connected `D` exists.

No Hessian, particle module, mass, coupling or Standard-Model target will be
used.  Only the new targeted verifier will be run; the full suite remains
excluded by user instruction.
