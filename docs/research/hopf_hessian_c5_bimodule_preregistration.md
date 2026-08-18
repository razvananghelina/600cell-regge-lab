# Preregistration: arbitrary-multiplicity `C^5` bimodule boundary

Date: 2026-08-11

Parent ten-state result: `09142e8`.

## Framing correction

The next step proposed in the ten-state verdict was to search for a larger
`A5`-equivariant `C^5` bimodule whose first-order support could contain the
Hessian family.  That search is premature.  Metric-dimension-zero
orientability may forbid the entire commutative arena independently of `D`,
exactly as it did for the six-fibration algebra `C^6`.

This audit therefore tests the cheaper universal obstruction before any
larger carrier or Hessian embedding is constructed.

## Complete hypotheses

Let `A=C^5` and let `pi` be a unital faithful representation on a nonzero
finite-dimensional complex Hilbert space `H`.  Assume:

1. order zero, so `H` is an `A-A^op` bimodule;
2. a grading `gamma` commuting with `pi(A)`;
3. KO6 reality `J gamma=-gamma J`;
4. metric-dimension-zero orientability,
   `gamma=sum_k pi(a_k) J pi(b_k) J^-1`;
5. the exact derived `A5` action normalizes `A`, preserves `gamma`, and is
   compatible with `J`, hence acts on the left/right character labels of the
   bimodule.

No hypothesis about a Dirac operator, first order, connectedness or the
Hessian selector is needed for the proposed obstruction.  If it holds, it
applies a fortiori when those extra gates are imposed.

## Blind exact tests

Before inspecting any larger Hessian realization:

1. reconstruct the exact 60-element group from the derived action;
2. construct the five-point action on the conjugates of an index-five `A4`;
3. enumerate all orbits on ordered pairs of the five points;
4. record whether every orbit is invariant under pair reversal;
5. solve the orbit-sign constraints

```text
epsilon_(g i,g j)=epsilon_(i,j),
epsilon_(j,i)=-epsilon_(i,j).
```

The full orbit and sign ledger is to be written before any use of
`Hhat_X`.

## Algebraic reason the sign test is decisive

Every finite `C^5` bimodule decomposes as

```text
H=direct_sum_(i,j) H_(i,j).
```

Every Hochschild zero-cycle acts on `H_(i,j)` as a scalar.  Therefore an
orienting `gamma` must equal `epsilon_(i,j) I` on the entire multiplicity
space, with `epsilon_(i,j)` equal to `+1` or `-1`.  Extra multiplicity cannot
alter this conclusion.

KO6 reality sends a block to its reversed block and requires opposite signs.
Equivariance makes the sign constant on each `A5` orbit.  An orbit containing
both `(i,j)` and `(j,i)` is therefore forbidden; a diagonal block is forbidden
by the same equation with `i=j`.

## Decision boundary

- If a nonempty orbit admits the sign constraints, enumerate the smallest
  corresponding bimodule carriers blindly and only later compare their
  first-order map space with the Hessian family.
- If no nonempty orbit admits them, record a **DERIVED FULL-ARENA NO-GO** for
  `A=C^5` with arbitrary multiplicities.  Do not perform a pointless larger
  carrier search.

## Scope limits

Even a full-arena negative would not cover:

- a noncommutative algebra;
- broken `A5` covariance;
- a real structure not compatible with the `A5` symmetry;
- `J gamma=+gamma J`;
- positive metric-dimension Hochschild orientability.

Those are changed hypotheses, not survivors of this census.
