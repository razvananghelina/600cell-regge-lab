# Preregistered correction: constrained-response synthetic block

Date: 2026-08-21

First failure preserved in commit `cddb3ca`.

Replace only the call that sends the frozen `2+2` synthetic matrix through the
production `12+10` block splitter.  Compute the same frozen control directly:

```text
C_q=q^T C q,
y=-(C_q)^(-1) q^T B^T p,
K=p^T(A p+B q y).
```

Continue to require exactly

```text
K=18,
n^T B^T(1,0)=3.
```

Use the existing linear-solve helper, not an explicit inverse.  Do not change
the production reduction, any input, pivot, threshold, class comparator,
outcome or scope statement.  Run only the corrected targeted verifier.

