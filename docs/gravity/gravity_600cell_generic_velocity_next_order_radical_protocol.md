# Radical-normalization protocol: next-order fixed-state census

Date: 2026-08-21

Scaled-jet implementation commit: `d8dde35`.

Preserved second-execution record:

```text
reproducible/gravity_600cell_generic_velocity_next_order_second_timeout.json
SHA-256 d17cd72341902f1a69dc111056179c964f3ac720b5fd6bf675d5e7772d8b6a6c
```

The second execution produced both exact next-order coefficients but left
positive polynomial square roots expanded.  The leading-zero structural
check did not reduce, and simplification of `a_C-a_P` was interrupted after
more than five minutes.  No common-root or composition verdict existed.

Before changing the verifier:

1. enumerate every square-root base in the two already-derived coefficients;
2. factor each base over `u=v^2`;
3. admit a replacement only when the factorization is an exact identity and
   every factor has a fixed positive sign for real `v`;
4. commit the complete list before applying it to either coefficient;
5. normalize the zeroth-order residuals and coefficients before constructing
   polynomials in `a`.

The seven factorizations already frozen in commit `f424f31` remain admissible
but are not assumed complete for the new derivative order.  No root,
velocity sample or outcome may be inspected while extending the list.
