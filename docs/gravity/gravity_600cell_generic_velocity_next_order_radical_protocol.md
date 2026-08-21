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

## Completed inventory before normalization

The coefficient-only inventory contains exactly eight square-root bases.  Set
`u=v^2>=0`.  Their exact factorizations and positive-branch replacements are

```text
u+4
  -> sqrt(u+4),

3u+8
  -> sqrt(3u+8),

3u^2+20u+32=(u+4)(3u+8)
  -> sqrt(u+4)*sqrt(3u+8),

9u^3+84u^2+256u+256=(u+4)(3u+8)^2
  -> (3u+8)*sqrt(u+4),

3u^3+32u^2+112u+128=(u+4)^2(3u+8)
  -> (u+4)*sqrt(3u+8),

9u^4+120u^3+592u^2+1280u+1024
  =(u+4)^2(3u+8)^2
  -> (u+4)(3u+8),

27u^5+432u^4+2736u^3+8576u^2+13312u+8192
  =(u+4)^2(3u+8)^3
  -> (u+4)(3u+8)*sqrt(3u+8),

27u^6+540u^5+4464u^4+19520u^3+47616u^2+61440u+32768
  =(u+4)^3(3u+8)^3
  -> (u+4)(3u+8)*sqrt(u+4)*sqrt(3u+8).
```

Every replacement is branch-safe for every real `v` because `u+4>0` and
`3u+8>0`.  The verifier must prove all eight squared identities before using
them and must confirm that the post-normalization inventory contains only the
two primitive positive radicals `sqrt(v^2+4)` and `sqrt(3v^2+8)`.
