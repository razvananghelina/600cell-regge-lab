# Preregistered precision correction: two-slab branch audit

Date: 2026-08-16

Original protocol commit: `29dcfa5`

First implementation/result commit: `1148dbf`

Status: **frozen after identifying an implementation limitation and before
re-evaluating the corrected branch certificates**.

## Limitation found by the post-result audit

The action values, central differences and momentum identities in
`1148dbf` were evaluated at 100 decimal digits.  However, the auxiliary Gram
signature audit converted each perturbation to binary64.  In particular,
`exp(1e-20)` rounds to exactly one in binary64.  Therefore the reported 902
branch-audit calls per parity were not 902 distinct represented geometries:
some of the smallest-step calls repeated the base geometry.

The passed binary64 `+/-1e-6` envelope, the minimum Gram-eigenvalue margins
and the `~1e-96` imaginary contamination make a branch change at the tiny
steps implausible, but they do not literally satisfy the preregistered
100-decimal contraction requirement.  Until this correction passes, the
outcome `TWO_SLAB_GLUING_CONTROL_PASSED` is **PROVISIONAL**.

## Frozen correction

Do not change any geometry, action, derivative step, tolerance, orbit map,
attempt count or momentum gate.  During every arbitrary-precision one-slab
and direct two-slab action evaluation already required by the protocol:

1. form each real symmetric `4 x 4` simplex Gram matrix at 100 decimal
   digits;
2. compute the leading principal minors `Delta_1,...,Delta_4` at the same
   precision;
3. require every minor to be nonzero and count the sign changes in
   `(1, Delta_1,...,Delta_4)` (Jacobi's signature criterion);
4. require exactly one sign change for every simplex, hence exactly one
   timelike Gram direction;
5. compute every complex angle argument at 100 decimal digits and retain the
   minimum modulus over every evaluated simplex;
6. require that minimum modulus to exceed the original `1e-6` gate and the
   maximum action imaginary contamination to remain below `1e-70`.

The existing binary64 envelope and minimum eigenvalue are retained as a
separate independent audit; they are no longer evidence for resolving the
tiny derivative steps.

## Decision boundary

- If all high-precision signatures and angle-argument gates pass in both
  parities while the frozen action/momentum gates reproduce, restore
  `TWO_SLAB_GLUING_CONTROL_PASSED` as **DERIVED CONTROL**.
- If any high-precision branch gate fails, report
  `TWO_SLAB_DERIVATIVE_CONTROL_FAILED`; do not alter a step or tolerance.

No outcome changes the claim boundary: this is an action-composition control,
not an evolving next frame.
