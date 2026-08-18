# Extended preregistration: 1,024-moment refined Whitney complexity

Date: 2026-08-11

## Provenance

This is an explicitly **post-result**, non-blind extension of protocol
`366fe4a` and result commit `8ac7b62`.  The first refined Krylov complexities
were `(34,60,70,1)`, and the first three exactly saturated half of their
frozen sequence lengths `(68,120,140,36)`.  Those are valid lower bounds but
are censored by the observation window.

The purpose of this extension is fixed before computing longer sequences:
distinguish modest growth just beyond the first bounds from substantially
larger refined spectral complexity.  No stopping rule based on an attractive
number is allowed.

## Frozen extension

Retain without change:

- the exact coarse and first-barycentric Whitney mass matrices;
- probe seed `60020260811 + form_degree`;
- primes `1000003`, `1000033`, `1000037`;
- the same integer sparse multiplication and Berlekamp--Massey code.

For form degrees 0, 1 and 2, compute exactly **1,024** scalar Krylov moments
for each prime.  For the scalar top-form control retain 36 moments; its
primitive mass is already certified as the identity.

Calibrate again on the coarse blocks: longer sequences must still return
exact complexities `(9,22,27,1)` for all three primes.

For each refined block report:

- all three extended complexities;
- their maximum as a rigorous minimal-degree lower bound;
- whether the value equals 512, the ceiling of a 1,024-term generic
  Berlekamp--Massey prefix;
- the corresponding inverse-degree lower bound (L-1);
- its ratio to the exact coarse inverse degree `(8,21,26)`.

## Labels

- **DERIVED EXTENDED LOWER BOUND:** every modular complexity is exact for its
  recorded sequence, regardless of whether it stabilizes.
- **CENSORED:** a value of 512 is a lower bound only; do not report it as the
  refined minimal degree.
- **STABILIZED CANDIDATE:** a common value below 512 across all three primes
  is still not an exact matrix degree until a whole-matrix annihilator is
  certified.  Label it candidate, not equality.

This extension does not prove divergence and does not justify a scaling fit
from two levels.  Only the targeted verifier will be run.
