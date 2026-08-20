# Correction protocol for the connection-norm assertion

Date: 2026-08-20

This correction is committed after the frozen `10/11` failure and before the
corrected verifier is executed.

## Permitted source change

Change only the `positive_ok` assertion.  It must require:

1. exact identity of the connection norm with
   `3*((lambda-1)^2+3*tau^2)`;
2. presence of the positive quadratic in its numerator factors;
3. an explicit real-domain certificate: if `tau != 0`, then
   `(lambda-1)^2+3*tau^2 > 0`.

It must not require `tau` to be a factor of the connection norm.  The
independent `tau != 0` hypothesis remains enforced through the recorded
vertex-solve determinants and the advertised real domain.

No formula, special point, ideal, factor whitelist, outcome hierarchy or
other control may change.  The frozen failure artifact remains immutable.
Only the corrected verifier may be rerun; the full suite remains out of
scope.

