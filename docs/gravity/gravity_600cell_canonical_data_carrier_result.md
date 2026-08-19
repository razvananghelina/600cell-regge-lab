# Canonical-data vertex-carrier result

Date: 2026-08-19

## Complete hypotheses

On the fixed lower regular 600-cell and the complete variable-face flat
Lorentzian frustum compatibility system, test the exact 240-column candidate
with one radial/scale and one normal/lapse variable at each upper vertex.  Use
the two rational nonstatic representatives `(lambda,tau)=(2,5),(3,11)`, two
exact right-inverse graphs, all complete face rows, and the preregistered
wrong-carrier and homogeneous controls.  No action, Hessian, propagation, or
physical unit is part of this test.

## Reproducible result

The final targeted execution of
`reproducible/verify_gravity_600cell_canonical_data_carrier.py` passed 10/10
checks and assigned

```text
CANONICAL_DATA_VERTEX_CARRIER_REFUTED
```

The final JSON artifact has SHA-256
`4065950aaac4180ec1cdd0b82f7a8bc403b2969c50d26cf14cc28592085cb2c5`.
For both representatives, the baseline and alternate right-inverse graphs
reject the candidate on exactly 3600 complete face rows and give the same
first exact residual.  The candidate data map itself has exact rank 240, the
wrong endpoint-difference and deleted-lapse controls are rejected, and its two
constant columns reproduce the frozen homogeneous controls.

## Verdict

**DERIVED NEGATIVE.** The compatible canonical-data carrier is not
`Q^120 vertex scale + Q^120 vertex lapse`.  Equality of dimensions was not
evidence of inclusion, and exact inclusion fails.

**OPEN.** The target-blind nullity 240 is stable over the two frozen primes,
representatives, and convention attacks, but its rational dimension and exact
geometric carrier remain unproved.  No physical phase-space, tensor mode,
tick, speed, `G`, or Planck scale follows from this negative.

The original 8/10 artifact and the intermediate 9/10 artifact remain preserved
in Git history together with the two verdict-harness corrections; neither
correction changed a scientific matrix entry or residual.
