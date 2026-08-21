# Control correction: refined local curvature mass verifier

Date: 2026-08-21

First-failure commit: `ead1a1c`.

The registered first execution returned `13/15` and
`REFINED_LOCAL_CURVATURE_MASS_CONTROL_FAILED`.  The failed artifact remains
preserved.  This note also corrects one sentence in the first-failure report:
the 64-character hash was already correct in the frozen protocol; only the
verifier's copied constant omitted `e0`.

Before another scientific execution, make only these two changes:

1. restore the omitted `e0` in the verifier's feasibility-source SHA-256,
   giving the actual unchanged digest
   `36fba835048e6e0f0676b749192a9d882406932770a00ba1396929bbc4d04a32`;
2. store `tau0` as text and convert it to an `mp.mpf` only after entering the
   100-decimal context.

The identity, `1e-68` tolerance, curvature allocation, alternative-mass and
corruption controls, outcome hierarchy and all input contents remain frozen.
No exploratory number is used to loosen a gate.
