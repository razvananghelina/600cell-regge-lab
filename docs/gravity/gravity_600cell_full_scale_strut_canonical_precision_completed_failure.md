# First completed multiprecision run: control failure

Date: 2026-08-20  
Status: **DERIVED diagnostic; primary conclusion remains OPEN**

## Frozen execution

The first execution that completed both preregistered precision levels is frozen in
commit `bc343e117df1eccc921e3884ee4c2f488314fc4b`.

- verifier SHA-256:
  `8534dc54e81654035f2ce1a0a6c0f6dd75cad91f5768ad0f46b6a747d751e8d9`;
- artifact SHA-256:
  `3e6032b735a33eb44014ef11ff24a1c26566aa49cf655a98177c7c955e2bb142`;
- protocol SHA-256:
  `f044c0738fc7f507b89b1bc3658836ba5fa7a1d34f00f533bf821146663686b0`;
- execution result: `12/15` checks, outcome
  `FULL_SCALE_STRUT_CANONICAL_PRECISION_CONTROL_FAILED`.

The failed checks were the P160 geometry controls for both parities and the
aggregate precision hierarchy. The artifact did not serialize the individual
conjuncts of the geometry control. Therefore the exact failing conjunct is
**not known from this execution** and must not be inferred after the fact.

## What the failed artifact contains

These are diagnostics, not accepted conclusions because the encompassing
geometry control failed.

- **DERIVED diagnostic:** all nonhomogeneous interval Gram determinants exclude
  zero at P100 and P160; the binary classification does not disagree.
- **DERIVED diagnostic:** the smallest nonhomogeneous singular-value-to-ball-radius
  ratio is `2.8204006836328425e132`.
- **DERIVED diagnostic:** the largest nonhomogeneous P100/P160 relative change is
  `8.342450824032887e-35`.
- **DERIVED diagnostic:** the source/target reversal control has margin
  `8.247181898282127e137`.
- **OPEN:** in each parity the homogeneous D and K matrices have exactly one P160
  midpoint singular value below `1e-50`, while exhaustive drop-one minors certify
  ranks at least 9 and 14. The next singular values are approximately
  `9.8234441621e-8` for D and `8.0184489625e-8` for K.
- **OPEN:** the frozen P100 candidate evaluated without refitting at P160 has
  residual approximately `2.4898227e-42` in D and K. This does not pass the
  preregistered `< 1e-50` zero gate.

No statement about a physical or gauge mode, propagation, a tick, `c`, `G`, the
Planck scale, or particle mass follows from this failed run.

## Required next action

Add diagnostic serialization only: expose every existing geometry-control
conjunct and its raw value, without changing precision, matrices, thresholds,
candidates, classifications, or outcome logic. Freeze that execution before any
possible repair.
