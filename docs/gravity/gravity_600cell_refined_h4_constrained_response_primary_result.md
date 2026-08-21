# Corrected primary result: one constrained H4 response class

Date: 2026-08-21

Status: clean primary result; not accepted before mechanically independent
replication.

## Provenance and correction ledger

| stage | commit |
|---|---|
| prior-art gate, final clarification | `8ecbd2a` |
| frozen primary protocol | `be10390` |
| verifier registered | `69ace62` |
| first synthetic-control crash preserved | `cddb3ca` |
| synthetic repair preregistered / implemented | `2ed8157`, `9b2436b` |
| first `18/19` primary failure preserved | `ba828ec` |
| directional diagnostic protocol / registration | `0baeeee`, `42b3b88` |
| diagnostic implementation failures and repairs | `dd7e05a`--`62abdd6` |
| corrected diagnostic result | `956ccf2` |
| combined-envelope protocol / registration | `5cefff6`, `711105b` |
| combined-envelope result | `3589610` |
| primary directional correction, final amendment | `7f2bfa9` |
| corrected primary implementation | `e281c57` |

The failed artifact remains frozen at its original path and hash.  The
corrected verifier passed `19/19` twice and produced the byte-identical new
artifact

```text
reproducible/gravity_600cell_refined_h4_constrained_response_corrected.json
SHA-256 85adea23f6a19153f61f3ed066137a5e40ab77b8901d4cc81cfc4f864e0bc093.
```

No full suite or deferred nonlinear census was run.

## Primary result

For every one of the 24 staircase schedules:

- the complete `22 x 22` Hessian passes its two-step/two-precision envelope;
- the analytic product tangent reproduces the isolated internal null line;
- the frozen compatibility row is reproduced;
- the restricted `9 x 9` internal complement is positive, with global minimum
  eigenvalue `1.3780099e-5`;
- the constrained solves and full internal residuals pass;
- changing the boundary pivot changes the matrix by congruence, and changing
  the internal complement leaves the bilinear form unchanged.

The complete target-free class census gives

```text
class count = 1,
membership  = all 24 schedules,
time reversal covariance = true.
```

The two basis controls differ by at most `2.752e-135` and `5.751e-133`.
The 12 high-order complete-action checks have exact target-independent
step-halving factors and agree with the response inside their combined
action/Hessian envelopes; the maximum used fraction is `0.00066391091`.

The frozen outcome is

```text
REFINED_H4_CONSTRAINED_RESPONSE_SINGLE_SCHEDULE_CLASS.
```

## Meaning and limits

- **PRIMARY DERIVED COMPUTATIONAL:** the constrained linearized
  boundary-momentum bilinear form on `ker(c^T)`, modulo the conormal `c`, has
  one class across all 24 bare staircase schedules in the `H4` sector.
- **PRIMARY DERIVED STRUCTURAL:** this response is independent of the two
  frozen choices of boundary basis and internal null complement.
- **DERIVED NEGATIVE:** the result is not an unconstrained Schur complement
  and does not license a Moore--Penrose extension away from `ker(c^T)`.
- **OPEN:** mechanically independent replication; until it passes, the
  single-class result is not accepted.
- **OPEN:** integration of the linear hyperplane to a nonlinear admissible
  boundary family.
- **NOT ESTABLISHED:** nonhomogeneous propagation, a tick, a dispersion
  relation, `c`, `G`, Planck units or particle physics.

Even after replication, a positive result advances only to building the
nonhomogeneous quadratic operator.  The invariant sector alone has no spatial
wavenumber and cannot determine a propagation speed.

