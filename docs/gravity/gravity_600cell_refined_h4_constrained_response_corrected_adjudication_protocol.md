# Protocol: corrected adjudication of the direct constrained H4 audit

Date: 2026-08-21

Status: frozen before the required unchanged rerun and before constructing the
corrected adjudication artifact.

## 1. Why this is a separate gate

The first direct-action adversarial verifier passed every decisive response
test but formally returned `15/17 CONTROL_FAILED`.  A later preregistered
diagnostic established, twice with a byte-identical artifact, that:

1. the fourth-order stationarity probe was under-resolved; and
2. termwise reality of individual Lorentzian off-shell curvatures was an
   overstrong category error, while the complete action stayed real on a
   continuous safe branch.

The failed execution remains failed.  This gate may adjudicate its frozen
direct data using the independently frozen diagnostic; it may not rewrite the
historical outcome, regenerate matrix entries with altered tolerances, or
discard any decisive comparison.

## 2. Required unchanged reproducibility rerun

Before implementing this gate, execute the unchanged verifier

```text
reproducible/verify_gravity_600cell_refined_h4_constrained_response_adversarial.py
  78f6b52f6f019a150a86ddadcb819b67c3757244c015687ab67f4649784ac53d
```

once more.  It is expected to exit nonzero with the same `15/17` historical
control failure.  Require its overwritten JSON to be byte-identical to the
first run:

```text
reproducible/gravity_600cell_refined_h4_constrained_response_adversarial.json
  a23ef4cc23d08ad8768f1df66789aa900cdb95a7f3529486df80697a53b1fe81.
```

Any different hash kills corrected adjudication as a reproducibility control
failure.  This second run checks the complete 24-schedule, 210-direction data,
not only the two 220-digit spot checks.

## 3. Frozen diagnostic evidence

Require exact hashes and scoped outcomes for

```text
reproducible/verify_gravity_600cell_refined_h4_constrained_response_auxiliary_diagnostic.py
  70beeffe19cf4b6e90a613d3936f9c30bd98021e0a7b6ae6b7e93d60c01c0bc4
reproducible/gravity_600cell_refined_h4_constrained_response_auxiliary_diagnostic.json
  f66177326afc3b3457a60b544745b739cbaa6b6d6e7f367b57d60f31eeeddeb7
docs/gravity/gravity_600cell_refined_h4_constrained_response_auxiliary_diagnostic_protocol.md
  2f6d9d72e04c4baf1dc385425ef7f26ba0a55f6249d6505db111aa21e0836405
docs/gravity/gravity_600cell_refined_h4_constrained_response_auxiliary_diagnostic_result.md
  afd3bb5bcec476bded8ea003c5749a83fa46c488b0d1c7da3d693da93fbe9423
```

Require diagnostic outcome

```text
REFINED_H4_CONSTRAINED_RESPONSE_AUXILIARY_FAILURES_RESOLVED,
tests 13/13,
```

all 240 stationarity zero gates at both precisions, all 240 precision gates,
all `5040/5040` curvature parity pairs, the two halving ratios in
`[1.99,2.01]`, and the complete-action analytic branch gates.  Do not replace
these with summaries from prose.

## 4. Independent re-adjudication of frozen matrices

Load all 24 stored direct `11 x 11` response matrices and their direct
envelopes from the reproducible failed artifact.  Also require every stored
`20 x 20` restricted second variation, `9 x 11` lift, nine internal
eigenvalues, positive internal minimum, direct solve residual, action
displacement, two 220-digit repeat gates, polynomial control and corruption
control.

Reconstruct `P=E(c,3)` from the accepted compatibility covector and derive
`T_R` from `R P=P T_R`; do not trust the stored class labels.  Recompute:

```text
K_o = T_R^T K_rev(o) T_R,
K_can,o = K_o                         if o <= rev(o),
          T_R^T K_o T_R               otherwise,
```

with the exact max-entry envelope propagation frozen in the direct protocol.
Rebuild the complete lexicographic class census using the sum of the two
direct envelopes as the equality gate.

Only after this target-free census, load the frozen primary artifact

```text
reproducible/gravity_600cell_refined_h4_constrained_response_corrected.json
  85adea23f6a19153f61f3ed066137a5e40ab77b8901d4cc81cfc4f864e0bc093
```

and independently recompute all 24 cross-method comparisons with

```text
||K_direct-K_primary||max
 <= 10(e_direct+e_primary).
```

Repeat the one-entry corruption check from the stored uncorrupted matrix and
require it both to split a two-matrix class census and to fail the primary
gate.  No stored boolean or stored class count is sufficient evidence for any
of these recomputations.

## 5. Framing attack

- This gate is a corrected adjudication, not a third mechanically independent
  construction.  Mechanical independence is supplied by the frozen direct
  scalar-action reconstruction; reproducibility is supplied by its required
  unchanged rerun.
- The auxiliary diagnostic changes only the status of two falsified controls.
  It cannot rescue multiple direct classes, failed time reversal, failed
  precision repeats or disagreement with the primary Hessian route.
- A corroborated result establishes schedule-independence only in the finite
  homogeneous `H4` constrained linear sector.  It does not establish a
  nonlinear tick or gravitational-wave propagation.

## 6. Frozen outcomes

Use the first applicable outcome:

1. `CORRECTED_ADJUDICATION_REFINED_H4_CONSTRAINED_RESPONSE_CONTROL_FAILED` if
   provenance, unchanged-rerun hash, diagnostic, dimensions, precision or
   corruption controls fail.
2. `CORRECTED_ADJUDICATION_REFINED_H4_CONSTRAINED_RESPONSE_DISAGREEMENT` if
   controls pass but direct time reversal fails, the recomputed direct census
   has more than one class, or any direct matrix fails its primary comparison.
3. `CORRECTED_ADJUDICATION_REFINED_H4_CONSTRAINED_RESPONSE_CORROBORATED` only
   if all controls pass, all 24 direct matrices form one class, reversal is
   covariant and every direct matrix matches the primary route.

Run the corrected adjudicator twice and require byte-identical JSON.  Run no
full suite, root search or deferred nonlinear census.

