# Receipt: unchanged direct-action adversarial rerun

Date: 2026-08-21

Corrected-adjudication protocol commit: `2a10336`.

After freezing that protocol, the unchanged registered verifier

```text
/home/razvan/science/.venv/bin/python \
  reproducible/verify_gravity_600cell_refined_h4_constrained_response_adversarial.py
```

was executed to completion a second time.  As preregistered, it exited nonzero
with the unchanged historical summary

```text
tests 15/17
direct classes 1
primary matches 24/24
outcome ADVERSARIAL_REFINED_H4_CONSTRAINED_RESPONSE_CONTROL_FAILED.
```

The regenerated artifact is byte-identical to the first complete run:

```text
reproducible/gravity_600cell_refined_h4_constrained_response_adversarial.json
SHA-256 a23ef4cc23d08ad8768f1df66789aa900cdb95a7f3529486df80697a53b1fe81.
```

The second execution also reproduced the displayed direct response envelope
`1.3984e-54`, internal minimum eigenvalue `1.3780099e-5`, normalized solve
residual `1.7954044e-181`, 220-digit fraction `6.8529575e-91`, cross-method
match count `24/24`, maximum cross-method fraction `6.6666667e-5`, and
corruption factor `1.9084791e24`.

This receipt establishes deterministic reproducibility of the full frozen
direct dataset.  It does not change the failed verifier's outcome and is not
an independent derivation.

