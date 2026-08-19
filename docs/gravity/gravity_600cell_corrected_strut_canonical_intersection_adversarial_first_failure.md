# First adversarial intersection execution: JSON boolean type failure

Date: 2026-08-19

The preregistered adversarial verifier from commit `b5739c8` completed every
scientific calculation and printed:

```text
even: all seven QR/Frobenius lower bounds certify zero intersection  PASS
even: synthetic, basis and convention controls pass                 PASS
odd:  all seven QR/Frobenius lower bounds certify zero intersection  PASS
odd:  synthetic, basis and convention controls pass                 PASS
role corruption changes an adversarial QR lower bound                PASS
```

It then failed before writing an artifact because `json.dumps` received a
NumPy boolean scalar in the nested payload:

```text
TypeError: Object of type bool is not JSON serializable
```

The process exited `1`.  No adversarial JSON artifact was created and no
scientific result is accepted from this execution.

Authorized repair: cast the uncertainty's binary and total terms to builtin
`float`, and cast stored comparison/control flags to builtin `bool`.  Also
cast the per-parity displayed minimum ratio to numeric `float` before taking
its minimum; the first run's displayed `2.3895e12` used lexicographic string
ordering and was diagnostic text only.  No input, matrix, QR decomposition,
threshold, control or outcome rule may change.

