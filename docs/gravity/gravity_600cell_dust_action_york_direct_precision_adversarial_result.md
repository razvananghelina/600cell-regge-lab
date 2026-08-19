# Adversarial result: a scale-dependent classifier contradicts itself

Date: 2026-08-19

The mechanically independent audit passed `10/10` checks and, under its
frozen hierarchy, returned

```text
ADVERSARIAL_DIRECT_REFUTATION_REFUTED.
```

That label is preserved exactly in commit `33e0fdc`.  It is not silently
rewritten after seeing the result.

## What the two observables actually did

For all sixteen schedule/sector/variant cells:

```text
r_span = ||(I-Q_BL Q_BL*) A L||       NONZERO_RESOLVED
r_comm = ||(I-L L*) B^-1 A L||        ZERO_CONSISTENT
rank[B L, A L] = 24                    (not 15)
```

The operational-primary values are essentially constant:

| quantity | value | empirical absolute error |
|---|---:|---:|
| `r_span` | `2.9633155952e-5` | `3.1364061730e-8` |
| sixteenth augmented singular value | `2.9633155950e-5` | rank threshold derived from the same floor |
| `r_comm` | `5.4834883103e-8` | `3.1364061730e-8` |

Thus the independent QR/span construction confirms that `A L` is not
contained in `B L`: the augmented image has rank 24 and a resolved sixteenth
singular direction.  The `B^-1 A` leakage is also nonzero at its midpoint, but
the protocol compared it to the same dimensionful absolute error used before
the inverse rescaled the operator.  It consequently fell below
`10*error`.

## Framing failure

Exact arithmetic gives the equivalence

```text
A L subset B L    iff    B^-1 A L subset L
```

because `B` is invertible.  The two numerical tests therefore cannot
legitimately support opposite mathematical verdicts.  The contradiction is
in the error normalization, not in the linear algebra.

The preregistered hierarchy nevertheless declared the primary refutation
"refuted" whenever either new residual was zero-consistent.  That implication
was too strong: zero-consistent under a loose or badly scaled bound is not
evidence of exact zero, especially when an equivalent residual and an
augmented-rank witness are resolved nonzero.

Status:

- the literal adversarial outcome is DERIVED COMPUTATIONAL and preserved;
- its use as a refutation of the primary result is REFUTED by the audit's own
  mutually incompatible observables;
- the independent span/rank evidence for pseudo-longitudinality is DERIVED
  COMPUTATIONAL / STRUCTURAL;
- a scale-invariant adversarial confirmation remains OPEN until preregistered
  and run.

The correction must normalize each residual by the norm of the operator it
tests and use a dimensionless backward-error floor.  It may not delete or
replace this failed protocol or artifact.
