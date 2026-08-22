# Frozen reporting correction: explicit Arb endpoints

Date: 2026-08-22.

First primary artifact commit: `7a40539`.

Status: frozen after the first primary run and before modifying or rerunning
the verifier.

## Observed problem

The first run returned `14/14` and every load-bearing comparison used Arb's
strict interval ordering.  However, the artifact serialized balls with
`str(ball)`.  Arb's low-relative-accuracy pretty printer may replace a
positive interval by a coarse display such as

```text
[+/- 1.03]
```

even when the internal `lower()` endpoint is strictly positive.  The display
is a valid outer enclosure but erases the sign evidence a reader needs to
audit the artifact.  The first artifact is preserved with SHA-256

```text
9eb252f18bd11361e6b4dffd4870feb40ad8a6c8283c876cc33ccf17877ed135.
```

## Frozen correction

Do not change the domain, formulas, precision, input artifacts, interval
objects, comparisons, controls or outcome hierarchy.  Change only artifact
reporting:

1. for every load-bearing Arb ball, store `lower()` and `upper()` separately;
2. store whether the ball contains zero;
3. for strict-positive gates, require both the existing comparison and an
   explicit `lower()>0` check;
4. add the already used implication that the accepted unimodality threshold
   lies above `q=1` when excluding negative physical roots;
5. preserve the original pretty string only as a diagnostic.

The corrected rerun must reproduce all mathematical values and the fifth-root
control.  Any sign disagreement changes the outcome to `OPEN`; no tolerance or
threshold may be adjusted.
