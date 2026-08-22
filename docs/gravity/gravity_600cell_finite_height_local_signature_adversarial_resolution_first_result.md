# Monotone-factor resolution first result: one auxiliary control remains OPEN

Date: 2026-08-22.

Verifier commit: `1aa915f`.

Status: **8/10 OPEN; ROOT AND TREE DISAGREEMENT RESOLVED, AUXILIARY IDENTITY WIDTH UNRESOLVED**.

## Frozen result

```text
RESULT: 8/10 checks passed
OUTCOME: LOCAL_SIGNATURE_ADVERSARIAL_DISAGREEMENT_OPEN
```

Artifact SHA-256:

```text
70448c78be2156ef84fbaa986c543c6063bcca8ca4395ee77bdbf657ab2760d1.
```

## What passed

- exact protocol and input provenance;
- independent reconstruction of the accepted `K` factor in `p'`;
- positive and negative polynomial controls;
- all five fixed-bracket gravity states;
- complete stationary and all-real root counts;
- all non-diagonal physical and endpoint gates;
- the exact ordered `DEAD` and `ENTERED_D` tree;
- comparison with both the preserved first-adversarial endpoint signs and all
  primary root intervals and gate signs;
- rejection of the hostile `m*q>126` entry claim.

Thus the first wide-interval dependency problem was resolved without interval
Newton, subdivisions or discovery roots.

## Sole failure

All five evaluated balls for

```text
r-(1+h*q)
```

contain zero. Their widths are approximately

```text
3.23e-125,
5.25e-125,
3.19e-120,
1.06e-113,
5.63e-106.
```

The preregistered auxiliary threshold was `1e-110`; only the last recursively
propagated edge exceeds it. Consequently the identity check and final outcome
had to fail even though no ball excludes zero.

## Permissible resolution

The width threshold cannot be loosened after seeing the result. A separate
protocol may instead certify the exact algebraic identity

```text
r-(1+h*q)=-E/(2*pi*mu(q)).
```

At every already certified unique root `E=0`, this proves the endpoint identity
without a numerical width requirement. The resolver must also verify
`mu(q)>0` and include a hostile false identity that is strictly rejected.

Until that separate gate passes, the primary local theorem remains
unconsolidated despite agreement of the complete roots and tree.

