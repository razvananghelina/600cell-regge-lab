# Threshold correction: adversarial H4 null coupling

Date: 2026-08-21

The formal disagreement and diagnosis are preserved in commit `f89d2d0`.

For primary schedule record `j`, replace only

```text
|c_adversarial-c_primary,j| < 1e-68
```

by

```text
|c_adversarial-c_primary,j| <= envelope_j,
```

where `envelope_j` is read directly from the frozen primary artifact and was
computed before the adversarial reconstruction.  Require all 24 envelopes to
be positive and the maximum error/envelope ratio to be `<1`.

Keep the actual-incidence boundary-vector comparison at `<1e-68`.  Do not
change the derived `1/8`, inputs, component order, vertical identity,
rank/reversal criteria, corruption thresholds, scope or outcome hierarchy.

Run only the corrected adversarial verifier twice and require a byte-identical
artifact.

