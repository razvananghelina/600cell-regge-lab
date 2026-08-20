# Preregistered correction: nonlinear time-reversal anchors

Date: 2026-08-20

First failure commit: `d2c018a`.

The first execution stopped at `10/11` before every scientific solve because
the two displaced anchors `A1,A2` left the real-curvature Lorentzian branch:
their minimum angle-argument moduli remained nonzero and their angle identities
passed, but their maximum imaginary curvatures were order one. Both members
of every time-reversal pair failed identically. The control incorrectly
required both anchors to remain branch-valid but chose displacements far
larger than the thin local branch around the induced geometry.

## Frozen narrow repair

Keep `A0=0`. Replace only the two diagnostic anchors by

```text
A1' = 1e-4 A1
A2' = 1e-4 A2.
```

Thus

```text
A1'=(1e-6,-1e-6,2e-6,-2e-6,3e-6,-3e-6,
     -1e-4,-1e-4,-1e-4,-1e-4),
A2'=(-3e-6,2e-6,-1e-6,1e-6,-2e-6,3e-6,
     -5e-5,-1e-4,-1.5e-4,-2e-4).
```

This deterministic scaling is chosen from continuity at the already certified
branch-valid base point, not by searching a passing amplitude. At 80 decimal
digits it remains a nonzero nonlinear substitution. Require the original
`1e-60` equality threshold and identical valid branch diagnostics.

Add per-anchor diagnostics to the JSON so another invalid anchor cannot be
hidden behind an infinity.

No scientific seed, bound, solver option, acceptance gate, class count,
look-elsewhere denominator or outcome changes. If either corrected anchor is
still invalid, preserve a second control failure; do not shrink it again in
the same execution.
