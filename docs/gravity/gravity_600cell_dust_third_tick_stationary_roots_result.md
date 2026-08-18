# Third-tick stationary roots — NO target comparison

Date: 2026-08-16

## Provenance

- prior-art gate: `7b9a676`;
- frozen protocol: `5d980b1`;
- registered implementation before evaluation: `addbac6`;
- targeted verifier:
  `reproducible/verify_gravity_600cell_dust_third_tick_stationary_roots.py`;
- artifact:
  `reproducible/gravity_600cell_dust_third_tick_stationary_roots.json`;
- artifact SHA-256:
  `02d4589a7df0851c67a31fc0a41c5ef8851a82c758214c1c5e8729afddfe479f`.

Only this targeted verifier was run.  It returned **5/5**.  The full suite was
not run.

**No accepted second-tick post-momentum was JSON-parsed or compared.**  The
verifier only byte-hashed the accepted second-tick artifact and serialized
`target_parsed=false`.

## Mechanical verdict

```text
THIRD_TICK_STATIONARY_ROOTS_ENUMERATED
N=2
```

On the frozen 257-point main grid

```text
x=C/abs(A1) in [-8,8], spacing=1/16,
```

there is one standalone sign-changing bracket and one near-zero node cluster.
The sign bracket touching the node cluster was merged exactly as
preregistered.

## Complete root multiset

### Root 0 — forward contracting candidate

```text
x = -6.00004539832334081494468046598,
C = -1.869649892404722217034061625667e-5,
bracket width = 1.6109650e-31,
G = -1.6853524e-33,
P = +0.00454045110614718387402634920808,
D = -2.4739258487162612e-6,
epsilon_D = 1.3416802e-28.
```

### Root 1 — exact time reversal of tick two

```text
x = -1,
C = A1 = -3.11605957669450169173470644420e-6,
G = -1.9652587e-48,
P = -0.00454044949003343733024916509405,
D = -2.4739258487462894e-6,
epsilon_D = 2.4434509e-29.
```

Both roots pass the complete 35 internal equations, derivative gates,
Lorentzian branch and independent schedule-parity comparison.  The two
determinants have resolved negative sign in both parities.

## Preregistered prediction audit

- predicted sampled candidate count `N=2`: **survives**;
- structural time-reversal root at `C=A1`: **exactly recovered**;
- pattern contracting root near `C=6*A1`: **survives**, with
  `C/abs(A1)=-6.0000453983`.

This agreement did not use the canonical momentum target.  It is stronger
than fitting a root after seeing the seam, but remains **PATTERN**, not a
recurrence theorem.

## Scope ledger

- **DERIVED COMPUTATIONAL ON THE FROZEN DOMAIN:** exactly two sampled
  sign/node candidates and their complete `(C,G,P,D)` data.
- **DERIVED:** root 1 is the time reversal of the accepted second slab.
- **DERIVED COMPUTATIONAL:** root 0 is locally regular rather than a fold;
  `|D|/epsilon_D` is about `1.84e22`.
- **OPEN:** even-multiplicity tangential roots between grid nodes.
- **OPEN:** roots outside `x in [-8,8]`; the sentinels are diagnostics only.
- **NOT SELECTED:** no root is yet declared the third physical tick.
- **NOT PARSED:** the desired third-tick momentum target.

## Firewall handoff

This note and artifact must be committed before any target comparison.  A
later verifier must compare **both** roots with the mapped accepted
second-tick post-momentum and report `hits/2`.  It may not modify either root
or silently discard the time-reversal candidate.
