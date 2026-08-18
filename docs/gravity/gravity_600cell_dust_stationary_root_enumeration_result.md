# Stationary-root enumeration result — NO target comparison

Date: 2026-08-16

## Provenance

- prior-art gate: `eecc80e`;
- frozen enumeration protocol: `07083cc`;
- initial implementation: `0c48e22`;
- counter-reporting failure: `7d27a20`;
- reporting-only fix: `8c31ec8`;
- corrected artifact SHA-256:
  `0ec5ba520ea25b39dd6cfd3c349d49fe480df2abee359854e1316b5af4d9fa2f`.

Only the targeted verifier was rerun.  It returned **5/5**.  The full suite
was not run.

**No desired second-tick momentum was parsed or compared in this step.**

## Mechanical verdict

```text
STATIONARY_ROOTS_ENUMERATED
```

On the frozen 257-point main grid

```text
x=b/|a1| in [-8,8],  spacing=1/16,
```

the target-independent equation `G(b,r1)=0` has exactly two enumerated
candidates:

```text
one sign-changing bracket,
one near-zero grid-node cluster,
N_candidates=2.
```

The 18 dyadic sentinels through `x=+-4096` remain Lorentzian and have positive
`G`; they are diagnostic only and do not strengthen the finite-domain root
count into an analytic theorem.

## Complete root multiset

### Root 0 — contracting stationary branch

```text
x = -3.00001691005308401262091769995,
b = -9.34823142281635992512450895990e-6,
bracket width = 1.6109650e-31,
G = -1.1975861e-34,
P = +0.00272430041185858338621233110987,
D = -2.4739522649094586e-6,
epsilon_D = 8.4385063e-29.
```

### Root 1 — time-reversal branch

```text
x = 0,
b = 0,
G = +4.2574745e-31,
P = -0.00272429879575512949578836966717,
D = -2.4739522649274752e-6,
epsilon_D = 1.0779799e-29.
```

Both roots pass the complete 35 internal equations, Lorentzian branch,
resolved nonzero `G_b`, resolved determinant sign and independent even/odd
derivative comparison.

## Status ledger

- **DERIVED COMPUTATIONAL ON THE FROZEN DOMAIN:** exactly two sampled
  sign/near-zero roots and their full characters `(b,G,P,D)`.
- **STRUCTURAL:** the `b=0` solution is the time reversal of the accepted first
  slab and carries the opposite pre-momentum sign.
- **DERIVED COMPUTATIONAL:** the contracting root is locally regular, not a
  fold; its determinant is about 22 orders of magnitude above its calibrated
  error.
- **OPEN:** tangential even-multiplicity roots between grid nodes and roots
  outside `x in [-8,8]`.
- **PATTERN:** positive sentinel values through `|x|=4096` suggest no further
  outer sign-changing roots, but do not prove it.
- **NOT SELECTED:** neither root is declared the physical branch in this step.

## Firewall handoff

This note and the full root artifact must be committed before the desired
second-tick momentum is loaded.  The later comparison may test both roots; it
may not discard the time-reversal root or adjust the contracting root.
