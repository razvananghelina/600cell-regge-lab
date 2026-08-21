# Result: the surviving history has a unique physical fourth slab

Date: 2026-08-21.

## Provenance

```text
fourth-slab prior-art gate and protocol            54c4554
primary verifier registered before first run       a58c46e
accepted primary artifact                          7601c8f
dual adversarial protocol                          530b2f0
adversarial verifier registered before first run   e217807
accepted adversarial artifact                      2a6b78b
```

Targeted verifiers:

```text
reproducible/verify_gravity_600cell_finite_height_fourth_slab.py
  7/7 PASS

reproducible/verify_gravity_600cell_finite_height_fourth_slab_adversarial.py
  9/9 PASS
```

Accepted artifacts:

```text
reproducible/gravity_600cell_finite_height_fourth_slab.json
SHA-256 cf322cf0d60668d8f3f58e251425c9ad6bf43b112f22f9f3aebbc28f86212468

reproducible/gravity_600cell_finite_height_fourth_slab_adversarial.json
SHA-256 ac1ed19fd72549cf7cd054107d921e2819580704391fbc294a55e106a8f7a1bd
```

Both targeted verifiers were rerun together and reproduced these hashes byte
for byte.  No full-suite run was performed.

A static AST audit of `reproducible/run_all.py` found 415 verifier files on
disk, 413 distinct registrations, two deliberate exclusions with reasons,
zero duplicate registrations, zero unregistered files and zero registered
files missing from disk.  This audit executed no verifier.

## Headline

```text
SURVIVING_HISTORY_HAS_UNIQUE_FOURTH_SLAB_
ADVERSARIALLY_CORROBORATED
```

> **DERIVED COMPUTATIONAL, FOUR-SLAB SCOPED / STRUCTURAL / ADVERSARIALLY
> CORROBORATED:** the unique physical three-slab history selected by one-step
> future extendibility has exactly one physical fourth-slab continuation.
> Finite-horizon uniqueness therefore persists for one additional step.

This remains a finite statement.  It neither proves an infinite history nor
turns the dimensionless state-dependent intervals into a universal or
fundamental unit of time.

## Complete hypotheses

Use the fixed homogeneous tetrahedral-frustum 600-cell action at zero
cosmological constant, conserved global dust, positive proper heights,
positive endpoint scales and the committed pre/post momentum convention.
The frozen history is

```text
v=3/2 -> unique first slab -> second branch B -> unique third slab.
```

No conclusion covers other initial states, nonhomogeneous data, refinement,
another action, a cosmological constant or a quantum boundary condition.

## Fourth incoming state

Direct reconstruction of all three preceding slabs gives

```text
m3  = 0.395744374852478801317435698526509...,
pi3 = -171.686225162332253486221608644291....
```

The recurrence

```text
m_next=m/r,
pi_next=p_post/r^2=p(q)+2*pi*h*mu(q)/r
```

agrees directly through slab three.  In the primary 170-digit calculation,
the third recurrence discrepancy rounded to exact zero at the working
precision.  Wrong momentum scaling, reversed sign and mass reset all change
the fourth incoming state.

## Complete fourth-root census

Both proofs find three real algebraic roots and exactly one physical root:

| `q4` | `h4` | `L4/L3` | Physical failure |
|---:|---:|---:|---|
| `-0.2137436416...` | `8.439306735...` | `-0.803848154...` | negative endpoint scale |
| `99.627601694...` | `-0.006876956864...` | `0.314865281...` | negative height |
| `316.698862258...` | `0.006876957360...` | `3.177924572...` | none |

The physical root passes the complete-action lapse, incoming-momentum and
shared-slice junction residuals below `7e-131`.  Its local determinant is
strictly positive by the exact determinant theorem.

## Independent proof structure

The primary route partitions the real line by solutions of

```text
p(q)=pi3
```

using `E4'(q)=p(q)-pi3`.  The adversarial route independently reconstructs
the entire history at 110 and 180 digits, solves the constraint first and
partitions by

```text
mu(q)=m3,
R4'(q)=4*pi[m3-mu(q)]/q^2.
```

It treats both one-sided zero limits, both infinite tails and `q=0`, and
requires all roots and labels to nest beyond 60 digits before reading the
primary artifact.  The two routes then agree beyond the stored precision.

## Status ledger

| Claim | Status |
|---|---|
| The frozen three-slab history has a fourth continuation | **DERIVED COMPUTATIONAL / ADVERSARIALLY CORROBORATED** |
| That fourth continuation is unique | **DERIVED COMPUTATIONAL, four-slab scoped** |
| Finite-horizon uniqueness persists from slab three to slab four | **DERIVED STRUCTURAL** |
| The history has a fifth slab | **OPEN** |
| The history extends indefinitely | **OPEN** |
| The apparent large-`q` regime is asymptotically self-similar | **PATTERN / OPEN theorem** |
| The absolute proper slab height converges or is selected | **PATTERN / OPEN theorem** |
| The result persists for other first states | **OPEN** |
| A deterministic fundamental tick is derived | **NO** |

## Interpretation and next gate

The successive physical slopes and scale ratios now show a conspicuous
large-`q` pattern:

```text
q2=31.279...,
q3=99.628...,
q4=316.699...,

r2=3.156...,
r3=3.176...,
r4=3.178....
```

This is **PATTERN**, not yet a theorem.  Brute-force addition of one slab at a
time is no longer the most informative next step.  The next target-free gate
is to derive the large-`q`, small-mass canonical map, identify any fixed point
in scale-free variables, and prove or refute that it controls an infinite
unique branch.  That calculation must also determine whether physical proper
heights converge only relative to the initial scale or supply any new time
scale.  Global scale covariance forbids the latter without an additional
dimensionful input.
