# Preregistration: stationary-root enumeration at the inherited lapse

Date: 2026-08-16

Prior-art gate: `eecc80e`.

Nonmonotone base result: `f9d7ada`.

Status: **frozen before evaluating the enumeration grid**.

## 1. Target firewall and fixed equation

Retain the accepted artifact SHA-256 firewall and the committed geometry
constants `(a1,r1)` from `ed1cd6a`.  Do not parse any momentum target.

At fixed lower log `a1` and lapse log `r1`, enumerate zeros of

```text
G(b)=mean(five complete pole equations).
```

Record `P(b)=mean(thirty pre-momenta)` only after roots are isolated; it is not
a root criterion.

## 2. Frozen enumeration domain and count

Set `x=b/abs(a1)` and evaluate exactly

```text
x_k = -8+k/16,  k=0,...,256.
```

Thus the main grid contains exactly `N=257` points on `[-8,8]` with spacing
`1/16` in `x`.  Also evaluate, for diagnostics only, the 18 dyadic sentinels

```text
x = +/-2^n,  n=4,...,12.
```

Sentinels do not enter the root count or isolation brackets.

Every point must record `b`, `G`, `P`, branch flags and branch margins.  A
branch failure is an explicit outcome; no point is skipped.

## 3. Root candidates and exact count

Before refinement form:

- one node-root candidate for every main-grid node with `|G|<1e-25`;
- one sign-bracket candidate for every adjacent pair whose endpoints both
  have `|G|>=1e-25` and resolved opposite signs.

Consecutive near-zero nodes form one node-root cluster, represented by their
middle grid node.  Sign brackets adjacent to a near-zero cluster are excluded.

The exact candidate counts and full sampled multiset must be written before
refinement.  More than eight total candidates is an enumeration-open outcome.

## 4. Refinement

For every sign bracket perform exactly 80 bisections, retaining the half with
opposite signs.  Require final `b` width `<1e-29` and midpoint `|G|<1e-25`.

For every node-root use the committed grid node unchanged and require
`|G|<1e-25`; do not Newton-correct it.

No alternate grid, target, Newton root search, interpolation or seed is
allowed.

## 5. Complete root certificates

At every distinct candidate root require:

```text
Lorentzian/complex-angle branch passes,
max abs(30 diagonal equations) < 1e-60,
max abs(5 pole equations)      < 1e-25,
within-type spreads            < 1e-60.
```

Compute the four calibrated Jacobians of `(G,P)` in `(b,r)` with the same
steps and entrywise gates as `ed1cd6a`.  Require `G_b` resolved nonzero; record
the determinant and whether its sign is resolved, but do not require the full
Jacobian to have rank two because a root may itself lie at a fold.

Substitute every even root into the odd schedule and require `G`, `P`, all 35
residuals and all derivative entries to agree within their combined calibrated
errors.

## 6. Mechanical outcomes

Assign the first applicable outcome:

1. `STATIONARY_ROOT_ENUMERATION_CONTROL_FAILED`;
2. `STATIONARY_ROOT_ENUMERATION_BRANCH_FAILED`;
3. `STATIONARY_ROOT_ENUMERATION_TOO_MANY_CANDIDATES`;
4. `STATIONARY_ROOT_ENUMERATION_REFINEMENT_FAILED`;
5. `STATIONARY_ROOT_ENUMERATION_DERIVATIVE_OPEN`;
6. `STATIONARY_ROOT_ENUMERATION_PARITY_FAILED`;
7. `STATIONARY_ROOTS_ENUMERATED` only if every frozen gate passes.

## 7. Interpretation boundary

The result enumerates sign-changing and sampled near-zero roots on the frozen
main grid.  It does not rule out tangential roots between nodes or roots
outside `x in [-8,8]`; the sentinels are **PATTERN** only.

No root is selected as the physical forward branch in this step.  The complete
root list and momentum multiset must be committed before any branch trace or
target comparison.

Only the new targeted verifier is run; the full suite remains excluded.
