# Adversarial protocol: direct-quotient half-strip certificate

Date: 2026-08-22.

Parent protocol commit: `559ca7a`.

Primary corrected artifact commit: `f5d8d0f`.

Status: frozen before the adversarial verifier exists or is executed.

## 1. Unchanged theorem and scope

Test exactly the same post-hoc half-strip

```text
0<m<=2/5,
x>=125
```

under the same homogeneous 600-cell action, zero cosmological constant,
conserved global dust, momentum convention and physical inequalities.  No
threshold, branch or outcome may change.

## 2. Independence from the primary certificate

The primary proof encloses removable quotients through exact integral Taylor
remainders on one full `u` interval.  The adversarial route must not use those
load-bearing integral enclosures.

Instead it must:

1. redifferentiate the complete one-slab action with new symbols;
2. use the unit-square coordinates

   ```text
   a=(5*m/2)^2,
   b=125/x,
   ```

   while evaluating the regular scalar

   ```text
   u=m^2/x^2=4*a*b^2/390625;
   ```

3. partition the full rational `u` interval into exactly 64 equal rational
   subintervals;
4. on the 63 subintervals separated from zero, evaluate `W`, `Bbar` and
   `-C'` by their direct quotient formulas, with no mean-value replacement;
5. on the first interval only, use explicit Lagrange Taylor formulas at
   `u=0`, with derivative remainders evaluated on that first interval;
6. repeat the complete certificate at 160 and 256 decimal digits;
7. build its own theorem verdict before reading the primary invariant-region
   artifact.

The primary artifact may be compared only after both precision records exist.

## 3. Direct quotient formulas

With `M`, `P`, `P0` and

```text
N=P-P0+2*pi*u*M,
```

evaluate away from zero

```text
W=(P-P0)/u,
Bbar=N/u^2,
W'=[u*P'-(P-P0)]/u^2,
-C'=-W'-4*pi*M'.
```

At the first interval, use

```text
W = P'(0)+P''(0)*u/2 + P'''(xi)*u^2/6,

Bbar = N''(0)/2+N'''(0)*u/6+N''''(xi)*u^2/24.
```

For `W'`, enclose

```text
P''(0)/2+P'''(0)*u/3
```

with absolute remainder at most

```text
(5/24)*sup|P''''|*u^2.
```

All endpoints and coefficients must be Arb balls with explicit lower and
upper endpoints in the artifact.

## 4. Frozen gates

At each precision independently require:

- positive lower bounds for `M`, `Bbar`, `C` and `-C'` on all 64 leaves;
- `U<1` throughout the half-strip;
- `y_plus/z>0`;
- `partial_z Y>0`;
- the exact normalized same-`x` gap is positive;
- the global `R'(q)` root argument excludes every other physical root;
- the accepted branch-B seed belongs to the half-strip;
- direct redifferentiation of the action reproduces its fifth successor.

The minimum lower bounds at 160 and 256 digits must overlap and have the same
sign.  No fitted agreement tolerance is used; comparison is by Arb overlap and
the already frozen serialized precision of the fifth artifact.

## 5. Hostile controls

- The first direct quotient interval must be unresolved if its zero endpoint
  is naively divided rather than treated by the frozen Taylor formula.
- Reversing the outgoing momentum sign must fail the direct shared-slice
  action equation.
- Omitting the boost term must change `Bbar(0)`.
- At `x=60`, the compactified positive-height gate must fail exactly.
- Diagnostic threshold perturbations `x>=124` and `x>=126` may be reported,
  but they cannot alter acceptance and cannot be used as fallback domains.

## 6. Outcomes

### `INVARIANT_HALF_STRIP_ADVERSARIALLY_CORROBORATED`

Both independent precision records pass, their rigorous bounds overlap, the
direct action agrees, and the later comparison matches the primary artifact.
Only then accept by induction a unique physical successor at every later
finite step for the branch-B history.

### `INVARIANT_HALF_STRIP_ADVERSARIAL_REFUTED`

A strict wrong sign or a direct-action disagreement is certified.  Preserve
the witness and retract the primary theorem.

### `INVARIANT_HALF_STRIP_ADVERSARIAL_OPEN`

Use for any unresolved leaf, non-overlap, axis problem, root-partition gap or
comparison ambiguity.  Do not shrink the domain.

Even the positive outcome remains homogeneous, representative-seed scoped and
dimensionless.  It does not make infinite extendibility a local law, classify
all original `v`, introduce inhomogeneous degrees of freedom or derive an
absolute tick.
