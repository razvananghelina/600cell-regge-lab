# Adversarial protocol: branch-selector audit

Date: 2026-08-21.

Primary selector protocol commit: `e3c8927`.

Primary selector verifier registration commit: `70129cf`.

Primary selector artifact commit: `b0c2fe9`.

Status: frozen before reading the primary selector artifact.  The two
composition roots were already known, so this is an independent-mechanism
replication, not blind discovery.

## 1. Independent inputs

Use the previously accepted adversarial composition artifact

```text
reproducible/gravity_600cell_finite_height_composition_adversarial.json
```

whose roots were obtained from direct solves of the complete differentiated
action at 80, 120 and 180 decimal digits without the primary scalar
elimination equation.

Do not read

```text
reproducible/gravity_600cell_finite_height_selector_audit.json
```

until the independent result is complete.

## 2. Independent causal reconstruction

For each 180-digit direct-action root reconstruct a representative
same-vertex strut in the certified Minkowski embedding:

```text
Delta_R = phi*h*q,
T       = sqrt(h^2+Delta_R^2),
s^2     = -T^2+Delta_R^2.
```

Require, without using the primary closed `beta` formula,

```text
s^2=-h^2<0,
T>0,
|Delta_R/T|<1,
1+h*q>0.
```

Also evaluate the direct cellular angle arguments and require both roots to
remain on the same real analytic branch:

```text
1/3 <= (q^2+2)/(2(q^2+3)) < 1/2,
q/sqrt(8(q^2+3)) real and finite.
```

No cubical central-height inequality may be imported.

## 3. Independent local regularity

Use the direct full-action Jacobians already reconstructed at all three
precisions.  Require each root's Jacobian to be nonzero at every precision
and stable beyond 60 digits.  Only after that compare the 180-digit direct
Jacobian with the coefficient identity printed by the primary selector
artifact.

This confirms local regularity by a different calculation.  It does not turn
local inversion into global injectivity.

## 4. Hostile controls

- Treat `T=h` as if `h` were central coordinate height.  This must make the
  high-`q` root appear superluminal and demonstrates the convention trap.
- Replace the timelike Minkowski sign by a Euclidean sign.  The strut
  certificate must fail.
- Confirm that the two direct roots solve the same normalized incoming
  `(m1,pi1)` data and remain distinct.

## 5. Outcome hierarchy

Use

```text
STANDARD_CANONICAL_SELECTORS_DO_NOT_RESOLVE_BRANCH_
ADVERSARIALLY_CORROBORATED
```

only if both direct roots are future oriented, causal in the independent
embedding, locally regular and on the same real action branch, while both
hostile conventions fail.  Label the conclusion **DERIVED NEGATIVE,
selector-scoped / adversarially corroborated**.

Use `SELECTOR_AUDIT_DISAGREEMENT` if any physical classification or local
regularity result disagrees with the primary artifact.  Such disagreement
stops interpretation.

Use `SELECTOR_AUDIT_ADVERSARIAL_OPEN` for any incomplete precision or
provenance gate.

Only the targeted adversarial verifier may be run.
