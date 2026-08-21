# Adversarial protocol: fourth-slab extendibility

Date: 2026-08-21.

Primary protocol commit: `54c4554`.

Primary verifier registration commit: `a58c46e`.

Accepted primary artifact commit: `7601c8f`.

Accepted primary artifact SHA-256:

```text
cf322cf0d60668d8f3f58e251425c9ad6bf43b112f22f9f3aebbc28f86212468
```

Status: frozen before reading the primary fourth-slab artifact.  Its reported
one physical root is known; this is independent replication, not blind
discovery.

## 1. Independent history reconstruction

At 110 and 180 decimal digits, reconstruct from the complete differentiated
action:

1. the first slab at `v=3/2`;
2. second branch B;
3. its unique physical third slab;
4. the outgoing state `(m3,pi3)` using `p_post/r3^2` only.

Do not import any fourth incoming value or root from the primary artifact.

## 2. Dual fourth-root proof

Do not use the primary function `E4` or its derivative.  Solve the constraint
first and classify

```text
h_C(q)=2[m3-mu(q)]/[q*mu(q)],
R4(q)=p(q)-pi3+4*pi[mu(q)-m3]/q,
R4'(q)=4*pi[m3-mu(q)]/q^2.
```

Enumerate all equal-`mu` stationary points, both infinite tails and both
one-sided limits at zero.  Treat `q=0` in the original constraint.  No finite
root box may support the count.

Require reduced equations for every algebraic root and direct full-action
equations for every root satisfying `h>0` and `1+h*q>0`.

## 3. Precision and hostile controls

- Counts, physical labels and roots must nest beyond 60 digits at 110 and
  180 digits.
- Wrong `p_post/r`, reversed post sign and reset mass must change the fourth
  incoming state at both precisions.
- Read the primary artifact only after the independent census; then compare
  all serialized roots and the outgoing state.

## 4. Outcome hierarchy

Use the primary outcome with suffix `_ADVERSARIALLY_CORROBORATED` only if the
dual proof gives the same complete real and physical count and every direct
gate passes.

Use `FOURTH_SLAB_EXTENDIBILITY_DISAGREEMENT` for a different count or
physical label.  Use `FOURTH_SLAB_EXTENDIBILITY_ADVERSARIAL_OPEN` for an
incomplete tail, exceptional level, precision or provenance gate.

Even a confirmed unique fourth slab remains **DERIVED COMPUTATIONAL,
four-slab scoped / STRUCTURAL**.  Indefinite continuation and any asymptotic
self-similarity remain **OPEN** until separately preregistered.
