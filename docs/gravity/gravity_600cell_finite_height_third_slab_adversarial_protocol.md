# Adversarial protocol: third-slab future extendibility

Date: 2026-08-21.

Primary protocol commit: `0f31fe8`.

First primary failure preserved: `564dca9`.

Bisection correction frozen: `9c40419`.

Corrected primary implementation: `cceb0cc`.

Accepted primary artifact commit: `a0abf0a`.

Accepted primary artifact SHA-256:

```text
6b0e92d031aa891fdc3e1b2045c35bd135a955bb1374c92f015dcd5727d3d8fc
```

Status: frozen before reading the primary third-slab artifact.  The reported
physical counts `A:0, B:1` are known, so this is a mechanically independent
replication, not blind discovery.

## 1. Independent reconstruction

Reconstruct the complete cellular action, first slab and both second slabs
directly at 100 and 160 decimal digits.  Use the already frozen rational
seeds but do not import either outgoing state or any third root from the
primary artifact.

Compute each second branch's outgoing state only as

```text
m2=m1/r2,
pi2=p_post,2/r2^2
```

from the redifferentiated action.

## 2. Dual elimination

Do not use the primary elimination function

```text
E(q)=4*pi[mu(q)-m]+q[p(q)-pi]
```

or its decisive derivative `E'=p-pi` to count roots.

For `q!=0`, solve the constraint first:

```text
h_C(q)=2[m-mu(q)]/[q*mu(q)].
```

Insert this into the momentum residual and classify zeros of

```text
R(q)=p(q)-pi+4*pi[mu(q)-m]/q.
```

Using the independently established state-derivative identity, derive

```text
R'(q)=4*pi[m-mu(q)]/q^2.
```

Thus this proof is controlled by equal-`mu` points, whereas the primary proof
is controlled by equal-`p` points.

Treat `q=0` directly in the original two equations; it is not covered by the
division above.

## 3. Complete real-line count

For each branch, enumerate every solution of `mu(q)=m2` from the exact shape
of `mu`:

- `mu` is even;
- it rises from `mu(0)` to its unique maximum at `v_star`;
- it decreases strictly to zero afterwards.

Use these points to partition

```text
(-infinity,0) union (0,+infinity)
```

into monotone intervals of `R`.  Certify the signs at all finite stationary
points and the four one-sided limits:

```text
q -> -infinity,
q -> 0-,
q -> 0+,
q -> +infinity.
```

No finite root box may be used as evidence.

For every counted zero reconstruct `h_C`, require the original direct
full-action constraint and momentum residuals below `1e-90`, and classify it
with the unchanged physical conditions `h>0` and `1+h*q>0`.

## 4. Precision and hostile controls

- Repeat the complete construction at 100 and 160 digits and require all
  roots and signs to nest beyond 60 digits.
- The wrong `p_post/r` normalization, reversed post sign and reset mass must
  change both outgoing states.
- Explicitly verify that `q=0` is not a simultaneous root for either frozen
  outgoing state.

Only after all independent classifications pass may the primary artifact be
read.  Compare branch counts and root values at their serialized precision.

## 5. Outcome hierarchy

Use

```text
ONE_SECOND_BRANCH_EXTENDS_UNIQUELY_ADVERSARIALLY_CORROBORATED
```

only if the equal-`mu` proof gives `A:0, B:1`, the B root passes the complete
action at both precisions and every hostile control passes.

Use `THIRD_SLAB_EXTENDIBILITY_DISAGREEMENT` for any different all-real or
physical count.  Such disagreement stops interpretation.

Use `THIRD_SLAB_EXTENDIBILITY_ADVERSARIAL_OPEN` for incomplete tails,
exceptional strata, precision or provenance.

Even a confirmed result is only **DERIVED COMPUTATIONAL, three-slab scoped /
STRUCTURAL selection candidate**.  It does not establish that future
extendibility is a fundamental physical axiom or that the B branch extends
indefinitely.
