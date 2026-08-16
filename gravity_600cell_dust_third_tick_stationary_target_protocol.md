# Frozen specification after diagnostic disclosure: third-tick target comparison

Date: 2026-08-16

Prior-art gate: `7b9a676`.

Target-independent root multiset committed first: `3401137`.

Status: **not a clean preregistration of the arithmetic comparison**.  After
the root multiset was safely committed at `3401137`, a local read-only shell
diagnostic parsed the second-tick target and printed the two scalar residuals
before this file was committed.  No root, grid, label, bound or later solver
state was changed.  The deterministic complete comparison below is frozen
before its registered verifier is run, but its result direction is already
known and must be labelled accordingly.

This process breach does **not** invalidate the target-independent root
enumeration: its commit and byte-level target firewall preceded the diagnostic.
It also does not authorize tuning a later correction.  Any target-dependent
Newton solve must receive its own clean pre-evaluation protocol and commit.

## Frozen inputs

Require exact SHA-256 values

```text
third stationary roots:
02d4589a7df0851c67a31fc0a41c5ef8851a82c758214c1c5e8729afddfe479f,

accepted second tick:
936984bc84a714140ce16917ee559b346b3c0d4a5ba92d8fb723398a120f8e70,

gluing map:
a5a22d219b71e49c154c1ef80ed9da93b1aef0b93cd2d6ed22f041b71f62db77.
```

Require `target_parsed=false`, outcome
`THIRD_TICK_STATIONARY_ROOTS_ENUMERATED`, `N=2` and 5/5 in the first artifact;
`SECOND_HOMOTHETIC_TICK_ACCEPTED` and 6/6 in the second.

## Target construction and exhaustive comparison

For each parity, map all 30 accepted second-tick post-momenta into the next
old-boundary ordering:

```text
target[i] = second_post[old_to_final_orbit_map[i]].
```

Require the stored map to be a permutation; apply it even if it is the
identity.  Compare **both** committed roots without changing `C`, `R2`, `P` or
any edge.  For each root and parity record all 30 values

```text
residual[i] = root_pre[i]-target[i],
```

their Euclidean norm, maximum absolute component and component spread.

A root is a hit only if both parities satisfy

```text
norm_2(residual) <= inherited second-tick junction_bound.
```

Report `hits/2`.  No relative or visually small mismatch counts as a hit.

## Structural labels

Assign labels before considering target proximity:

- `CONTRACTING` iff `C<B2`;
- `TIME_REVERSAL` iff `C=A1` within `1e-25` and the root is the committed node
  cluster.

Do not select or discard a candidate based on its residual.

## Mechanical outcomes

Assign exactly one:

1. `THIRD_TICK_TARGET_CONTROL_FAILED`;
2. `STATIONARY_THIRD_TICK_HIT` if at least one unmodified root passes both
   parity bounds;
3. `STATIONARY_THIRD_TICK_NO_HIT` otherwise.

A no-hit result only says that the inherited lapse `R2` does not glue exactly.
It may seed one separately preregistered local `(C,R)` correction on the
geometrically contracting root, but no alternate root or solver endpoint.

Only the new targeted verifier will be run.  The full suite will not be run.
