# Preregistration: stationary roots versus the second canonical target

Date: 2026-08-16

Prior-art gate: `fcc4d7c`.

Status: frozen after the target-independent root multiset was committed at
`caaf1f1` and before either root is compared with the desired momentum.

## 1. Frozen inputs

Load exactly:

```text
reproducible/gravity_600cell_dust_stationary_root_enumeration.json
reproducible/gravity_600cell_dust_homothetic_canonical_lapse.json
reproducible/gravity_600cell_dust_two_slab_gluing.json
```

Require the first file to have outcome `STATIONARY_ROOTS_ENUMERATED`, two
roots, `target_parsed=false`, artifact SHA-256
`0ec5ba520ea25b39dd6cfd3c349d49fe480df2abee359854e1316b5af4d9fa2f`
and provenance commit `07083cc`.  Require the accepted first-tick artifact to
have SHA-256
`4b1c59c0518eec11b88b140cdecdf558d762c0d70b4826a758f67544e14ac5b9`,
outcome `HOMOTHETIC_CANONICAL_LAPSE_SELECTED`, root-result commit `b788258`
and a passing complete gate for both parities.  The artifact/result commit
containing that accepted evaluation is `46a7361`; it is external Git
provenance rather than a JSON field.

The total dust mass remains the fixed original mass.  No mass or scale is
readjusted in this comparison.

## 2. Target construction

For each parity, take the 30-component `post_momentum` of the accepted first
tick and map it into the next slab's old-boundary ordering with the independently
derived `old_to_final_orbit_map` from the gluing artifact:

```text
target[i] = first_post[old_to_final_orbit_map[i]].
```

Require each map to be a permutation of `0,...,29`; do not replace it by the
identity even if the stored permutation happens to be the identity.

## 3. Exhaustive frozen comparison

Compare **both** committed roots, in their committed order, without changing
`b`, `r`, `P` or any internal variable.  For root `j` and parity `p`, define

```text
residual[j,p,i] = root_pre[j,i]-target[p,i].
```

Record the complete residual vector, Euclidean norm, maximum absolute
component and component spread.  A root is an exact canonical hit only if,
for both parities,

```text
norm_2(residual) <= junction_bound,
```

where `junction_bound` is inherited verbatim from the accepted first-tick
artifact.  No looser visual, relative or fitted tolerance is allowed.

Report the hit fraction `hits/2`.  A small residual outside the bound remains
a failure; it may only seed a separately preregistered later correction.

## 4. Structural labels

Independently of target values:

- call a root `CONTRACTING` only if its upper log `b` is strictly below the
  committed lower log `a1`;
- call it `TIME_REVERSAL` only if `b=0` within the root-enumeration tolerance
  and the enumeration already labels it as such;
- do not choose a root by whichever has the smaller target residual.

## 5. Mechanical outcomes

Assign exactly one:

1. `SECOND_TICK_TARGET_CONTROL_FAILED` if input/provenance/map controls fail;
2. `STATIONARY_SECOND_TICK_HIT` if at least one unmodified root passes both
   parity bounds;
3. `STATIONARY_SECOND_TICK_NO_HIT` if controls pass and neither root does.

Outcome 2 would select an exact second tick without a target-dependent solve.
Outcome 3 is not a no-go for local canonical correction, because each stored
root is stationary only at the inherited first-tick lapse.  It merely proves
that holding that lapse exactly fixed does not glue.

Only the new targeted verifier will be run.  The full suite will not be run.
