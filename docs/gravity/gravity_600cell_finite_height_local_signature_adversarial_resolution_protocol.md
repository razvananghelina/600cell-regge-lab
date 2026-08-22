# Resolution protocol: monotone-factor direct bisection

Date: 2026-08-22.

First adversarial `OPEN` commit: `3ee8dcb`.

Status: frozen before the resolution verifier exists or any gravity bracket is
evaluated by the resolution method.

## 1. Disagreement to resolve

The first direct-bisection route obtained strict opposite endpoint signs on
the initial stationary bracket `[5,6]` and root bracket `[9,10]`, but raw Arb
evaluation of `p'` and `p-pi` over a whole unit interval contained zero due to
dependency. It therefore proved existence but not uniqueness and correctly
returned `LOCAL_SIGNATURE_ADVERSARIAL_OPEN`.

The resolution must retain every frozen integer bracket, use no interval
Newton, discovery seed or primary root ball, and must not subdivide until a
desired sign appears.

## 2. Independent accepted monotonicity theorem

Use the already accepted and adversarially corroborated one-slab
classification identity

```text
p'(q)=-720*K(q^2)/[
  (q^2+4)^(3/2)*(q^2+3)*sqrt(3*q^2+8)
],
```

where `K(x)` has exactly one positive zero `x_star` and

```text
v_star=sqrt(x_star).
```

Consequently, for real `q`, `p` is strictly:

```text
increasing on (-infinity,-v_star),
decreasing on (-v_star,v_star),
increasing on (v_star,+infinity).
```

Recheck the derivative factorization symbolically and verify the hash and
outcome of the accepted classification artifact. Do not numerically infer
monotonicity from the local tree.

## 3. Frozen brackets and endpoint proof

Use exactly the stationary and root brackets in
`gravity_600cell_finite_height_local_signature_adversarial_protocol.md`.

For a stationary bracket of `p(q)-pi`:

1. prove that the bracket lies in one of the three strict monotonicity regions;
2. certify strict opposite Arb signs at the rational endpoints;
3. conclude exactly one stationary point by continuity and monotonicity.

For a root bracket of `E`:

1. certify strict opposite Arb signs of `E` at the rational endpoints;
2. partition the bracket conceptually only at `-v_star` or `v_star` if either
   turning point lies inside it;
3. evaluate `p(q)-pi` at both rational endpoints and at every contained
   turning point;
4. require all these values to have the same strict sign;
5. use monotonicity of `p` on each resulting segment to conclude that
   `E_q=p(q)-pi` has that sign everywhere on the complete root bracket.

This proves existence and uniqueness without interval evaluation of a raw
derivative over a wide ball.

## 4. Direct bisection and recursion

After uniqueness is proved, bisect each exact rational sign-changing bracket
for at most 420 dyadic steps. Retain the sign-changing half. If a midpoint Arb
ball first contains zero after at least step 240, retain the current bracket;
earlier ambiguity is `OPEN`.

Propagate the entire final root interval through the exact state update. Repeat
the stationary-count, tail, origin, complete-root, strict-gate and terminal
checks from the first adversarial protocol. Require the endpoint identity
width below `1e-110`, exactly one strict `D` entry, and rejection of `m*q>126`.

## 5. Controls and comparison ordering

Before gravity evaluation:

- `[1,2]` must certify the unique root of `q^2-2` using direct monotonicity;
- `[2,3]` must fail the endpoint sign gate.

After the complete resolution tree exists in memory:

- compare it with the preserved first-adversarial endpoint signs;
- compare ordered root enclosures, counts and strict gate signs with the
  primary artifact;
- require all comparisons to agree without repairing the resolution result.

## 6. Outcomes

### `LOCAL_SIGNATURE_ADVERSARIAL_DISAGREEMENT_RESOLVED`

The exact monotonicity factor resolves every wide-bracket dependency; all
fixed brackets, root counts, recursive states, terminal gates and controls
pass; the primary and resolution trees agree.

This counts as the mechanically different adversarial corroboration required
by Rule 4.

### `LOCAL_SIGNATURE_ADVERSARIAL_DISAGREEMENT_OPEN`

Use for any factorization/provenance failure, bracket crossing not handled by
the frozen turning points, inconsistent sign, early midpoint ambiguity,
terminal disagreement or failed control. Preserve the artifact and do not
accept the primary theorem.

## 7. Scope

Resolution would prove only existence of an unspecified local neighbourhood
on the special incoming curve. It does not compute its radius, classify the
global basin, derive `v=3/2`, make extendibility local or add nonhomogeneous
physics.

