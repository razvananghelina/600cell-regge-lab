# Protocol: exact reconciliation of the two 119-dimensional carriers

Date: 2026-08-19

Prior-art/framing commit: `fe9fc18`.

This target-disclosed protocol is frozen before constructing the coordinate
intertwiner.  It tests a literal equality, not a desired spectrum.

## Frozen provenance

| input | SHA-256 |
|---|---|
| reconciliation prior-art gate | `46134f34b396ecc7dd844c85a83a561da29c8db0896d33afea8865fecc00718e` |
| prism-shift gluing protocol | `f76a2216f20ea00aa66c31891d81c1f72a4ee86fe519458b6527118a4fca2251` |
| prism-shift gluing result | `2cefd7a24a6ac132da34cbe450210446b1b3bd0b39a4920dcc36176fc4a68e1a` |
| prism-shift verifier | `0faa50e20f3efd89b8828426d83aba5d92401bc59e72a1091653761c4ab23519` |
| prism-shift artifact | `1ab6654ae57c83a49dd4f427154b891c0b8ae613631773ab6733a1227b9999fa` |
| static-gradient protocol | `c689202b94abfe9436bdc8f5db6f79fa00d6a7bc743d79bc3ae34d391890c17b` |
| static-gradient adversarial protocol | `53c094c120411e0b63cafd9b8b5e4b60880c7f1339eabec25eed7f75ccdf2805` |
| static-gradient verifier | `6974ef8c85ee62b32daa2277ba221ee1cfa96f1c7cc7a92a8ec91fad576b124f` |
| static-gradient artifact | `ce018db5c66c78e89e4ca32360385955ea520b9ac8e42955b110d190432239c0` |
| adversarial static-gradient verifier | `ac5774353d516ee33355e7850f1f32fc6a678d28d0979bb59141607405690c5f` |
| adversarial static-gradient artifact | `682ede43b77ada62ce8f7badb3c1672468ad14636e4cf7581af3e0b1a1a92632` |

All three artifacts must retain their passing outcomes and dimensions before
the reconciliation is interpreted.

## Carrier and coordinate maps

Rebuild the 600-cell exactly from the golden-field vertex set and clique
incidence.  Require

```text
f=(120,720,1200,600),
two tetrahedra per face.
```

For every sorted tetrahedron assign the canonical local points

```text
p0=( 1, 1, 1), p1=( 1,-1,-1),
p2=(-1, 1,-1), p3=(-1,-1, 1).
```

Put

```text
D=[p1-p0,p2-p0,p3-p0],
Q_T=D^T.
```

The direct sum `Q=direct-sum_T Q_T` maps a local Cartesian translation to
its three covector evaluations.  Require exactly

```text
det(Q_T)=-16,
rank(Q)=1800.
```

## Potential embeddings

Construct independently:

```text
B : Q^120 -> Q^1800,
(B phi)_(T,i)=phi(v_i)-phi(v_0),
```

which is the old prism-shift coordinate map, and

```text
G : Q^120 -> Q^1800,
(G phi)_T=(D^T)^(-1)(B phi)_T,
```

which is the Cartesian translation part of the new static Poincare map.
Require the entrywise rational identity

```text
Q G=B.
```

Require `rank(B)=rank(G)=119` over both primes `1000003,1000033`, with the
constant vector as their exact common kernel.

## Face operators

For every sorted shared face `(r,u,v)`, construct the old two prism-shift
rows by equating the incident covectors on oriented edges `(r,u)` and
`(r,v)`.

Independently construct the reflected Euclidean transition `R_f` from the
target local tetrahedron to the source local tetrahedron.  Construct the new
two Cartesian rows

```text
t^T (s_source-R_f s_target)=0,
t in {p_u-p_r,p_v-p_r} in the source frame.
```

With `C_old` and `C_new` ordered by the same faces and oriented edge pair,
require entrywise over the rationals

```text
C_new=C_old Q.
```

This stronger identity implies equal row spaces and, because `Q` is
invertible,

```text
Q ker(C_new)=ker(C_old).
```

Reproduce ranks `1681` and nullities `119` for both matrices over both
disclosed primes.

## Falsification controls

1. Inspect only the derived local spatial transports and select the
   lexicographically first face on which at least one of the two shared
   tangents has non-identity target transport.  On that face replace the
   target transition by the identity while leaving its distinct target frame
   unchanged.  The local row identity must fail.  The first execution showed
   that the lexicographically first face itself has identity tangential
   transport and is therefore an inert control; that `10/11` artifact is
   preserved in commit `a442429`.
2. Flip one Cartesian axis in exactly one cell of `Q` without modifying `G`
   or the face transitions.  Both `QG=B` and at least one incident face
   identity must fail.
3. Apply the odd canonical relabelling `(0 1)` globally, rebuild every map and
   require the valid intertwining identities and ranks to persist.

## Outcome hierarchy

1. `RECONCILIATION_CONTROL_FAILED` if provenance, incidence, invertibility,
   rank or negative controls fail.
2. `STATIC_GRADIENT_IS_PRISM_SHIFT_EXACTLY` if all exact potential and face
   intertwining identities pass.
3. `STATIC_GRADIENT_ONLY_ABSTRACTLY_ISOMORPHIC` if ranks agree but either
   literal intertwiner fails while controls pass.
4. `STATIC_GRADIENT_PRISM_SHIFT_RECONCILIATION_REFUTED` if a frozen upstream
   theorem or the common-potential relation fails.

## Interpretation firewall

A positive result proves that the two project labels denote one carrier in
different local coordinates.  All previously derived facts transfer:

- its equal-scale action Hessian is a positive multiple of the 600-cell graph
  Laplacian;
- unequal homogeneous scale removes it as an independent common-strut mode;
- the complete canonical equations eliminate all 119 relative directions.

This is **DERIVED RECONCILIATION**, not a second physical discovery.  No full
suite, time, wave speed, graviton, mass, Planck scale or `G` calculation is
authorized.
