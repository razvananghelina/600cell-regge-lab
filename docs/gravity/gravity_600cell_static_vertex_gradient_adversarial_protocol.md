# Adversarial protocol: static rotation/translation block decomposition

Date: 2026-08-19

This protocol is committed after the primary rational exhaustion theorem and
before constructing an independent spatial operator.  It must not use the
four-dimensional lateral-face reflection, Poincare stabilizer, raw `12 x 13`
face block or its algebraic elimination.

## Frozen inputs

| input | SHA-256 |
|---|---|
| prior-art gate | `9541383a435e069be13ed9c2175674036a9cef5cd4e17ee455f524bd1c1c6a7d` |
| primary protocol | `c689202b94abfe9436bdc8f5db6f79fa00d6a7bc743d79bc3ae34d391890c17b` |
| primary verifier | `6974ef8c85ee62b32daa2277ba221ee1cfa96f1c7cc7a92a8ec91fad576b124f` |
| primary artifact | `ce018db5c66c78e89e4ca32360385955ea520b9ac8e42955b110d190432239c0` |
| target-blind modular artifact | `61cebd1cd67fcdc56de088855b1fc7b805d0f70f9f9b3029d4a61209d7a53944` |

The primary artifact must retain `10/10` and
`STATIC_KERNEL_EXACTLY_VERTEX_GRADIENTS`.

## Independent static operator

At `lambda=1`, each local six-vector decomposes directly as

```text
(omega_T, v_T) in R^3 rotations + R^3 translations.
```

Rebuild only the exact **spatial** 600-cell face transitions.  For adjacent
tetrahedra `T <- T'`, construct the unique Euclidean isometry mapping their
three shared canonical face vertices and reflecting the opposite apex.  Let
its linear part be `R_f`.

### Rotation block

Compute the `3 x 3` adjoint `Ad_f` by conjugating the three explicit
antisymmetric rotation matrices with `R_f`; do not use an axial-vector sign
formula.  Impose

```text
omega_T - Ad_f omega_T' = 0
```

for all 1200 faces.  The resulting exact matrix has shape `3600 x 1800`.
The disclosed prediction is full column rank 1800: spatial curvature permits
no nonzero globally parallel rotation field.

As a flat negative control, replace every `Ad_f` by the identity.  Connected
graph incidence must then have rank 1797 and the three constant rotations as
kernel.

### Translation block

Let `t_1,t_2` span the source shared triangle.  Impose only tangential
continuity

```text
t_a . (v_T - R_f v_T') = 0,  a=1,2.
```

The normal difference is the derived hinge mode and is deliberately not
constrained.  The exact matrix has shape `2400 x 1800`.  The disclosed
prediction is

```text
rank 1681, nullity 119.
```

On a breadth-first spanning tree require rotation rank 1797 and translation
rank 1198, hence combined rank/nullity `2995/605`.

## Independent gradient certificate

Construct the continuous P1 gradient from canonical spatial barycentric
coordinates and one scalar per global vertex.  Require exact rational
annihilation by the translation block, modular rank 119, and failure after a
deliberate shared-vertex mismatch on the first face.

The combined block-diagonal static operator must therefore have

```text
rank = 1800+1681 = 3481,
nullity = 119,
kernel = zero rotations direct-sum P1 gradients.
```

All ranks are evaluated exactly modulo `1000003` and `1000033`.  The explicit
rational gradient image plus the modular lower bounds closes the rational
rank argument as in the primary proof.

## Convention attack

Apply the odd canonical relabelling `(0 1)` and rebuild the spatial
transitions, rotation adjoints, tangential rows and gradients.  Every rank and
inclusion must persist.

## Outcome hierarchy

1. `ADVERSARIAL_STATIC_DECOMPOSITION_CONTROL_FAILED` if provenance,
   incidence, spatial isometry, flat control, tree control, discontinuity or
   relabelling fails.
2. `ADVERSARIAL_STATIC_KERNEL_IS_P1_GRADIENTS` if the rotation sector is full
   rank, the translation sector is exhausted by the 119 P1 gradients and the
   combined rank is 3481 at both primes.
3. `ADVERSARIAL_STATIC_EXTRA_ROTATION_OR_TRANSLATION` if an explicit rational
   mode outside the gradient image is found.
4. `ADVERSARIAL_STATIC_DECOMPOSITION_OPEN` otherwise.

## Interpretation firewall

Corroboration identifies a discrete-gradient kinematic kernel.  It still does
not prove that the scalar is canonical gauge or physical time.  That requires
the action/Hessian or a constraint algebra.

No dynamics or full suite is authorized here.
