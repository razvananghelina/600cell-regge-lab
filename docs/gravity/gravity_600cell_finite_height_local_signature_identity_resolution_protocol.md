# Protocol: exact resolution of the adversarial endpoint identity

Date: 2026-08-22.

Monotone-factor `OPEN` artifact commit: `91bb6db`.

Status: frozen before the identity-resolution verifier exists.

## 1. Exact disagreement

The monotone-factor adversarial route certified all roots, states, strict gates
and terminals but returned `8/10 OPEN` because the fifth interval evaluation of

```text
r-(1+h*q)
```

had width `5.63e-106`, above the preregistered auxiliary threshold `1e-110`.
All five identity balls contained zero. No root, sign or terminal disagreed
with the primary certificate.

The numerical width threshold must not be changed after the result.

## 2. Exact algebraic resolution

Treat `mu`, `p`, `m`, `pi`, and `q` as independent real symbols with
`mu!=0`. Define

```text
E=4*pi*(mu-m)+q*(p-pi),
h=(p-pi)/(2*pi*mu),
r=2*m/mu-1.
```

Verify symbolically

```text
r-(1+h*q)=-E/(2*pi*mu).
```

The monotone-factor artifact already certifies that each physical bracket
contains exactly one root of `E`. Verify directly with Arb that `mu(q)>0` on
every physical root interval. The displayed identity then proves
`r=1+h*q` at the actual root, independently of the width of an interval
evaluation away from that root.

## 3. Frozen artifact audit

Verify hashes and outcomes of:

- the primary `10/10` artifact;
- the first direct-wide-interval `4/11 OPEN` artifact;
- the monotone-factor `8/10 OPEN` artifact.

From the last artifact, independently require:

- exactly five state paths with counts `(2,1)`, `(3,2)`, `(2,0)`, `(3,1)`,
  `(3,1)`;
- every non-diagonal stationary and root bisection marked certified and
  monotonicity-certified;
- every primary-comparison row passed;
- exactly five physical transitions and one strict `D` entry;
- every stored numerical endpoint-identity ball contains zero;
- the preserved first-`OPEN` endpoint comparison passed;
- the hostile `m*q>126` claim was rejected.

Do not rerun or refit any root bracket in this resolver.

## 4. Negative control

On every physical root interval evaluate the deliberately false relation

```text
r=1+2*h*q.
```

Require `r-(1+2*h*q)` to exclude zero on every edge. If any false-identity
ball contains zero, the resolver is `OPEN`.

## 5. Outcomes

### `LOCAL_SIGNATURE_ENDPOINT_IDENTITY_EXACTLY_RESOLVED`

The exact symbolic factorization, positive `mu`, complete artifact audit and
false-identity control all pass. Combined with the monotone-factor root/tree
certificate, this resolves the sole remaining adversarial failure and permits
consolidation of the local theorem.

### `LOCAL_SIGNATURE_ENDPOINT_IDENTITY_OPEN`

Use for any symbolic, positivity, artifact, root-existence, tree or negative
control failure. Preserve the artifact and leave the primary theorem
unconsolidated.

## 6. Scope

This resolver removes an interval-width artifact. It adds no explicit
neighbourhood radius, global basin, local physical selection rule or
nonhomogeneous degree of freedom.

