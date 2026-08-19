# Protocol: exact certification of the 119 static vertex-gradient modes

Date: 2026-08-19

This protocol is committed after the target-blind modular census and before
constructing any candidate static null vector.

## Frozen provenance

| input | SHA-256 |
|---|---|
| vertex-gradient prior-art gate | `9541383a435e069be13ed9c2175674036a9cef5cd4e17ee455f524bd1c1c6a7d` |
| global closure protocol | `ed79c6a15ade377ae09854b3cad3028eb1c0f43cf8e85789d46993fe25ed1b49` |
| global closure verifier | `ec44be8e4d82634e30944739d10d3f80fbb9f6fee0883ec1f612690c38d90ab6` |
| target-blind modular artifact | `61cebd1cd67fcdc56de088855b1fc7b805d0f70f9f9b3029d4a61209d7a53944` |
| local variable-face theorem | `2db55cb87ec1c01d537cdbc11010bc9ea740762c598108e4c2de0f3acca72cc8` |

The artifact must retain `12/12`, expanding ranks 3600 and static ranks 3481
at both primes.  This protocol may not change the static closure operator.

## Exact local gradient

For the canonical spatial tetrahedron vertices `p_0,...,p_3`, let

```text
D=[p_1-p_0, p_2-p_0, p_3-p_0]
```

using their first three coordinates.  For local scalar values
`(u_0,u_1,u_2,u_3)`, define

```text
grad_T u = (D^T)^(-1) (u_1-u_0,u_2-u_0,u_3-u_0)^T.
```

This gives an exact rational `3 x 4` matrix `H`.  Require

```text
H (1,1,1,1)^T = 0,
rank(H)=3,
(grad_T u).(p_i-p_0)=u_i-u_0.
```

No vector from the observed global modular kernel may be used to define `H`.

## Global map

Rebuild the exact 600-cell incidence and the same static reduced face blocks
as the frozen global verifier.  For every sorted tetrahedron `T`, map its four
global vertex values through `H` and place the result in the local static
translation coordinates:

```text
z_T(u)=(0,0,0, grad_T u).
```

Assemble

```text
G : Q^120 -> Q^(6*600),
shape 3600 x 120.
```

Require exactly:

1. `G 1=0`;
2. rank 119 modulo both `1000003` and `1000033`;
3. every one of all 1200 rational face blocks satisfies its local identity

   ```text
   C_f [G_source;G_target]=0;
   ```

4. consequently `C_static G=0` without assembling a fitted nullspace.

The explicit constant kernel gives `rank_Q(G)<=119`; a modular rank 119 gives
`rank_Q(G)>=119`.  Together they certify exact rational rank 119.

## Exhaustion proof

Independently reproduce the static modular ranks

```text
rank(C_static mod 1000003)=3481,
rank(C_static mod 1000033)=3481.
```

Hence `rank_Q(C_static)>=3481` and

```text
dim ker_Q(C_static) <= 3600-3481=119.
```

The rational 119-dimensional image of `G` lies in that kernel, so require the
verifier to assign

```text
ker_Q(C_static)=im_Q(G),
dim ker_Q(C_static)=119.
```

This conclusion uses rank inequalities and exact inclusions, not an assumed
equality of modular and rational nullspaces.

## Falsification controls

1. On the lexicographically first shared face, deliberately swap two target
   shared-vertex scalar columns while leaving the source columns unchanged.
   The resulting discontinuous local map must not be annihilated by that face
   block.
2. Apply the odd canonical relabelling `(0 1)`, rebuild both `C` and `G`, and
   require the same exact inclusion and ranks.
3. Repeat the face identities under `eta -> -eta`.

## Outcome hierarchy

1. `STATIC_VERTEX_GRADIENT_CONTROL_FAILED` if provenance, incidence, local
   gradient, modular rank, discontinuity, relabelling or metric-sign control
   fails.
2. `STATIC_KERNEL_EXACTLY_VERTEX_GRADIENTS` if all controls pass, `G` has
   exact rational rank 119, `C G=0`, and the exhaustion inequalities close.
3. `STATIC_VERTEX_GRADIENTS_PROPER_SUBSPACE` if a rational 119-dimensional
   gradient image is certified but an explicit additional rational kernel
   vector exists.
4. `STATIC_VERTEX_GRADIENT_HYPOTHESIS_REFUTED` if any correctly assembled
   face fails to annihilate `G` while construction controls pass.
5. `STATIC_VERTEX_GRADIENT_OPEN` otherwise.

## Interpretation firewall

The accepted result, if positive, is a finite-element/de Rham identification.
It does not name the scalar as lapse, gauge, time or matter.  Such a name
requires an action or canonical constraint with the same nullspace.

No Hessian, dynamics, wave speed or full suite is authorized here.
