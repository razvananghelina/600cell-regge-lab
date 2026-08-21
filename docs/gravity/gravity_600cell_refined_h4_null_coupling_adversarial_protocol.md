# Adversarial protocol: reconstruct null coupling from spatial hinge curvature

Date: 2026-08-21

Primary-result commit: `fd0e5f3`.

This protocol is frozen before writing or running the adversarial verifier.
It must not import or execute the primary null-coupling verifier, the full
Hessian verifier or any Lorentzian action evaluator.

## 1. Frozen inputs

Require:

```text
reproducible/gravity_600cell_refined_h4_null_coupling.json
  6b6fbd95b07f365b3fcac332fa3546021e8d756a510af0184bc974e52d5efa79
docs/gravity/gravity_600cell_refined_h4_null_coupling_primary_result.md
  5dfefd1f0b2fdcae02cede4a9e7d069a5e7e3d0c29b6f1d324368a4cdbe8803a
reproducible/gravity_600cell_refined_local_curvature_mass_adversarial.json
  c59890d12bf929c4677dffed1b932ad8c05ab0ac00980be15ba780e62744c28e
reproducible/gravity_600cell_refined_boundary_cotangent_adversarial.json
  19c888a43bdba9d57166d6e3595c6d5b51dd019ebf616efdbf1189e25078f808
docs/gravity/gravity_600cell_refined_local_curvature_mass_result.md
  ef6e29fc1e4c89d893a40ee2b5efb3ab6c833e0d73ec232bdbd41033bc4f0f94
docs/gravity/gravity_600cell_refined_boundary_cotangent_result.md
  391a317b9f8823a5479f450dde43a43177e210a2d81192aedc938e90fc8006d1
```

Require all four accepted upstream outcomes and the primary `16/16`
null-coupling outcome.

## 2. Independent product-hinge derivation

Use SymPy and Heron's squared-area polynomial for one lower-boundary product
hinge with squared data

```text
(x,y,z)=(l^2,-tau^2,l^2-tau^2).
```

Derive, rather than enter as constants,

```text
A=i*l*tau/2,
dA/dlog(l^2)=i*l*tau/4.
```

With the repository's `-i A epsilon` convention this gives the boundary
action-gradient contribution

```text
g_boundary,e=tau*l_e*epsilon_e/4.
```

Differentiate this expression with respect to `u=log(tau^2)` and require

```text
d g_boundary,e/du=tau*l_e*epsilon_e/8.            (1)
```

No primary Hessian quantity may enter this derivation.

## 3. Actual-incidence reconstruction

Read the six actual-incidence curvature totals

```text
C_rs=sum_(actual edges of pair rs) l_e epsilon_e
```

only from the already accepted adversarial boundary artifact.  Construct in
pair order `(01,02,03,12,13,23)`

```text
c_hinge=(tau0*C_rs/8 for six pairs),
c_adv=(c_hinge,c_hinge).                          (2)
```

Write (2) before reading the primary compatibility rows.  Then compare it
componentwise with all 24 primary rows; require maximum error `<1e-68`.

As an independent vertical-null consistency check, use the four actual
rank-curvature totals `K_r` and verify algebraically

```text
m_r=K_r/(8*pi),
K_r/2-4*pi*m_r=0,                                 (3)
```

so the product vertical residual
`tau*(K_r/2-4*pi*m_r)` and its log-lapse derivative vanish for arbitrary
positive `tau`.  The six cross components of the null image are not rebuilt
by this verifier; their evidence remains the primary 72-point finite-family
control.  This limitation must be explicit in the artifact.

## 4. Controls

- Reconstruct the accepted adversarial boundary `pre/post` vectors from
  `C_rs` and require exact agreement inside `1e-68`.
- Replacing `1/8` in (2) by `1/4`, reversing the sign, dropping the largest
  `C_rs`, and swapping two unequal pair totals must each fail by `>1e-6`.
- Require the repeated old/new row in (2) to have explicit rank one and to be
  invariant under the fixed layer swap.
- Do not execute an action evaluator, numerical derivative, eigensolver,
  full Hessian, pseudoinverse, root search, spectrum or physical target.

## 5. Frozen outcomes

1. `ADVERSARIAL_REFINED_H4_NULL_COUPLING_CONTROL_FAILED` for provenance,
   symbolic, actual-incidence, vertical-balance, corruption or scope failure;
2. `ADVERSARIAL_REFINED_H4_NULL_COUPLING_DISAGREEMENT` if (2) disagrees with
   any primary compatibility component or rank/reversal statement;
3. `ADVERSARIAL_REFINED_H4_NULL_COUPLING_CORROBORATED` otherwise.

Outcome 3 accepts the mechanically independent statement that the internal
product-lapse line couples nontrivially to one schedule-independent boundary
compatibility covector.  It does not compute the constrained effective
Hessian or derive a physical tick, `c`, `G` or Planck scale.

Register before execution, run twice with byte-identical JSON, and perform
only the static registry audit.  Do not run the full suite.

