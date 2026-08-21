# Primary result: the refined action selects a renormalized boundary covector

Date: 2026-08-21

Status: primary result; not accepted before adversarial replication.

## Provenance

| stage | commit |
|---|---|
| prior-art gate | `e7a1545` |
| frozen protocol | `5ee4f3e` |
| verifier registered before execution | `56b4ea6` |
| first `15/16` control failure preserved | `5ed87f6` |
| narrow control correction | `8991447` |
| corrected implementation | `a0edc31` |

The corrected targeted verifier passed `16/16` twice and wrote byte-identical
artifact

```text
reproducible/gravity_600cell_refined_boundary_cotangent.json
SHA-256 4e7bf0beb0327a3ee1bddbec13126fbef99380970e62cecf74eb24ce8d6dafaa.
```

No full suite, nested root census, Hessian or spectrum was run.

## Complete six-component result

At the curvature-matched on-shell static seed, in orbit-total logarithmic
squared-edge coordinates `(01,02,03,12,13,23)`, every one of the 24 staircase
schedules gives

```text
P_pre =
(-0.0365507943663027407026...,
 -0.00243609423423285173839...,
 -0.0000655036546510496875794...,
 -0.07335464323435491928399...,
 -0.000989979010535868633659...,
 -0.03822370401923157424126...).

P_post=-P_pre.
```

The maximum component spread over schedules is `1.91e-96`.  All internal
cross residuals are below `1.87e-96`, all curvature-mass vertical residuals
below `5.91e-77`, and independent Richardson boundary derivatives agree to
`1.73e-42` relative.

Thus the action selects one point in the five-dimensional inverse fiber left
open by geometry and symplecticity alone.  A deliberately kernel-shifted
six-vector has the same scalar pullback and proves that the component result,
not the scalar sum, carries this information.

## Coarse/fine distinction

The homothetic pullback satisfies

```text
p_pre,fine=-tau0*K_fine/2=-4*pi*tau0*M_fine,
p_post,fine=+tau0*K_fine/2=+4*pi*tau0*M_fine
```

to maximum absolute error `4.32e-76`.

At fixed unit volume radius, however,

```text
p_pre,fine/p_pre,coarse
 = K_fine/K_coarse
 = M_fine/M_coarse
 = 0.984190377388521915998...,
```

so the bare raw momentum differs by about `1.58096%`.  After dividing each
momentum by its own curvature-selected mass,

```text
(p_fine/M_fine)/(p_coarse/M_coarse)=1
```

within the numerical envelope.

The frozen outcome is

```text
REFINED_BOUNDARY_COTANGENT_SELECTED_RENORMALIZED.
```

## Interpretation before replication

- **PRIMARY DERIVED COMPUTATIONAL / STRUCTURAL:** the on-shell refined action
  selects a complete schedule-independent boundary covector.
- **DERIVED NEGATIVE:** the bare fixed-radius coarse and fine covectors are
  not equal.
- **DERIVED PATTERN/IDENTITY CANDIDATE:** the mismatch is exactly the spatial
  curvature/mass renormalization, while momentum per selected mass is
  invariant.
- **NOT ESTABLISHED:** a perfect action, physical coarse/fine equivalence,
  refinement convergence, a tick, `c`, `G`, Planck units or particle masses.

The adversarial gate must derive the six components from actual spatial hinge
curvatures rather than reevaluate the same Lorentzian action code.
