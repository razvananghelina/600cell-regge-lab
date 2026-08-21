# Adversarial protocol: derive the boundary covector from spatial hinges

Date: 2026-08-21

Primary-result commit: `f805557`.

This protocol is frozen before writing or executing the adversarial verifier.
It must not import or execute the primary boundary-cotangent verifier or the
Lorentzian stationary-fill action evaluator.

## 1. Frozen inputs

Require:

```text
reproducible/verify_gravity_600cell_refined_boundary_cotangent.py
  ababad0e8e667e31c290b9e8bbf61005308faed20af4d09ae7affbc32b3509d7
reproducible/gravity_600cell_refined_boundary_cotangent.json
  4e7bf0beb0327a3ee1bddbec13126fbef99380970e62cecf74eb24ce8d6dafaa
docs/gravity/gravity_600cell_refined_boundary_cotangent_primary_result.md
  ec69cd7c6521b3cff3ded777d80a4c740b9065bcce73bd52f53c0591433c9074
reproducible/gravity_600cell_refined_local_curvature_mass_adversarial.json
  c59890d12bf929c4677dffed1b932ad8c05ab0ac00980be15ba780e62744c28e
reproducible/gravity_600cell_dust_regular_lapse_identity.json
  5079428fade247f730ebc07e5e2eae388b48045cd5201e84afb3186bfc248a51
commons/cell600.py
  ea5bce4b6c52e0834539ca4b1df9c6a67a3a5ed4da32f4e0298a493fc5315c7f
```

Require the primary `16/16`
`REFINED_BOUNDARY_COTANGENT_SELECTED_RENORMALIZED` outcome and the accepted
actual-incidence curvature artifact.

## 2. Independent local product identity

For one spatial edge of length `l` and product proper time `tau`, the lower
vertical hinge triangle has squared edge data

```text
(x,y,z)=(l^2,-tau^2,l^2-tau^2).
```

Starting from Heron's squared-area polynomial

```text
A^2=[2(xy+xz+yz)-x^2-y^2-z^2]/16,
```

derive symbolically, while holding `y,z` fixed under the old-boundary
variation,

```text
A=i*l*tau/2,
dA/dlog(x)=i*l*tau/4.
```

The Schlaefli identity removes angle derivatives.  With the repository's
`-i A epsilon` Lorentzian convention, a lower boundary edge therefore
contributes

```text
G_old,e=+tau*l_e*epsilon_e/4,
P_pre,e=-tau*l_e*epsilon_e/4.                    (1)
```

The upper product triangle gives the opposite canonical sign,
`P_post,e=+tau*l_e*epsilon_e/4`.

As a known-answer control, (1) on the unrefined regular 600-cell must
reproduce the frozen exact per-edge formula
`P_pre=-L*epsilon3*tau/4`.

## 3. Actual-incidence six-vector

Independently rebuild the 600-cell, its barycentric flag complex and all
17,040 actual refined edges.  Compute every edge's tetrahedron incidence,
exact rank-pair length and deficit, then group

```text
C_rs=sum_(edge rank pair rs) l_e epsilon_e.
```

Without reading the primary vector during construction, form

```text
P_pre,rs=-tau0*C_rs/4,
P_post,rs=+tau0*C_rs/4.                          (2)
```

Only after writing the six values compare (2) componentwise with the frozen
primary vector.  Require maximum absolute error below `1e-68`.

Require also:

- `2 sum P_pre=-tau0 sum C_rs/2`;
- the raw coarse/fine ratio from (2) equals the independently reconstructed
  curvature ratio within `1e-68`;
- momentum per curvature-selected mass equals `-4*pi*tau0` within `1e-68`.

## 4. Controls and scope

- Dropping one actual edge from one `C_rs` must change the corresponding
  component by more than `1e-8` relative to the complete vector.
- Replacing the factor `1/4` by `1/2` must fail both the primary-vector and
  pullback comparisons.
- Reversing the canonical sign must fail the pre-vector comparison.
- A kernel perturbation must retain the scalar pullback while disagreeing
  with the six action-selected components.

Record that no primary function, Lorentzian action evaluator, root search,
nested census, Hessian, spectrum or physical target is executed.

## 5. Frozen outcomes

Use the first applicable outcome:

1. `ADVERSARIAL_REFINED_BOUNDARY_COTANGENT_CONTROL_FAILED` for provenance,
   topology, symbolic-identity, known-answer, corruption or scope failure;
2. `ADVERSARIAL_REFINED_BOUNDARY_COTANGENT_DISAGREEMENT` if the independent
   six-vector or normalized comparison disagrees with the primary;
3. `ADVERSARIAL_REFINED_BOUNDARY_COTANGENT_CORROBORATED` otherwise.

Outcome 3 accepts **DERIVED COMPUTATIONAL / STRUCTURAL** selection of the
renormalized refined covector.  It does not establish an exact fixed-radius
perfect action, refinement convergence, a physical tick, `c`, `G`, Planck
units or particle masses.

## 6. Provenance warning and deliverables

The component formula (1) was recognized after the primary vector was known.
Its evidential role is an analytic/mechanically independent explanation, not
a blind prediction.  Register the verifier before execution, run it twice
with byte-identical JSON and perform only the static registry audit.  Do not
run the full suite.
