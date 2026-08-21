# Preregistered correction: primary constrained directional control

Date: 2026-08-21

The formal first primary failure is preserved in commit `ba828ec`.  The
directional diagnosis and its implementation failures are preserved through
commit `956ccf2`.  The frozen combined-envelope audit is accepted in commit
`3589610`.

## Frozen new inputs

Add exact provenance checks for

```text
reproducible/gravity_600cell_refined_h4_constrained_combined_envelope.json
  34e2d598a6f608c9436217024138b32f0095e5df8e32d3ff91df2b182843aa0d
docs/gravity/gravity_600cell_refined_h4_constrained_combined_envelope_result.md
  7439a9707df4531b506f736e6b59b7dc292d939e9bacd43d347616a0063536d1
```

Require outcome `REFINED_H4_COMBINED_ENVELOPE_CORROBORATED`, `10/10`, and
`12/12` combined matches.

## Allowed code change

Replace only the primary verifier's direct-action directional block and its
serialized records.  Keep the same schedules `0,1,22,23`, coefficient
directions, computed primary boundary basis, lifts and quadratic responses.

At both 140 and 180 digits evaluate centred complete-action differences at

```text
h_j=1e-10/2^j, j=0,...,4,
```

and form

```text
R_j=(4D_(j+1)-D_j)/3,
X_j=(16R_(j+1)-R_j)/15,
Y_j=(64X_(j+1)-X_j)/63.
```

Require, independently of the quadratic target, both successive `R`
difference ratios in `[8,32]` and the `X` difference ratio in `[32,128]` for
all twelve directions.

For each direction set

```text
e_action=100*max(|Y_0,180-Y_1,180|,
                 |Y_1,140-Y_1,180|)
         +1e-50*max(1,|Y_1,180|),
e_hessian=||y||_1^2 e_K,
e_total=e_action+e_hessian.
```

Require

```text
|Y_1,180-y^T K y| <= e_total.
```

This exact componentwise propagation replaces the impossible fixed relative
`1e-28` gate.  Store all ladders, ratios and the two uncertainty contributions.

Do not change any Hessian step, matrix envelope, basis, lift, solve, reversal,
class comparator, corruption control, outcome hierarchy or scope.  Do not
load a target or select/average a schedule.

Preserve the frozen failed artifact at
`gravity_600cell_refined_h4_constrained_response.json`, because both correction
audits require its exact hash.  The repaired primary verifier must write
`gravity_600cell_refined_h4_constrained_response_corrected.json` instead.
This is a provenance-preserving output rename, not a scientific change.

## Execution

Commit the implementation before rerunning.  Run the complete targeted
primary verifier twice and require a byte-identical artifact.  Do not run the
full suite or deferred nonlinear census.  A clean single-class primary result
still requires a mechanically different adversarial audit before acceptance.
