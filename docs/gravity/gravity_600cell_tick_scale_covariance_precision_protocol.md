# Precision-adjudication protocol for the direct tick covariance audit

Date: 2026-08-21

Upstream protocol: `fd1f8a8`.

Recorded first adversarial failure: commit `9951c7a`.

Status: frozen before any arbitrary-precision direct 2400-simplex evaluation.

Only the new targeted verifier and the static registry audit may be run.  The
failed binary64 verifier and its artifact must remain unchanged.

## 1. Exact discrepancy to adjudicate

Both the 100-decimal orbit evaluator and the binary64 direct 2400-simplex
evaluator separately obeyed the frozen scaling identity for both parities and
both scale factors.  The sole blocking gate is their base-action discrepancy:

```text
even: 1.648e-6,
odd : 1.338e-6,
```

against an inherited `2e-8` cross-implementation threshold.  The direct
evaluator's separately frozen conditioning-aware bounds were of order `1e-6`.
No threshold in the failed verifier may be changed.

## 2. Mechanically distinct precision construction

At 80 decimal digits, reconstruct the same off-shell state directly from the
published formulas and frozen perturbation.  Do not call either the primary
orbit action or the binary64 `full_evaluation` function.

Reimplement the following locally:

1. signed simplex volumes from Cayley/Gram determinants;
2. Lorentzian dihedral sine, cosine and branch logarithm;
3. one curvature accumulator keyed by every individual triangle;
4. a literal loop over all 2400 four-simplices;
5. a literal sum over all individual hinges, with no orbit multiplicity;
6. the point-dust action from the independently reconstructed mass.

The carrier and edge lookup may be loaded from the certified slab
reconstruction.  Reusing its combinatorics is declared shared upstream, not
mechanical independence of the geometry.

## 3. Frozen evaluations

For each of the even and odd parities evaluate exactly three states:

```text
alpha = 1,
alpha = 3/5,
alpha = 7/4,
```

with `q -> alpha^2 q` and `M -> alpha M`.  Construct and evaluate all six
states before reading either covariance artifact.

Every state must have exactly 2400 Lorentzian simplices with one negative Gram
direction, minimum nonzero leading-minor modulus above `1e-20`, minimum angle
argument modulus above `1e-6`, and maximum imaginary total-action part below
`1e-60`.

## 4. Frozen acceptance gates

Within each parity require

```text
relative_error(S_alpha, alpha^2 S_1) < 1e-55
```

for both nontrivial scales.  Require the two direct high-precision base actions
to agree with each other within `1e-55`.

Only after those values are complete, load the stored artifacts.  Require:

1. the primary outcome is `TICK_SCALE_COVARIANCE_PRIMARY_CONFIRMED`;
2. the first adversarial outcome is
   `TICK_SCALE_COVARIANCE_IMPLEMENTATIONS_DISAGREE` with `10/12`;
3. each high-precision direct base action agrees with the corresponding stored
   100-decimal orbit action within `1e-45` (the stored string has 50 digits);
4. its discrepancy from the stored binary64 direct action reproduces the
   recorded `1e-6` disagreement within a factor of two;
5. high precision improves the orbit/direct agreement by at least `10^30`
   relative to the binary64 discrepancy.

The factor-of-two gate concerns reproduction of a recorded binary64 error, not
a physics equality.  The `10^30` improvement requirement prevents a marginal
threshold change from resolving the disagreement.

## 5. Outcome hierarchy

Assign exactly one:

1. `TICK_SCALE_COVARIANCE_PRECISION_CONTROL_FAILED` if provenance, state,
   branch, artifact or error-reproduction controls fail;
2. `TICK_SCALE_COVARIANCE_HIGH_PRECISION_DISAGREES` if the direct high-precision
   action fails scaling or fails to match the orbit action;
3. `ABSOLUTE_CLASSICAL_TICK_NO_GO_PRECISION_ADJUDICATED` only if every gate
   passes.

Outcome 3 resolves the earlier disagreement as binary64 absolute-action loss
and, together with the exact degree proof and independent raw-gradient audit,
licenses the conditional statement:

> **DERIVED EXACT / ADVERSARIALLY CORROBORATED:** the stated classical
> zero-cosmological-constant Regge-plus-dust theory cannot select an absolute
> nonzero tick when all geometrized masses scale with the geometry.

This still does not exclude dimensionless or relational time, nor time measured
relative to a separately justified dimensionful input.

