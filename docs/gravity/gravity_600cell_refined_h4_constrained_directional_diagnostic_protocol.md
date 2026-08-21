# Protocol: constrained-response directional truncation diagnostic

Date: 2026-08-21

This correction diagnostic is frozen after the formal `18/19` primary control
failure in commit `ba828ec` and before evaluating any additional action
difference.

## 1. Frozen question

The primary constrained Hessian/action comparison used one fourth-order
Richardson estimate at steps `1e-10` and `5e-11` and missed its copied
relative `1e-28` threshold, with maximum error `1.3085e-20`.  The lifted
directions have maximum components up to approximately `1.545e5`, so the
largest coordinate displacement was approximately `1.545e-5`.

Determine whether the discrepancy is the expected finite-step truncation of
the centred action difference or a genuine mismatch with the stored
constrained quadratic response.  Do not recompute or modify any Hessian,
lift, response matrix or schedule class.

## 2. Frozen inputs

Require exact hashes for:

```text
reproducible/verify_gravity_600cell_refined_h4_stationary_fill.py
  89aab727792e20a81e7577e0425f8fa4b1e84e2a7ae66caa9e79a4aebf3581e7
reproducible/gravity_600cell_refined_local_curvature_mass.json
  180010a79177ba16620ebea9847443c57a7a6d2d8a3df71ad6ecb83f454ef091
reproducible/gravity_600cell_refined_h4_constrained_response.json
  f029260c9ee6e3b763293d237aae27e6ff7c1256eb8bc19c35725084ff385888
docs/gravity/gravity_600cell_refined_h4_constrained_response_primary_first_result.md
  633a57f3d2b4a054cce20d08544d409dac8fdaf53c39bae72ab2e9fceb4e83eb
commons/cell600.py
  ea5bce4b6c52e0834539ca4b1df9c6a67a3a5ed4da32f4e0298a493fc5315c7f
```

Require the frozen primary outcome `CONTROL_FAILED`, `18/19`, one tentative
schedule class and the same twelve directional records.  Load only action
function definitions; do not import or execute either Hessian verifier.

## 3. Multi-level action ladder

For each of the same 12 directions (schedule indices `0,1,22,23`, times the
first-basis, all-ones and alternating coefficient vectors), reconstruct the
full 22-component tangent only from the stored boundary basis and stored
primary lift.  Let its stored quadratic value be `q`.

At both `140` and `180` decimal digits evaluate

```text
D_j=[S(+h_j v)-2S(0)+S(-h_j v)]/h_j^2,
h_j=1e-10/2^j,  j=0,...,4.
```

Build

```text
R_j=(4 D_{j+1}-D_j)/3,            j=0,...,3,
X_j=(16 R_{j+1}-R_j)/15,          j=0,...,2,
Y_j=(64 X_{j+1}-X_j)/63,          j=0,1.          (1)
```

For an analytic action, `D`, `R`, `X`, `Y` successively remove the `h^2`,
`h^4` and `h^6` errors.  Use `Y_1` at 180 digits as the final action estimate.

Define separately for each direction

```text
e_dir=100*max(|Y_0,180-Y_1,180|,
              |Y_1,140-Y_1,180|)
      +1e-50*max(1,|Y_1,180|).                   (2)
```

Require `|Y_1,180-q|<=e_dir`.  No fixed replacement tolerance is allowed.

## 4. Truncation and precision gates

For each direction require:

1. `|R_0-q|>|R_1-q|>|R_2-q|>|R_3-q|`;
2. the three successive ratios `|R_j-q|/|R_(j+1)-q|`, `j=0,1,2`, lie in
   `[8,32]`, bracketing the theoretical fourth-order factor `16`;
3. `|X_0-q|>|X_1-q|>|X_2-q|`;
4. raw imaginary parts of all estimates are below (2);
5. the 140/180 difference is included explicitly in (2).

These gates must pass for all twelve directions; a peak or median is
insufficient.

## 5. Independent controls

1. Apply the same ladder to the exact even polynomial
   `f(x)=7*x^2/2+11*x^4+13*x^6+17*x^8+19*x^10`; recover the exact second
   derivative `7` within its own (2).
2. For every real direction replace `q` by
   `q_bad=q+1e-12*max(1,|q|)`.  Require
   `|Y_1-q_bad|>10^6 e_dir`.
3. Require every maximum coordinate displacement
   `h_0*||v||max` below `2e-5` and every finest displacement
   `h_4*||v||max` below `2e-6`; this records the regime actually tested without
   changing it.
4. No Hessian, solve, schedule census, root search, spectrum, continuum target
   or physical constant is computed.

## 6. Frozen outcomes

Use the first applicable outcome:

1. `REFINED_H4_DIRECTIONAL_DIAGNOSTIC_CONTROL_FAILED` for provenance,
   reconstruction, polynomial, precision, imaginary, corruption or scope
   failure;
2. `REFINED_H4_DIRECTIONAL_DIAGNOSTIC_NONASYMPTOTIC` if any monotonicity or
   ratio gate fails;
3. `REFINED_H4_DIRECTIONAL_HESSIAN_ACTION_MISMATCH` if the asymptotic gates
   pass but any final estimate fails (2);
4. `REFINED_H4_DIRECTIONAL_TRUNCATION_CONFIRMED` otherwise.

Outcome 4 licenses only a separately preregistered replacement of the primary
directional check by (1)--(2), followed by a complete rerun of that targeted
verifier.  It does not retroactively turn the first failed artifact into an
accepted result.

## 7. Execution

Fill the committed first-result hash and this protocol hash, register before
execution, run twice and require a byte-identical artifact.  Do not run the
full suite or deferred nonlinear census.
