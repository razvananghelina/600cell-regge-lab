# Protocol: scale covariance and the absolute-tick no-go

Date: 2026-08-21

Prior-art gate:
`docs/gravity/gravity_600cell_tick_scale_covariance_prior_art.md`.

Status: frozen before evaluating any rescaled complete-action state.

Only the new targeted verifiers and the static registry audit may be run.  Do
not run the full suite or any deferred nonlinear census.

## 1. Frozen model and transformation

Use the two already-derived order-24 staircase parities and the complete
100-decimal 95-variable Lorentzian Regge-plus-dust evaluator underlying the
canonical continuation.  Its inputs are 30 old-boundary squared lengths, 35
internal magnitudes (30 cross-edge squares and five positive pole squares), 30
new-boundary squared lengths, and the geometrized total dust mass.

For positive `alpha`, define `r=alpha^2` and transform

```text
(q_old, x, q_new, M) -> (r q_old, r x, r q_new, alpha M).
```

Use exactly `alpha in {3/5, 7/4}`.  No scale may be selected from output.

## 2. Off-shell state fixed before evaluation

Start from the unrounded published state.  Make it nonstationary and
nonhomogeneous using the deterministic dimensionless factors

```text
q_old[i] *= exp(1e-6*((i mod 7)-3)),             i=0..29,
x[i]     *= exp(1e-6*((i mod 5)-2)),             i=0..34,
q_new[i] *= exp(1e-6*((i mod 11)-5)),            i=0..29.
```

Do not solve any equation.  Require the base and every transformed state to
remain on the same certified Lorentzian branch.  The perturbation prevents a
stationary zero from making derivative covariance vacuous.

## 3. Exact algebraic certificates

Using exact SymPy algebra, certify:

1. the triangle area-square polynomial has degree two in squared lengths;
2. a 4-simplex Gram determinant, facet volume square, hinge area square and
   the dihedral sine/cosine formulas have scale degrees `4,3,2,0` in `r` as
   required for angle invariance;
3. the gravitational hinge term has length degree two;
4. `m sqrt(rho)` has length degree two only when `m -> alpha m`;
5. the coarse mass `M=(90/pi) epsilon3 L` and the refined local rule
   `m_v=K_v/(8*pi)` both have length degree one.

The refined rule may be certified from the definition of three-dimensional
Regge scalar curvature `K_v=sum_(edges at v) length*deficit`, whose angles are
scale invariant.  No stored numerical match is sufficient for this gate.

## 4. Primary complete-action test

At 100 decimal digits, for both parities and both frozen `alpha` values,
evaluate the complete action and all 95 logarithmic squared-length derivatives
at the base and transformed states.  Require

```text
S_scaled = r S_base,
g_scaled = r g_base
```

componentwise with normalized error at most `1e-65`.  Require at least 20
internal components and at least 20 boundary components of `g_base` to exceed
`1e-20`, so the test cannot pass on an accidentally sparse stationary vector.

Record branch data, maximum action error, maximum component error, component
counts, both scale factors, and both parities.

## 5. Mandatory hostile control

Repeat each transformed evaluation while holding `M` fixed.  The geometry is
still scaled, but the dust term now has degree one in `alpha`, not degree two.
Require the claimed total-action covariance to fail by more than `1e-8` and at
least one pole derivative covariance error to exceed `1e-8` for every parity
and scale.  If it does not, the test lacks falsification power and the outcome
is a control failure.

Also certify symbolically that the fixed-mass defect is

```text
(alpha-alpha^2) S_dust_base,
```

up to the sign convention used to compare scaled and predicted actions.

## 6. Mechanically different adversarial audit

A separately registered verifier must not call the 100-decimal orbit action.
It must construct the direct 2400-simplex binary64 action using the existing
full-simplex evaluator, add the dust term independently, and use derivatives
with respect to raw squared-length variables.  It must independently construct
all scaled states before reading the primary result.

For both parities and both scales require, within propagated binary64 bounds,

```text
S_scaled = r S_base,
dS_scaled/dq = dS_base/dq.
```

It must execute the same fixed-mass negative control.  Agreement of orbit and
direct results is required before accepting the scientific statement.

## 7. Outcome hierarchy

The primary verifier assigns exactly one:

1. `TICK_SCALE_COVARIANCE_CONTROL_FAILED` if provenance, algebra, branch,
   nonzero-support or hostile-control gates fail;
2. `TICK_SCALE_COVARIANCE_REFUTED` if any simultaneous geometry-and-mass
   scaling identity fails;
3. `TICK_SCALE_COVARIANCE_PRIMARY_CONFIRMED` otherwise.

The adversarial verifier assigns exactly one:

1. `TICK_SCALE_COVARIANCE_ADVERSARIAL_CONTROL_FAILED` for a failed independent
   construction or hostile control;
2. `TICK_SCALE_COVARIANCE_IMPLEMENTATIONS_DISAGREE` if its covariance result
   conflicts with the primary artifact;
3. `ABSOLUTE_CLASSICAL_TICK_NO_GO_CORROBORATED` only if every independent gate
   passes.

## 8. Interpretation boundary

Outcome 3 of the adversarial audit establishes the following conditional
statement:

> **DERIVED EXACT / ADVERSARIALLY CORROBORATED:** on the fixed carrier, for the
> stated zero-cosmological-constant classical Regge-plus-dust action with all
> geometrized masses scaled with the geometry, stationary and canonical
> solutions occur in global scale families.  The equations cannot select an
> absolute nonzero tick.

It does not exclude selection of `tau/L`, `tau_next/tau0`, relational dust
time, or a duration relative to an externally fixed mass or cosmological
scale.  It does not derive `c`, `G`, `hbar`, Planck units, or a quantum scale.

