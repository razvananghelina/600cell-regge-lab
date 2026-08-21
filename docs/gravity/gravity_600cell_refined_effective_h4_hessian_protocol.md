# Protocol: on-shell effective H4 boundary Hessian across all schedules

Date: 2026-08-21

Prior-art commit: `f7bf3c1`.

This protocol is frozen before constructing or inspecting any new full or
effective Hessian on the curvature-matched matter branch.

## 1. Frozen inputs

Require these exact files and hashes:

```text
reproducible/verify_gravity_600cell_refined_h4_stationary_fill.py
  89aab727792e20a81e7577e0425f8fa4b1e84e2a7ae66caa9e79a4aebf3581e7
reproducible/gravity_600cell_refined_local_curvature_mass.json
  180010a79177ba16620ebea9847443c57a7a6d2d8a3df71ad6ecb83f454ef091
reproducible/gravity_600cell_refined_local_curvature_mass_adversarial.json
  c59890d12bf929c4677dffed1b932ad8c05ab0ac00980be15ba780e62744c28e
reproducible/gravity_600cell_refined_boundary_cotangent.json
  4e7bf0beb0327a3ee1bddbec13126fbef99380970e62cecf74eb24ce8d6dafaa
docs/gravity/gravity_600cell_refined_boundary_cotangent_result.md
  391a317b9f8823a5479f450dde43a43177e210a2d81192aedc938e90fc8006d1
docs/gravity/gravity_600cell_refined_effective_h4_hessian_prior_art.md
  d111b896265ccbd0534ec50fec184c81067133665b97a99ac0a69df834877934
commons/cell600.py
  ea5bce4b6c52e0834539ca4b1df9c6a67a3a5ed4da32f4e0298a493fc5315c7f
```

Require the accepted curvature-mass and boundary-cotangent outcomes.  Load
only function definitions from the frozen action source; executing its
top-level P1 calculation is forbidden.

## 2. Coordinates and exact background

Use total-orbit logarithmic squared-edge coordinates in fixed order

```text
B=(old_01,old_02,old_03,old_12,old_13,old_23,
   new_01,new_02,new_03,new_12,new_13,new_23),
I=(cross_01,cross_02,cross_03,cross_12,cross_13,cross_23,
   rho_0,rho_1,rho_2,rho_3).
```

Rebuild `K0=P(sd K_600)` and all 24 schedules directly.  Use the exact static
product coordinates at `tau0=0.0102`.  Set the action evaluator's legacy P1
mass to zero.  Add the already frozen selected masses only through

```text
H_dust[rho_r,rho_r]=-2*pi*m_r*tau0.               (1)
```

No other Hessian entry receives a dust term.  Before forming a Hessian,
recheck the ten total internal gradients by adding
`-4*pi*m_r*tau0` to the four lapse rows; require maximum residual `<1e-60`
for every schedule.

## 3. High-precision Hessian construction

At `100` and `140` decimal digits, differentiate the analytic 22-component
gravitational log-gradient by centred log-coordinate differences.  Freeze

```text
h0=1e-10, h1=5e-11, h2=2.5e-11.
```

At 100 digits construct

```text
H100a=(4*D(h1)-D(h0))/3,
H100b=(4*D(h2)-D(h1))/3,
```

and at 140 digits construct `H140b` from `(h1,h2)`.  Add (1) independently
to all three matrices, then use the real symmetric part of `H140b` as the
reported Hessian.

For each schedule define

```text
e_H = 100*max(||H100a-H100b||max,
              ||H100b-H140b||max)
      + 1e-50*max(1,||H140b||max).                (2)
```

Require raw imaginary part and raw antisymmetry no larger than `e_H`.

## 4. Internal elimination and certified envelopes

Let `lambda_j` be the eigenvalues of the real symmetric `H_ii`.  Define the
conservative spectral envelope `e_lambda=100*e_H`.  Do not invert unless

```text
min_j |lambda_j| > e_lambda.                      (3)
```

If (3) holds, solve linear systems rather than constructing a numerical
pseudoinverse and form

```text
K=H_bb-H_bi solve(H_ii,H_ib).                     (4)
```

Apply (4) separately to `H100a`, `H100b` and `H140b`, and define

```text
e_K = 100*max(||K100a-K100b||max,
              ||K100b-K140b||max)
      + 1e-45*max(1,||K140b||max).                (5)
```

Require relative linear-solve residual `<1e-60` and raw Schur antisymmetry no
larger than `e_K`.  Store all 24 full `12x12` matrices, internal eigenvalues,
inertias, condition numbers and envelopes.

Two matrices are equal only when their maximum entrywise difference is at
most the larger of their two `e_K` values.  This relation, applied in fixed
lexicographic schedule order, defines the reported class census.  Print the
complete class membership before assigning the outcome.

## 5. Time reversal and controls

Let `R` be the fixed `12x12` permutation exchanging the six old and six new
coordinates.  For each schedule `o`, require

```text
K_o = R^T K_reverse(o) R
```

inside the same pairwise envelope.  This is the only allowed identification;
no rank-colour relabelling is permitted.

Controls:

1. A synthetic block matrix
   `Hbb=[[5,1],[1,4]]`, `Hbi=[[1],[2]]`, `Hii=[[2]]` must give exactly
   `K=[[4.5,0],[0,2]]`.
2. The difference between the selected-matter and gravity-only full Hessians
   must be zero off the four lapse diagonal entries and agree with (1) on
   them within `1e-80`.
3. For schedules with lexicographic indices `0,1,22,23`, use boundary
   directions `old_01`, common old/new scale, and old-minus-new scale.  Lift
   each direction by `u=-H_ii^(-1)H_ib v`.  Centred complete-action second
   differences at `1e-10` and `5e-11`, Richardson combined, must reproduce
   `v^T K v` to relative error `<1e-28`.
4. Adding `1e-6*max(1,||K||max)` to one component of a copied matrix must be
   detected as inequivalent by the frozen class comparator.
5. The calculation must not execute a root search, nested census, spectrum,
   continuum comparison or physical-constant extraction.

## 6. Frozen outcomes

Use the first applicable outcome:

1. `REFINED_EFFECTIVE_H4_HESSIAN_CONTROL_FAILED` for provenance, topology,
   on-shell, precision, reality, symmetry, solve, directional or corruption
   control failure;
2. `REFINED_EFFECTIVE_H4_HESSIAN_INTERNAL_SINGULAR` if any schedule fails
   (3);
3. `REFINED_EFFECTIVE_H4_HESSIAN_TIME_REVERSAL_FAILED` if the fixed reversal
   covariance fails;
4. `REFINED_EFFECTIVE_H4_HESSIAN_MULTIPLE_SCHEDULE_CLASSES` if more than one
   effective class remains;
5. `REFINED_EFFECTIVE_H4_HESSIAN_SINGLE_SCHEDULE_CLASS` otherwise.

Outcome 4 is a **DERIVED COMPUTATIONAL NEGATIVE** for a canonical bare
staircase evolution in the `H4` sector, subject to adversarial replication.
Outcome 5 advances only to a nonhomogeneous quadratic operator; it does not
derive propagation or `c`.  Any result with physical weight remains
unaccepted until a mechanically different adversarial test.

## 7. Execution boundary

Register the verifier before its first run.  Run it twice and require a
byte-identical artifact.  Perform the static registry audit only.  Do not run
the full suite or the deferred 12-case nonlinear census.

