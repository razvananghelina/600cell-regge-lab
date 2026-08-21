# Protocol: product-lapse null line and boundary compatibility covector

Date: 2026-08-21

Prior-art commit: `32854ba`.

This protocol is frozen before evaluating any new null image or compatibility
covector.

## 1. Frozen inputs

Require:

```text
reproducible/verify_gravity_600cell_refined_h4_stationary_fill.py
  89aab727792e20a81e7577e0425f8fa4b1e84e2a7ae66caa9e79a4aebf3581e7
reproducible/gravity_600cell_refined_local_curvature_mass.json
  180010a79177ba16620ebea9847443c57a7a6d2d8a3df71ad6ecb83f454ef091
reproducible/gravity_600cell_refined_local_curvature_mass_adversarial.json
  c59890d12bf929c4677dffed1b932ad8c05ab0ac00980be15ba780e62744c28e
reproducible/gravity_600cell_refined_boundary_cotangent.json
  4e7bf0beb0327a3ee1bddbec13126fbef99380970e62cecf74eb24ce8d6dafaa
reproducible/gravity_600cell_refined_effective_h4_hessian.json
  56e08db9a840b95e686fadb2763e89400b09220e88b80e9d35c17c1e73eef0a3
docs/gravity/gravity_600cell_refined_effective_h4_hessian_first_result.md
  f8bf5679e153fcca8a076064bc5b98e881d91c2add9aaffa4c6858247538f1b8
docs/gravity/gravity_600cell_refined_h4_null_coupling_prior_art.md
  4a6b55689535feeb729db3a71a9a28d2a7a8ccbeb6554a79371f1073b6d669dd
commons/cell600.py
  ea5bce4b6c52e0834539ca4b1df9c6a67a3a5ed4da32f4e0298a493fc5315c7f
```

Require all accepted upstream outcomes and the primary Hessian census
`24 x (9,1,0)`.  Load only function definitions from the frozen action
source; do not execute its P1 top level or import the full-Hessian verifier.

## 2. Exact tangent and product-family control

In internal order

```text
(cross_01,cross_02,cross_03,cross_12,cross_13,cross_23,
 rho_0,rho_1,rho_2,rho_3),
```

construct at the exact static product

```text
n=(-tau0^2/q_cross,01,...,-tau0^2/q_cross,23,1,1,1,1).  (1)
```

Do not obtain `n` from an eigensolver.  At `tau/tau0` equal to `1/2`, `1`
and `2`, rebuild the exact product cross and lapse coordinates while keeping
the spatial boundaries and selected masses fixed.  Require every one of the
ten internal gradients below `1e-60` for all 24 schedules and all three
ratios.  This tests the finite family whose tangent is (1), rather than only
one infinitesimal cancellation.

## 3. Directional Hessian image

Set the legacy P1 mass to zero.  At `100` and `140` decimal digits,
differentiate the complete 22-component gravitational log-gradient along
the fixed internal direction (1), using exponential log-coordinate shifts
and

```text
h0=1e-10, h1=5e-11, h2=2.5e-11,
w100a=(4*D(h1)-D(h0))/3,
w100b=(4*D(h2)-D(h1))/3,
w140b=(4*D(h2)-D(h1))/3 at 140 digits.
```

Add the selected dust Hessian image only on the four lapse rows:

```text
(H_dust n)_rho,r=-2*pi*m_r*tau0.                  (2)
```

Define

```text
e_w=100*max(||w100a-w100b||max,
            ||w100b-w140b||max)
    +1e-50*max(1,||w140b||max).                   (3)
```

Split `w140b=(c,H_ii n)` into 12 boundary and ten internal components.
The analytic null prediction passes only if

```text
||H_ii n||max <= e_w.                             (4)
```

Since the frozen internal nullity is exactly one, (4) then identifies the
entire internal null line without fitting an eigenvector.

## 4. Boundary prediction fixed before comparison

For each schedule evaluate its base boundary action-gradient vector

```text
g_B=(g_old01,...,g_old23,g_new01,...,g_new23).
```

Print the computed `c` before loading the frozen boundary-cotangent vector.
Then require componentwise

```text
c=(1/2)*g_B                                           (5)
```

inside (3).  Independently reconstruct

```text
g_B=(-P_pre,P_post)
```

from the frozen primary boundary artifact and require the same comparison.
Require `||c||max > 10^6 e_w`, so a numerical zero cannot pass as coupling.

The 24 unnormalised compatibility rows must have rank one in the explicit
sense that all are equal within pairwise envelopes and the common row is
nonzero.  For every order `o`, require the fixed old/new layer-swap relation
between `c_o` and `c_reverse(o)`.  No rank-colour permutation is allowed.

## 5. Controls and scope

1. A synthetic affine gradient `g(z)=A z+b` with a fixed rational `22x22`
   symmetric `A` must reproduce `A*(0,n)` under the same directional ladder
   to error `<1e-80`.
2. Replacing `1/2` in (5) by `1` and by `-1/2` must fail by more than
   `10^6 e_w`.
3. Dropping the largest component of the predicted compatibility row must
   fail the componentwise comparison by more than `10^6 e_w`.
4. Raw imaginary components must be below (3).
5. No full Hessian, pseudoinverse, Schur complement, root search, nested
   census, spectrum, continuum target or physical constant is computed.

## 6. Frozen outcomes

Use the first applicable outcome:

1. `REFINED_H4_NULL_COUPLING_CONTROL_FAILED` for provenance, product-family,
   precision, reality, synthetic, corruption or scope failure;
2. `REFINED_H4_PRODUCT_TANGENT_NOT_NULL` if (4) fails;
3. `REFINED_H4_PRODUCT_NULL_DECOUPLED` if (4) passes and `c` is zero-compatible;
4. `REFINED_H4_NULL_COUPLING_FORMULA_REFUTED` if `c` is nonzero but (5) fails;
5. `REFINED_H4_NULL_COUPLING_SCHEDULE_DEPENDENT` if (5) passes separately but
   the 24 rows or time-reversal comparisons disagree;
6. `REFINED_H4_NULL_COUPLING_COMPATIBILITY_CONFIRMED` otherwise.

Outcome 6 means the lapse line is internally flat but not a full-Hessian
gauge direction; it supplies one schedule-independent linear boundary
compatibility condition.  It does not yet produce a constrained effective
Hessian, propagation, a tick or `c`.

## 7. Execution and acceptance boundary

Register before the first execution.  Run twice with byte-identical output,
then perform only the static registry audit.  No physically weighted outcome
is accepted before a mechanically independent boundary-hinge reconstruction.
Do not run the full suite.

