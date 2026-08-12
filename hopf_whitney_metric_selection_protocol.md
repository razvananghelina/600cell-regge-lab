# Protocol: does the geometry select the round or fixed-Regge Whitney metric?

Date: 2026-08-12

## Provenance

This is a post-recognition hostile audit.  Before registration it was known
that the repository uses the round metric as an auxiliary norm while keeping
the pushed-forward piecewise-flat metric as the exact Whitney target.  It was
also known that the new Hopf curvature response was derived at the round
metric.  The unresolved question frozen here is whether the already accepted
geometric criteria distinguish those metrics.  This is not a blind search.

No heat spectrum, gravitational target or phenomenological number is to be
inspected before the metric-family census is complete.

## Complete carrier and hypotheses

1. `P` is the regular unit-circumradius 600-cell, `K=boundary(P)`, and

   ```text
   R:K -> S3,   R(x)=x/|x|
   ```

   is the radial identification.

2. On every tetrahedral facet there are two already present metrics on the
   same transported Whitney spaces:

   ```text
   g_R=(R^-1)^* g_flat,
   g_0=the unit round metric on S3.
   ```

   Equivalently on the flat facet, `R^*g_R=g_flat`, while

   ```text
   (R^*g_0)_x(v,w)
      =(v.w)/r^2-(x.v)(x.w)/r^4.
   ```

3. The refinement carrier is the already selected rank-edgewise tower

   ```text
   K_n=Esd_(2^n)(sd K).
   ```

   The differential and Whitney degrees of freedom are metric-independent.
   Inner products are exact `L2(g)` integrals; mass lumping and quadrature
   approximations are forbidden in the canonicity conclusion.

4. The admissible selection criteria being audited are exactly:

   - equivariance under every orthogonal symmetry of the 600-cell;
   - construction from the existing radial map and ambient metrics, with no
     vertex labels;
   - positivity and tangential compatibility across shared faces;
   - use of the same Whitney subcomplex and exact refinement inclusions;
   - level-independent uniform equivalence to the existing metric.

5. No new action, curvature minimization, cutoff function, field equation or
   measured target may be used to choose a metric.  Such input could be a
   future selector but is outside this audit.

## Frozen tests

1. Re-derive the exact facet distance

   ```text
   a^2=(7+3 sqrt(5))/16
   ```

   and the pullback-round eigenvalues relative to the flat metric:

   ```text
   a^2/r^4, 1/r^2, 1/r^2.
   ```

2. At a fixed interior point on a facet, prove exactly that `g_R` and `g_0`
   are not proportional.  Endpoint or floating-point evidence does not
   count.

3. Prove that both metrics are equivariant under the full orthogonal
   symmetry group and have matching tangential restrictions on shared
   faces.

4. Enumerate, before any action comparison, the affine family

   ```text
   g_u=(1-u)g_R+u g_0,   0<=u<=1.
   ```

   Prove positivity, injectivity in `u`, symmetry equivariance and the same
   uniform equivalence bounds.  This family is evidence of insufficiency of
   the audited criteria, not a claim that arbitrary `u` is selected.

5. On an independently built tetrahedral refinement, use a generic positive
   diagonal metric `diag(x,y,z)` and exact Whitney integration to verify in
   every degree

   ```text
   P_p^T M_f,p(g) P_p=M_c,p(g).
   ```

   This must establish that refinement isometry does not distinguish a
   metric; checking only one Euclidean metric is insufficient.

6. Verify that the two endpoint local geometries are physically different:
   `g_0` has scalar curvature six, while `g_R` is flat in every facet
   interior and carries its curvature on the Regge skeleton.  Norm
   equivalence must not be used to transfer a heat coefficient or Hessian.

7. Audit the accepted Whitney result and the new Hopf result for their stated
   metric scopes.

## Decision boundary

- **SELECTED ROUND TRANSFER:** only `g_0` passes all frozen geometric and
  refinement criteria.
- **SELECTED REGGE TRANSFER:** only `g_R` passes.
- **DERIVED METRIC-SELECTION NO-GO:** both endpoints, or a nontrivial
  continuum between them, pass all frozen criteria.  Then the round Hopf
  heat Hessian cannot be assigned to the exact fixed-Regge Whitney theory by
  symmetry/refinement alone.
- **STRUCTURAL OPENING:** even under the no-go, a canonical Hopf tangent
  tensor may exist over each chosen baseline.  That does not select the
  baseline or action.

The result must say whether this closes emergent gravity or only the present
metric-selection argument.  A future independently derived action is not
excluded by failure of these kinematic criteria.
