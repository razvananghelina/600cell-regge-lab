# Protocol: does the continuum Kähler--Dirac heat trace see the Hopf spin-two carrier?

Date: 2026-08-12

## Provenance

This is a post-recognition protocol.  Before registration it was already
known from the standard heat-kernel formula that the ordinary (not
supertrace) de Rham heat trace can contain a scalar-curvature term.  The
calculation below is therefore a hostile derivation and scope audit, not a
blind discovery.  No measured gravitational constant, mass, speed or
four-dimensional target will be inspected.

## Complete hypotheses

1. The continuum carrier is the unit round three-sphere
   `S3=SU(2)`, with a left-invariant orthonormal frame satisfying

   ```text
   [e_a,e_b]=2 epsilon_abc e_c.
   ```

2. The perturbation space is exactly the five-dimensional homogeneous Hopf
   carrier already certified in `hopf_spin2_tensor_carrier_result.md`:
   constant-frame symmetric tracefree tensors

   ```text
   H in Sym^2_0(R^3).
   ```

   The left-handed space is tested explicitly.  The right-handed result may
   follow only from inversion symmetry of the round metric; it is not to be
   counted as ten independent local fibre components.

3. A positive-definite left-invariant metric is represented by its Gram
   matrix `G` in the fixed Lie-algebra frame.  Its Levi--Civita connection,
   curvature and scalar curvature are to be derived directly from the Koszul
   formula.  No Einstein equation or Regge action is inserted.

4. Remove the pure volume direction with the scale-invariant, fixed-volume
   Einstein functional

   ```text
   Y(G)=R(G) det(G)^(1/3).
   ```

   This is used only to read the universal scalar-curvature response of the
   heat coefficient.  Calling `Y` a theory-selected physical action is
   forbidden.

5. The spectral operator is the continuum de Rham Kähler--Dirac operator

   ```text
   D_g=d+d*_g
   ```

   on the complete exterior algebra in three dimensions.  For its ordinary
   heat trace, derive the local `a2` coefficient by summing the
   Weitzenböck curvature traces over degrees `p=0,1,2,3`.  The finite moments
   `(2640,14880,55920)` are not Seeley--DeWitt coefficients and may not be
   used as substitutes.

6. The current repository operators retain their recorded scopes:

   - `verify_kahler_dirac.py` uses the incidence inner product and has no
     continuous metric variable;
   - the Whitney construction has exact metric-dependent mass matrices but
     presently fixes the piecewise-flat Regge metric;
   - no cutoff function, heat scale, Newton normalization or variable-metric
     spectral action is selected.

7. No transfer from the round metric to the fixed Regge metric, no
   Lorentzian time, no diffeomorphism quotient, no stress-tensor source and no
   Planck scale are assumed.

## Frozen calculations

1. Derive `Gamma`, `Riemann`, `Ricci` and `R(G)` from the structure constants
   and a generic symmetric `3 x 3` Gram matrix.
2. Verify the round normalization `R(I)=6`.
3. Compute the exact Hessian of `Y(G)` at `G=I`, restrict it to an explicit
   orthogonal basis of `Sym^2_0(R^3)`, and report its rank and inertia.  Check
   separately that the scale direction is null.
4. Derive degree by degree

   ```text
   tr(E_p-curvature contribution)
   ```

   for the Hodge Laplacian and hence the scalar-curvature coefficient in

   ```text
   Tr exp(-t D_g^2).
   ```

5. Combine the two exact results.  Along determinant-one Hopf metric paths,
   the volume coefficient is constant; determine whether the first
   curvature coefficient has a nondegenerate quadratic response on all five
   tensor directions.
6. Audit the authoritative constructors for the three absences in Hypothesis
   6.  This is a coverage statement, not evidence for a general no-go beyond
   the repository.

## Decision boundary

- **DERIVED CONTINUUM KINEMATICS:** the exact curvature and heat-coefficient
  identities under the complete smooth round hypotheses.
- **STRUCTURAL INDUCED-GRAVITY ADVANCE:** the theory's own continuum operator
  has a nonzero, nondegenerate scalar-curvature response on the Hopf carrier,
  with no relative tensor coefficient fitted.
- **DERIVED NEGATIVE:** the response vanishes or has a kernel on the
  five-dimensional carrier.
- **NO PHYSICAL PROMOTION:** even a positive result is not emergent gravity
  unless the finite/refined theory selects the variable metric, spectral
  functional and normalization and later supplies Lorentzian/gauge/source
  gates.

The result must state plainly whether the computation advances only a
continuum compatibility mechanism or an already selected sector of the
finite theory.
