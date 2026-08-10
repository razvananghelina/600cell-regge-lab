# Preregistration: local-tensor extension to barycentric modes

Date: 2026-08-10

## Definition fixed before refined spectra

The STEP-1 blind base enumeration found one and only one Hopf fiber edge in
every coarse tetrahedron, for all six fibrations.  The following extension is
fixed before inspecting any first-refinement gap ratio.

For a parent tetrahedron `T`, let `e_T` be its unique fiber-edge vector, `V_T`
its Euclidean volume and `w_T` the equal-edge coefficient in its local scalar
Whitney stiffness.  Define the positive tangent tensors

`Q_f(T)=(w_T/V_T) e_T tensor e_T`,

`Q_c(T)=P_tangent(T)-Q_f(T)`.

For every one of the 24 affine barycentric children `t` of `T`, define

`K_f(t)_(ij)=integral_t grad(lambda_i).Q_f(T).grad(lambda_j)`,

and analogously for `K_c`.  No new fine edge is labelled fiber or cross.

This rule is geometrically local, positive, contains no level coefficient and
restricts the same parent tensor to all children.  It also has a hard audit:
the consistent Whitney inclusion must compress the fine mass, fiber form and
cross form back to their coarse counterparts.

## Blind output

`reproducible/verify_whitney_hopf_refinement_blind.py` will write:

- tensor ranks/eigenvalues;
- all three compression residuals;
- kernel multiplicities reached before the first positive mode;
- low generalized spectra and first positive gaps at 120 and 2640 scalar
  degrees of freedom;
- the two observed cross/fiber ratios.

The verifier contains no bootstrap integer or proposed physical speed.  Its
definition is committed before execution and target comparison.

## Acceptance and kill criteria reserved for the later comparison

- **Acceptance candidate:** the refined dynamical ratio is stable without a
  new coefficient and the relevant low-mode sector remains identifiable.
- **Kill candidate:** new refinement modes change the first-gap ratio, or the
  fiber form develops an unresolved/extensive kernel which makes a single
  propagation speed ill-defined.

No verdict is made in this preregistration file.
