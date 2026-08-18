# Prior-art gate: shifted principal-function precision

Date: 2026-08-18

## Exact object

The precision-open persistence result traces its shifted centered coefficients
through

```text
high-precision slab Hessian
  -> canonical tangent T_2
  -> binary64 midpoint/radius archive
  -> re-enclosed tangent blocks [A B; C D]
  -> principal-function block S_2,10 = C-D B^-1 A
  -> shifted K_- = S_2,10
  -> centered M,V
  -> restricted shape stiffness.
```

The proposed audit does not recompute an eigenvalue.  It asks whether the
mandatory half-ULP enclosure of the serialized binary64 tangent, after the
ill-conditioned principal-function reconstruction, dominates the final
restricted-form error.  A direct high-precision reconstruction is admissible
only if this audit selects serialization as the dominant source.

## Primary literature checked after the open result

1. Squire and Trapp, [*Using Complex Variables to Estimate Derivatives of
   Real Functions*](https://doi.org/10.1137/S003614459631241X), establish the
   cancellation-avoiding complex-step principle used by the local derivative
   variants.  It supports recomputing derivatives at high precision; it does
   not justify discarding serialization error after conversion to binary64.
2. Martins, Sturdza and Alonso, [*The Complex-Step Derivative
   Approximation*](https://doi.org/10.1145/838250.838251), develop accurate
   and robust implementations and connect complex-step differentiation with
   algorithmic differentiation.  It is a numerical-method control, not a
   result about Regge spectra.
3. Fike and Alonso, [*The Development of Hyper-Dual Numbers for Exact
   Second-Derivative Calculations*](https://doi.org/10.2514/6.2011-886), give
   an alternative route to second derivatives without a differencing step.
   Replacing the current derivative construction with hyper-duals would be a
   new implementation and is not silently introduced here.
4. Dittrich, Freidel and Speziale, [*Linearized dynamics from the 4-simplex
   Regge action*](https://arxiv.org/abs/0707.4513), demonstrate that the Regge
   Hessian is the appropriate linearized object.  They do not address the
   binary serialization of a full 600-cell canonical tangent.

## Ledger

- **KNOWN:** a binary64 midpoint differs from the unrounded high-precision
  value by at most half an ULP per real component; a rigorous re-enclosure
  must retain that uncertainty.  Solving through `B^-1` may amplify it.
- **CONTROL:** the stored Flint radii of the `T_2` tangent blocks are about
  `1e-61`, while the downstream shifted `K_-` radii are about `1e-6`.  This
  scale comparison was observed after the persistence result and is not yet a
  certified attribution.
- **OPEN:** the separate half-ULP, stored-ball, inversion and downstream
  carrier contributions; whether direct reconstruction from the original
  high-precision Hessian reduces the complete error enough to decide a sign.
- **FORBIDDEN:** deleting the half-ULP while retaining only the binary
  midpoint; shrinking a radius until a desired eigenvalue becomes negative;
  changing derivative steps in response to the target sign.
- **OPEN external novelty:** no source found performs this exact precision
  audit for the three accepted 600-cell dust-Regge slabs.  The search is not
  proof of novelty.

The audit is numerical hygiene.  Even a successful precision recovery would
certify only a finite action-selected stiffness carrier, not a graviton or a
physical instability.
