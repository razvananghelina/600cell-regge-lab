# Final trace-step correction after the disclosed 39/41 run

Date: 2026-08-12

Resource correction commit: `4a28c41`

The resource-safe run completed as `39/41`.  Both failures were the new direct
trace controls.  Their real traces agreed with the gradient-derived Hessian
at relative errors `4.864e-5` and `4.678e-5`, inside the frozen `2e-4`
tolerance.  Their imaginary residuals were `1.21e-2` and `1.01e-2`, outside
the frozen `2e-4` tolerance.  This is consistent with cancellation of three
order-one complex actions divided by `(2e-5)^2`; it is not accepted as a
passing independent trace.

Freeze the following change before evaluating it:

1. restore the originally proposed direct-trace step `h=5e-4`;
2. use that same `h` for the gradient-derived Hessian and the direct-action
   trace, so both reuse the same 70 perturbations and no extra evaluations;
3. retain rank thresholds `1e-7,1e-9,1e-11` and all parity separation bounds;
4. require Hessian symmetry below `2e-5`, direct real trace agreement below
   `2e-5`, and direct trace imaginary residual below `2e-4`;
5. compare the new singular spectra with the already disclosed `h=2e-5`
   spectra; require relative change below `2e-4` for each parity;
6. retain the `38/39` and `39/41` failed-run provenance in the result.

The decisive parity separation remains valid only if it again exceeds `1e-3`.
No root search is introduced.
