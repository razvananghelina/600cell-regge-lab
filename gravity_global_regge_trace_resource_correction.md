# Resource correction: reuse the preregistered R0 action perturbations

Date: 2026-08-12

Parity correction commit: `08b638c`

The first corrected run retained every original passing check and both
Hessians, then the process terminated before producing the newly requested
direct Hessian trace.  Therefore no direct-trace value at the proposed
`5e-4` step was observed.

The script already evaluates the complete reduced action at

```text
R0 +/- 2e-5 e_j,   j=0,...,34,
```

as part of the preregistered direct-gradient/Hessian control.  Freeze the
following resource-preserving replacement before inspecting its result:

1. retain every plus and minus action from those 70 evaluations;
2. compute

   ```text
   sum_j [S(x+h e_j)-2S(x)+S(x-h e_j)]/h^2,
   h=2e-5;
   ```

3. compare it with the trace of the independently gradient-differentiated
   Hessian;
4. require relative agreement `2e-4` and imaginary residual `2e-4`;
5. perform no additional action evaluation for this trace.

The larger tolerance accounts for subtracting three order-one action values
at `h^2=4e-10`.  This is still much smaller than the disclosed `1.62e-2`
parity spectral separation.  All other thresholds and the decisive parity
gap `(equal <2e-5, separated >1e-3)` remain unchanged.

This is an implementation/resource correction, not a response to a measured
trace: the attempted `5e-4` calculation produced none.
