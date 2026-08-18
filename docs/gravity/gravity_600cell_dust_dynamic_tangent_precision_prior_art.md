# Prior-art gate: precision audit of the dynamic tangent spectrum

Date: 2026-08-17

## Question

Can the `6.86e-11` matched-eigenvalue separation between the two blind dynamic
shape maps be resolved after accounting for nonnormal eigenvalue sensitivity,
or was the frozen binary64 classifier overconfident?

This is a numerical correction audit of the already committed blind artifact,
not a new physical spectrum search.  No continuum eigenvalue, desired
degeneracy, speed, or experimental number is a target.

## Primary-source check

1. F. L. Bauer and C. T. Fike, *Norms and exclusion theorems*, Numerische
   Mathematik 2 (1960), 137--141,
   <https://eudml.org/doc/131452>.  For a diagonalizable nonnormal matrix, a
   matrix perturbation can move eigenvalues by the perturbation norm multiplied
   by an eigenvector-basis condition factor.  Therefore a raw eigensolver floor
   proportional only to machine precision and spectral radius is not a complete
   error model here.
2. L. N. Trefethen, *Pseudospectra of Linear Operators*, SIAM Review 39 (1997),
   <https://doi.org/10.1137/S0036144595295284>.  Strongly nonnormal spectra can
   be much more sensitive than norms or finite-time propagation, so eigenvalues
   alone must not be given a stability interpretation without a sensitivity
   audit.
3. LAPACK Users' Guide, *Error Bounds for the Nonsymmetric Eigenproblem*,
   <https://www.netlib.org/lapack/lug/> and the `xGEEVX` test documentation,
   <https://www.netlib.org/lapack/explore-html/de/d5d/cget23_8f_abe009f4d30645a5bc7c97ab01712081d.html>.
   The standard nonsymmetric workflow explicitly distinguishes eigenvalue and
   eigenvector condition estimates and validates left/right residuals.

No source found in this gate computes the present order-24 600-cell dust
tangent maps.  The relevant prior art supplies the error analysis, not the
answer.

## Consequence for this project

The blind verifier used `numpy.linalg.eig`, then estimated the condition of the
whole eigenvector matrix, but its frozen `epsilon_eig` did not propagate that
condition factor.  The shape estimates (`9.15e4` and `5.64e5`) make the omitted
effect potentially comparable with or larger than the reported schedule
separation.

The correction must preserve the original artifact and verdict.  It will use
two independent checks:

- high-precision, optimally matched eigenvalues with residuals and a
  Bauer--Fike-style calibrated uncertainty;
- normalized power traces, which are basis-independent spectral invariants and
  avoid eigenvector conditioning.

The precise invariant set, normalization and mechanical verdicts must be
committed before evaluating them.

## Status

- **DERIVED:** the frozen binary64 schedule classifier omitted nonnormal
  sensitivity from its uncertainty.
- **STRUCTURAL:** Bauer--Fike is an upper-bound framework, not evidence that the
  two spectra are equal.
- **OPEN:** whether the two committed maps have a calibration-resolved spectral
  difference.
