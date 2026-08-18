# High-precision audit of the dynamic tangent spectrum

Date: 2026-08-17

## Verdict

> **DERIVED COMPUTATIONAL:** the previously reported schedule dependence was a
> binary64 eigensolver false positive.  At 160 decimal digits, both the complete
> 58-member power-trace family and the optimally matched eigenvalue multisets are
> consistent between the two schedules by margins enormously larger than the
> frozen calibration requires.

Frozen outcome:

```text
SCHEDULE_SPECTRUM_NOT_RESOLVED
```

This supersedes the scientific interpretation of the blind artifact's
mechanical label `DYNAMIC_SHAPE_TANGENT_SCHEDULE_DEPENDENT`; it does not edit or
delete that preregistered result.

## Provenance

- blind dynamic census commit: `da5e5ce`
- precision prior-art gate: `16f4310`
- precision protocol, before comparison: `0e8f561`
- registered verifier implementation, before first execution: `fa401fb`

Verifier:
`reproducible/verify_gravity_600cell_dust_dynamic_tangent_precision.py`

Artifact:
`reproducible/gravity_600cell_dust_dynamic_tangent_precision.json`

Artifact SHA-256:
`c460f2573b32eb0cf05a5905839dd755e7718f78e5ff11a3a060f9e83f229104`

The targeted verifier passed `8/8`.  The full suite was not run, following the
explicit instruction for the current work.

## Numerical result

The original binary64 comparison reported a maximum matched eigenvalue
distance

```text
6.8628065e-11.
```

The independent 160-digit calculation instead gives

```text
maximum matched distance       = 4.107115179e-41
calibrated eigen uncertainty   = 6.474596258e-12
distance / uncertainty         = 6.343430564e-30.
```

The normalized eigensystem residuals are `9.75e-162` and `1.05e-161`.  The
Frobenius eigenvector-basis condition factors are approximately `1.920e5` and
`1.912e5`; propagating them is exactly what the blind classifier omitted.

Independently, all power sums

```text
Tr((T_shape/2^20)^k), k=1,...,58,
```

pass the frozen consistency threshold.  The largest difference is at `k=1`:

```text
difference                    = 1.462300618e-56
calibrated uncertainty        = 1.869123491e-21
difference / uncertainty      = 7.823456426e-36.
```

Because the first 58 power sums determine the characteristic polynomial of a
58 by 58 matrix, this is a basis-independent check of the entire finite
characteristic spectrum, not a favorable subset of eigenvalues.

## Correct interpretation

- **DERIVED:** the two finite order-24 dynamic shape maps are spectrally
  indistinguishable at the committed calibration.
- **DERIVED:** the earlier `6.86e-11` separation came from finite-precision
  diagonalization of nonnormal matrices, not from the committed high-precision
  maps.
- **PATTERN:** agreement to roughly forty eigenvalue digits and fifty-six trace
  digits strongly suggests a structural isospectrality.
- **OPEN:** exact isospectrality has not been proved.  The inputs are calibrated
  numerical Hessians, not exact algebraic matrices.
- **STRUCTURAL:** spectral agreement does not imply equality of the two maps,
  identical singular vectors, or schedule-independent nonlinear evolution.
- **OPEN:** no graviton identification, continuum limit, refinement stability,
  causal cone, limiting speed, or physical tick follows from this audit.

The shape spectral radius `45.3745` is still not labelled an instability.  The
maps are highly nonnormal (largest singular value about `4.495e5`), and no
physical norm or continuum mode identification has yet been derived.

## Post-result primary-source audit

The high-precision result does not contradict the canonical discrete-gravity
framework of Dittrich and Hoehn,
<https://arxiv.org/abs/0912.1817>: their formalism allows discretization-broken
symmetries and background-dependent pseudo-constraints, so it does not supply a
general schedule-independence theorem for this curved dust slab.

Product-eigenvalue theory treats cyclic reorderings of operator factors as a
natural source of shared spectra; see D. S. Watkins, *Product Eigenvalue
Problems*, <https://doi.org/10.1137/S0036144504443110>.  In the elementary
invertible two-factor case, `AB` and `BA` are similar.  This is a plausible
explanation only if the two geometrically derived tick maps admit the required
factorization or carrier conjugacy.  None has yet been derived here.

Therefore the next honest question is not to fit an arbitrary similarity
matrix--any two diagonalizable isospectral matrices would make that nearly
content-free.  It is to preregister the finite set of conjugacies supplied by
the 600-cell carrier, schedule reversal, boundary permutation and canonical
symplectic structure, then test only those.
