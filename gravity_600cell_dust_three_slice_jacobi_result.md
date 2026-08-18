# Result: a regular action-derived three-slice Jacobi operator

Date: 2026-08-18

## Provenance

```text
prior-art gate                 aca7971
blind protocol                 a3b50fc
registered implementation     a92ab24
blind passing artifact         6bf9392
```

Targeted verifier:
`reproducible/verify_gravity_600cell_dust_three_slice_jacobi.py`.

Artifacts:

```text
reproducible/gravity_600cell_dust_three_slice_jacobi.json
SHA-256 514e01937d621e82c240ea5cad621fb2bc699d09c4940b9be46fa1498152d90c

reproducible/gravity_600cell_dust_three_slice_jacobi.npz
SHA-256 63d95e79c11b25cada660f9a2422654eb92180263dad64e1cbf0ecc30b67d7f8
```

The first targeted execution completed in `109.76 s` with

```text
8/8 PASS
THREE_SLICE_JACOBI_CERTIFIED.
```

A second complete targeted execution took `139.49 s` and reproduced the same
`8/8` outcome and byte-identical JSON/NPZ SHA-256 values.

No spatial spectrum, desired degeneracy, continuum mode or limiting speed was
loaded.  The full suite was not run.

## What was constructed

The already certified one-step maps were partitioned as

```text
T_i = [ A_i  B_i ]
      [ C_i  D_i ]
```

on all `720+720` boundary position/momentum variables.  The committed Flint
balls were re-enclosed after binary serialization.  In every minimal symmetry
sector, the verifier reconstructed the quadratic Hamilton principal function
of each slab and formed the linearized seam equation

```text
K_- delta q_0 + K_0 delta q_1 + K_+ delta q_2 = 0,
```

where each full `q_n` has `720` logarithmic signed-squared spatial edge
variables.  No internal edge remains in this equation.

The normalized recurrence is

```text
delta q_2 = P delta q_1 + Q delta q_0.
```

**DERIVED COMPUTATIONAL:** this is the complete position-space linearized
discrete Euler--Lagrange equation along the first two dynamically solved dust
slabs.  It is equivalent to, but more directly wave-equation-like than, the
canonical phase-space cocycle.

## The boundary twist is safely regular

All

```text
2 schedules * 7 sectors * 4 derivative variants * 2 slabs = 112
```

Flint determinant balls of `B_i` exclude zero.  The operational midpoint
singular data span

```text
minimum singular value        0.00180250 ... 0.00184120,
maximum singular value        0.0152371  ... 0.0549024,
condition number              8.2757     ... 29.8188,
singular uncertainty          2.22e-15   ... 2.23e-15.
```

The weakest minimum singular value exceeds its calibrated uncertainty by
more than `8.08e11`.

**DERIVED COMPUTATIONAL:** eliminating canonical momentum to obtain a
three-position recurrence is not marginal or singular.  This also localizes
the earlier enormous tangent condition numbers: they are not caused by a
failure of the boundary discrete Legendre twist.

## Variational identities

For each slab the reconstructed blocks satisfy, entrywise in Flint balls,

```text
S_00 = S_00*,
S_11 = S_11*,
S_10 = S_01*,
```

and reconstruct the original `A,C,D` tangent blocks.  The three-slice equation
then reproduces both upper blocks of the independently committed product
`T_2 T_1`, in both implicit and solved forms.

Across all residuals, the largest midpoint Frobenius norm is `3.01e-10` and
the largest propagated radius envelope is `5.83e-5`; every entry contains
zero.  The enclosure is conservative rather than marginal: the maximum
midpoint-to-radius diagnostic is below `8.49e-4`.

**DERIVED COMPUTATIONAL:** the momentum signs, boundary order, canonical seam
and Hessian reconstruction are mutually consistent.  The recurrence is not
an independently fitted matrix that merely happens to share a spectrum.

## Schedule robustness

All fourteen target-free comparisons pass:

```text
7 sectors * ([K_- K_0 K_+] and [P Q]) = 14 SCHEDULE_ROBUST.
```

The even/odd ordered-singular-spectrum distances are

```text
K families: 9.09e-13 ... 4.32e-12,
P,Q families: 2.66e-15 ... 4.71e-14,
```

versus calibrated errors from `2.72e-7` to `6.53e-6` and `1.64e-8` to
`5.71e-7`, respectively.

This is robustness of invariant singular data under the two independently
derived staircase schedules, not an assertion that their displayed block
matrices are entrywise identical.

## The first exact wave-operator form

After the blind artifact was committed, the three coefficients expose the
identity

```text
M = (K_- + K_+)/2,
N = (K_+ - K_-)/2,
V =  K_- + K_0 + K_+,

M (delta q_2 - 2 delta q_1 + delta q_0)
+ N (delta q_2 - delta q_0)
+ V delta q_1 = 0.
```

This decomposition is algebraic and unique for the centered three-point
stencil.  It uses no fitted coefficient:

- `M` multiplies the centered second difference;
- `N` multiplies the centered first difference and records the changing
  background;
- `V` is the response to a perturbation held constant across the three
  identified slices.

The already blind background-asymmetry diagnostic

```text
||K_+ - K_-*|| / max(||K_+||,||K_-||)
```

lies between `0.22098` and `0.23738` in every sector.  Hence the two slabs are
not close enough to a time-translation-invariant stencil to discard the
first-difference term.

**DERIVED ALGEBRAIC:** the theory now supplies an exact finite equation of the
form “acceleration + background drift + stiffness” for anisotropic edge
perturbations.

**OPEN:** `V` has not yet been shown to be a spatial Laplacian or the tensor
operator of continuum gravity.  `M` has not yet been shown to be a positive
physical kinetic metric.  Therefore this is not yet a derived d'Alembertian,
dispersion relation or value of `c`.

## Framing correction

The existence of the recurrence is mathematically necessary for any regular
discrete variational system.  It is not, by itself, evidence that the model
describes gravity.  What is project-specific and nontrivial is that:

1. the complete `720`-position twist survives all interval and schedule
   controls;
2. it is exactly compatible with the independently constructed Regge
   phase-space cocycle;
3. all internal variables have been eliminated without fitting;
4. the resulting operator is now suitable for a falsifiable spatial-mode
   comparison.

The large `120+120` phase-space branch count must not be read as a wave count
before the kinetic and stiffness operators are classified.

## Post-result prior-art reconciliation

- [Marsden--West](https://doi.org/10.1017/S096249290100006X) establish the
  general equivalence of regular discrete Legendre maps and discrete
  Euler--Lagrange recurrences.
- [Dittrich--Hoehn](https://arxiv.org/abs/1108.1974) provide the corresponding
  action-generated canonical formalism for simplicial gravity.
- [Dittrich--Hoehn](https://arxiv.org/abs/0912.1817) explain the
  background-dependent linearized/pseudo-constraint structure.
- [Hoehn](https://arxiv.org/abs/1411.5672) shows why curvature and gauge
  classification are required before calling Regge perturbations gravitons.
- [Rostworowski](https://arxiv.org/abs/1902.05090) is a continuum control:
  FLRW metric perturbations reduce, after constraint handling, to wave-type
  master equations plus a matter transport mode.  It does not establish the
  same reduction on this finite carrier.

No located primary source prints this exact full dust 600-cell Jacobi
operator.  External novelty remains **OPEN**; a search is not proof.

## Next falsifiable gate

Preregister the centered operators `M,N,V` before reading their spectra.
The next gate must decide:

1. whether `M` is regular and admits a defensible positive kinetic form;
2. whether `N` is quantitatively resolved and schedule robust rather than a
   coordinate artefact;
3. whether the generalized stiffness `M^-1 V` has real, symmetry-stable
   spectral sectors;
4. only after committing that blind spectrum, whether it intertwines with an
   independently derived intrinsic 600-cell operator.

Failure of regularity, positivity or real generalized stiffness would block
the proposed wave interpretation immediately.  A later match of spectra
without an intertwiner would remain a **PATTERN**, not a derivation of `c`.
