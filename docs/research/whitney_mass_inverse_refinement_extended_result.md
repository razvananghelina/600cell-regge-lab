# Extended refined Whitney inverse lower bounds

Date: 2026-08-11

Extended preregistration commit: `14e187d`

This is explicitly a post-result extension of the shorter-window protocol
`366fe4a`; it is not presented as blind discovery.

Targeted verifier:
`reproducible/verify_whitney_mass_inverse_refinement.py`

Targeted result: **12/12 PASS**.  No full-suite run was performed.

## Result

With 1,024 exact modular Krylov moments at each of three preregistered primes,
the first-refinement minimal-polynomial lower bounds become

\[
\deg m_{\rm fine}\geq(117,512,512,1),
\]

and therefore

\[
\deg M_{\rm fine}^{-1}\geq(116,511,511,0).
\]

The comparison with the exact unrefined inverse degrees is

| form degree | coarse inverse degree | refined inverse lower bound | growth factor lower bound | window status |
|---:|---:|---:|---:|---|
| 0 | 8 | 116 | 14.50 | below ceiling; stabilized candidate |
| 1 | 21 | 511 | 24.33 | censored at ceiling |
| 2 | 26 | 511 | 19.65 | censored at ceiling |
| 3 | 0 | 0 | n/a | exact scalar control |

All three primes returned exactly the same complexities.  The same extended
procedure still returns the already certified coarse degrees `(9,22,27,1)`,
so the longer-window calibration passes.

## Exact versus candidate statements

Every number in `(117,512,512,1)` is a **DERIVED exact lower bound**.  This
requires no genericity assumption: a recurrence for the full rational matrix
would reduce to a recurrence for each finite modular probe sequence.

The interpretation of equality differs by block:

- Degree 0 returns 117 well below the 512 ceiling at all three primes.  This
  is a **STABILIZED CANDIDATE** for its full minimal degree, but equality is
  **OPEN** until a degree-117 whole-matrix annihilator is reconstructed and
  certified.
- Degrees 1 and 2 return exactly 512, the maximum informative complexity of
  the frozen 1,024-term windows.  They are **CENSORED**; their true degrees
  may be much larger.
- Degree 3 is the primitive identity matrix, so degree one is exact.

It would be false to print `(117,512,512,1)` as the refined minimal degrees.

## Consequence

The earlier fixed-complex statement—finite polynomial inverses of degrees
`(8,21,26,0)`—is mathematically correct but physically fragile.  At the very
first exact barycentric refinement, any inverse polynomial in the same mass
operator requires at least `(116,511,511,0)` degrees.

> **DERIVED NEGATIVE:** exact Whitney metric glue cannot use the base
> polynomial depths as a resolution-independent finite internal tick.

> **PATTERN:** inverse complexity grows by more than an order of magnitude at
> one refinement in every propagating form degree.

This is strong evidence against bounded-depth explicit mass inversion, but it
is not yet a divergence theorem.  Only two levels exist, and the upper bounds
for degrees 1 and 2 are unknown.  No scaling exponent is fitted.

The result is also specific to algorithms expressing the inverse as a
polynomial of the local mass.  A different reversible construction with
ancillas, multigrid structure or an enlarged first-order carrier could in
principle bypass this degree.  Such a construction must be explicit and
coefficient-free; its mere abstract existence would have no evidential
weight.

## Status ledger

- **DERIVED:** extended fine minimal-degree lower bounds
  `(117,512,512,1)`.
- **DERIVED:** extended fine inverse-degree lower bounds
  `(116,511,511,0)`.
- **DERIVED NEGATIVE:** resolution-independent base-depth Whitney glue.
- **PATTERN:** severe first-refinement complexity growth.
- **OPEN:** exact fine degree 0 polynomial certificate.
- **OPEN:** uncensored fine degrees 1 and 2.
- **OPEN:** repeated-refinement divergence or a bounded alternative dilation.
- **NOT CLAIMED:** physical infinity, time, inertia, mass or (c).

## Reproduction

```bash
/home/razvan/science/.venv/bin/python \
  reproducible/verify_whitney_mass_inverse_refinement.py
```

Expected result: `12/12`.
