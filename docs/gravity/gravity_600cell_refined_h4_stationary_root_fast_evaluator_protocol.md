# Preregistered performance correction: binary64 search evaluator

Date: 2026-08-20

Performance-interruption commit: `08a6c2e`.

The literal 50-decimal SciPy `3-point` search was interrupted after about 20
CPU-minutes without completing its first six-seed class. No scientific
endpoint was recorded. This correction changes only how the already frozen
search-phase residual is evaluated.

## Exact fast evaluator

Implement the same reduced action from the committed integer simplex and
triangle signatures using NumPy complex128 arithmetic:

- construct each `4x4` simplex Gram matrix from the same signed squared edges;
- use direct determinants and inverses for simplex/facet/hinge signed volume
  squares;
- use principal complex square roots;
- implement the committed `log_minus` convention explicitly, including
  `log(-x)-i*pi` on a resolved negative real argument;
- form the same complex dihedral angles, boundary/internal curvatures,
  triangle areas, Schlaefli log gradients and four dust terms;
- use the exact high-precision geometry rounded once to binary64.

No approximate geometric carrier, rounded incidence multiplicity or alternate
equation is allowed.

## Cross-evaluator gate before any solve

For all 12 class representatives compare the binary64 and 80-decimal action,
ten gradients and branch diagnostics at

```text
0, A1', A2', 0.25 p_sigma, 0.5 p_sigma, p_sigma,
```

where `A1',A2'` and `p_sigma` are already committed. Every point must be
branch-valid under the corrected physical-reality convention. Require

```text
max |S_fast-S_mp| / max(1,|S_mp|) < 5e-9,
max_i |G_fast_i-G_mp_i| / max(1,max_i|G_mp_i|) < 5e-9,
|minimum_argument_fast-minimum_argument_mp| < 5e-10,
maximum binary angle-identity residual < 5e-11,
maximum binary physical imaginary part < 5e-9.
```

Time-reversed fast rows must agree below `5e-9`. If any control fails, stop
before the search. Report all 72 cross-evaluator anchors.

## Frozen unchanged search

Use the fast evaluator inside the already preregistered SciPy call with the
same `method='trf'`, `jac='3-point'`, `diff_step=1e-5`, tolerances,
`max_nfev=1200`, six seeds, five boxes, 12 classes and 120-attempt denominator.

The endpoint eligibility gate remains `norm(H(0)^-1 G_fast)<1e-7`. Every
eligible endpoint still undergoes the original 100-decimal damped-Newton
refinement and 140-decimal independent action-derivative validation. Thus no
binary64 endpoint can become a scientific root by itself.

Add per-attempt elapsed time and checkpoint the deterministic JSON after every
completed class. The final artifact replaces the checkpoints only after the
frozen outcome is assigned.

No bound, seed, outcome, physical criterion or look-elsewhere count changes.
