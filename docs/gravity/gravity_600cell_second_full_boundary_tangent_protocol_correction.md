# Pre-evaluation correction: repository two-step precedent

Date: 2026-08-22

Original prior-art commit: `eae691a`.

Original protocol commit: `f80e174`.

Status: **CORRECTED BEFORE REGISTERING OR IMPLEMENTING THE NEW VERIFIER AND
BEFORE CONSTRUCTING ANY NEW SECOND-SLAB HESSIAN OR TANGENT.**

## Omission

The original prior-art gate searched the primary literature and the current
finite-height tangent chain, but failed to locate an older repository route
whose filename already states the relevant object:

```text
reproducible/verify_gravity_600cell_dust_two_step_full_tangent.py
docs/gravity/gravity_600cell_dust_two_step_full_tangent_result.md
```

This is a genuine prior-art-gate failure.  The new calculation had not yet
started, so no numerical result or threshold was available when the omission
was found.

## Frozen additional inputs

```text
docs/gravity/gravity_600cell_dust_two_step_full_tangent_prior_art.md
  e7d865b9e72a411eee61e0ce091cde0d912fd9d9f773708f00a9c7046a6785f9

docs/gravity/gravity_600cell_dust_two_step_full_tangent_protocol.md
  d5dd44ece724b65351b35fd18e6d334dbf4b68e9f2757484b76fdaf6c42fe0cf

docs/gravity/gravity_600cell_dust_two_step_full_tangent_result.md
  014a86460433e9e8ab72a2aae029bed774306c2095aef1a543522d780c783038

reproducible/verify_gravity_600cell_dust_two_step_full_tangent.py
  c1a3fb09146188c1932ab81629ab69817f2a2f19108fdf8d9e89d78b6de8f717

reproducible/gravity_600cell_dust_two_step_full_tangent.json
  f7fbf18535cc00dacec9a9ffa95f97f2d1847ac83073f27d39fcdb7968b0bafc

reproducible/gravity_600cell_dust_two_step_full_tangent.npz
  ce78ebf415584b1cdcf1d2cb07687135b624ad4939e0a4e54650653f7b384e6d
```

Require the old artifact to retain `16/16` and
`TWO_STEP_FULL_TANGENT_COCYCLE_CERTIFIED`.  It is a method and hostile-framing
control, not a numerical target for the new background.

## What was already done

The older route already constructed complete `1440`-dimensional first- and
second-slab canonical tangents, multiplied them, checked Flint-ball
symplectic identities and compared schedule-invariant spectra.  A later route
also reconstructed the equivalent three-slice Jacobi equation.

Therefore:

- composable full 600-cell dust tangents are **REPOSITORY-KNOWN**;
- the action-to-Jacobi machinery is **REPOSITORY-KNOWN**;
- the method itself cannot be claimed as a new discovery;
- the phrase "for the first time two perturbative ticks" is withdrawn.

## Why the new calculation is not identical

The old and new missions have different frozen backgrounds and different
acceptance gates.

1. The old route used the published/nearly-static chain with inherited
   `tau=0.0102` and an older fixed scale.  The new route uses the exact
   finite-height branch reconstructed from `v=3/2`, with
   `h1 approximately 0.20405`, `r1 approximately 2.96260`, followed by the
   accepted branch-B slab `h2 approximately 0.06893`,
   `r2 approximately 3.15620`.
2. The old route inherited an earlier full tangent.  The new route starts
   from the 180/140-decimal finite-height tangent that was replicated by a
   complete dense real-space construction.
3. The old route assembled its unequal-scale second slab directly in
   physical units but did not separately compare it with a normalized
   assembly through the exact degree-two scale lift.  That comparison is
   load-bearing here.
4. The old route compared even-even and odd-odd spectral invariants.  The new
   route compares direct maps in literally common bases and all four
   `(first schedule, second schedule)` products.
5. The old route computed eigenvalue and singular-value censuses.  The new
   blind mission explicitly forbids them and classifies only regularity,
   canonicality, physical-unit conjugacy and direct schedule robustness.

Thus a positive result would be a **new-background replication and stronger
unit-consistency control**, not a new two-step formalism.  A negative result
would be more informative: it would show that the older near-static cocycle
does not survive on the newly selected finite-height branch.

## Protocol amendment

All clauses of `gravity_600cell_second_full_boundary_tangent_protocol.md`
remain frozen, with these additions:

1. the six hashes above join the mandatory provenance set;
2. the new verifier may reuse audited helper definitions from the old source,
   but must not execute its bottom-level scientific calculation or parse its
   spectra before the new labels are frozen;
3. the old tangent, product, eigenvalues and singular values are forbidden as
   numerical targets;
4. the final interpretation must say explicitly that two-step tangent and
   Jacobi machinery were already present on a different background;
5. external novelty of the method is **REFUTED**; external novelty of the
   coefficient-level new-background result remains **OPEN**.

No derivative step, bracket, matrix gate, 10/100 classifier, hostile control
or outcome hierarchy is changed.

