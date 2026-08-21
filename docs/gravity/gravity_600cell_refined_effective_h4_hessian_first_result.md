# First result: the curvature-matched internal Hessian is singular

Date: 2026-08-21

Status: primary computational result, not accepted before null-space analysis
and adversarial replication.

## Provenance and reproducibility

| stage | commit |
|---|---|
| prior-art gate | `f7bf3c1` |
| frozen protocol | `eed7891` |
| verifier registered before execution | `f99b8f4` |

The targeted verifier was run twice.  Both runs returned `16/16` with the
same outcome and byte-identical artifact

```text
reproducible/gravity_600cell_refined_effective_h4_hessian.json
SHA-256 6760ce4dbceb43dab88f3c3052e634f8e5e4d05c69f8fc6fcd02bb3fa6fdc254.
```

No full suite or deferred nonlinear census was run.  The static registry
audit reports `375` registered entries, `375` distinct names, zero duplicates,
zero unregistered files and zero missing files, plus two reasoned exclusions.

## Primary result

The background controls pass for every one of the 24 schedules:

```text
maximum internal residual       5.91e-77
maximum branch identity error   4.44e-139
maximum imaginary curvature     9.39e-139
maximum full-Hessian envelope   2.07e-29
```

Every `10x10` internal Hessian has certified inertia

```text
(positive, zero-compatible, negative)=(9,1,0).
```

For the first schedule, its smallest displayed eigenvalue is
`-8.92e-39`, inside the `2.01e-27` spectral envelope.  The next eigenvalue is
`1.19521e-5`, far outside it.  The same isolated one-dimensional null sector
occurs for all 24 schedules; it is not a rank decision set by a marginal
tolerance.

The frozen outcome is

```text
REFINED_EFFECTIVE_H4_HESSIAN_INTERNAL_SINGULAR.
```

## What this does and does not establish

- **PRIMARY DERIVED COMPUTATIONAL:** the curvature-selected matter branch
  changes the off-shell P1 internal inertia from the earlier full-rank
  `(8,0,2)` result to an on-shell `(9,1,0)` result for every schedule.
- **DERIVED PROCEDURAL STOP:** the ordinary inverse in
  `K=H_bb-H_bi H_ii^(-1)H_ib` is forbidden.  No effective boundary Hessian,
  schedule class, dispersion relation or speed was obtained.
- **OPEN:** whether the one-dimensional null line is the exact product-lapse
  tangent, whether it is a gauge direction of the full Hessian, and whether
  it decouples from all boundary variations.
- **OPEN:** if it couples to the boundary, which compatibility constraint it
  imposes and whether a reduced principal form exists only on that constrained
  boundary subspace.

An arbitrary Moore--Penrose inverse would silently choose a representative
unless the null coupling is proved to vanish.  It is therefore forbidden.

## Reporting issue found after the run

Because the singular gate correctly prevented construction of `K`, three
downstream branches contained zero samples: directional effective-action
checks, time-reversal/class enumeration and the synthetic corruption of an
effective matrix.  The implementation represented those conditional skips
as boolean successes, producing misleading lines such as `classes=0` and a
zero corruption difference under a `PASS` label.

Those lines cannot falsify anything and are **not evidence**.  They do not
enter the scientific outcome, which was already fixed by the singular gate,
but the reporting must be corrected before this verifier is treated as
accepted.  The allowed correction is limited to explicit conditional-skip
labels and assertions that no downstream effective claim was made.

The correction was preregistered in `ebc32a2` and implemented in `1ab4333`.
The corrected verifier again passed `16/16` twice.  It now prints and
serializes the unavailable quantities as
`NOT_COMPUTED_INTERNAL_SINGULAR`/`null`; the corrected byte-identical artifact
has SHA-256

```text
56e08db9a840b95e686fadb2763e89400b09220e88b80e9d35c17c1e73eef0a3.
```

All eigenvalues, envelopes and the frozen singular outcome are unchanged.

## Next exact question

Compute, for all schedules and independently of an arbitrary eigensolver
sign:

1. the overlap of the null line with the analytically induced product-lapse
   tangent;
2. the norm of `H_bi n` relative to a propagated numerical envelope;
3. the image/rank of the resulting boundary compatibility covector;
4. whether that covector is schedule-independent after time reversal.

Only an independently verified `H_bi n=0` licenses gauge-independent
elimination.  A nonzero coupling means the internal equations constrain the
boundary perturbation and requires a constrained, not ordinary, Hamilton
principal form.
