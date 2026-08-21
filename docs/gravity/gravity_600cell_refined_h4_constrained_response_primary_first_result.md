# Primary first result: one tentative class, directional control failed

Date: 2026-08-21

## Provenance

| stage | commit |
|---|---|
| prior-art gate, amended | `8ecbd2a` |
| frozen protocol, final envelope amendment | `be10390` |
| verifier registered before execution | `69ace62` |
| first pre-Hessian implementation failure | `cddb3ca` |
| synthetic correction preregistered | `2ed8157` |
| synthetic correction implemented | `9b2436b` |

The corrected targeted verifier completed with

```text
18/19 PASS
REFINED_H4_CONSTRAINED_RESPONSE_CONTROL_FAILED
```

and wrote

```text
reproducible/gravity_600cell_refined_h4_constrained_response.json
SHA-256 f029260c9ee6e3b763293d237aae27e6ff7c1256eb8bc19c35725084ff385888.
```

No full suite or deferred nonlinear census was run.

## What passed

- all 24 full Hessian ladders passed precision, reality and symmetry;
- all 24 rebuilt the analytic internal null line and frozen nonzero coupling;
- all 24 restricted internal complements were positive, with minimum
  eigenvalue `1.3780099e-5`;
- all constrained solves and full compatibility residuals passed;
- both frozen boundary and internal basis changes preserved the form, with
  maximum differences `2.752e-135` and `5.751e-133`;
- fixed time reversal passed;
- the complete target-free census put all 24 schedules in one class;
- the matrix-corruption and scope controls passed.

These are **TENTATIVE PRIMARY COMPUTATIONAL** observations only.  The frozen
outcome is control failure, so the single-class statement is not accepted.

## What failed

The direct complete-action second-difference check required relative error
below `1e-28`; the maximum was

```text
1.3084747685858131e-20.
```

The largest discrepancy occurs for the all-ones coefficient direction.  The
quadratic response and action estimate are approximately

```text
-11973729.9554093456978641863658
-11973729.9554093456980208596012.
```

The associated lifted tangent has maximum component about `1.5449e5`, so the
largest frozen displacement at `h=1e-10` is about `1.54e-5`.  The copied
absolute accuracy target was not calibrated to these large constrained
lifts.  This suggests fourth-order truncation after one Richardson step, but
that explanation is **OPEN** until a frozen step-halving test distinguishes
it from a real Hessian/action mismatch.

## Required resolution

Preserve this formal failure.  Before changing the primary verifier, register
a separate directional diagnostic using a multi-level centred/Richardson
ladder, an independent precision repeat, asymptotic convergence tests and a
corrupted-quadratic control.  Only a resolved extrapolation to the stored
quadratic response can license a preregistered replacement of the inadequate
fixed `1e-28` gate.

