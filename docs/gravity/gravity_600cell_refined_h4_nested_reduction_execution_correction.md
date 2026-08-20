# Execution correction: resumable nested H4 census

Date: 2026-08-20

## Status before the correction

The preregistered nested `6+3+1` verifier was committed as `5dcdcac`
before its first scientific execution.  That execution passed all eight
pre-search controls.  The user directly observed it reach `12/12` after
approximately 2.5 hours and report no root before the WSL / tool session
failed.  The recoverable tool transcript available after the failure ended
at `2/12`, and no JSON artifact survived.  Therefore the completed-run report
is recorded, but its scalar values, class signs, exact final outcome string
and test count cannot be independently recovered from the filesystem.

**DERIVED OPERATIONAL NEGATIVE:** the all-at-end write strategy is not robust
enough for this calculation.  Approximately 2.5 hours of completed CPU work
and its artifact were lost.

**OPEN / USER-WITNESSED RESULT:** the completed run found no root.  This is
useful recovery evidence but is not promoted to `DERIVED COMPUTATIONAL`
because the exact artifact and a reproducing run are absent.

## Frozen correction before the next execution

Add an atomic checkpoint after every completed schedule class.  A checkpoint
is reusable only when all of the following match exactly:

- the verifier's own SHA-256;
- every frozen input SHA-256 already required by the protocol;
- prior-art commit `7714933`;
- protocol commit `b284aa1`;
- the exact 19-point continuation order.

On a matching restart, reuse only whole completed classes and continue with
the next class in numerical order.  Never checkpoint or reuse a partial
class.  Preserve the exact class-local continuation, two inner seeds,
fallback rule, tolerances, bisection, high-precision validation and outcome
hierarchy from the committed protocol.  The correction changes persistence,
not the mathematical search.

Provide `H4_NESTED_FORCE_RECOMPUTE=1` as a replication mode which ignores an
otherwise valid checkpoint and recomputes all twelve classes.  The final
accepted result still requires a complete forced recomputation to reproduce
the scientific fields of the first complete artifact.  Merely replaying a
checkpoint is not an independent rerun.

## Epistemic scope

- **STRUCTURAL:** checkpoint/resume does not enlarge or shrink the search.
- **OPEN:** the user-witnessed no-root outcome awaits an artifact-producing
  reproduction before it can be accepted as a computational result.
- This correction is committed before the next scientific execution and
  before inspecting any class scalar values.
