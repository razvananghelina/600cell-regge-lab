# Refuted carrier with one remaining harness failure

Date: 2026-08-19

The first execution after commit `cbb91ed` assigned the scientifically correct
outcome `CANONICAL_DATA_VERTEX_CARRIER_REFUTED`.  It preserved every exact
residual from the original failure and changed the alternate-graph check from
FAIL to PASS because both graphs agree that the candidate is not included.

The resulting artifact
`reproducible/gravity_600cell_canonical_data_carrier.json` has SHA-256
`4a9fbbfea5287de10e4a396a55dafbe6c31a58d2c44856e068bad92a2fef6347`.
It records 9/10 checks:

```text
OUTCOME: CANONICAL_DATA_VERTEX_CARRIER_REFUTED
(lambda,tau)=(2,5): 3600 nonzero inclusion rows
(lambda,tau)=(3,11): 3600 nonzero inclusion rows
baseline/alternate decision: agrees
```

**DERIVED NEGATIVE.** The candidate is refuted exactly and independently of
the right-inverse graph.

**STRUCTURAL SOFTWARE DEFECT.** The remaining red check asks the target
proposition itself to be true.  Thus a preregistered negative scientific
outcome cannot exit successfully.  The separate
`gravity_600cell_canonical_data_carrier_outcome_check_correction.md` freezes
the replacement by an accounting check before another execution.  It changes
no scientific value.
