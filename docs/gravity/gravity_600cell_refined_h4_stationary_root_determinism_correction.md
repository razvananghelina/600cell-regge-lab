# Preregistered correction: deterministic root-search artifact

Date: 2026-08-20

First complete result commit: `a57418b`.

The fast-evaluator protocol simultaneously required per-attempt elapsed time
and byte-identical deterministic JSON. Wall-clock measurements make exact
artifact reproduction impossible and do not contribute to the scientific
classification.

Remove only `elapsed_seconds` from serialized attempt records. The verifier
may print total runtime to the console, but must not include clocks,
timestamps beyond the fixed date, process identifiers or unordered runtime
metadata in the JSON.

Do not change any equation, evaluator, cross-check, seed, bound, solver
option, endpoint, residual, diagnostic, threshold, attempt count, hit count or
outcome. Run the corrected targeted verifier twice and require identical
SHA-256 hashes. If any scientific field differs from the first complete
artifact beyond deleting the 120 elapsed-time fields, stop and record a
scientific reproducibility failure.
