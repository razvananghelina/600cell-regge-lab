# Corrected precision run: serialization-only failure

Date: 2026-08-20

After the conditioning-aware correction committed in `dd30486`, the targeted
verifier printed:

```text
[PASS] both carrier matrices reproduce the frozen binary64 diagnostics
[PASS] row-order reversal preserves both direct spectra
[PASS] deleting the first pole coefficient strictly weakens both carriers
[PASS] both high-precision comparisons are evaluated and classified
       decisive=True
[PASS] the preregistered precision hierarchy assigns one outcome
       FULL_SCALE_STRUT_PRECISION_RESOLVED
```

It then raised

```text
TypeError: Object of type bool is not JSON serializable
```

inside `json.dumps(payload, ...)`.  No output artifact was written.  Under
NumPy 2, `numpy.bool_` reports the class name `bool` in this exception but is
not a Python built-in `bool` accepted by the standard JSON encoder.

**DERIVED HARNESS FAILURE.** The mathematical and numerical computation
completed, but its result is not accepted or frozen because serialization
failed.

The authorized repair is limited to a strict JSON fallback converting only
`numpy.bool_` to Python `bool` and raising on every other unknown type.  No
matrix construction, precision, threshold, control, outcome rule or payload
field may change.  The preserved `6/8` first-failure artifact must continue
to pass its hash check before the rerun.

