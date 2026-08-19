# Preserved first execution: upstream `SystemExit(0)` ended the wrapper

Date: 2026-08-19

The first execution of
`reproducible/verify_gravity_600cell_corrected_strut_alignment.py` printed the
new provenance pass and began the frozen old-alignment reconstruction, then
returned process code zero without executing a corrected comparison or
writing a new artifact.

Cause: the audited upstream verifier ends with

```python
sys.exit(0 if passed == tests else 1)
```

and `runpy.run_path(...)` propagated its `SystemExit(0)` through the wrapper.
The operating-system success code was therefore not evidence that the new
verifier completed.  The absence of
`gravity_600cell_corrected_strut_alignment.json` is the decisive control.

Before rerunning, the implementation may catch `SystemExit` around exactly
this frozen upstream call, accept only exit codes `None` or `0`, and re-raise
every nonzero code.  No carrier, target matrix, comparison, tolerance,
look-elsewhere count or outcome rule may change.

No corrected target angle was evaluated in this failed run.

