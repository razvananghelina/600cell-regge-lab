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

Because `runpy` does not return the executed namespace after an exception,
the implementation may temporarily replace `sys.exit` around exactly this
frozen upstream call by an audited function that returns only for codes
`None` or `0` and raises `SystemExit` for every nonzero code.  It must restore
the original function in `finally`.  No carrier, target matrix, comparison,
tolerance, look-elsewhere count or outcome rule may change.

No corrected target angle was evaluated in this failed run.
