# First execution failure: refined local curvature mass verifier

Date: 2026-08-21

The first registered execution of
`verify_gravity_600cell_refined_local_curvature_mass.py` stopped with the
frozen control outcome

```text
REFINED_LOCAL_CURVATURE_MASS_CONTROL_FAILED
13/15 PASS
```

The preserved artifact has SHA-256

```text
43cdfd58a4bea192255c9ca92975fe08754a012506ffd8756f9676bc1a451321.
```

No scientific result is accepted from this execution.

## Failure 1: malformed frozen source hash

The verifier copied the feasibility-source hash as a 62-character string:

```text
36fba835048e6f0676b749192a9d882406932770a00ba1396929bbc4d04a32
```

The actual unchanged SHA-256, correctly frozen in the new protocol as well as
the earlier stationary-fill protocol, is

```text
36fba835048e6e0f0676b749192a9d882406932770a00ba1396929bbc4d04a32.
```

The transcription omitted `e0`; it did not expose a changed scientific input.
The provenance guard correctly failed.

## Failure 2: premature low-precision conversion of `tau0`

The verifier constructed

```python
TAU = mp.mpf("0.0102")
```

at module import, before entering the 100-decimal working context.  Mpmath
therefore stored only its default-precision approximation.  Reusing that
rounded object inside the high-precision block produced a maximum identity
error `2.7445033e-18`, rather than testing the exact decimal input required by
the protocol.

This is a numerical implementation defect, not evidence against the identity.
The geometry, positivity, schedule, rank, mass-conservation and all negative
controls passed.

## Allowed correction

Before rerunning:

1. correct the malformed hash in the verifier;
2. retain `tau0` as the string `"0.0102"` and construct `mp.mpf(TAU_TEXT)`
   only inside each high-precision context;
3. keep every equation, tolerance, input file, alternative mass vector,
   corruption control and outcome unchanged.

The failed artifact remains part of the record.
