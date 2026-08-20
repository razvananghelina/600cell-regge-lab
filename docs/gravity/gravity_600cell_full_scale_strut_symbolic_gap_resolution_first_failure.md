# First coordinate-free gap-resolution failure

Date: 2026-08-20

## Frozen run

- registered source commit: `260e89b`
- failure-artifact commit: `3d43113`
- source SHA-256:
  `fb4a526f2f168c415ceb147829c679c038854516391cd9500f6f089c6a297c03`
- frozen first-failure artifact: `reproducible/gravity_600cell_full_scale_strut_symbolic_gap_resolution_first_failure.json`
- failure-artifact SHA-256:
  `a56ddcf20c0b6ec4aeee1b4e8c076e5c4c955a0780441684fc319ec7c05f6254`
- result: `10/11`, `FULL_SCALE_STRUT_GAP_CONTROL_FAILED`

## Diagnosis

**DERIVED SOFTWARE NEGATIVE.**  The failing assertion required the numerator
of the Euclidean squared component norm of the connection vector to contain
both `tau` and

```text
(lambda-1)^2+3 tau^2.
```

The exact computed norm was instead

```text
3*((lambda-1)^2+3 tau^2).
```

This does not violate the committed protocol.  The protocol requires the
positive quadratic to be factored from the norm and proved nonzero on the
real domain with `tau != 0`; it never requires `tau` itself to divide that
norm.  Requiring it would be geometrically wrong: the connection line can
remain nonzero at `tau=0`, even though `tau=0` is independently excluded by
the four vertex solves.

All result-bearing hostile checks passed in the failed run:

- the pivot-free wedge ideal equalled the disclosed endpoint formula;
- fixed-frame gluing still gave the unit ideal;
- the `D+1` corruption left 96 nonzero wedge residuals;
- no unallowed actual denominator or rank factor remained;
- both exact real points on `lambda+tau-1=0` and `lambda-tau-1=0` retained
  full local ranks, a one-dimensional connection, and the disclosed ideal.

No scientific conclusion is accepted from this failed run.
