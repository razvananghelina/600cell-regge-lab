# Regular-lapse identity: first exactness-control failure

Date: 2026-08-16

Prior-art gate: `898d123`

Frozen protocol: `cf492b9`

Implementation registered before evaluation: `55d953c`

First artifact SHA-256:
`2e796705fba06803b6b7d2c4964b2ad3ecdf0860b893626fe42db9b0c3dc5344`

Status: **IMPLEMENTATION EXACTNESS FAILURE; SCIENTIFIC OUTCOME NOT CLOSED**.

The first targeted run returned `10/13` and mechanically classified the
result as `REGULAR_LAPSE_PATTERN_ONLY`.  All twelve preregistered 100-digit
evaluations passed.  Their largest errors were

```text
internal gradient absolute error  6.354e-95
pre-momentum relative error        1.121e-90
post-momentum relative error       5.774e-94
action relative error              1.831e-92
maximum imaginary contamination   9.253e-97
```

The exact action, dust, pole, diagonal and pre/post formula flags also all
passed.  The three failed checks were confined to the symbolic Gram/angle,
log-cut and angle-product certificates.

## Cause

The symbolic angle implementation formed its Gram derivative with

```python
(integer expression)/2
```

where the numerator was a Python integer.  This inserted binary floats into
an otherwise exact SymPy matrix.  Consequently exact identities appeared as
`1.000...` rather than rationals and could not satisfy equality-to-zero
gates.  The same construction is harmless in the upstream arbitrary-
precision numerical evaluator because its matrices are `mpmath` matrices;
it is invalid in the new exact route.

## Allowed correction

Replace only that division by multiplication with `sympy.Rational(1,2)`.
No carrier, formula, domain, control point, tolerance, branch classification
or outcome rule changes.  Preserve this first artifact before rerunning.

