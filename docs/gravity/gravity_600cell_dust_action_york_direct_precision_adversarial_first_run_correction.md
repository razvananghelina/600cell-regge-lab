# First-run correction: missing sparse-helper dependency

Date: 2026-08-19

The first execution of
`verify_gravity_600cell_dust_action_york_direct_precision_adversarial.py`
stopped after the provenance, primary-outcome and imported-geometry controls,
before constructing a rigidity matrix, carrier or adversarial residual.

The AST-extracted audited helper `incidence_data` calls
`scipy.sparse.csr_matrix`.  The standalone wrapper omitted the corresponding
`import scipy.sparse as sp`, producing:

```text
NameError: name 'sp' is not defined
```

No JSON artifact and no scientific label were produced.  The frozen repair is
to add only that missing import.  No input hash, matrix, rank rule, error rule,
selected sector or outcome criterion changes.
