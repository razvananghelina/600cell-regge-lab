# Preserved first run: structural SymPy equality rejected an exact identity

Date: 2026-08-19

The first execution of the preregistered target-blind corrected-strut carrier
verifier returned

```text
CORRECTED_STRUT_CARRIER_CONTROL_FAILED
TOTAL: 12/13 tests PASSED
```

Its JSON artifact has SHA-256

```text
dee2581931b682fd70ea0b2828ff6d97ea417455373e7b467ddb9a451fbacfe1
```

All geometric, rational-block, incidence, equivariance, collective,
corruption and 100-digit/binary64 controls passed.  The only failed predicate
printed

```text
differential residual             [0,0]
natural Jacobian determinant      2*(lambda-1)
static natural rank               1
```

but compared the determinant with another construction of
`2*(lambda-1)` using SymPy's structural `==`.  The two printed expressions
are mathematically equal but have different internal construction trees:

```text
generic_det == 2*(lambda-1)                         False
simplify(generic_det-2*(lambda-1)) == 0              True
```

The disclosed mathematical gate is equality after exact simplification, not
equality of SymPy expression trees.  Before rerunning, the implementation may
replace only this structural comparison by the simplified-zero comparison.
No carrier coefficient, matrix row, spectrum, tolerance, outcome hierarchy
or target firewall may change.

The failed artifact is intentionally retained.  It already contains the
complete target-blind census but cannot be accepted while the frozen control
harness exits nonzero.

