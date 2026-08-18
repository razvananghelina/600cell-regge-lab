# Result: arbitrary-precision solver repair stopped at calibration

Date: 2026-08-14

Prior-art gate: `6b7f9e4`

Protocol commit: `4b6b10c`

Implementation commit: `f7a70fa`

Targeted verifier:
`reproducible/verify_gravity_600cell_dust_action_solver_repair.py`

Machine-readable result:
`reproducible/gravity_600cell_dust_action_solver_repair.json`

## Status ledger

| Item | Status | Evidence |
|---|---|---|
| Imported complete-action control | **DERIVED CONTROL** | `14/14` retained |
| New verifier implementation gates | **DERIVED COMPUTATIONAL** | `5/5` |
| Target rows evaluated | **DERIVED NEGATIVE COUNT** | exactly `0` |
| Branch audits | **DERIVED CONTROL** | both pairs pass in both parities |
| Imaginary contamination | **DERIVED CONTROL** | below `3.16e-78`, threshold `1e-70` |
| Operational/validation self-agreement | **DERIVED CONTROL** | component difference `2.186e-37`; transverse difference `1.173e-36` |
| Agreement with old 60-decimal reference | **FAILED CALIBRATION GATE** | `3.544337829e-9 > 1e-10`, both pairs and parities |
| Complete-action transverse roots | **OPEN NUMERICALLY** | forbidden to evaluate after calibration failure |
| Reduced scalar / physical tick | **OPEN** | no target root and no scalar classification |

Global outcome:

```text
DERIVATIVE_CALIBRATION_FAILED
```

This is a calibration failure, not evidence that a transverse root is absent.

## 1. Exact outcome

For both even and odd schedule parities, the operational and validation pairs
passed every gate except comparison with the old control derivative:

```text
operational reference error = 3.5443378290216933e-9
validation  reference error = 3.5443378290216933e-9
frozen threshold            = 1e-10
```

The new calculations agreed with each other much more closely:

```text
maximum component difference     = 2.1863876919233215e-37
component stability proxy        = 2.7329846146308530e-27
transverse difference norm       = 1.1728023260836656e-36
transverse stability-proxy sum   = 1.4660029074579816e-26
```

At the symmetric control the operational transverse norm was
`1.4660e-37`; the validation norm was `1.3194e-36`.  These are controls only,
not target results.

Maximum imaginary contamination was:

```text
even: 2.95e-82
odd : 3.15e-78
```

All 71 branch audits per row passed.  The minimum Gram margins were
`1.0404e-4` (even) and `5.2020e-5` (odd); the minimum angle-argument modulus
was `0.99534`.

The verifier obeyed the preregistered stop rule and wrote
`target_rows_evaluated = 0`.

## 2. Framing attack: what actually failed

The frozen `1e-10` reference threshold was not justified by the precision of
the reference.  The old 60-decimal control uses a second-order additive
central difference at relative step `3e-6`, and its own registered comparison
with the binary analytic gradient is only:

```text
even: 1.105642964e-8
odd : 1.261789290e-8
```

with a frozen upstream acceptance threshold `5e-8`.  A new derivative cannot
be required to agree with that object at `1e-10` and simultaneously treat the
old object as the reference.  The preregistration incorrectly inferred a
small truncation constant merely from `(3e-6)^2`; the coefficient multiplying
that power was unknown.

There is also a coordinate distinction.  The old control differences the
action additively in `u` and is later multiplied by `u/24`.  The new method
differences symmetrically in `log(u)`.  Both converge to the same derivative,
but their finite-step `h^2` terms are not identical.

Therefore:

- **DERIVED:** the preregistered calibration gate failed;
- **DERIVED:** the two new disjoint derivative pairs are internally
  consistent far below their stability proxies;
- **STRUCTURAL:** the failure is attributable to comparing against a
  lower-accuracy reference with an unjustified tighter threshold;
- **OPEN:** the new target solver, roots and scalar remain untested.

This does not license deleting the failed run or silently relaxing its gate.

## 3. Post-result prior-art update

Centered finite-difference accuracy is an asymptotic statement with a
derivative-dependent coefficient, not a universal numerical bound obtained by
writing `h^2`.  Standard finite-difference theory explicitly states orders of
accuracy while retaining function-dependent truncation terms.  A relevant
primary source is B. Fornberg, *Generation of Finite Difference Formulas on
Arbitrarily Spaced Grids*, Mathematics of Computation 51 (1988), 699--706,
<https://doi.org/10.1090/S0025-5718-1988-0935077-0>.

Nothing in that literature supplies the repository-specific root or scalar.
External novelty remains **OPEN**.

## 4. Legitimate continuation

A separate protocol may replace only the invalid cross-reference gate by the
already registered upstream accuracy class `5e-8`, while retaining unchanged:

- both new derivative pairs;
- their self-agreement gate;
- all target states and starting points;
- the solver, active-merit switch, damping list and iteration limit;
- independent validation, scalar thresholds and look-elsewhere count.

This repair is target-blind: no deformed target row, root or scalar was
evaluated in the failed experiment.  The observed `3.544e-9` must be reported
but must not become the new threshold.

