# First adversarial outcome: formal null-coupling disagreement from an impossible threshold

Date: 2026-08-21

The registered adversarial verifier returned

```text
ADVERSARIAL_REFINED_H4_NULL_COUPLING_DISAGREEMENT
11/11 PASS
```

with artifact SHA-256

```text
4b7f4866d6a1480f5006e8484325e9bf5d578106325fb19ad353c834a60062aa.
```

All symbolic, incidence, curvature-balance and corruption controls passed.
The independently derived row is

```text
c_rs=tau0*C_rs/8
```

and its maximum difference from every primary row is `6.493e-44`.

The frozen adversarial protocol incorrectly required `<1e-68`, even though
the primary finite-difference artifact had already certified per-schedule
envelopes between roughly `6.14e-41` and `9.74e-41`.  The strict requirement
therefore asked the adversarial exact formula to agree more accurately than
the comparison data were known.  The observed difference is approximately
three orders of magnitude *inside* every relevant primary envelope.

This remains the formal first adversarial outcome and is not overwritten.
It is not evidence for a physical disagreement: the test criterion was
incompatible with the frozen uncertainty of its input.

The only admissible correction is to compare each adversarial row with the
explicit `envelope` stored alongside that primary schedule.  The independent
boundary pre/post reconstruction retains its `<1e-68` gate because both of
those inputs support that precision.  No formula, geometry, component,
factor, sign or corruption threshold may change.

