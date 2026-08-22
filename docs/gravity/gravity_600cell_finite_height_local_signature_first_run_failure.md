# First local-signature run: unresolved symbolic identity and control crash

Date: 2026-08-22.

Verifier commit: `241a7d4`.

Status: **PRESERVED FIRST RUN; LOCAL THEOREM OPEN**.

The first execution was made only after the prior-art note, protocol and
registered verifier had been committed and pushed. It reported:

```text
[PASS] the local protocol and all accepted inputs are frozen
[FAIL] the derivative and persistent diagonal identities hold symbolically
[PASS] the exact initial diagonal tangency is nondegenerate
[PASS] every visited state has a complete stationary-and-tail root certificate
       states=5
[PASS] every non-diagonal real root has strict physical-gate signs
[PASS] the ordered recursive tree is one, then two, with DEAD and ENTERED_D terminals
       D_entries=1
[PASS] the certified roots reproduce the delayed discovery seeds
```

The run then raised a `ValueError` while converting the bracketed string form
of an Arb midpoint to `mpmath` for the delayed full-action control. No JSON
artifact was written.

The midpoint conversion is an infrastructure defect and can be repaired by
using an explicit Arb decimal export. The failed symbolic identity is
load-bearing and is not classified as infrastructure. Before any corrected
run, the implementation must separately determine whether:

1. `partial_q E=p(q)-pi` is exact but the selected SymPy simplifier failed;
2. the displayed scalar functions use a mismatched normalization; or
3. the stationary partition used by the discovery census is mathematically
   wrong.

Case 3 would invalidate the discovery completeness claim and is the primary
result until repaired. No tolerance, root bracket or numerical agreement may
override the symbolic disagreement.

