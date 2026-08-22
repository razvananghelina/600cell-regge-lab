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

## Resolution frozen before the second run

The disagreement is case 1, with a specific auditable cause. SymPy reduces

```text
partial_q E-(p-pi)
```

to a rational prefactor times

```text
sqrt(9*q^6+84*q^4+256*q^2+256)
-(3*q^2+8)*sqrt(q^2+4).
```

The exact polynomial identity

```text
9*q^6+84*q^4+256*q^2+256
=(3*q^2+8)^2*(q^2+4)
```

and the strict positivity of both `3*q^2+8` and `sqrt(q^2+4)` for real `q`
prove that this remainder is zero. The corrected verifier records the
factorization and positivity queries separately before replacing the radical;
it does not accept numerical agreement as the identity proof.

The control crash is corrected only by replacing `str(arb.mid())`, which
includes brackets, with Arb's explicit fixed-digit decimal export. Neither
correction changes a scientific criterion, root bracket, precision or
expected terminal.

## Second run: remaining Arb export failure

Verifier commit: `a9c68dd`.

The second execution passed the symbolic factorization, the nondegenerate
diagonal check, all five complete root censuses, every strict physical gate,
the exact terminal tree and the delayed seed control. It then reached the same
full-action control and failed because `arb.mid().str(140)` still emitted a
bracketed ball at the verifier's 220-digit context.

This failure occurs after the primary certificate is assembled but before the
control and artifact are completed. A pre-commit infrastructure check showed
that even `lower().str()` can retain Arb's rounding-radius notation. The
narrow correction therefore asks Arb for the 140-digit midpoint rendering and
extracts its explicit centre token from `[centre +/- radius]`. This changes
only the decimal representative supplied to the non-load-bearing `mpmath`
action residual control; the primary certificate continues to use complete
outward-rounded balls.
