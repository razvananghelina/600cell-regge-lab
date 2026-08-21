# First adversarial execution failure: refined boundary cotangent

Date: 2026-08-21

The registered adversarial verifier returned

```text
ADVERSARIAL_REFINED_BOUNDARY_COTANGENT_CONTROL_FAILED
11/12 PASS
```

The failed artifact is preserved with SHA-256

```text
7254260e86813f1bae10fbfacd01d6e0eab0bcc80b0ace53f3fc3adc9585bf8b.
```

The sole failure repeats the already diagnosed mixed-precision synthetic
comparison: the artifact records zero kernel pullback change and distance
`0.0000010000...`, but the boolean compared the 100-decimal distance with a
new default-precision `mpf("1e-6")` outside its working context.

All independent scientific comparisons passed, including the symbolic
product-hinge derivative, all 17,040 actual edge incidences, all twelve
boundary components at maximum error `4.05e-77`, the regular coarse control,
and exact curvature/mass-normalized pullbacks.  None is accepted from this
failed run.

The only allowed correction is to construct both the `1e-6` target and
`1e-90` envelope inside a 100-decimal context and freeze the boolean there.
No equation, input, physical threshold or outcome may change.
