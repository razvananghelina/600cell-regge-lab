# First symbolic adversarial run: factor-classifier failure

Date: 2026-08-20

The verifier registered in commit `1fb14f3` passed provenance, the exact
regular geometry, the one-cell ideal and the face-connection controls.  It
then printed

```text
agrees=True,
basis=[
  (B lambda^2-2 B lambda+B-2 lambda^2+4 lambda-2 tau^2-2)/(lambda-1)^2,
  (D lambda-D-lambda)/(lambda-1)
]
```

for the generic shared-face ideal.  Before any convention or corruption
control completed, it raised

```text
CoercionFailed: expected Rational object, got -3*D - 1
```

while collecting exceptional `(lambda,tau)` factors.  No JSON artifact was
written.

The failure is mechanical: `factor_set` factored expressions containing
`A,B,C,D` but then attempted to coerce every factor into `Q[lambda,tau]`.
The preliminary ideal equality is not accepted because the protocol did not
finish.

The authorized repair is exact and target-independent:

1. factor each numerator and denominator over
   `Q[A,B,C,D,lambda,tau]`;
2. retain only nonconstant factors whose free symbols are a subset of
   `{lambda,tau}`;
3. normalize only those retained factors in `Q[lambda,tau]`;
4. change no geometry, residual, ideal, allowed-factor list, convention
   control or outcome rule.

