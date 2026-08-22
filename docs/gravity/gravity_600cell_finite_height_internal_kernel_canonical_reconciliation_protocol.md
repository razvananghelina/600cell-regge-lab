# Protocol: exact canonical reconciliation of the internal kernel

Date: 2026-08-22.

Prior-art/framing commit: `5c30454`.

Status: **FROZEN AFTER DISCLOSURE OF THE EXPLORATORY RATIO AND BEFORE THE
REGISTERED SYMBOLIC/PROJECTOR CALCULATION.**

## 1. Frozen inputs

Reject any hash mismatch:

```text
docs/gravity/gravity_600cell_finite_height_internal_kernel_canonical_reconciliation_prior_art.md
  9e474113a426c61e44fcce000ff6e4a0262bc00c951a09dd72eb3cca347d0c4b

reproducible/gravity_600cell_finite_height_internal_carrier_rank.json
  513fdea33f6b868efa6d6f2b2526bade7ce615ea949f955588916a8d0baee0c8

reproducible/gravity_600cell_finite_height_internal_carrier_rank_matrices.npz
  97f5b8318be2b3ccf843db87e678ac1ac6ce402db262023c6bbc63a7b647321b

reproducible/gravity_600cell_finite_height_internal_carrier_rank_adversarial.json
  ddd4704b7d1deb6360e752b2ebfe5cc0b66d03819c9f8df7b74e24373aa98fb5

reproducible/gravity_600cell_finite_height_internal_carrier_rank_adversarial_matrices.npz
  45ee09642485dc0e18c6378b9454882414bb07d7e5ad9dbd0c3a8896fd8a7f74

reproducible/gravity_600cell_finite_height_selector_audit.json
  956cd655b8b3a5106029fb852df74b85bb59f922a4984542bc2e089f54799676

reproducible/verify_gravity_600cell_finite_height_selector_audit.py
  aca44fa0cf0f6a464ca1a8eaa61356941e8408950ae3801085b85ce314741503

docs/gravity/gravity_600cell_finite_height_carrier_quadratic_protocol.md
  f73ee892258e33d43991fc8c74bc6f44e6c7f2ae57be56f057050c86ff646fad

reproducible/gravity_600cell_finite_height_carrier_quadratic.json
  0ec142bfc68d04498992a6cdba7437933560b860244573d187cb6e018ece78f9
```

Require the primary rank outcome `25/25`, the adversarial outcome `19/19`,
their full nullities `1`, diagonal nullities `121`, and their stored parity
kernel agreements.  Require the selector audit to retain its exact determinant
certificate and `10/10` outcome.

## 2. Exact symbolic calculation

Before reading any stored kernel vector, use independent SymPy symbols
`h,q,u,u',p'` and reconstruct

```text
C_h =4*pi*q*u,
C_q =8*pi*u' +4*pi*h*(u+q*u'),
P_h =-2*pi*u,
P_q =p' -2*pi*h*u'.
```

Prove

```text
det partial(C,P)/partial(h,q)
 =8*pi^2*h*u^2 +4*pi*u*(4*pi*u'+q*p').
```

Then substitute the exact state identity

```text
4*pi*u'+q*p'=0
```

and require the determinant to become `8*pi^2*h*u^2`.

For the carrier coordinates

```text
sigma=delta lambda,
c=delta log rho,
lambda=1+h*q,
rho=h^2,
```

derive rather than assume

```text
partial(h,q)/partial(sigma,c)
 = [[0,h/2],[1/h,-q/2]],
det = -1/2.
```

Therefore prove

```text
det partial(C,P)/partial(sigma,c)=-4*pi^2*h*u^2 !=0
```

for `h>0,u>0`.

The exact `C`-tangent ratio is

```text
c/sigma = -2*C_q/[h*(h*C_h-q*C_q)].
```

Require its denominator to be nonzero at the frozen state.  Prove directly
that the vector `(1,c/sigma)` annihilates `dC` and not `dP`.

## 3. Independent numerical bridge

Reconstruct from their definitions

```text
epsilon(q)=2*pi-5*acos((q^2+2)/(2(q^2+3))),
mu(q)=180*epsilon(q)/(pi*sqrt(q^2+4)),
p(q)=180*q*epsilon(q)/sqrt(q^2+4)
     -600*sqrt(3)*asinh(q/sqrt(8(q^2+3))).
```

Evaluate derivatives and the exact ratio independently at 120 and 180
decimal digits.  Require relative agreement below `1e-100`, `dC` residual
below `1e-100`, the determinant identity below `1e-100` relatively, and
`|dP|>1e-3` for the representative `sigma=1` tangent.

The exploratory value

```text
c/sigma = 0.4589898592210...
```

is disclosed and is not a blind target.

## 4. Four stored projector comparisons

Only after the symbolic and independent numerical ratio are fixed, load the
primary and adversarial rank matrix artifacts.  Form the normalized physical
240-vector

```text
(1 repeated 120 times, (c/sigma) repeated 120 times)
```

and its rank-one projector.  Compare it with the even and odd projectors from
both constructions.

Use the maximum stored projector uncertainty from the two rank artifacts.
Agreement requires distance at most ten times that uncertainty; dependence
requires more than one hundred times it; otherwise the result is open.

## 5. Hostile coordinate controls

Repeat the analytic-projector construction with each wrong convention:

1. `sigma=delta log lambda` instead of `delta lambda`;
2. `c=delta rho` instead of `delta log rho`.

Each wrong projector must differ from the accepted adversarial projector by
more than one hundred times the frozen uncertainty.  The symbolic coordinate
determinant must also reject omission of the factor `1/2` in
`delta h=h*c/2`.

## 6. Logical closure

If the complete internal kernel is exactly the `C` tangent and `dP` is
nonzero on it, conclude

```text
ker R_p intersect ker dP = {0}
```

for both parities.  Because every internal survivor is already known to lie
on that one homogeneous line, no untested nonhomogeneous momentum derivative
is needed for this fixed-input intersection statement.

This establishes local isolation of the output at fixed incoming canonical
data.  It does not say that the forced derivative with respect to varying
incoming data is zero.

## 7. Outcome hierarchy

1. `CANONICAL_RECONCILIATION_CONTROL_FAILED` for provenance, exact formula,
   positivity, branch or hostile-control failure.
2. `CANONICAL_RECONCILIATION_OPEN` for a projector comparison in its open
   band or inadequate precision separation.
3. `INTERNAL_KERNEL_NOT_THE_HOMOGENEOUS_CONSTRAINT_TANGENT` for a resolved
   projector disagreement with all controls passing.
4. `INTERNAL_KERNEL_IS_LAPSE_CONSTRAINT_TANGENT_FIXED_INPUT_REMOVES_IT` only
   if every exact, numerical, projector and hostile gate passes.

Outcome 4 is **DERIVED EXACT/COMPUTATIONAL** and a bounded negative for reading
the internal line itself as fixed-input evolution.  The next physically
meaningful object would be a forced canonical/Jacobi response to varying
incoming data.

Run only the new verifier and static registry checks.  Do not run the full
suite.
