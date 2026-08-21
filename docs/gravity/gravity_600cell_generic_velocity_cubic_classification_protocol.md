# Frozen adjudication: global classification of the cubic mismatch

Date: 2026-08-21

Original cubic protocol commit: `d2efdf4`.

Unexecuted verifier registration commit: `628b769`.

First completed execution commit: `fcd5af5`.

First artifact SHA-256:
`b594a924fd2ebcfeb2a59c7734158b50063228603e6170890ec96f37889c0cc0`.

Status: frozen after the first exact run returned
`GENERIC_CUBIC_OPEN`, before changing its classification logic or outcome.

## 1. Preserved first-run observation

The first verifier derived both on-shell coefficients independently of any
classification rule.  They are affine in `c`, their slopes reproduce the
accepted first-correction slopes, and their exact cross-resultant reduced to

```text
Delta(v)
=C2(v,0) coefficient_c(P2)-P2(v,0) coefficient_c(C2)
=129600 [2*pi-5*theta(v)]^2/(v^2+4),

theta(v)=acos((v^2+2)/(2(v^2+3))).
```

All preregistered direct arbitrary-precision coefficient controls passed.
The first implementation deliberately left the zero set unresolved and
therefore labelled the result `OPEN`.

## 2. Frozen exact positivity proof

Put `x=v^2>=0` and

```text
epsilon(x)=2*pi-5*acos((x+2)/(2(x+3))).
```

Use the already-certified exact identities

```text
epsilon'(x)=5/[(x+3)sqrt(x+4)sqrt(3x+8)]>0,
cos(2*pi/5)=(sqrt(5)-1)/4<1/3.
```

Because `acos` is strictly decreasing on `[-1,1]`, the second identity gives

```text
acos(1/3)<2*pi/5,
epsilon(0)>0.
```

Strict monotonicity then gives `epsilon(x)>0` for every `x>=0`.  Since
`x+4>0`, the exact cross-resultant is strictly positive everywhere on the
real velocity line.

## 3. Consequence for the registered domain

On the registered domain

```text
v!=0,
K(v^2)!=0,
```

the lapse slope is

```text
coefficient_c(C2)
=1440 v K(v^2)/[
  sqrt(v^2+4)sqrt(3v^2+8)(v^2+3)(v^2+4)
],
```

and is nonzero.  Therefore the lapse equation has exactly one `c`.  On that
root the momentum residual is the nonzero cross-resultant divided by the
nonzero lapse slope, up to the registered sign convention.  A zero or
degree-drop of the momentum slope cannot rescue compatibility.

Hence there is no common real `c` anywhere on the complete registered
domain.  The excluded pair `K=0` has no quadratic endpoint jet and cannot be
reintroduced at cubic order; `v=0` belongs to the separately certified
turning-point stratum.

## 4. Permitted correction and fixed interpretation

The correction may only:

1. certify the exact factorization of the first-run cross-resultant;
2. certify the positivity proof above;
3. replace the unresolved classification with
   `NO_COMMON_C_ON_COMPLETE_REGISTERED_DOMAIN`;
4. assign the preregistered outcome
   `GENERIC_CUBIC_FIXED_STATE_OBSTRUCTION`.

It may not change the action, physical state, endpoint ansatz, coefficient
extraction, sample points, numerical tolerances or hostile controls.

The result is a **DERIVED NEGATIVE, scoped**: the fixed finite 600-cell action
does not possess the registered arbitrary-small fixed-state endpoint jet at
cubic order.  It does not select an isolated positive duration and therefore
does not derive a tick.  A mechanically independent adversarial derivation
is still required before consolidation.
