# The full off-diagonal Hessian produces the exact Hopf selector cubic

Date: 2026-08-11

Protocol commits:

- `3767638`: complete spectral census preregistered;
- `21a988e`: STEP 1 enumeration committed with no target comparison.

STEP 2 verifier: `reproducible/verify_hopf_full_hessian_spectral_target.py`.
Targeted result: `19/19`.

## Verdict

There is a **STRUCTURAL ADVANCE**:

```text
Tr_(1^perp)(H_X^3) = -23328 C_box(X),
```

where `H_X` is the complete six-label Hessian of the extended cubic and every
off-diagonal label channel is retained.  On `Tr(X^2)=7200`, this single
spectral moment has exactly the six `+Box_i` as its global minima.

This resolves the algebraic action problem that survived the earlier
single-trace audits: the equal-weight projector cubic is now an ordinary
single trace of a completely specified operator.  No diagonal conditional
expectation, label superselection, fitted combination or preferred fibration
is used.

It is not yet a physical vacuum theorem.  The repository has not derived the
extended Hessian as a physical fluctuation operator or a positive coefficient
for its cubic moment in the action.

## Provenance and complete census

Before target comparison, the five characteristic coefficients of the
physical restriction `Hhat_X=H_X|_(1^perp)` were committed in full:

```text
det(lambda I-Hhat_X)
 =lambda^5-e1 lambda^4+e2 lambda^3-e3 lambda^2+e4 lambda-e5.
```

The blind result was

```text
N=5,
e1=0,
e2=-9331200 q,
nonconstant normalized-sphere characters: e3,e4,e5.
```

Power sums were recorded only as Newton-identity checks, not counted as
additional attempts.  Thus the exact look-elsewhere denominator was fixed
before comparison.

Preregistration does not make the investigation fully blind: the spectrum at
one `Box_i` was known from the earlier Hessian audit, and the operator family
itself was chosen because of the Hopf problem.  The commit ordering proves the
narrower, checkable statement that the primitive list was not enlarged or
trimmed after target comparison.

## Exact cubic identity

Let

```text
s_i(X)=Tr(X Box_i),
C_box(X)=sum_i s_i(X)^3.
```

Coefficientwise comparison in the two-dimensional invariant-cubic space
gives

```text
e3(X) = 0*Tr(X^3) - 7776*C_box(X).
```

Because `e1=Tr(Hhat_X)=0`, the third Newton identity gives

```text
Tr(Hhat_X^3)=3e3(X)=-23328*C_box(X).
```

The zero coefficient of the old cubic is load-bearing: the complete Hessian
does not merely reproduce the previously failed mixture.  It lands exactly
on the missing equal-weight invariant line.

## Exact global minimum proof

The overlap map is an isomorphism from the five-dimensional field space to
the zero-sum hyperplane in `R^6`, with

```text
sum_i s_i=0,
sum_i s_i^2=8640 q.
```

At a constrained stationary point of `sum_i s_i^3`, every `s_i` solves the
same quadratic Lagrange equation.  Hence there are at most two distinct
values.  Enumerating their multiplicity `k=1,...,5` is exhaustive.  On
`q=7200`, the exact cubic values are

| `k` | `C_box` |
|---:|---:|
| 1 | `358318080000` |
| 2 | `44789760000 sqrt(10)` |
| 3 | `0` |
| 4 | `-44789760000 sqrt(10)` |
| 5 | `-358318080000` |

The unique maximum class `k=1` has overlap vector

```text
(7200,-1440,-1440,-1440,-1440,-1440)
```

up to permutation.  Injectivity of the overlap map identifies these six
vectors exactly with `+Box_i`; the `k=5` minima are `-Box_i`.

Since the spectral cubic is `-23328 C_box`, its global minima are precisely
the six positive Hopf operators.  This is an exact stationary classification,
not a numerical optimization.

## The other primitive coefficients

The quartic reduces exactly to

```text
e4=15672832819200 q^2+139968 sum_i s_i^4.
```

An exhaustive stationary-root census for the fourth moment gives the sharp
bound

```text
sum_i s_i^4 <= (7/10)(sum_i s_i^2)^2.
```

Equality occurs only for one-versus-five values.  Therefore the twelve
signed points `+/-Box_i` are exactly the global maxima of `e4`.  Using them as
minima would require changing the sign of this primitive and is not silently
counted as the same physical claim.

The determinant coefficient `e5` is independent of

```text
q Tr(X^3), q C_box(X), sum_i s_i^5.
```

All `+Box_i` are exact strict local maxima of `e5`, and all `-Box_i` strict
local minima.  A global certificate has not been obtained.

**OPEN:** global uniqueness for `e5`.  Numerical reconnaissance and an exact
solver returning `unknown` are not evidence and are not used in the verdict.

## Look-elsewhere accounting

- Primitive spectral coefficients fixed in STEP 1: `N=5`.
- Nonconstant coefficients on the normalized sphere: `3`.
- Certified Hopf extremal loci: `e3,e4`, hence `2/5` overall and `2/3` among
  nonconstant primitives.
- Coefficients whose ordinary positive-sign minimization selects exactly the
  six `+Box_i`: only `e3`, hence `1/5` overall and `1/3` among nonconstant
  primitives.
- `e5` is excluded from the hit count until its global problem is certified.

The small fixed denominator makes the `e3` identity more than a large-search
coincidence, but the non-blind provenance still requires the label
**STRUCTURAL ADVANCE**, not a prediction.

## Gaussian and physical-action audit

At every `+Box_i`, the physical Hessian has exact signature

```text
(3 positive, 2 negative, 0 zero).
```

The omitted constant label direction is the separate universal zero mode.
Thus a real bosonic Euclidean Gaussian with quadratic operator `Hhat_X`
diverges.  A fermionic determinant, an absolute determinant, or replacing
`Hhat_X` by `Hhat_X^2` would be new input.

The cubic trace itself is perfectly well-defined as finite spectral algebra.
What remains **OPEN** is why the theory's action contains
`+Tr(Hhat_X^3)`.  It is not derived from a convergent bosonic one-loop
integration merely because it is a spectral invariant.

## Status ledger

- **DERIVED:** `e3=-7776 C_box` and
  `Tr(Hhat_X^3)=-23328 C_box` coefficientwise.
- **DERIVED:** the positive-sign cubic spectral moment has global minima
  exactly at the six `+Box_i` on `q=7200`.
- **DERIVED:** `e4` has global maxima exactly at the twelve signed vertices.
- **OPEN:** the global extremal classification of the independent quintic
  `e5`.
- **DERIVED NEGATIVE:** `Hhat_Box` is indefinite, so the selector is not
  generated by a convergent real bosonic Gaussian.
- **STRUCTURAL ADVANCE:** the full nondiagonal Hessian supplies the previously
  missing single-operator realization of the selector, bypassing rather than
  deriving label superselection.
- **OPEN PHYSICS:** derive `Hhat_X` as a licensed fluctuation operator and
  derive the presence and positive coefficient of its cubic spectral moment.

## Consequence

The earlier statement that off-diagonal Hessian channels merely obstruct the
selector is too strong and must be corrected.  They obstruct identifying the
Hessian itself with the diagonal Gram operator `D_aux`, but their complete
cubic spectral invariant produces the correct selector exactly.

The next gate is no longer algebraic recognition.  It is dynamical: determine
whether the existing spectral-action or fluctuation machinery produces
`+Tr(Hhat_X^3)` without adding a cutoff, sign or field after seeing this
result.
