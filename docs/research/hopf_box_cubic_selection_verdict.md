# The existing Hopf--Box cubic does not select the six fibrations

Date: 2026-08-10

Protocol commit: `2b754de`.  Registered verifier:
`reproducible/verify_hopf_box_cubic_selection.py`.  Targeted result: `14/14`.

## Complete hypotheses of the tested claim

Use exactly the already derived real operator space

```text
W = span_R{Box_i : i=0,...,5},
Box_i = 6 A_f,i - A,
sum_i Box_i = 0,
dim_R W = 5.
```

No new weights, action terms or normalization conditions are added.  The two
existing invariant constraints are

```text
q(X) = Tr(X^2) = 7200,
f(X) = Tr(X^3) = 14400.
```

The question frozen before the calculation was whether their simultaneous
real level set consists of exactly the six derived `Box_i`.

## Exact coordinate forms

In the basis `E_a=Box_a-Box_5`, `a=0,...,4`, write
`X=sum_a u_a E_a`.  The exact Gram matrix has diagonal `17280` and
off-diagonal `8640`, hence

```text
q(u)/17280 = sum_a u_a^2 + sum_{a<b} u_a u_b.
```

The exact cubic has 15 nonzero monomials:

```text
f(u)/(-155520) =
    u0^2 u1 + u0^2 u2 + u0 u1^2
  + 2 u0 u1 u2 + 2 u0 u1 u3 + u0 u2^2 + 2 u0 u2 u4
  + u1^2 u3 + u1 u3^2 + 2 u1 u3 u4
  + u2^2 u4 + 2 u2 u3 u4 + u2 u4^2
  + u3^2 u4 + u3 u4^2.
```

These are coefficientwise identities obtained from integer incidence
matrices, not fitted polynomial models.

## The desired six points really are stationary

The coordinates of `Box_i`, for `i<5`, are `5/6` in position `i` and
`-1/6` elsewhere; `Box_5` has all five coordinates `-1/6`.  Direct exact
substitution gives, for every `i`,

```text
q(Box_i)=7200,
f(Box_i)=14400,
grad f(Box_i)=3 grad q(Box_i).
```

**DERIVED:** all six fibration operators are stationary points of the cubic
on the fixed sphere.  Equality of their values is therefore not an accidental
failure of the stationarity condition.

## An exact additional solution

Consider the rational projective line

```text
u(t)=(-1, 2t-1, 0, 0, 0).
```

On it,

```text
q(t)=17280 (4t^2-6t+3),
f(t)=311040 (t-1)(2t-1).
```

After the positive normalization `s=sqrt(7200/q(t))`, the requirement
`f(sX)=+14400` is equivalent to

```text
1800 f(t)^2 = q(t)^3
```

with `f(t)>0`.  The residual factors exactly as

```text
1800 f(t)^2-q(t)^3 = -1289945088000 P(t),

P(t) = 256t^6 - 1152t^5 + 1764t^4 - 972t^3
       - 27t^2 + 162t - 27.
```

Exact rational endpoint values are

```text
P(1/10) = -185576/15625,
P(1/4)  = 157/64.
```

A Sturm count gives exactly one real root `alpha` in `(1/10,1/4)`; for
orientation only, `alpha` is approximately `0.212498536026`.  The polynomial
is square-free.  Throughout this interval `q(t)>0` and `f(t)>0`.  Consequently

```text
X_extra = sqrt(7200/q(alpha)) X(u(alpha))
```

satisfies exactly

```text
q(X_extra)=7200,
f(X_extra)=14400.
```

It is not a `Box_i`: its last three `E` coordinates vanish, whereas every
`Box_i` coordinate vector has nonzero entries there.

## It is a continuum, not one stray algebraic point

The root of `P` is simple, so the projective line meets the homogeneous level
hypersurface transversely.  If `df` were proportional to `dq` at this point,
Euler homogeneity for degrees three and two would force

```text
d(1800 f^2-q^3)=0
```

in every direction, contradicting the nonzero line derivative.  Thus `dq`
and `df` are independent at the normalized point.  The real implicit-function
theorem then gives a local common level set of dimension

```text
dim_R W - 2 = 3.
```

**DERIVED NEGATIVE:** the fixed quadratic and cubic constraints contain a
local real three-dimensional continuum in addition to the six desired
vertices.  They do not select the Hopf fibrations.

This is a mathematical continuum inside a five-dimensional real coefficient
space.  It is not evidence that physical spacetime, the carrier or the number
of degrees of freedom is infinite; the underlying matrices remain `120 x 120`.

## Protocol and scope audit

The preregistration requested a complete critical-scheme enumeration if the
six vertices were stationary, and independently declared that one exact extra
level-set point was sufficient to kill selection.  The extra point reaches
that explicit kill boundary, so the full critical scheme was not enumerated:
it cannot restore uniqueness.  Its detailed orbit stratification remains
**OPEN** but is no longer load-bearing for this route.

The negative is scoped precisely.  It closes selection by the pair
`(Tr(X^2),Tr(X^3))` on the derived Hopf--Box space.  It does not refute the
separate equal-weight projector cubic, which was proved to select the six
projector vertices.  What remains missing is a derivation of that particular
cubic from a licensed theory operator or action.

## Status ledger

- **DERIVED:** the six `Box_i` are stationary at the common levels
  `(q,f)=(7200,14400)` with Lagrange multiplier `3`.
- **DERIVED:** an exact extra algebraic point occurs on the rational line
  above; Sturm counting supplies existence, uniqueness in its isolating
  interval and simplicity.
- **DERIVED NEGATIVE:** a local real three-dimensional continuum shares the
  two levels, so the existing cubic is not a six-fibration selector.
- **STRUCTURAL:** treating arbitrary real combinations in `W` as physical
  fields; the negative is stronger because it defeats selection even under
  that permissive reading.
- **OPEN:** a licensed construction whose action produces the distinct
  equal-weight projector cubic with a derived coefficient and sign.
- **OPEN, non-load-bearing:** the complete critical-orbit decomposition of
  the operator cubic on the sphere.

## Consequence

The cheapest internal shortcut is closed.  The fact that all six desired
operators have the same cubic moment is not a selection theorem.  Progress
now requires an additional, independently derived structure that distinguishes
the projector cubic from `Tr(X^3)`; appending it because it has the desired
vacua would be fitting.
