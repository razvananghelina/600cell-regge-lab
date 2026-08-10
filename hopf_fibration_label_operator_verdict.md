# A canonical auxiliary operator recognizes exactly the signed Hopf simplex

Date: 2026-08-10

Protocol commit: `8dcc164`.  Registered verifier:
`reproducible/verify_hopf_fibration_label_operator.py`.  Targeted result:
`17/17`.

## Verdict

There is a genuine **STRUCTURAL ADVANCE**:

1. the unordered set of six derived Hopf fibrations supplies a canonical
   six-dimensional permutation carrier;
2. the overlap map from the Hopf--Box order-parameter space into that carrier
   is the unique `A5` intertwiner up to scale;
3. its sharp norm fixes the scale;
4. a canonical sign-neutral positive operator has kernel exactly at the
   twelve points `+/-Box_i`;
5. intersecting this kernel locus with the already derived positive cubic
   condition leaves exactly the six actual `+Box_i`.

This is not yet a physical vacuum theorem.  The missing step is why a physical
field must saturate this auxiliary spectral bound, or equivalently why the
founding nontrivial-kernel rule applies to this new carrier.

## Canonical carrier and uniqueness

Let `F` be the six-element set of certified Hopf fibrations and set

```text
H_F=R^F.
```

Quaternionic conjugation reconstructs the full 60-element `A5` action.  Its
cycle census on `F` is

```text
1^6:       1
1^2 2^2: 15
3^2:      20
1 5:      24.
```

The permutation module decomposes as

```text
H_F = 1 + W,
dim_R W=5.
```

The exact character inner products are

```text
<W,W>=1,
<W,1>=0,
<W,H_F>=1.
```

Therefore

```text
dim Hom_A5(W,H_F)=1.
```

At matrix level every icosahedral rotation permutes the six `Box_i` exactly,
so the Hilbert--Schmidt overlap map is a nonzero intertwiner and hence spans
this unique line.

**DERIVED:** no Schur coefficients are hidden in the map.  Once the
six-fibration carrier is admitted, equivariance leaves only an overall scale.

## Sharp normalized analysis operator

On the fixed sphere `Tr(X^2)=7200`, define

```text
r_i(X)=Tr(X Box_i)/7200,
Phi(X)=diag(r_0(X),...,r_5(X)).
```

The denominator is the exact sharp Cauchy bound

```text
||X||_HS ||Box_i||_HS = sqrt(7200)sqrt(7200)=7200.
```

The tight-frame identities give

```text
sum_i r_i(X)=0,
sum_i r_i(X)^2=Tr(X^2)/6000=6/5.
```

The analysis map has rank five and maps `W` isomorphically onto the zero-sum
hyperplane in `R^6`.  Every coordinate has the sharp range

```text
-1 <= r_i(X) <= 1.
```

The exact equality cases are

```text
r_i(X)=+1  iff X=+Box_i,
r_i(X)=-1  iff X=-Box_i.
```

This is an exhaustive Cauchy equality classification, not numerical
optimization.

## The sign ambiguity and its exact count

For an affine diagonal slack

```text
D_{a,b}(X)=aI_6+bPhi(X),
```

positivity over the full sphere is equivalent to `a>=|b|`.  A nonzero sharp
operator acquires a kernel only when equality holds.  Up to positive scale,
there are exactly two candidates:

```text
I-Phi: kernel locus {+Box_i},
I+Phi: kernel locus {-Box_i}.
```

Thus the desired-set hit fraction is `1/2`.  Choosing the minus sign after
comparison would be a look-elsewhere choice.  Equivariance and positivity do
not by themselves distinguish the two signed simplex orbits.

## Sign-neutral operator and exact kernel

The unique sharp quadratic slack that is even under `Phi -> -Phi` is

```text
K(X)=I_6-Phi(X)^2 >= 0.
```

It is singular exactly when some `|r_i|=1`; hence

```text
ker K(X) != 0  iff X in {+/-Box_0,...,+/-Box_5}.
```

There are exactly twelve points and no additional curve or component.  At a
signed vertex the spectrum is

```text
{0, 24/25,24/25,24/25,24/25,24/25}.
```

Moreover, `K` is the Schur complement of the operator linear in `X`

```text
D_aux(X) = [ I    Phi(X) ]
           [ Phi(X)  I   ].
```

A fixed sheet-Hadamard transform diagonalizes it as

```text
(I+Phi) direct-sum (I-Phi).
```

Thus `D_aux>=0` on the whole sphere and has the same twelve-point kernel
locus.  At every signed vertex it has exactly one zero mode and spectrum

```text
{0,2,(4/5)x5,(6/5)x5}.
```

**DERIVED:** the sign-neutral selector can be represented by a canonical
positive `12 x 12` operator linear in the order parameter.  It is not merely
an arbitrarily chosen nonlinear polynomial.

## The missing multi-trace becomes one ordinary trace

On the label carrier,

```text
Tr_HF(Phi(X)^3)
    = sum_i (Tr(X Box_i)/7200)^3
    = C_box(X)/7200^3.
```

At the signed vertices its values are exactly

```text
Tr_HF(Phi(+Box_i)^3)=+24/25,
Tr_HF(Phi(-Box_i)^3)=-24/25.
```

**DERIVED:** the previously missing canonical multi-trace on the 120-vertex
carrier is an ordinary single trace on the geometry-selected six-label
carrier.  This resolves the algebraic representation problem without weights
or fitted coefficients.

It does not fix why the physical action contains this odd trace or its sign.

## Conditional six-point selection

The sign-neutral kernel condition first leaves all twelve signed vertices.
The founding positive cubic condition already used by the `Box` construction
is

```text
Tr(X^3)=N^2=14400.
```

Direct exact moments give

```text
Tr((+Box_i)^3)=+14400,
Tr((-Box_i)^3)=-14400.
```

Therefore

```text
Tr(X^2)=7200,
ker D_aux(X) != 0,
Tr(X^3)=14400
```

have exactly six solutions: `X=+Box_i`.

**DERIVED CONDITIONAL:** reusing both the sharp-kernel bootstrap and the old
positive cubic constraint selects the six actual Hopf operators exactly.  The
earlier three-dimensional continuum is removed.

## Hostile framing audit

This result must not be oversold.

1. `Phi` is built from the six already known `Box_i`.  It recognizes the
   frame vertices; it does not independently predict that the six fibrations
   exist.
2. A generic `X` on the sphere gives a strictly positive `D_aux`.  Geometry
   alone does not force its smallest eigenvalue to close.
3. The repository has used nontrivial kernels to select the coefficient
   `c=6` in the original `Box_F(c)`, but it has not stated a universal axiom
   that every auxiliary correlation operator must be singular.
4. Extending `Tr(Box^3)=N^2` from the original one-parameter wave operator to
   an arbitrary order-parameter field is likewise a conditional reuse of the
   founding bootstrap.
5. `H_F` is a commutative six-label carrier.  It supplies no non-abelian gauge
   sector and is not a constructed real finite spectral triple.

Thus “the theory dynamically derives six vacua” remains **OPEN**.  What is now
closed is the finite algebraic bridge that was missing: there is a unique
equivariant auxiliary operator whose sharp kernel, together with the old
positive cubic, has exactly those six solutions.

## Relation to `a1=5`

The carrier decomposition is now an exact representation-theoretic statement:

```text
R^6 = 1 + W_5,
dim W_5=5=a1.
```

The five-dimensional equality and construction are **DERIVED**.  The claim
that the bootstrap integer `a1` physically means the number of order-parameter
components remains **STRUCTURAL**.

## Status ledger

- **DERIVED:** `dim Hom_A5(W,H_F)=1`.
- **DERIVED:** the sharp normalized analysis is a tight-frame isomorphism.
- **DERIVED:** exactly two affine PSD sharp signs exist; desired hit `1/2`.
- **DERIVED:** the sign-neutral doubled operator is PSD and has exactly twelve
  signed kernel points.
- **DERIVED:** `Tr_HF(Phi^3)=C_box/7200^3`.
- **DERIVED CONDITIONAL:** the auxiliary kernel plus the old positive cubic
  leaves exactly six `+Box_i`.
- **STRUCTURAL:** imposing auxiliary kernel saturation as a physical axiom.
- **OPEN:** derive `D_aux` as a Hessian, Schur complement or fluctuation of an
  already licensed theory operator, rather than merely a canonical recognizer.
- **OPEN:** a physical action or stability principle fixing saturation and
  the positive sign.

## Next investigation

The next load-bearing question is no longer whether a canonical selector
operator exists—it does.  It is whether `D_aux` arises from the second
variation or Schur complement of the already derived `Box` variational
bootstrap.  A positive result would turn kernel saturation from a new axiom
into a consequence; a negative would delimit this construction as a useful
but auxiliary recognizer.
