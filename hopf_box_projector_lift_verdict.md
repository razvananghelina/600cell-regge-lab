# The six Hopf `Box_F` operators realize the 5D simplex, but `Tr(Box^4)` does not select it

Date: 2026-08-10

Protocol commit: `818738a`.  Registered verifier:
`reproducible/verify_hopf_box_projector_lift.py`.  Targeted result: `11/11`.

## Fixed construction

For each of the six already certified discrete Hopf fibrations, let `A_f,i`
be its fibre-edge adjacency and use the derived wave operator

```text
Box_i = 6 A_f,i - A.
```

No weights were introduced.  All moment calculations after the combinatorial
fibrations were identified use integer matrices.

## Exact edge partition and zero centre

Each fibration contains 120 of the 720 undirected 600-cell edges.  The six
sets partition the edge set exactly:

```text
sum_i A_f,i = A.
```

It follows immediately that

```text
sum_i Box_i = 6 sum_i A_f,i - 6 A = 0.
```

Thus the affine centre `Box_bar=(1/6)sum_i Box_i` is exactly zero.  This is a
combinatorial identity, not a cancellation within numerical tolerance.

## Exact operator simplex

The Hilbert--Schmidt Gram matrix is

```text
Tr(Box_i^2)     =  7200,
Tr(Box_i Box_j) = -1440  for i != j,
rank Gram       = 5.
```

The normalized cross inner product is again `-1/5`.  Therefore the six
theory-defined `Box_i` form a centered regular 5-simplex in their operator
span.  The unique simplex map

```text
T_i=P_i-I/3  ->  Box_i
```

has fixed squared scale

```text
7200/(2/3) = (-1440)/(-2/15) = 10800.
```

**DERIVED:** the five-dimensional projector order parameter has a canonical
realization inside the theory's own six-fibration wave-operator family.  This
is substantially stronger than the numerical equality `dim=5=a1`, although
it still does not identify this span with a licensed inner-one-form space.

## Fourth-moment kill

Write the canonical affine interpolation as

```text
Box(Q) = Box_bar + X(Q) = X(Q).
```

The homogeneous cubic part of the fourth moment is normally
`4 Tr(Box_bar X^3)`.  Here it vanishes coefficient by coefficient because
`Box_bar=0`.  In fact

```text
Tr(Box(Q)^4) = Tr(X(Q)^4) = Tr(Box(-Q)^4).
```

At all six vertices,

```text
Tr(Box_i^3) =  14400,
Tr(Box_i^4) = 756000.
```

The fourth value is identical at `-Box_i`.

**DERIVED NEGATIVE:** the canonical `Box` realization does not rescue the
fourth-moment action gate.  Its fourth moment cannot distinguish the six
positive fibration vertices from the six negative, non-projector vertices.
There is no sign to read: the cubic coefficient is exactly zero.

This reaches the preregistered kill condition for the `Tr(Box^4)` shortcut.

## The third moment is a different open route

On the five-dimensional span, `Tr(X^3)` is nonzero.  The complete exact
polynomial comparison shows that it is not proportional to the equal-weight
projector cubic

```text
C3(Q)=sum_i Tr(QT_i)^3.
```

Therefore the existing identity `Tr(Box_i^3)=N^2` cannot be silently
relabelled as the six-axis selector.  A five-dimensional `A5` module admits
more than the single cubic form used in the projector construction.

**OPEN:** whether the constraint or extremization of `Tr(X^3)`, together with
all already derived normalization conditions, selects exactly the six
`Box_i`.  Equality at the desired vertices is not enough; the full solution
or critical set must be exhausted before making that claim.

## Revised next boundary

The fourth moment on the canonical `Box` interpolation is closed negative.
The cheapest remaining internal question is now sharply posed:

```text
On Tr(X^2)=7200 in span{Box_i}, what is the complete exact critical set of
Tr(X^3), and does Tr(X^3)=14400 select only the six Box_i?
```

If additional points or continua satisfy the same conditions, the existing
variational cubic is not a selector.  If the six vertices are uniquely
selected by the already derived constraints, the route advances without
adding a coupling by hand.
