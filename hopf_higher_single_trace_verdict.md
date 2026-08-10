# Higher single traces do not derive the Hopf selector

Date: 2026-08-10

Protocol commits: `aea13cb`, followed by the pre-full-functional framing
correction `35566e8`.  Registered verifier:
`reproducible/verify_hopf_higher_single_trace.py`.  Targeted result: `23/23`.

## Complete hypotheses

Use only the fixed 120-vertex family

```text
D(epsilon,X)=A+epsilon X,
X in W=span_R{Box_i},
Tr(X^2)=7200.
```

No fibration weights, auxiliary projectors or fitted combinations of trace
words are admitted.  For a moment `p`, its action-selected cubic is

```text
K_p(X)=[epsilon^3]Tr((A+epsilon X)^p).
```

The sixth moment is privileged because it is the next even term in the formal
finite heat-trace Taylor series.  The repository still has no licensed
fluctuated Dirac operator of this form and no fixed heat scale.

## Why the complete word search is finite

The adjacency satisfies exactly

```text
m_A(z)=z(z-12)(z-3)(z+2)(z+3)
       (z^2-6z-36)(z^2-4z-16).
```

Integer Horner evaluation gives `m_A(A)=0` coefficientwise.  An exact nonzero
`9 x 9` minor of the flattened matrices `I,A,...,A^8` proves that no lower
degree polynomial annihilates `A`.  Hence this is the minimal polynomial.

Every single-trace cubic word therefore reduces to

```text
T_abc(X)=Tr(A^a X A^b X A^c X),
0 <= a <= b <= c <= 8.
```

Cyclicity and transpose reversal generate every permutation of `(a,b,c)`.
The exact number of reduced individual words is

```text
N=binomial(11,3)=165.
```

## Exact enumeration and look-elsewhere count

The full invariant cubic space was independently checked to have dimension
two.  The basis used for comparison was

```text
G_0(X)=Tr(X^3),
C_box(X)=sum_i Tr(X Box_i)^3.
```

Two fixed evaluation points give the exact matrix

```text
                         G_0          C_box
Box_0                   14400      358318080000
X(-1,1,1,-1,1)         622080                 0
```

whose determinant is `-222902511206400000`, so these two evaluations are
coordinates on the complete invariant space.

All 330 word values were reconstructed by CRT using six declared primes.  The
modulus is

```text
1000292032458727685153601621373570283.
```

For both evaluation points it exceeds twice the rigorous norm bound

```text
120 * 12^24 * ||X||_infinity^3,
```

so signed reconstruction is unique.  Every value also passes an unused-prime
check.

The complete result is

```text
individual words                       165
distinct projective cubic lines         24
span rank                                2
words proportional to C_box              0
hit fraction                          0/165
words on the old G_0 line                 2
old-line exponent triples       (0,0,0), (0,0,1)
```

**DERIVED NEGATIVE:** no individual adjacency-polynomial single-trace word
is the Hopf selector.  This is a complete finite enumeration, not failure to
find one at low exponent.

The 165 words collectively span both invariant directions, so a linear
combination can manufacture `C_box`.  That does not advance the theory:
choosing coefficients after target comparison is precisely fitted
Schur/invariant data.

## The action-selected moment cubics

The exact noncommutative word sum gives, in the basis `(G_0,C_box)`,

```text
K_3 =      G_0,
K_4 =   -8 G_0,
K_5 =   40 G_0 + C_box/1244160,
K_6 = -280 G_0 + C_box/124416.
```

Thus the sixth moment opens the second invariant direction but does not isolate
it.  Its cubic is a fixed mixture.  Cancelling `G_0` would require an
additional moment coefficient.

After the comparison one notices

```text
K_6-35K_4=C_box/124416.
```

The numerical equality `35=dim Sym^3(R^5)` is exact, but no principle in the
repository derives this subtraction.  It is therefore **PATTERN**, disclosed
after target comparison, not evidence for an action.

The same cancellation occurs in the order-six truncated heat series at
`t=3/35`, whose angular part is proportional to

```text
S_4-S_6/35.
```

This value of `t` was obtained by cancelling the unwanted target component;
it was not independently derived.

## The complete sixth moment fails dynamically

Purity of the cubic is not logically necessary if the complete fixed action
selects the desired points.  This was corrected in commit `35566e8` before
the full sixth moment was evaluated.

For

```text
S_6(X)=Tr((A+X)^6),
```

all signed simplex vertices are stationary:

```text
point       multiplier         S_6
+Box_i        38880         111974400
-Box_i        69360         200678400
```

Let `Y=X(-1,1,1,-1,1)` and `s=sqrt(5)/6`, so that
`Tr((sY)^2)=7200`.  Both signs are exact constrained stationary points, with

```text
S_6(+sY)=165542400-38592000sqrt(5),
S_6(-sY)=165542400+38592000sqrt(5).
```

The first is strictly below `S_6(+Box_i)`, so a positive sixth-moment action
does not select the six.  The second is strictly above it, so the formal heat
coefficient `-S_6` does not select them either.

**DERIVED NEGATIVE:** exact stationary competitors defeat the desired orbit
for both signs motivated by the sixth term.

The post-hoc truncated-heat cancellation fails even more cheaply.  Although
its cubic part has the desired `-C_box` sign,

```text
(S_4-S_6/35)(+Box_i) = -15863040/7,
(S_4-S_6/35)(-Box_i) = -31991040/7.
```

The negative vertices are lower by `2304000`.  Cancelling one cubic component
does not control the other odd-degree pieces of the complete action.

## Status ledger

- **DERIVED:** the adjacency minimal polynomial has degree nine, reducing the
  complete word space to 165 candidates.
- **DERIVED:** those words occupy 24 projective lines and span the full
  two-dimensional invariant cubic space.
- **DERIVED NEGATIVE:** selector hit fraction `0/165` for individual words.
- **DERIVED:** `K_6=-280G_0+C_box/124416`.
- **DERIVED NEGATIVE:** the full sixth moment and its formal heat sign both
  prefer exact stationary competitors over the six Hopf vertices.
- **PATTERN:** the post-comparison cancellation coefficient `35` and
  `t=3/35`; the complete truncated functional fails anyway.
- **STRUCTURAL:** treating `A+X` as a fluctuated Dirac operator.
- **OPEN:** a derivation of the canonical multi-trace `C_box` from integrating
  out a licensed operator, with regulator and coefficient fixed beforehand.
- **OPEN:** a valid `A5`-equivariant finite spectral triple whose inner
  one-forms realize `W`.

## Consequence

The canonical fourth and sixth single-trace moments are now closed negative.
The full adjacency-polynomial word algebra contains enough freedom to fit the
selector, but it does not select a word or combination.  Testing successively
higher moments after seeing each failure would create an unbounded
look-elsewhere search unless a moment or spectral function is fixed by an
independent axiom.

The next honest route is no longer “try another trace power.”  It is to derive
why the action should contain the already canonical multi-trace `C_box`, for
example from a specified auxiliary-field or integration mechanism, and to
preregister that mechanism before computing its induced coefficients.
