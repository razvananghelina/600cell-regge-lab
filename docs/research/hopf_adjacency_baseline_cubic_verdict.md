# The canonical adjacency baseline does not generate the Hopf selector

Date: 2026-08-10

Protocol commit: `48bbb00`.  Registered verifier:
`reproducible/verify_hopf_adjacency_baseline_cubic.py`.  Targeted result:
`15/15`.

## Fixed hypotheses

The tested operator family was frozen as

```text
W = span_R{Box_i},                       dim_R W=5,
Box_i = 6 A_f,i-A,                      sum_i Box_i=0,
D_A(X)=A+X,
q(X)=Tr(X^2)=7200.
```

The baseline `A` is the already derived 600-cell vertex adjacency.  No
fibration weights, matrix insertions or fitted linear combinations of trace
words were allowed.  At a desired vertex, `A+Box_i=6A_f,i`.

This family is canonical on the 120-vertex carrier but is not already a
licensed fluctuated finite spectral-triple Dirac operator.  The polynomial
test therefore could at most have produced a **STRUCTURAL** action bridge.

## Symmetry permits exactly two cubic tensors

Quaternionic conjugation reconstructs all 60 rotations and their permutation
action on the six Hopf fibrations.  Its exact cycle distribution is

```text
cycle type             number
1^6                         1
1^2 2^2                    15
3^2                        20
1 5                        24
```

The five-dimensional zero-sum module `W` consequently has character
distribution

```text
chi_W = 5, 1, -1, 0
count = 1,15, 20,24.
```

For every group element,

```text
chi_Sym3(g) = (chi(g)^3 + 3 chi(g)chi(g^2) + 2 chi(g^3))/6.
```

The exact group average is `120/60=2`.  Hence

```text
dim Sym^3(W*)^A5 = 2.
```

**DERIVED:** icosahedral symmetry does not uniquely determine a cubic on this
order-parameter space.  There are exactly two invariant directions.  Thus
calling a cubic “the `A5` invariant” would hide one real coefficient of
freedom.

## The adjacency baseline remains on the wrong invariant line

The zero-baseline cubic and the homogeneous cubic part of the fixed fourth
moment are

```text
G_0(X) = Tr(X^3),
G_A(X) = 4 Tr(A X^3).
```

Exact coefficient comparison in a basis of `W` gives the stronger identity

```text
G_A(X) = -8 G_0(X)
```

for every `X in W`.  Adding `A` therefore changes the coefficient and sign
but does not open the second invariant direction.

The canonical operator-simplex cubic

```text
C_box(X) = sum_i Tr(X Box_i)^3
```

is independent of that line.  The coefficient ranks are

```text
rank(G_0,G_A)       = 1,
rank(G_0,G_A,C_box) = 2.
```

It has value `358318080000` on each `Box_i` and is the operator-scaled version
of the previously derived equal-weight projector selector.

**DERIVED NEGATIVE:** the cubic part of `Tr((A+X)^4)` is not the Hopf
selector.  The complete two-dimensional invariant space was accounted for,
so this is not a failed numerical proportionality test inside an unknown
larger space.

## The complete fourth moment also fails dynamically

The full fixed functional is

```text
S_4(X)=Tr((A+X)^4),             q(X)=7200.
```

All twelve signed simplex vertices are exact constrained stationary points:

```text
point       multiplier      S_4
+Box_i          216         933120
-Box_i          264        1163520
```

Thus a positive fourth-moment coefficient does have the desired sign when
only `+Box_i` and `-Box_i` are compared.  Reading that pairwise comparison as
selection would nevertheless be wrong.

Take the exact direction in `E_a=Box_a-Box_5` coordinates

```text
v=(-1,1,1,-1,1),
Tr(X(v)^2)=51840,
s=sqrt(5)/6.
```

Then `X_*=sX(v)` lies on the required sphere and

```text
S_4(X_*) = 1048320 - 115200 sqrt(5),

S_4(X_*)-S_4(Box_i) = 115200(1-sqrt(5)) < 0.
```

It is itself an exact stationary point, with multiplier
`240-24sqrt(5)`.  In zero-sum six-vertex coordinates its `A5` orbit has size
ten and stabilizer order six; its negative lies in the other ten-point orbit.

**DERIVED NEGATIVE:** none of the six desired vertices is a global minimum of
the complete fixed fourth moment.  The failure is not caused by dropping its
quadratic or quartic pieces.

The ten-point orbit has the same group-theoretic cardinality and stabilizer
order as the ten unoriented icosahedral face axes.  Identifying the two
geometrically without an explicit intertwiner is recorded only as a
**PATTERN**, not used in the negative proof.

## Framing correction

“Canonical” and “`A5` invariant” are insufficient here.  The geometry admits
two invariant cubic lines.  The standard single-trace construction with
baselines `0` and `A` samples only one of them, whereas the six-axis selector
is the other.  Moving between them requires an additional operation, not a
normalization choice.

The equal-weight multi-trace `C_box` is itself canonical as a polynomial once
the unordered six-fibration set is given.  What is absent is a reason for the
theory's action to contain this multi-trace rather than the single trace.  To
append it now because it selects the desired orbit would be fitting.

## Status ledger

- **DERIVED:** `dim Sym^3(W*)^A5=2`.
- **DERIVED:** `4Tr(AX^3)=-8Tr(X^3)` on all of `W`.
- **DERIVED:** `C_box` spans the second invariant cubic line and retains the
  exact six simplex extrema.
- **DERIVED NEGATIVE:** the adjacency-baseline fourth-moment cubic is not
  proportional to `C_box`.
- **DERIVED NEGATIVE:** the complete `S_4` has an exact lower stationary
  ten-point orbit, so it does not dynamically select the six fibrations.
- **PATTERN:** the lower orbit has the face-axis orbit type.
- **STRUCTURAL:** regarding `A+X` as a fluctuated Dirac operator.
- **OPEN:** whether a higher already motivated single-trace moment reaches
  the second invariant line, and whether its complete fixed action has the
  correct extrema.
- **OPEN:** a licensed spectral triple realizing `W` as inner one-forms.

## Next finite test

The next non-fitted calculation is the cubic-in-`X` coefficient of the fixed
higher moments

```text
Tr((A+epsilon X)^p),
```

especially the even sixth moment already singled out by the heat-trace
expansion.  Because the entire invariant cubic space is two-dimensional, an
exact rank test will decide whether these higher single traces ever reach
`C_box`.  If they remain on the `Tr(X^3)` line, the whole adjacency-polynomial
single-trace route closes; if the sixth moment reaches the selector line, its
full sign and competing extrema must then be tested.
