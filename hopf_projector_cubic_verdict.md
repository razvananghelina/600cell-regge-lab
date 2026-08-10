# The Hopf choice as a five-dimensional projector order parameter

Date: 2026-08-10

Protocol commit: `402de35` (the expected simplex idea was declared there and
was not blind).  Registered verifier:
`reproducible/verify_hopf_projector_cubic.py`.  Targeted result: `9/9`.

## Exact result

For the six already certified fivefold-axis projectors `P_i`, define

```text
T_i = P_i - I/3  in Sym^2_0(R^3).
```

Exact arithmetic gives

```text
sum_i T_i = 0,
Tr(T_i^2) = 2/3,
Tr(T_i T_j) = -2/15  for i != j,
rank span{T_i} = 5,
sum_i |T_i><T_i| = (4/5) I_5.
```

Thus the six centered Hopf projectors form a regular 5-simplex in the
five-real-dimensional space `Sym^2_0(R^3)`.  This is **DERIVED**, not a
dimension fit.

The transitive six-axis permutation representation already certified in the
symmetry audit has one invariant line.  Its zero-sum complement is exactly
the five-dimensional space realized by the `T_i`.  In `A5` language this is
the nontrivial five-dimensional irreducible component of `1+5`.

## The cubic selector

The equal-weight cubic

```text
C3(Q) = sum_i Tr(Q T_i)^3
```

is canonical once the six-axis orbit and Frobenius metric are fixed.  On the
sphere `Tr(Q^2)=2/3`, write `s_i=Tr(QT_i)`.  The tight-frame identities give

```text
sum_i s_i = 0,
sum_i s_i^2 = 8/15.
```

At a constrained stationary point of `sum_i s_i^3`, every `s_i` solves the
same quadratic Lagrange equation, so there are at most two distinct values.
Enumerating their multiplicity `k=1,...,5` is exhaustive.  The exact classes
are

| multiplicity `k` | number | `C3` |
|---:|---:|---:|
| 1 | 6 | `64/225` |
| 2 | 15 | `8 sqrt(10)/225` |
| 3 | 20 | `0` |
| 4 | 15 | `-8 sqrt(10)/225` |
| 5 | 6 | `-64/225` |

Therefore the only global maxima are the six `T_i`, and the only global
minima are the six `-T_i`.  No sampling or fitted coefficients enter.

## Exact relation to the old sixth-order form

For an unconstrained vector `n`, set

```text
Q(n) = n n^T - (n.n) I/3.
```

The verifier proves

```text
C3(Q(n)) = S6(n) - (34/45)(n.n)^3.
```

The difference is radial.  Hence the previously derived degree-six angular
selector is exactly the pullback of a cubic on the natural unoriented-axis
order parameter.

**DERIVED CORRECTION:** “the selector first occurs at degree six” is true in
the vector `n`; it is false as an unqualified statement.  In the projector
variable `Q` it occurs at degree three.

## Relation to `a1=5`

The repository's `a1=5` is defined as the unique integer level selected by
the Fibonacci/`SU(2)` quantum-dimension bootstrap.  It is not defined as an
order-parameter dimension.

Starting from that level, the established chain selects `2I`, whose six C10
axes produce the canonical space

```text
span{P_i-I/3} = Sym^2_0(R^3),  dimension = 5 = a1.
```

The equality and the construction of the five-dimensional space are
**DERIVED**.  Identifying the semantic meaning of the bootstrap level `a1`
with the number of projector-order-parameter components is still
**STRUCTURAL**: no functor or physical axiom in the repository equates those
two roles.

## Consequence for the action gate

The earlier degree-ceiling argument assumed `D_A` was linear in the vector
`n`.  It correctly ruled out a degree-six term from a fourth moment under
that hypothesis.  It does not rule out `C3(Q)` if a certified fluctuation is
linear in `Q`: expansion of `Tr(D+A(Q))^4` can contain terms cubic in `Q`.

The current action verdict remains negative because the repository has not
constructed such an `A5`-equivariant five-real-dimensional one-form channel,
nor computed a nonzero cubic coefficient or its sign.  But `D^6` is no longer
the unique minimal continuation.  The cheaper exact gate is now:

1. find a licensed fluctuation subspace isomorphic to the centered six-axis
   `A5` module;
2. restrict the already available fourth spectral moment to it;
3. compute whether its cubic part is nonzero and whether its sign selects
   `T_i` rather than `-T_i`.

Appending `C3` by hand would still be fitting.  Deriving it from the
theory's operator would be a genuine advance.

The first canonical 120-vertex attempt is now closed negative in
`hopf_adjacency_baseline_cubic_verdict.md`: the invariant cubic space has
dimension two, while both `Tr(X^3)` and `4Tr(A X^3)` lie on the wrong line.
The full adjacency-baseline fourth moment also has a lower exact stationary
ten-point orbit.
