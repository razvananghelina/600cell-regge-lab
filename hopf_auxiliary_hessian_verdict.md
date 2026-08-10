# The auxiliary selector is a Gram recognizer, not the Hessian of the original bootstrap

Date: 2026-08-10

Protocol commit: `90627f7`.  Registered verifier:
`reproducible/verify_hopf_auxiliary_hessian.py`.  Targeted result: `15/15`.

## Verdict

- **DERIVED NEGATIVE for direct dynamics:** the actual one-parameter
  variational bootstrap has identically zero Hessian and cannot generate the
  twelve-dimensional auxiliary selector.
- **DERIVED CONDITIONAL bridge:** after extending the cubic to the full
  Hopf--Box field space, the diagonal of its six-label Hessian is exactly the
  overlap operator `Phi`.
- **DERIVED OBSTRUCTION:** the complete Hessian contains nonzero off-diagonal
  information of rank five and cannot be simultaneously diagonalized.  Keeping
  only its diagonal is a new fibration-label locality/superselection
  hypothesis.
- **DERIVED:** `D_aux` is exactly a normalized pair-Gram operator.  Its Schur
  complement detects Cauchy equality, not an equation of motion inherited
  from the old bootstrap.

Thus the prior six-point conditional selection remains mathematically exact,
but kernel saturation has not been derived dynamically.

## The actual bootstrap Hessian is zero

For every one of the six fibrations, the original operator family is

```text
B_i(c)=cA_f,i-A.
```

The complete integer trace polynomial is

```text
Tr(B_i(c)^3)=3600(c-2).
```

Equivalently, its coefficient vector in descending powers is

```text
(alpha_3,alpha_2,alpha_1,alpha_0)=(0,0,3600,-7200).
```

Therefore

```text
d^2/dc^2 Tr(B_i(c)^3)=0
```

identically for all six fibrations.

**DERIVED NEGATIVE:** neither the Hessian nor a quadratic stability operator
of the actual one-dimensional bootstrap can be `D_aux`.  The linearity that
made `c=6` uniquely solvable also removes precisely the second variation that
this continuation hoped to use.

## The extended cubic has the disclosed diagonal

On the already derived field space `W`, let

```text
V(X)=Tr(X^3),
H_X(i,j)=3Tr(X(Box_i Box_j+Box_j Box_i)).
```

This is the full Hessian in the six distinguished synthesis directions.  The
constant label vector is always in its kernel because `sum_i Box_i=0`.

Exact coefficient comparison gives

```text
H_X(i,i)=12Tr(X Box_i),
diag(H_X)=86400 Phi(X).
```

Hence

```text
Phi(X)=diag(H_X)/(12*7200).
```

**DERIVED CONDITIONAL:** the auxiliary overlap field is exactly the normalized
diagonal response of the extended cubic action.  This is a real bridge and
uses no fitted coefficient.

It is not yet the full Hessian response.

## The discarded sector contains all five field directions

The off-diagonal map

```text
X -> {H_X(i,j):i<j}
```

has exact rank five, the same as `W`.  Already for one basis direction its
nonzero values include

```text
-155520,-103680,-51840,51840,103680,155520.
```

At `X=Box_0`, the complete Hessian is

```text
[ 86400 -17280 -17280 -17280 -17280 -17280]
[-17280 -17280 -69120 -69120  86400  86400]
[-17280 -69120 -17280  86400 -69120  86400]
[-17280 -69120  86400 -17280  86400 -69120]
[-17280  86400 -69120  86400 -17280 -69120]
[-17280  86400  86400 -69120 -69120 -17280].
```

Its spectrum is

```text
0                                      multiplicity 1
103680                                 multiplicity 1
-25920-77760sqrt(5)                    multiplicity 2
-25920+77760sqrt(5)                    multiplicity 2.
```

The universal zero is only the redundant constant synthesis direction.  On
the physical zero-sum sector the signature is `(3 positive,2 negative)` and
there is no zero.  Thus the full Hessian's kernel condition is either automatic
and content-free on `R^6`, or absent after removing the redundancy.

Moreover,

```text
rank[H_Box0,H_Box1]=4,
max absolute commutator entry=26873856000.
```

Two real symmetric matrices admit a common orthogonal eigenbasis only if they
commute.  Therefore no fixed change of label basis converts the full Hessian
family into the diagonal `Phi` family.

**DERIVED OBSTRUCTION:** the diagonal is not the full Hessian written in a
better basis.  The omitted sector is rank five and load-bearing.

## What remains canonical

Both the full Hessian and the diagonal conditional expectation

```text
E_diag:H -> diag(H_00,...,H_55)
```

are exactly `A5`-equivariant under permutation of fibration labels.  Therefore
the diagonal projection requires no weights or preferred fibration.

But covariance is weaker than derivation.  Selecting `E_diag(H_X)` asserts
that different global fibration labels do not mix.  The original vertex
bootstrap contains no such locality or superselection axiom.

Status:

- **DERIVED:** the conditional expectation is canonical once the commutative
  label algebra is adopted;
- **STRUCTURAL:** adopting label locality and discarding the full off-diagonal
  Hessian.

## Exact Gram and Schur interpretation

For normalized `X` and each `Box_i`, the Hilbert--Schmidt Gram matrix is

```text
G_i(X)=[1    r_i(X)]
       [r_i(X)  1  ].
```

Consequently

```text
D_aux(X)=direct_sum_i G_i(X),
Schur(G_i)=1-r_i(X)^2,
K(X)=I-Phi(X)^2.
```

This fully derives positivity and the twelve signed singular points.  It also
clarifies their meaning:

```text
det G_i(X)=0  iff X and Box_i are linearly dependent.
```

Every normalized pair of vectors has exactly this Gram property.  Therefore
the singularity recognizes that `X` is a signed frame vertex; it does not
explain why a dynamical field must become collinear with one.

Calling `K` “a Schur complement derived from `Box`” without this distinction
would be circular.  It is the Schur complement of the newly constructed Gram
operator, not of the original `120 x 120` wave operator.

## Scoped source audit

The authoritative original constructors

```text
reproducible/verify_variational_bootstrap.py
reproducible/verify_hopf_fibration_invariants.py
```

contain only the vertex adjacency, fibre adjacency and their `120 x 120`
operators.  They define no state-space map or block coupling from `C^120` to
the six-label carrier.

This is an absence result for the specified construction, not a theorem that
no future coupling can exist.

## Status ledger

- **DERIVED NEGATIVE:** original bootstrap Hessian `=0` for all six
  fibrations.
- **DERIVED CONDITIONAL:** `diag(H_X)/86400=Phi(X)` for the extended cubic.
- **DERIVED:** full and off-diagonal Hessian maps both have rank five.
- **DERIVED:** the Hessian family is not simultaneously diagonalizable.
- **DERIVED:** diagonal expectation remains `A5`-equivariant.
- **STRUCTURAL:** fibration-label locality/superselection.
- **DERIVED:** `D_aux` is the normalized direct sum of pair-Gram matrices.
- **DERIVED NEGATIVE:** Gram singularity is a collinearity recognizer, not a
  consequence of the original equation of motion.
- **OPEN:** an actual coupling making fibration labels dynamical and local.
- **OPEN:** a stability or integration mechanism that suppresses the
  off-diagonal Hessian rather than deleting it.

## Consequence and next boundary

The direct Hessian/Schur route is closed negative.  The useful surviving fact
is narrower: `Phi` is the canonical *local diagonal response* of the extended
cubic.

This continuation is now closed in `hopf_label_superselection_verdict.md`.
The minimal `C^6` representation enforces diagonality only by killing all
one-forms and connectedness.  A canonical pair bimodule passes order zero,
first order, connectedness, KO6 and `A5` with 360 off-diagonal channels.
Adding metric-zero orientability kills every nonzero `A5`-equivariant KO6
`C^6` bimodule, even with arbitrary multiplicity.  Thus the commutative label
algebra does not license the selector under the current axioms.
