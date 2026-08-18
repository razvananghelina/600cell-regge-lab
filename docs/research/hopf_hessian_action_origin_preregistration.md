# Preregistration: can the certified spectral action generate the Hessian selector?

Date: 2026-08-11

## Complete question and hypotheses

Let `Hhat_X` be the complete full-label Hessian restricted to the canonical
five-dimensional subspace `1^perp`.  The preceding blind/target protocol has
already proved

```text
Tr(Hhat_X)=0,
Tr(Hhat_X^2)=constant on q(X)=7200,
Tr(Hhat_X^3)=-23328 C_box(X),
```

and that its positive cubic moment has global minima exactly at `+Box_i`.

This audit separates two claims:

1. **existing-action origin:** the already certified Kähler--Dirac spectral
   action contains this label-Hessian block and its cubic response;
2. **minimal structural extension:** the uniquely equivariant affine baseline
   on the five-dimensional label module produces the same selector in a
   complete even fourth moment, without tuning a moment combination.

No claim may move from (2) to (1) merely because the matrix construction is
canonical.

## Frozen test A: current certified action

Audit the authoritative constructors and record exactly:

- the carrier and operator used by `verify_spectral_action.py`;
- whether a map `X -> Hhat_X` or any six-label block is present;
- whether the certified functional uses `D`, `D^2`, or both;
- the exact moments presently certified.

For any functional `Tr f(Hhat_X^2)`, prove coefficientwise that its Taylor
series is even under `X -> -X` and therefore contains no
`Tr(Hhat_X^3)`.  For any grading-odd Dirac operator, independently verify that
all odd full traces vanish by spectral symmetry.

**Kill for existing-action origin:** the current action has no `Hhat_X`
coupling, and its even/graded parity forbids the required cubic if `Hhat_X` is
inserted without a nonzero baseline.

## Frozen test B: exhaustive equivariant affine baseline

The physical label module is the irreducible real `A5` module `W_5`.  Compute

```text
dim End_A5(W_5).
```

If it is one, every equivariant constant baseline is `b I`.  This makes the
affine family exhaustive:

```text
B_(b,c)(X)=b I_5+c Hhat_X.
```

Place it in the minimal grading-odd self-adjoint double

```text
D_(b,c)(X) = [[0,B_(b,c)(X)],
              [B_(b,c)(X),0]].
```

Compute the **complete** fourth moment `S4=Tr(D_(b,c)^4)`, not only its cubic
coefficient.  On `q=7200`, classify its exact global minima for all parameter
branches:

```text
b*c^3 > 0,
b*c^3 < 0,
b*c^3 = 0.
```

No value of `|b/c|` may be selected after target comparison.  Report whether
the extremal locus is independent of every positive magnitude ratio and give
the sign look-elsewhere fraction.

## Physical gates for the minimal extension

Record before promotion:

1. whether this doubled operator already occurs in the repository;
2. whether an algebra representation, `J`, order zero, first order,
   orientability and nonzero inner one-forms have been constructed for it;
3. whether the relative sign of `b` and `c` is fixed by an existing axiom;
4. whether `Hhat_X` is a licensed fluctuation rather than a Hessian covariant
   adopted as a new field.

Passing only self-adjointness, grading and `A5` covariance is
**STRUCTURAL**, not a finite spectral triple or a physical action.

## Decision boundary

- **Physical action advance:** the already certified action contains the
  block/coupling and fixes the selector coefficient and sign.
- **Structural action advance:** the exhaustive equivariant affine double has
  a complete fourth moment whose global minima are exactly one signed Hopf
  orbit for an entire sign branch, with no magnitude tuning.
- **Pattern:** only a selected ratio, truncation or post-comparison sign works.
- **Kill:** parity removes the cubic, or the complete fourth moment has lower
  exact competitors for both sign branches.

## Provenance limitation

The identity `Tr(Hhat_X^3)=-23328 C_box` is already known, so this protocol is
not blind to the desired cubic coefficient.  Its falsifiable new content is
whether the complete fourth moment preserves the target for all magnitude
ratios, and whether the existing repository action actually contains the
required operator.  Those outcomes are frozen here before calculation.
