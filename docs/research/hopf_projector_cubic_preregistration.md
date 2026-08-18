# Preregistration: projector-order-parameter cubic audit

Date: 2026-08-10

This protocol is registered before writing or running the exhaustive verifier.
It is **not blind**: the centered-projector idea and the expected regular
simplex relation were recognized algebraically before this commit.  The test
is therefore a derivation/correction audit, not independent evidence that the
six-vacuum target was discovered without looking.

## Fixed input

Use exactly the six rank-one projectors `P_i` onto the fivefold axes already
certified by `verify_hopf_symmetry_selector.py`.  No alternative axes,
weights, or normalization may be searched.

Define

```text
T_i = P_i - I_3/3  in Sym^2_0(R^3),
C3(Q) = sum_i Tr(Q T_i)^3.
```

Use the Frobenius inner product.  Coordinates and all comparisons remain in
`Q(sqrt(5))` or exact rational arithmetic.

## Claims to try to falsify

1. `Sym^2_0(R^3)` has real dimension five and the six `T_i` span it.
2. The `T_i` have equal norm and constant normalized cross inner product
   `-1/5`, hence form a regular 5-simplex centered at zero.
3. Their frame operator on `Sym^2_0(R^3)` is `(4/5) I_5` in an orthonormal
   basis.
4. On the fixed Frobenius sphere, `C3` has exactly six global maxima at the
   positive simplex vertices and six global minima at their negatives, after
   normalization.  This must be proved by an exhaustive stationary-point or
   constrained-coordinate argument, not sampling.
5. For `Q(n)=n n^T-(n.n)I/3`, the pullback is exactly

   ```text
   C3(Q(n)) = S6(n) - (34/45)(n.n)^3.
   ```

   Thus it has the same angular anisotropy as the previously certified `S6`.

## Falsifiers

The route fails if the centered projectors do not span dimension five, are
not a regular simplex, the cubic has additional or continuous global extrema,
or the exact pullback identity fails.

## Interpretation boundary

Even if all claims pass, the result establishes only a canonical five-real-
dimensional order-parameter geometry and a cubic invariant.  It does **not**
establish that `A1=5` denotes this space, that a certified inner fluctuation
contains this `Q`, that the spectral action generates `C3`, or that its sign
is fixed.  Those remain separate gates.
