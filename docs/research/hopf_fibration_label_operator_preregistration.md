# Preregistration: canonical six-fibration label operator

Date: 2026-08-10

## Provenance and anti-circularity warning

The idea is known before this protocol: the six derived `Box_i` form a tight
simplex frame, so their Hilbert--Schmidt overlaps can be placed on the
diagonal of the canonical six-point permutation module.  This is not a blind
discovery.

The construction uses the six fibrations as input.  Therefore even a perfect
kernel locus would be a selector/recognizer inside their span, not an
independent derivation that the six fibrations exist.  That limitation must
remain explicit.

## Fixed data

Use exactly

```text
W = span_R{Box_i},                 dim_R W=5,
q(X)=Tr(X^2)=7200,
H_F=R^{six fibrations}.
```

The group `A5` acts on `H_F` by permuting the six fibrations.  Define the
normalized analysis map

```text
r_i(X)=Tr(X Box_i)/7200,
Phi(X)=diag(r_0(X),...,r_5(X)).
```

The denominator is frozen as the sharp Cauchy--Schwarz bound on the sphere:
`||X|| ||Box_i||=7200`.  No coordinate weights or post-comparison rescaling
are allowed.

## Frozen exact tests

1. Reconstruct the 60-element `A5` action on the six fibrations and its
   zero-sum module `H_F^0`.
2. Compute exactly

   ```text
   dim Hom_A5(W,H_F)=1.
   ```

   Verify that the overlap map is a nonzero intertwiner and hence unique up
   to scale; verify that the sharp norm normalization above fixes that scale.
3. Prove coefficientwise that

   ```text
   sum_i r_i(X)=0,
   sum_i r_i(X)^2=6/5
   ```

   on `q(X)=7200`, and that the normalized analysis map is an isomorphism
   from `W` onto the zero-sum hyperplane.
4. Determine exhaustively the affine diagonal positive-semidefinite slack
   operators

   ```text
   D_{a,b}(X)=a I_6+b Phi(X)
   ```

   that are positive for every `X` on the sphere and acquire a kernel
   somewhere.  Count them up to positive overall scale.
5. For every surviving affine operator, compute its complete kernel locus.
   Report the hit fraction for `+Box_i`, `-Box_i` and their union.  Do not
   choose a sign after seeing which one hits the desired set.
6. Test the sign-neutral sharp slack

   ```text
   K(X)=I_6-Phi(X)^2.
   ```

   Prove positivity and determine its complete kernel locus.
7. Only after that, intersect the sign-neutral kernel locus with the already
   derived founding conditions

   ```text
   q(X)=7200,
   Tr(X^3)=N^2=14400.
   ```

   Record whether this removes the negative vertices without introducing a
   fitted coefficient.

## Exact extremal criterion

No numerical optimization may prove a kernel locus.  On the zero-sum sphere

```text
sum_i r_i=0,             sum_i r_i^2=6/5,
```

use the exact equality case of Cauchy--Schwarz/Lagrange multipliers to classify
all points with `r_i=+1` or `r_i=-1`.

## Decision boundary

- **Operator bridge advance:** the overlap map is the unique equivariant map,
  the sign-neutral slack is canonical, and its kernel locus is exactly the
  twelve signed simplex vertices.
- **Conditional selection advance:** intersecting with the previously derived
  positive cubic condition leaves exactly the six actual `Box_i`.
- **Kill:** non-unique intertwiners, additional kernel points/continua, a
  normalization not fixed by the sharp norm, or failure of the old cubic to
  remove the negative set.

Even an advance is not yet a physical vacuum theorem.  The remaining
load-bearing hypothesis would be:

> the founding nontrivial-kernel bootstrap is also an axiom for the auxiliary
> six-fibration slack operator.

That reuse is **STRUCTURAL** unless independently justified.  Without it,
`K(X)` is a canonical recognizer whose determinant happens to vanish at its
own frame vertices; it does not dynamically force a field to saturate the
bound.
