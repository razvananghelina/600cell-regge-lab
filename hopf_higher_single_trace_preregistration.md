# Preregistration: higher single-trace Hopf selector audit

Date: 2026-08-10

## Provenance

The target cubic `C_box` is known, and the fourth-moment baseline test is
closed negative.  This is not a blind target discovery.  The protocol fixes
the entire finite single-trace search space and the privileged sixth-moment
test before computing either.

## Fixed operator family

Use exactly

```text
D(epsilon,X)=A+epsilon X,
X in W=span_R{Box_i},
q(X)=Tr(X^2)=7200.
```

No baseline besides `A`, no fibration weights, no inserted projector and no
chosen linear combination of trace words is allowed.

For each `p>=3`, define the action-selected cubic coefficient

```text
K_p(X) = [epsilon^3] Tr((A+epsilon X)^p).
```

The primary candidate is `K_6`, because the heat-trace Taylor expansion
already singles out the sixth moment with coefficient `-t^3/6` for `t>0`.
The formal Taylor coefficient does not license a heat-kernel interpretation
or a physical scale on the finite carrier.

## Complete finite word space

The registered adjacency spectrum has nine distinct eigenvalues.  Verify
coefficientwise that its minimal polynomial is

```text
m_A(z)=z(z-12)(z-3)(z+2)(z+3)(z^2-6z-36)(z^2-4z-16).
```

Hence every power of `A` reduces to degrees `0,...,8`.  Cyclicity of trace
and transpose reversal identify exponent triples up to all permutations.
The complete individual-word list is therefore

```text
T_abc(X)=Tr(A^a X A^b X A^c X),
0 <= a <= b <= c <= 8,
N_words = binomial(11,3) = 165.
```

This exhausts the cubic-in-`X` tensors obtainable from one trace with
arbitrary polynomial functions of `A` between the three `X` insertions.  It
does not license an arbitrary linear combination of the 165 words as an
action.

## Frozen exact tests

1. Reconstruct `A`, the six `Box_i`, `W`, `C_box` and the already proved
   two-dimensional `A5`-invariant cubic space.
2. Verify `m_A(A)=0` and that no proper square-free factor annihilates `A`.
3. Enumerate all 165 exponent triples before target comparison.  Use exact
   integer arithmetic or modular non-vanishing certificates followed by exact
   checks for every modular survivor.
4. Record:

   - the rank of the complete 165-word cubic span;
   - the number of individual words proportional to `C_box`;
   - their exponent triples, if any;
   - the number of distinct projective cubic lines among the 165 words.

5. Independently compute `K_p` for `p=3,4,5,6` from the full noncommutative
   word sum, not from a commuting binomial shortcut.  Record each coordinate
   in the exact basis `(Tr(X^3),C_box)`.
6. Compare `K_6` with `C_box`.  If it is proportional, use the fixed
   heat-trace sign `-K_6` and determine whether it favours `+Box_i`; then test
   the complete sixth moment `Tr((A+X)^6)` on the fixed sphere.  Desired
   vertices must first be exact stationary points.  An exact lower competitor
   is sufficient to kill global selection.
7. If `K_6` is not proportional, still record whether it opens the second
   invariant line.  A mixed cubic is not a selector without a separately
   derived cancellation of its other component.

## Evaluation certificate

Because the full invariant cubic space has dimension two and
`(Tr(X^3),C_box)` is a basis, two fixed evaluation points may be used as exact
coordinates only after verifying that their `2 x 2` evaluation determinant
is nonzero.  Modular arithmetic may rule out proportionality: a nonzero
determinant modulo a declared prime proves that the corresponding integer
determinant is nonzero.  Every modular zero must be checked exactly; it may
not be counted as a hit.

## Decision boundary

- **Sixth-moment advance:** `K_6` is a nonzero multiple of `C_box`, has the
  required fixed sign, and the complete sixth moment selects exactly the six
  desired vertices.
- **Sixth-moment kill:** `K_6` is on the old line, is a nontrivial mixture,
  has the wrong sign, makes the desired vertices nonstationary, or has an
  exact lower competitor.
- **Complete single-trace algebraic kill:** all 165 words remain on the old
  `Tr(X^3)` line.
- **Look-elsewhere result:** if isolated words hit `C_box` but `K_6` does not,
  report the hit fraction out of 165 as **PATTERN**.  Choosing such a word or
  a fitted linear combination after comparison is forbidden.

Even a positive polynomial result is **STRUCTURAL** until `D=A+X` is licensed
as a fluctuated Dirac operator satisfying the finite spectral-triple gates.
