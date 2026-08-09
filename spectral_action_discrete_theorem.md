# Discrete Spectral-Action Theorem

Working note for the spectral-action repair, before folding the result back
into `one_integer_paper.tex`.

## Exact discrete statement

Let `X` be the simplicial 600-cell and let

- `D = d + d*`
- `H = R^V \oplus R^E \oplus R^F \oplus R^C`
- `D^2 = Delta_0 \oplus Delta_1 \oplus Delta_2 \oplus Delta_3`

with `f`-vector

- `V = 120`
- `E = 720`
- `F = 1200`
- `C = 600`.

Then the finite spectral moments are exactly

- `c_0 = Tr(I) = 2640`
- `c_1 = Tr(D^2) = 14880`
- `c_2 = (1/2) Tr(D^4) = 55920`.

After dividing by `2N = 240`, one gets the reduced triple

- `A_0 = 11`
- `A_1 = 62`
- `A_2 = 233`

which satisfies the exact identity

- `2 A_1^2 + 1 = 3 A_0 A_2`.

In the paper's polynomial parametrization,

- `A_0 = 2 a_1 + 1`
- `A_1 = 2 a_1^2 + 2 a_1 + 2`
- `A_2 = 6 a_1^2 + 15 a_1 + 8`

and the deficit factors as

- `2 A_1^2 + 1 - 3 A_0 A_2 = (a_1 - 5)(8 a_1^3 + 20 a_1^2 + 16 a_1 + 3)`.

So the finite-moment Diophantine identity singles out `a_1 = 5` among
positive integers.

These numbers are not Seeley--DeWitt coefficients.  On the finite Hilbert
space the exact heat trace is analytic at the origin,

`Tr exp(-t D^2) = c_0 - c_1 t + c_2 t^2 + O(t^3)`,

and tends to `2640`, not to a negative power of `t`.  Therefore the triple
contains no heat-kernel dimension.  It remains a genuine arithmetic identity
among finite matrix moments.

## Exact gauge-side statement

From the `A_5` action on the 12-vertex icosahedral figure:

- `12 = 1 + 3 + 3' + 5`
- `ad(SU(3))|_{A_5} = 3 + 5`

hence the exact gauge-dimension decomposition is

- `12 = 1 + 3 + 8`.

This is the discrete bosonic channel split that is rigorously available from
the current framework.

## What is not fixed by this theorem

The following stronger claims are still not established by the discrete proof
above:

- a canonical 4D almost-commutative spectral triple of Connes type
- a controlled continuum limit from the 3D simplicial model to the full 4D
  Chamseddine-Connes action
- any inference of a geometric or spacetime dimension from the finite moment
  triple
- a rigorous derivation of the continuum gauge prefactors
  `(8/15, 1/3, 2/15)` as opposed to the exact discrete split `1:3:8`

So the safe statement is:

- the finite 600-cell operator fixes the exact moment triple
  `(2640, 14880, 55920)` and the exact discrete gauge-channel skeleton
  `1 + 3 + 8`;
- matching this to the full continuum Standard Model bosonic Lagrangian remains
  an interpretation step, not a theorem.
