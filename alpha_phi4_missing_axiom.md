# Alpha tree-level coefficient: minimal missing axiom

Date: 2026-07-22

## Decision

The identification

`1/alpha_0 = 4 a_1 phi^4 = 20 phi^4`

remains **STRUCTURAL**. The separate algebraic/spectral identities entering
the right-hand side are **DERIVED**, but no unique variational principle in
the current discrete theory forces their product to be the physical inverse
electromagnetic coupling.

## Derived content

- **DERIVED:** `4 a_1 = N/b_1 = 20`.
- **DERIVED:** for the decagonal Hopf fiber,
  `lambda_f = 2-phi = phi^{-2}`, hence `lambda_f^{-2} = phi^4`.
- **DERIVED:** the gauge-projected second moment satisfies
  `Tr(Box_gauge^2)/n_base = N` in the normalization currently used.
- **DERIVED:** therefore
  `Tr(Box_gauge^2)/(n_base b_1 lambda_f^2) = 4 a_1 phi^4` as a finite spectral
  identity.
- **DERIVED:** the Hopf-fiber holonomy supplies the geometric angle `2 pi`.
- **NORMALIZATION:** the constant term `1` in the alpha quadratic.

## Why this is not yet a coupling derivation

Let a discrete electromagnetic field be represented by `A`, with curvature
linearized as `F(A)`. Replacing the physical identification by
`A_phys = s A`, for any `s > 0`, rescales the quadratic action coefficient by
`s^{-2}` while leaving `Box`, its spectrum, the Hopf fibers, and all identities
listed above unchanged. The existing discrete variational problem therefore
does not fix the absolute field normalization and cannot uniquely select
`20 phi^4` as `1/alpha_0`.

This rescaling freedom is the precise obstruction. Searching additional heat
moments of the same operator cannot remove it: spectral numbers do not define
which normalized discrete cochain is the physical electromagnetic potential.

## Minimal missing axiom

The minimal extra input is a **STRUCTURAL field-normalization axiom**:

> The physical `U(1)` potential is the image of a specified normalized
> discrete 1-cochain under a discrete-to-continuum map, and that map is an
> isometry between a stated discrete inner product and the continuum
> `L^2` inner product on one full Hopf fiber/base cell.

Together with a fixed `U(1)` charge unit, this axiom would make the coefficient
of `||F||^2` meaningful and test whether the spectral expression above is
forced. Merely declaring the spectral expression to be `1/alpha_0` is an
equivalent normalization ansatz, not a proof.

## Status of the quadratic

- `2 pi`: **DERIVED** as Hopf holonomy; its use in the physical coupling
  equation remains part of the **STRUCTURAL** matching.
- `4 a_1 phi^4`: **STRUCTURAL** as a physical coupling coefficient; its
  internal spectral equality is **DERIVED**.
- `1`: **NORMALIZATION**.
- The physical alpha equation as a whole: **STRUCTURAL**.
- A variational principle that fixes the normalization without the added
  axiom: **OPEN**.
