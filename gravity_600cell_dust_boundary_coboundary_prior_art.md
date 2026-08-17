# Prior-art gate: boundary-coboundary equivalence of the two dust schedules

Date: 2026-08-17

Status: completed before evaluating any new on-shell action difference or
mixed rectangle.

## 1. Exact object, carrier and hypotheses

For each already derived ordered schedule

```text
p in {even, odd},
```

let

```text
S_p(o,x,n)
```

be the same Lorentzian Regge--dust slab action used by the certified canonical
calculation.  Here `o`, `x`, and `n` are respectively the logarithms of the 30
old-boundary squared edge-orbit lengths, 35 internal squared-length variables,
and 30 new-boundary squared edge-orbit lengths.  The carrier is the fixed
order-24 quotient of the 600-cell slab: 30 old, 35 internal and 30 new
variables.  The dust mass, accepted dynamic background, Lorentzian logarithm
branch and unique already selected physical-edge permutation `P` are frozen.

For boundary data in the local branch, solve only

```text
partial_x S_p(o,x,n) = 0
```

and define Hamilton's principal function

```text
W_p(o,n) = S_p(o,x_p(o,n),n).
```

The comparison object is

```text
Delta(o,n) = W_odd(P o,P n) - W_even(o,n).
```

This is a Dirichlet/on-shell action comparison.  It is not the earlier
initial-value comparison at common `(o,p_pre)`, whose two solutions generally
have different final boundaries.

The narrow question is whether, on the tested local branch, the two principal
functions can differ only by endpoint terms,

```text
Delta(o,n) = B_old(o) + B_new(n) + constant.                 (1)
```

Such terms shift pre/post momenta but do not change the mixed Lagrangian
two-form.  The more restrictive time-homogeneous coboundary

```text
Delta(o,n) = F(n) - F(o) + constant                         (2)
```

is not assumed.  Equation (2) would require additional equality between the
two endpoint functions after identifying consecutive slices.

## 2. Exact mathematical criterion

For any four boundary points in a product neighbourhood define the rectangle

```text
R_Delta(o0,o1;n0,n1)
 = Delta(o1,n1) - Delta(o1,n0)
 - Delta(o0,n1) + Delta(o0,n0).                             (3)
```

- **DERIVED:** equation (1) implies `R_Delta = 0` identically.
- **DERIVED:** for a twice differentiable `Delta`, equation (1) on a connected
  product neighbourhood is equivalent to

  ```text
  partial_o partial_n Delta = 0
  ```

  there.  Integrating the vanishing mixed derivative yields the two endpoint
  functions.
- **DERIVED:** one resolved nonzero rectangle falsifies (1) for the present
  action, carrier, branch and boundary identification.
- **STRUCTURAL:** finitely many zero rectangles do not prove (1) on a
  neighbourhood.

The criterion is target-blind: it does not use a desired continuum equation,
speed, chirality sign or experimental number.

## 3. Primary literature

- Dittrich and Hoehn, [*Canonical simplicial
  gravity*](https://arxiv.org/abs/1108.1974), use Hamilton's principal
  function as the generator of discrete pre/post canonical evolution.  This
  is the direct framework for `W_p`.
- Dittrich and Hoehn, [*From covariant to canonical formulations of discrete
  gravity*](https://arxiv.org/abs/0912.1817), derive the canonical dynamics
  from the covariant action and show that nonlinear Regge terms can break
  linearized symmetries and produce pseudo-constraints.  Hence equality of
  tangent maps does not imply equality of nonlinear principal functions.
- Bahr and Dittrich, [*Improved and Perfect Actions in Discrete
  Gravity*](https://arxiv.org/abs/0907.4323), define improved actions by
  solving refined bulk equations at fixed boundary data and evaluating the
  action on that solution.  This motivates comparing the on-shell boundary
  functions rather than raw off-shell slab formulas.
- Marsden, Patrick and Shkoller, [*Multisymplectic Geometry, Variational
  Integrators, and Nonlinear
  PDEs*](https://authors.library.caltech.edu/records/74nvs-vb440), derive the
  discrete symplectic structure from derivatives of an action function.  An
  endpoint addition changes boundary one-forms but has no mixed derivative.

No located primary source compares these two ordered dust-filled 600-cell
schedules.  External novelty is **OPEN**.

## 4. Repository controls

- **DERIVED:** both schedules have the same accepted homothetic dynamic
  boundary values and the same 35 internal values after the physical
  identification.
- **DERIVED:** the complete 65 by 65 pre-Legendre Jacobian is locally regular
  for each schedule.
- **DERIVED:** their complete canonical tangent maps agree under the unique
  physical-edge permutation within the registered calibration.
- **DERIVED NEGATIVE:** their nonlinear canonical maps disagree on all 32
  frozen anisotropic initial-value rays, with a quadratic-compatible defect.
- **CONTROL:** the four zero-sum unit boundary directions and `ETA=1e-4` were
  fixed in an older target-blind nonlinear protocol.
- **OPEN:** whether the nonlinear disagreement is entirely an endpoint
  momentum convention of form (1).

## 5. Framing attack and scope

A nonzero rectangle would remove the cheapest explanation of the observed
nonlinear discrepancy, but it would not by itself prove new continuum
physics.  A nonlinear field redefinition, a more general canonical
transformation mixing endpoints, refinement, or the full 720-edge carrier
could still alter the comparison.  Conversely, zero on a finite direction
set would be only a **PATTERN**, not proof of a perfect action or schedule
equivalence.

The test must therefore use common Dirichlet boundary data, solve all 35 bulk
equations separately, retain the certified Lorentzian branch, compare
operational and validation solvers, and report every tested ordered direction
pair.  No result may be promoted beyond this fixed quotient and local branch.

## 6. Status before calculation

- **KNOWN:** action-generated discrete evolution and endpoint-term ambiguity.
- **CONTROL:** action, branch, carrier, dynamic base, `P`, four directions and
  scale `ETA`.
- **OPEN:** whether any target-blind mixed rectangle of `Delta` is nonzero.
- **NOT TESTED:** all 29 shape directions, a neighbourhood theorem, nonlinear
  boundary coordinate changes, refinement and the full carrier.
