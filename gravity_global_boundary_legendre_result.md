# Global Regge boundary-Legendre rank

Date: 2026-08-13

Frozen protocol: `gravity_global_boundary_legendre_protocol.md`

Protocol commit: `8c2482b`

Registered verifier:
`reproducible/verify_gravity_global_boundary_legendre.py`

Machine-readable result:
`reproducible/gravity_global_boundary_legendre.json`

Targeted run: **33/33 passed**.  The full repository suite was deliberately
not run at the user's request.

## Headline

> **DERIVED COMPUTATIONAL:** for either ordered-schedule parity, the complete
> invariant Regge action has 65 variables: 35 internal and 30 final-boundary.
> The internal-equation Jacobian has stable ranks
> `rank(J)=rank(J_internal)=35` and `rank(J_final)=30` at all three frozen
> thresholds.  The old and final boundary derivatives are real, independently
> reproduced by complete-action differences, and obey regular time reversal.

This establishes a valid off-shell boundary-Legendre formulation.  It does
not establish an on-shell tick: no stationary point was searched in this
audit.

## 1. Fixed construction and independent controls

For each of the two phase-order parity representatives the variables are

```text
30 internal staircase diagonals
 5 internal pole magnitudes
30 final-boundary edge lengths
--------------------------------
65 variables.
```

Old-boundary edges stay fixed at the regular value.  Their 30 action
derivatives are nevertheless reconstructed as pre-momentum orbits.  Only the
first 35 derivatives are equations of motion; final-boundary derivatives are
post-momenta and are not set to zero.

The verifier reconstructs the orbit structure:

```text
old boundary:    720 edges = 30 orbits x 24
final boundary:  720 edges = 30 orbits x 24
boundary hinges: 2400 triangles = 100 orbits x 24.
```

At all three frozen controls `B0`, `B1`, `B2`, both parity representatives
remain on the real Lorentzian branch.  The 100-simplex-orbit evaluator agrees
with the direct sum over all 2400 four-simplices in the action, every one of
the 65 variable derivatives, all 30 old-boundary derivatives and every hinge
curvature.  Relative discrepancies are at most order `1e-14`.

Centered differences of the complete action independently reproduce all 65
derivatives at `B0` and `B1`; the largest reported relative error is
`1.20e-8`.  This rules out an orbit multiplicity or boundary-sign convention
as the source of the result.

## 2. Frozen rank result

At the regular, nonstationary control `B0`, the Jacobian is

```text
J = [J_internal | J_final],
shape(J)          = 35 x 65,
shape(J_internal) = 35 x 35,
shape(J_final)    = 35 x 30.
```

The ranks are identical at relative thresholds `1e-7`, `1e-9`, `1e-11`:

| parity | `rank(J)` | `rank(J_internal)` | `rank(J_final)` |
|---|---:|---:|---:|
| even | 35 | 35 | 30 |
| odd  | 35 | 35 | 30 |

The smallest singular values are well separated from the frozen thresholds:

| parity | combined | internal | final |
|---|---:|---:|---:|
| even | 0.176396 | 0.0375126 | 0.103418 |
| odd  | 0.180616 | 0.0330180 | 0.103608 |

The internal-block condition numbers are `71.89` and `81.78`.  Thus the rank
statement is not being created by a barely resolved numerical mode.

The final-boundary block has its maximum possible rank 30, but it maps into
only a 30-dimensional subspace of the 35 internal equations.  Its
five-dimensional left-null complement is therefore unavoidable.  At `B0`,
varying only the final boundary leaves residual norms

```text
even: 2.3393416551
odd:  2.3333840992
```

from a common starting norm `3.1380483927`.  Hence a final-boundary change
alone cannot remove the regular internal residual even at linear order.

Allowing internal and final variables together solves the linearized residual
to `8.21e-15` (even) and `7.51e-15` (odd).  The corresponding minimum-norm
steps have norms about six, so this is a tangent-space controllability result,
not evidence for a nearby nonlinear root.

## 3. Boundary momenta

With the standard generating-function convention,

```text
p_pre  = -dS/dq_old,
p_post =  dS/dq_final.
```

Both 30-component momentum vectors are real.  Their multisets transform into
one another under complete reversal of the regular staircase slab with
relative error `2.04e-15`.  Direct complete-action differences reproduce the
old and final derivatives with errors between `4.85e-9` and `5.79e-9`.

These are **DERIVED COMPUTATIONAL OFF SHELL** momenta.  Since the internal
gradient at `B0` is nonzero, they do not yet define physical on-shell initial
and final data.

## 4. Attack on the framing and prior art

The broad construction is established Regge calculus, not a new invention.
The parallel tent-move evolution, nonadjacent vertex batches and staircase
triangulation appear in Barrett et al.:
[A Parallelizable Implicit Evolution Scheme for Regge Calculus](https://arxiv.org/abs/gr-qc/9411008).
That paper's claimed four classes of 30 were corrected to five classes of 24
by De Felice and Fabri:
[The Friedmann universe of dust by Regge Calculus](https://arxiv.org/abs/gr-qc/0009093)
and
[Singularities of the closed RW metric in Regge Calculus](https://arxiv.org/abs/gr-qc/0106077).
Canonical boundary momenta and evolution generated by a discrete action are
also standard:
[Canonical simplicial gravity](https://arxiv.org/abs/1108.1974).

The 2000 calculation already recognizes three inequivalent old-edge types
between each pair of color classes.  For the newly created slant edges it then
imposes additional equalities, fixes each pole as lapse, and solves five
successive systems of four equations in four unknowns.  It includes dust and
checks residual equations that were omitted under its lapse/shift treatment.
Later closed-universe work also uses a positive cosmological constant.

The current calculation does **not** impose those within-class-pair edge
equalities: the three stabilizer orbits remain independent.  It checks all 35
internal orbit derivatives of the complete 2400-simplex action and all 30
independent final-boundary responses.  But it uses a pure curvature action:
no dust and no cosmological-constant volume term.  Consequently it is neither
a bibliographically new tent construction nor a reproduction of the known
Friedmann tick.  Whether this exact 35/65-orbit audit has appeared elsewhere
remains **OPEN**.

This distinction also weakens the original motivation.  A regular closed
Friedmann universe is not a nontrivial vacuum solution with zero cosmological
constant.  An asymmetric compact vacuum solution is not ruled out merely by
that observation, but searching for one is a much less canonical physical
target than testing the established dust or positive-`Lambda` controls.

## 5. Status ledger

| Claim | Status |
|---|---|
| The 65-variable evaluator is the exact invariant restriction of the complete slab action | **DERIVED COMPUTATIONAL** |
| Both parity branches are real and Lorentzian at all controls | **DERIVED COMPUTATIONAL** |
| The combined/internal/final ranks are 35/35/30 | **DERIVED COMPUTATIONAL** |
| Final-boundary variations alone cancel the regular residual | **DERIVED NEGATIVE at linear order** |
| Combined variations cancel every linearized residual direction | **DERIVED COMPUTATIONAL at `B0`** |
| Real pre/post momentum orbits obey regular time reversal | **DERIVED COMPUTATIONAL OFF SHELL** |
| A nonlinear stationary tick exists | **OPEN; not searched here** |
| A stationary tick is near `B0` | **NOT SUPPORTED; the linear step is large** |
| Tent moves or five-phase 600-cell evolution are new | **REFUTED by prior art** |
| The present 35/65-variable full-action audit is externally novel | **OPEN pending dedicated review** |
| The pure zero-`Lambda` vacuum is the physically preferred next target | **OPEN / weakly motivated** |

## 6. Correct next test

Do not launch another broad blind vacuum search merely because the rectangular
Jacobian has full row rank.  The rank makes such a search numerically possible;
it does not make its physical target canonical.

The next protocol should first add two preregistered external controls on the
same exact slab and orbit reduction:

1. a positive cosmological-constant volume term, for comparison with closed
   vacuum `Lambda`-FLRW Regge models;
2. the homogeneous dust term used in the published 600-cell evolution.

Only after reproducing at least one known physical evolution should the pure
zero-`Lambda` asymmetric vacuum continuation be interpreted.  Failure of the
controls would expose an action/branch/boundary mismatch; success would make
the remaining vacuum negative or positive scientifically legible.
