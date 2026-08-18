# Prior-art gate: complete dynamics of the negative-stiffness shape sectors

Date: 2026-08-18

## Exact object, carrier and hypotheses

This is a **TARGET-DISCLOSED FOLLOW-UP**, not a second blind discovery test.
The preceding preregistered census selected minimal binary-tetrahedral sectors
`4` and `5` because each has the resolved Hermitian inertia

```text
15 negative + 10 positive
```

on its `25`-dimensional action-relative shape block.  Both irreps have
dimension one, so together they give the `30` certified negative directions
on the full `600`-position carrier.

In either selected sector let `W` be the already certified shape basis and

```text
Gamma_S = W* Gamma W,
Omega_S = W* Omega W.
```

The exact centered recurrence is

```text
q_2 - 2q_1 + q_0 + Gamma_S(q_2-q_0) + Omega_S q_1 = 0.
```

If `I+Gamma_S` is regular, it gives the unique doubled position-state map

```text
[q_0]       [ 0                         I                    ] [q_0]
[q_1]  ->   [-(I+Gamma_S)^-1(I-Gamma_S)  (I+Gamma_S)^-1(2I-Omega_S)] [q_1].
```

Call this finite `50 x 50` companion `C_S`.  It is one non-autonomous update
constructed from the first two accepted slabs; it is not assumed to be a
time-translation-invariant map that may be iterated forever.

Independently, let `E_-` be the `15`-dimensional negative eigenspace of the
Hermitian restricted stiffness form

```text
A=-W*[(V+V*)/2]W.
```

The first question is whether `E_-` is invariant under both `Gamma_S` and
`Omega_S`.  Only if it is may `E_- direct-sum E_-` be called an autonomous
negative-stiffness state carrier and receive its own `30 x 30` companion.

Complete hypotheses:

- the fixed regular 600-cell and `720` logarithmic signed-squared boundary
  edge coordinates;
- the first two accepted nonstationary fixed-total-mass dust-Regge slabs;
- the literal adjacent-slice identification of the centered Jacobi equation;
- both staircase schedules and all four derivative variants;
- the canonical conformal map and action-selected shape complement;
- the two sectors selected by the already committed negative-stiffness result;
- no exact constraint quotient, independent dust perturbation, proper-time
  normalization, refinement or continuum tensor target.

## Constraint-framing correction

The fixed carrier does not admit a canonical constraint quotient at present.
The already certified complete pre-Legendre Jacobian has

```text
rank 1560/1560,
error-consistent nullity 0.
```

Its `120` weak directions are resolved nonzero
`PSEUDOCONSTRAINT_CANDIDATE`s.  They are not a kernel.  Removing them by a
numerical threshold would be an arbitrary modification of the finite theory.

This is consistent with the primary formalism:

- Dittrich--Hoehn derive pre/post constraints and reduced phase spaces for
  variational discrete systems from the degenerate directions of the
  Lagrangian two-form: <https://arxiv.org/abs/1303.4294>.
- Hoehn's quadratic classification likewise uses null vectors of that form:
  <https://arxiv.org/abs/1407.6641>.
- Bahr--Dittrich show that curvature generically breaks exact Regge gauge
  symmetry and replaces constraints by background-dependent
  pseudo-constraints: <https://arxiv.org/abs/0905.1670>.

Therefore the fixed model must be tested as a regular dynamical system.  A
reduced quotient could become legitimate only after an exact symmetry,
extended matter carrier or controlled refinement limit is independently
derived.

## Primary prior art

- Marsden--West give the general variational origin of discrete
  Euler--Lagrange recurrences and their regular Legendre maps:
  <https://doi.org/10.1017/S096249290100006X>.
- Dittrich--Hoehn formulate the corresponding canonical evolution in
  simplicial gravity: <https://arxiv.org/abs/1108.1974>.
- Hoehn warns that lattice gravitons require gauge/curvature identification,
  not merely a raw tangent eigenvalue: <https://arxiv.org/abs/1411.5672>.
- Rostworowski's continuum FLRW control separates two gravitational master
  wave variables from a matter transport scalar after constraint reduction:
  <https://arxiv.org/abs/1902.05090>.

Constructing a companion matrix from a regular second-order recurrence and
testing invariant subspaces is **KNOWN LINEAR ALGEBRA**.  Negative stiffness
alone need not decide stability in the presence of first-difference coupling
and mode mixing.

No located primary source publishes these two dust-600-cell sector
companions, their negative-subspace leakage or their complete controlled
spectrum.  Search failure is not evidence of novelty; external novelty is
**OPEN**.

## KNOWN / CONTROL / OPEN

- **KNOWN:** the companion formula follows uniquely from the centered
  recurrence when `I+Gamma_S` is regular.
- **CONTROL:** the complete shape factor is reducing for both `Gamma` and
  `Omega` under the committed error model.
- **CONTROL:** sectors `4` and `5` and their `15+10` inertia are fixed by the
  preceding result, not selected from the new companion spectrum.
- **CONTROL:** the Hermitian kinetic form is definite and far from singular.
- **OPEN:** whether `E_-` is reducing for the actual recurrence.
- **OPEN:** regularity and non-normal conditioning of `I+Gamma_S` and `C_S`.
- **OPEN:** resolved contracting/expanding/unit/open counts of the complete
  companion under calibrated errors.
- **OPEN:** whether any expansion is a physical scalar/tensor instability,
  a frozen-dust artifact or a finite-discretization pseudo-constraint effect.

## Framing attack

1. The follow-up is conditional on a post-result target.  It can characterize
   the obstruction but cannot increase the evidential status of discovering
   the two sectors.
2. Eigenvalues of one non-autonomous companion are not Lyapunov exponents.
   Long-time growth requires further dynamically solved slabs and a cocycle
   norm.
3. Negative stiffness can mix with positive stiffness through `Gamma_S` and
   `Omega_S`.  If `E_-` is not invariant, projecting the recurrence onto it
   would choose forbidden Schur coefficients.
4. Even a resolved `|mu|>1` multiplier is local hyperbolicity in the declared
   finite coordinates, not automatically a ghost or a continuum graviton.
5. No comparison with desired `S^3` harmonics, two polarizations, `c` or
   Planck units is licensed by this mission.

The meaningful outcomes are therefore: isolated negative dynamics exists;
the negative stiffness is dynamically mixed and only the full sector has a
meaning; or the finite classification remains open.  All three constrain the
next step without fitting.
