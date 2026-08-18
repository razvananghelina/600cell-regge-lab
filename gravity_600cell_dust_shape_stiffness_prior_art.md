# Prior-art gate: action-relative shape stiffness census

Date: 2026-08-18

## Exact object, carrier and hypotheses

This mission starts only from two already certified finite constructions:

1. the centered three-slice dust-Regge Jacobi equation

   ```text
   M (q_2 - 2 q_1 + q_0) + N (q_2 - q_0) + V q_1 = 0;
   ```

2. the action-selected reducing split, in every frozen minimal symmetry
   sector,

   ```text
   K = im C,
   S_H = ker(C* H_M),
   H_M = (M + M*)/2.
   ```

Here `C` is the literal unsigned vertex--edge incidence map

```text
(C sigma)_uv = sigma_u + sigma_v.
```

On the complete carrier, `dim K=120` and `dim S_H=600`.  The normalized
centered operators `Gamma=M^-1 N` and `Omega=M^-1 V` preserve both factors
within the committed numerical/ball envelope.

For an orthonormal basis `W` of `S_H`, define two distinct finite objects:

```text
M_S       = W* H_M W,
V_S       = W* H_V W,       H_V=(V+V*)/2,
Omega_S   = W* Omega W.
```

The Hermitian generalized pencil is

```text
V_S x = lambda M_S x.
```

The invariant normalized block `Omega_S` is the actual coefficient in the
shape recurrence.  They coincide only when the relevant adjointness relations
hold; this mission must measure their disagreement rather than assume it away.

Complete hypotheses:

- the fixed regular 600-cell and its literal `720` spatial edge coordinates;
- the first two accepted nonstationary fixed-mass dust-Regge slabs;
- the same adjacent-slice edge identification used by the certified Jacobi
  equation;
- all seven minimal binary-tetrahedral sectors, both staircase schedules and
  all four frozen derivative variants;
- the certified conformal/shape carrier and its stored uncertainty model;
- no independent dust perturbation, gauge quotient, continuum refinement or
  proper-time normalization.

## Primary literature

- Marsden and West's discrete variational mechanics gives the general setting
  in which a regular discrete Lagrangian produces a three-point linearized
  recurrence and a canonical update:
  <https://doi.org/10.1017/S096249290100006X>.
- Dittrich and Hoehn derive action-generated canonical evolution for
  simplicial gravity and explain the role of pre/post constraints:
  <https://arxiv.org/abs/1108.1974>.
- Bahr and Dittrich show that curvature generically breaks exact discrete
  diffeomorphism symmetry and replaces constraints by background-dependent
  pseudo-constraints: <https://arxiv.org/abs/0905.1670>.
- Hoehn identifies lattice gravitons only after gauge and curvature degrees of
  freedom have been separated in linearized Regge calculus:
  <https://arxiv.org/abs/1411.5672>.
- Hartle, Miller and Williams show that the Lund--Regge supermetric has an
  indefinite, triangulation- and point-dependent signature; extra negative
  directions need not be gauge: <https://arxiv.org/abs/gr-qc/9609028>.
- Christiansen relates the quadratic three-dimensional Euclidean Regge action
  to `curl^T curl` and proves an eigenpair convergence result in that stated
  setting: <https://arxiv.org/abs/1106.4266>.  It does not identify the present
  four-dimensional Lorentzian dust pencil.
- Rostworowski's continuum FLRW control obtains two gravitational master
  scalars satisfying wave equations and one matter transport scalar only
  after solving the continuum perturbation constraints:
  <https://arxiv.org/abs/1902.05090>.
- De Felice and Fabri evolve a generalized 600-cell cosmology, but do not give
  this complete action-relative shape pencil or its spectrum:
  <https://arxiv.org/abs/gr-qc/0106077>.

The use of a quadratic-action kinetic/stiffness pencil is therefore
**KNOWN**.  A spectrum of the unreduced metric Hessian is not by itself a
graviton spectrum.

No located primary source computes the exact `600`-position carrier, the
present dust background, or the comparison between its Hermitian pencil and
its invariant normalized recurrence block.  Search failure is not proof;
external novelty remains **OPEN**.

## KNOWN / CONTROL / OPEN

- **KNOWN ALGEBRA:** if `M_S` is definite, the Hermitian generalized
  eigenproblem has a real spectrum.
- **CONTROL:** the source Jacobi coefficients are regular, variational and
  schedule robust under the committed ball model.
- **CONTROL:** `S_H` is selected by the action kinetic form, is complementary
  to the canonical conformal image and is reducing for `Gamma` and `Omega`.
- **CONTROL:** same-dimensional Euclidean spectral and Fourier splits do not
  reduce the recurrence, so closure is not a dimension-count artifact.
- **OPEN:** the sign census of the generalized shape stiffness.
- **OPEN:** whether `Omega_S` is self-adjoint in the restricted kinetic form
  strongly enough to agree with the Hermitian pencil.
- **OPEN:** a constraint-derived physical quotient of `S_H`.
- **OPEN:** identification with transverse-traceless tensor modes, a spatial
  Lichnerowicz operator, dispersion, refinement or an effective speed.

## Framing attack

The words `shape`, `stiffness` and `mode` are only structural at this stage.
In particular:

1. `600` shape coordinates are not `600` physical gravitons.  Scalar,
   longitudinal, pseudo-constraint, matter-coupled and discretization modes
   may remain.
2. A positive generalized stiffness is only a necessary diagnostic for the
   simplest undamped oscillator reading.  The recurrence also contains the
   background-drift term `Gamma`, and the background is nonstationary.
3. A resolved negative value would refute a blanket positive-stiffness claim,
   but would not alone prove a ghost or kill Regge gravity.  Cosmological and
   constraint sectors can grow without being propagating negative-norm
   particles.
4. Symmetrizing `M` and `V` is an action diagnostic, not a license to replace
   the actual recurrence.  Agreement with `Omega_S` must be tested.
5. Comparing the resulting numbers with a desired `S^3` Laplacian spectrum,
   desired degeneracies or a desired value of `c` before deriving a physical
   carrier map would be fitting and is forbidden here.

The mission is therefore a blind finite-operator census.  It can close a
naive stable-wave shortcut or establish a clean necessary condition.  It
cannot by itself open the physical graviton gate.
