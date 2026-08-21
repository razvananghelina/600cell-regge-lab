# Prior-art gate: generic-velocity next-order lapse selection

Date: 2026-08-21

Status: written after the adversarially corroborated leading generic-velocity
theorem and its post-result search, before evaluating any next-order action,
constraint or momentum coefficient.

## 1. Exact object and complete hypotheses

Retain the certified homogeneous cellular 600-cell Regge-plus-dust action on
the positive Lorentzian branch, with zero cosmological constant and one
conserved physical dust mass `M`.  Set the incoming spatial scale to one by
the already-proved global scale covariance.  Exclude the singular
time-symmetric stratum `v=0`, which has its own completed theorem.

For a fixed real `v!=0`, freeze the same physical incoming state for every
candidate duration:

```text
M  =mu(v),
p0 =p(v),

mu(v)=180 epsilon_v/[pi sqrt(v^2+4)],
p(v)=180 v epsilon_v/sqrt(v^2+4)-600 sqrt(3) eta(v),

epsilon_v=2*pi-5*acos((v^2+2)/(2(v^2+3))),
eta(v)=asinh(v/sqrt(8(v^2+3))).
```

Let `h=sqrt(rho)>0` be the actual proper height of one slab, not an external
coordinate label.  The local endpoint ansatz is

```text
L_plus=exp(v h+a h^2+O(h^3)).
```

No coefficient of `M` or `p0` may depend on whether a coarse or half slab is
being tested.  Allowing `M(h)` or `p0(h)` would compare different states and
repeat the defect already exposed in the turning-point `lambda` hierarchy.

## 2. Two questions that must not be conflated

### Lapse selection

For fixed `(L_minus,p0,M)`, the exact one-slab equations are

```text
F=rho partial_rho S=0,
p_pre=-L_minus partial_(L_minus)S/2=p0.
```

They may admit a local one-parameter family in `h`, no nonzero local branch,
or isolated finite positive roots.  Only the last possibility is a candidate
state-dependent relational duration.  Even then it is not a fundamental
time quantum without a separate principle selecting the carrier and showing
stability under refinement.

### Composition/refinement

One coarse slab and two half slabs are compared by gluing at a stationary
intermediate boundary and matching canonical momenta.  Equality tests whether
the discrete principal function is exact/perfect to the registered order.
It does not by itself derive a tick.  Exact equality for every `h` instead
supports reparametrization and divisibility; failure can be ordinary
discretization error.

This separation corrects the earlier shorthand that treated step doubling as
if it were itself the tick-selection criterion.

## 3. Primary prior art

- Bahr and Dittrich show that curved finite Regge discretizations generically
  break continuum gauge symmetry and replace constraints by lapse-dependent
  pseudo-constraints: [(Broken) Gauge Symmetries and Constraints in Regge
  Calculus](https://arxiv.org/abs/0905.1670).
- Dittrich and Hoehn formulate simplicial evolution through Hamilton's
  principal function, pre/post momenta and stationary gluing:
  [Canonical simplicial gravity](https://arxiv.org/abs/1108.1974) and
  [From covariant to canonical formulations of discrete
  gravity](https://arxiv.org/abs/0912.1817).
- Marrero, Martin de Diego and Martinez construct the exact discrete
  Lagrangian associated with continuous Hamiltonian flow; a generic discrete
  Lagrangian is an approximation to this object:
  [On the exact discrete Lagrangian function for variational integrators](https://arxiv.org/abs/1608.01586).
- Bahr and Dittrich construct improved/perfect actions designed to reproduce
  continuum dynamics and restore the corresponding gauge symmetry:
  [Improved and Perfect Actions in Discrete Gravity](https://arxiv.org/abs/0907.4323).

The literature therefore predicts the logical alternatives but not the
coefficients of the present 600-cell action.  The search found no primary
source printing this generic-velocity next-order branch census.  That search
does not prove novelty; external novelty remains **OPEN**.

## 4. KNOWN / CONTROL / OPEN

### KNOWN

- The exact action is globally scale covariant, so no absolute classical
  duration can emerge from it alone.
- At generic velocity, the leading action, lapse constraint and pre-momentum
  are exactly independent of the interval factor.
- At the time-symmetric scaling, the same-state analytic half-step is absent.
- Stationary gluing, not endpoint addition alone, is the correct two-slab
  composition law.

### CONTROL

- Differentiate the complete exact action before or consistently with the
  series expansion; retain the boundary term and fixed-mass chain rule.
- Use the same exact `M=mu(v)` and `p0=p(v)` in the coarse and fine histories.
- Enumerate the complete common-root set of the first nonzero lapse and
  momentum coefficients before any sampled velocity.
- If a one-slab jet exists, impose both second-slab lapse stationarity and
  intermediate momentum matching before comparing final endpoints.
- Keep lapse selection and composition outcomes in separate fields.

### OPEN

- Whether one next-order endpoint coefficient satisfies both fixed-state
  one-slab equations for symbolic `v!=0`.
- Whether exceptional nonzero velocities exist.
- Whether the two-half-slab stationary history exists and matches the coarse
  history to the same order.
- Whether the exact finite equations have isolated positive lapse roots.
- Whether any such root is stable under carrier refinement or is only a
  pseudo-constraint artifact.

## 5. Framing attack and decision boundary

A common next-order endpoint coefficient for every generic `v` means that
the duration remains locally free to this order: **no tick selected**.

No common coefficient means only that the exactly leading continuum state
does not lie on an arbitrary-small-step discrete solution family.  It is a
finite-discretization obstruction, not yet a derived finite tick.  An exact
finite positive-root census would then be required.

An isolated positive root is only a **candidate relational lapse**.  It may be
promoted further only if it is unique on a stated physical branch, stable
under perturbations/refinement and selected without changing the incoming
state.  Seconds or Planck time still require one dimensionful physical scale
or a separately derived dimensional-transmutation mechanism.

