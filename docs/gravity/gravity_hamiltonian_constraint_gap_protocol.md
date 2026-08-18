# Preregistered protocol: Hamiltonian-constraint coverage of the current theory

Date: 2026-08-12

Status at registration: **PROTOCOL ONLY -- NO NEW DYNAMICAL CONSTRUCTION**

## 1. Question and framing

The negative conformal Hessian of an Euclidean curvature functional is not by
itself a theorem that Lorentzian gravity is unstable. In canonical general
relativity scalar metric data are constrained, and lapse/shift act as
multipliers rather than propagating graviton components. The relevant primary
benchmark is the ADM formulation:
<https://arxiv.org/abs/gr-qc/0405109>.

Therefore the current question is deliberately narrower:

> Does any already certified dynamical construction in this repository
> provide a metric phase space and a physical first-class Hamiltonian/
> diffeomorphism constraint capable of removing or controlling the conformal
> metric directions found in the `A2` Hessians?

This is an audit of the **current repository**, not a no-go theorem against a
future discrete gravity theory.

## 2. Frozen inventory

The authoritative current candidates are frozen before the audit:

1. the signed Kähler--Dirac local unitary tick and its refinement audit;
2. the canonical constrained Whitney action and its Dirac--Bergmann
   classification;
3. the minimal linear first-class conversion;
4. the weak first-class local-Hamiltonian audit;
5. the newly certified smooth and finite `A2` saddles.

No partition, matter, flavor, Hopf-label or phenomenological script can count
as a gravity constraint merely because it uses the word “constraint”.

## 3. Necessary coverage fields

For this audit, a construction counts as a present Hamiltonian-gravity route
only if its certificate specifies all of:

```text
Q_metric   : dynamical spatial metric variables
P_metric   : conjugate metric momenta / extrinsic-curvature data
N, N^a     : lapse/shift or a derived discrete replacement
H, H_a     : physical first-class scalar/vector constraints
closure    : preservation/constraint algebra
evolution  : the same Hamiltonian advances Q_metric and P_metric
```

On the finite 600-cell the natural configuration **arena**, not yet a theory,
is the open set of 720 admissible edge lengths. Merely possessing a unitary
operator on 2640 cochains or 14880 incidence arcs does not make those edge
lengths dynamical.

A constraint that only enforces equality of duplicated cochain copies also
does not qualify. It must act on metric phase-space data and be physically
first class, rather than having redundancy only in multiplier coordinates.

## 4. Frozen mechanical checks

The verifier must read the committed machine certificates rather than rerun
their expensive calculations and must check:

1. the local tick acts on the fixed cochain/incidence carrier and explicitly
   leaves Lorentzian/physical selection open;
2. the refinement audit says this tick is not the accepted metric Whitney
   dynamics;
3. both certified Whitney levels have physical first-class count zero;
4. their constraints are second class and reduction reproduces the assembled
   cochain evolution;
5. minimal conversion introduces first-class auxiliaries but the unique
   physical dressing and Hamiltonian contain `G^-1`, so it supplies no local
   tick;
6. even weak first-classness fails the frozen endpoint-local quadratic
   Hamiltonian gate;
7. the smooth `l=2` conformal direction is non-gauge under spatial
   diffeomorphisms, and the finite scale quotient has no extra Hessian nulls;
8. none of the authoritative dynamics sources declares lapse, shift,
   extrinsic curvature or edge-length canonical momenta.

The last check is a scoped source-coverage guard, not a mathematical proof
that no differently named variable exists anywhere in history. The verdict
must rely primarily on the positive contents and scopes of the structured
certificates above.

## 5. Preregistered decision boundary

- **DERIVED CURRENT HAMILTONIAN-CONSTRAINT GAP:** every frozen certificate
  passes the inventory, but none covers the six necessary fields. The present
  theory therefore cannot dismiss its conformal negative directions as a
  derived Hamiltonian constraint or gauge mode.
- **REFUTED GAP:** an already committed candidate covers all six fields with a
  physical first-class constraint acting on metric data.
- **OPEN/INCOMPLETE AUDIT:** the frozen machine certificates do not support a
  complete inventory decision.

A gap result does **not** say the conformal modes are experimentally physical.
It says their removal has not been constructed. It also does not invalidate
the mathematical cochain tick, Whitney reduction, `A2` Hessians or smooth
`A4` coefficient.

## 6. Next-construction boundary

If the gap is confirmed, the next valid project must add, without fitting:

1. a canonical metric phase-space carrier across spatial slices;
2. a selected local symplectic/Poisson structure;
3. lapse/shift multipliers or a derived discrete alternative;
4. first-class constraints with a checked algebra;
5. a Hamiltonian that evolves the metric and has a refinement/causal limit.

Choosing whether lapse lives on vertices, tetrahedra, dual cells or chambers
after inspecting which choice removes the 150 negative modes is forbidden.
Those carrier choices must first be enumerated and preregistered.

Only the targeted audit verifier and a static registry check will run. No full
suite and no PDF build.
