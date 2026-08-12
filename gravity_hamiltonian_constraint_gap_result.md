# The current theory has no Hamiltonian metric constraint

Date: 2026-08-12

Preregistered protocol commit: `6c151e1`

Registered verifier:
`reproducible/verify_gravity_hamiltonian_constraint_gap.py`

Machine-readable result:
`reproducible/gravity_hamiltonian_constraint_gap.json`

## Headline

> **DERIVED CURRENT HAMILTONIAN-CONSTRAINT GAP.** None of the currently
> certified canonical dynamical constructions supplies a phase space for the
> 600-cell metric together with lapse/shift, a physical first-class
> Hamiltonian/diffeomorphism constraint and a Hamiltonian that evolves that
> metric. Therefore the theory cannot currently dismiss the negative
> conformal `A2` directions as constraint or gauge modes.

This is an inventory theorem about the present repository, not a no-go
against a future discrete gravity construction. It also does **not** prove
that the conformal modes are experimentally propagating.

The targeted verifier passes `15/15`. No full-suite run was performed.

## 1. Why the framing matters

The Euclidean conformal-factor sign is not, by itself, a proof of Lorentzian
instability. Canonical general relativity has spatial metric variables and
conjugate momenta, while lapse and shift impose scalar and vector constraints.
The benchmark structure is the original
[ADM Hamiltonian formulation](https://arxiv.org/abs/gr-qc/0405109).

Accordingly there are two different claims:

1. `A2` is an indefinite Euclidean functional — now derived both smoothly and
   at the finite Regge point;
2. the negative scalar directions survive as physical Lorentzian degrees of
   freedom — **not derived**.

Claim 2 could fail if the theory produced the right first-class constraint.
The present audit asks whether it already has one.

## 2. Required coverage

A current construction was counted as a metric Hamiltonian route only if it
contained all six items:

| field | required content |
|---|---|
| `Q_metric` | dynamical spatial metric, here naturally the 720 admissible edge lengths |
| `P_metric` | conjugate momenta/extrinsic-curvature data |
| `N,N^a` | lapse/shift or a derived discrete replacement |
| `H,H_a` | physical first-class scalar/vector constraints |
| closure | preservation and constraint algebra |
| evolution | the same Hamiltonian advances `Q_metric,P_metric` |

No currently certified candidate covers all six.

## 3. What the local tick actually evolves

The exact local unitary construction acts on

```text
2640 fixed cochains
14880 directed incidence arcs.
```

It has a strict cone of one Hasse incidence per micro-tick and is a genuine
mathematical evolution. But its discriminant is the normalized signed
incidence Kähler--Dirac operator of one fixed complex. Its own certificate
leaves Lorentzian physical selection open.

The refinement audit is more decisive: the unweighted tick is not a lift of
the accepted metric Whitney codifferential on an anisotropic barycentric
tetrahedron. Thus it neither evolves the 720 edge lengths nor supplies their
conjugate momenta. Calling its flip permutation an ADM “shift” would be a word
collision, not a physical derivation.

## 4. What the Whitney constraints actually constrain

The canonical Whitney constraints identify duplicated local cochain values.
For positive local mass matrix `M`, their bracket Gram matrix is

```text
G=C M^-1 C*.
```

The exact classification gives:

| level | real second-class | physical first-class |
|---|---:|---:|
| base 600-cell | 12,720 | 0 |
| first barycentric | 306,240 | 0 |

Their row-cycle redundancy acts only on multiplier coordinates. Dirac
reduction reproduces the already assembled cochain evolution; it does not
create metric evolution or a Gauss/ADM gauge sector.

The minimal auxiliary conversion does make

```text
Phi=C u+eta
```

first class, but the unique gauge-invariant coordinate and Hamiltonian contain
`G^-1`. The global solve is relocated rather than removed. The later weak
first-class audit is stronger: even under the correct weak condition, the
load-bearing quadratic Hamiltonian block has exact support at maximum
600-cell distance and fails endpoint locality.

These are useful constraint theorems. They simply concern cochain assembly,
not the metric Hamiltonian constraint required here.

## 5. Relation to the saddle results

The smooth `l=2` conformal direction is already proven not to be a spatial
diffeomorphism or scale variation: its scalar-curvature variation is
`delta R=40f`. The finite edge Hessian has no quotient null directions beyond
global scale and has inertia `(569,0,150)`.

Those facts prevent two shortcuts:

- the negative modes cannot be relabelled as already established gauge zeros;
- the absence of Hessian nullity cannot substitute for a Hamiltonian
  constraint, which is a phase-space/dynamical statement.

Smooth `A4` supplies a positive higher-derivative conformal term, but no
selected finite heat time. It is not a constraint and is UV-suppressed.

## 6. Status ledger

| Claim | Status |
|---|---|
| Fixed incidence tick is a local unitary evolution | **DERIVED KINEMATIC** |
| That tick evolves the 600-cell metric | **REFUTED AS CURRENTLY DEFINED** |
| Whitney copy constraints are physically first class | **REFUTED** |
| Their physical first-class count is zero at two levels | **DERIVED** |
| Minimal conversion creates a local metric gauge theory | **REFUTED IN THE STATED CLASS** |
| Weak first-classness restores endpoint locality | **REFUTED IN THE STATED CLASS** |
| Current theory contains a Hamiltonian metric constraint | **DERIVED ABSENT** |
| Current work removes the conformal modes as gauge/constraint | **REFUTED** |
| Conformal modes are experimentally propagating instabilities | **OPEN / NOT CLAIMED** |
| A future selected first-class metric extension is impossible | **NOT CLAIMED** |
| Lorentzian gravity, `G` or Planck units are derived | **OPEN** |

## 7. The precise missing construction

The next gravity step is no longer vague. It must supply, without choosing a
carrier after seeing which option removes the 150 negative modes:

1. a canonical metric phase space across slices;
2. a selected local symplectic structure;
3. lapse/shift multipliers or a derived alternative;
4. physical first-class constraints with a checked algebra;
5. a Hamiltonian evolving the metric with a causal/refinement limit.

There is an immediate canonicity problem: lapse could be placed on vertices,
tetrahedra, dual cells or chambers, and a 4D slab/product structure is not yet
selected. Those alternatives must be enumerated before testing whether one
repairs the conformal sector.

The static theory therefore has a promising spatial operator and a kinematic
clock, but still no dynamical spacetime metric. This is presently the central
gravity gap.

No PDF was built.
