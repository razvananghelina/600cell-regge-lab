# Prior-art gate: conserved inhomogeneous dust-mass response of one 600-cell tick

Date: 2026-08-18

Status: completed before preregistration and before evaluating a mass-response
matrix.  This is not a novelty proof; external novelty remains **OPEN**.

## 1. Exact proposed object and complete hypotheses

Use the already certified, non-static, homogeneous dust solution on each of the
two five-stage 600-cell schedules.  Keep unchanged:

- the complete 2,400-simplex Lorentzian Regge slab;
- all 720 old, 840 internal and 720 new logarithmic squared-edge variables;
- the accepted angle branch and total background mass `M`;
- the 120 comoving pole world-lines, one per 600-cell vertex;
- the fixed incoming canonical data `(q_old,p_old)`.

Replace only the homogeneous world-line weights `m_v=M/120` by conserved
parameters

```text
m_v = M/120 + delta m_v,
sum_v delta m_v = 0
```

at linear order.  The declared matter action is therefore still

```text
S_dust = -8 pi sum_v m_v sqrt(rho_v).
```

The zero-sum condition holds total mass fixed.  The uniform mass variation is
retained only as a control.  No particle position, velocity, dust proper-time
field or conjugate dust momentum is added.

Let `J` be the already certified complete pre-Legendre Jacobian from the 840
internal plus 720 new-boundary variations to the 840 internal equations plus
720 old momenta.  At the homogeneous background the mixed source derivative is
fixed analytically:

```text
B_m = partial(internal equations, old momenta) / partial m_v,
```

with one nonzero entry per pole equation.  Since `J` is regular, the forced
slab response at fixed incoming canonical data is

```text
Y_m = -J^-1 B_m.
```

The proposed new object is the induced outgoing boundary-phase response

```text
R_m : delta m in R^120 -> (delta q_new, delta p_new) in R^1440,
```

with the new boundary returned to the frozen canonical vertex/edge ordering.
The physically contentful questions are the rank and orientation of the
zero-total-mass image, its internal-curvature response, and its relation to the
already committed expanding/contracting tangent subspaces.  The rank of the
unprojected solved-variable map `Y_m` is an algebraic control, not evidence: it
follows from injectivity of `B_m` and invertibility of `J`.

## 2. KNOWN

De Felice and Fabri use point-dust world-line terms in the Sorkin evolution of
the 600-cell and keep the total dust mass fixed:

- A. De Felice and E. Fabri, *The Friedmann universe of dust by Regge
  Calculus: study of its ending point*, arXiv:`gr-qc/0009093`;
- A. De Felice and E. Fabri, *Singularities of the closed RW metric in Regge
  Calculus: a generalized evolution of the 600-cell*,
  arXiv:`gr-qc/0106077`.

Inhomogeneous mass perturbations in closed lattice cosmology are direct prior
art.  Liu and Williams perturb one mass in a different Collins--Williams
reduction and find a well-behaved response:

- R. G. Liu and R. M. Williams, *Regge calculus models of closed lattice
  universes*, arXiv:`1502.03000`, DOI
  `10.1103/PhysRevD.93.023502`.

Thus no broad claim that inhomogeneous Regge dust cosmology is new is allowed.
The located paper does not compute the complete `120 -> 1440` forced canonical
response of the present Sorkin slab or compare it with its full anisotropic
tangent invariant spaces.

The distinction between a conserved particle weight and dynamical dust clock
fields is load-bearing.  Brown and Kuchar introduce dust proper time and
comoving labels as canonical fields, with their own conjugate momenta:

- J. D. Brown and K. V. Kuchar, *Dust as a Standard of Space and Time in
  Canonical Quantum Gravity*, arXiv:`gr-qc/9409001`, DOI
  `10.1103/PhysRevD.51.5600`.

Gauge-invariant perturbations of such reference-field models contain extra
matter degrees of freedom and coupled scalar equations; they are not recovered
merely by differentiating fixed world-line masses:

- K. Giesel, B.-F. Li and P. Singh, *Relating dust reference models to
  conventional systems in manifestly gauge invariant perturbation theory*,
  arXiv:`2012.14443`.

## 3. CONTROL

- The background total mass and all 120 equal weights reproduce the accepted
  homogeneous dust action.
- `B_m` is derived analytically from the same action, not by finite
  differencing or fitting.
- Its nonzero rows are exactly the 120 pole equations selected independently
  by edge incidence.
- The full `J` blocks, tangent blocks and binary-tetrahedral bases must
  reproduce their frozen artifacts before a new comparison is read.
- A uniform `delta m` must agree with an independent finite difference of the
  homogeneous source term.
- A zero source must give exactly zero forced response.

## 4. OPEN

- whether the zero-sum mass response has rank 119 after projection to outgoing
  boundary phase;
- whether it is identified with, separated from, or numerically unresolved
  against either 119-dimensional strong tangent branch;
- whether its four-dimensional Regge-curvature response is injective;
- whether the result is schedule robust;
- a dynamical dust clock, dust velocity/momentum perturbations, a
  gauge-invariant scalar/vector/tensor split, refinement and a continuum limit;
- external novelty of this exact finite response map.

## 5. Attack on the framing

The proposed calculation does **not** add matter degrees of freedom to phase
space.  The `m_v` are conserved source parameters.  Consequently a successful
alignment with a strong tangent branch would identify that branch as the
geometric sensitivity to comoving density inhomogeneity; it would not by
itself prove a propagating dust mode, a clock, or a scalar cosmological
perturbation.

Conversely, separation would refute only this simplest density-source
explanation.  It would not turn the strong modes into gravitons: constraint
violation, longitudinal geometry and discretization artifacts would remain.

The zero-sum restriction is fixed before calculation by conservation of total
mass, not selected to obtain dimension `119`.  The numerical equality with the
already known `119` strong pairs is disclosed in advance, so any comparison is
confirmatory and must report all schedule/sector attempts.  No continuum
spectrum, desired speed, polarization count or particle target may be loaded.

## 6. Decision

Proceed to a target-disclosed, fully enumerated protocol.  This route is
preferable to imposing another spectral quotient: its source and response are
selected directly by the already declared action.  It remains a fixed-carrier
source-sensitivity experiment, not the refinement test ultimately required to
recover continuum gauge symmetry.
