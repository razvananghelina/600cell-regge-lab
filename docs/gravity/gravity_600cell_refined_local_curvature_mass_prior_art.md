# Prior-art and disclosure gate: local curvature mass on the refined static slab

Date: 2026-08-21

Status: written before a registered verifier, but after the exploratory
four-orbit numerical comparison disclosed below.

## Exact object and complete hypotheses

Use the fixed projected barycentric carrier

```text
K0=P(sd K_600),  f=(2640,17040,28800,14400),
```

its exact rank-derived chordal metric, the induced static Minkowski-product
fill with supplied positive proper time `tau0=0.0102`, all 24 colour-ordered
staircase triangulations, and the corrected Lorentzian Regge action with
boundary terms.  Restrict only by the common spatial `H4` action, which gives
four vertical-edge rank orbits.

For every spatial edge `e`, let

```text
C_e = l_e epsilon_e
```

be its three-dimensional Regge curvature contribution.  Split that hinge
contribution equally between its two endpoints and define

```text
K_v = (1/2) sum_(e incident on v) C_e,
K_r = sum_(rank(v)=r) K_v,
K   = sum_e C_e = sum_r K_r.
```

The already selected total mass is `M=K/(8*pi)`.  Replace the conditional
`P1` rank masses temporarily by four unspecified conserved rank totals
`mu_r`, with point-particle action

```text
S_dust = -8*pi sum_r mu_r sqrt(rho_r).
```

The narrow question is whether, at the induced product fill, the gravitational
log-lapse derivatives satisfy

```text
dS_grav/dlog(rho_r) = tau0*K_r/2.                 (1)
```

If (1) holds, the four internal equations are diagonal in the four masses and
select uniquely

```text
mu_r = K_r/(8*pi),                                (2)
```

whose sum is automatically the older global mass `M`.  This mission does not
solve a new internal geometry, vary the boundary, compute a Hessian or claim
that the selected rank density is homogeneous.

## Exploratory-result disclosure

This was not a blind discovery.  Before this gate was written, the four
committed `P1` lapse residuals were stripped of their known dust derivative
and compared with the endpoint-half spatial curvature totals.  The fractions
agreed to the precision available in the manually copied residuals:

```text
K_r/K =
(0.1287831657723389984...,
 0.3657000761313201399...,
 0.3759856918014127686...,
 0.1295310662949280930...).
```

The largest displayed discrepancy was approximately `2.2e-31`.  Therefore a
registered positive result is a formal validation and structural
reconciliation, not preregistered evidence discovered target-blind.  The
verifier must retain a corruption control and an independent action-gradient
route before the equality is accepted.

## KNOWN from primary literature

- Regge's action localizes the integrated scalar curvature on codimension-two
  hinges as hinge measure times deficit angle: T. Regge, *General Relativity
  Without Coordinates*, DOI `10.1007/BF02733251`.
- Vertex scalar curvature in Regge calculus requires a rule for distributing
  hinge curvature and a vertex volume; circumcentric-dual constructions are
  developed by McDonald and Miller, *A geometric construction of the Riemann
  scalar curvature in Regge calculus*, arXiv:`0805.2411`.  This warns that an
  endpoint-half split is a declared localization convention, not the unique
  scalar-curvature density on every irregular mesh.
- Local versus global variation and the initial-value/Hamiltonian constraint
  at a moment of time symmetry are standard issues in closed Regge
  cosmology: Liu and Williams, *Regge calculus models of the closed vacuum
  Lambda-FLRW universe*, arXiv:`1501.07614`.
- Conserved simplicial dust particles contribute minus mass times proper
  worldline length: Dittrich, Gielen and Schander, *Lorentzian quantum
  cosmology goes simplicial*, arXiv:`2109.00875`, DOI
  `10.1088/1361-6382/ac42ad`.

These sources motivate the objects in (1)--(2).  None of them proves the
four-rank identity for this projected 600-cell carrier or selects the present
endpoint-half localization as fundamental.  Search absence is not novelty
evidence; external novelty is **OPEN**.

## Repository controls

- **DERIVED:** the global static product identity fixes `M=K/(8*pi)` and
  leaves the common lapse arbitrary.
- **DERIVED COMPUTATIONAL:** with conditional `P1` masses `M/4`, every one of
  the 24 schedules has zero-compatible cross equations, four nonzero
  rank-lapse equations and a vanishing common-lapse sum.
- **DERIVED CONDITIONAL:** vertex-only affine-exact `P1` quadrature puts
  exactly one quarter of total mass in each rank; choosing the `P1` matter
  ansatz is **STRUCTURAL**.
- **CONTROL:** `sum_r K_r=K`, positivity of every `K_r`, rank four of the mass
  response, agreement of all 24 schedules and failure after a deliberate
  residual corruption.

## Framing attack

A positive result would not show that barycentric refinement preserves a
homogeneous dust density.  The `P1` dual volume is exactly `1/4` of the total
in each rank, whereas (2) is visibly nonuniform relative to that volume.
Equation (2) would instead be the unique rank-mass distribution that makes
this particular time-symmetric discrete geometry satisfy its four local
lapse equations.

That is legitimate discrete initial-data solving, analogous to imposing a
Hamiltonian constraint, but it changes the matter data.  It must not be
described as deriving ordinary homogeneous dust, a physical clock, a tick,
`c`, `G` or a Planck scale.  The strongest defensible outcome is:

- **DERIVED / STRUCTURAL:** local stationarity selects curvature-matched
  conserved rank masses on the fixed refined static geometry;
- **DERIVED NEGATIVE:** the earlier `P1` homogeneous-density ansatz is not
  locally on shell there;
- **OPEN:** whether the curvature-matched density converges toward a uniform
  continuum density under further refinement.

## Next admissible calculation

Freeze the exact rank geometry and the committed high-precision action
artifact.  Reconstruct `K_r` without using its lapse residuals, reconstruct
the gravitational lapse derivative by removing the analytically known `P1`
dust term, and test (1) with an error envelope inherited from the committed
100/140-digit calculation.  Print the four fractions before interpreting
them.  A separate adversarial verifier must reconstruct spatial incidences
and a dust-free action derivative by a mechanically different path.
