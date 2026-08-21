# Prior-art gate: cubic generic-velocity formal integrability

Date: 2026-08-21

Status: written after the adversarially corroborated first-correction theorem
and its post-result search, before evaluating any cubic endpoint coefficient.

## 1. Exact object and hypotheses

Retain the certified homogeneous cellular 600-cell Regge-plus-conserved-dust
action, positive Lorentzian branch and zero cosmological constant.  Set the
incoming scale to one.  Use exactly the accepted incoming state

```text
M=mu(v),
p0=p(v),
```

with real

```text
v!=0,
K(v^2)!=0.
```

The two already-classified points `v=+-v_star` are excluded because no
quadratic endpoint jet exists there.  Freeze the accepted coefficient

```text
a(v)=-B(v^2)/K(v^2)
```

and extend only the endpoint:

```text
L_plus=exp(v h+a(v)h^2+c h^3),
rho=h^2,
h>0.
```

Do not change `M`, `p0`, `v`, `a(v)` or the carrier while solving for `c`.

## 2. Exact question

The accepted branch cancels the leading and first-correction lapse and
momentum residuals.  Extract the next coefficients

```text
C2(v,c)=lim_(h->0+) [2F/h]/h^2,
P2(v,c)=lim_(h->0+) [p_pre-p0]/h^2.
```

Classify every real common root in `c` on the complete domain
`v!=0, K(v^2)!=0`, including all degree-drop or denominator loci.

This is a formal-integrability test.  It is not an exact finite-height root
census and cannot by itself derive a tick.

## 3. Prior art

- [Bahr--Dittrich](https://arxiv.org/abs/0905.1670) show that nonlinear
  corrections in curved Regge calculus replace exact gauge constraints by
  lapse-dependent pseudo-constraints.
- [Dittrich--Hoehn](https://arxiv.org/abs/0912.1817) show how higher-order
  broken symmetries appear as consistency conditions on lower-order gauge
  parameters in canonical discrete gravity.
- [Marrero--Martin de Diego--Martinez](https://arxiv.org/abs/1608.01586)
  relate a discrete Lagrangian's local error to the exact discrete
  Lagrangian.
- [Schmitt--Leok](https://arxiv.org/abs/1609.02309) connect resonances of
  variational integrators with ill-posed boundary-value problems for their
  generating functions.
- [Bahr--Dittrich perfect actions](https://arxiv.org/abs/0907.4323) provide
  the refinement/improvement route if a finite action fails formal
  composition or gauge consistency.

No located primary source prints the present `K`, `B` or cubic 600-cell
coefficient.  The search does not prove novelty; external novelty remains
**OPEN**.

## 4. KNOWN / CONTROL / OPEN

### KNOWN

- Leading generic duration is a reparametrization.
- The quadratic endpoint coefficient is uniquely `a=-B/K` away from exactly
  two velocities.
- At those two velocities no quadratic endpoint jet exists.
- A global scale theorem forbids an absolute classical time unit.

### CONTROL

- Reconstruct the previous `C1`, `P1`, `K`, `B` identities before extracting
  `C2`, `P2`.
- Use exact second total derivatives of the scaled expressions and direct
  unexpanded arbitrary-precision controls.
- Record both coefficients of `c`, the cross-resultant and every exceptional
  real locus before sampling.
- Treat the intervals `0<|v|<v_star` and `|v|>v_star` separately when proving
  signs or denominators.

### OPEN

- Whether `C2` and `P2` have a common root generically.
- Whether new isolated degree-drop velocities appear.
- Whether formal duration freedom continues to all orders.
- Whether any finite positive duration is selected by the exact equations.
- Stability under spatial/carrier refinement and physical interpretation of
  any exceptional value.

## 5. Framing boundary

A common `c(v)` is **DERIVED formal integrability through cubic endpoint
order**, not an exact flow and not a tick.  An incompatible pair is a
**DERIVED NEGATIVE local obstruction**, not a finite tick.  A degree-drop
value is a discretization singularity until refinement and inhomogeneous
propagation say otherwise.  No result in this gate yields seconds, `c`, `G`
or Planck units.

