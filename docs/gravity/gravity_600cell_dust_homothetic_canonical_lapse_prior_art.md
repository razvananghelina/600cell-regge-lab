# Prior-art gate: homothetic canonical lapse selection

Date: 2026-08-16

Status: **completed after the fixed-lapse negative and before evaluating a
two-variable Jacobian or solve**.

Upstream fixed-lapse result: `b788258`.

This is a targeted primary-source map, not proof of external novelty.

## 1. Exact object and complete hypotheses

Keep the published lower regular 600-cell boundary, fixed total dust mass,
both derived order-24 staircase schedules, complete Lorentzian Regge+dust
action, and the preceding slab's already committed post-momentum.

The next homothetic slab has exactly two unknown logarithmic coordinates:

```text
s = log(L_+/L0),
z = log(rho_next/rho0),
rho0=(0.0102)^2.
```

Its geometry is fixed, not fitted:

```text
q_old       = L0^2,
q_new       = exp(2s)*L0^2,
pole        = rho0*exp(z),
diagonal    = exp(s)*L0^2-rho0*exp(z).
```

The two invariant equations are

```text
F0(s,z) = mean of the five complete pole equations = 0,
F1(s,z) = mean[p_pre(s,z)-P p_post(static)]         = 0.
```

Every candidate must then be substituted into all 35 internal equations and
all 30 momentum equations.  A reduced two-scalar zero alone does not pass.

The mass and preceding canonical datum remain fixed.  Only the next slab's
positive lapse magnitude and upper spatial scale vary.

## 2. KNOWN

The conceptual mechanism is not new.  In consistent discretizations,
variables that are Lagrange multipliers in the continuum can be determined
by the discrete canonical equations, producing an unconstrained canonical
transformation:

- Gambini and Pullin, *Consistent discretization and canonical classical and
  quantum Regge calculus*, <https://arxiv.org/abs/gr-qc/0511096>;
- Gambini and Pullin, *Consistent discretizations for classical and quantum
  general relativity*, <https://arxiv.org/abs/gr-qc/0108062>.

In canonical simplicial gravity, pre/post momenta and evolution arise from
Hamilton's principal function, while a priori free data can become fixed a
posteriori:

- Dittrich and Hoehn, *Canonical simplicial gravity*,
  <https://arxiv.org/abs/1108.1974>.

Curved Regge backgrounds generally break exact vertex-displacement gauge and
replace constraints by background-dependent pseudo-constraints:

- Bahr and Dittrich, <https://arxiv.org/abs/0905.1670>;
- Dittrich and Hoehn, <https://arxiv.org/abs/0912.1817>.

De Felice and Fabri instead choose the 600-cell lapse explicitly and evolve
the remaining geometry:

- <https://arxiv.org/abs/gr-qc/0009093>.

Thus a lapse-like variable determined by finite discrete consistency is
**KNOWN STRUCTURE**, not by itself a discovery or proof of emergent time.

## 3. CONTROL

The new calculation must recover all committed fixed-lapse data before
moving in `z`:

- the nonzero root
  `s=-3.1160706675973036032799885512889e-6`;
- maximum internal residual `4.318e-32`;
- uniform per-component canonical mismatch
  `+1.61609399353141e-9`;
- exact even/odd agreement;
- the target orbit permutation and uncertainty from the two-slab artifact.

The solve begins only from `(s_root,z=0)`.  No alternate seed or root search
is allowed.

## 4. OPEN difference

The following are open:

1. whether the two-equation Jacobian has calibrated rank two;
2. whether a locally unique positive-lapse solution exists;
3. whether the reduced solution passes all 65 internal/momentum components;
4. whether the two schedule parities select the same pair;
5. whether the selected lapse differs significantly from `0.0102`;
6. whether any selection survives refinement or is a discretization
   pseudo-constraint;
7. whether the resulting slab can be iterated to another frame;
8. whether any absolute time scale is derived rather than inherited from
   the initial `tau0`.

No located primary source gives this exact two-variable solve for the
repository carrier.  External novelty remains **OPEN**.

## 5. Framing verdict before calculation

Even a successful solve cannot establish a fundamental clock.  The original
proper interval `tau0=0.0102` entered the preceding sandwich by hand and
sets the scale against which `rho_next` is measured.  A passing result would
show only that the discrete canonical map selects the *next* lapse relative
to already supplied canonical data.

The physical headline would still be meaningful but narrower: first locally
unique, canonically glued homogeneous next slab on this fixed finite carrier.
Calling it emergent time requires refinement stability and a rule that no
longer inherits an initial time unit.
