# Local conserved-mass homothetic slab: result

Date: 2026-08-16

## 1. Provenance

- prior-art gate: `8865346`;
- frozen protocol: `428330e`;
- pre-evaluation outcome clarification: `ff8d352`;
- registered implementation before evaluation: `7b28d2f`;
- targeted verifier:
  `reproducible/verify_gravity_600cell_dust_homothetic_mass_conservation.py`;
- artifact:
  `reproducible/gravity_600cell_dust_homothetic_mass_conservation.json`;
- artifact SHA-256:
  `72225b1ca17de18f6d77aac43972f4fdca18e24575c8640e8be5e5636316fad0`.

Only the targeted verifier was run.  It returns **10/10**.  The full suite
was not run, following the active instruction.

## 2. Exact geometric reduction

**DERIVED.**  For unit 600-cell vertices, adjacent dot product `phi/2`
gives edge chord `1/phi`; hence a regular 600-cell of edge `L` has embedding
circumradius `phi*L`.  Two homothetic slices with lower and upper edges
`L_-`, `L_+` and same-vertex proper separation `tau` have staircase
cross-edge square

```text
d^2 = L_-*L_+ - tau^2.
```

This is an exact consequence of the five-dimensional Lorentzian embedding,
not a fitted internal length.  The verifier rederives it symbolically.

## 3. Frozen derivative result

Set

```text
L_- = L0,
L_+ = exp(s)*L0,
tau = 0.0102,
M   = (90/pi)*(2*pi-5*acos(1/3))*L0,
d^2 = exp(s)*L0^2-tau^2.
```

At `s=0`, fresh 100-decimal evaluations reproduce the exact static-family
theorem.  The largest local residual is `2.09e-99` for the even schedule and
`1.89e-97` for the odd schedule; the exact pre/post momenta are reproduced
below `2e-97` relative error.

At the six preregistered nonzero scale displacements, every action and
gradient evaluation remains Lorentzian and away from the angle branch cut.
The independent complete-action derivative in `log(tau^2)` reproduces the
chain-rule lapse equation at all twelve parity/state points.

The fourth-order Richardson derivatives are identical between schedules
within the frozen error bands:

```text
||d R_local/ds|| = 1.218350019332e-2
proxy norm        = 7.839e-9

d E_lapse/ds      = 2.724313463615e-2
proxy             = 1.753e-8.
```

The mechanical outcome is

```text
LOCAL_HOMOTHETIC_STATIC_ONLY_GLOBAL_AND_LOCAL.
```

**DERIVED COMPUTATIONAL LOCAL.**  Since the exact static solution line is
`s=0` for every admitted positive lapse, the resolved transverse derivative
isolates that line in some neighborhood of the published point by the
implicit-function theorem.  This statement is local and does not exclude a
second nearby root.

## 4. Stronger structure learned after the frozen verdict

The 35-vector derivative has a special support:

```text
30 diagonal components: zero within about 1e-92 or better,
5 pole components:       0.005448626927230728... each.
```

More strongly, at all twelve non-static parity/state points the complete
diagonal residual norm is between about `6e-99` and `2e-96`, while the pole
residual alone accounts for the displayed `R_local` norm.  Thus the tested
homothetic geometry satisfies every one of the 30 independent staircase
diagonal equations to 96 decimal orders, not merely their average.

**PATTERN, not yet an exact theorem.**  The evidence says that on the
homothetic embedding all 30 diagonal equations vanish identically and the
complete 35-equation system reduces to five equal pole equations.  An exact
angle/Schlaefli derivation is still required before calling this an identity.

## 5. Preregistered data predict a second, contractant root

The sampled lapse equation is positive at both frozen negative scale points

```text
E_lapse(-2.5e-5) = 4.782867351584e-6,
E_lapse(-5.0e-5) = 2.049140388683e-5,
```

while its resolved derivative at zero is positive.  Smoothness on the
certified branch therefore gives **DERIVED COMPUTATIONAL EVIDENCE** for a
nonzero negative root between zero and `-2.5e-5`.  Quadratic interpolation
of only the preregistered values predicts

```text
s_root approximately -3.116121868e-6,
L_+/L_- approximately exp(s_root),
```

a contraction of about 3.1 parts per million.  The numerical value is
**PATTERN** until an independently frozen root solve is performed.

The boundary-momentum derivatives are uniform across all 30 components:

```text
mean d p_pre/ds  = -582.8542908882002...,
mean d p_post/ds = -582.8533827838162....
```

Combining the first number with the predicted root changes
`p_pre=-0.0009081044...` toward `+0.0009081044...`, which is precisely the
preceding slab's post-momentum.  This is a striking **PATTERN / preregistered
prediction**, not yet a canonical-junction result.  No root or momentum at
the predicted point was evaluated in this mission.

## 6. Physical interpretation and framing correction

If the predicted root passes all 35 equations and the independent junction,
it will be the first nondegenerate homogeneous forward slab found by the
repository.  It would be a discrete Friedmann tick: a chosen positive
proper interval `tau=0.0102` produces a derived scale contraction.

It would **not** mean that mass conservation selected the duration of the
tick.  This mission held `tau` fixed.  The stronger Claude proposal

```text
mass conservation selects tau and supplies the missing clock
```

is therefore not established.  At most, fixed mass plus the constraint may
select the scale change for a chosen lapse.  A physical clock, absolute time
unit, `c` and Planck scales remain **OPEN**.

The root would also be disconnected from the already proved positive-lapse
static branch of the artificial momentum homotopy.  Its selection must come
from the exact homothetic geometry and forward canonical junction, not from
choosing an unrelated numerical root after inspection.

## 7. Post-result primary-source audit

The learned distinction is known in broad form.  Liu and Williams explicitly
compare global and local variation in closed Regge--FLRW models and find that
local variation does not generally give the same viable model:
<https://arxiv.org/abs/1501.07614>.  De Felice and Fabri evolve a fixed-mass
dust 600-cell with a chosen lapse and solve for the next geometry:
<https://arxiv.org/abs/gr-qc/0009093>.  Broken discrete gauge symmetry and
pseudo-constraints are standard Regge phenomena:
<https://arxiv.org/abs/0905.1670>.

No located source prints the present exact homothetic diagonal relation,
complete 35-orbit cancellation pattern and canonical-root prediction for
this order-24 staircase pair.  External novelty remains **OPEN**; a search
does not prove novelty.

## 8. Next falsification test

Before another evaluation, freeze a bracket entirely inside
`(-2.5e-5,0)`, a deterministic root method, full Lorentzian gates and these
two independent acceptance tests:

1. all 35 internal orbit residuals vanish at the root, not only
   `E_lapse`;
2. `p_pre(root)` equals the previously committed `P p_post(static)` in all
   30 components.

If the pole equation has a root but either test fails, the candidate is a
minisuperspace pattern and not a physical tick.  If both pass, expansion or
contraction is derived kinematically for one slab, while clock selection and
multi-tick evolution remain open.
