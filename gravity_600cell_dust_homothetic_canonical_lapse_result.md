# Homothetic canonical-lapse solve: precision-corrected result

Date: 2026-08-16

## Provenance

- prior-art gate: `c7f3e29`;
- original frozen protocol: `ded77c5`;
- unresolved coarse result: `d30854c`;
- precision-correction protocol: `3c34a59`;
- corrected implementation committed before evaluation: `e4784ad`;
- corrected artifact SHA-256:
  `4b1c59c0518eec11b88b140cdecdf558d762c0d70b4826a758f67544e14ac5b9`.

Only the registered targeted verifier was run.  It returned **7/7**.  The
full suite was not run, following the user's explicit instruction.

## Mechanical verdict

```text
HOMOTHETIC_CANONICAL_LAPSE_SELECTED
```

Both independently derived staircase parities converge in two undamped
Newton steps to

```text
s = -3.1160595766945016917347064441986394e-6,
z = -3.5592531351706334372503053396391740e-6,

L_+/L0          = 0.9999968839452782140983335686805886,
rho_next/rho0   = 0.9999964407531989632917428351804216,
tau_next/tau0   = 0.9999982203750159491039423415664135.
```

For the inherited `tau0=0.0102`, this is

```text
tau_next = 0.0101999818478251626808602118839774.
```

The scale contracts by `3.11605` parts per million and the next proper lapse
is smaller by `1.77962` parts per million.

## Numerical falsification gates

At the endpoint, for both parities:

```text
||F||_infinity                 = 1.0247254e-26,
||p_pre-p_target||_2           = 5.6126524e-26,
allowed junction bound         = 3.6513654e-21,
max |five pole equations|      = 4.2574746e-31,
max |thirty diagonal equations|< 1.7e-96.
```

The endpoint Jacobian has

```text
singular values = (582.8522169181..., 4.2445618171e-9),
epsilon         = 1.3301400039e-22.
```

Thus the weak singular value is about `3.19e13` times its calibrated error.
The rank-two conclusion is resolved, but the condition number remains about
`1.37e11`: the selected direction is physically/numerically soft even though
it is not zero on this carrier.

Even/odd differences are below `3e-83` in `(s,z)` and below `3e-86` in the
post-momenta.  The result is therefore not a staircase-parity artifact at the
tested precision.

## What selected the lapse

The strongest version of the motivating claim is refuted.  Conserved mass
alone does **not** select `tau`: the complete pole equation is one equation in
the two unknowns `(s,z)` and therefore generically defines a curve.

The local pair is selected only by imposing both

```text
mean(complete five pole equations) = 0,
mean(p_pre - committed p_post)     = 0.
```

The first is the conserved-mass Regge equation; the second is canonical seam
consistency with already supplied boundary data.  Their calibrated Jacobian
has rank two.  This distinction is load-bearing: calling the result "time
from mass conservation" would overstate it.

## Interpretation ledger

- **DERIVED COMPUTATIONAL LOCAL:** on the fixed regular 600-cell carrier, at
  fixed total dust mass and fixed incoming canonical datum, the homothetic
  two-variable problem has the recorded locally unique non-static root.
- **DERIVED COMPUTATIONAL:** the candidate satisfies all 35 internal
  equations and all 30 canonical matching components, not only the two
  reduced means.
- **STRUCTURAL:** the selected next slab is contracting for the chosen
  canonical orientation.  It is the first canonically glued non-static
  homothetic step in this calculation, not yet a cosmological trajectory.
- **STRUCTURAL / candidate pseudo-constraint:** the relative lapse is selected
  by a very soft finite-triangulation direction.  Curved Regge calculus is
  known to replace exact gauge constraints by background-dependent
  pseudo-constraints.
- **OPEN:** full 65-variable uniqueness, stability against anisotropic
  perturbations, iteration to a second independently solved step, time-reverse
  expanding partner, refinement stability, continuum convergence and whether
  the weak lapse selection survives rather than tending to gauge.
- **NOT DERIVED:** an absolute time unit, the speed of light, Planck time,
  emergent time, inflation or the history of the observed universe.  The
  initial `tau0` was inherited.

## Post-result primary-source gate

The refined search terms were `canonical Regge momentum matching`, `lapse
selected by discrete consistency`, `600-cell dust canonical evolution` and
`curved Regge pseudo-constraint`.

The mechanism remains known: consistent discretization can determine
continuum multiplier-like variables and implement dynamics as a canonical
transformation (Gambini--Pullin,
<https://arxiv.org/abs/gr-qc/0511096>); simplicial pre/post data can fix
a-priori free data in later moves (Dittrich--Hoehn,
<https://arxiv.org/abs/1108.1974>); and curved Regge backgrounds replace exact
gauge constraints by pseudo-constraints (Bahr--Dittrich,
<https://arxiv.org/abs/0905.1670>).  Published 600-cell work evolves enlarged
variable sets with Sorkin's scheme but does not, in the located material,
report this exact two-scalar fixed-mass/canonical-seam root
(De Felice--Fabri, <https://arxiv.org/abs/gr-qc/0106077>).

Therefore the general physical mechanism is **KNOWN**.  No located primary
source gives this exact carrier, action, target and root, but external novelty
remains **OPEN** until a dedicated review.

## Next falsification

Do not infer a full evolution law from one local root.  The next decisive and
minimal test is to use this accepted slab's own post-momentum as the target
for a second homothetic solve, with mass still fixed and no new seed search.
Before evaluation, preregister three outcomes: continued contraction,
time-reversed/branch failure, or collapse of the weak Jacobian toward gauge.
Only successful iteration begins to justify the word "dynamics".
