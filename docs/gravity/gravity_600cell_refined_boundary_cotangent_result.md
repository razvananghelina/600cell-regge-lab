# Consolidated result: the refined action selects a renormalized boundary cotangent

Date: 2026-08-21

Status: accepted within the frozen static-product scope after a mechanically
independent adversarial reconstruction.

## Complete hypotheses and scope

The carrier is the projected barycentric subdivision

```text
K0=P(sd K_600),  f=(2640,17040,28800,14400),
```

at the supplied product lapse `tau0=0.0102`.  The matter branch is the already
accepted curvature-selected one,

```text
m_v=K_v/(8*pi),
```

which makes the internal static equations stationary for all 24 staircase
schedules.  Boundary coordinates are the six orbit totals of logarithmic
squared-edge variables, ordered `(01,02,03,12,13,23)`, with canonical signs
`P_pre=-dS/dq_pre` and `P_post=+dS/dq_post`.

No unequal-boundary solution, Hessian, spectrum, nested root census, physical
clock, `c`, `G`, Planck unit or particle observable is computed here.

## Provenance ledger

| stage | commit |
|---|---|
| prior-art gate | `e7a1545` |
| primary protocol frozen | `5ee4f3e` |
| primary verifier registered before execution | `56b4ea6` |
| primary first `15/16` control failure preserved | `5ed87f6` |
| primary narrow correction preregistered | `8991447` |
| primary implementation correction | `a0edc31` |
| primary result | `f805557` |
| adversarial protocol frozen | `fbe6613` |
| adversarial verifier registered before execution | `2199d5e` |
| adversarial first `11/12` control failure preserved | `99940ad` |
| adversarial narrow correction preregistered | `3c35385` |
| adversarial implementation correction | `516d7d0` |

Both disclosed first failures were mixed-precision failures of a synthetic
kernel control.  Their scientific comparisons had passed, but no result was
accepted from either failed run.  The only corrections moved each boolean
comparison into its declared 100-decimal context; no geometric input,
equation, threshold or outcome hierarchy changed.

## Primary computation

The registered primary verifier

```text
reproducible/verify_gravity_600cell_refined_boundary_cotangent.py
```

passed `16/16` twice and produced byte-identical artifact

```text
reproducible/gravity_600cell_refined_boundary_cotangent.json
SHA-256 4e7bf0beb0327a3ee1bddbec13126fbef99380970e62cecf74eb24ce8d6dafaa.
```

Every staircase schedule selects

```text
P_pre =
(-0.0365507943663027407026225087...,
 -0.00243609423423285173839481667...,
 -0.0000655036546510496875793500932...,
 -0.0733546432343549192839929500...,
 -0.000989979010535868633659292999...,
 -0.0382237040192315742412603313...),

P_post=-P_pre.
```

The maximum pre/post schedule spreads are respectively `1.44e-96` and
`1.91e-96`.  The maximum internal cross residual is `1.87e-96`, the maximum
curvature-mass vertical residual is `5.91e-77`, and a Richardson derivative
of the complete action agrees to `1.73e-42` relative.

## Mechanically independent reconstruction

The adversarial verifier

```text
reproducible/verify_gravity_600cell_refined_boundary_cotangent_adversarial.py
```

does not import the primary functions and does not evaluate the Lorentzian
action.  It starts from the product-hinge Heron polynomial.  For one lower
boundary spatial edge,

```text
A^2=-l^2*tau^2/4,
dA/dlog(l^2)=i*l*tau/4,
P_pre,e=-tau*l_e*epsilon_e/4.                    (1)
```

It independently rebuilds all 17,040 actual refined edges and groups their
hinge curvatures `C_rs=sum l_e epsilon_e`.  Equation (1) then gives

```text
P_pre,rs=-tau0*C_rs/4,
P_post,rs=+tau0*C_rs/4.
```

All twelve components reproduce the primary result with maximum absolute
error `4.05e-77`; the unrefined regular 600-cell control is exact.  Dropping
one edge, changing the factor by two, reversing the sign and moving within
the five-dimensional pullback kernel all behave as preregistered.

The corrected verifier passed `12/12` twice and wrote a byte-identical
artifact:

```text
reproducible/gravity_600cell_refined_boundary_cotangent_adversarial.json
SHA-256 19c888a43bdba9d57166d6e3595c6d5b51dd019ebf616efdbf1189e25078f808.
```

Its frozen outcome is

```text
ADVERSARIAL_REFINED_BOUNDARY_COTANGENT_CORROBORATED.
```

## What is selected and what is not

The homothetic pullback row is `(2,2,2,2,2,2)`.  The action therefore gives

```text
p_pre,fine=-tau0*K_fine/2=-4*pi*tau0*M_fine,
p_post,fine=+tau0*K_fine/2=+4*pi*tau0*M_fine.
```

This selects one definite six-component covector inside the five-dimensional
inverse fiber left open by geometry and symplecticity alone.  The selection
is schedule-independent at first derivative order.

At fixed unit-volume radius, however,

```text
p_fine/p_coarse
 = K_fine/K_coarse
 = M_fine/M_coarse
 = 0.984190377388521915998... .                  (2)
```

Thus the bare refined and coarse momenta differ by about `1.58096%`.  Only the
mass-normalized quantity is invariant:

```text
(p_fine/M_fine)/(p_coarse/M_coarse)=1.
```

## Status ledger

- **DERIVED COMPUTATIONAL / STRUCTURAL, adversarially corroborated:** the
  curvature-matched on-shell action selects the complete refined boundary
  cotangent; it is not an arbitrary lift of the scalar momentum.
- **DERIVED:** its six components and homothetic pullback are independent of
  all 24 staircase schedules at first derivative order.
- **DERIVED NEGATIVE:** exact bare fixed-radius coarse/fine momentum transport
  fails by `1.58096%`; this is not yet a perfect action.
- **DERIVED identity within this product family:** momentum divided by the
  independently selected curvature mass is unchanged.  This follows from the
  shared curvature factor and is not by itself evidence for continuum
  universality.
- **OPEN:** whether the raw ratio tends to one under further controlled
  refinements, whether a perfect/improved action removes the mismatch, and
  whether the quadratic dynamics is schedule-independent.
- **OPEN / not computed:** a physical tick, causal speed, Newton constant,
  Planck units, particle masses or continuum gravity.

Primary literature establishes boundary momenta from a discrete action and
the need for improved/perfect actions under coarse graining; it does not, in
the sources located so far, contain this six-orbit projected-600-cell result.
Search absence is not a novelty proof, so external novelty remains **OPEN**.

## Registry and execution scope

The static coverage audit after registration reports `374` entries, `374`
distinct names, zero duplicates, zero unregistered verifiers and zero missing
files, with two explicitly reasoned exclusions.  Per the user's instruction,
no full suite and no deferred 12-case nested census was run.

## Next falsifiable gate

The first derivative is now selected, but propagation lives at second order.
Before extracting any wave speed, compute the reduced boundary Hessian at
this on-shell curvature-matched seed for all 24 schedules and compare the
physical quadratic forms, not their raw coordinate matrices.  If legitimate
staircase schedules give inequivalent reduced forms, this carrier has no
canonical quadratic evolution and cannot derive a unique `c`.  If they agree,
the common generalized eigenproblem against the spatial Laplacian is the
first admissible route to a dispersion relation.

