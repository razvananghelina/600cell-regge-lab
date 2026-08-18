# Anchored endpoint counterterms are refuted on all pure-momentum rays

Date: 2026-08-17

Only the lightweight targeted verifier was run.  No action was re-evaluated
and the full suite was not run.

## 1. Provenance

- prior-art gate: `20e9a26`
- preregistered corollary protocol: `b007bd3`
- verifier registered before component parsing: `65a07e6`
- frozen nonlinear input SHA-256:
  `a1e00071fa41f986dfaee84ea6e7689a14c50823f6c87d76889e6cb9346a7e3f`
- result artifact:
  `reproducible/gravity_600cell_dust_anchored_counterterm.json`
- result SHA-256:
  `1ee391c736a6acd119e5639ad56e26612d674edebc523914de531730f3d573b1`

This is explicitly a post-result corollary of nonlinear data committed in an
older blind protocol.  The underlying action outputs were not newly
preregistered; the component criterion and complete hypotheses were committed
before their configuration-only defects were parsed.

## 2. Complete theorem hypothesis

After the unique physical boundary identification, let

```text
Phi_p:(q,p_pre)->(Q,p_post),  p=even,odd,
```

be the two regular local canonical maps.  Suppose they differ only by
endpoint-generated vertical momentum translations and suppose both endpoint
functions are anchored so that their gradients vanish at the common base
configuration.  Anchoring preserves the already equal canonical momenta of
the base solution.

At fixed `q=q0`, the old-end translation is then the identity for every nearby
momentum, and the new-end translation changes only `p_post`.  Therefore the
final configurations must obey

```text
Q_odd(q0,p)=Q_even(q0,p)                                (1)
```

throughout the common branch.  This statement uses the regular initial-value
map and does not assume the invalid unconstrained Dirichlet function from the
preceding OPEN test.

## 3. Mechanical result

The verifier passed `8/8` controls and found

```text
ANCHORED_REFUTED: 16 / 16
OUTCOME: ANCHORED_ENDPOINT_COUNTERTERM_REFUTED_ON_FROZEN_RAYS.
```

Every selected ray keeps the old configuration exactly at `q0` and perturbs
only `p_pre`.  Nevertheless, the even and odd schedules give resolved
different final configurations.

| direction | configuration-defect range | defect / uncertainty range | order range |
|---:|---:|---:|---:|
| 1 | `1.6364e-19 .. 6.5457e-19` | `2.9049e30 .. 1.8900e33` | `2.000000062 .. 2.000000065` |
| 2 | `2.9795e-30 .. 1.1923e-29` | `1.1848e22 .. 3.9683e22` | `1.999814314 .. 2.000185645` |
| 3 | `3.0761e-30 .. 1.2308e-29` | `2.3092e21 .. 6.3205e23` | `1.999868252 .. 2.000131744` |
| 4 | `1.7020e-31 .. 6.8094e-31` | `2.5310e19 .. 5.7953e21` | `1.999897822 .. 2.000102189` |

All eight half/full pairs are `QUADRATIC_COMPATIBLE`.  Even the weakest case
is more than nineteen decimal orders above the conservative empirical
uncertainty.  The conclusion is not set by a marginal tolerance.

## 4. Interpretation

- **DERIVED COMPUTATIONAL COROLLARY:** under the complete anchoring
  hypotheses, no pair of independent endpoint functions accounts for the two
  schedule maps on any of the 16 frozen pure-momentum rays.
- **DERIVED:** the obstruction is in the final geometry `Q`, not merely the
  convention-dependent outgoing momentum.
- **STRUCTURAL:** the schedule maps have the same value and first derivative
  at the background but distinct quadratic configuration response.
- **NOT A FORCE LAW:** quadratic amplitude scaling diagnoses the first broken
  order; it is not a derived interaction or dispersion relation.
- **OPEN:** an unanchored endpoint translation that redefines the base
  momenta, a general canonical transformation mixing `q,p`, a nonlinear field
  redefinition, an improved/perfect action, refinement and the full carrier.

## 5. Framing boundary

Anchoring is load-bearing.  If arbitrary nonzero gradients at the base are
allowed, a source endpoint term translates all momenta at fixed `q0`; equation
(1) no longer follows.  Such a proposal changes the canonical interpretation
of the published base momenta and introduces up to 30 free first derivatives.
Without a geometric derivation it would recreate the fitting freedom this
programme is designed to avoid.  This result therefore does not chase an
unanchored fit.

Nor does the result select the even or odd schedule.  It says the bare
schedule dependence cannot be dismissed as the simplest base-preserving
boundary convention.  Both may still be discretization approximants rather
than physical alternatives.

## 6. Post-result prior-art status

The mechanism is standard variational/canonical structure:

- [*Canonical simplicial gravity*](https://arxiv.org/abs/1108.1974) derives
  pre/post evolution from the action;
- [*Constraint analysis for variational discrete
  systems*](https://arxiv.org/abs/1303.4294) supplies the constrained discrete
  setting;
- [*Multisymplectic Geometry, Variational Integrators, and Nonlinear
  PDEs*](https://authors.library.caltech.edu/records/74nvs-vb440) derives
  canonical boundary forms from action variation;
- [*On the exact discrete Lagrangian function for variational integrators:
  theory and applications*](https://arxiv.org/abs/1608.01586) emphasizes that
  an exact discrete Lagrangian is the object controlling comparison with the
  continuum trajectory.

No primary source located this exact 16-ray corollary for the dust 600-cell
schedule pair.  A search cannot prove novelty; external novelty is **OPEN**.

## 7. Consequence and next route

The cheap explanations are now bounded as follows:

```text
ordinary common Dirichlet principal function: invalid without constraints;
anchored endpoint counterterm:               refuted 16/16;
unanchored fitted momentum translation:      mathematically open, physically
                                               unselected;
general canonical equivalence:               open;
improved/perfect action or refinement:        open and now the clean route.
```

The next useful calculation should not fit 30 counterterm derivatives.  It
should ask whether a geometry-defined coarse-graining, schedule sum, or common
refinement supplies a unique effective action and whether its schedule defect
decreases.  Until that succeeds, the quotient does not define a unique
nonlinear physical tick.
