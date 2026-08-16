# Second canonical homothetic tick: direct-predictor result

Date: 2026-08-16

## Provenance

- prior-art gate: `1865e13`;
- frozen protocol: `fb119d1`;
- implementation and registration before evaluation: `1cc52a4`;
- artifact SHA-256:
  `51a674ad2ade10ade2cfccc6a593805d7554b49612789fc7443cc46587414768`.

Only the targeted verifier was run.  It returned **8/8**.  The full suite was
not run.

## Mechanical verdict

```text
SECOND_TICK_NEWTON_OPEN
```

This is neither an accepted second tick nor a no-go theorem.  The sole frozen
constant-increment predictor stays Lorentzian and every calibrated Jacobian
has resolved rank two, but the deterministic eight-iteration Newton solve
does not reach its residual gate.

## Seed and trajectory

The predictor begins at

```text
u = -3.1160595767e-6,
v = 0,
F = (-3.3956483e-8, -1.8162085e-3).
```

Thus carrying the first lapse forward is already far from the inherited
canonical momentum.  The first two undamped corrections move to approximately

```text
(u,v)=(-4.36250e-6,-0.800007),
(u,v)=(-3.32889e-6,-0.999243).
```

The Armijo rule then accepts six damped steps and stops at the frozen
iteration limit with

```text
u = -3.820550117259453e-6,
v = -0.768315571799095,

L2/L1       = 0.999996179457181,
rho2/rho1   = 0.463793637802952,
tau2/tau1   = 0.681023962723011,

F = (1.7524869e-10, -7.2955682e-5),
||p_pre,2-target||_2 = 3.9959473e-4.
```

The allowed junction bound is `3.6513654e-21`, so the endpoint misses the
canonical seam by about seventeen orders of magnitude.  It is not an
approximate accepted tick.

Both schedule parities reproduce the same trajectory to the stored precision.

## Weak direction

The smallest calibrated singular value changes along the attempted path:

```text
seed              2.1222817e-8,
after step 1       4.1960991e-9,
minimum observed  7.8539492e-10,
final attempt      1.1567290e-9.
```

All are resolved against calibrated errors of order `1e-22`; no exact rank
loss was observed.  The large lapse motion and softening nevertheless support
the **STRUCTURAL / candidate pseudo-constraint** warning.  They do not prove
that the lapse tends to zero or becomes gauge.

## Interpretation ledger

- **DERIVED COMPUTATIONAL NEGATIVE:** the preregistered direct
  constant-increment Newton procedure does not produce a second tick.
- **DERIVED COMPUTATIONAL:** the failure is not caused by the Lorentzian
  branch, schedule parity or unresolved Jacobian rank.
- **PATTERN:** the inherited momentum appears to require a substantially
  shorter second lapse.  No physical claim is licensed because no root was
  reached.
- **OPEN:** existence of a connected second-tick root.  Eight Newton
  iterations are an algorithmic boundary, not a mathematical nonexistence
  proof.
- **NOT DERIVED:** continued contraction, a recollapse history, emergent time
  or any absolute clock.

## Post-result primary-source check

The refined searches used `canonical Regge evolution lapse collapse`,
`pseudo-constraint next time step` and `600-cell dust numerical evolution`.
The located literature continues to support only the general interpretation:
canonical simplicial moves propagate pre/post data
(Dittrich--Hoehn, <https://arxiv.org/abs/1108.1974>), while curved Regge
backgrounds can replace exact lapse gauge by pseudo-constraints
(Bahr--Dittrich, <https://arxiv.org/abs/0905.1670>).  Earlier 600-cell
evolutions encounter finite-step stopping/causality difficulties
(De Felice--Fabri, <https://arxiv.org/abs/gr-qc/0009093>,
<https://arxiv.org/abs/gr-qc/0106077>), but this does not explain or prove the
present local solver failure.

No located primary source supplies this exact second-tick equation or root.
External novelty remains **OPEN**.

## Next admissible attack

Do not simply increase the iteration limit after seeing the miss.  A separate
preregistered target homotopy can test connected existence without using the
failed endpoint: continuously interpolate the lower boundary and incoming
canonical target from the already accepted first-tick problem to the desired
second-tick problem, starting at the exact accepted root.  Fixed homotopy
steps and the same calibrated Newton gates can distinguish a connected branch
from another direct-solver accident.
