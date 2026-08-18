# Second homothetic dust tick: accepted local canonical correction

Date: 2026-08-16

## Provenance

- target-independent stationary roots: `caaf1f1`;
- prior-art and fixed-mass framing gate: `fcc4d7c`;
- stationary target comparison: `1d32334`;
- frozen local-correction protocol: `54dd336`;
- registered implementation before evaluation: `7b3552d`;
- targeted verifier:
  `reproducible/verify_gravity_600cell_dust_second_tick_local_correction.py`;
- artifact:
  `reproducible/gravity_600cell_dust_second_tick_local_correction.json`;
- artifact SHA-256:
  `936984bc84a714140ce16917ee559b346b3c0d4a5ba92d8fb723398a120f8e70`.

Only this targeted verifier was run.  It returned **6/6**.  The full suite was
not run.

## Verdict

```text
SECOND_HOMOTHETIC_TICK_ACCEPTED
```

**DERIVED COMPUTATIONAL:** at fixed conserved dust mass, the contracting
stationary root has a unique resolved local correction of absolute upper
scale and lapse that simultaneously satisfies the pole equation and the
canonical momentum seam.

Both independently built order-24 schedule parities give

```text
b2 = log(L2/L0)
   = -9.3481870589058271363382229926575337302742819401e-6,

r2 = log(rho2/rho0)
   = -1.4237027552009802996130054524247481533837837066e-5.
```

Relative to the accepted first output,

```text
u2 = log(L2/L1)
   = -6.2321274822113254446035165484588942890526300145e-6,

v2 = log(rho2/rho1)
   = -1.0677774416839169558879749184608307568118093612e-5,

L2/L1       = 0.9999937678919374548098562766940532526983,
rho2/rho1   = 0.9999893222825903911752645241121734947739,
tau2/tau1   = 0.9999946611270434133643461998608157349719.
```

The scale continues to contract.  The proper duration changes only by about
5.34 parts per million; it does not collapse toward zero as in the failed
target homotopy.

## Full gates

For each parity:

```text
Newton iterations                  3,
reduced residual infinity norm     6.047e-44,
canonical junction norm            3.312e-43,
inherited junction bound            3.651e-21,
maximum diagonal residual          1.41e-99 even / 2.17e-96 odd,
maximum pole residual              1.97e-48.
```

The endpoint Jacobian has

```text
singular values = (5.8284599359e2, 4.2445618156e-9),
epsilon         = 1.3300997012e-22.
```

Thus the weak singular direction remains ill-conditioned but is separated
from the calibrated error by about `3.19e13`.  Both parities agree to at least
90 decimal orders in `(b2,r2)` and pass the complete Lorentzian branch tests.

## The new sequence and its status

Let the first accepted increments be

```text
u1=a1=-3.1160595766945e-6,
v1=r1=-3.5592531351706e-6.
```

The second increments satisfy

```text
u2/u1 = 2.00000267287005,
v2/v1 = 3.00000421754978,
b2/u1 = 3.00000267287005.
```

**PATTERN:** the absolute scale logs now look like the first triangular
numbers `(1,3)` times `u1`, while the absolute lapse logs look like the first
squares `(1,4)` times `v1`.  With only two accepted ticks this is not a law,
not an exact integer identity and not evidence of cosmological acceleration.
The third tick is the first calculation capable of falsifying this pattern.

## What is physically established

- **DERIVED COMPUTATIONAL:** two consecutive non-static homogeneous slabs can
  be glued canonically on the fixed carrier with one conserved dust mass.
- **DERIVED COMPUTATIONAL:** the geometry contracts on both steps; the second
  scale decrement is approximately twice the first.
- **DERIVED:** the exact static all-lapse cancellation is broken by unequal
  boundary scales at conserved mass, and the joint pole/seam equations select
  relative lapse corrections.
- **STRUCTURAL:** this is a discrete homogeneous/isotropic Regge evolution
  witness, i.e. a two-step minisuperspace trajectory.
- **OPEN:** arbitrary-tick recurrence, stability under anisotropic
  perturbations, refinement convergence, continuum Friedmann agreement,
  global branch uniqueness and the ultimate causal endpoint.
- **OPEN:** an absolute clock.  The starting `tau0=0.0102` remains externally
  supplied, so neither a time unit, `c`, Planck time nor Planck mass is derived.

## Post-result primary-source audit

Multiple-step 600-cell Regge evolution is not new in broad form.  Barrett et
al. give an implicit vertex-evolution scheme and illustrate it on 600-cell
Friedmann cosmology: <https://arxiv.org/abs/gr-qc/9411008>.  De Felice and
Fabri evolve a dust-filled 600-cell and study its causality-breaking endpoint:
<https://arxiv.org/abs/gr-qc/0009093>.  Dittrich, Gielen and Schander analyze
Lorentzian Regge shells of closed universes, including dust, and emphasize
comparison across discretizations: <https://arxiv.org/abs/2109.00875>.

No located primary source prints this particular fixed order-24 staircase,
target-independent root firewall, action-generated canonical seam and the
two numerical ratios above.  External novelty is **OPEN**; a search cannot
prove novelty.

## Next falsification test

Repeat the same provenance order for tick three:

1. at fixed lower geometry `(b2,r2)`, enumerate the complete stationary-root
   multiset without reading `p_post,2`;
2. commit that multiset;
3. compare every root with the mapped canonical target;
4. only then apply one frozen local correction to the geometrically forward
   contracting branch.

Preregister, but do not fit, the two-tick pattern predictions

```text
u3/u1 approximately 3,
absolute log(L3/L0)/u1 approximately 6,
v3/v1 approximately 5,
absolute r3/v1 approximately 9.
```

Failure of any prediction kills the integer-sequence pattern while leaving
the already accepted two-tick evolution intact.
