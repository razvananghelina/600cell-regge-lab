# Result: the action-selected conformal/shape split is dynamically closed

Date: 2026-08-18

## Headline

**DERIVED COMPUTATIONAL, with post-result falsifiability controls.**  On the
fixed regular 600-cell, along the first two accepted dust-Regge slabs, the
complete `720`-dimensional position recurrence has an action-selected reducing
decomposition

```text
720 edge-metric perturbations
= 120 vertex-conformal perturbations
  direct-sum
  600 kinetic-orthogonal shape perturbations.
```

Both normalized centered operators `Gamma=M^-1 N` and `Omega=M^-1 V`
preserve both factors.  All `224/224` primary residuals are
`ZERO_CONSISTENT`.  In the same arithmetic, all `448/448` declared alternative
split controls are `NONZERO_RESOLVED`, so the primary zero is selective rather
than a consequence of broad tolerances or dimension counting.

This produces an autonomous finite shape recurrence.  It does **not** yet
identify its `600` position directions with transverse-traceless gravitons or
reduce them to two physical polarizations.

## Provenance ledger

| stage | commit |
|---|---|
| prior-art and framing gate | `42313fb` |
| blind primary protocol | `4039244` |
| registered primary verifier before first execution | `3edfa46` |
| preserved first `224/224` artifact | `9bab7a3` |
| post-result power-control preregistration | `6a8ce90` |
| power-control implementation before execution | `4ed3009` |
| powered result artifact | `a232dbb` |

Targeted verifier:

```text
reproducible/verify_gravity_600cell_dust_conformal_shape_dynamics.py
```

Certified artifact:

```text
reproducible/gravity_600cell_dust_conformal_shape_dynamics.json
SHA-256 c5bbeaa2a64d07688061bc5098a33361dc2f5300d637e44a10b6cccbbd1bb162
```

Two complete powered runs returned byte-identical artifacts and

```text
12/12 PASS
CONFORMAL_SHAPE_DYNAMICS_DECOUPLED_POWER_CERTIFIED.
```

Only this mission-specific verifier and its direct imported geometry controls
were run.  The full suite was not run.

## Complete hypotheses

The claim is conditional on all of the following:

1. the fixed regular 600-cell spatial carrier and its literal `720` edge
   coordinates `delta log ell_e^2`;
2. the accepted first two nonstationary fixed-mass dust-Regge slabs and their
   certified Lorentzian branch;
3. the exact adjacent-slice edge identification used by the committed
   three-slice Jacobi equation;
4. the vertex-conformal incidence map

   ```text
   (C sigma)_uv = sigma_u + sigma_v;
   ```

5. the Hermitian kinetic form `H=(M+M*)/2` and its action-relative complement

   ```text
   K = im C,
   S_H = ker(C*H);
   ```

6. both frozen staircase schedules, all seven minimal binary-tetrahedral
   sectors and all four derivative variants;
7. the stored ball/roundoff envelopes and the preregistered `10/100`
   classification bands.

It is not a theorem for arbitrary triangulations, arbitrary Regge solutions,
later ticks, nonlinear perturbations or a continuum limit.

## Carrier controls

The calculation independently retains:

- `120` vertices, `720` edges, degree `12`, and `1,200` triangles;
- exact equivariance of `C` under all `24` frozen group actions;
- exact conformal rank `120`, proved by connectedness and an odd cycle;
- sector conformal dimensions `5d` for
  `d=3,2,2,2,1,1,1`;
- `56` resolved complements of dimensions `25d`;
- positive kinetic restriction on the conformal factor and negative
  restriction on the shape factor, in the frozen action convention;
- direct-sum condition numbers only `1.1881 ... 1.8065`;
- shape-row singular gaps `45.5355 ... 124.8248`;
- no open carrier rank or subspace.

Changing the overall sign of the action swaps the words positive and negative
but does not alter the decomposition or closure statement.

## Primary closure census

The normalized recurrence is

```text
delta^2 q + Gamma delta q + Omega q = 0.
```

For each normalized operator, conformal invariance was tested by

```text
||(I-P_K) X U_K||_2,
```

and shape invariance by

```text
||C* H X U_S||_2.
```

The complete ledger is:

| operator | carrier | zero-consistent | open | resolved nonzero |
|---|---|---:|---:|---:|
| `Gamma` | conformal | 56 | 0 | 0 |
| `Gamma` | shape | 56 | 0 | 0 |
| `Omega` | conformal | 56 | 0 | 0 |
| `Omega` | shape | 56 | 0 | 0 |
| **total** | — | **224** | **0** | **0** |

Raw norms and their ratios to the complete errors are:

| operator/carrier | norm range | error-unit range |
|---|---:|---:|
| `Gamma`, conformal | `4.61e-8 ... 1.71e-7` | `0.0586 ... 4.0228` |
| `Gamma`, shape | `3.50e-6 ... 2.55e-5` | `0.0484 ... 1.9988` |
| `Omega`, conformal | `6.15e-8 ... 2.28e-7` | `0.0568 ... 3.6577` |
| `Omega`, shape | `4.67e-6 ... 3.40e-5` | `0.0462 ... 1.8351` |

The zero threshold is `10` error units.  Therefore the weakest primary margin
is not an exact symbolic zero, but every audit remains on the zero side of the
unchanged calibrated rule.

All `112/112` even/odd residual-norm comparisons are `SCHEDULE_ROBUST`; their
largest distance/error ratio is `2.03e-9`.

## Why the positive is not automatic

The first primary positive run was intentionally not accepted on faith.  After
preserving it, two alternative splits of identical dimensions were frozen:

1. the Euclidean positive/negative spectral split of `H`, already known to be
   distinct from the conformal carrier;
2. a deterministic dense Fourier coordinate split.

Both factors of both controls were tested for both normalized operators in
all `112` schedule/sector/variant/operator cells:

| control factor | resolved nonzero | minimum error units | maximum error units |
|---|---:|---:|---:|
| spectral positive | 112/112 | `7.26e4` | `4.92e6` |
| spectral negative | 112/112 | `8.39e4` | `2.66e6` |
| Fourier low | 112/112 | `1.38e5` | `1.42e7` |
| Fourier high | 112/112 | `1.55e5` | `1.40e7` |

Thus every one of the `112` power cells is `POWER_HIT`.  The same classifier
which assigns the physical split at most `4.03` error units assigns every
control at least `72,602` error units.  Neither group equivariance nor equal
dimension can explain the primary zero: all these subspaces live inside the
same symmetry-reduced blocks, where each irrep has multiplicity.

The control was preregistered only after the first positive was known, so it
is labelled honestly as a post-result falsifiability audit.  It strengthens
the numerical selectivity but does not retroactively make the original result
cognitively blind.

## Algebraic consequence

Within the calibrated finite calculation,

```text
Gamma K subset K,       Omega K subset K,
Gamma S_H subset S_H,   Omega S_H subset S_H.
```

Therefore, in the action-selected basis `[K,S_H]`, both normalized operators
are block diagonal and the centered equation separates into two autonomous
equations.  On the doubled first-order phase carrier this induces invariant
dimensions

```text
240 conformal phase directions
1200 shape phase directions.
```

This is stronger than the previous invariant global-scale plane.  The old
order-24 quotient separated only one uniform scale coordinate from `29`
zero-sum edge-orbit coordinates.  It did not cover the other four conformal
vertex-orbit directions or the full `720`-edge carrier.

## Reconciliation with the failed rigidity/self-stress route

There is no contradiction with
`RIGIDITY_YORK_DECOUPLING_REFUTED`.

The failed carrier was

```text
ker R*,
```

the Euclidean self-stress complement of the embedded bar-and-joint rigidity
map.  It is selected by geometry alone and has dimension `250`.  Hundreds of
its Regge dynamic cross blocks are resolved nonzero.

The successful carrier is

```text
ker(C*H),
```

the complement of vertex conformal variations selected jointly by geometry
and the action-derived kinetic bilinear form.  It has dimension `600` and is
not a framework self-stress space.

This is exactly the conceptual distinction the two hostile tests were meant
to expose: a plausible geometry-only complement is not enough; the dynamics
selects the relevant orthogonality.

## Physical meaning and limits

What is now supported:

- **DERIVED COMPUTATIONAL:** one conformal plus five shape directions per
  vertex in the kinetic signature;
- **DERIVED COMPUTATIONAL:** the corresponding complete conformal/shape
  position split is dynamically reducing for the centered recurrence;
- **DERIVED COMPUTATIONAL:** the result is schedule robust and selectively
  distinguishes four alternative control factors per audit;
- **STRUCTURAL:** this is a finite DeWitt-like conformal/shape decomposition
  selected by the declared Regge action and time-fibre identification.

What is not supported:

- `600` shape directions are not `600` gravitons;
- no Hamiltonian/momentum constraint quotient has produced two tensor
  polarizations per spatial point;
- the shape factor may still contain longitudinal, scalar, vector,
  constraint-violating and discretization modes;
- frozen dust has not been promoted to independent perturbation variables;
- there is no refinement theorem, continuum dispersion, physical clock,
  limiting speed, Planck scale or particle-mass formula.

The positive therefore opens a **well-defined subsystem**, not the complete
physics gate.

## Post-result prior-art check

The second search used `simplicial conformal mode dynamics`, `Regge conformal
traceless decoupling`, `Lund-Regge invariant subspace` and `lattice conformal
gauge`.

- Glickenstein's discrete conformal map is established prior art:
  <https://arxiv.org/abs/0906.1560>.
- Brown treats continuum conformal--traceless variables together with their
  evolution and gauge conditions: <https://arxiv.org/abs/gr-qc/0501092>.
- Hamber and Williams discuss lattice conformal/diffeomorphism separation,
  mainly in weak-field and two-dimensional settings:
  <https://arxiv.org/abs/hep-th/9607153>.
- Catterall--Mottola--Bhattacharya study conformal-mode dynamics in simplicial
  quantum gravity, not this Lorentzian Regge transfer operator:
  <https://arxiv.org/abs/hep-lat/9809114>.
- Bahr--Dittrich's curved-background warning remains applicable:
  <https://arxiv.org/abs/0905.1670>.

No located primary source gives the present fixed dust-600-cell reducing
subspace, the `224/224` closure census or its powered controls.  Search failure
does not prove novelty; external novelty remains **OPEN**.

## Next load-bearing gate

The next calculation should remain entirely inside the now certified
`600`-dimensional shape recurrence.  It must not impose a gauge quotient by
hand.  The correct order is:

1. construct the action-derived weak/pseudo-constraint response restricted to
   the shape factor;
2. determine whether it selects a smaller invariant null/weak carrier;
3. test curvature response and two-step transport of that carrier;
4. only then compare the surviving generalized stiffness with independently
   derived spatial tensor harmonics;
5. require the comparison and low-mode dispersion to survive a declared
   refinement family before discussing an effective `c`.

The present result makes that sequence finite and non-arbitrary.  It does not
prejudge whether the remaining physical quotient exists.
