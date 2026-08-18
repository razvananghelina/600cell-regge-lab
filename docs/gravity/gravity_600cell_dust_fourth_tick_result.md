# Out-of-sample fourth weak-lapse dust tick

Date: 2026-08-16

## Provenance

- prior-art gate: `40d77c7`;
- frozen prediction protocol: `9eae63b`;
- registered implementation, committed before the first fourth-slab
  evaluation: `0e650aa`;
- targeted verifier:
  `reproducible/verify_gravity_600cell_dust_fourth_tick.py`;
- artifact: `reproducible/gravity_600cell_dust_fourth_tick.json`;
- artifact SHA-256:
  `4d8d03957675a6f454c1ad05102ffd1711f48c2e5a19f09b2898a60c9f07020d`.

Only this targeted verifier was run.  It returned **5/5**.  The full suite was
not run.

The calculation did not evaluate a `lambda=1` fourth state.  Its only fourth
step targets were the five integers committed in the protocol before the new
states existed.

## Verdict

```text
FOURTH_TICK_WEAK_LAPSE_PREDICTION_CONFIRMED
```

**DERIVED COMPUTATIONAL:** the fourth canonical seam exists on all three
weak-lapse controls and both independently derived schedule parities.  All six
new solves pass the complete action, momentum, Lorentzian-branch and calibrated
Jacobian gates.

The five predictions inferred from `n<=3` are confirmed at the new iteration
index:

```text
u4/u1          -> 4,
a4/u1          -> 10,
v4/v1          -> 7,
r4/v1          -> 16,
p_post,4/k     -> 9.
```

Every error decreases under both lapse halvings, all ten observed convergence
orders are within `[1.999993,2.000004]`, and every fine Richardson intercept
passes both its internal consistency gate and the fixed training band

```text
B_train = 4.6222921056804246599831556548181231e-10.
```

The band was committed before tick four and was not enlarged after seeing the
result.

## Numerical certificate

The value columns correspond to `lambda=(1/2,1/4,1/8)`; even and odd schedules
agree far inside the frozen tolerances.

| observable | target | computed values | observed orders | fine Richardson intercept |
|---|---:|---|---|---:|
| `u4/u1` | 4 | `4.000006682172`, `4.000001670541`, `4.000000417635` | `2.0000019`, `2.0000005` | `3.9999999999998145` |
| `a4/u1` | 10 | `10.000010023255`, `10.000002505811`, `10.000000626453` | `2.0000016`, `2.0000004` | `9.9999999999997639` |
| `v4/v1` | 7 | `7.000014761441`, `7.000003690352`, `7.000000922588` | `2.0000031`, `2.0000008` | `6.9999999999993388` |
| `r4/v1` | 16 | `16.000021087763`, `16.000005271931`, `16.000001317982` | `2.0000026`, `2.0000006` | `15.9999999999992700` |
| `p_post,4/k` | 9 | `8.999879824473`, `8.999969955989`, `8.999992488989` | `1.9999938`, `1.9999984` | `8.9999999999892271` |

The accepted even-parity absolute states are

```text
lambda=1/2:
  a4 = -7.790157397351537991742103038533878603853e-6,
  r4 = -1.423703946988178813286660848777531984116e-5;

lambda=1/4:
  a4 = -1.947537925747164658479963798680486534415e-6,
  r4 = -3.559256859498406998682651740151786502448e-6;

lambda=1/8:
  a4 = -4.868843924625024332471111287329187825080e-7,
  r4 = -8.898140268768103332891158857465487331030e-7.
```

Newton required at most three undamped accepted steps.  Endpoint reduced
residuals range from `3.85e-28` to `8.14e-42`; every calibrated weak Jacobian
direction is resolved.

## Consolidated four-step law

**DERIVED COMPUTATIONAL:** through the out-of-sample index `n=4`, the canonical
map has the leading weak-lapse recurrence

```text
p_post,n/k       = 2n+1       + O(lambda^2),
u_n/u1           = n          + O(lambda^2),
v_n/v1           = 2n-1       + O(lambda^2),
a_n/u1           = n(n+1)/2   + O(lambda^2),
r_n/v1           = n^2        + O(lambda^2).
```

The fourth step is genuine out-of-sample evidence in `n`: it was not used to
identify the recurrence or its integers.

**STRUCTURAL:** this is the discrete constant-leading-force / constant-leading-
acceleration pattern expected of a coherent local variational map.  The result
shows that the accepted slabs are one dynamics rather than three isolated
solutions.  It is not, by itself, a new law of gravity.

**OPEN:** an analytic derivation from the reduced action, validity for arbitrary
`n`, spatial refinement, anisotropic perturbations, continuum Einstein
convergence and emergent time.

## Post-result primary-source audit

The new search again located the established frameworks for canonical Regge
dynamics and variational-order analysis, but no primary source giving these
five fourth-step limits for the present fixed action:

- <https://arxiv.org/abs/gr-qc/0511096>;
- <https://arxiv.org/abs/1108.1974>;
- <https://arxiv.org/abs/gr-qc/0009093>;
- <https://arxiv.org/abs/gr-qc/0106077>;
- <https://arxiv.org/abs/1102.2685>.

External novelty remains **OPEN**; search is not proof.

The lapse warning is load-bearing.  Curved Regge configurations generically
break exact diffeomorphism symmetry and replace constraints by
pseudo-constraints (<https://arxiv.org/abs/0905.1670> and
<https://arxiv.org/abs/0909.5688>).  Therefore a lapse selected after symmetry
breaking is not automatically a physical fundamental clock; it may be a
discretization-dependent consistency condition.

## Correct next physical mission

Do not accumulate more homogeneous ticks.  The homogeneous recurrence is now
validated far enough to serve as a background.  The next discriminating object
is the linearized anisotropic canonical map around that background.

Before calculating it, preregister:

1. the full boundary perturbation space and its scale, gauge/pseudo-gauge and
   shape decomposition;
2. the Hessian/symplectic operator whose eigenmodes actually evolve;
3. the 600-cell symmetry irreducible sectors and expected degeneracies;
4. stability criteria and a comparison with tensor harmonics on continuous
   `S^3`;
5. a refinement test separating carrier artifacts from persistent physics.

A continuous compact `S^3` already has a discrete spectrum, a gap and symmetry
degeneracies.  Merely finding those features on the 600-cell would not be new.
The potentially informative signal is a symmetry-controlled deviation with a
documented refinement law, or an instability/refutation.

Likewise, a limiting speed cannot be defined merely as spatial edge length
divided by the selected lapse.  It would require a derived causal propagation
cone, invariance under lapse convention and a stable refinement limit.  Until
those tests pass, a physical tick and `c` remain **OPEN**.
