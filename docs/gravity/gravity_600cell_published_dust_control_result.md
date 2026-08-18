# Published 600-cell dust sandwich: complete-action reproduction

Date: 2026-08-13

Prior-art gate: `gravity_600cell_evolution_prior_art.md`, commit `e7d8bd5`

Frozen protocol: `gravity_600cell_published_dust_control_protocol.md`, commit
`cc0902b`

Recorded binary64 failure and frozen precision correction: commit `3056c7d`

Registered verifier:
`reproducible/verify_gravity_600cell_published_dust_control.py`

Machine-readable result:
`reproducible/gravity_600cell_published_dust_control.json`

Targeted run: **14/14 implementation checks passed**.  The full repository
suite was not run at the user's request.

## Headline

> **DERIVED EXTERNAL CONTROL:** the unrounded De Felice--Fabri
> time-symmetric 600-cell dust sandwich satisfies every one of the 35 internal
> orbit equations of the repository's complete 2400-simplex one-slab action,
> for both ordered-schedule parities.

This validates the action normalization, Lorentzian angle branch, dust source
sign and the dynamics engine against published physics.  It is a successful
reproduction, not a new physical tick: the source paper already reports that
its unused slanting/horizontal equations have tiny residuals.

## 1. Input fixed by the publication

The verifier reconstructs, without fitting,

```text
M_star = 10
zeta   = (pi^2 sqrt(2)/50)^(1/3)
R0     = 4 M_star/(3 pi)
l0     = zeta R0
M      = (90/pi)[2 pi-5 acos(1/3)] l0
tau    = 0.0102
d^2    = l0^2-tau^2.
```

Numerically,

```text
M    = 10.202069074391675
l0^2 = 7.69379990138305
d^2  = 7.69369586138305.
```

The maximum discrepancy from the decimals printed in the 2000 paper is
`7.64e-14`.  Both boundary 600-cells use `l0^2`; all five pole orbits use
`-tau^2`; all 30 staircase-diagonal orbits use `d^2`.

The action convention is the paper's action multiplied by `8*pi`:

```text
S_total = S_gravity - (8*pi*M/120) sum_(120 poles) tau_i.
```

For one 24-pole orbit this gives

```text
dS_dust/d(rho_k) = -4*pi*M/(5 sqrt(rho_k)).
```

The independent pole-length form reproduces the published right-hand side
`pi*M/15` with relative errors at most `4.31e-12` in binary64.

## 2. Complete carrier checks

At the published point, every one of the 2400 four-simplices is
nondegenerate Lorentzian.  The minimum absolute Gram eigenvalue is

```text
even: 1.040e-4
odd:  5.202e-5,
```

and the minimum angle-argument modulus is `0.9955`.  Thus stationarity is not
being obtained at a causal or branch singularity.

For each parity, the 100-simplex-orbit gravitational action and all 95
derivatives (65 variable plus 30 old-boundary) agree with the direct
2400-simplex evaluation.  The largest reported gradient discrepancies are
`4.62e-11` and `4.94e-11`; action discrepancies are `3.51e-9` and `3.77e-9`.

## 3. Why the first independent derivative failed

The original preregistered binary64 centered-difference checks failed and are
retained in the output:

| parity | relative derivative error | imaginary derivative |
|---|---:|---:|
| even | `1.271e-3` | `9.933e-6` |
| odd | `1.000` | `5.356e-2` |

This was recorded before correction in
`gravity_600cell_published_dust_control_float_failure.md`.  At
`rho=tau^2=0.00010404`, the frozen relative step `3e-6` is only
`3.1212e-10`.  Gravitational and dust orbit derivatives of order `10^3`
cancel to a stationary value; binary64 action noise divided by that small
step overwhelms the result.

The correction, committed before its first evaluation, retained the same
geometry, source, action, branch, step and thresholds and reimplemented the
action only at 60 decimal digits.  It did not reuse the analytic gradient.

## 4. Corrected independent result

The 60-decimal base action agrees with the certified binary64 orbit action to
`1.90e-11` (even) and `4.59e-10` (odd).  Its centered differences give:

| quantity | even | odd |
|---|---:|---:|
| error versus analytic total gradient | `1.106e-8` | `1.262e-8` |
| largest imaginary derivative | `1.04e-50` | `3.53e-47` |
| largest direct per-edge residual | `4.607e-10` | `4.607e-10` |

The frozen stationarity threshold is `1e-7`, so both parities pass by more
than two orders of magnitude.  The residual dust/gravity cancellation ratios
on the five poles range from `7.79e-17` to `2.94e-16`.

The direct 60-decimal residual summaries are identical for the two parities:

```text
maximum pole residual       3.0846e-14
maximum diagonal residual   4.6068e-10
35-component residual norm  2.4504e-9.
```

The diagonal value is the centered-difference truncation at the frozen step;
the separately certified analytic residuals are `1e-14`--`1e-13`.  All ten
phase-pair sums pass, but the stronger fact is that every individual orbit
passes.

## 5. Post-result prior-art correction

The source itself states that unused slanting and horizontal edge equations
had residuals below `2e-15`, and separately follows the unused timelike-edge
equation as a Bianchi/source diagnostic:
[The Friedmann universe of dust by Regge Calculus](https://arxiv.org/abs/gr-qc/0009093).
The earlier Barrett et al. paper had explicitly left verification of its
ignored full equations as future work:
[A Parallelizable Implicit Evolution Scheme for Regge Calculus](https://arxiv.org/abs/gr-qc/9411008).

Therefore the result is not evidence that we discovered a previously unknown
dust solution.  What is added is an independent, registered reconstruction
on the corrected five-phase carrier, with a complete 35-orbit table, both
schedule parities, full/reduced controls and high-precision action
differences.  External novelty of this packaging is not a physics claim.

## 6. What this establishes physically

- **DERIVED EXTERNAL CONTROL:** the code can reproduce a known classical
  Regge evolution datum with matter.
- **DERIVED:** the regular time-symmetric sandwich is insensitive to the two
  schedule-parity classes at the level of the complete stationary equations.
- **DERIVED:** the earlier zero-`Lambda`, matter-free regular no-go was not an
  action bug; adding the published dust term supplies exactly the missing
  pole balance.
- **NOT DERIVED:** matter, its mass or the tick duration from the 600-cell
  geometry.  All are external inputs in this control.
- **NOT DERIVED:** subsequent multi-tick evolution; only the initial
  time-symmetric sandwich is reproduced here.
- **NOT NEW PHYSICS:** the dust sandwich and small unused-equation residuals
  were already published.

## 7. Status ledger

| Claim | Status |
|---|---|
| Published constants reconstructed without fitting | **DERIVED** |
| Complete slab is Lorentzian and off all branch boundaries | **DERIVED COMPUTATIONAL** |
| Full 2400-simplex and 100-orbit actions/derivatives agree | **DERIVED COMPUTATIONAL** |
| Binary64 total-action differences certify stationarity | **DERIVED NUMERICAL FAILURE, retained** |
| Corrected 60-decimal action differences certify stationarity | **DERIVED COMPUTATIONAL** |
| All 35 individual equations pass in both parities | **DERIVED EXTERNAL CONTROL** |
| This dust tick is a new discovery | **REFUTED BY PRIOR ART** |
| The code is now calibrated on known 600-cell dynamics | **DERIVED** |
| The geometry selects `M`, `tau`, dust or a clock rate | **OPEN / NOT CLAIMED** |
| A generalized nonhomogeneous tick exists | **OPEN** |

## 8. Next valid step

The engine has now earned a controlled continuation.  The next protocol
should hold the published source and lapse fixed and perturb the final
boundary in a target-free basis selected before solving.  Use the rank-30
boundary response and solve all 35 internal equations, comparing the result
with the homogeneous published branch.

This asks a genuinely sharper question than repeating Friedmann evolution:
does the exact 600-cell action propagate small nonhomogeneous boundary data,
and how many physical versus gauge directions survive?  Before computing,
perform a dedicated prior-art search for linearized/nonhomogeneous
perturbations of the 600-cell Sorkin evolution and preregister the perturbation
basis.  Do not choose a direction because it happens to produce a root.
