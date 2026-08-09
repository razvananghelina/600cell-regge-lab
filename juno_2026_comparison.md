# JUNO consistency check and the `1/45` data-driven scoping decision

Date: 2026-07-22

## Provenance correction (2026-07-24)

This document originally framed the exercise as a prediction confronting new
data.  That was wrong.  JUNO arXiv:2511.14593 was submitted on 2025-11-18
and NuFIT 6.0 arXiv:2410.05380 on 2024-10-07.  The repository's neutrino
formulae were fixed only in February 2026; `prompt_neutrino_masses.md`
explicitly supplies the already measured splittings.  Consequently every
agreement below involving `sin^2(theta_12)`, `sin^2(theta_13)`,
`delta_CP`, or a mass splitting is a **RETRODICTION / CONSISTENCY CHECK**,
not a successful prediction.  Variant I remains a **PATTERN**, and its
selection is explicitly data-driven.

## Primary experimental inputs

The inputs were independently checked online on 2026-07-22 and rechecked
against the linked arXiv records on 2026-07-24.  JUNO, DESI, the
dynamical-dark-energy analysis, and KATRIN match the pinned values.  The
NuFIT arXiv record was reachable, but the official parameter-table PDF was
not retrievable during the 2026-07-24 audit; its table entries therefore
remain pinned to the earlier check rather than newly independently verified.
Values are
pinned in `reproducible/verify_juno_comparison.py`.

- JUNO Collaboration, [First measurement of reactor neutrino oscillations at
  JUNO](https://arxiv.org/abs/2511.14593), first 59.1 days:
  `sin^2(theta_12)=0.3092+-0.0087` and
  `Delta m^2_21=(7.50+-0.12) 10^-5 eV^2` for normal ordering.
- NuFIT 6.0, [arXiv:2410.05380](https://arxiv.org/abs/2410.05380) and the
  [official parameter table](https://www.nu-fit.org/sites/default/files/v60.tbl-parameters.pdf),
  normal ordering:
  - without SK atmospheric data,
    `Delta m^2_31=2.534^(+0.025)_(-0.023) 10^-3 eV^2`;
  - with SK atmospheric data,
    `Delta m^2_31=2.513^(+0.021)_(-0.019) 10^-3 eV^2`,
    `sin^2(theta_13)=0.02215^(+0.00056)_(-0.00058)`,
    `theta_13=8.56+-0.11 deg`, and
    `delta_CP=212^(+26)_(-41) deg`.
- DESI Collaboration, [arXiv:2503.14744](https://arxiv.org/abs/2503.14744),
  DR2 BAO plus DR1 full shape in LambdaCDM:
  `sum m_nu < 0.0642 eV` at 95% confidence.
- The model-dependent dynamical-DE analysis
  [arXiv:2507.16589](https://arxiv.org/abs/2507.16589), combining DESI DR2,
  CMB, DESY5 and DESY1, reports
  `sum m_nu=0.098^(+0.016)_(-0.037) eV` and a 2.7-sigma positive preference
  in `w0 wa CDM`.
- KATRIN Collaboration, [arXiv:2406.13516](https://arxiv.org/abs/2406.13516)
  and Science 388 (2025), gives `m_beta<0.45 eV` at 90% confidence.

**EXTERNAL discrepancy found:** JUNO, DESI and KATRIN agree with the mission
inputs.  The paper's old atmospheric entry
`Delta m^2_32=(2.453+-0.033) 10^-3 eV^2` is not the official NuFIT 6.0
normal-ordering table entry.  NuFIT reports `Delta m^2_31`, with the two
values above.  The comparison below uses both official fits.  Where
`Delta m^2_32` is displayed, its experimental central value is inferred as
`Delta m^2_31-Delta m^2_21`; errors are combined in quadrature without
correlations and are explicitly secondary to the direct `Delta m^2_31`
comparison.

## Theory definitions and variants

The exact framework inputs are

`phi=(1+sqrt(5))/2`,

`alpha=(20 phi^4-sqrt((20 phi^4)^2-8 pi))/(4 pi)`,

`r=alpha phi^3=0.030912038750648...`, and `m_1=0`.

The bare and corrected heavy masses are

`m3_b=2m_e/phi^35=49.531627 meV`,

`m3_c=2m_e/phi^(35-1/45)=50.064140 meV`,

with `m2_b=m3_b sqrt(r)=8.708561 meV` and
`m2_c=m3_c sqrt(r)=8.802187 meV`.

- Variant I: `(m2,m3)=(m2_b,m3_c)`; eigenvalue-local correction.
- Variant II: `(m2,m3)=(m2_c,m3_c)`; reimpose the ratio after correction.
- Variant III: `(m2,m3)=(m2_b,m3_b)`; bare package.

## Complete numerical comparison

Signed sigma values show `(framework retrodiction-experiment)/sigma`.  For asymmetric
NuFIT errors, the error on the prediction's side of the central value is
used.  Limits have no Gaussian sigma and are reported as margins.  The
dynamical-DE sigma uses its lower `0.037 eV` uncertainty because all
predictions lie below its central value.

| Observable | Variant I | Variant II | Variant III | Experimental comparison |
|---|---:|---:|---:|---|
| `m2` [meV] | 8.708561 | 8.802187 | 8.708561 | no direct measurement |
| `m3` [meV] | 50.064140 | 50.064140 | 49.531627 | no direct measurement |
| `Delta m2_21` [`eV^2`] | `7.583904e-5` (`+0.699 sigma`) | `7.747850e-5` (`+2.065 sigma`) | `7.583904e-5` (`+0.699 sigma`) | JUNO `7.50+-0.12 e-5` |
| `Delta m2_31` [`eV^2`] | `2.506418e-3` (`-0.346 sigma` SK; `-1.199 sigma` no-SK) | same | `2.453382e-3` (`-3.138 sigma` SK; `-3.505 sigma` no-SK) | NuFIT 6.0 official fits |
| `Delta m2_32` [`eV^2`] | `2.430579e-3` (`-0.390 sigma` SK; `-1.234 sigma` no-SK) | `2.428940e-3` (`-0.476 sigma`; `-1.305 sigma`) | `2.377543e-3` (`-3.176 sigma`; `-3.537 sigma`) | inferred using JUNO; correlations neglected |
| `sum m_nu` [meV] | 58.772702 | 58.866327 | 58.240189 | all below 64.2 meV LCDM limit by 5.43, 5.33, 5.96 meV |
| deviation from dynamical-DE preference | `-1.060 sigma` | `-1.058 sigma` | `-1.075 sigma` | model-dependent `98^(+16)_(-37)` meV |
| `m_beta` [meV] | 8.837860 | 8.865225 | 8.770929 | all far below KATRIN 450 meV bound; sigma N/A |
| `m_betabeta` [meV] | 3.103043 | 3.129064 | 3.095782 | no measurement; sigma N/A |

The two-Gaussian diagnostic using JUNO `Delta m^2_21` and NuFIT-with-SK
`Delta m^2_31` is

| Variant | `chi^2` for two inputs |
|---|---:|
| I | 0.6089 |
| II | 4.3859 |
| III | 10.3346 |

This diagnostic ignores correlations between experiments.  It ranks the
variants but is not a model-selection significance.

## Mixing observables common to all variants

| Observable | Theory | Experiment | Signed deviation |
|---|---:|---:|---:|
| `sin^2(theta_12)` | `2/(phi+5)=0.302204553` | JUNO `0.3092+-0.0087` | `-0.804 sigma` |
| `sin^2(theta_13)` | `1/45=0.022222222` | NuFIT `0.02215^(+0.00056)_(-0.00058)` | `+0.129 sigma` |
| `theta_13` | `8.573105 deg` | NuFIT `8.56+-0.11 deg` | `+0.119 sigma` |
| `delta_CP` | `3 atan(sqrt(5))=197.715472 deg` | NuFIT `212^(+26)_(-41) deg` | `-0.348 sigma` |

## Scoping decision

### What is derived

- **DERIVED:** the bare formula `m3_b=2m_e/phi^35`.
- **DERIVED within the stated framework chain:** the bare spectral ratio
  `Delta m^2_21/Delta m^2_31=alpha phi^3`.
- **PATTERN:** `m3_b -> m3_c=m3_b phi^(1/45)`.  The repository explicitly
  says that its detailed spectral-action mechanism has not been computed.

### What the derivation does not decide

The available kernel-zeta/inner-fluctuation argument derives the ratio for the
blind bare package.  It does not compute a perturbation operator for the
`nu_3-nu_e` correction and does not prove whether that perturbation commutes
with, precedes, or follows the ratio construction.

Therefore:

- applying the correction only to `m3` (I) is not ruled out by the ratio
  theorem, if the theorem is read as fixing the bare spectrum before a
  state-specific mixing correction;
- propagating it to `m2` (II) is not forced, but is the consistent choice if
  the ratio is postulated to hold for the corrected physical eigenvalues;
- omitting it (III) retains the fully bare package but now has a roughly
  `3.1 sigma` atmospheric tension against the official NuFIT-with-SK fit.

### Decision

- **Most physically motivated:** Variant I.  A leading `nu_3-nu_e` mixing
  correction is naturally local to the `m3` eigenvalue and has no computed
  mechanism acting on `m2`.
- **Most consistent with the pinned splitting data:** Variant I, with
  `chi^2=0.609` for the two displayed Gaussian inputs.
- **PATTERN:** this scoping choice remains PATTERN, because the finite
  perturbation has not been derived.  Data preference cannot upgrade it.
- **DERIVED empirical statement:** JUNO now disfavors the propagated Variant
  II relative to I/III in the solar channel (`2.065 sigma` versus
  `0.699 sigma`).  NuFIT atmospheric data independently favor corrected `m3`
  (I/II) over bare III.

## Blind falsification ledger

- A robust LambdaCDM cosmological bound below `58.24 meV` excludes all three
  variants under `m1=0`; a DESI DR3 bound below about `58 meV` therefore kills
  the package in LambdaCDM.  Cosmological-model dependence must be stated.
- Future splitting measurements can reject formula variants, but these are
  not blind predictive assets because their targets predated construction.
- KATRIN's present `0.45 eV` bound is not discriminating.  A direct bound or
  measurement below the predicted `8.77--8.87 meV` range would exclude the
  package, but that is far beyond current KATRIN sensitivity.
- A reliable positive `m_betabeta` measurement incompatible with roughly
  `3.10--3.13 meV`, under the assumed Majorana phases and nuclear matrix
  elements, rejects that phase-specific prediction.
- A nonzero lightest mass excludes `m1=0`; an inverted-ordering
  determination excludes the strict normal-ordering claim.

## Classification summary

- **DERIVED:** exact three-variant observables and deviations from pinned
  inputs; bare package; current empirical ranking.
- **STRUCTURAL:** reading the ratio theorem as a bare spectral relation before
  later eigenstate mixing.
- **PATTERN:** the `1/45` mass correction and the Variant-I versus Variant-II
  scoping.
- **OPEN:** the actual finite perturbation operator and whether it preserves
  the ratio relation on corrected masses.
