# Prediction-provenance ledger

Date: 2026-07-24

## Rule and chronology

This ledger records **predictive provenance**, not the separate question of
whether a formula follows from assumptions inside the framework.  A formula
can be internally DERIVED and still be a **RETRODICTION** if the comparison
value was public when the formula was chosen.  In every unresolved case the
classification is against the theory:

- **RETRODICTION**: a relevant measured value was public before the formula's
  earliest evidenced fixation;
- **BLIND**: only a bound, or no measurement, was public at fixation and the
  stated point outcome remains unmeasured;
- **AMBIGUOUS**: the repository does not establish the order of fixation and
  publication.

The decisive neutrino record is unusually clear.  The older
`600cell_theory_v3.6_backup.tex` says neutrino masses were not addressed.
`masa_neutrinilor.md` (filesystem/repository date 2026-02-11) explores a
different exponent 34.  `prompt_neutrino_masses.md` (2026-02-23) begins by
listing the already measured splittings and already contains the PMNS
formulae.  The present exponent-35/splitting package was therefore fixed no
earlier than February 2026.  JUNO arXiv:2511.14593 was submitted
2025-11-18 and NuFIT 6.0 arXiv:2410.05380 on 2024-10-07.  The repository
contains no immutable pre-result registration.

Primary publication records used below: [JUNO
2511.14593](https://arxiv.org/abs/2511.14593), [NuFIT 6.0
2410.05380](https://arxiv.org/abs/2410.05380), [DESI
2503.14744](https://arxiv.org/abs/2503.14744), [KATRIN
2406.13516](https://arxiv.org/abs/2406.13516), and [PDG 2024
CKM review](https://pdg.lbl.gov/2024/reviews/rpp2024-rev-ckm-matrix.pdf).
Dates are the first arXiv submission dates, not later journal dates.

## Observable ledger

The formula date is an upper bound on novelty: where only the dated February
work log survives, that date is used even if the formula may have been tried
earlier in an unrecorded session.

| Observable | Formula / point claim | Earliest repo evidence fixed | Experimental value used; first public date | Class | Fitting-risk annotation |
|---|---|---|---|---|---|
| `alpha_s(M_Z)` | `1/(2 phi^3)=0.118034` | `experiments_log.md`, EXP-043/044, 2026-02-04 | PDG values long predate 2026; paper uses `0.1179+-0.0009` | RETRODICTION | Target known; golden expression selected in an exploratory log. |
| `sin^2 theta_W(M_Z)` | `6/26=3/13` | `experiments_log.md`, EXP-045, 2026-02-04 | precision electroweak value public before 2026 | RETRODICTION | Scale and rational form were chosen with target known. |
| `alpha^-1` | smaller root of `2 pi alpha^2-20 phi^4 alpha+1=0` | February 2026 paper/log chain; exact first day unclear | CODATA `137.035999084(21)`, public before 2026 | RETRODICTION | Very high precision agreement cannot be treated as forecast; coefficient/normalisation search occurred after CODATA. |
| charged masses `mu,tau,u,c,t,d,s,b` | `m_e phi^n` plus norm-log/radiative corrections | `experiments_log.md`, 2026-02-04 onward | PDG 2024 masses, published 2024-08-01 (and much earlier measurements) | RETRODICTION | Explicit correction construction and global chi-square used known masses; eight targets and multiple discrete/correction choices create a large trials factor. |
| Koide `Q` | corrected masses give `Q=0.66668` | supplementary February 2026 history; exact day unclear | Koide relation published 1981 | RETRODICTION | Known target was used to assess corrections. |
| CKM `sin theta_12` / `V_us` | golden K-matrix exponent 3 plus correction | February 2026 manuscript/log history; exact day unclear | PDG 2024 `|V_us|=0.22431+-0.00085`, 2024-08-01 | RETRODICTION | Exponents/correction choices selected with CKM hierarchy known. |
| CKM `sin theta_23` / `V_cb` | exponent 7 plus correction | same | PDG 2024 `|V_cb|=(41.1+-1.2)10^-3`, 2024-08-01 | RETRODICTION | Same look-elsewhere family. |
| CKM `sin theta_13` / `V_ub` | exponent 12 plus correction | same | PDG 2024 CKM review, 2024-08-01 | RETRODICTION | Same look-elsewhere family. |
| CKM phase | `atan(sqrt(5))` | February 2026 manuscript history; exact day unclear | CKM global fits public before 2026 | RETRODICTION | Golden angle chosen after the phase was known. |
| Jarlskog `J` | CKM formula, `3.12 10^-5` | same | global CKM fits public before 2026 | RETRODICTION | Correlated consequence of already-retrodictive CKM inputs, not an independent test. |
| `m_Z` | `m_e phi^25` times measured running ratio | EXP-065--068 era, 2026-02-05 or earlier | LEP precision mass public decades earlier | RETRODICTION | Uses an experimental running input and a known target. |
| `m_H/m_W` | `sqrt(phi^2-16 alpha phi)` | `experiments_log.md`, EXP-070, 2026-02-05 | Higgs and W masses public before 2026 | RETRODICTION | Log calls it a discovery while comparing the known Higgs mass. FCC-ee can discriminate two post-data formula variants, but that future discrimination does not restore blind origin. |
| `m_e/m_P` | `alpha^(4 phi^2)` | February 2026 manuscript; exact day unclear | both constants public before 2026 | RETRODICTION | Dimensionless hierarchy target known. |
| `sin^2 theta_12` | `2/(phi+5)=0.302205` | present in `prompt_neutrino_masses.md`, 2026-02-23 | NuFIT 6.0 (2024-10-07); JUNO `0.3092+-0.0087` (2025-11-18) | RETRODICTION | Formula fixed after both targets; JUNO comparison is consistency only. |
| `sin^2 theta_13` | `1/(5*9)=1/45` | present in `prompt_neutrino_masses.md`, 2026-02-23 | NuFIT 6.0 `0.02215(+0.00056,-0.00058)`, 2024-10-07 | RETRODICTION | Paper admits `1/45` is essentially optimal and numerical minimum is `0.0224`; direct data-driven selection risk. |
| `sin^2 theta_23` | `4/7` | present in `prompt_neutrino_masses.md`, 2026-02-23 | global-fit value public before 2026 | RETRODICTION | Target and octant information known at construction. |
| `delta_CP` (PMNS) | `3 atan(sqrt(5))=197.715 deg` | present in `prompt_neutrino_masses.md`, 2026-02-23 | NuFIT 6.0 best fit `212(+26,-41) deg`, 2024-10-07 | RETRODICTION | Multiplier and golden angle selected with a broad known target. |
| `Delta m^2_21` | bare `m3^2 alpha phi^3` | exponent-35 package fixed no earlier than 2026-02-23 | JUNO `(7.50+-0.12)10^-5 eV^2`, 2025-11-18 | RETRODICTION | Prompt explicitly supplies measured splitting before derivation. |
| `Delta m^2_31` / `Delta m^2_32` | `m3^2` / difference; optional `phi^(1/45)` correction | no earlier than 2026-02-23; scoped correction documented by 2026-07-22 | NuFIT 6.0 values, 2024-10-07 | RETRODICTION | Correction and Variant-I scope were assessed against known splittings. Variant I is correctly PATTERN and data-driven. |
| ratio `Delta m^2_21/Delta m^2_31` | `alpha phi^3` | no earlier than 2026-02-23 | both splittings measured before 2026 | RETRODICTION | Selecting a simple golden/coupling expression with the ratio known carries a substantial trials factor. |
| `m2` and `m3` point values | `m3=2m_e/phi^35`, `m2=m3 sqrt(alpha phi^3)` | no earlier than 2026-02-23 | not directly measured, but constructed to reproduce already-measured splittings | RETRODICTION | They are algebraic reparameterisations of known oscillation splittings, so they are not blind absolute-mass forecasts. |
| `m1` | exactly zero (rank two) | no earlier than 2026-02-23 | no measurement; only absolute-mass constraints | BLIND | Clean point claim: any established `m1>0` falsifies it. |
| strict ordering | normal, with `m1=0` | no earlier than 2026-02-23 | ordering not settled at fixation or as of ledger date | BLIND | JUNO/combined oscillation programme provides the direct test; inverted ordering falsifies. |
| `sum m_nu` | `58.2--58.9 meV` by variant | no earlier than 2026-02-23 | DESI `<64.2 meV` (95%, LambdaCDM), 2025-03-19; no model-independent measurement | BLIND | Bound is not a measurement. A robust same-model upper bound below `58.24 meV`, or incompatible positive determination, falsifies. |
| `m_beta` | `8.77--8.87 meV` | no earlier than 2026-02-23 | KATRIN `<0.45 eV` (90%), arXiv 2024-06-19 | BLIND | Far below current reach; direct result excluding this interval falsifies. |
| `m_betabeta` | phase-specific `3.10--3.5 meV` | no earlier than 2026-02-23 | no measurement, only experiment/NME-dependent limits | BLIND | Positive measurement incompatible with the stated phase-specific interval falsifies; NME/phase assumptions must be retained. |
| Majorana phases | `alpha_1,2=+-4pi/5` | February 2026 manuscript; exact first day unclear | unmeasured | BLIND | Test only through phase-sensitive observables; no current direct value. |
| no fourth generation | `N_gen=3` / fourth forbidden | February 2026 framework | collider constraints existed before construction, but exact universal absence is not directly measurable as one scalar | AMBIGUOUS | Both old exclusions and a broader structural claim are involved; no clean publication/formula ordering certificate. |
| `a_mu(NP)=0` | no new-physics contribution | February 2026 manuscript | anomaly interpretation was already public and theory-dependent | AMBIGUOUS | A null BSM component is not an isolated observable without an SM calculation; provenance and future falsifier are model-dependent. |
| `theta_QCD=0`, no axion | exact Galois invariance | February 2026 manuscript | neutron-EDM bounds long predate construction; exact value unmeasured | BLIND | A nonzero strong-CP angle or confirmed QCD axion falsifies the claim; present bounds are not a zero measurement. |
| inflation `r` | `12/60^2=0.003` | February 2026 manuscript | only bounds existed | BLIND | A reliable primordial result outside the stated window (`r>0.01` or `<0.001`) falsifies this registered point claim. |
| scalar tilt `n_s` | Starobinsky value `0.9654` | February 2026 manuscript | Planck 2018 value public | RETRODICTION | Standard Starobinsky consistency, selected after Planck. |
| scalar amplitude `A_s` | patterned `f4=N^4/(2pi)` gives about `2.1e-9` | February 2026 manuscript | Planck amplitude public | RETRODICTION | Explicit pattern with free test-function moment chosen against a known value. |
| cosmological constant | `Lambda_P=alpha^(57-alpha_s)` | February 2026 manuscript | observed order `10^-122` long known | RETRODICTION | Exponent decomposition searched with target known; paper already labels mechanism PATTERN. |
| `Omega_DM/Omega_b` | `7-phi` | February 2026 manuscript | cosmological ratio public before 2026 | RETRODICTION | Simple-expression match selected post data. |
| dark-sector particle spectrum | nine stable EM-dark partners, stated masses/degeneracy | February 2026 manuscript | no particles measured | BLIND | Collider/astrophysical detection inconsistent with the spectrum, stability or darkness falsifies the scoped claim; production abundance is not derived. |
| no Galois-sector stochastic GW background | null | February 2026 manuscript | no sector-specific measurement | BLIND | A background causally identified with that sector falsifies; attribution makes the window presently open-ended. |
| graviton mass | about `10^-34 eV` after setting `R=R_Hubble` | February 2026 manuscript | only much weaker bounds | BLIND | A nonzero measurement/bound excluding the point falsifies, conditional on the external radius identification. |

## Look-elsewhere estimate (documented, not exculpatory)

This is an illustrative trials count, **not** a reconstruction of every
expression inspected in the sessions.  Define a deliberately narrow family

`F = { (p/q) phi^k, (p+q phi)/r : 1<=p,q,r<=10, -40<=k<=40 }`

plus reciprocals where finite.  It already contains 8,100 power/rational
forms and 1,000 linear golden forms before duplicate removal.  For each
retrodiction, `verify_prediction_provenance.py` reports how many distinct
members land inside the pinned one-sigma interval after using the natural
dimensionless normalization stated in the script.  This family omits
products with `alpha`, arctangents, powers of `a1,b1,N_eig`, correction
terms, scale choices, and correlated variant choices, so it is a lower-bound
illustration of flexibility, not a p-value.  Conversely, its members were
not demonstrably all tried, so the hit count is not a literal effective
trials factor.

The honest conclusion is:

> Consistency was achieved with values known at construction time.
> Evidential weight rests on the BLIND set and on internal derivational
> rigidity, not on retrodictive agreement.

## The genuinely blind asset set and falsification windows

The registered blind neutrino assets, fixed no earlier than 2026-02-23, are
`m1=0`, strict normal ordering, `sum m_nu=58.2--58.9 meV`, `m_beta=8.77--8.87
meV`, and phase-specific `m_betabeta=3.10--3.5 meV`.  Their windows are,
respectively: a future absolute-mass reconstruction; JUNO plus global
ordering data; DESI/Euclid/CMB analyses with the cosmological model stated;
future sub-10-meV beta spectroscopy; and next-generation neutrinoless
double-beta searches with nuclear-matrix and phase assumptions stated.

Other blind, more model-dependent registrations are the Majorana phases,
exact strong-CP null/no axion, `r=0.003`, the dark-partner spectrum and
degeneracy, the sector-specific GW null, and the conditional graviton-mass
point.  These are the paper's only genuine predictive assets in this ledger.
No measured oscillation parameter belongs to that set.

