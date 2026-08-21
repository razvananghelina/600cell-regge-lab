# Accepted constrained H4 boundary response

Date: 2026-08-21

Status: **DERIVED COMPUTATIONAL, adversarially corroborated** in the frozen
homogeneous `H4` sector.

## Exact object and hypotheses

On the curvature-matched internally stationary Lorentzian product over

```text
K0=P(sd K_600), f=(2640,17040,28800,14400), tau0=0.0102,
```

use the corrected Regge action and boundary terms, the selected rank masses,
the accepted product-duration null line `n`, and its nonzero compatibility
covector `c=H_bi n`.  In the fixed total-orbit log-squared-edge coordinates,
the object is the constrained linearized boundary-momentum bilinear form on

```text
ker(c^T) subset R^12,
```

after stationary elimination on the nine-dimensional internal slice
`ker(n^T)`.  It is not an unconstrained Schur complement and is not assumed to
be the intrinsic Hessian of an unconstructed nonlinear boundary surface.

The tested question is target-free: after the fixed temporal-orientation
identification, how many response classes occur among all 24 colour-ordered
staircase triangulations?

## Provenance ledger

| stage | commit |
|---|---|
| prior-art gate, final clarification | `8ecbd2a` |
| primary protocol, final version | `be10390` |
| primary verifier registration | `69ace62` |
| primary failure/correction chain | `cddb3ca`--`7f2bfa9` |
| corrected primary implementation/result | `e281c57`, `d9f3abc` |
| direct adversarial protocol, final dependency freeze | `8c16996`--`a9ee74b` |
| direct adversarial verifier registration | `7de8f9e` |
| first direct `15/17` failure frozen | `136de6d` |
| auxiliary diagnostic protocol, final precision floor | `af8aa91`, `9e44dd7` |
| diagnostic registration/result | `c991173`, `516a721` |
| corrected adjudication protocol | `2a10336` |
| unchanged full rerun receipt | `ab81207` |
| corrected adjudicator registration | `defb04c` |

No full suite or deferred nonlinear root census was run.

## Primary route

The primary route differentiated the analytic 22-component action gradient,
formed the full Hessian, restricted to algebraic bases

```text
P=E(c,3), Q=E(n,9),
```

and solved the nonsingular `9 x 9` internal block.  Its corrected verifier
passed `19/19` twice with byte-identical artifact

```text
reproducible/gravity_600cell_refined_h4_constrained_response_corrected.json
SHA-256 85adea23f6a19153f61f3ed066137a5e40ab77b8901d4cc81cfc4f864e0bc093.
```

It found one class across all 24 schedules, time-reversal covariance, internal
minimum eigenvalue `1.3780099e-5`, and basis invariance.  At this stage the
claim remained primary-only and therefore **OPEN**.

## Mechanically independent direct-action route

The adversarial route never read a primary Hessian block or internal lift.
For every schedule it evaluated the complete scalar action along all

```text
20 basis directions + 190 pair sums = 210 directions,
```

used five centred step sizes and eighth-order extrapolation, reconstructed the
entire `20 x 20` restricted second variation by polarization, then independently
solved its `9 x 9` internal block.  The complete audit comprised 55,610
scalar-action evaluations and stored all 24 second variations, lifts and
responses.

The direct response envelopes lie between `1.3984066e-54` and
`1.3984072e-54`.  The two complete 220-digit repeats agree with the 180-digit
route with maximum used fraction `6.8530e-91`.  A deliberate one-entry
corruption is resolved by `1.9085e24` primary-comparison gates.

The unchanged full verifier was executed twice and regenerated the same
artifact byte-for-byte:

```text
reproducible/gravity_600cell_refined_h4_constrained_response_adversarial.json
SHA-256 a23ef4cc23d08ad8768f1df66789aa900cdb95a7f3529486df80697a53b1fe81.
```

Its historical outcome remains honestly `15/17 CONTROL_FAILED`; the two
auxiliary failures were not erased.

## Resolution of the two failed controls

The preregistered diagnostic passed `13/13` twice and produced byte-identical

```text
reproducible/gravity_600cell_refined_h4_constrained_response_auxiliary_diagnostic.json
SHA-256 f66177326afc3b3457a60b544745b739cbaa6b6d6e7f367b57d60f31eeeddeb7.
```

### Stationarity

Across all ten internal directions, all 24 schedules and both 180/220 decimal
precisions,

```text
legacy fourth-order maximum  5.1661660e-54,
tenth-order maximum          5.9097062e-77.
```

All zero and cross-precision gates pass.  A known nonstationary linear term
`1e-20` is recovered and rejected, so the test cannot make a nonzero gradient
look stationary.

**DERIVED COMPUTATIONAL:** the historical stationarity failure was
finite-difference truncation, not a resolved physical residual.

### Lorentzian off-shell terms

For all `24*210=5040` schedule-direction pairs, the individual imaginary
curvature maxima at successive half-steps are

```text
2.0750292e-8, 1.0375146e-8, 5.1875730e-9,
ratios 2.0000000, 2.0000000.
```

All `5040/5040` finest-step pairs have smooth odd-leading behaviour, with
maximum even/odd ratio `1.6772e-10`.  At the same time the relative imaginary
complete action is `9.1262e-177`, the angle identity residual is
`8.6673e-178`, and the minimum branch argument is `0.91222028`.  A synthetic
discontinuous branch gives halving ratios `(1,1)` and is rejected.

**DERIVED COMPUTATIONAL / STRUCTURAL:** termwise off-shell Lorentzian
curvatures may be complex while vanishing smoothly at the real background and
combining into a real action.  Termwise reality was the wrong branch gate;
analytic continuity and reality of the complete action are the relevant gates
for this calculation.

## Corrected adjudication and accepted result

The corrected adjudicator did not trust stored class labels.  It reconstructed
the temporal reversal action from the accepted compatibility row, reloaded all
24 direct matrices, rebuilt the target-free class census, and only then loaded
the primary matrices.  It passed `12/12` twice with byte-identical artifact

```text
reproducible/gravity_600cell_refined_h4_constrained_response_corrected_adjudication.json
SHA-256 6b4194c278548d06a9a59a74867c55c80e97ac5c6fcb86b9ee65b4db0905ff5f.
```

The frozen outcome is

```text
CORRECTED_ADJUDICATION_REFINED_H4_CONSTRAINED_RESPONSE_CORROBORATED.
```

The recomputed results are:

- direct class count `1`, containing all `24/24` schedules;
- time reversal covariant, maximum difference `1.03755e-59`, using at most
  `2.23905e-6` of its direct envelope;
- direct/primary matches `24/24`, maximum difference `1.16741e-28`, using at
  most `6.66667e-5` of the combined cross-method gate;
- internal minimum eigenvalue `1.37800993565e-5`;
- maximum normalized internal solve residual `1.79540e-181`.

**DERIVED COMPUTATIONAL, ADVERSARIALLY CORROBORATED:** on this fixed on-shell
seed and in the complete homogeneous `H4` invariant sector, the constrained
linearized boundary-momentum response is independent, within the frozen
direct numerical envelopes, of all 24 bare staircase schedules (after the
fixed temporal-orientation identification).

This is stronger than a dimension match or equality inside the much wider
primary Hessian envelope: the independent scalar-action route resolves the
class directly at approximately `1e-54`.

## Post-result primary-literature check

The learned concepts are not new in general:

- Dittrich and Hoehn use Hamilton's principal function to generate canonical
  simplicial evolution and implement Pachner moves in
  [Canonical simplicial gravity, arXiv:1108.1974](https://arxiv.org/abs/1108.1974).
- Dittrich and Steinhaus analyze triangulation independence and factorized
  linearized Regge actions under Pachner moves in
  [Path integral measure and triangulation independence in discrete gravity,
  arXiv:1110.6866](https://arxiv.org/abs/1110.6866).
- Bahr and Dittrich explain improved/perfect actions as the route to exact
  discretization independence in
  [Improved and Perfect Actions in Discrete Gravity,
  arXiv:0907.4323](https://arxiv.org/abs/0907.4323).
- Bahr and Dittrich show that curvature generically breaks exact discrete
  gauge symmetry and yields pseudo-constraints in
  [(Broken) Gauge Symmetries and Constraints in Regge Calculus,
  arXiv:0905.1670](https://arxiv.org/abs/0905.1670).
- Asante, Dittrich and Padua-Arguelles study analytic continuation and branch
  cuts of the complex Lorentzian Regge action in
  [Complex actions and causality violations,
  arXiv:2112.15387](https://arxiv.org/abs/2112.15387).
- De Felice and Fabri explicitly evolve a 600-cell with a Sorkin algorithm in
  [Singularities of the closed RW metric in Regge Calculus,
  arXiv:gr-qc/0106077](https://arxiv.org/abs/gr-qc/0106077).

These sources make the general interpretation unsurprising and sharply warn
against promoting the result to a gauge-symmetry theorem on a curved lattice.
The search did not identify the same projected barycentric carrier, selected
rank-matter seed, constrained null reduction and complete 24-schedule `H4`
census.  Search absence is not proof: external novelty remains **OPEN** pending
a dedicated literature review.

## Limits and next gate

- **DERIVED NEGATIVE:** this does not restore an unconstrained Schur
  complement or license a Moore--Penrose extension off `ker(c^T)`.
- **STRUCTURAL:** schedule independence is established only in the homogeneous
  invariant sector; high symmetry may be exactly why the schedules collapse.
- **OPEN:** integration of the compatibility hyperplane to a nonlinear
  admissible boundary family.
- **NOT ESTABLISHED:** a physical tick, nonhomogeneous propagation, a wave
  equation, effective `c`, Newton's `G`, Planck units or particle masses.

The next load-bearing gate is therefore not another homogeneous rerun.  It is
the target-free nonhomogeneous quadratic operator: decompose boundary and
internal perturbations into symmetry sectors, perform the same constrained
stationary elimination in every licensed sector, and test schedule covariance,
stability and the presence of propagating mode pairs before comparing with any
continuum dispersion law.

