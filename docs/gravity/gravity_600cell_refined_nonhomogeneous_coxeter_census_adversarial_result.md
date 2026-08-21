# Adversarial result: all 24 internal kernels are the duration line

Date: 2026-08-21

Status: **41/41 adversarial checks passed.**  The complete internal-kernel
claim is computationally corroborated for all 24 labelled schedules, subject
to the shared frozen Hessian definition and the declared numerical error
model.

## Frozen provenance

- accepted schedule-0 adversarial result: `66e47a7`;
- frozen primary all-schedule result: `d0adbac`;
- adversarial all-schedule protocol: `3e8e797`;
- registered adversarial implementation before first execution: `c3d10a4`;
- adversarial artifact:

  ```text
  reproducible/gravity_600cell_refined_nonhomogeneous_coxeter_census_adversarial.json
  SHA-256 787d09d75f810a88d0ce09f33a5dc40a9c9b924ecfbba8e0319feb3395bf8e9b
  ```

Only the active verifier was run.  Neither the full suite nor the old failed
sparse census was rerun.

## Claim and independence boundary

For each of the twelve time-reversal representatives, with the complete
hypothesis list frozen in `2eef0e1`, define

```text
K_s = [[C_s, n_s],
       [n_s^T, 0]].
```

Schedule 0 was already independently checked.  This run attacked the other
eleven representatives with ambient `R4` reflections, geometrically induced
edge permutations, explicit sparse isometries `Q_k`, sparse products
`Q_k^* C_s Q_k`, SciPy's `evd` driver and eigenvector residuals lifted to the
original sparse matrix.  No primary permutation, block or diagonalization
function was loaded.

The primary spectra were not parsed until all eleven adversarial spectra had
been constructed.  Atomic blind-phase checkpoints contained no primary
comparison.

This is a mechanically independent spectral replication, not an independent
derivation of `C_s`.  Both routes share the frozen Regge action, mass rule,
coordinates and binary64 source-matrix assembly.  A common upstream physics
or assembly error is therefore outside the independence claim.

## Controls

The geometric route independently reproduces:

```text
H4 pair orders              (3,3,5)
Coxeter and reverse orders  30
vertex matching residual    1.12038e-10  < 5e-8
edge cycles per schedule    656 x 30
```

For every new representative:

- the source matrix reproduces its preregistered CSR digest;
- forward and reverse edge actions are inverse;
- all thirty `Q_k` matrices are isometries, with maximum defect
  `2.22045e-16`;
- all sixteen independent blocks and the weighted dimension `19,681` are
  complete;
- lifted residual, trace and Frobenius/Parseval controls pass;
- a `1e-4` diagonal corruption violates covariance;
- no zero-compatible candidate occurs.

The single-reflection control has order 2 rather than 30.  The wrong-phase
control on the first new representative fails the unrelabeled same-sector
test as preregistered.

## Adversarial margins and primary comparison

| schedule | weakest sector | eigenvalue | gate | margin | primary difference |
|---:|---:|---:|---:|---:|---:|
| 1 | 1 | `1.455649e-9` | `2.4569e-10` | 5.9247 | `7.105e-15` |
| 2 | 1 | `1.455649e-9` | `2.5654e-10` | 5.6741 | `6.661e-15` |
| 3 | 1 | `1.455649e-9` | `2.8921e-10` | **5.0332** | `5.329e-15` |
| 4 | 1 | `1.455649e-9` | `2.7847e-10` | 5.2273 | `5.329e-15` |
| 5 | 1 | `1.455649e-9` | `2.8913e-10` | 5.0346 | `7.105e-15` |
| 6 | 1 | `1.455649e-9` | `2.4581e-10` | 5.9218 | `1.021e-14` |
| 7 | 11 | `1.455649e-9` | `2.4575e-10` | 5.9234 | `1.510e-14` |
| 8 | 1 | `1.455649e-9` | `2.5888e-10` | 5.6228 | `1.155e-14` |
| 10 | 11 | `1.455649e-9` | `2.7870e-10` | 5.2230 | `7.105e-15` |
| 12 | 1 | `1.455649e-9` | `2.4847e-10` | 5.8585 | `1.110e-14` |
| 14 | 1 | `1.455649e-9` | `2.4834e-10` | 5.8616 | `5.329e-15` |

The least adversarial zero-exclusion margin is `5.033170`, at schedule 3,
sector 1.  The worst complete-spectrum disagreement is `1.50990e-14`, at
schedule 7, under a `4.91262e-10` gate.  These are resolved margins, not
near-threshold classifications.

## Verdict

- **DERIVED COMPUTATIONAL / ADVERSARIALLY CORROBORATED:** for every one of
  the twelve representatives, `K_s` is numerically nonsingular and hence
  `ker(C_s) = span(n_s)` under the complete frozen hypotheses.
- **DERIVED COMPUTATIONAL NEGATIVE / ADVERSARIALLY CORROBORATED:** using the
  explicitly checked time-reversal congruences, none of the 24 labelled
  schedules has an additional nonhomogeneous internal zero mode.
- **STRUCTURAL:** the duration tangent is the unique surviving internal
  degeneracy; all transverse internal directions are stiff at the fixed
  configuration.
- **OPEN:** a symbolic exact certificate and an independent derivation of the
  source Hessian from the physical action.
- **OPEN / NOT TESTED:** boundary evolution, propagation, a graviton law,
  `c`, `G`, or Planck units.

The computational internal-kernel question is now closed for the frozen
matrices.  This is not yet a derived tick.  The value `tau0=0.0102` is still
an input, and the surviving direction is exactly the freedom that changes
the product duration.

## What is required to derive a tick

A boundary-to-boundary variational problem must treat the temporal lengths
or duration as variables:

```text
S(q_lower, q_upper, tau; matter).
```

With physical data fixed on the lower boundary, the internal equations and
`partial S / partial tau = 0` must select an isolated nonzero solution
`(q_upper, tau_*)`.  A candidate tick is accepted only if it is not the
zero-lapse copy, is not removable by gauge, is stable under perturbations and
schedule changes, and survives refinement.

Conserved dust proper time is the cleanest existing clock candidate.  If
mass conservation and the dust clock still leave `tau` flat, the honest
verdict is that duration is gauge and this geometry has no intrinsic tick.
If only the dimensionless ratio `tau_*/L` is selected, conversion to seconds
still requires a physical scale; a dimensionful Planck tick cannot arise
from scale-free combinatorics alone.
