# Canonical Legendre rank of the 600-cell dust slab

Date: 2026-08-16

Status: **DERIVED COMPUTATIONAL LOCAL — canonical rank 65/65 in both
order-24 schedule carriers.**  This accepts a separate nonlinear canonical
inversion protocol.  It does not yet construct a next spatial frame.

## 1. Provenance ledger

The numerical spectrum was not inspected before the following commits:

| Item | Commit |
|---|---:|
| prior-art gate | `31b1690` |
| main preregistration | `4837a16` |
| derivative-tolerance clarification | `be145e4` |
| pseudoconstraint/verifier-outcome clarification | `92182e1` |
| calibrated null-projection rule | `9684918` |
| registered implementation before evaluation | `9b34e48` |

The first complete numerical evaluation produced the same rank result below
but then failed while serializing a `numpy.bool_`.  Commit `7459b83` casts
only certificate booleans and two displayed errors to native JSON types.  It
does not change an action, derivative, matrix, tolerance, singular value or
scientific outcome.  The verifier was then rerun from the beginning.

Certified artifact:

```text
reproducible/gravity_600cell_dust_canonical_legendre_rank.json
SHA256 a3e820a2d9b96f5f6388f34522e9f55264cfb4cfbd4ecb022de844216a32b25c
```

The verifier is registered exactly once in `reproducible/run_all.py`.  The
registry contains 220 entries and 220 distinct names.  In accordance with the
active instruction, only this verifier was run, not the full suite.

## 2. Complete hypotheses and object

The claim below is conditional on all of the following:

1. the fixed published De Felice--Fabri dust slab at `tau=0.0102`;
2. one of the two derived even/odd order-24 schedule carriers;
3. the complete Lorentzian Regge curvature plus dust action and its certified
   complex-angle branch;
4. the 95 logarithmic coordinates
   `(q_old[30], x_internal[35], q_new[30])`;
5. the 65 canonical unknowns `(x_internal[35], q_new[30])` at fixed old
   geometry and fixed pre-momentum;
6. the pre-Legendre Jacobian

```text
J_can = [[ K_XX,  K_XN],
         [-K_OX, -K_ON]],
```

   where `K` is the complete action Hessian divided by the common orbit size
   24;
7. the fixed 100-decimal central-difference calibration and frozen global
   singular-value error bands in the preregistration.

This is not a statement about the unreduced 720-edge boundary phase space,
about arbitrary triangulations, or about a continuum limit.

## 3. Controls

The targeted verifier returns **17/17**.

For both parities:

- all 761 base/perturbed gradients remain on the one-timelike-direction
  Lorentzian branch;
- all 9,025 Hessian entries pass the operational/validation calibration;
- the complete 95-dimensional Hessian is reciprocal;
- the independently rebuilt internal and final-boundary blocks reproduce the
  stored precision-corrected blocks below the frozen `1e-6` relative
  Frobenius gate.

The last control is not hidden: its discrepancies are approximately
`7.083e-7` for `K_XX` and `6.833e-7` for `K_XN`, hence pass but are close to
the frozen gate.  The new rank is nevertheless separated from the calibrated
arbitrary-precision error by thirteen orders beyond that comparison:

| Quantity | even | odd |
|---|---:|---:|
| `s_max` | 2961.1840490 | 2963.5578384 |
| `s_min` | 4.2445618120e-9 | 4.2445618120e-9 |
| `epsilon_global` | 6.80219e-23 | 6.80242e-23 |
| `s_min/epsilon_global` | 6.23999e13 | 6.23979e13 |
| condition number | 6.97642e11 | 6.98201e11 |
| resolved rank | 65 | 65 |
| error-consistent nullity | 0 | 0 |
| numerically open modes | 0 | 0 |

The binary64 spectrum is only a control; its normalized discrepancy from the
100-decimal spectrum is below `4.55e-14` and is not used for the verdict.

## 4. Result and interpretation

**DERIVED COMPUTATIONAL LOCAL.**  In both parity carriers, `J_can` has
calibrated numerical rank 65.  Fixing the old geometry and pre-momentum lifts
the exact one-parameter degeneracy previously seen in the Dirichlet internal
problem.  Consequently the canonical implicit system is locally invertible
at the published slab within this reduced carrier.

This is the missing permission for a high-precision Newton solve; it is not
the solve itself.  No new frame, expansion factor, physical tick duration,
inflation, value of `c`, Planck scale or particle mass is derived here.

Exactly five singular directions satisfy the preregistered diagnostic
`RESOLVED_PSEUDOCONSTRAINTS` criterion `s/s_max < 1e-6`.  They form an almost
degenerate cluster at `4.244561812e-9`, separated from the next singular value
near 43.08 by a factor about `1.015e10`:

- one is `0.999999999982` on the analytic collective internal-lapse sector;
- four are almost entirely internal-transverse;
- their final-boundary components are tiny, but none is a calibrated null.

**STRUCTURAL.**  Calling these five weak directions pseudo-constraints is
consistent with their nonzero lifting on a curved, matter-filled discrete
background and with the literature, but this calculation has not yet shown
their scaling with curvature, matter density or refinement.  They must not be
called gauge modes.

The two parity matrices are not entrywise the same (relative Frobenius
difference about 0.630), and their strong spectra differ by up to about 2.55%
under the displayed normalization.  Their identical weak cluster is therefore
not merely a byte-for-byte duplicate calculation.

## 5. Post-result prior-art audit

**KNOWN.**  Bahr and Dittrich show that curved Regge solutions generally lose
exact discrete gauge symmetries and replace constraints by pseudo-constraints:
<https://arxiv.org/abs/0905.1670>.  Dittrich and Hoehn derive the corresponding
canonical formalism from the discrete action and explain the dependence of
pseudo-constraints on background data:
<https://arxiv.org/abs/0912.1817>.  Hoehn's flat-background analysis identifies
exact vertex-displacement generators and lattice gravitons only in the
linearized flat regime: <https://arxiv.org/abs/1411.5672>.

**KNOWN.**  Sorkin-type implicit evolution and 600-cell Friedmann models
predate this repository: <https://arxiv.org/abs/gr-qc/9411008>.  De Felice and
Fabri specifically evolve a dust 600-cell and discuss its causal endpoint:
<https://arxiv.org/abs/gr-qc/0009093> and
<https://arxiv.org/abs/gr-qc/0106077>.

**OPEN NOVELTY.**  The post-result search did not locate the present complete
`65 x 65` pre-Legendre rank census for the exact order-24 quotient, parity
pair, Lorentzian branch and dust action used here.  A search cannot prove
novelty; external novelty remains OPEN pending a dedicated review.

## 6. What this says about the nonzero seam

The earlier two-slab cusp lies, to numerical precision, entirely on the
homogeneous boundary-scale direction.  Together with the present full-rank
canonical Jacobian this makes a first local continuation mathematically
well-posed in the reduced carrier.  It still does not prove that the next
solution expands: symmetry forces the cusp direction at the published point,
and the ill-conditioned five-mode cluster can amplify small right-hand-side
components.

The statement “only one scalar Newton step is missing” is therefore too
strong.  A scalar solve is a useful discrete-Friedmann control, but the
acceptance test must solve all 65 variables and then measure whether the
result remains on the homogeneous-scale subspace.

## 7. Next frozen sequence

The next mission must be preregistered separately:

1. **reproduction control:** with published `q_old` and published `p_pre` as
   canonical data, a 65-variable high-precision solve must recover the
   published internal variables and `q_new`;
2. **candidate forward tick:** set the incoming canonical momentum to the
   correctly permuted published `p_post`, with no fitted coefficient, and
   solve the same complete 65 equations;
3. decompose the displacement into internal lapse/transverse and final
   scale/shape sectors before calling it expansion;
4. reject or mark OPEN any solve that leaves the certified Lorentzian branch,
   depends materially on precision/initialization, or is not unique in the
   calibrated local basin;
5. only after that, perturb the symmetric data and test linear stability.

The current result advances this sequence from **blocked by unknown canonical
rank** to **allowed to attempt the nonlinear reproduction control**.
