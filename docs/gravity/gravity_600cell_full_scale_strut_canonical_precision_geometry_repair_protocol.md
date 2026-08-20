# Preregistered precision-aware geometry-control repair

Date: 2026-08-20  
Status: **preregistered control repair; primary classification remains frozen**

## Motivation and provenance caveat

The diagnostic frozen in commit `2dc7d91` localized the only P160 geometry
failure to the absolute test `maximum_imaginary < 1e-130`. The observed residues
are of order `10^-119`, consistent with working precision `10^-160` amplified by
the smallest central-difference step `10^-40`.

This repair is necessarily designed after observing that control failure. It is
therefore not target-blind evidence for the primary intersection result. To avoid
merely loosening a threshold, it adds a new, higher-precision convergence
calculation and leaves every primary P100/P160 matrix and classification gate
unchanged.

## Frozen object and hypotheses

The carrier, action kernel, parities, sectors, source/target convention, P100 and
P160 matrices, interval determinants, candidate vector, no-refit gate `1e-50`,
and outcome hierarchy remain unchanged.

The exact mathematical action kernel is expected to be real. The nonzero
imaginary part produced by its complex Lorentzian intermediate formulae is
classified as arithmetic contamination only if an independently recomputed
higher-precision value converges toward zero while its real entries remain
stable.

## New P200G geometry-only audit

For each parity, rebuild the geometry and all four representative action kernels
at 200 decimal digits. Use exactly the P160 finite-difference steps:

- operational primary `1e-40`;
- operational shadow `1e-30`;
- validation primary `3e-40`;
- validation shadow `3e-30`.

Do not project sectors, build intersection matrices, refit candidates, or repeat
the scientific classification at P200G.

The P160 imaginary-residue test is replaced by the conjunction below:

1. all six non-imaginary geometry controls pass at P200G;
2. the key support of every P200G kernel equals its P160 counterpart exactly;
3. `max_imaginary(P200G) < 1e-30 * max_imaginary(P160)`;
4. `max_imaginary(P200G) < 1e-150`;
5. the maximum absolute P160/P200G difference over every real kernel entry is
   below `1e-110`.

The predicted arithmetic scale is `10^-200 / 10^-40 = 10^-160`; the decay and
absolute gates leave ten decimal orders of slack while remaining at least sixty
orders below the frozen homogeneous no-refit gate.

P100 retains its original absolute `1e-70` control. At P160, the six
non-imaginary controls must pass locally and the new P200G conjunction must pass
before `all_level_controls` can be true.

## Preregistered outcomes

- If either parity fails any P200G condition, the outcome remains
  `FULL_SCALE_STRUT_CANONICAL_PRECISION_CONTROL_FAILED`.
- If both pass, the unchanged preregistered hierarchy assigns the scientific
  outcome from the already frozen P100/P160 calculations.
- Any accepted material claim still requires the separate adversarial replication
  mandated by Rule 4.
