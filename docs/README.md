# Documentation index

The public documentation is now restricted to the active 600-cell Regge
programme. Legacy fitted particle-physics notes are absent from the current
tree and remain recoverable only through Git history.

## Authoritative status

Read [`gravity/CURRENT_STATUS.md`](gravity/CURRENT_STATUS.md) before starting
or interpreting another calculation. It records the accepted checkpoint,
complete hypotheses, exact artifact hashes, negative results and next open
gate.

## Current finite-height chain

Read in this order:

1. [finite-height prior-art gate](gravity/gravity_600cell_finite_height_prior_art.md);
2. [exact finite-height classification](gravity/gravity_600cell_finite_height_classification_result.md);
3. [two-slab composition obstruction](gravity/gravity_600cell_finite_height_composition_result.md);
4. [canonical selector audit](gravity/gravity_600cell_finite_height_selector_result.md);
5. [third-slab branch distinction](gravity/gravity_600cell_finite_height_third_slab_result.md);
6. [fourth-slab continuation](gravity/gravity_600cell_finite_height_fourth_slab_result.md);
7. [scale-free map and fifth-slab continuation](gravity/gravity_600cell_finite_height_asymptotic_map_result.md);
8. [invariant half-strip and finite-step forward continuation](gravity/gravity_600cell_finite_height_invariant_region_result.md).
9. [incoming-state basin prior-art gate](gravity/gravity_600cell_finite_height_incoming_basin_prior_art.md);
10. [frozen incoming-state discovery protocol](gravity/gravity_600cell_finite_height_incoming_basin_discovery_protocol.md);
11. [incoming-state candidate skeleton](gravity/gravity_600cell_finite_height_incoming_basin_discovery_result.md).
12. [local-stability prior-art gate](gravity/gravity_600cell_finite_height_local_signature_prior_art.md);
13. [primary local-stability protocol](gravity/gravity_600cell_finite_height_local_signature_protocol.md);
14. [adversarially corroborated local branch theorem](gravity/gravity_600cell_finite_height_local_signature_result.md).
15. [finite-height quadratic carrier canonicity](gravity/gravity_600cell_finite_height_carrier_quadratic_result.md).
16. [finite-height internal-carrier rank](gravity/gravity_600cell_finite_height_internal_carrier_rank_result.md).

The chain proves the exact continuous limiting family and, by induction in a
rigorous invariant region, a unique physical successor at every later finite
step for the representative branch.  It does **not** prove convergence,
infinite total proper duration or completeness.  The complete-domain finite
census finds 36 distinct depth-four signatures and leaves 1080 diagnostic
inputs with a branch live outside that region.  Generic incoming-state
evolution therefore remains **OPEN**; local physics and a fundamental tick
remain not derived.

The local theorem proves that the representative `DEAD+ENTERED_D` tree is
stable on some unspecified nonzero neighbourhood of `v=3/2`. It removes an
isolated-point objection but does not select `v=3/2`, measure a basin or make
future extendibility a local law.

The latest nonhomogeneous gate is mixed.  Two mechanically different
calculations find no staircase-parity dependence in the one-sided quadratic
action.  A second preregistered pair then finds that internal stationarity
kills every nonhomogeneous direction of the exact rank-240 carrier and leaves
one common homogeneous line.  This is **DERIVED COMPUTATIONAL** for the
internal constraint only.  The fixed incoming canonical momentum has not yet
been imposed, so physical boundary evolution remains **OPEN**.

## Continuum-control chain

- [exact homogeneous all-tick result](gravity/gravity_600cell_cellular_weak_lapse_all_n_result.md);
- [blind projected-refinement coefficients](gravity/gravity_600cell_projected_refinement_acceleration_blind_result.md);
- [post-commit Friedmann comparison](gravity/gravity_600cell_projected_refinement_acceleration_comparison_result.md).

## Executable evidence

Registered verifiers and frozen artifacts are in
[`../reproducible/`](../reproducible/). Shared exact 600-cell construction is
in [`../commons/`](../commons/).
