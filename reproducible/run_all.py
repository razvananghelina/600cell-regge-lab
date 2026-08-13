"""
Run all verification scripts for the 600-cell framework paper.
"""
import os
import subprocess
import sys
import time

# Always run from the directory containing this script
os.chdir(os.path.dirname(os.path.abspath(__file__)))

scripts = [
    "verify_coupling_constants.py",
    "verify_spectrum_600cell.py",
    "verify_masses_and_mixing.py",
    "verify_berry_phase.py",
    "verify_spectral_action.py",
    "verify_discrete_scalar_response.py",
    "verify_mckay_chirality.py",
    "verify_galois_kernel.py",
    "verify_neutrino_masses.py",
    "verify_juno_comparison.py",
    "verify_prediction_provenance.py",
    "verify_muon_gminus2_null.py",
    "verify_polytope_uniqueness.py",
    "verify_tqft_entanglement.py",
    "verify_modular_tqft_layer.py",
    "verify_higgs_mass.py",
    "verify_ckm_structural.py",
    "verify_derived_quantities.py",
    "verify_e8_decomposition_ckm.py",
    "verify_uniqueness_quantized.py",
    "verify_uniqueness_simple_edges.py",
    "verify_minimal_edge_lifts.py",
    "verify_branch_identification.py",
    "verify_global_uniqueness_constructive.py",
    "verify_vacuum_selection.py",
    "verify_universality_class.py",
    "verify_edge_gauge_spectrum.py",
    "verify_gauge_continuum_map.py",
    "verify_low_mode_sampling_intertwiner.py",
    "verify_vector_sampling_intertwiner.py",
    "verify_a5_equivariant_brackets.py",
    "verify_matter_trace_indices.py",
    "verify_matter_module.py",
    "verify_matter_functor.py",
    "verify_galois_doubling.py",
    "verify_bimodule_krajewski.py",
    "verify_edge_matter.py",
    "verify_preprojective_matter.py",
    "verify_bratteli_tower.py",
    "verify_segregation_theorem.py",
    "verify_dirac_selection.py",
    "verify_c3_dynamical_selection.py",
    "verify_variational_bootstrap.py",
    "verify_neutral_vacuum_scale.py",
    "verify_kahler_dirac.py",
    "verify_kahler_dirac_local_tick.py",
    "verify_kahler_dirac_tick_vertex_isotropy.py",
    "verify_kahler_dirac_tick_refinement.py",
    "verify_weighted_szegedy_metric_nogo.py",
    "verify_incidence_polynomial_metric_nogo.py",
    "verify_whitney_element_local_assembly.py",
    "verify_whitney_mass_inverse_polynomial.py",
    "verify_whitney_mass_inverse_refinement.py",
    "verify_whitney_local_kkt_pencil.py",
    "verify_whitney_neighbour_constraints.py",
    "verify_whitney_constraint_dirac_bergmann.py",
    "verify_whitney_positive_metric_completion.py",
    "verify_whitney_first_class_conversion.py",
    "verify_whitney_finite_stiffness.py",
    "verify_whitney_stiffness_refinement.py",
    "verify_whitney_trace_stiffness.py",
    "verify_whitney_trace_stiffness_full.py",
    "verify_whitney_trace_stiffness_second_refinement.py",
    "verify_whitney_barycentric_shape.py",
    "verify_whitney_rank_edgewise_refinement.py",
    "verify_whitney_rank_edgewise_stiffness.py",
    "verify_whitney_trace_penalty_role.py",
    "verify_whitney_finite_stiffness_circle.py",
    "verify_whitney_central_flux_circle.py",
    "verify_whitney_3d_flux_canonicity.py",
    "verify_whitney_broken_feec.py",
    "verify_whitney_broken_feec_universality.py",
    "verify_whitney_projection_selection.py",
    "verify_whitney_commuting_recovery.py",
    "verify_whitney_dual_constraint_resolution.py",
    "verify_whitney_dual_resolution_edgewise.py",
    "verify_whitney_dual_resolution_all_k.py",
    "verify_whitney_reducible_poisson_brst.py",
    "verify_whitney_weak_first_class_hamiltonian.py",
    "verify_whitney_edgewise_continuum_dynamics.py",
    "verify_whitney_regge_continuum_transfer.py",
    "verify_tetrahedral_dirac_walk_bridge.py",
    "verify_tetrahedral_dirac_walk_robust.py",
    "verify_tetrahedral_dirac_walk_connectivity.py",
    "verify_h4_three_bond_walk.py",
    "verify_h4_three_bond_local_isotropy.py",
    "verify_h4_local_tight_frame.py",
    "verify_h4_matrix_coin_selection.py",
    "verify_invariant_spectrum.py",
    "verify_holographic_dimension.py",
    "verify_tower_spacetime.py",
    "verify_warped_spacetime.py",
    "verify_inner_fluctuations.py",
    "verify_primal_dual_triple.py",
    "verify_free_orbifold_arenas.py",
    "verify_multiplicity_mixing_J.py",
    "verify_canonical_bimodule_arena.py",
    "verify_q8_transplant_2I.py",
    "verify_orbifold_groupoid_pd.py",
    "verify_canonical_chamber_algebras.py",
    "verify_oriented_chamber_double.py",
    "verify_alpha_spectral.py",
    "verify_gravity.py",
    "verify_gravity_box4_full_hessian.py",
    "verify_hopf_spin2_tensor_carrier.py",
    "verify_hopf_kahler_induced_gravity.py",
    "verify_hopf_whitney_metric_selection.py",
    "verify_hopf_spectral_metric_selector.py",
    "verify_regge_de_rham_cone_selector.py",
    "verify_round_regge_interior_a2.py",
    "verify_round_regge_a2_interval.py",
    "verify_round_regge_spectral_action_sign.py",
    "verify_round_a2_transverse_hessian.py",
    "verify_finite_regge_a2_hessian.py",
    "verify_smooth_derham_a4_stabilization.py",
    "verify_gravity_hamiltonian_constraint_gap.py",
    "verify_gravity_metric_phase_space_canonicity.py",
    "verify_gravity_time_slab_canonicity.py",
    "verify_gravity_tent_move_regge.py",
    "verify_gravity_lorentzian_tent.py",
    "verify_gravity_lorentzian_asymmetric_tent.py",
    "verify_gravity_lorentzian_tent_regular_evolution.py",
    "verify_gravity_lorentzian_tent_legendre.py",
    "verify_gravity_global_tent_schedule.py",
    "verify_gravity_global_regge_orbits.py",
    "verify_gravity_global_regge_roots.py",
    "verify_gravity_global_boundary_legendre.py",
    "verify_gravity_lorentzian_volume_selection.py",
    "verify_rg_bootstrap.py",
    "verify_hopf_fibration_invariants.py",
    "verify_polytope_selection_intrinsic.py",
    "verify_edge_endomorphism_type.py",
    "verify_q8_abelian.py",
    "verify_s08_edge_fibration_uniformity.py",
    "verify_s08b_bridge_nogo.py",
    "verify_fibonacci_nonbinary_dynamics.py",
    "verify_exceptional_nonbinary_selection.py",
    "verify_nonnormal_c10_selection.py",
    "verify_noncentral_context_J.py",
    "verify_chamber_symmetry_sat.py",
    "verify_incidence_operator_enumeration.py",
    "verify_spectral_dimension_flow.py",
    "verify_dimension_reconciliation.py",
    "verify_chamber_b1_counterexample.py",
    "verify_chamber_b1_embedding_orbits.py",
    "verify_chamber_noncomm_no_go_refutation.py",
    "verify_chamber_rigidity_audit.py",
    "verify_math_to_physics_bridge.py",
    "verify_orbifold_incidence_route.py",
    "verify_hd600_holonomy.py",
    "verify_hd600_connection_space.py",
    "verify_hd600_heat_kernel_state.py",
    "verify_missing_link_selection.py",
    "verify_missing_link_modular_state.py",
    "verify_missing_link_refinement_scaling.py",
    "verify_inductive_spectral_dynamics.py",
    "verify_inductive_relativistic_gate.py",
    "verify_local_refinement_dynamics_gate.py",
    "verify_whitney_kahler_induction.py",
    "verify_whitney_circle_calibration.py",
    "verify_barycentric_shape_regular_gate.py",
    "verify_whitney_hopf_blind_enumeration.py",
    "verify_whitney_hopf_target_comparison.py",
    "verify_whitney_hopf_refinement_blind.py",
    "verify_whitney_hopf_refinement_comparison.py",
    "verify_hopf_vertical_infinity.py",
    "verify_smooth_hopf_refinement_blind.py",
    "verify_smooth_hopf_refinement_comparison.py",
    "verify_smooth_hopf_red_refinement_blind.py",
    "verify_smooth_hopf_red_refinement_comparison.py",
    "verify_hopf_symmetry_selector.py",
    "verify_hopf_sixth_order_selector.py",
    "verify_hopf_selector_action_gate.py",
    "verify_hopf_projector_cubic.py",
    "verify_hopf_box_projector_lift.py",
    "verify_hopf_box_cubic_selection.py",
    "verify_hopf_adjacency_baseline_cubic.py",
    "verify_hopf_higher_single_trace.py",
    "verify_hopf_fibration_label_operator.py",
    "verify_hopf_auxiliary_hessian.py",
    "verify_hopf_label_superselection.py",
    "verify_hopf_label_crossed_product.py",
    "verify_hopf_dynamical_superselection.py",
    "verify_hopf_full_hessian_spectral_enumeration.py",
    "verify_hopf_full_hessian_spectral_target.py",
    "verify_hopf_hessian_action_origin.py",
    "verify_hopf_hessian_ten_state_triple.py",
    "verify_hopf_hessian_c5_bimodule_nogo.py",
    "verify_hopf_hessian_c5_crossed_product_blind.py",
    "verify_hopf_hessian_crossed_krajewski_blind.py",
    "verify_hopf_hessian_crossed_real_form.py",
    "verify_hopf_hessian_binary_crossed_real.py",
    "verify_hopf_six_crossed_real_galois.py",
    "verify_hopf_six_galois_spectral_split.py",
    "verify_hopf_six_spectral_krajewski.py",
    "verify_hopf_six_equivariant_dirac.py",
    "verify_hopf_six_existing_operator_lift.py",
    "verify_hopf_six_equivariant_dirac_universal.py",
    "verify_hopf_six_w5_yukawa_enumeration.py",
    "verify_hopf_six_w5_yukawa_vacuum_gate.py",
    "verify_hopf_residual_symmetry_breaking.py",
    "verify_hopf_axis_orientation.py",
    "verify_hopf_chamber_crossed_product.py",
]

# Coverage guard.  On 2026-07-28 ten verifier files were sitting on disk
# unregistered -- among them verify_oriented_chamber_double.py (the chamber
# spectral triple), verify_alpha_spectral.py, verify_gravity.py and
# verify_rg_bootstrap.py, which underwrite whole sections of
# logic_chain_map.md.  A green "N/N" therefore did not mean what it looked
# like: those files could have been broken for weeks without anyone noticing.
# Fail loudly instead of silently under-reporting coverage.  If a verifier is
# deliberately excluded (too slow, experimental), add it to DELIBERATELY_SKIPPED
# with a reason rather than dropping it from the list.
DELIBERATELY_SKIPPED: dict[str, str] = {
    "verify_chamber_global_minimum.py": (
        "bounded branch-and-bound experiment; exits INCOMPLETE on its "
        "default resource limit"
    ),
    "verify_chamber_global_smt.py": (
        "bounded global SMT attempt; exits INCOMPLETE when its solver "
        "timeout is reached"
    ),
}

on_disk = {
    name for name in os.listdir(".")
    if name.startswith("verify_") and name.endswith(".py")
}
registered = set(scripts)
duplicates = sorted({name for name in scripts if scripts.count(name) > 1})
unregistered = sorted(on_disk - registered - set(DELIBERATELY_SKIPPED))
missing_file = sorted(registered - on_disk)

if unregistered or missing_file or duplicates:
    print("=" * 70)
    print("COVERAGE ERROR -- the suite does not cover the repository")
    print("=" * 70)
    for name in unregistered:
        print(f"  on disk but NOT registered : {name}")
    for name in missing_file:
        print(f"  registered but MISSING     : {name}")
    for name in duplicates:
        print(f"  registered more than once  : {name}")
    print("\nAdd them to `scripts`, or to DELIBERATELY_SKIPPED with a reason.")
    sys.exit(2)

results = {}
total_time = 0

# The chamber verifier contains a declared exhaustive C5 contraction census.
# A standalone run on 2026-08-09 crossed 300 s while remaining CPU-active and
# exited before 600 s.  Give that registered exhaustive
# calculation its measured class of allowance; retain 300 s for every other
# verifier so a general hang is still caught.
SCRIPT_TIMEOUTS = {
    "verify_chamber_symmetry_sat.py": 600,
    "verify_gravity_global_regge_roots.py": 1800,
    "verify_gravity_global_boundary_legendre.py": 1800,
}

for script in scripts:
    print("\n" + "=" * 70)
    print(f"RUNNING: {script}")
    print("=" * 70)
    t0 = time.time()
    try:
        result = subprocess.run(
            [sys.executable, script],
            capture_output=False,
            timeout=SCRIPT_TIMEOUTS.get(script, 300)
        )
        elapsed = time.time() - t0
        total_time += elapsed
        results[script] = "PASS" if result.returncode == 0 else "FAIL"
        print(f"\n  Completed in {elapsed:.1f}s -- {results[script]}")
    except subprocess.TimeoutExpired:
        results[script] = "TIMEOUT"
        print(f"\n  TIMEOUT after {SCRIPT_TIMEOUTS.get(script, 300)}s")
    except Exception as e:
        results[script] = f"ERROR: {e}"
        print(f"\n  ERROR: {e}")

print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
for script, status in results.items():
    print(f"  {script:40s} {status}")
print(f"\nTotal time: {total_time:.1f}s")

n_pass = sum(1 for s in results.values() if s == "PASS")
print(f"Result: {n_pass}/{len(scripts)} scripts completed successfully.")
if n_pass != len(scripts):
    sys.exit(1)
