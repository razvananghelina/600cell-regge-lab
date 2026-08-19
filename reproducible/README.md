# Reproducible Verification Code

**Paper:** "One Integer, Three Generations: Deriving Particle Physics from a_1 = 5"
**Author:** Razvan-Constantin Anghelina
**Version:** 3.9 (February 2026)

## Requirements

- Python 3.8+
- NumPy
- SciPy

```
pip install numpy scipy
```

## Scripts

Each script is directly executable and prints PASS/FAIL for every claim.  Some
later verifiers import and recheck an explicitly identified upstream verifier
before extending its certified construction.

| Script | What it verifies | Runtime |
|--------|-----------------|---------|
| `verify_coupling_constants.py` | alpha, alpha_s, sin^2(theta_W) from a_1=5 | <1s |
| `verify_spectrum_600cell.py` | 600-cell Laplacian spectrum, 9 eigenvalues in Z[phi], localization of 3/3' irreps in Galois-conjugate eigenspaces, spectral gap 4*sqrt(5) | ~30s |
| `verify_masses_and_mixing.py` | All 9 fermion masses (bare + norm-log corrected to 0.11% RMS), CKM angles, PMNS angles, CP phases | <1s |
| `verify_berry_phase.py` | Berry phase = 1/phi^4 over all 1200 faces, face-transitivity, 5 fiber values | ~60s |
| `verify_spectral_action.py` | Simplicial complex counts, boundary operators, Hodge decomposition 119+601=720, Betti numbers, exact discrete spectral coefficients `(c_0,c_1,c_2)=(2640,14880,55920)`, reduced triple `(11,62,233)`, Seeley-DeWitt Diophantine identity, gauge group from `A_5` | ~60s |
| `verify_discrete_scalar_response.py` | Exact static scalar-response theorem on the 600-cell: `B^+ d_0 = d_0 Delta_0^+`, every point source gives `h=d_0 Phi`, zero coexact content `d_1 h = 0`, and full recovery of `Phi` from `h` | <1s |
| `verify_mckay_chirality.py` | McKay graph = affine E8 (from 2I character table), bipartite chirality gamma_F = (-1)^{2j}, Casimir sum rules (sum C_2 = 26), Wilson line masses (9/9), fermionic action dimensions (H+ = H- = 39600), mass quantum number sum rules, rep ring Z/2 grading (45/45) | <1s |
| `verify_galois_kernel.py` | Galois kernel theorem: ker(b1*A_fiber - A) = rho_0+rho_1+rho_8 (dim=9), b1=6 unique, stability across Hopf fibrations, alpha*alpha'=1/(2pi) | ~60s |
| `verify_neutrino_masses.py` | Neutrino masses m_3=2*m_e/phi^35, splitting ratio r=alpha*phi^3, seesaw exponent n=35 (8 identities), PMNS angles, cosmological bounds | <1s |
| `verify_juno_comparison.py` | Pinned primary-source JUNO 2026 confrontation of three `1/45` correction scopes, exact masses/splittings/effective masses, NuFIT-with/without-SK deviations, DESI/KATRIN screens, and two-input chi-square ranking | <1s |
| `verify_muon_gminus2_null.py` | Dark-sector no-go for photon-coupled BSM loops: the dark electromagnetic coupling has no real root, `alpha_s' < 0`, hence the framework prediction is `a_mu(new physics) = 0` | <1s |
| `verify_polytope_uniqueness.py` | All 6 regular 4D polytopes tested against 7 SM criteria. 600-cell: 7/7, 120-cell (dual): 6/7, all independent alternatives: 0/7. Uniqueness of a1=5, N_gen from ring structure, Weinberg angle, alpha, mass hierarchy, mixing angles, anomaly cancellation | <1s |
| `verify_uniqueness_quantized.py` | Strengthened bounded-search uniqueness route `C1 + Q + C2`: irreducible node values in `Z[phi]`, edge Galois quantization, flatness sum, unique neutral assignment in `|t|<=15` | <1s |
| `verify_uniqueness_simple_edges.py` | Stronger bounded-search uniqueness route `C1 + S`: simple edge Wilson lines `dz in {0, +/-phi^r, +/-2phi^r}` already select the unique neutral assignment in `|t|<=15` | <1s |
| `verify_minimal_edge_lifts.py` | Constructive local replacement for `S`: each edge jump `Delta n` is lifted by the unique minimal-`L^1` solution of `5 da + 6 db = Delta n`, reconstructing the neutral assignment and making `C1`, `Q`, and `C2` emergent | <1s |
| `verify_branch_identification.py` | Resolves the neutral branch labels using affine-E8 leg geometry, bipartite McKay chirality, and the generation-2 up/down exponent formulas: length-2 BLACK endpoint = top, short-leg WHITE endpoint = bottom | <1s |
| `verify_global_uniqueness_constructive.py` | Full bounded-search-free constructive chain for the physical assignment: main-chain minimal lifts, branch identification, top from the generation-2 up-type branch jump, and bottom forced by the prime-sector Galois relation `z_b = phi*sigma(z_t)` | <1s |
| `verify_neutral_vacuum_scale.py` | The neutral electroweak selector `n=25`: unique `ker(A)` block of multiplicity `25`, vanishing generator character, IR projector `Tr(e^{-5A^2})≈25`, and `Box` refinement `12^(5) ⊕ (6/phi)^(10) ⊕ (-6phi)^(10)` | <1s |
| `verify_edge_gauge_spectrum.py` | Edge-space gauge spectrum from `Box_1`: `ker(Box_1)=13=rho_0+2rho_5`, 12 gauge modes localized on fiber edges, exact Hodge eigenvalues `C=5`, `B=16/5`, and the `A_5` decomposition `12=1+3+3'+5 => 1+3+8` | ~60s |
| `verify_gauge_continuum_map.py` | Controlled gauge continuum map on the Hopf base: the 12 fibers form an icosahedron, the sampled scalar harmonics reproduce `S^2` exactly up to `l=2` (`1+3+5`), the remaining `3'` is the first alias sector, and a quadratic polynomial in the base Laplacian reproduces the continuum Casimir `l(l+1)` on that low-mode sector | <1s |
| `verify_low_mode_sampling_intertwiner.py` | Exact band-limited bridge theorem: evaluation of `H_0+H_1+H_2` on the 12-point Hopf base is injective, has exact reconstruction, intertwines the round-sphere scalar Casimir with a quadratic polynomial in the base Laplacian, and separates the 3-dimensional alias sector | <1s |
| `verify_vector_sampling_intertwiner.py` | Vector bridge on the Hopf base: geodesic edge integration samples exact/coexact `l=1,2` vector harmonics into a rank-16 subspace, a cubic polynomial in the edge Hodge Laplacian exactly intertwines the continuum 1-form spectrum, the 14-dimensional alias complement is resolved, and the four norm factors test the isometry obstruction | <1s |
| `verify_a5_equivariant_brackets.py` | Complete real classification of `A5`-equivariant antisymmetric brackets on `3'+5`: Jacobi variety, compact/split/degenerate classes, explicit `3'` embedding, Killing normalization, canonical edge-metric constants `(20,16)`, and metric-compatible branch classification | <1s |
| `verify_matter_trace_indices.py` | Exact gauge-index target `8:5:2`, one-generation SM benchmark, finite inventory of existing candidate spaces, and negative test of the only derived common module `(0,2,3)` | <1s |
| `verify_matter_module.py` | Exhaustive 15/16-dimensional nonabelian block counts, exact SM trace/anomaly benchmark, and no-go tests for McKay exponent, C10 residue, and `Z[phi]` unit gradings as hypercharge | <1s |
| `verify_matter_functor.py` | Route-specific chiral-functor audit: McKay `16+14` KO6 obstruction, Galois grading compatibility, C10 winding alias/no-selection theorem, exact M15/M16 anomaly factorizations, and the exact negative character test of even-16 versus diagonal M16 | <1s |
| `verify_galois_doubling.py` | Exact A5/2I outer-versus-Galois character action, real-color conjugation negative, doubled KO sign variants, scoped node-algebra first-order failure, and the 240-dimensional bidirected-Hom decomposition | <1s |
| `verify_bimodule_krajewski.py` | Exhaustive 5x5 Krajewski multiplicity-matrix obstruction on the 30/60-dimensional McKay node spaces, exact doubled Schur commutants, the node-space M3 no-go, and the edge/Hom abstract-embedding check | <1s |
| `verify_edge_matter.py` | Complete bidirected affine-E8 edge census, exact two-sided commutant and tensor factorizations, canonical KO6 orientation/adjoint signs, and the first-order-plus-oddness Dirac no-go | <1s |
| `verify_bratteli_tower.py` | Exact Fibonacci and rooted affine-E8 McKay floors through level 12, Smith/stable dimension-group certificates, canonical trace ranges, parity-separated matter seeds, and the consecutive-floor first-order witness | <1s |
| `verify_warped_spacetime.py` | Exact warped mode decomposition, frozen `N=8,16,24` spectral-dimension scan for golden/McKay/no warp, Dirichlet/Neumann sensitivity, conformal-boundary recovery, and exact R4/H4 heat-kernel benchmarks | ~30s |
| `verify_segregation_theorem.py` | Unique involution and complete odd-subgroup classification in `SL(2,5)`, exact `C3/C5/C1` odd-operator dimensions on `W`, and canonical cyclic Krajewski escape tests | <1s |
| `verify_dirac_selection.py` | Conjugacy of all `C3` escapes, the complete legal `16+14` Dirac-block count `d0=132`, KO6 reality count, generic gauge/scale moduli dimension `122`, polynomial critical-family obstruction, and exact `C3` isotypic census | <1s |
| `verify_missing_link_selection.py` | Exact rank-44 projective cometric split and spectral-functional reduction `Tr f(K_h+tQ)=Tr_h f(K_h)+44f(t)`; positive-moment, heat-cutoff, and coefficient-ratio selection tests | <1s |
| `verify_missing_link_modular_state.py` | Exact finite modular-state audit: trivial Haar/trace flow, arbitrary density-ratio clock for nontracial states, and the faithful `M2` GNS image/commutant obstruction to selecting `M4` | <1s |
| `verify_missing_link_refinement_scaling.py` | Exact full-barycentric f-vector growth, dominant factor 24, and conditional level-weight spectral law `d=log(24)/log(b)`, showing that compact resolvent/convergence does not select dimension | <1s |
| `verify_inductive_spectral_dynamics.py` | Actual 600-root/24-branch inductive path dynamics: root dual graph, exact conditional-expectation Laplacian, Markov heat flow, volume-selected dimension-three Dirac, critical KMS balance, and ordinary-versus-twisted Cuntz commutator gate | ~5s |
| `verify_inductive_relativistic_gate.py` | Exact inertia--mass--causality control plus the 24-child barycentric incidence audit: the hierarchical generator couples all 276 child pairs, including 240 non-neighbours, and has no spatial dispersion band | <1s |
| `verify_local_refinement_dynamics_gate.py` | Exact one-tetrahedron face-local refinement gate: Galerkin compression uniquely fixes Laplacian scale 4, while directional vertical leakage disproves exact operator intertwining | <1s |
| `verify_whitney_kahler_induction.py` | Exact all-degree Whitney/FEEC induction on a tetrahedron and its barycentric subdivision: `d` commutation, `L2` isometry, metric Kähler--Dirac compression, adjoint leakage, and finite generalized spectra | ~5s |
| `verify_whitney_circle_calibration.py` | Known-answer unit-circle calibration of consistent Whitney Dirac induction and second-order spectral convergence; exposes the exact-induction versus finite-speed tradeoff under mass lumping | ~2s |
| `verify_barycentric_shape_regular_gate.py` | Exact repeated-flag witness that tetrahedral barycentric refinement has volume ratio `24^-n` but unbounded affine condition at least `2^n`, blocking standard shape-regular FEEC convergence | <1s |
| `verify_whitney_hopf_blind_enumeration.py` | STEP-1 target-blind enumeration of all raw/generalized Whitney scalar fiber/cross spectra for all six discrete Hopf fibrations; writes the committed JSON before any bootstrap/speed comparison | ~5s |
| `verify_whitney_hopf_target_comparison.py` | STEP-3 comparison against independently derived `a1=5`: common four-dimensional generalized gap eigenspace, exact mass cancellation, nine-dimensional weighted-difference kernel, and preservation-versus-selection audit | ~2s |
| `verify_whitney_hopf_refinement_blind.py` | Target-blind first-refinement extension of the unique-per-tetrahedron Hopf edge to positive local tensors; exact Galerkin compression and low generalized spectra on 120/2640 scalar carriers | ~2m |
| `verify_whitney_hopf_refinement_comparison.py` | Post-preregistration kill test: coarse Whitney gap ratio `5.0000000005` drifts to `5.3388401713` under the coefficient-free local-tensor refinement despite `1e-15` Galerkin residuals | <1s |
| `verify_hopf_vertical_infinity.py` | Exact round-`S3` Hopf spectrum audit: the separated vertical operator has infinite kernel/positive multiplicities, while every positive combined horizontal+vertical operator has only the constant zero mode | <1s |
| `verify_smooth_hopf_refinement_blind.py` | Target-free projected-barycentric `P1` discretization of the smooth Hopf field with true rank-1/rank-2 orthogonal projectors; records canonical charged/base Ritz modes and low spectra before continuum comparison | ~2m |
| `verify_smooth_hopf_refinement_comparison.py` | Post-commit continuum calibration: smooth Hopf charged modes converge to `(1,2,3)`, base pullbacks to `(0,8,8)`, and the combined operator retains one zero; separates the genuine Hopf split from the old combinatorial factor five | <1s |
| `verify_smooth_hopf_red_refinement_blind.py` | Target-free projected tetrahedral `1-to-8` red refinement through `120/840/6480` nodes; shortest central-octahedron diagonal, element-quality audit, smooth-Hopf canonical Ritz modes and low spectra | ~1m |
| `verify_smooth_hopf_red_refinement_comparison.py` | Post-commit two-level gate: exact closed topology, stable element quality, canonical smooth-Hopf Ritz errors decreasing twice, one combined zero mode; second-order rate and `3/5` low bands labelled pattern | <1s |
| `verify_hopf_symmetry_selector.py` | Exact `Q(sqrt(5))` enumeration of both handed `qH/Hq` C10 Hopf-fibration classes; conjugation orbit, invariant coefficient space, equiangular projector Gram matrix and tight-frame isotropy | <5s |
| `verify_hopf_sixth_order_selector.py` | Exact icosahedral moment and Gröbner audit: radial equal-weight moments through degree four, exhaustive degree-six C10/C4/C6 critical orbits, and the conditional six-vacuum Hopf selector | <10s |
| `verify_hopf_spin2_tensor_carrier.py` | Exact six-projector tight frame for `Sym^2_0(R^3)` and its left/right homogeneous TT lifts on round `S3` | <1s |
| `verify_hopf_kahler_induced_gravity.py` | Direct `SU(2)` curvature formula and ordinary full-de-Rham heat `A2` Hessian `(8/3)I5` on the Hopf spin-two carrier | <5s |
| `verify_hopf_whitney_metric_selection.py` | Exact round-versus-fixed-Regge distinction, full admissible metric continuum, and all-degree Whitney refinement-isometry no-go | <5s |
| `verify_hopf_spectral_metric_selector.py` | Global fixed-volume round selection by smooth de Rham `A2`, exact scale/time covariance, heat-order reversal, `31/11` cutoff audit, and common-endpoint scope gate | <1s |
| `verify_regge_de_rham_cone_selector.py` | Exact full-de Rham cone coefficient on the fixed 600-cell Regge metric; link/domain audit and equal-volume round endpoint comparison | <1s |
| `verify_round_regge_interior_a2.py` | Preregistered full affine round--Regge path: exact bulk/face/edge formulas, five-order Duffy convergence, both endpoint controls, and target-free 201-point monotonicity audit (PATTERN, not an interval theorem) | ~25s |
| `verify_round_regge_a2_interval.py` | Preregistered Arb certificate for the continuous round--Regge path: uniform complex analyticity bounds, degree-18 Taylor cover, validated spatial orders 16/20/24, endpoint stationarity, and independent radial/Duffy control | ~45m (4 workers) |
| `verify_round_regge_spectral_action_sign.py` | Exact Mellin-moment audit: every standard positive cutoff gives a positive `A2` weight, while a Gaussian cutoff family and an exact heat-order reversal prove that its magnitude, scale and full finite-action minimum remain unselected | <1s |
| `verify_round_a2_transverse_hessian.py` | Preregistered exact conformal-Hessian audit on smooth round `S3`: positive Hopf directions coexist with a negative non-gauge `l=2` direction, so ordinary de Rham `A2` is a saddle outside the homogeneous sector | <1s |
| `verify_finite_regge_a2_hessian.py` | Preregistered full 720-edge automatic-differentiation Hessian of the exact equal-volume conical de Rham `A2` at the equilateral 600-cell, including quotient inertia, edge-stabilizer probes and the complete vertex-quadratic conformal carrier | ~15s |
| `verify_smooth_derham_a4_stabilization.py` | Preregistered exact smooth full-de Rham `A4` trace reduction and conformal Hessian: positive on every non-gauge scalar mode, but suppressed by `Lambda^-2` relative to the unstable `A2`; the finite Gaussian truncation threshold remains unselected | <2s |
| `verify_gravity_hamiltonian_constraint_gap.py` | Preregistered structured-certificate audit of the current dynamical candidates: the unitary tick is fixed-background, Whitney copy constraints are second class, and existing first-class conversions provide no local metric Hamiltonian constraint | <1s |
| `verify_gravity_metric_phase_space_canonicity.py` | Exact full-H4 orbital audit on the 720-edge Regge cotangent phase space: 47 symmetric kinetic parameters, with four surviving nearest-neighbour locality and three surviving common-tetrahedron locality | <2s |
| `verify_gravity_time_slab_canonicity.py` | Exact cone/product/Pachner audit: canonical two-boundary CW cylinder, H4 obstruction to a vertex-only simplicial product, 115,200 barycentric four-flags, and exhaustive 8! negative comparison with the existing robust walk | <3s |
| `verify_gravity_tent_move_regge.py` | Exact local icosahedral tent carrier and full Euclidean Regge boundary-action audit: fivefold incidence gives the conditional static root `t/a=phi^-1`, while unrestricted final length, volume terms and overlapping move orbits block a physical tick reading | <3s |
| `verify_gravity_lorentzian_tent.py` | Exact symmetric Lorentzian tent no-go: every admissible timelike-hinge angle is below the fivefold target, so no zero-volume stationary pole exists in the one-orbit final-boundary family | <3s |
| `verify_gravity_lorentzian_volume_selection.py` | Frozen six-class census and spectral-action audit: the current repository does not select a Lorentzian tent volume coefficient, while Gaussian cutoff freedom can fit every positive value | <2s |
| `verify_gravity_lorentzian_asymmetric_tent.py` | Preregistered Arb certificate of a strictly admissible asymmetric zero-volume stationary pole, including independent Minkowski normals and the target-found/non-selected boundary-data audit | <5s |
| `verify_gravity_lorentzian_tent_regular_evolution.py` | Preregistered Arb proof that the asymmetric internal pole equation has nonzero pole Hessian and defines a unique local `rho(q)` for nearby boundary data; this certificate feeds the later local Legendre-map audit | <3s |
| `verify_gravity_lorentzian_tent_legendre.py` | Corrected plus-branch complex Regge action on all 92 tent hinges, real pre/post star momenta, and independently confirmed full-rank `12x12` on-shell mixed Legendre block; common link and adjacent moves remain open | ~12s |
| `verify_gravity_global_tent_schedule.py` | Exact five-phase global tent census: `alpha=24`, `chi=5`, 25 maximal 24-cells, ten `H4`-equivalent covers, two ordered parity orbits, all 1200 staircase schedules, and a coherent but nonstationary Lorentzian product slab | ~12s |
| `verify_gravity_global_regge_orbits.py` | Corrected full-versus-orbit Lorentzian Regge audit: 840 internal edges reduce honestly to 35 orbits, 100 simplex orbits reproduce the 2400-simplex action and gradients, the regular Hessian has rank 35, and the two phase parities have distinct quadratic spectra | ~15s |
| `verify_gravity_global_regge_roots.py` | Preregistered equal-boundary global search: six identical starts in each 35-variable phase-parity sector, zero validated roots, seven artificial-box contacts, one causal/branch contact, four iteration limits, and exact parallel reproduction of all terminal residuals | ~10min (8 workers) |
| `verify_gravity_global_boundary_legendre.py` | Preregistered variable-final-boundary audit: exact 100-orbit/full-2400-simplex agreement for a 65-variable Regge action, stable ranks `rank(J)=rank(J_internal)=35`, `rank(J_final)=30`, and independently checked real pre/post boundary momenta | ~10min (8 workers) |
| `verify_gravity_600cell_dust_full_anisotropic_legendre_rank.py` | Preregistered full `2280`-edge, `2T`-resolved canonical Legendre-rank census at the accepted dust tick; four independently stepped Hessians, direct action controls, seven frozen representation blocks, and target-free regular/open/degenerate classification | potentially hours (8 workers) |
| `verify_gravity_600cell_dust_full_lapse_schur.py` | High-precision `2T` Schur audit of all 120 geometry-selected pole/lapse directions: 100-decimal representative kernels, 80-decimal Flint ball solves, calibrated rank classification, and frozen vertex-lapse subspace comparison | potentially hours |
| `verify_gravity_600cell_dust_lapse_stiffness_origin.py` | Disclosed post-result decomposition of the regular lapse Schur scale into its exact affine dust shift and gravitational remainder; quantifies scalar consistency and the near cancellation without promoting it to a Hamiltonian constraint | <1s |
| `verify_gravity_600cell_published_dust_control.py` | External reproduction of the published De Felice--Fabri time-symmetric dust sandwich: both schedule parities solve all 35 complete one-slab orbit equations; includes recorded binary64 cancellation failure and corrected 60-decimal action differences | ~6min (8 workers) |
| `verify_gravity_600cell_dust_implicit_jacobian.py` | Preregistered three-step Jacobian/Hessian audit at the published dust sandwich, with branch, symmetry, cross-reciprocity and independent 60-decimal weakest-curvature controls; classifies local regularity without treating full rank as a PASS target | ~15min (8 workers) |
| `verify_gravity_600cell_dust_lapse_schur.py` | Preregistered 80-decimal correction restricted to the five unresolved pole/lapse directions: canonical Schur lifts, 15 action-only curvatures, collective/relative split and a three-point lapse-family control | ~12min (8 workers) |
| `verify_gravity_600cell_dust_exact_lapse_path.py` | Preregistered 100-decimal direct audit of the published path `rho=tau^2 exp(t)`, `q=l0^2-rho`: eleven full residuals, exact-path action constancy/curvature, and the frozen four relative pseudo-constraint modes | ~7min (8 workers) |
| `verify_gravity_600cell_dust_gauge_quotient.py` | Preregistered gauge-fixed linear response: exact lapse-tangent quotient, 90-decimal mixed-action reconstruction of the boundary compatibility row, and the 34-coordinate response to all compatible final-boundary directions | ~12min (8 workers) |
| `verify_hopf_six_galois_spectral_split.py` | Exhaustive two-handed fibre-incidence audit: all 1,440 internal edges select one D5 inverse class, exact golden spectral projectors separate the two M12 blocks, and the Galois automorphism maps edges to distance-two chords | <5s |
| `verify_hopf_six_spectral_krajewski.py` | Eight-reading audit of the incidence-labelled four-node real algebra: 936-state off-diagonal enveloping carrier, exhaustive matrix-unit order zero, explicit KO6 orientation cycle, unimodular Poincare form, and scoped legal-Dirac-position census | <10s |
| `verify_hopf_six_equivariant_dirac.py` | Exact induced-A5 Hom Gram matrix and 4/8 eligibility census; all 32 normalized equivariant rook operators pass KO6/first-order/nonzero-form gates but fail connectedness with commutant dimension 109 or 141 | <10s |
| `verify_hopf_six_existing_operator_lift.py` | Frozen seven-family provenance audit: five existing operators have no faithful map to the 936-state carrier, while all left/right crossed-product-generated candidates preserve central cells and have zero odd projection | <1s |
| `verify_hopf_selector_action_gate.py` | Current-action kill gate: exact A5 stabilizer of the valid chamber B1 embedding, its two-doublet `M2` one-form block, and the `D^4` polynomial-degree ceiling excluding the sixth-order Hopf selector | <10s |
| `verify_chamber_b1_counterexample.py` | Exact noncommutative `M2(C)+C^3` counterexample to B1 on the fixed chamber `D,J,gamma`: full-basis order zero/first order, nonzero forms, explicit orientability cycle, nondegenerate KO6 pairing, and connectedness | ~30s |
| `verify_chamber_b1_embedding_orbits.py` | Exact finite non-selection certificate for the B1 embedding: reconstructs `Aut(S)=A5`, verifies two disjoint free support orbits, rebuilds every B1 gate on both, and proves the residual 9-real-dimensional internal `M2` ambiguity; explicitly does not claim a complete census | ~10s |
| `verify_gravity_600cell_cellular_frustum_relative_poincare.py` | Exact Poincare decomposition of the six cellular-frustum flexes: direct-kernel equality, static-versus-expanding Lorentz projection, symbolic translation determinant, and origin/frame/labelling covariance controls | <10s |
| `verify_gravity_600cell_cellular_frustum_relative_poincare_covariance_correction.py` | Preserved-failure correction using exact Poincare intertwiners and observer-normal stabilizers; distinguishes invariant Lorentz covariance from the noninvariant coordinate rotation/boost split | <10s |
| `verify_gravity_600cell_cellular_frustum_relative_poincare_adversarial.py` | Independent irregular-frustum replication using a redundant 20-parameter affine-Lorentz system, direct polynomial kernel equality, a different rational boost and metric-sign controls | <10s |
| `verify_gravity_600cell_two_frustum_face_gluing.py` | Exact two-frustum fixed-frame gluing: positive one-dimensional full-Poincare face-stabilizer control versus the constrained relative mode inside the accepted six-flex kernels | <10s |
| `verify_gravity_600cell_two_frustum_face_gluing_adversarial.py` | Independent irregular five-vertex union audit: direct 14-by-20 squared-length Jacobian, redundant full-Poincare face control, metric-sign and exact-boost transport | <10s |
| `verify_gravity_600cell_global_flex_holonomy.py` | Exact candidate global-rigidity gate: golden-field fivefold edge incidence, rational Regge deficit rotations, affine Poincare adjoints and common fixed spaces of two nonparallel hinge holonomies | <30s |
| `verify_gravity_600cell_global_flex_holonomy_adversarial.py` | Complete exact dual-complex audit: 2400 vertex-matched face transitions, all 720 actual five-tetrahedron hinge loops, and six-loop closure of the local flex seed with reversal/relabel controls | potentially minutes |
| `verify_gravity_600cell_variable_face_connection.py` | Exact two-frustum correction gate allowing the uniquely derived lower-triangle Poincare stabilizer to vary the face transition, with fixed-frame, Lorentz-covariance and metric-sign controls | <10s |
| `verify_gravity_600cell_variable_face_connection_adversarial.py` | Independent direct polynomial audit on two irregular reflected carriers: 52-variable gluing Jacobian with redundant Lorentz matrix entries, frozen-transition and pointwise-face controls | <30s |
| `verify_gravity_600cell_global_variable_face_closure.py` | Complete exact variable-transition closure: full Lorentzian lateral-face maps, 1200 rank-five body-hinge blocks and dual-prime finite-field certificates on the 6000-by-3600 global matrix | potentially minutes |
| `verify_gravity_600cell_canonical_data_admissibility.py` | Exact forced variable-face audit of the 840 upper-edge/strut tangent data: local rank-ten right inverses, complete augmented compatibility ranks, two rational homothetic witnesses and graph/convention attacks | potentially minutes |
| `verify_gravity_600cell_canonical_data_carrier.py` | Target-disclosed exact falsification of the proposed 240-dimensional vertexwise radial-scale plus normal-lapse carrier, with alternate-graph and wrong-carrier controls | potentially minutes |
| `verify_gravity_600cell_canonical_data_projection.py` | Target-blind modular census of the compatible-data kernel's intersections with and projections onto the 720 upper-edge and 120 strut coordinates, with synthetic quotient, graph and convention attacks | potentially minutes |
| `verify_gravity_600cell_static_vertex_gradient.py` | Exact rational exhaustion of the 119-dimensional static closure kernel by gradients of continuous P1 vertex scalars modulo constants, with discontinuity, relabelling and metric-sign attacks | potentially minutes |
| `verify_gravity_600cell_static_vertex_gradient_adversarial.py` | Independent static block decomposition using only spatial face transports: full-rank curved rotations plus tangentially continuous translations exhausted by 119 P1 gradients | potentially minutes |
| `verify_gravity_600cell_static_gradient_prism_shift_reconciliation.py` | Exact cellwise and facewise intertwiner proving that the newly derived 119 static-gradient translations are the previously certified prism-shift carrier in Cartesian rather than covector coordinates | potentially minutes |
| `verify_kahler_dirac_local_tick.py` | Target-blind signed Grover--Szegedy lift of the 2640-state Kähler--Dirac incidence: exact local reflections, a 14880-state unitary tick, one-Hasse-edge cone, exact spectral map, ballistic cycle calibration, and the 9600-dimensional extra-carrier audit | ~2s |
| `verify_kahler_dirac_tick_refinement.py` | Exact refinement gate for the local tick: dyadic circle quasienergies scale perfectly, the barycentric top block remains metric-uniform, but lower Whitney adjoints acquire 10 and 12 off-incidence entries, restricting the walk to combinatorial kinematics | <1s |
| `verify_weighted_szegedy_metric_nogo.py` | Universal support no-go for all rank-one weighted/phased Szegedy coins on the original Hasse arcs: 22 exact Whitney entries cannot be produced and all require incidence depth at least three | <1s |
| `verify_incidence_polynomial_metric_nogo.py` | Exact symmetry no-go for the algebra generated by signed incidence and cochain-degree projectors: all words retain S4 symmetry, while the barycentric Whitney operator has trivial stabilizer; also certifies complete cancellation of naive three-step paths | <1s |
| `verify_tetrahedral_dirac_walk_bridge.py` | Exact H4 chamber-carrier audit for the published tetrahedral Dirac walk: canonical 14,400-node four-coloured bipartite geometry, but a source-count theorem shows the economical four-amplitude shift is two-to-one under natural chamber handedness | <1s |
| `verify_tetrahedral_dirac_walk_robust.py` | Literal Appendix-B three-swap audit on the doubled 115,200-state H4 chamber carrier: all stages compose to a local permutation, while the separately printed published expansion differs in 28,800 positions and is non-bijective | <1s |
| `verify_tetrahedral_dirac_walk_connectivity.py` | Universal connectivity no-go for the literal robust transplant: colours 2 and 3 form 1,440 decagons, and every word with arbitrary chamber-local coins preserves their exact orbit projectors | <1s |
| `verify_h4_three_bond_walk.py` | Exhaustive six-order support audit for the three consecutive H4 Coxeter bonds: all local templates are permutations and all 6/6 periodic schedules are strongly connected on the 172,800-state time-expanded carrier | ~2s |
| `verify_h4_three_bond_local_isotropy.py` | Calibrated one-period tangent-moment gate for the connected three-bond walk: exact fixed coin and maximally mixed spin input give 0/6 isotropic schedules; the designated covariance eigenvalue ratio is 11.69 | ~1s |
| `test_holographic_rg.py` | Exploratory go/no-go test for emergent 4D holography from the 600-cell boundary spectrum using counting-function scaling and heat-flow effective dimension; returns `CONTINUE`, `STOP`, or `INCONCLUSIVE` | ~10s |

## Running

```bash
python verify_coupling_constants.py
python verify_spectrum_600cell.py
python verify_masses_and_mixing.py
python verify_berry_phase.py
python verify_spectral_action.py
python verify_discrete_scalar_response.py
python verify_mckay_chirality.py
python verify_galois_kernel.py
python verify_neutrino_masses.py
python verify_muon_gminus2_null.py
python verify_polytope_uniqueness.py
python verify_uniqueness_quantized.py
python verify_uniqueness_simple_edges.py
python verify_minimal_edge_lifts.py
python verify_branch_identification.py
python verify_global_uniqueness_constructive.py
python verify_edge_gauge_spectrum.py
python verify_gauge_continuum_map.py
python verify_low_mode_sampling_intertwiner.py
python verify_vector_sampling_intertwiner.py
python verify_a5_equivariant_brackets.py
python verify_neutral_vacuum_scale.py
```

Or run all at once:
```bash
python run_all.py
```

## Expected Output

Each script prints a summary table ending with:
```
TOTAL: X/Y tests PASSED
```

All tests should pass. Any failure indicates either a bug in the verification code or a discrepancy with the paper's claims.

## What These Scripts Do NOT Verify

- The cosmological constant formula (pattern, not derivation)
- The dark matter abundance ratio (speculative)
- The 3D->4D continuum limit of the spectral triple
- The full vector-harmonic / nonabelian gauge continuum completion
- Gauge-sector trace normalizations `(8/15, 1/3, 2/15)`
- The physical field-normalization axiom needed to identify `4*a1*phi^4` with `1/alpha_0`
- A canonical vector isometry: edge integration has unequal norm factors on the four resolved `l=1,2` exact/coexact blocks

## License

MIT License. Use freely for verification purposes.
