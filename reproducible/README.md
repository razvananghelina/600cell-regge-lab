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

Each script is self-contained (no cross-imports) and prints PASS/FAIL for every claim.

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
| `verify_chamber_b1_counterexample.py` | Exact noncommutative `M2(C)+C^3` counterexample to B1 on the fixed chamber `D,J,gamma`: full-basis order zero/first order, nonzero forms, explicit orientability cycle, nondegenerate KO6 pairing, and connectedness | ~30s |
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
