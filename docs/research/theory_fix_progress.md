# Theory Fix Progress

## Session 2026-07-27, fifth session: final canonical bimodule boundary

- Rechecked `/tmp/science-python-deps`: `numpy 2.5.1`, `scipy 1.18.0`,
  and `sympy 1.14.0` were present.
- Before project writes, re-established the registered baseline:
  `PYTHONPATH=/tmp/science-python-deps python3 reproducible/run_all.py`
  completed in `297.6 s` with
  `Result: 52/52 scripts completed successfully.`
- Added `canonical_bimodule_arena.md` and
  `reproducible/verify_canonical_bimodule_arena.py`; registered the verifier
  as script 53.
- **REFUTED setup identification:** for non-simple
  `A=C[2I]`, the abstract `A tensor A^op` has dimension 14,400 but is not
  `End_C(A)=M120(C)`.  Two-sided multiplication on `A` has exact image
  dimension `sum d_i^4=2628` and kernel dimension 11,772.
- **DERIVED bimodule data:** flip-star has `J^2=+1` and exchanges the two
  commuting actions, so order zero is exact.
- **DERIVED adjoint decomposition:** in McKay order,
  `A_Ad=(9,0,7,0,9,0,6,0,7)` and
  `(A tensor A^op)_Ad=(296,0,736,0,1192,0,932,0,736)`.
  The central element acts trivially, so the segregation theorem's
  opposite-central-character hypothesis is absent.
- **DERIVED adjacency identification:** the twelve neighbors of the identity
  are the inverse-closed trace-`phi` conjugacy class.  Their sum `c` is
  self-adjoint and central; the graph Laplacian is `12-c`.
- **DERIVED candidate verdict:** the minus two-factor lift has `JD=-DJ`; the
  plus lift has `JD=+DJ`.  Class-sum/Laplacian lifts and the equivariant
  Hopf vertex-Box lifts satisfy order zero and first order but commute with
  the represented left algebra, so all inner one-forms vanish.
- **DERIVED grading obstruction:** adjoint and McKay central parities commute
  with the candidates.  The explicit triangle proves the 600-cell graph is
  nonbipartite, excluding a vertex-sign grading; form parity has no derived
  transport to dimension 14,400.
- **DERIVED SM-corner scope:** both Galois choices have support rank
  `14*120=1680`, not 14,400.  Restriction leaves the zero-form and grading
  failures unchanged and remains nonunital.
- **CLOSING BOUNDARY:** no derived `(D,gamma,J)` row passes KO6 oddness with
  nonzero fluctuations.  Per the binding stopping rule, no larger arena is
  proposed.  Existing multiplicity-mixing `J` varieties and specified
  twisted orbifold fibers remain **OPEN** on already authorized arenas.
- The physics gate remains closed; no `Y`, anomaly, multiplet, generation,
  Yukawa, or frozen mass-target trial was activated.
- Standalone verifier completed `21/21` checks.
- Final registered suite:
  `PYTHONPATH=/tmp/science-python-deps python3 reproducible/run_all.py`
  completed in `301.8 s` with
  `Result: 53/53 scripts completed successfully.`
- No PDF build was attempted.

## Session 2026-07-27, fourth session: multiplicity-mixing J boundary

- Rechecked `/tmp/science-python-deps`: `numpy 2.5.1`, `scipy 1.18.0`,
  and `sympy 1.14.0` were present.
- Before project writes, re-established the registered baseline:
  `PYTHONPATH=/tmp/science-python-deps python3 reproducible/run_all.py`
  completed in `631.0 s` with
  `Result: 51/51 scripts completed successfully.`
- Added `multiplicity_mixing_J.md` and
  `reproducible/verify_multiplicity_mixing_J.py`; registered the verifier as
  script 52.
- **DERIVED classification correction:** order-zero embeddings are classified
  by nonnegative `9 x 9` bimodule matrices `k_ij` satisfying both weighted
  margins `sum_j d_j k_ij=d_i m` and
  `sum_i d_i k_ij=d_j m`.  One nine-entry global multiplicity vector does
  not classify them.
- **DERIVED moduli dimensions:** a type `K` has commutant-conjugacy orbit
  dimension
  `sum_i ((d_i m)^2-sum_j k_ij^2)`; its fixed-image implementing-unitary
  torsor has dimension `120m^2`.  For the diagonal types the orbit dimensions
  are 53,724 (`m=22`) and 214,896 (`m=44`).
- **DERIVED necessary reality condition:** because all complex `2I` irreps
  are self-conjugate, a real structure requires a transpose-symmetric
  bimodule type, with additional real/symplectic FS constraints.
- **DERIVED geometric-stratum verdict:** every repository-derived
  coefficient/Galois/inversion/Hodge/star composition still fails at least
  one of KO6 grading, order zero, first order, or `JD=+/-DJ`.
- **STRUCTURAL/OPEN:** the arbitrary continuous `U` varieties have not been
  exhaustively solved.  Neither existence nor a 2I-specific universal no-go
  is claimed.  The Q8 result forbids promoting the geometric negative to an
  arbitrary-`J` theorem.
- **DERIVED algebra scope:** the two SM-type block choices have corner rank
  14 per regular copy, not 120.  They are nonunital on the full arena until
  an action on the other six Wedderburn sectors is supplied.
- Added explicit Q8/multiplicity-mixing scope corrections to
  `bimodule_krajewski_result.md`, `edge_matter_krajewski.md`,
  `preprojective_matter.md`, and `kahler_dirac_matter.md`.
- The physics gate remains closed: no `Y`, Route C, multiplet, generation,
  Yukawa, or frozen `Z[phi]` comparison was activated; new look-elsewhere
  count is zero.
- Standalone structural verifier completed `22/22` checks.
- Final registered suite:
  `PYTHONPATH=/tmp/science-python-deps python3 reproducible/run_all.py`
  completed in `300.1 s` with
  `Result: 52/52 scripts completed successfully.`
- No PDF build was attempted.

## Session 2026-07-27, third session: free claim refuted and orbifold arena

- Re-established the untouched registered baseline before project writes:
  `PYTHONPATH=/tmp/science-python-deps python3 reproducible/run_all.py`
  completed in `280.7 s` with
  `Result: 50/50 scripts completed successfully.`
- Added `free_arena_nogo_theorem.md`, `orbifold_arena.md`, and
  `reproducible/verify_free_orbifold_arenas.py`; registered the verifier as
  script 51.
- **REFUTED proposed theorem:** for
  `A=R(C[G]) tensor I_m`, the commutant is
  `L(C[G]) tensor M_m(C)`, not merely `L(C[G])`.  Order zero does not force
  a factor-preserving side exchange when `m>1`; the proper-subalgebra
  commutant is larger still.
- **DERIVED exact counterexample:** on the free arena
  `C[Q8] tensor C^16`, a factor-swap/inversion antiunitary and a noncentral
  right-convolution Dirac give KO6 signs `(+,+,-)`, exact order zero, exact
  first order, full free-left equivariance, and nonzero inner one-forms.
  This is finite real-even data for the listed axioms, not a completed
  physical spectral triple.
- **DERIVED corrected theorem:** if `J` is additionally factor-preserving,
  then reality puts equivariant `D` in
  `(R(C[G]) intersect L(C[G])) tensor M_m
   =Z(C[G]) tensor M_m`; canonical right inner one-forms then vanish.
  This exactly scopes the previous primal/primal--dual split obstruction.
- **DERIVED orbifold complex:** the oriented icosahedral cochains have
  f-vector `(12,30,20)`, exact ranks `(11,19)`, Betti numbers `(1,0,1)`,
  and nonzero self-adjoint form-odd `D=d+d*`.
- **DERIVED stabilizers:** A5 stabilizers `(C5,C2,C3)` lift to the exact 2I
  stabilizers `(C10,C4,C6)` on vertices, edges, faces; their odd parts are
  `(C5,C1,C3)`.
- **DERIVED isotypic rows:** in A5 order `(1,3,3',4,5)`, the three layers
  have multiplicities
  `(1,1,1,0,1)`, `(0,2,2,2,2)`, `(1,1,1,2,1)`.
  The central involution acts trivially, so the module is non-regular and
  contains no spinorial/quaternionic sector.
- **DERIVED scalar-arena gate:** `A_cell=C^62` with coefficient conjugation
  has signs `(+,+,+)`, passes order zero, fails first order, and has exact
  one-form dimension 240.  Global scalars pass first order but fluctuate
  trivially.
- **SCOPE:** the full stabilizer algebras are
  `C[C10],C[C4],C[C6]`; they are not faithfully represented on scalar fixed
  cells.  Twisted/projective fibers, their algebra/opposite action, and
  stabilizer-changing Dirac maps remain **STRUCTURAL/OPEN**.
- **DERIVED/PATTERN split:** the odd stabilizer parts
  `{C1,C3,C5}` exactly equal the segregation theorem's escape list.
  Interpreting this as physical matter localization is **PATTERN/OPEN**.
- The orbifold physics gate does not open, so `Y`, Route C, multiplets,
  generations, Yukawa blocks, and frozen mass comparisons are not activated.
- Final registered suite:
  `PYTHONPATH=/tmp/science-python-deps python3 reproducible/run_all.py`
  completed in `297.7 s` with
  `Result: 51/51 scripts completed successfully.`  The new standalone
  verifier completed `27/27` exact checks.
- No PDF build was attempted.

## Session 2026-07-27, second session: primal--dual real-structure gate

- Re-established the untouched registered baseline before project writes:
  `PYTHONPATH=/tmp/science-python-deps python3 reproducible/run_all.py`
  completed in `252.7 s` with
  `Result: 49/49 scripts completed successfully.`
- Added `primal_dual_triple.md` and
  `reproducible/verify_primal_dual_triple.py`; registered the verifier as
  script 50.
- **DERIVED:** the oriented cellular dual has f-vector
  `(600,1200,720,120)` and incidence convention
  `q_j=d_(2-j)^T`.  All primal/dual `d^2=0` and transpose identities hold
  exactly over `Z`.
- **STRENGTHENED TO EXACT:** GF2 lower ranks plus chain-complex upper bounds
  prove primal ranks `(119,601,599)`, dual ranks `(599,601,119)`, and Betti
  numbers `(1,0,0,1)` on both sides.
- **DERIVED:** the dual `2I` action is free with orbit counts `(5,10,6,1)`;
  hence `C_dual=22 Reg` and the doubled arena is
  `H=C[2I] tensor C^44`, dimension 5280.
- **DERIVED:** the four cellular-star variants have sign tables
  `(+,+,-)`, `(-,+,-)`, `(+,-,-)`, and `(-,-,-)` for
  `(J^2,JD,Jgamma)`.  The first is KO6 and the second KO2.  Antiunitary
  phases change none of the signs.
- **DERIVED gate failure:** every pure-star variant preserves right
  multiplication, so the noncommutative canonical right algebra fails order
  zero and the exhibited first-order double commutator.  This includes both
  simultaneous Galois block choices and the KO6 variant.
- **DERIVED complementary failure:** composing with orbitwise inversion
  sends right to left and makes order zero and first order exact, but all
  four variants have neither `JD=DJ` nor `JD=-DJ`.  Residual nonzero counts
  are `(22000,40760)` or the exchanged pair.
- **DERIVED candidate census:** exact coefficient-span ranks
  `(1,4,9,12,12,12,12,4,9)` give
  `dim_C Omega_D^1=1191`.  The star-stable self-adjoint real dimension is
  1191; all directions are degree-off-diagonal and trace-unimodular.  The
  full unimodular gauge Lie algebra has real dimension 119.
- **SCOPE:** these are candidate field dimensions, not physical gauge or
  Yukawa fields.  No enumerated variant constructs a real spectral triple.
  Therefore `Y`, Route C, the `M15/M16` census, generations, Yukawa/mass
  blocks, and the frozen mass-exponent comparison are not activated.
- **OPEN:** a geometrically derived 44-orbit mixing unitary beyond the
  enumerated star/orientation/inversion/Galois family.
- Final registered suite:
  `PYTHONPATH=/tmp/science-python-deps python3 reproducible/run_all.py`
  completed in `296.9 s` with
  `Result: 50/50 scripts completed successfully.`  The new standalone
  verifier completed `64/64` checks.
- No PDF build was attempted.

## Session 2026-07-27: free-cell convolution and inner fluctuations

- Confirmed that the interrupted earlier attempt had created neither
  `inner_fluctuation_dichotomy.md` nor
  `reproducible/verify_inner_fluctuations.py`.  The workspace has no Git
  metadata, so no version-control cleanliness claim is made.
- Rechecked `/tmp/science-python-deps`: `numpy 2.5.1`, `scipy 1.18.0`, and
  `sympy 1.14.0` were present.  Before any project write, the untouched
  registered suite completed in `421.8 s` with
  `Result: 48/48 scripts completed successfully.`
- Added `inner_fluctuation_dichotomy.md` and
  `reproducible/verify_inner_fluctuations.py`; registered the verifier as
  script 49.
- **DERIVED:** left `2I` acts freely on all four oriented-cell layers.  The
  exact orbit counts are `(1,6,10,5)`, and 22 explicit simplex
  representatives give a signed-basis isomorphism
  `C=C[2I] tensor C^22`.
- **DERIVED:** every equivariant operator has unique form
  `sum R_(w_alpha,beta) tensor E_(alpha,beta)`.  The verifier extracts `D`
  as 124 integer group coefficients in 112 blocks and exactly reconstructs
  all 14,880 nonzero signed incidences without diagonalization.
- **DERIVED:** `C[2I]` has Wedderburn type
  `M1+M2^2+M3^2+M4^2+M5+M6`.  Exact FS indicators are
  `(+,-,+,-,+,-,+,-,+)`.  The SM-type corner choices are
  `(rho0,rho1,rho8)` and its single simultaneous Galois flip
  `(rho0,rho7,rho2)`; a concrete quaternionic matrix form and full unital
  complement allocation remain structural choices.
- **DERIVED negative:** `[D,L_g]=0` for the full 120-element spanning set,
  so all inner one-forms for the canonically left-placed group algebra
  vanish.  This placement produces no gauge or Yukawa fluctuation.
- **REFUTED predicted dichotomy:** an exact right multiplication `R_s` is
  fully `2I`-equivariant, acts on the group-algebra factor rather than only
  `C^22`, and satisfies `[D,R_s]!=0`.  Right convolution is a third horn.
- **REFUTED inverse-problem premise:** for explicit orbitwise inversion
  `J`, both left and right group algebras satisfy order zero and first order,
  but are incomparable; any algebra containing both fails order zero.
  Hence there is no unique greatest algebra and no honest single Wedderburn
  answer to the unqualified “the maximal algebra” question.  Each specified
  left/right ansatz has the Wedderburn type above; left fluctuations vanish
  and right fluctuations do not.
- **DERIVED inverse ansatz solve:** the 22-label coefficient-support graph of
  `D` is connected.  In the explicit equivariant ansatz
  `R(C[2I])` times a coordinate-diagonal multiplicity algebra, first order
  therefore reduces the latter to scalars, so the maximal algebra in that
  ansatz is exactly `R(C[2I])`.  Classification of arbitrary non-diagonal
  subalgebras of the 58,080-dimensional equivariant ambient remains open.
- **DERIVED real-structure negative:** orbitwise inversion has `J^2=+1` and
  commutes with form parity, but has neither `JD=DJ` nor `JD=-DJ`
  (`11000` and `20380` nonzeros in the two residuals).  Coefficient
  conjugation has signs `(+,+,+)` but fails order zero for a noncommutative
  left algebra.  A primal Hodge-star endomorphism does not exist; Galois
  doubling has conditional KO6 signs only on the different doubled-node
  arena, where the SM representation/order conditions remain undefined.
- **SCOPE:** the right-convolution result reopens the algebraic order-one
  screen but does not construct a real spectral triple.  A valid derived
  `J`, matter representation, `Y`, anomaly forcing, color orientation,
  generations, gauge fields, and Yukawa data remain **OPEN**.
- Final registered suite:
  `PYTHONPATH=/tmp/science-python-deps python3 reproducible/run_all.py`
  completed in `256.6 s` with
  `Result: 49/49 scripts completed successfully.`  The new standalone
  verifier completed `45/45` exact finite checks.
- No PDF build was attempted.

## Session 2026-07-24: prediction-provenance honesty correction

- Confirmed the untouched registered baseline at `47/47 PASS` in `97.8 s`
  with `PYTHONPATH=/tmp/science-python-deps`.
- Corrected a material error in our own 2026-07-22 record: JUNO's first
  result (arXiv:2511.14593, 2025-11-18) and NuFIT 6.0 (arXiv:2410.05380,
  2024-10-07) predate the February 2026 neutrino formulae.  The earlier
  “post-prediction data” / “prediction confronting new data” framing was
  false.  These comparisons are retrodictions and consistency checks.
- Added `prediction_provenance_ledger.md`, separating internal formula status
  from predictive provenance and auditing the paper's coupling, mass,
  CKM/PMNS, Higgs, and cosmological summary claims conservatively.
- Reclassified every measured oscillation parameter as RETRODICTION,
  including `sin^2(theta_12)`, `sin^2(theta_13)`, `sin^2(theta_23)`,
  `delta_CP`, and both mass splittings.  The `1/45` correction and Variant-I
  scope remain PATTERN and are explicitly data-driven.
- Restricted the genuine blind neutrino set to `m1=0`, strict normal
  ordering, `sum(m_nu)`, `m_beta`, and phase-specific `m_betabeta`; recorded
  their bounds, assumptions, and falsification windows.
- Added a documented illustrative golden-expression trials census.  It is
  neither a p-value nor proof of innocence; the binding conclusion is that
  consistency was achieved with known values and evidential weight rests on
  blind claims and internal rigidity.
- Corrected `one_integer_paper.tex`, `one_integer_supplementary.tex`, and
  `juno_2026_comparison.md`; no PDF build was attempted.
- Added `reproducible/verify_prediction_provenance.py` and registered it as
  test 48.  It pins formula/publication dates and sources, asserts all
  measured neutrino quantities are RETRODICTION, and checks that every
  neutrino BLIND item is described as unmeasured rather than merely below a
  bound.
- Final registered suite completed `48/48 PASS` in `109.5 s` with
  `PYTHONPATH=/tmp/science-python-deps`.  No PDF build was attempted.

## Session 2026-07-24, evening: C3 dynamical-selection audit

- Read the full 2026-07-22/23/24 repair ledger and treated the hostile
  adversarial audit as binding: the object studied is a
  residual-equivariant Krajewski-legal **candidate**, not a constructed
  spectral-triple Dirac.
- Dependencies were intact at `numpy 2.5.1`, `scipy 1.18.0`, and
  `sympy 1.14.0`.  The prescribed untouched baseline completed
  `46/46 PASS` in `64.9 s`.
- Added `c3_dynamical_selection.md` and
  `reproducible/verify_c3_dynamical_selection.py`; registered the verifier as
  test 47.
- **STEP 0 RECONFIRMED:** the unrestricted residual-C3 odd arena is 148 real
  dimensional; the complete legal census is `2+16+0+48=66` complex, hence
  `d0=132` real; the generic gauge-and-scale quotient is 122 dimensional.
- **DERIVED-CONSTRAINT:** for the least-assumption direct sum
  `D_tot=D_g direct-sum D_m`, the exact moments are
  `2670`, `14880+q_m`, and `55920+h_m`, with
  `q_m=Tr(D_m^2)` and `h_m=Tr(D_m^4)/2`.  Parameterizing the matter
  contributions is bookkeeping and cuts no moduli; the residual generic
  dimension remains 122.
- **EXACT CHOICE result:** demanding the old total quadratic moment remain
  exactly 14880 forces `Tr(D_m^2)=0` and hence `D_m=0`; the nonzero
  projective moduli is empty and the quartic condition is redundant.  This
  is not promoted to a closure theorem because unchanged enlarged totals
  are an extra choice, not an already-derived identity.
- **DERIVED scope negative:** the verified `Tr(Box^3)=N^2` bootstrap depends
  on the 120-vertex Hopf-fiber adjacency and face/triangle census.  The
  30-dimensional candidate has no derived analogue, so no bootstrap cut was
  imposed.
- **DERIVED scope negative:** the paired mass law
  `z_b=phi sigma(z_t)` does not define arithmetic sigma covariance of
  arbitrary Dirac coefficients.  No coefficient lattice or paired block law
  exists, so no Galois cut was imposed.
- **DERIVED conditional on the candidate model:** on the exact quartic
  critical circle,
  `T_AC(theta)=diag(r,r exp(i theta))`; the single-sheet spectrum is
  `{+r x2,-r x2,0 x26}`, with
  `Tr(D_m^2)=4r^2` and `Tr(D_m^4)/2=2r^4`.  All allowed moment cuts are
  constant on the circle.
- **PATTERN rejected:** neither registered phase
  `arctan(sqrt(5))` nor `3 arctan(sqrt(5))` is selected.  Quartic trace,
  canonical metric, and allowed moment constraints give six null
  target/diagnostic comparisons; the arithmetic Galois overlap is undefined
  and was skipped.
- **OPEN:** the defensible conditions leave 122 generic dimensions.  A
  constructed real finite triple, derived matter/geometric coupling,
  nonzero matter moments or spectral functional, arithmetic block lattice,
  masses, `Y`, hypercharge, and anomaly forcing remain open.
- Final registered suite:
  `PYTHONPATH=/tmp/science-python-deps python3 reproducible/run_all.py`
  completed in `63.1 s` with
  `Result: 47/47 scripts completed successfully.`
- No PDF build was attempted.

## Session 2026-07-24: paper and supplementary consolidation

- Read `adversarial_audit_report.md` first and treated its weakened/broken
  findings as binding.  No new theory claim was introduced.
- Updated `one_integer_paper.tex` with a binding
  `Audited results, scope, and boundaries` section, including:
  exact scalar/vector sampling theorems; the conditional canonical-metric
  wording for compact color; conditional KO6 Galois doubling and the negative
  color-conjugation result; the central-parity segregation/escape
  classification; the C3 candidate (not a constructed Dirac) and open
  122-dimensional moduli; the three-certificate matter boundary theorem;
  Kähler--Dirac carrier and invariant-spectrum results; the ordered Fibonacci
  tower; calibrated static/tower/warped dimension negatives; and audited JUNO
  scoping.
- Added the explicitly positive `Scope and boundaries` presentation: knowing
  that canonical nodes, edges, and smooth preprojective fibers do not yield SM
  matter is a derived boundary theorem.  Hypercharge is forced by anomaly
  factorization only if an appropriate module exists.
- Downgraded the gauge prefactors to **PATTERN** with the trace-index no-go
  (`5:3:3` for a real SM generation, not `8:5:2`); made the physical
  `4 a1 phi^4 -> 1/alpha0` step **STRUCTURAL** with the missing normalized
  cochain-isometry/charge axiom; flagged static `d_ST=4` as **OPEN**; and made
  the `1/45` correction and Variant-I scope **PATTERN** (empirically favored).
- Synchronized `one_integer_supplementary.tex` with the same master ledger,
  corrected the stale constructed-`D_F`/zero-Yukawa-parameter claim, corrected
  the Berry mean from exact to numerical, and added a complete 46-script
  reproducibility appendix.
- Added `consolidation_summary.md` with the section-by-section before/after
  ledger.  No PDF build was attempted.
- Baseline before manuscript edits:
  `PYTHONPATH=/tmp/science-python-deps python3 reproducible/run_all.py`
  completed with `Result: 46/46 scripts completed successfully.` in 250.0 s.
- Final post-edit suite:
  `PYTHONPATH=/tmp/science-python-deps python3 reproducible/run_all.py`
  completed in `260.5 s` with
  `Result: 46/46 scripts completed successfully.`

## Session 2026-07-24: warped inflation/AdS tower

- Dependencies were intact at `numpy 2.5.1`, `scipy 1.18.0`, and
  `sympy 1.14.0`; the untouched baseline completed `45/45 PASS` in `247.5 s`.
- Froze `N=8,16,24`, proper spacing `ell=log(phi)`, Dirichlet conditions at
  both finite ends as primary, and Neumann at both ends as sensitivity.
  Plateau criteria are unchanged from `verify_holographic_dimension.py`.
- Added `warped_spacetime.md` and
  `reproducible/verify_warped_spacetime.py`; registered it as test 46.
- **DERIVED:** the full warped spectrum is exactly the union of the
  `T_lambda` Jacobi spectra with the 52 invariant spatial levels and their
  multiplicities. A directly assembled 5-cell/N=4 control agrees with the
  union at zero reported error; all production multiplicities close at
  `2640 N` without rediagonalizing the 2640-complex.
- **DERIVED negative:** w1 `phi^(-2n)` and w2 `2^(-2n)` have no qualifying
  4D counting or heat plateau at any registered N. The result is identical
  for Dirichlet and Neumann endpoints. The w3 no-warp control is also
  negative. The single floor retains its width-`0.795`, `d_N=3.0688` 3D
  counting plateau and no 4D plateau.
- **VERDICT:** `NO_N_STABLE_4D_PLATEAU_WARPED_NEGATIVE`. This closes the
  stated finite warped-Jacobi tower route under the frozen protocol, not
  every possible continuum limit or altered radial operator.
- Exact R4 and H4 heat kernels are implemented and checked, including the H4
  `9/4` curvature gap with radius `log(phi)`. Because no 4D window appeared,
  geometry-shape fitting was not activated and no H4/Euclidean-AdS4 claim is
  made.
- **STRUCTURAL boundary:** conformal rescaling exactly recovers all 52
  600-cell `D3^2` levels and multiplicities on every floor. **INDEX OPEN:**
  the old vertex index `-4` is not an index of this warped `D4^2`.
- The standalone warped verifier completed in `27.7 s`; the post-change
  registered suite completed `46/46 PASS` in `246.8 s`. No PDF was built.

## Session 2026-07-24: inflation tower as dynamical fourth direction

- The prescribed untouched baseline completed `44/44 PASS` with
  `PYTHONPATH=/tmp/science-python-deps`; dependencies were intact at
  `numpy 2.5.1`, `scipy 1.18.0`, and `sympy 1.14.0`.
- Added `tower_spacetime.md` and
  `reproducible/verify_tower_spacetime.py`; registered the verifier as test
  45.  It imports the frozen plateau machinery unchanged from
  `verify_holographic_dimension.py`.
- **DERIVED:** form parity anticommutes with `D3`, and an exact symbolic
  certificate proves
  `D4^2=D3^2 tensor I + I tensor D_tower^2`.  Product spectra are computed
  only as eigenvalue pair sums, with multiplicity closure `2640 N`.
- **SCOPE FIX:** the Fibonacci PF ratio `phi` and McKay PF ratio `2` are
  derived tower data, but the earlier Bratteli audit explicitly left a Dirac
  coefficient sequence open.  Treating `phi^-n` or `2^-n` as Jacobi hops is
  therefore the mission-registered **STRUCTURAL** scale dictionary, not a
  newly derived AF Dirac.
- **DERIVED negative:** neither w1 (`phi^-n`) nor w2 (`2^-n`) has a
  qualifying 4D counting or heat plateau at `N=8,12,16`.  The single floor
  retains its width-`0.795`, `d_N=3.0688` counting plateau and no 4D
  plateau.
- **CONTROL negative:** the uniform path has both 4D product plateaus only
  at `N=8` (`d_N=3.7132`, `d_s=4.3352`); both disappear at `N=12,16`.
  Hence it is not an `N`-stable positive product control.
- **VERDICT:** `NO_N_STABLE_4D_PLATEAU_FINITE_SIZE_NEGATIVE`.  The registered
  derived-weight decision rule fails.  Geometry-type fitting was not
  activated, and no flat-R4, hyperbolic-H4, or Euclidean-AdS4 claim is made.
- **INDEX OPEN:** the vertex `Box` index `-4` belongs to a different operator
  hierarchy and gains no derived interpretation as an index of `D4` or a
  boundary index.  It remains unexplained as a spacetime dimension.
- **ANTI-NUMEROLOGY:** the frozen `{5,6,25,35}` exponent search found zero
  hits in 134 distinct positive tower-eigenvalue ratios.  w1 hop separations
  5 and 6 are tautological members of all allowed integer separations, not
  distinguished features.
- **Verification:** the registered suite completed `45/45 PASS` in `215.9 s`
  with `PYTHONPATH=/tmp/science-python-deps`.
- No PDF build was attempted.

## Session 2026-07-24: calibrated holographic-dimension controls

- The prescribed pre-change suite completed `43/43 PASS` with
  `PYTHONPATH=/tmp/science-python-deps`; dependency versions were
  `numpy 2.5.1`, `scipy 1.18.0`, and `sympy 1.14.0`.
- Added `holographic_dimension.md` and
  `reproducible/verify_holographic_dimension.py`; registered the verifier.
  Plateau thresholds and the complete control roster are frozen at the top
  of the verifier before spectral construction.
- **DERIVED convention fix:** because the eigenvalues are those of `D^2`,
  both estimators include the same factor two:
  `d_N=2 d log N/d log Lambda` and
  `d_s=-2 d log Tr exp(-tD^2)/d log t`.
- **DERIVED negative:** the earlier selected `d_N=3.9951` fit is not a
  plateau under the frozen full-curve rule.  The 600-cell has no accepted 4D
  counting or heat plateau.  Its sole accepted weighted interval is counting
  3D: width `0.795` decade, `d_N=3.0688`, log-RMSE `0.0659`, and local
  standard deviation `0.2303`.
- **DERIVED controls:** identical estimators were applied to the 5-, 16-, and
  triangulated 24-cell boundaries, periodic `4^3 T^3`, and genuinely 4D
  periodic `3^4 T^4`.  Fourier blocks keep the latter tractable.  Kernel
  dimensions close at 2 for the S3 controls, 8 for T3, and 16 for T4.
- **ARTIFACT / INCONCLUSIVE-CALIBRATION:** no S3 control manufactures a 4D
  plateau, but the small T4 control also fails to show 4D in both estimators.
  The old 4D window claim is therefore rejected; the present small controls
  are not a positive plateau calibration.
- **DERIVED degeneracy specificity:** stripping all shell degeneracies removes
  even the accepted 600-cell 3D counting interval and produces no 4D
  interval.  Since the weighted 4D anomaly itself fails, “symmetry
  degeneracy mimics an extra dimension” is not supported.
- **THEORY CONFRONTATION:** the static triangulated complex supplies no
  spectral derivation of `d_ST=4`.  The registered vertex-index argument is
  retained but flagged in precise tension; a fourth direction must be
  dynamical (RG scale, inflation/Bratteli tower, or time) if the theory is to
  remain four-dimensional.
- **OPEN:** larger-control finite-size sequences and a preregistered scaling
  extrapolation.  No PDF build was attempted.
- **Verification:** the registered suite completed `44/44 PASS` in `182.2 s`
  with `PYTHONPATH=/tmp/science-python-deps`.

## Session 2026-07-24: invariant Kähler--Dirac spectrum

- Froze the anti-numerology registry in
  `reproducible/verify_invariant_spectrum.py` before spectral computation.
  The untouched suite first completed `42/42 PASS`.
- Added `invariant_spectrum.md` and registered the new verifier.
- **DERIVED:** the carrier-free `2I`-isotypic spectrum has 52 distinct
  positive `D^2` levels, total dimension 2640, and kernel dimension two.
  A central class sum labels irreps without choosing multiplicity carriers.
- **DERIVED:** the Witten index is zero, nonzero even/odd spectra pair, and
  every quaternionic-sector degeneracy is even.
- **DERIVED:** the first energy is `phi^-4`; its SUSY multiplet is
  `rho1:C2(4)+C3(4)`, so `D=+/-phi^-2` with multiplicity four per sign.
- **DERIVED:** this same operator gives
  `Tr(1)=2640`, `Tr(D^2)=14880`, and
  `(1/2)Tr(D^4)=55920`, hence reduced triple `(11,62,233)`.
- **DERIVED negative:** across all 1326 distinct-level ratios, no registered
  mass exponent occurs.  Ten exact golden-power ratios occur (eight
  `phi^2`, one `phi^4`, one `phi^6`) and are reported with the full
  look-elsewhere count.
- **DERIVED negative:** the full gap is `phi^-4`, not
  `1/(2 phi^2)`; their ratio is exactly `3-sqrt(5)`.
- Zeta values at `s=3/4,1,2` and theta values at `t=5,phi` produce no
  registered identity.
- **INCONCLUSIVE:** the holographic protocol has a passing 4D counting
  window but a passing 3D heat window; it does not establish dimension flow.
- **PATTERN rejected:** finite shell counts do not support an asymptotic
  Cardy law, eta quotient, or McKay--Thompson identification.  Ihara poles
  share the vertex block only tautologically.
- Numerical unresolved roots are honestly limited to symmetric-double
  precision (about 13 reliable digits); exact minimal polynomials and
  50+-digit refinement remain **OPEN** rather than being overstated.
- No PDF build was attempted.

## Session 2026-07-23: C3 Dirac-selection audit

- Added `dirac_selection.md` and
  `reproducible/verify_dirac_selection.py`; registered and documented the
  verifier.
- Prescribed dependencies were intact.  The untouched baseline completed
  `40/40 PASS`.
- **DERIVED:** the ten `C3` subgroups of `2I=SL(2,5)` form one conjugacy
  orbit; the normalizer of a concrete `C3` has order 12.
- **DERIVED:** for the full witness
  `H+=(2,2)+2(3bar,2)`, `H-=(2,1bar)+2(3bar,2)`, the legal upper-block
  dimensions are exactly `2+16+0+48=66` complex.  Thus the self-adjoint
  first-order space is `d0=132` real, a codimension-16 subspace of the
  148-dimensional equivariant odd arena.
- **DERIVED conditional:** on the established KO6 double, `JD=DJ` determines
  the second-sheet Dirac from the first and leaves 132 real parameters.
- **DERIVED scope:** outer/Galois sheet covariance is already contained in
  that reality equation.  No arithmetic sigma action on arbitrary complex
  Dirac coefficients is derived.  A bare `sigma(D)=phi D` law forces `D=0`.
- **DERIVED:** the generic bimodule-gauge orbit has dimension 9 and the
  gauge-plus-scale moduli has generic dimension `132-9-1=122`.
- **DERIVED negative:** positive quadratic/quartic spectral actions select
  only `D=0`; the symmetry-breaking quartic has an exact gauge-inequivalent
  critical `S1`.  Polynomial criticality therefore does not yield a finite
  selection without a derived polynomial and further constraints.
- **DERIVED:** the single sheet restricts as `10 Reg(C3)`; the KO6 double is
  `20 Reg(C3)`, with `10 Reg(C3)` in each doubled chirality.
- **PATTERN, not derived:** a cyclic three-sector generation reading requires
  a noncanonical matching of the three 10-dimensional character
  multiplicity spaces and is not an SM `M15/M16` generation.
- **DERIVED negative on the generic stratum:** the surviving commutant is
  only the common scalar `u(1)`; anomaly forcing makes its common charge
  trivial, not the SM hypercharge tuple.
- **OPEN:** the selection principle remains open with residual generic
  moduli dimension 122; no physical mass/Yukawa data is predicted.
- **Verification:** the new verifier passes `20/20` internal checks and the
  registered suite completes `41/41 PASS` with the prescribed dependency
  path.
- PDF build intentionally not attempted.

## Session 2026-07-23: central-parity segregation theorem

- Added `segregation_theorem.md` and
  `reproducible/verify_segregation_theorem.py`; registered the verifier.
- Prescribed temporary dependencies were reinstalled after `/tmp` was wiped.
  The untouched baseline completed `39/39 PASS` before changes.
- **DERIVED:** for `z=-1`, every real- or complex-linear `2I`-equivariant
  operator commutes with `rho(z)`.  Hence every equivariant gamma-odd
  operator is zero.  The proof includes real and quaternionic irreducibles.
- **DERIVED, scoped:** the lemma reproduces the zero equivariant odd sector
  on the node module, source-graded edge/Hom module, preprojective regular
  fiber, and consecutive Bratteli floors.  The older endpoint-first-order and
  missing-support results remain independent stronger obstructions.
- **DERIVED:** exact enumeration in `SL(2,5)` proves that `-I` is the unique
  involution.  Odd-order subgroup types are exactly `C1,C3,C5`; there are ten
  `C3` and six `C5` subgroups.  An exact normalizer computation excludes
  order 15.
- **DERIVED:** a nonzero parity-odd Dirac has maximal residual symmetry at
  most `C3` or `C5`.
- **DERIVED:** on `W`, parity-odd equivariant operator dimensions are
  `88` for `C5`, `148` for `C3`, and `448` for `C1` (complex dimension for
  all odd operators, equivalently real dimension for self-adjoint odd
  operators).
- **DERIVED positive:** the canonical seed restrictions to `C3` admit the
  explicit Krajewski witness
  `H+=(2,2)+2(3bar,2)`, `H-=(2,1bar)+2(3bar,2)`, with a nonzero legal odd
  first-order block.
- **DERIVED negative:** the corresponding exact `C5` grading-wise reblocking
  has zero solutions despite its 88-dimensional odd operator space.  The
  obstruction is cyclic character content.
- **OPEN:** the `C3` result does not select `D`; the loss of uniqueness is
  quantified by the 148-dimensional self-adjoint odd space.  Its
  gamma-preserving commutant admits a 10-real-dimensional self-adjoint family
  before `D`, so no canonical `Y` is selected and Route C cannot yet start.
- **PATTERN:** `C3` matches `N_gen=3` and `C5` matches `a1=5`; no principle
  yet derives which subgroup is physical or why.
- **Verification:** the new verifier passes `17/17` internal checks and the
  registered suite completes `40/40 PASS` with the prescribed dependency
  path.
- No PDF build was attempted.

Working log for theory-side repairs before updating `one_integer_paper.tex`.

## 2026-03-28

- Confirmed the Kirchhoff bug: `tau(icosahedron) = 5,184,000`, so the
  matrix-tree justification for `L(3)L(5)L(3') = 120` must be removed.  The
  algebraic identity itself remains valid.

- Added adversarial audit:
  `reproducible/referee_check.py`

- Strengthened bounded-search uniqueness from the paper's `C1 + C2` to:
  `C1 + Q + C2`
  via:
  `reproducible/verify_uniqueness_quantized.py`

- Found a still stronger bounded-search criterion:
  `C1 + S`
  where
  `S: dz in {0} U {+/-phi^r} U {+/-2phi^r}`
  via:
  `reproducible/verify_uniqueness_simple_edges.py`

- Replaced the simple-edge condition by a more primitive constructive rule:
  minimal `L^1` edge lifts of
  `5 da + 6 db = Delta n`
  on the neutral branch ordering.
  Files:
  `reproducible/verify_minimal_edge_lifts.py`
  `minimal_edge_lift_theorem.md`

- Verified that the constructive edge-lift rule reconstructs the full neutral
  assignment and makes `C1`, `Q`, and `C2` emergent.

- Added all three new uniqueness verifiers to:
  `reproducible/run_all.py`
  and documented them in:
  `reproducible/README.md`

- Current suite status:
  `python reproducible\\run_all.py`
  -> `17/17 PASS`

- Branch-label issue isolated:
  the arithmetic uniqueness scripts use neutral endpoints `brA/brB`, while the
  physical McKay graph distinguishes the length-2 BLACK endpoint from the
  short-leg WHITE endpoint.

- Added branch-identification verifier:
  `reproducible/verify_branch_identification.py`
  Goal: resolve `top/bottom` from McKay geometry + chirality + generation
  formulas, not by manual relabeling.

- Branch identification now closed:
  the affine-E8 leg geometry plus bipartite McKay chirality and the internal
  generation formulas force
  `length-2 BLACK endpoint -> top (n=26)`
  and
  `short-leg WHITE endpoint -> bottom (n=19)`.
  This is verified in:
  `reproducible/verify_branch_identification.py`

- Added a bounded-search-free constructive chain for the full physical
  assignment:
  `reproducible/verify_global_uniqueness_constructive.py`
  It reconstructs
  `e,u,d,s,mu,c,tau,t`
  from main-chain and branch minimal lifts, then forces
  `b`
  from the prime-sector Galois relation
  `z_b = phi * sigma(z_t)`.

- Strengthened the spectral-action side on the exact discrete layer:
  `reproducible/verify_spectral_action.py`
  now verifies
  `c0 = 2640`,
  `c1 = 14880`,
  `c2 = 55920 = (1/2) Tr(D^4)`,
  the reduced triple
  `(A0, A1, A2) = (11, 62, 233)`,
  and the exact Diophantine identity
  `2 A1^2 + 1 = 3 A0 A2`.

- Added working note:
  `spectral_action_discrete_theorem.md`
  to separate the exact discrete theorem from the still-open continuum
  identification with the full Chamseddine-Connes Standard Model action.

- Added an exact verifier for the static scalar-response sector:
  `reproducible/verify_discrete_scalar_response.py`
  It proves numerically on the full 600-cell that
  `B^+ d_0 = d_0 Delta_0^+`,
  every point source gives `h = d_0 Phi`,
  `d_1 h = 0`,
  and `d_0^+ h = Phi`.

- Added working note:
  `discrete_scalar_response_theorem.md`
  to separate the exact discrete theorem from the stronger continuum PPN
  interpretation.

- Updated `one_integer_paper.tex` discussion/summary wording to match the
  repaired theory state:
  - constructive `(a,b)` uniqueness is now marked as derived, not merely
    bounded-search-backed;
  - the spectral-action discrete coefficients
    `(2640, 14880, 55920)` and reduced triple `(11, 62, 233)` are stated as
    exact;
  - the gauge prefactors `(8/15, 1/3, 2/15)` are downgraded to an open
    structural normalization / continuum-identification issue;
  - the gravity wording now explicitly distinguishes exact
    `gamma_disc = 1` from any continuum PPN claim.

- Attempted local PDF rebuild of `one_integer_paper.tex`, but verification is
  currently blocked by the local MiKTeX installation state
  (`fresh installation` / `Access is denied` under the user profile setup
  directory).

## Session 2: 2026-03-28 (wording hardening pass)

Changes made:

### 1. Gauge prefactors (8/15, 1/3, 2/15) - DEFINITIVELY DOWNGRADED
- Main paper eq (lagrangian_explicit): underbrace now says "Yang--Mills (pattern)"
- Main paper line after the equation: explicitly separates DERIVED (1+3+8 skeleton)
  from PATTERN (specific fractions), with explanation that the continuum trace
  computation is the missing step.
- Summary table already had Pattern; no change needed there.

### 2. Stale wording cleanup
- Supplement master table: `gamma_PPN = 1` -> `gamma_disc = 1 (discrete scalar
  response)` with status `Derived^ddagger` and new footnote explaining the
  discrete-only scope.
- Supplement uniqueness proof: replaced old bounded-search proof with reference
  to the constructive derivation as primary, bounded search as cross-check.
  Added analytic bound argument (leading term -19t^2/25 ensures |N|>19 for |t|>5).
- Supplement footnote p: updated from "exhaustive search" to "constructive
  derivation; exhaustive search confirms as cross-check."

### 3. Alpha derivation audit
- Main paper: added "Geometric identification" paragraph explaining why phi^4
  enters (spectral diffusion factor from Cayley graph), with explicit caveat
  that the identification 1/alpha_0 = 4*a1*phi^4 is STRUCTURAL (motivated by
  Kaluza-Klein), not a theorem from the spectral action alone.
- Main paper derivation status: 2pi = DERIVED, 4*a1*phi^4 = STRUCTURAL,
  1 = NORMALIZATION. Three distinct categories, not a blanket "all derived."
- Main paper summary table: alpha separated from alpha_s and sin^2(tW), with
  status STRUCTURAL (not Derived).
- Supplement master table: alpha status changed to Structural^a.

### 4. Gravity wording hardening
- Main paper: already properly hardened (gamma_disc, explicit disclaimers).
  No changes needed.
- Supplement: gamma_PPN -> gamma_disc in master table (the only remaining issue).

### 5. Supplement synchronization
- Master table: gamma_disc, alpha structural, uniqueness footnote updated.
- Uniqueness proof: constructive primary, bounded search secondary.
- ddagger footnote added for discrete scalar response scope.

### What was strengthened
- Nothing. This session was purely about honesty hardening.

### What was downgraded
- Gauge prefactors: from implicit "arise from" to explicit "pattern"
- Alpha: from "Derived" to "Structural" (the phi^4 identification)
- gamma in supplement: from PPN to disc

### What remains open
- Continuum gauge normalization (8/15, 1/3, 2/15): definitively a pattern.
- Alpha phi^4 factor: structural identification, not theorem.
- Continuum PPN interpretation: not claimed.
- TeX rebuild still blocked by local MiKTeX installation.

### Verification
- `python reproducible/run_all.py` -> `20/20 PASS`

## Current status

- Fixed on the theory side:
  - bounded-search uniqueness is no longer just a brute-force `C1 + C2` claim;
  - there is now a constructive local edge principle behind the neutral
    assignment;
  - the neutral branch labels are now resolved structurally to `top/bottom`.
  - the full physical assignment now has a bounded-search-free constructive
    derivation chain.
  - the spectral-action section now has an exact discrete core:
    coefficient triple `(2640, 14880, 55920)`,
    reduced triple `(11, 62, 233)`,
    and a Diophantine identity that again singles out `a1 = 5`.
  - the gravity section now has an exact discrete scalar-response core:
    point sources induce only exact edge fields, with no independent static
    coexact scalar component on the tested sector.
  - gauge prefactors definitively downgraded to pattern.
  - alpha derivation honestly separated into derived / structural / normalization.
  - supplement synchronized with main paper on all five topics.

- Still open:
  - once TeX is usable locally, rebuild `one_integer_paper.tex` and check the
    rendered discussion/summary section for overflow or caption issues.

## Session: 2026-07-22

### Thread 1: gauge prefactors

- Added `gauge_prefactor_no_go.md`.
- **DERIVED:** the gauge kernel and its `1+3+3'+5` decomposition, together
  with the scalar restrictions `C=5 I`, `B=(16/5) I`, and
  `Delta_1=(41/5) I` on the 12-dimensional gauge sector.
- **PATTERN:** `(8/15, 1/3, 2/15)` remains un-derived. The present symmetry
  permits independent invariant quadratic-form normalizations on the
  `1`, `3`, and `8` sectors, and the available scalar operator restrictions
  do not fix them.
- **OPEN:** a derived Lie bracket on `3'+5`, a common matter representation,
  and a single normalized trace (including `U(1)` charge normalization).

### Thread 2: alpha structural step

- Added `alpha_phi4_missing_axiom.md`.
- **DERIVED:** the finite identity
  `Tr(Box_gauge^2)/(n_base*b1*lambda_f^2)=4*a1*phi^4` and the separate Hopf
  holonomy `2*pi`.
- **STRUCTURAL:** identification of that spectral number with physical
  `1/alpha_0`; an arbitrary rescaling of the discrete-to-physical `U(1)` field
  leaves all present spectral data fixed but changes the action coefficient.
- Minimal missing axiom: a normalized discrete-to-continuum `U(1)` cochain
  map that is an isometry for specified discrete and continuum inner products,
  together with a fixed charge unit.
- Corrected stale overclaims in `reproducible/verify_coupling_constants.py`,
  `one_integer_paper.tex`, and `logic_chain_map.md` without weakening any
  finite spectral identity.

### Thread 3: exact scalar bridge

- Added `low_mode_sampling_intertwiner_theorem.md` and
  `reproducible/verify_low_mode_sampling_intertwiner.py`.
- **DERIVED:** evaluation of `H_0+H_1+H_2` on the 12-point Hopf base is an
  injective rank-9 sampling map with exact left reconstruction. A quadratic
  polynomial in the icosahedral graph Laplacian exactly intertwines the
  round-sphere scalar Casimir values `0,2,6`, and reconstruction annihilates
  the remaining 3-dimensional alias sector.
- **OPEN:** refinement convergence, vector fields, local transport, a Lie
  bracket, and the nonabelian continuum connection.

### Reproducibility

- Registered the new verifier in `reproducible/run_all.py` and documented it
  in `reproducible/README.md`.
- Initial system-`python3` baseline: `7/27 PASS`; 20 scripts failed at import
  because NumPy was absent, not because verifier assertions failed.
- Installed NumPy/SciPy only into `/tmp/science-python-deps` for this session.
- Full suite with that temporary dependency path: `28/28 PASS`.
- PDF build intentionally not attempted, per the session rules.

## Session: 2026-07-22 — vector lift

- Added `vector_sampling_intertwiner_theorem.md` and
  `reproducible/verify_vector_sampling_intertwiner.py`.
- Registered the verifier in `reproducible/run_all.py` and documented it in
  `reproducible/README.md`.
- **DERIVED:** geodesic edge integration is the canonical de Rham target for
  continuum 1-forms; it commutes exactly with `d` on sampled scalar modes.
- **DERIVED:** gradient and curl vector harmonics for `l=1,2` give an injective
  rank-16 edge-cochain space with exact/coexact split `8+8`.
- **DERIVED:** the raw edge Hodge eigenvalues are
  `(5-sqrt(5),6,3-sqrt(5),2)`, and an explicit cubic polynomial in `Delta_1`
  exactly reproduces the continuum 1-form eigenvalues `(2,6,2,6)`.
- **DERIVED:** the alias complement has dimension 14, split as 3 exact and 11
  coexact modes.
- **DERIVED (negative):** the four discrete/continuum norm factors are unequal.
  Since the icosahedron has one edge orbit, no local diagonal invariant edge
  weight makes the full low vector sampling map an isometry.
- **STRUCTURAL:** a mode-dependent spectral metric can force isometry, but its
  four weights and its extension to the alias complement are extra choices.
- **OPEN:** a principle selecting that metric; therefore the physical `U(1)`
  normalization needed for alpha remains open and `1/alpha_0=20*phi^4`
  remains structural.
- **DERIVED (negative):** `d_0` maps the 12 gauge amplitudes with rank 11,
  kills the constant sector, resolves `3+5`, and leaves `3'` as an exact alias.
  It therefore supplies no common trace on `1+3+8` and does not derive the
  gauge prefactors, which remain pattern.
- PDF build not attempted.

## Session: 2026-07-22 — `A5` bracket classification

- Added `a5_equivariant_bracket_theorem.md` and
  `reproducible/verify_a5_equivariant_brackets.py`.
- Registered and documented the new verifier.
- **DERIVED:**
  `Lambda^2(3'+5)=2(3)+3(3')+2(4)+5`, giving three equivariant maps to
  `3'` and one to `5`.
- **DERIVED:** the complete real Jacobi variety is `b=0`, with either
  `c=d=0` or `a=c`.
- **DERIVED:** its branches are compact `su(3)`, split `sl(3,R)`, semidirect
  and direct degenerations, a two-step nilpotent algebra, and the abelian
  point.
- **DERIVED:** the compact-simple class is unique up to scale and
  `A5`-equivariant isomorphism. Its apparent continuous negative-sign family
  is only a block-rescaling orbit.
- **STRUCTURAL:** equivariance and Jacobi do not select compact simplicity
  over the other branches, so physical color selection retains an explicit
  compactness/simple-algebra criterion.
- **DERIVED:** an explicit `3'` embedding with `Tr(AB)=phi'` generates 60
  rotations and gives `ad(su(3))=3'+5` in the edge-kernel convention.
- **DERIVED:** the compact Killing form gives the same Frobenius normalization
  on `3'` and `5`, closing the internal color-block normalization freedom.
- **OPEN:** relative trace normalization between `su(3)`, `su(2)`, and the
  abelian factor; a common matter representation and `U(1)` charge unit remain
  necessary. Gauge prefactors remain **PATTERN**.
- **DERIVED:** the `3` bracket is the unique nonzero cross product up to scale;
  the `1` bracket is necessarily abelian.
- Synchronized stale `3+5`/`3'+5` wording with the edge-kernel convention in
  `reproducible/verify_spectral_action.py` and `one_integer_paper.tex`.
- PDF build not attempted.
## Session 2026-07-22: canonical edge metric and compact color

- **DERIVED:** identified the 12 gauge modes with the alternating edge mode on
  each Hopf fiber's ten-cycle (not the constant lift), giving Gram matrix
  `10 I_12` in fiber-amplitude coordinates.
- **DERIVED:** exact icosahedral sampling moments give the canonical color
  metric `20 Frob` on `3'` and `16 Frob` on `5`, with ratio `5/4`.
- **DERIVED:** ad-invariance of this actual metric is exactly `b=0` and
  `16c+20d=0`. Its intersection with the Jacobi variety contains only the
  abelian bracket, `so(3)+R^5` with central `R^5`, and the compact-simple
  `su(3)` family `(a,0,a,-4a/5)`.
- **DERIVED:** the split `sl(3,R)`, semidirect, and nilpotent branches fail
  invariance of the canonical positive metric.
- **DERIVED (conditional):** within the metric-compatible list, trivial center
  uniquely selects `su(3)` up to overall bracket scale.
- **STRUCTURAL:** the remaining color input is the axiom that the bracket is a
  center-free metric Lie algebra for the canonical edge inner product.
- Updated `reproducible/verify_a5_equivariant_brackets.py` and
  `reproducible/verify_edge_gauge_spectrum.py`; updated the bracket theorem and
  gauge-prefactor no-go notes. No PDF build was attempted.
## Session 2026-07-22: finite matter trace-index audit

- Added `matter_trace_index_no_go.md` and
  `reproducible/verify_matter_trace_indices.py`; registered and documented the
  verifier.
- **DERIVED:** the prefactor target is exactly `(T1:T2:T3)=(8:5:2)`, or
  `5T1=8T2` and `2T2=5T3`.
- **DERIVED:** one Standard Model Weyl generation has ordinary-hypercharge
  indices `(10/3,2,2)`, ratio `5:3:3`; GUT-normalized hypercharge gives
  `(2,2,2)`. Neither is the target.
- **DERIVED conditional on the compact bracket choice:** the `1+3+8` gauge
  algebra acting on itself has adjoint indices `(0,2,3)` and fails the target.
- **DERIVED (negative):** no unconditionally derived matter module has all
  three actions; the conditionally available adjoint module does not have
  index ratio `8:5:2`.
- **OPEN:** the 120-vertex regular sectors, affine-E8 McKay nodes, `(a,b)`
  fermion slots, Hopf/spectral candidates outside the adjoint construction do
  not carry a derived simultaneous `u(1)+su(2)+su(3)` action. Their trace
  indices cannot yet be formed.
- **PATTERN:** the gauge prefactors remain non-derived. The missing object is
  an explicit common matter module with a derived rational `U(1)` generator.
- No PDF build was attempted.

## Session 2026-07-22: exact inflation/Bratteli towers

- Added `bratteli_inflation.md` and
  `reproducible/verify_bratteli_tower.py`; registered and documented the
  verifier.
- Baseline with the prescribed dependency path completed `38/38 PASS` before
  changes.
- **DERIVED:** the Fibonacci fusion and rooted affine-E8 McKay towers were
  constructed exactly through level 12, with total McKay floor dimension
  `2^n` and exact bipartite support.
- **DERIVED:** Fibonacci ordered `K0` is the golden dimension group
  `(Z[phi], Z[phi] intersect R_+ ,1)` with its canonical PF state.
- **DERIVED:** after removing unreachable parity vertices and telescoping,
  rooted McKay `K0` is `lim(Z^4,M)` for the explicit determinant-4 matrix in
  the results note.  Its canonical trace range is `Z[1/2]`, with a rank-3
  infinitesimal kernel and exact sequence
  `0 -> Z^3 -> K0 -> Z[1/2] -> 0`; no unsupported splitting is claimed.
- **DERIVED pattern check:** the theory's `Z[phi]` lattice is canonically the
  Fibonacci trace lattice but cannot embed trace-preservingly in the rational
  McKay trace range.  A physical scale identification remains STRUCTURAL.
- **DERIVED:** McKay block multiplicities grow and remove the old block-size
  obstruction.
- **DERIVED all-level negative:** rooted parity prevents the scalar/color
  seed nodes and quaternionic weak seed node from coexisting on any floor.
- **DERIVED scoped negative:** after the STRUCTURAL change to a cumulative
  endpoint model, every nonzero consecutive-floor shift is gamma-odd but has
  an exact nonzero first-order double-commutator witness.
- **OPEN:** a different derived AF representation/opposite action, real
  structure, color orientation, generation blocks, and `Y`.  Route C remains
  unavailable.  The coefficient sequence `c_n` is NOT DERIVED.
- **Verification:** `verify_bratteli_tower.py` passes `25/25` internal checks;
  the registered suite completed `39/39 PASS` in 12.8 s with
  `PYTHONPATH=/tmp/science-python-deps`.
- No PDF build was attempted.

## Session 2026-07-22, session 2: matter-module construction audit

- Added `matter_module_no_go.md` and
  `reproducible/verify_matter_module.py`; registered the verifier.
- **DERIVED:** exhaustive dimension equations give 42 decompositions at
  dimension 15 and 48 at dimension 16 from block sizes `6,3,2,1`.
- **DERIVED:** imposing the SM inventory gives the expected nonabelian
  multiplicities `(1,2,1,1)` or `(1,2,1,2)`.
- **STRUCTURAL:** the defining `2I/SU(2)` and `A5/SU(3)` actions permit an
  abstract external tensor-product module, but no discrete matter space is
  identified with it and the `A5` restriction does not select `3` versus
  `bar(3)` chirality.
- **DERIVED:** the SM benchmark has indices `(10/3,2,2)`, ratio `5:3:3`, all
  four local anomaly sums zero, and even Witten parity.
- **DERIVED (negative):** the McKay exponent pairs on quark doublets are
  `(3,5),(16,11),(26,19)` and the C10 residues are
  `(3,5),(6,1),(6,9)`.  Both fail the necessary condition that a commuting
  `u(1)` generator be scalar on each weak doublet.
- **DERIVED (negative):** the `Z[phi]` unit exponent is not defined on all
  nine slots; their norms include `0`, `5`, `19`, and `-19`.
- **OPEN:** a derived chiral matter functor, color-conjugation choice,
  neutrino placement, and everywhere-defined integer commutant generator.
  Hypercharge and anomaly cancellation remain un-derived.
- **PATTERN:** the old `8:5:2` target remains a pattern; the physical target
  is the verified `5:3:3` benchmark.
- **Verification:**
  `PYTHONPATH=/tmp/science-python-deps python3 reproducible/run_all.py`
  completed with `32/32 PASS` (184.0 s), including
  `verify_matter_module.py`.
- No PDF build was attempted.

## Session 2026-07-22, session 3: derived chiral matter-functor routes

- Added `matter_functor_route_no_go.md` and
  `reproducible/verify_matter_functor.py`; registered the verifier.
- Baseline reconfirmed with the prescribed dependency path: `32/32 PASS`.
- **DERIVED (Route A):** the 30-dimensional McKay sum has bipartite grading
  dimensions `16+14`.  Therefore no invertible KO6-type real structure can
  anticommute with this grading.
- **DERIVED (Route A):** the Galois node involution preserves dimensions and
  commutes with the bipartite grading; it is not yet an antiunitary
  endomorphism with a derived Hilbert metric and opposite-algebra action.
- **DERIVED (negative, Route A):** the repository supplies McKay Dirac
  topology and `gamma_F`, but not matrices for a representation of
  `C+H+M3(C)`, its opposite action, or `J`.  Hence order-zero, first-order,
  and unimodularity equations are presently undefined rather than failed.
- **DERIVED (Route B):** the SM benchmark reduces to C10 residues
  `(1,6,2,7,6,0)` and passes the commutant gate.
- **DERIVED (negative, Route B):** six unconstrained block characters leave
  `10^6` assignments, and C10 cannot distinguish integer lifts differing by
  10; notably `u^c` and `e^c` both have residue 6.
- **PATTERN:** the displayed SM winding tuple is consistent but unselected.
- **DERIVED (Route C):** the complete M15 anomaly system factorizes as
  `18q(2q-u)(4q+u)=0`.  For nonzero primitive `q`, it fixes the SM charge
  tuple up to overall sign/scale and exchange of the two colored singlets.
- **DERIVED (Route C):** for M16 the factor is
  `18q(-n+2q-u)(-n+4q+u)`, leaving the extra-singlet charge parameter; anomaly
  cancellation alone does not force `nu^c` neutrality.
- **OPEN:** a derived real associative-algebra representation or an
  equivariant matter-section functor with canonical integral winding lift,
  plus chiral/color orientation and generation-blind block assignment.
- **Verification:**
  `PYTHONPATH=/tmp/science-python-deps python3 reproducible/run_all.py`
  completed with `33/33 PASS` (205.7 s), including
  `verify_matter_functor.py` (`17/17` internal checks).
- No PDF build was attempted.

## Session 2026-07-22, session 3 follow-up: the `16=M16` test

- Extended `matter_functor_route_no_go.md` and
  `reproducible/verify_matter_functor.py` with an exact 2I character test.
- **DERIVED:** the even half is
  `rho_1+rho_4+rho_5+rho_6+rho_8`, dimension 16, containing the existing
  slots `e,d,b,tau,mu`.  The odd half is
  `rho_2+rho_3+rho_7+rho_9`, dimension 14, containing `u,t,s,c`.
- **DERIVED:** in class order
  `(1A,2A,4A,6A,3A,10A,5A,5B,10B)`, their characters are
  `(16,16,0,1,1,1,1,1,1)` and `(14,-14,0,1,-1,1,-1,-1,1)`.
- **DERIVED:** the correct diagonal restriction uses color `rho_5=3'` and
  weak `rho_2=2`; the McKay rule gives `rho_5 tensor rho_2=rho_9`.  Hence
  `M16|2I=rho_9+2rho_5+rho_2+2rho_1`.
- **DERIVED (negative):** the even half is not diagonal M16.  Already at the
  central element, its character is 16 while the M16 character is 0.
  Therefore `16=dim(M16)` is a numerical coincidence for this action.
- **DERIVED (negative, scoped):** the odd half matches neither diagonal M15
  with one singlet removed nor gauge-12 content plus the defining spinor.
- **OPEN:** no physical interpretation of the odd half; other functors remain
  possible but must still construct color orientation and `Y`.
- **Verification:** the extended `verify_matter_functor.py` passes `26/26`
  internal checks, and the registered suite completed `33/33 PASS` in
  157.4 s with `PYTHONPATH=/tmp/science-python-deps`.
- No PDF build was attempted.

## Session 2026-07-22: JUNO 2026 confrontation and `1/45` scope

- Added `juno_2026_comparison.md` and
  `reproducible/verify_juno_comparison.py`; registered and documented the
  verifier.
- Independently checked primary sources online: JUNO arXiv:2511.14593,
  NuFIT~6.0 arXiv:2410.05380 plus its official parameter table, DESI
  arXiv:2503.14744, dynamical-DE arXiv:2507.16589, and KATRIN
  arXiv:2406.13516 / Science 388 (2025).
- **EXTERNAL discrepancy:** the paper's stale
  `Delta m^2_32=(2.453+-0.033) 10^-3 eV^2` is not the official NuFIT~6.0
  normal-ordering entry.  Official `Delta m^2_31` is
  `2.534^(+0.025)_(-0.023) 10^-3` without SK atmospheric and
  `2.513^(+0.021)_(-0.019) 10^-3` with SK.
- **DERIVED recomputation:** JUNO solar deviations are `+0.699 sigma` for
  Variants I/III and `+2.065 sigma` for Variant II.  NuFIT-with-SK atmospheric
  deviations are `-0.346 sigma` for corrected `m3` (I/II) and
  `-3.138 sigma` for bare Variant III.
- **DERIVED empirical diagnostic:** the two-input chi-squares are
  `0.6089`, `4.3859`, and `10.3346` for I, II, and III.  This ranks the
  variants but is not a model-selection significance.
- **STRUCTURAL:** treating the derived ratio as a bare spectral relation
  before a later eigenstate-local correction motivates Variant I.
- **PATTERN:** the `1/45` mass correction and its scope remain pattern because
  no finite perturbation operator determines whether the ratio is reimposed.
- **EXTERNAL:** all three sums, `58.24--58.87 meV`, pass the DESI LambdaCDM
  `64.2 meV` bound and lie about `1.06--1.08 sigma` below the model-dependent
  dynamical-DE positive preference.  KATRIN is non-discriminating.
- Updated only stale experimental comparison values and the scoping paragraph
  in `one_integer_paper.tex`; no theory formula was changed.
- **Verification:** `verify_juno_comparison.py` passes `17/17` internal checks;
  the registered suite completed `34/34 PASS` in 210.5 s with
  `PYTHONPATH=/tmp/science-python-deps`.
- No PDF build was attempted.

## Session 2026-07-22: Galois-doubled McKay triple

- Added `galois_doubling_triple.md` and
  `reproducible/verify_galois_doubling.py`; registered the verifier.
- **DERIVED:** the (A_5) outer automorphism equals
  `sqrt(5) -> -sqrt(5)` on characters, and its lift to (2I) acts as
  `2 <-> 2'`, `3 <-> 3'`, fixing `1,4,4_s,5,6`.
- **DERIVED (linchpin negative):** complex conjugation fixes the real derived
  `3'` embedding, whereas the outer twist changes it to inequivalent `3`.
  Thus Galois twisting is not finite-subgroup color `3 <-> bar3`.
- **DERIVED, conditional:** on `W + sigma(W)`, opposite sheet chirality and
  compatible twisted adjacency realize the KO6 signs
  `J^2=+1`, `JD=DJ`, `J gamma=-gamma J`; the old `16 != 14` obstruction is
  removed because `J` pairs two 30-dimensional sheets.
- **STRUCTURAL/OPEN:** arithmetic Galois conjugation is not itself
  antiunitary; isometric anti-linear intertwiners and the doubled sheet must
  still be selected by the geometry.
- **OPEN:** no derived unital action of `C + H + M3(C)` on the node space is
  available, so target order zero, first order, and the surviving `Y` cannot
  be evaluated without inserting the missing bimodule.
- **DERIVED (scoped negative):** diagonal `C^9` vertex actions satisfy order
  zero, but McKay adjacency fails first order for independent endpoint
  projectors; only the constant test is automatic.
- **DERIVED:** the bidirected quiver Hom space has dimension 240 and diagonal
  decomposition `16 rho_2 + 6 rho_3 + 16 rho_7 + 22 rho_9`; it contains no
  integer-spin color/singlet irreps and exposes no direct M15/M16 block.
- **Verification:** `verify_galois_doubling.py` passes `34/34` internal checks;
  the registered suite completed `35/35 PASS` in 197.2 s with
  `PYTHONPATH=/tmp/science-python-deps`.
- No PDF build was attempted.

## Session 2026-07-22, session with Claude: derived bimodule decisive audit

- Added `bimodule_krajewski_result.md` and
  `reproducible/verify_bimodule_krajewski.py`; registered and documented the
  verifier.
- Baseline with the prescribed dependency path completed `35/35 PASS` before
  changes.
- **DERIVED:** exhaustive generating-function enumeration gives 188,908,396
  dimension-only `5 x 5` Krajewski multiplicity matrices at dimension 30 and
  1,362,811,872,984 at dimension 60.  Dimension arithmetic is not the
  obstruction.
- **DERIVED (negative):** requiring the derived diagonal restrictions
  `1 -> rho_1`, weak `2 -> rho_2`, and color `3,bar3 -> rho_5=3'` leaves
  untwisted block support only on
  `{rho_1,rho_2,rho_4,rho_5,rho_8,rho_9}`.  The exhaustive compatible search
  for `W` has zero solutions.
- **DERIVED (negative):** allowing all outer-twisted blocks on the Galois
  sheet still omits `rho_6,rho_7`; the exhaustive search for
  `W+W^sigma` also has zero solutions.  The killing constraint is therefore
  `2I` compatibility, before McKay-Dirac legality is reached.
- **DERIVED:** `End_2I(W+W^sigma)=M2(C)^9`.  Opposite-sheet gamma reduces it
  to `C^18`; the conditional sheet-swap `J` gives `M2(R)^9`, and imposing
  both gives `C^9` as a real algebra.
- **DERIVED (negative):** no product of `M2(C)` factors contains a unital
  `M3(C)` subalgebra.  Hence neither the maximal doubled-node commutant nor
  any `J`/gamma/Dirac-selected subalgebra can contain
  `C+H+M3(C)`.  This is a fundamental node-matter no-go.
- **DERIVED:** the 240-dimensional edge/Hom decomposition instead has
  commutant `M16(C)+M6(C)+M16(C)+M22(C)`.
- **STRUCTURAL:** that edge commutant admits an abstract faithful
  `C+H+M3(C)` representation (`H` in an `M16` factor and `M3` in `M6`), but
  no derived embedding, grading, real structure, or first-order Dirac selects
  it.
- **OPEN:** the edge/Hom Krajewski construction, color orientation,
  generation-blind commutant generator, and hence hypercharge.  Route C
  anomaly factorization remains available only after these objects exist.
- **Verification:** `verify_bimodule_krajewski.py` completed `20/20 PASS`; the
  registered suite completed `36/36 PASS` in 127.3 s with
  `PYTHONPATH=/tmp/science-python-deps`.
- PDF build was intentionally not attempted.

## Session 2026-07-22: edge-space Krajewski construction

- Added `edge_matter_krajewski.md` and
  `reproducible/verify_edge_matter.py`; registered the verifier.
- Baseline with the prescribed dependency path completed `36/36 PASS` before
  changes.
- **DERIVED:** the 16 oriented affine-E8 arrows are distinct irreducible
  `2I x 2I` outer products and have total dimension 240.
- **DERIVED:** their full two-sided commutant is `C^16`.  Thus restoring the
  missed two-sided action does not select the abstract `M3` and `H` available
  in the much larger diagonal-`2I` commutant.
- **DERIVED:** `2 tensor 3'=6`, whereas `2 tensor 3=2+4s` and
  `2 tensor 4=6+2'`.  Hence only the eight arrows over `0--1`, `1--2`,
  `2--3`, and `5--8` have ambient actions on both endpoints from the supplied
  seeds, symmetric powers, and the multiplicity-one product factorization.
- **DERIVED (negative):** the real Galois pair `3,3'` still does not encode
  color `3,bar3`; complex conjugation fixes each real embedding.
- **DERIVED:** orientation grading and Hilbert--Schmidt adjoint reversal give
  `J^2=+1` and `J gamma=-gamma J` canonically.
- **DERIVED decisive no-go:** for the canonical endpoint bimodule,
  first-order-legal Dirac blocks share a source or target and therefore
  preserve orientation.  Odd blocks reverse orientation and share neither.
  Thus first order plus evenness forces `D=0`.
- **DERIVED:** length-two path multiplication has zero projection back to the
  edge space, while the preprojective moment map is quadratic and
  vertex-valued rather than a linear Dirac operator.
- **OPEN:** a separately derived smaller endpoint algebra might change the
  first-order test, but no `C+H+M3(C)` action, generation census, or
  generation-blind `Y` is selected.  Route C therefore cannot start.
- **Verification:** `verify_edge_matter.py` passes `21/21` internal checks;
  the registered suite completed `37/37 PASS` in 27.7 s with
  `PYTHONPATH=/tmp/science-python-deps`.
- No PDF build was attempted.

## Session 2026-07-22: preprojective smooth-fiber last-room audit

- Added `preprojective_matter.md` and
  `reproducible/verify_preprojective_matter.py`; registered the verifier.
- Baseline with the prescribed dependency path completed `37/37 PASS` before
  changes.
- **DERIVED:** a smooth free-orbit fiber is the regular 120-dimensional
  `2I` module, with commutant factors
  `M1+M2+M3+M4+M5+M6+M4+M2+M3` and balanced regular grading `60+60`.
- **DERIVED:** the FS indicators in McKay-chain order are
  `(+,-,+,-,+,-,+,-,+)`: integer-spin irreps are real and spinors are
  quaternionic.
- **DERIVED / STRUCTURAL distinction:** the two weak `M2` factors and two
  color `M3` factors are Galois pairs, but selecting concrete embeddings and
  allocating the other seven factors to obtain a unital
  `C+H+M3(C)` representation is not canonical.
- **DERIVED decisive negative:** inversion on `C[2I]` has `J^2=+1` and
  commutes with isotypic parity.  Galois composition preserves parity too, so
  the canonical proposal has `J gamma=+gamma J`, not the KO6 minus sign.
- **DERIVED negative:** for the maximal multiplicity algebra the regular
  bimodule has only diagonal Krajewski vertices; first order admits no
  off-diagonal odd block, hence `D=0`.
- **DERIVED scope:** coordinate/preprojective arrows are parity-odd and evade
  the edge orientation argument, while existing equivariant vertex operators
  are parity-even.  A smaller hand-picked SM algebra could change legality,
  but no canonical unital restriction exists, so that equation is undefined
  rather than failed.
- **DERIVED boundary theorem:** SM matter is not derivable from the canonical
  node, edge, or smooth preprojective-fiber modules of this geometry, with
  three separate certificates recorded in `preprojective_matter.md`.
- **OPEN:** any noncanonical smaller allocation, a parity-reversing real
  structure, color orientation, generation blocks, and `Y`; Route C cannot
  start without them.
- **Verification:** `verify_preprojective_matter.py` passes `22/22` internal
  checks; the registered suite completed `38/38 PASS` with
  `PYTHONPATH=/tmp/science-python-deps`.
- No PDF build was attempted.
## Session 2026-07-24: full-complex Kaehler--Dirac matter audit

- Reinstalled `numpy`, `scipy`, and `sympy` under
  `/tmp/science-python-deps`; the prescribed pre-change baseline completed
  `41/41 PASS` in 100.1 s.
- Added `kahler_dirac_matter.md` and
  `reproducible/verify_kahler_dirac.py`; registered the verifier.
- **DERIVED:** the exact sparse oriented complex has f-vector
  `(120,720,1200,600)`, integral `d^2=0`, numerical ranks `(119,601,599)`,
  and Betti numbers `(1,0,0,1)`.
- **DERIVED:** `D=d+d*` is nonzero, self-adjoint, fully `2I`-equivariant,
  odd for form parity, and even for the independent central/spin parity.
  Therefore the segregation theorem is not violated and does not kill it.
- **DERIVED:** `c0=2640` is exactly `Tr(I_C)=dim C` for this same cochain
  Hilbert space, not a coincidence.
- **DERIVED:** exact signed cell characters give
  `C^0=Reg`, `C^1=6Reg`, `C^2=10Reg`, `C^3=5Reg`, hence `C=22Reg`.
  Total McKay-chain multiplicities are
  `(22,44,66,88,110,132,88,44,66)`.
- **DERIVED negative at the matter gate:** weak/color-sized factors are
  abundant but no invariant carrier/allocation selects a unital
  `C+H+M3(C)` action. Full multiplicity `U(m)` naturality leaves only
  block scalars, so canonical order-zero/first-order SM data are undefined;
  choosing subspaces would insert the missing bimodule by hand.
- **DERIVED correction:** a Hodge star maps primal to dual cells and cannot
  be an invertible `C^k -> C^(3-k)` endomorphism here because the paired
  primal dimensions differ. Coefficient conjugation and antipodal
  conjugation both have sign table `(+,+,+)`.
- **DERIVED taste facts:** `ker D` has dimension 2 and the rest is 1319
  positive/negative pairs. **PATTERN rejected:** no uniform threefold taste,
  generation count, or distinguished `C3` is selected.
- **OPEN:** a derived SM subalgebra, justified primal-dual real structure,
  `Y`, Route C, and the multiplet census. No PDF build was attempted.

## Session 2026-07-24: hostile adversarial audit

- Confirmed the prescribed pre-audit baseline:
  `PYTHONPATH=/tmp/science-python-deps python3 reproducible/run_all.py`
  completed `46/46 PASS` in 243.8 s.
- Added `adversarial_audit_report.md` with explicit
  `ATTACK -> METHOD -> OUTCOME` logs for all eight priority targets.
- **BROKEN:** the holographic/tower/warped “preregistration” provenance is
  not auditable: the workspace has no version-control history or immutable
  external registry, and the older protocol lacks the numerical thresholds.
  Notes now say fixed-for-rerun with unverified pre-result provenance.
- **WEAKENED:** the C3 result is a residual-equivariant
  Krajewski-legal block candidate, not a constructed real spectral triple or
  matrix-level first-order Dirac.  Added an independent common-character
  check.
- **WEAKENED:** the A5 edge-metric representative is conditional on the
  displayed sampling embedding because the verifier does not compare it to
  the actual 720-edge kernel intertwiner.
- **WEAKENED:** `gap=phi^-4` and the mass-ratio negative are numerical
  double-precision certificates, not exact-algebraic theorems.
- **STRENGTHENED:** invariant spectral moments now have an independent exact
  integer-matrix path giving `(14880,55920)`.
- **STRENGTHENED:** the Krajewski verifier now asserts the advertised
  dimension-only counts `188,908,396` and `1,362,811,872,984` exactly.
- Rechecked cited primary arXiv records for JUNO, DESI, the dynamical-DE
  analysis, and KATRIN.  The NuFIT arXiv record was reachable but its official
  table PDF could not be re-fetched; the note records that limitation.
- Standalone corrected checks passed: segregation, A5 brackets, Bratteli,
  bimodule Krajewski, JUNO, Kähler--Dirac, and invariant spectrum.
- Final post-correction registered suite completed `46/46 PASS` in 236.9 s
  with `PYTHONPATH=/tmp/science-python-deps`.
- No PDF build was attempted.
## Session 2026-07-27, sixth session: modular/TQFT foundational audit

- Read the binding 2026-07-22--27 audit/provenance record, the S02 closure
  note, the paper/supplement TQFT passages, and the Fibonacci-tower audit.
  Dependencies were intact.  The untouched registered baseline completed
  `53/53 PASS` in `147.3 s`.
- Added `modular_tqft_layer.md` and
  `reproducible/verify_modular_tqft_layer.py`; registered the latter as test
  54.  It gives 25 exact symbolic checks.
- **DERIVED / SCOPE FIX:** the requested original bootstrap uses the
  fundamental dimension `d_(1/2)(n)=2 cos(pi/n)`, whereas the later S02 note
  silently used the spin-one dimension
  `d_1(n)=1+2 cos(2pi/n)`.  The equation
  `2 cos(pi/n)=(1+sqrt(n))/2` has the unique integer solution `n=5`:
  `n>=9` is excluded by bounds and `3<=n<=8` is checked exactly.
- **VERDICT:** uniqueness conditional on the chosen matching is DERIVED; the
  matching prescription is STRUCTURAL; its content is the pentagon identity,
  not a dynamical vacuum-selection principle.  Physical selection remains
  OPEN.
- **DERIVED modular data:** computed fusion, dimensions, total dimensions,
  `S`, twists/`T`, and central charges for Fib=`(G2)_1`, `(F4)_1`, `(E8)_1`,
  and `SU(2)_3`.  The full `SU(2)_3` Verlinde ring has rank four; only its
  even Fib subring is `Z[phi]`.
- **DERIVED new bridge:** verified
  `(G2)_1 x (F4)_1 subset (E8)_1`, including `14/5+26/5=8`, branching
  vector `(1,0,0,1)`, exact `S` invariance, and `T` compatibility from
  `2/5+3/5=1`.  Connecting this extension to 600-cell operators is
  STRUCTURAL/OPEN.
- **DERIVED negative:** the Fib Galois conjugate is non-unitary Yang--Lee,
  not the complementary unitary `(F4)_1` factor; their twists and categorical
  dimensions differ.
- **DERIVED, right for the wrong reason:** the formal global-dimension ratio
  gives exactly `log(D_Fib/D_YL)=log(phi)`.  Calling it physical
  visible/dark topological entanglement is OPEN because no non-unitary entropy
  prescription or dark ground state/cut is constructed.
- **DERIVED arithmetic / PATTERN physics:** `||N_(1/2)||_F^2=6` for
  `SU(2)_3`, while the Fib matrix gives 3.  The `S` prefactor and norm do not
  derive `c_eff=sqrt(2/5)(1-6 alpha^2)` without a propagator, contraction,
  and fermion map.  The radiative interpretation is PATTERN.
- **PATTERN:** `N_gen=3` and `k=3` are an equality of unrelated counts; a
  categorical generation functor remains OPEN.
- **DERIVED positive:** the Fib based fusion ring and the stationary tower
  `K0` are canonically the same ordered ring (explicit basis swap, same PF
  positive cone).  Physical mass/scale interpretation remains STRUCTURAL.
- Corrected the abstract, introduction/bootstrap, radiative, dark/entropy,
  conclusion, and supplementary TQFT wording; updated S02's binding closure
  status.  No PDF build was attempted.
- Final registered suite:
  `PYTHONPATH=/tmp/science-python-deps python3 reproducible/run_all.py`
  completed in `167.8 s` with
  `Result: 54/54 scripts completed successfully.`

## Session 2026-08-09--10: dimension reconciliation

- Repaired the verifier registry in commit `12952c9`: the parent contained
  79 entries but only 78 distinct names, with
  `verify_incidence_operator_enumeration.py` as the sole duplicate.  The
  coverage guard now rejects duplicates as well as missing registrations and
  missing files.
- Added the independent registered reconstruction in commit `7246e9d`.
  **DERIVED:** the used carrier is the three-dimensional boundary of the
  ambient four-dimensional 600-cell, with `f=(120,720,1200,600)`, Betti
  numbers `(1,0,0,1)`, and exactly two harmonic Kahler--Dirac modes.
- **DERIVED NEGATIVE:** the complete fixed `D^2` heat flow has global maximum
  `3.295663771`, so it contains no four-dimensional diffusion regime.
  **PATTERN:** the exhaustive target-free 2% shoulder is `3.274268` over only
  `0.376` decade and does not meet the registered half-decade plateau gate.
- **DERIVED NEGATIVE:** `9-13+1-1=-4` is an alternating nullity of four
  operators that do not intertwine the simplicial differential, not a
  cohomological or spacetime index.  The actual Euler characteristic is zero.
- **DERIVED NEGATIVE:** `(2640,14880,55920)` are exact finite heat moments.
  The finite heat trace is analytic at zero, so these are not
  Seeley--DeWitt coefficients and carry no dimensional exponent.
- Corrected `logic_chain_map.md` and the paper/supplement labels.  No audited
  probe on the fixed carrier supports a fourth static dimension; a fourth
  direction, if retained, requires new dynamical structure.
- Commit `f7adc8c` gives the exhaustive chamber census its measured timeout.
  Final unique registered suite:
  `/home/razvan/science/.venv/bin/python reproducible/run_all.py` completed
  `79/79 PASS` in `901.0 s` with process exit `0`.
- No PDF build was attempted.

## Session 2026-08-10: exact refutation of chamber conjecture B1

- **DERIVED REFUTATION:** on the unchanged 120-chamber `D,J,gamma`, the
  algebra `A=M2(C)+C^3` has an exact faithful unital noncommutative
  representation satisfying order zero, first order, grading compatibility,
  nonzero one-forms, metric-dimension-zero orientability, nondegenerate KO6
  intersection form, and connectedness.
- The positive-sheet multiplicities are
  `mu_01=2, mu_12=25, mu_31=12, mu_23=19`, with cell dimensions
  `(4,25,12,19)`. An exact 60-vertex colouring sends all 90 edges of
  `S=(D J)|H+` into composable cell pairs. The `M2` node is a pure source, so
  every shared first-order factor is scalar; full order zero and first order
  are nevertheless checked on all seven complex basis elements.
- **DERIVED:** the orientability cycle equals `gamma` entrywise; the
  intersection matrix has rank four, Pfaffian `38`, and determinant `1444`;
  the exact commutator map has rank six and one-dimensional scalar kernel.
- **DERIVED framing corrections:** total faithfulness does not imply the old
  `k<=7` bound (one-sheet support permits `k=60`). Under orientability, first
  order, and invertibility of fixed `S`, each occupied cell has dimension at
  most 30, giving the corrected general bound `k<=30`. Odd summand parity
  remains a valid necessary PD obstruction.
- **DERIVED NEGATIVE:** orientability and connectedness are not functions of
  the Wedderburn type. The old failing witness and the new all-gate witness
  have the same type `M2(C)+C^3`; type-only enumeration cannot decide B1.
- Corrected the previous “first-order detector”: it enforced only the
  Krajewski support mask and omitted identity-intertwiner constraints on
  non-scalar shared factors. Its exact census remains a cell-mask screen and
  its numerical output remains only a PATTERN.
- **OPEN:** B1 strengthened to demand an integrally unimodular intersection
  form; geometry-selection of the algebra; a color/matter/Standard-Model
  sector. The exact witness proves possibility, not physical selection.
- Added and registered `verify_chamber_b1_counterexample.py` (`17/17`) in
  commit `4c05edb`. Final unique registered suite completed `80/80 PASS` in
  `803.2 s`, process exit `0`.
- No PDF build was attempted.
