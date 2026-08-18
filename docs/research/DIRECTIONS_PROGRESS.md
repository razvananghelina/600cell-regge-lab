# Directions Progress - Post-Review Action Items
# Source: directii_16_02_2026.txt + 16_feb_2026_review.txt
# Created: 2026-02-16

## STATUS LEGEND
- [ ] NOT STARTED
- [~] IN PROGRESS
- [x] DONE
- [-] SKIPPED (with reason)

---

## I. MUST-FIX BEFORE SUBMISSION

### 1. [x] Spectral action -> dim=12 -> Theorem 1 chain -- exp330
**Gap**: Theorem 1 assumes compact Lie group with dim=12. Why dim=12?
**Fix**: dim=12 = vertex degree = 2(a_1+1), DERIVED from a_1=5. Permutation rep
12 = 1+3+3'+5 under A5. Adjoint identification 1+3'+3+5 = 1+3+8. Uniqueness of
partition 11=3+8. Three DISTINCT YM coefficients 8:5:2 confirm 3 gauge sectors.
**Result**: exp330 verified complete chain. dim(G)=12 is NOT assumed but DERIVED.
**Location**: Section 4, add Remark after Theorem 1
**Status**: DONE (computation complete, paper edit pending)

### 2. [x] Distler-Garibaldi disclaimer -- exp330
**Gap**: Reviewer may think framework uses E8 embedding a la Lisi
**Fix**: Clarified 4 points: (1) E8 via McKay graph, not gauge group. (2) Fermions
from 2I irreps, not E8 reps. (3) Gauge group from A5 decomposition (Thm 1). (4) E8
enters only through McKay topology, edge decomposition, Seeley-DeWitt coefficients.
**Location**: Section 4 or Discussion
**Status**: DONE (argument formulated, paper edit pending)

### 3. [x] m_W consistency (experimental vs derived) -- exp330
**Gap**: Paper uses m_W experimental for m_H prediction sometimes, derived other times
**Result**: m_H = 124.6 GeV (from derived m_W=79.90, err 0.39%) and m_H = 125.36 GeV
(from exp m_W=80.377, err 0.20% = 2.3 sigma vs ATLAS 125.11). Both within tolerance.
**Location**: Section 7 (boson masses), present BOTH values
**Status**: DONE (computation complete, paper edit pending)

### 4. [x] Update experimental values (JUNO + ATLAS) -- exp330
**Gap**: Paper uses m_H = 125.25 (old), sin2_th12 from pre-JUNO
**Result**: All updated values verified. Key changes:
  - m_H: 125.25+-0.17 -> 125.11+-0.11 (ATLAS 2023, arXiv:2308.04775)
  - sin2_th12: 0.303+-0.012 -> 0.3092+-0.0087 (JUNO 2025, arXiv:2511.14593)
  - sin2_th13: 0.02219+-0.00062 -> 0.02225+-0.00056 (NuFIT 6.0)
  - All predictions remain consistent. sin2_th12 IMPROVES (0.8 sigma).
**Location**: Table 1, Table 2, relevant equations
**Status**: DONE (values verified, paper edit pending)

---

## II. HIGH-PRIORITY EXPERIMENTS

### 5. [x] D_F eigenvalues on McKay graph (Directia 2) -- exp329
**Goal**: Construct finite Dirac operator D_F on E8-hat McKay graph with weights
from Theorem 3 ([3,2,6,0,5] main chain). Compute eigenvalues. Check if proportional
to phi^n_f (mass exponents).
**Result**: CONFIRMED NEGATIVE. 6 schemes tested, all fail. 14 positive eigenvalues
(need 9), max dynamic range 630 (need 271443). Best corr=0.975 but structurally wrong.
D_F = selection rules, Z[phi] = scale. COMPLEMENTARY structures.
**Status**: DONE

### 6. [ ] Remark: CC exponent N/2 = |A5| via spectral balance
**Insight**: N/2 = |A5| = |2I/Z2|, where Z2 cancellation is the spectral balance
(Proposition spectral_balance). So 57 = |A5| - N_gen.
**Fix**: Add Remark in CC section (Section 12)
**Effort**: Small (1 paragraph)
**Impact**: Structural motivation for exponent 57 (currently unmotivated)

### 7. [x] Alpha derivation from Hopf fiber spectral action -- exp331
**Gap**: Factor pi/a_1 in self-energy is ad-hoc
**Result**: pi/a_1 DERIVED from spectral gap of Hopf fiber C_{2*a_1} Laplacian:
  lambda_1 = 1/phi^2 => theta = arccos(1-lambda_1/2) = arccos(phi/2) = pi/a_1.
  KEY IDENTITY: cos(pi/5) = phi/2, SPECIFIC to a_1=5.
  Hopf fibration verified: 12 fibers * 10 vertices, all C_10 cycles.
  Edge angles: ALL exactly pi/5 (verified numerically). Total holonomy = 2*pi.
  Fiber-Cayley DUALITY: SUM L(3)+L(3')=2*a_1 counts fiber edges (topological),
  PRODUCT L(3)*L(3')=4*a_1 gives tree coupling (metric). Same eigenvalue pair!
  120 fiber edges + 600 cross edges = 720 total. Ratio 600/120 = a_1 = 5.
**Status**: DONE (pi/a_1 upgraded from ad-hoc to DERIVED)

### 5b. [x] Product Dirac D_F x Delta_0 on 3600-dim space (Direction A) -- exp332
**Goal**: Construct D = D_ext (x) I_30 + Gamma (x) D_F on 3600-dim Hilbert space.
Gamma = Hopf T3 grading introduces sector-dependent scaling. Cross-terms between
D_F (which knows fermion identity) and Delta_0 (which knows 600-cell geometry)
could generate mass hierarchy.
**Result**: CONFIRMED NEGATIVE. 6 forms tested:
  - Best dynamic range: phi^14.6 (Form 4: Adj+T3), still far from phi^26
  - Best RMS: 2.80 (Form 4), WORSE than random baseline (mean 1.79 +- 0.76)
  - Z-score: -1.33 (NOT significant), 94% of random targets fit better
  - Cross-term: 16.1% of D_Hopf^2 (nonzero but not dominant)
  - D_ext dominates spectrum (98% of Frobenius norm in D^2)
  - T3 grading introduces ratio phi:1:1/phi between chirality levels (interesting but insufficient)
**Structural insight**: D_ext (120-dim, eigenvalues in Z[phi]) provides scale but is
INDIFFERENT to fermion identity. D_F (30-dim) knows fermion identity but has tiny
dynamic range. Product structure cannot amplify D_F's range through D_ext.
**Status**: DONE

### 8. [ ] m_H/m_W correction: why 8 = rank(E8)?
**Gap**: -8*alpha correction identified as rank(E8) but could be dim(SU(3))
**Fix**: Derive from spectral action on icosahedron: multiplicity of lambda=3
eigenspace on 600-cell is 16, but adj representation is 8-dim. Show which.
**Effort**: Medium
**Impact**: Strengthens Higgs mass derivation

---

## III. MEDIUM-PRIORITY (strengthen paper)

### 9. [ ] Mention Furey-Hughes 2024 + triality-D4 connection
**What**: Cite Phys. Lett. B 2024. Note complementarity: they derive N_gen=3
from division algebra triality, we derive from Fibonacci/Galois. D4 node
of E8 Dynkin = triality symmetry S3.
**Location**: Discussion section
**Effort**: 1 paragraph + reference

### 10. [ ] Mention Singh 2025 (J3(O_C) mass ratios)
**What**: Cite arXiv:2508.10131. Note: first Singh relation sqrt(m_tau/m_mu) =
sqrt(m_s/m_d) = phi^3 is EXACTLY satisfied in our framework (trivially,
both exponent differences = b_1 = 6). Second relation 1:2:3 is NOT satisfied
(we give 1:phi^(3/2):phi^(5/2) = 1:2.058:3.330).
**Location**: Discussion section
**Effort**: 1 paragraph + reference

### 11. [ ] Neutrino mass factor 2: derive from |2I/A5| = 2
**Gap**: Factor 2 in m_3 = 2*m_e/phi^35 motivated as "Majorana = 2 DOF"
**Fix**: Derive from |2I/A5| = 2 (binary extension quotient) or from eta-invariant
double counting (symmetric eigenvalues).
**Effort**: Small (rewrite motivation)

### 12. [ ] Mention Boyle-Mygdalas spacetime quasicrystals (2026)
**What**: arXiv:2601.07769. 600-cell as 4D quasicrystal prototype (Elser-Sloane 1987).
**Location**: Discussion
**Effort**: 1 sentence + reference

---

## IV. LONG-TERM DIRECTIONS (post-submission)

### 13. [ ] Random NCG on 600-cell (Directia 1)
**What**: Path-integral over Dirac operators on 600-cell simplicial complex.
Fix f_2*Lambda^2 at phase transition critical point.
**Test**: <Tr(D^2)> in random ensemble = 14880 (= c_1)?
**Effort**: Large (new program)

### 14. [x] Cohn-Kumar optimality -> m_e scale (Directia 3 / Direction F) -- exp333
**What**: Equate 600-cell gravitational minimum energy to SM Casimir energy on S^3.
Solve for R, which fixes the scale and hence m_e.
**Result**: CONFIRMED NEGATIVE (structural no-go).
  - S_2 = 5395 = 5*13*83 (RATIONAL, Galois-invariant, phi cancels exactly)
  - Energy balance E_grav~1/R^2 vs E_Casimir~1/R gives R=O(1) pure number
  - Need R~10^22 for m_e/m_P~10^-23. Polynomial ratios CANNOT produce this.
  - The hierarchy requires EXPONENTIAL mechanism (alpha^z), not power-law balance.
  - Also tested Direction E (theta function on Z[phi]): NEGATIVE (linear form not modular,
    indefinite norm, no correlation between norm and mass hierarchy)
  - POSITIVE: Lambda^2 = c_1/(2*c_0) = 31/11 = (a1^2+b1)/(a1+b1) (clean expression)
**Status**: DONE

### 15. [ ] Triality D4 -> 3 legs of E8 (Directia 4)
**What**: Map E8 legs (dim 11,6,9) to triality factorizations (vector, S+, S-)
**Effort**: Medium (algebraic)

### 16. [ ] Beta functions from spectral action cutoff variation (Directia 5)
**What**: Compute dc_k/dLambda on 600-cell. Compare with one-loop SM beta functions.
**Partial**: exp326 already derived b_1, b_2, b_3 from N_gen=3, N_H=1.
**Remaining**: Derive from spectral action directly.
**Effort**: Medium

---

## COMPLETED (from earlier work, relevant to directions)

- [x] exp324: Gauge Group Identification Theorem (partial fix for item 1)
- [x] exp325: Chirality (48+48 split)
- [x] exp326: Beta function coefficients DERIVED
- [x] exp327: Electron mass CONFIRMED NEGATIVE (7 routes)
- [x] exp328: CC spectral action CONFIRMED NEGATIVE (6 routes)
- [x] exp328: Proposition spectral_balance added to paper
- [x] Singh relation tested: first EXACT (phi^3), second NOT (phi^(3/2) != 2)
- [x] exp332: Product Dirac D_F x Delta_0 CONFIRMED NEGATIVE (Direction A from alte_directii.txt)
