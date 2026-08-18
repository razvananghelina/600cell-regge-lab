# PROMPT: Derive Neutrino Masses from a₁ = 5

## Why this matters
If the framework predicts neutrino masses BEFORE they are precisely measured,
it's a genuine prediction — not post-hoc fitting. This is the single most
important thing we can do to counter the "Texas Sharpshooter" critique.

Current experimental status:
- Mass splittings measured: Δm²₂₁ = 7.53 × 10⁻⁵ eV², |Δm²₃₂| = 2.453 × 10⁻³ eV²
- Absolute masses UNKNOWN (only upper bounds: Σmᵢ < 0.12 eV from cosmology)
- Mass ordering UNKNOWN (normal vs inverted hierarchy)
- KATRIN direct limit: m_νe < 0.45 eV (2024)

So we CAN predict: absolute masses, mass ordering, and lightest mass.
These are testable by KATRIN, JUNO, DESI, and Euclid within ~5 years.

## What we have from the framework

### Constants
- a₁ = 5, b₁ = 6, φ = (1+√5)/2, N = 120, N_gen = 3, N_eig = 9
- m_e = 0.51099895 MeV (dimensional anchor)
- Mass formula: m_f = m_e · φ^(5a + 6b) with norm-log corrections

### PMNS angles (already derived and verified)
- sin²(θ₂₃) = 4/7 = (a₁-1)/(a₁+2)
- sin²(θ₁₂) = 2/(φ+a₁) = 2/(φ+5)
- sin²(θ₁₃) = 1/(a₁·N_eig) = 1/45
- δ_PMNS = 3·arctan(√5) ≈ 199.1°

### Charged lepton assignments
- e:   (a,b) = (0,0),  n = 0
- μ:   (a,b) = (1,1),  n = 11
- τ:   (a,b) = (1,2),  n = 17

### Key constraint: neutrinos are MUCH lighter than charged leptons
- m_ν < 0.1 eV vs m_e = 0.511 MeV → ratio > 5 × 10⁶
- In phi units: φ^n where n ~ 32-35 below electron
- So neutrino exponents must be LARGE AND NEGATIVE

## Strategy: BLIND DERIVATION (no peeking at experimental masses)

### THE CORRECT PROCEDURE

**Step 1: Topological deduction (a priori, no data)**
Use the SAME geometric rules that selected charged fermion (a,b) to
determine neutrino (a,b). No scanning, no fitting, no searching.

Primary route: GALOIS CONJUGATE MECHANISM
- Charged leptons use z = a + b·φ with mass m ∝ φ^n
- Neutrinos use the Galois image: z' = a + b·φ' with mass m ∝ |φ'|^n = φ^(-n)
- SAME quantum numbers as partner lepton, different Galois sector

This gives IMMEDIATELY (no search needed):
  ν_e partner of e:   (a,b) = (0,0), n = 0 → m_νe = m_e · |φ'|^0 = m_e ??? 
  
  NO — the mechanism must be more subtle. The Galois conjugate acts on the
  MASS FORMULA, not trivially. Investigate:
  
  Option A: n_ν = -n_charged (Galois flip of exponent)
    ν_e: n=0, ν_μ: n=-11, ν_τ: n=-17
    m_νe = m_e, m_νμ = m_e·φ^(-11), m_ντ = m_e·φ^(-17)
    These are HUGE (m_νe = 0.511 MeV). Wrong.
    
  Option B: Seesaw via Galois — m_ν = m_e²/M where M = m_e·φ^N_Galois
    The Galois structure provides the heavy scale M naturally.
    
  Option C: Neutrino exponent = Galois conjugate of CHARGED exponent
    n_charged uses φ-sector eigenvalues. Neutrino uses φ'-sector.
    From Galois kernel theorem: ker has ρ₀⊕ρ₁⊕ρ₈
    ρ₁ is φ-sector (charged), ρ₈ is φ'-sector (neutral?)
    Eigenvalue ratio: λ(ρ₁)/λ(ρ₈) = 6φ/(6-6φ) = φ/φ' = -φ² 
    
  Option D: Use L(ρ₈)/L(ρ₁) ratio (Laplacian eigenvalues)
    L(ρ₁) = 12-6φ = 2.292, L(ρ₈) = 12-(6-6φ) = 6+6φ = 15.708
    Ratio L₈/L₁ = (6+6φ)/(12-6φ) = ... = φ'⁻² ??? Compute this.
    
  TASK: Compute all Galois-related ratios and find which one provides
  the neutrino mass suppression relative to charged leptons.

**Step 2: PMNS constraint (fix mass ratios)**
The PMNS angles are ALREADY DERIVED:
  sin²(θ₁₂) = 2/(φ+5), sin²(θ₂₃) = 4/7, sin²(θ₁₃) = 1/45

These angles diagonalize the neutrino mass matrix:
  M_ν = U_PMNS · diag(m₁, m₂, m₃) · U_PMNS^T

The mixing angles constrain the mass ratios. Use the PMNS matrix
to determine m₁:m₂:m₃ from the geometric structure, BEFORE looking
at any experimental mass splitting.

Specifically: if neutrino masses come from a matrix M_ν that is
diagonalized by U_PMNS, the structure of U_PMNS (which we derived
from A₅) constrains which mass patterns are allowed.

**Step 3: BLIND comparison (the supreme test)**
Only AFTER completing Steps 1 and 2 (with zero experimental input 
on neutrino masses), compute:
  - Δm²₂₁ = m₂² - m₁²
  - Δm²₃₂ = m₃² - m₂²
  - Σmᵢ = m₁ + m₂ + m₃
  - m_β = sqrt(Σ |U_ei|² mᵢ²) (for KATRIN)
  - m_ββ = |Σ U_ei² mᵢ| (for 0νββ)

Then compare with experiment. The comparison is meaningful ONLY
because we did NOT use the experimental values in Steps 1-2.

## Computational tasks

### Task 1: Galois suppression mechanism (Step 1)
Compute ALL Galois-related ratios on the 600-cell that could provide
the neutrino mass scale:
- L(ρ₈)/L(ρ₁) = (6+6φ)/(12-6φ) = ?
- λ(ρ₁)/λ(ρ₈) = 6φ/(6-6φ) = ?
- Eigenvalue products in φ'-sector vs φ-sector
- Galois kernel theorem: what does ρ₈ contribute vs ρ₁?
- The Hopf fiber Galois structure: fiber eigenvalue φ vs φ'

Find the UNIQUE ratio R such that m_ν = m_charged · R^(something).
This ratio must come from geometry, not fitting.

### Task 2: Assign neutrino quantum numbers (Step 1 continued)
Using the Galois mechanism from Task 1:
- Same (a,b) as partner charged lepton? Or Galois-conjugated (a,b)?
- Compute the three neutrino masses from pure geometry
- Record these BEFORE looking at experimental splittings

### Task 3: PMNS mass ratio constraint (Step 2)
Build the neutrino mass matrix from our derived PMNS angles:
- Construct U_PMNS from sin²θ₁₂=2/(φ+5), sin²θ₂₃=4/7, sin²θ₁₃=1/45
- Determine what mass ratios m₁:m₂:m₃ are compatible with this U
- Check consistency with Task 2 masses

### Task 4: Blind comparison (Step 3)
ONLY after Tasks 1-3, compute:
- Δm²₂₁ = m₂² - m₁²  → compare with (7.53 ± 0.18) × 10⁻⁵ eV²
- |Δm²₃₂| = |m₃² - m₂²| → compare with (2.453 ± 0.033) × 10⁻³ eV²
- Σmᵢ → compare with < 0.12 eV (Planck) or < 0.072 eV (DESI)
- Mass ordering: normal (m₁<m₂<m₃) or inverted?
- m_β for KATRIN, m_ββ for 0νββ

### Task 5: Document the result
Three possible outcomes (all valuable):
A) Framework predicts masses consistent with experiment → MAJOR RESULT
B) Framework predicts masses outside current bounds → FALSIFIABLE PREDICTION  
C) Framework cannot uniquely determine masses → HONEST NEGATIVE, document why

## What counts as success
- A UNIQUE or small set of (a,b) triplets that simultaneously:
  1. Give masses in the sub-eV range
  2. Match Δm²₂₁ and |Δm²₃₂| within experimental errors
  3. Are consistent with the framework's PMNS angles
  4. Follow the same Z[φ] algebra as charged fermions

- Bonus: predict mass ordering (normal/inverted) and absolute scale

## What counts as failure (document honestly)
- No (a,b) triplet matches both splittings simultaneously
- Multiple triplets match equally well (no predictive power)
- The framework requires ad hoc modifications for neutrinos
- Mass ordering contradicts experimental hints

## Framework constants for reference
a1 = 5, b1 = 6, phi = 1.6180339887, phi' = -0.6180339887
N = 120, N_gen = 3, N_eig = 9, degree = 12, h = 30
C = 2/13, c_ell = C·φ³/4, sin²θ_W = 6/26
alpha = 1/137.036, alpha_s = 1/(2φ³)
m_e = 0.51099895 MeV = 511.0 keV = 0.511 × 10⁶ eV

## Experimental constraints (PDG 2024 + NuFIT 5.3)
Δm²₂₁ = (7.53 ± 0.18) × 10⁻⁵ eV²
Δm²₃₂ = (2.453 ± 0.033) × 10⁻³ eV² (normal ordering)
Δm²₃₂ = -(2.536 ± 0.034) × 10⁻³ eV² (inverted ordering)
Σmᵢ < 0.12 eV (Planck 2018, 95% CL)
Σmᵢ < 0.072 eV (DESI + CMB, 2024, preliminary)
m_β < 0.45 eV (KATRIN 2024, 90% CL)
