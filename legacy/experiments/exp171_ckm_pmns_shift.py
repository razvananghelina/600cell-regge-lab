"""
EXP-171: CKM/PMNS Exponent Relationship
========================================
CKM mixing exponents: {3, 7, 12} (for theta_12, theta_23, theta_13)
PMNS mixing exponents: {1, 0, 4} (for theta_12, theta_23, theta_13)

Known geometric origins:
  CKM:  3 = dim(Im(H)), 7 = dim(Im(O)), 12 = degree(600-cell)
  PMNS: 0 = identity,   1 = edge/radius, 4 = dim(R^4)

Questions:
1. Is there a transformation CKM_i = f(PMNS_i)?
2. Is the shift related to color (SU(3))?
3. Can we find a unified formula?
4. What do sums, products, differences reveal?

Category: RESEARCH (searching for PATTERN)
"""

import numpy as np
from itertools import product, combinations
from collections import Counter

PHI = (1 + np.sqrt(5)) / 2
ALPHA = 1/137.036

print("=" * 70)
print("EXP-171: CKM/PMNS Exponent Relationship")
print("=" * 70)

# ============================================================
# Step 1: Basic data
# ============================================================
print("\n--- Step 1: Mixing Exponents ---")

# Format: angle_label: (CKM_n, PMNS_n, CKM_deg, PMNS_deg)
mixing_data = {
    'theta_12': {'ckm_n': 3, 'pmns_n': 1,
                 'ckm_deg': np.degrees(np.arctan(PHI**(-3))),
                 'pmns_deg': np.degrees(np.arctan(PHI**(-1)))},
    'theta_23': {'ckm_n': 7, 'pmns_n': 0,
                 'ckm_deg': np.degrees(np.arctan(PHI**(-7))),
                 'pmns_deg': np.degrees(np.arctan(PHI**(0)))},
    'theta_13': {'ckm_n': 12, 'pmns_n': 4,
                 'ckm_deg': np.degrees(np.arctan(PHI**(-12))),
                 'pmns_deg': np.degrees(np.arctan(PHI**(-4)))},
}

print(f"\n  {'angle':>10} {'CKM_n':>6} {'PMNS_n':>7} {'CKM_deg':>9} {'PMNS_deg':>10} {'diff_n':>7} {'ratio_n':>8}")
print("  " + "-" * 65)
for angle, data in mixing_data.items():
    diff = data['ckm_n'] - data['pmns_n']
    ratio = data['ckm_n'] / data['pmns_n'] if data['pmns_n'] != 0 else float('inf')
    print(f"  {angle:>10} {data['ckm_n']:>6} {data['pmns_n']:>7} {data['ckm_deg']:>9.3f} {data['pmns_deg']:>10.3f} {diff:>7} {ratio:>8.3f}")

ckm = [3, 7, 12]
pmns = [1, 0, 4]
diffs = [c - p for c, p in zip(ckm, pmns)]
print(f"\n  CKM exponents:  {ckm}")
print(f"  PMNS exponents: {pmns}")
print(f"  Differences:    {diffs}")

# ============================================================
# Step 2: Arithmetic relationships
# ============================================================
print("\n" + "=" * 70)
print("Step 2: Arithmetic Relationships")
print("=" * 70)

print(f"\n  --- Sums ---")
print(f"  CKM sum:  3+7+12 = {sum(ckm)}")
print(f"  PMNS sum: 1+0+4  = {sum(pmns)}")
print(f"  Total:    {sum(ckm)+sum(pmns)}")
print(f"  Ratio:    {sum(ckm)}/{sum(pmns)} = {sum(ckm)/sum(pmns):.4f}")
print(f"  phi^3 = {PHI**3:.4f}")
print(f"  Sum ratio vs phi^3: {abs(sum(ckm)/sum(pmns) - PHI**3)/PHI**3*100:.1f}%")

print(f"\n  --- Products ---")
print(f"  CKM product:  3*7*12 = {np.prod(ckm)}")
print(f"  PMNS product: 1*0*4  = {np.prod(pmns)} (degenerate due to 0)")

print(f"\n  --- Differences (CKM - PMNS for each angle) ---")
print(f"  diffs = {diffs}")
print(f"  Sum of diffs: {sum(diffs)}")
print(f"  Product of diffs: {np.prod(diffs)}")
print(f"  2 * 7 * 8 = {2*7*8}")
print(f"  Compare: dim(SU(2))=3, dim(SU(3))=8, diff_sum=17")

print(f"\n  --- Ratios ---")
print(f"  12/3 = {12/3} (CKM: theta_13/theta_12 exponents)")
print(f"  12/7 = {12/7:.4f}")
print(f"  7/3  = {7/3:.4f}")
print(f"  4/1  = {4/1} (PMNS: theta_13/theta_12 exponents)")
print(f"  Compare: 12/3 = 4 = 4/1 (same ratio!)")

# ============================================================
# Step 3: Test specific hypotheses
# ============================================================
print("\n" + "=" * 70)
print("Step 3: Specific Hypotheses")
print("=" * 70)

# Hypothesis A: CKM = PMNS + color_correction
print(f"\n  --- Hypothesis A: CKM_i = PMNS_i + color_shift_i ---")
for angle, data in mixing_data.items():
    shift = data['ckm_n'] - data['pmns_n']
    print(f"  {angle}: shift = {shift}")
print(f"  Shifts = {diffs}")
print(f"  Not constant -> simple additive shift FAILS")

# Hypothesis B: CKM_i = k * PMNS_i + offset
print(f"\n  --- Hypothesis B: CKM_i = k * PMNS_i + c ---")
# Solve: 3 = k*1 + c, 12 = k*4 + c (skip theta_23 since PMNS_n=0)
# From these: 12 - 3 = k*(4-1) => k = 9/3 = 3, c = 3 - 3 = 0
# Predict theta_23 CKM: 3*0 + 0 = 0 (actual: 7) -> FAILS
k_12_13 = (12 - 3) / (4 - 1)
c_12_13 = 3 - k_12_13 * 1
pred_23 = k_12_13 * 0 + c_12_13
print(f"  From theta_12 and theta_13: k={k_12_13:.1f}, c={c_12_13:.1f}")
print(f"  Predicted CKM theta_23: {pred_23:.1f} (actual: 7)")
print(f"  Error: {abs(pred_23 - 7):.1f} -> Linear model FAILS")

# Hypothesis C: CKM_i = PMNS_i + SU(3)_dim_related
print(f"\n  --- Hypothesis C: Color-related shifts ---")
print(f"  shift_12 = 2 (= 2 = dim(U(1))*2?)")
print(f"  shift_23 = 7 (= 7 = dim(Im(O)))")
print(f"  shift_13 = 8 (= 8 = dim(SU(3)))")
print(f"  Pattern: 2, 7, 8 -> no obvious single structure")

# Hypothesis D: CKM_i = f(PMNS_i) with phi
print(f"\n  --- Hypothesis D: CKM = phi-transform of PMNS ---")
tests = [
    ("PMNS*phi + c", lambda p: p*PHI),
    ("PMNS*phi^2 + c", lambda p: p*PHI**2),
    ("PMNS^2 + c", lambda p: p**2),
    ("Fibonacci", None),  # Special
]

# Try: is {3,7,12} a Fibonacci-like transform of {0,1,4}?
fib = [0, 1, 1, 2, 3, 5, 8, 13, 21]
print(f"  Fibonacci sequence: {fib}")
print(f"  CKM exponents in Fibonacci? 3={'YES' if 3 in fib else 'NO'}, 7={'YES' if 7 in fib else 'NO'}, 12={'YES' if 12 in fib else 'NO'}")
print(f"  PMNS exponents in Fibonacci? 0={'YES' if 0 in fib else 'NO'}, 1={'YES' if 1 in fib else 'NO'}, 4={'YES' if 4 in fib else 'NO'}")

# Hypothesis E: Both sets from a single sequence
print(f"\n  --- Hypothesis E: Single underlying sequence ---")
all_exponents = sorted(set(ckm + pmns))
print(f"  All exponents sorted: {all_exponents}")  # [0, 1, 3, 4, 7, 12]
print(f"  Gaps: {[all_exponents[i+1]-all_exponents[i] for i in range(len(all_exponents)-1)]}")

# Check if these are related to icosahedral numbers
print(f"\n  Icosahedral/geometric numbers:")
ico_numbers = {
    'edges_tet': 6, 'verts_tet': 4, 'faces_tet': 4,
    'edges_oct': 12, 'verts_oct': 6, 'faces_oct': 8,
    'edges_ico': 30, 'verts_ico': 12, 'faces_ico': 20,
    'degree_600': 12, 'cells_600': 600, 'verts_600': 120,
    'a_1': 5, 'b_1': 6,
    'dim_Im_H': 3, 'dim_Im_O': 7, 'dim_R4': 4,
}
for name, val in sorted(ico_numbers.items(), key=lambda x: x[1]):
    in_ckm = val in ckm
    in_pmns = val in pmns
    marker = ""
    if in_ckm: marker += " [CKM]"
    if in_pmns: marker += " [PMNS]"
    if marker:
        print(f"    {name}: {val}{marker}")

# ============================================================
# Step 4: The 12/3 = 4/1 coincidence
# ============================================================
print("\n" + "=" * 70)
print("Step 4: The Ratio Coincidence 12/3 = 4/1 = 4")
print("=" * 70)

print(f"""
  CKM theta_13/theta_12 exponent ratio: 12/3 = 4
  PMNS theta_13/theta_12 exponent ratio: 4/1 = 4

  SAME RATIO! The 13/12 mixing hierarchy is IDENTICAL in both sectors.

  This means: sin(theta_13)/sin(theta_12) has the same phi-power law
  for both quarks and leptons:
    sin(theta_13)/sin(theta_12) ~ phi^(-n_13)/phi^(-n_12) = phi^(n_12-n_13)

  For CKM:  phi^(3-12) = phi^(-9)
  For PMNS: phi^(1-4)  = phi^(-3)

  Ratio of ratios: phi^(-9)/phi^(-3) = phi^(-6)

  So the SECTOR DEPENDENCE is captured by a single factor phi^(-6)!
""")

print(f"  phi^(-6) = {PHI**(-6):.6f}")
print(f"  Also: phi^(-6) = (1/phi)^6 = phi_conj^6 (with appropriate sign)")
print(f"  And: 6 = b_1 = intersection number of 600-cell")
print(f"")

# Test: CKM_n = PMNS_n + 6 * something?
for angle, data in mixing_data.items():
    ratio = (data['ckm_n'] - data['pmns_n']) / 6 if data['pmns_n'] != data['ckm_n'] else 0
    print(f"  {angle}: (CKM-PMNS)/6 = ({data['ckm_n']}-{data['pmns_n']})/6 = {ratio:.4f}")

print(f"\n  Not exactly multiples of 6. But interesting pattern.")

# ============================================================
# Step 5: The CKM origins revisited
# ============================================================
print("\n" + "=" * 70)
print("Step 5: Geometric Origins Comparison")
print("=" * 70)

print(f"""
  PMNS (leptons see GEOMETRY):
    n=0: Identity (trivial representation)
    n=1: Edge/radius ratio = 1/phi (fundamental geometric ratio)
    n=4: dim(R^4) (spacetime dimension)

  CKM (quarks see ALGEBRA):
    n=3: dim(Im(H)) = 3 (imaginary quaternions)
    n=7: dim(Im(O)) = 7 (imaginary octonions)
    n=12: degree(600-cell) = 12 (graph connectivity)

  Key observation: PMNS uses {0, 1, 4} = simple geometry
                   CKM uses {3, 7, 12} = algebraic structures

  The ALGEBRA is deeper than the GEOMETRY.
  Quarks, having color charge (SU(3)), "see" the algebraic structure.
  Leptons, being color-neutral, "see" only the geometric structure.
""")

# Is there a unifying pattern?
print(f"  --- Unification attempts ---")

# attempt 1: All from division algebras
# R: dim=1 (PMNS theta_12? n=1)
# C: dim=2 (?)
# H: dim=4 (PMNS theta_13? n=4), dim(Im(H))=3 (CKM theta_12)
# O: dim=8 (dim(SU(3))=8), dim(Im(O))=7 (CKM theta_23)
print(f"\n  Division algebra dimensions: R(1), C(2), H(4), O(8)")
print(f"  Imaginary dimensions: Im(R)=0, Im(C)=1, Im(H)=3, Im(O)=7")
print(f"  PMNS uses: 0, 1, 4 -> dim(Im(R)), dim(Im(C)), dim(H)")
print(f"  CKM uses:  3, 7, 12 -> dim(Im(H)), dim(Im(O)), degree(600-cell)")
print(f"")
print(f"  PMNS shift: Im(R)->Im(C)->H (total algebra dimension)")
print(f"  CKM shift:  Im(H)->Im(O)->graph_degree")

# The pattern: PMNS uses the FULL algebra dimension (0,1,4)
# while CKM uses the IMAGINARY part dimension (3,7) plus graph degree
# For n=12: is there an algebra with Im dimension 12? No standard one.
# But: 12 = dim(SU(2)) * dim(H) = 3*4? Or 12 = degree of graph?

print(f"\n  --- Attempt: Division algebra hierarchy ---")
print(f"  Lepton sector (PMNS): exponents from total dim of algebras")
print(f"    R=0, C=1(?), H=4  -> {0, 1, 4}")
print(f"    Actually: 0=trivial, 1=dim(R), 4=dim(H)")
print(f"")
print(f"  Quark sector (CKM): exponents from imaginary part")
print(f"    Im(H)=3, Im(O)=7, ??? = 12")
print(f"    12 = 3+4+5 (triangle number)?")
print(f"    12 = dim(SU(2)xU(1)) + dim(SU(3)) = 4+8?")
print(f"    12 = degree of 600-cell vertex (GEOMETRIC, not algebraic!)")

# ============================================================
# Step 6: The theta_23 anomaly
# ============================================================
print("\n" + "=" * 70)
print("Step 6: The theta_23 Anomaly")
print("=" * 70)

print(f"""
  theta_23 is special in both sectors:

  PMNS: n=0 -> theta_23 = arctan(1) = 45 deg (MAXIMAL mixing!)
  CKM:  n=7 -> theta_23 = arctan(phi^-7) = 2.0 deg (very small)

  The PMNS theta_23 is the ONLY angle with n=0 (identity).
  The CKM theta_23 uses n=7 (octonion imaginary dimension).

  Is n=7 related to n=0 through color?
  7 = 0 + 7? Shift by 7 = dim(Im(O))

  Meanwhile:
  theta_12: CKM n=3 = PMNS n=1 + 2
  theta_13: CKM n=12 = PMNS n=4 + 8
  theta_23: CKM n=7 = PMNS n=0 + 7

  Shifts: 2, 7, 8

  2 + 7 + 8 = 17 = tau level (n=17 for tau!)

  Also: 2*7*8 = 112 = 16*7 (16 = SM multiplet states)
""")

# Is 17 meaningful?
print(f"  17 = n_tau = 5*1 + 6*2 (tau level in mass formula)")
print(f"  Coincidence? Or connection between mixing and mass?")
print(f"")

# Check: are {2, 7, 8} related to 600-cell geometry?
print(f"  The numbers 2, 7, 8:")
print(f"  2 = diameter of icosahedron vertex figure at distance 2?")
print(f"  7 = dim(Im(O)) = number of imaginary octonion units")
print(f"  8 = dim(SU(3)) = number of gluons")
print(f"  {2,7,8} in spectral structure: eigenvalue -2 has mult 36")

# ============================================================
# Step 7: Unified exponent formula
# ============================================================
print("\n" + "=" * 70)
print("Step 7: Unified Exponent Formula Attempts")
print("=" * 70)

# Attempt 1: n_CKM = n_PMNS + dim(Im(D_i))
# where D_1 = C (dim Im = 1), D_2 = O (dim Im = 7), D_3 = ??? (dim Im = 8?)
print(f"  Attempt 1: n_CKM(theta_ij) = n_PMNS(theta_ij) + s_ij")
print(f"    s_12 = 2 = dim(C) = 2")
print(f"    s_23 = 7 = dim(Im(O))")
print(f"    s_13 = 8 = dim(su(3)) = dim(o_imaginary) + 1?")
print(f"    Pattern: 2, 7, 8 = ?, Im(O), su(3)")
print(f"    NOT a clean pattern from division algebras alone")

# Attempt 2: All from graph theory
print(f"\n  Attempt 2: Both from 600-cell graph distances")
print(f"    600-cell has 9 distance classes (diameter 5)")
print(f"    Vertex counts at distance d from any vertex:")
print(f"    d=0: 1, d=1: 12, d=2: 42, d=3: 44, d=4: 20, d=5: 1")
print(f"    Cumulative: 1, 13, 55, 99, 119, 120")

# Are exponents related to cumulative counts?
cum = [1, 13, 55, 99, 119, 120]
for n in [0, 1, 3, 4, 7, 12]:
    matches = [(d, cum[d]) for d in range(6) if cum[d] == n or d == n]
    print(f"    n={n}: distance d={n} count = {[1,12,42,44,20,1][min(n,5)] if n<=5 else 'N/A'}")

# Attempt 3: Factorization
print(f"\n  Attempt 3: Prime factorization")
print(f"    CKM: 3=3, 7=7, 12=2^2*3")
print(f"    PMNS: 0=0, 1=1, 4=2^2")
print(f"    CKM without PMNS factor: 3/1=3, 7/0=inf, 12/4=3")
print(f"    12/4 = 3 = 3/1: consistent factor of 3 for theta_12 and theta_13!")
print(f"    theta_23 is the anomaly (0 in denominator)")

print(f"\n  *** KEY FINDING: For theta_12 and theta_13: CKM_n = 3 * PMNS_n ***")
print(f"    theta_12: 3 = 3 * 1  CHECK")
print(f"    theta_13: 12 = 3 * 4  CHECK")
print(f"    theta_23: 7 != 3 * 0  (anomalous)")
print(f"")
print(f"  The factor 3 = dim(Im(H)) = dim(SU(2)) = number of colors!")
print(f"  Quarks feel 3x the mixing of leptons (for 12 and 13 channels)!")
print(f"  theta_23 is special: PMNS has n=0 (maximal mixing),")
print(f"  so the quark version uses a different mechanism entirely.")

# ============================================================
# Step 8: Testing CKM = 3*PMNS hypothesis
# ============================================================
print("\n" + "=" * 70)
print("Step 8: Testing CKM_n = 3 * PMNS_n + correction")
print("=" * 70)

print(f"\n  For theta_23: if CKM_n = 3*PMNS_n + delta, then delta = 7")
print(f"  For theta_12 and theta_13: delta = 0")
print(f"")
print(f"  Alternative: CKM_n = 3*PMNS_n for 12,13; CKM_23 = 7 (Im(O)) independently")
print(f"  This is a 'mixed' formula: two from multiplication, one from octonions")

# Physical meaning of factor 3
print(f"\n  Physical meaning of factor 3:")
print(f"  - 3 = number of color charges in SU(3)")
print(f"  - 3 = dimension of imaginary quaternions (SU(2) generators)")
print(f"  - 3 = number of generations (!!)")
print(f"  - Quarks carry color => mixing 'amplified' by color factor?")

# Numerical check
print(f"\n  Numerical verification:")
for angle in ['theta_12', 'theta_13']:
    data = mixing_data[angle]
    pred_ckm = 3 * data['pmns_n']
    actual_ckm = data['ckm_n']
    print(f"  {angle}: 3 * {data['pmns_n']} = {pred_ckm}, actual = {actual_ckm}, match = {pred_ckm == actual_ckm}")

# Predicted CKM angles from PMNS * 3
print(f"\n  Predicted CKM angles from n_CKM = 3*n_PMNS:")
for angle, data in mixing_data.items():
    pred_n = 3 * data['pmns_n']
    pred_deg = np.degrees(np.arctan(PHI**(-pred_n)))
    actual_deg = data['ckm_deg']
    # Experimental
    exp_vals = {'theta_12': 13.04, 'theta_23': 2.38, 'theta_13': 0.201}
    exp = exp_vals[angle]
    err_pred = abs(pred_deg - exp) / exp * 100
    err_actual = abs(actual_deg - exp) / exp * 100
    print(f"  {angle}: n_pred={pred_n:>2} -> {pred_deg:>8.3f} deg (err={err_pred:>5.1f}%)"
          f"  vs actual n={data['ckm_n']:>2} -> {actual_deg:>8.3f} deg (err={err_actual:>5.1f}%)")

# ============================================================
# CONCLUSIONS
# ============================================================
print("\n" + "=" * 70)
print("CONCLUSIONS")
print("=" * 70)

print(f"""
KEY FINDINGS:

1. RATIO COINCIDENCE: CKM_13/CKM_12 = PMNS_13/PMNS_12 = 4
   The 13-to-12 hierarchy is IDENTICAL in both sectors. DERIVED.

2. FACTOR OF 3: For theta_12 and theta_13:
   CKM_n = 3 * PMNS_n  (3=1*3, 12=4*3)
   This is EXACT for 2 out of 3 angles.
   Factor 3 = number of colors = dim(SU(2)) = dim(Im(H)).
   Status: PATTERN (2/3 match, 1 anomalous)

3. THETA_23 ANOMALY: PMNS has n=0 (maximal, identity).
   CKM has n=7 = dim(Im(O)). This angle doesn't follow the factor-3 rule.
   It appears to use the OCTONION structure independently.

4. SHIFT SUM: sum of shifts (2+7+8) = 17 = n_tau. Coincidence or not?

5. GEOMETRIC vs ALGEBRAIC: PMNS exponents come from geometry (0,1,4).
   CKM exponents come from algebra (3,7,12). This reflects the color
   charge difference between leptons and quarks.

STATUS:
- Factor 3 rule: PATTERN (2/3 exact, theta_23 anomalous)
- Ratio 4 coincidence: DERIVED (follows from the exponents)
- Geometric origin of each exponent: DERIVED (from exp085)
- Unified formula: PARTIAL (need to explain theta_23)
""")
