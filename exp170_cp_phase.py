"""
EXP-170: CKM/PMNS CP Phase from Golden Ratio
==============================================
The CP-violating phase in CKM: delta = 1.144 +/- 0.027 rad = 65.5 +/- 1.5 deg (PDG 2022)
The CP-violating phase in PMNS: delta = ~197 +/- 40 deg (poorly measured)

Can we express delta_CKM and/or delta_PMNS as functions of phi?

Approach: systematic scan of phi-based expressions and compare to experiment.

Category: RESEARCH (searching for PATTERN)
"""

import numpy as np
from itertools import product as iterproduct

PHI = (1 + np.sqrt(5)) / 2

print("=" * 70)
print("EXP-170: CKM/PMNS CP Phase from Golden Ratio")
print("=" * 70)

# ============================================================
# Step 1: Experimental values
# ============================================================
print("\n--- Step 1: Experimental Values ---")

# CKM CP phase (PDG 2022)
delta_CKM_rad = 1.144  # +/- 0.027
delta_CKM_deg = np.degrees(delta_CKM_rad)
delta_CKM_err = 0.027  # rad
delta_CKM_err_deg = np.degrees(delta_CKM_err)

print(f"  CKM delta = {delta_CKM_rad:.3f} +/- {delta_CKM_err:.3f} rad")
print(f"           = {delta_CKM_deg:.1f} +/- {delta_CKM_err_deg:.1f} deg")

# Also: Wolfenstein parameters give gamma angle
# gamma = arg(-V_ud V_ub* / V_cd V_cb*) = 65.5 +/- 3 deg
# In standard parametrization, delta_13 = delta_CKM
gamma_deg = 65.5
gamma_err = 3.0

# PMNS CP phase (NOvA + T2K combined)
# Very uncertain: delta_PMNS ~ 197 +/- 40 deg (from global fits)
# Some fits give ~230 deg, others ~195 deg
delta_PMNS_deg = 197.0  # central value
delta_PMNS_err = 40.0
delta_PMNS_rad = np.radians(delta_PMNS_deg)

print(f"\n  PMNS delta = {delta_PMNS_rad:.2f} +/- {np.radians(delta_PMNS_err):.2f} rad")
print(f"            = {delta_PMNS_deg:.0f} +/- {delta_PMNS_err:.0f} deg")

# ============================================================
# Step 2: CKM CP phase candidates
# ============================================================
print("\n" + "=" * 70)
print("Step 2: CKM CP Phase Candidates")
print("=" * 70)

candidates_CKM = []

# Category A: Simple arctan expressions
for n in range(-5, 6):
    if n == 0:
        val_deg = 45.0
        formula = "arctan(1)"
    else:
        val_deg = np.degrees(np.arctan(PHI**n))
        formula = f"arctan(phi^{n})"
    err = abs(val_deg - delta_CKM_deg)
    err_pct = err / delta_CKM_deg * 100
    candidates_CKM.append((err_pct, val_deg, formula, 'arctan'))

# Category B: Multiples of arctan(phi^n)
for n in range(-4, 5):
    for m in range(2, 8):
        val_deg = m * np.degrees(np.arctan(PHI**n))
        if 0 < val_deg < 180:
            err_pct = abs(val_deg - delta_CKM_deg) / delta_CKM_deg * 100
            formula = f"{m}*arctan(phi^{n})"
            candidates_CKM.append((err_pct, val_deg, formula, 'multi-arctan'))

# Category C: pi/phi^n and multiples
for n in range(-3, 5):
    val_deg = np.degrees(np.pi / PHI**n)
    if 0 < val_deg < 360:
        err_pct = abs(val_deg - delta_CKM_deg) / delta_CKM_deg * 100
        formula = f"pi/phi^{n}"
        candidates_CKM.append((err_pct, val_deg, formula, 'pi-phi'))

# Category D: Combinations of Cabibbo angle
theta_C = np.degrees(np.arctan(PHI**(-3)))  # 13.28 deg
for m in range(2, 8):
    val_deg = m * theta_C
    err_pct = abs(val_deg - delta_CKM_deg) / delta_CKM_deg * 100
    formula = f"{m}*theta_C ({m}*arctan(phi^-3))"
    candidates_CKM.append((err_pct, val_deg, formula, 'cabibbo-mult'))

# Category E: arctan of simple phi expressions
phi_exprs = [
    (2*PHI, "2phi"),
    (PHI+1, "phi+1 = phi^2"),
    (PHI-1, "phi-1 = 1/phi"),
    (2*PHI-1, "2phi-1 = sqrt(5)"),
    (PHI**2-1, "phi^2-1 = phi"),
    (2/PHI, "2/phi"),
    (PHI/2, "phi/2"),
    (3*PHI-2, "3phi-2"),
    (PHI**2+PHI, "phi^2+phi = phi+1+phi = 2phi+1"),
]
for val, label in phi_exprs:
    val_deg = np.degrees(np.arctan(val))
    err_pct = abs(val_deg - delta_CKM_deg) / delta_CKM_deg * 100
    formula = f"arctan({label})"
    candidates_CKM.append((err_pct, val_deg, formula, 'phi-expr'))

# Category F: Differences and sums of arctan(phi^n)
for n1 in range(-3, 4):
    for n2 in range(n1+1, 5):
        # arctan(phi^n2) - arctan(phi^n1)
        val = np.degrees(np.arctan(PHI**n2) - np.arctan(PHI**n1))
        if val > 0:
            err_pct = abs(val - delta_CKM_deg) / delta_CKM_deg * 100
            formula = f"arctan(phi^{n2}) - arctan(phi^{n1})"
            candidates_CKM.append((err_pct, val, formula, 'diff-arctan'))
        # sum
        val = np.degrees(np.arctan(PHI**n2) + np.arctan(PHI**n1))
        if 0 < val < 180:
            err_pct = abs(val - delta_CKM_deg) / delta_CKM_deg * 100
            formula = f"arctan(phi^{n2}) + arctan(phi^{n1})"
            candidates_CKM.append((err_pct, val, formula, 'sum-arctan'))

# Category G: arcsin and arccos
for n in range(-4, 4):
    val_sin = PHI**n
    if -1 < val_sin < 1:
        val_deg = np.degrees(np.arcsin(val_sin))
        if val_deg > 0:
            err_pct = abs(val_deg - delta_CKM_deg) / delta_CKM_deg * 100
            formula = f"arcsin(phi^{n})"
            candidates_CKM.append((err_pct, val_deg, formula, 'arcsin'))

    val_cos = PHI**n
    if -1 < val_cos < 1:
        val_deg = np.degrees(np.arccos(val_cos))
        if val_deg > 0:
            err_pct = abs(val_deg - delta_CKM_deg) / delta_CKM_deg * 100
            formula = f"arccos(phi^{n})"
            candidates_CKM.append((err_pct, val_deg, formula, 'arccos'))

# Category H: Angles from regular polygons / icosahedral
special_angles = [
    (36.0, "pi/5 (pentagon)"),
    (72.0, "2pi/5 (pentagon)"),
    (108.0, "3pi/5"),
    (60.0, "pi/3 (hexagon)"),
    (120.0, "2pi/3"),
    (np.degrees(np.arctan(2)), "arctan(2)"),
    (np.degrees(np.arccos(1/3)), "arccos(1/3) (tetrahedron)"),
    (np.degrees(np.arccos(PHI/2)), "arccos(phi/2)"),
    (np.degrees(np.arccos(1/(2*PHI))), "arccos(1/(2phi))"),
    (np.degrees(2*np.arctan(PHI)), "2*arctan(phi)"),
]
for val_deg, label in special_angles:
    err_pct = abs(val_deg - delta_CKM_deg) / delta_CKM_deg * 100
    candidates_CKM.append((err_pct, val_deg, label, 'special'))

# Category I: Using SM parameters from theory
alpha_em = 1/137.036
sin2tw = 6/26
theta_W_deg = np.degrees(np.arcsin(np.sqrt(sin2tw)))

ckm_from_sm = [
    (np.degrees(np.pi * sin2tw), f"pi*sin^2(tW) = pi*6/26"),
    (np.degrees(np.arctan(PHI**3 * alpha_em)), f"arctan(phi^3 * alpha)"),
    (theta_W_deg, f"theta_W = arcsin(sqrt(6/26))"),
    (2*theta_W_deg, f"2*theta_W"),
    (90 - theta_W_deg, f"90 - theta_W"),
]
for val_deg, label in ckm_from_sm:
    err_pct = abs(val_deg - delta_CKM_deg) / delta_CKM_deg * 100
    candidates_CKM.append((err_pct, val_deg, label, 'SM-param'))

# Sort by error
candidates_CKM.sort(key=lambda x: x[0])

print(f"\n  Top 20 CKM CP phase candidates (target: {delta_CKM_deg:.1f} deg):")
print(f"  {'err%':>6} {'value':>8} {'formula':<50} {'type':<15}")
print("  " + "-" * 85)
for err_pct, val_deg, formula, cat in candidates_CKM[:20]:
    within = "< 1sigma" if err_pct < delta_CKM_err_deg/delta_CKM_deg*100 else ""
    print(f"  {err_pct:>5.2f}% {val_deg:>8.2f} {formula:<50} {cat:<15} {within}")

# Check best candidate within 1sigma
sigma_pct = delta_CKM_err_deg / delta_CKM_deg * 100
n_within_sigma = sum(1 for e,_,_,_ in candidates_CKM if e < sigma_pct)
print(f"\n  1sigma threshold: {sigma_pct:.2f}%")
print(f"  Candidates within 1sigma: {n_within_sigma}")

# ============================================================
# Step 3: PMNS CP phase candidates
# ============================================================
print("\n" + "=" * 70)
print("Step 3: PMNS CP Phase Candidates")
print("=" * 70)

candidates_PMNS = []

# Similar scan but for delta_PMNS ~ 197 deg
target = delta_PMNS_deg

# arctan-based (remembering delta > 180 means we're in 3rd/4th quadrant)
# delta_PMNS ~ 197 = 180 + 17
# So the "effective angle" from pi is about 17 degrees

# Category A: 180 + arctan(phi^n)
for n in range(-5, 5):
    val_deg = 180 + np.degrees(np.arctan(PHI**n))
    err_pct = abs(val_deg - target) / target * 100
    formula = f"180 + arctan(phi^{n})"
    candidates_PMNS.append((err_pct, val_deg, formula, '180+arctan'))

# Category B: k * pi/5 (icosahedral)
for k in range(1, 11):
    val_deg = k * 36.0
    err_pct = abs(val_deg - target) / target * 100
    formula = f"{k}*pi/5 = {k}*36"
    candidates_PMNS.append((err_pct, val_deg, formula, 'icosahedral'))

# Category C: n * arctan(phi^m)
for n in range(2, 16):
    for m in range(-3, 4):
        val_deg = n * np.degrees(np.arctan(PHI**m))
        if 150 < val_deg < 250:
            err_pct = abs(val_deg - target) / target * 100
            formula = f"{n}*arctan(phi^{m})"
            candidates_PMNS.append((err_pct, val_deg, formula, 'multi-arctan'))

# Category D: pi + specific phi expressions
phi_based = [
    (180 + np.degrees(np.arctan(PHI**(-3))), "180 + arctan(phi^-3) = 180 + theta_C"),
    (180 + 2*np.degrees(np.arctan(PHI**(-3))), "180 + 2*theta_C"),
    (360/PHI, "360/phi"),
    (360/PHI**2, "360/phi^2"),
    (180*PHI/PHI**2, "180/phi"),
    (360 - 360/PHI, "360*(1-1/phi)"),
    (180 + 180/PHI**3, "180 + 180/phi^3"),
    (np.degrees(PHI * np.pi), "phi*pi (in degrees)"),
    (np.degrees(np.pi + 2*np.arctan(PHI**(-3))), "pi + 2*arctan(phi^-3)"),
]
for val_deg, label in phi_based:
    if 100 < val_deg < 300:
        err_pct = abs(val_deg - target) / target * 100
        candidates_PMNS.append((err_pct, val_deg, label, 'phi-special'))

# Category E: Fractions of 360
for num in range(1, 20):
    for den in [PHI, PHI**2, PHI**3, 5, 2*PHI, 3, 7, 12]:
        val_deg = 360 * num / den
        if 150 < val_deg < 250:
            err_pct = abs(val_deg - target) / target * 100
            if isinstance(den, float):
                formula = f"360*{num}/phi^? "
            else:
                formula = f"360*{num}/{den}"
            candidates_PMNS.append((err_pct, val_deg, f"360*{num}/{den:.4f}", 'fraction'))

# Category F: Using PMNS angles
theta_23_pmns = np.degrees(np.arctan(1))  # 45 for n=0
theta_12_pmns = np.degrees(np.arctan(1/PHI))  # ~31.7 for n=1
theta_13_pmns = np.degrees(np.arctan(PHI**(-4)))  # ~8.3 for n=4

pmns_combos = [
    (theta_23_pmns + theta_12_pmns + theta_13_pmns + 180, "180 + theta_23 + theta_12 + theta_13"),
    (4*theta_23_pmns + theta_12_pmns, "4*theta_23 + theta_12"),
    (6*theta_12_pmns, "6*theta_12_PMNS"),
    (180 + theta_12_pmns/2, "180 + theta_12/2"),
    (theta_23_pmns + 3*theta_23_pmns + theta_12_pmns, "4*45 + 31.7"),
]
for val_deg, label in pmns_combos:
    if 150 < val_deg < 250:
        err_pct = abs(val_deg - target) / target * 100
        candidates_PMNS.append((err_pct, val_deg, label, 'pmns-combo'))

# Sort
candidates_PMNS.sort(key=lambda x: x[0])

print(f"\n  Top 20 PMNS CP phase candidates (target: {target:.0f} +/- {delta_PMNS_err:.0f} deg):")
print(f"  {'err%':>6} {'value':>8} {'formula':<55} {'type':<15}")
print("  " + "-" * 90)
for err_pct, val_deg, formula, cat in candidates_PMNS[:20]:
    within = "< 1sigma" if abs(val_deg - target) < delta_PMNS_err else ""
    print(f"  {err_pct:>5.2f}% {val_deg:>8.2f} {formula:<55} {cat:<15} {within}")

# ============================================================
# Step 4: Deep analysis of best CKM candidates
# ============================================================
print("\n" + "=" * 70)
print("Step 4: Analysis of Best CKM Candidates")
print("=" * 70)

# Focus on candidates with < 2% error
best_ckm = [(e,v,f,c) for e,v,f,c in candidates_CKM if e < 2.0]
print(f"\n  Candidates with < 2% error ({len(best_ckm)}):")
for err_pct, val_deg, formula, cat in best_ckm:
    val_rad = np.radians(val_deg)
    print(f"  {err_pct:>5.2f}%  {val_deg:>8.3f} deg  {val_rad:>6.4f} rad  {formula}")

# The most "natural" candidate
print(f"\n  --- Special analysis: 5*theta_Cabibbo ---")
val = 5 * np.degrees(np.arctan(PHI**(-3)))
err = abs(val - delta_CKM_deg)
err_pct = err / delta_CKM_deg * 100
print(f"  5*arctan(phi^-3) = {val:.4f} deg")
print(f"  Experimental: {delta_CKM_deg:.1f} +/- {delta_CKM_err_deg:.1f} deg")
print(f"  Error: {err:.2f} deg = {err_pct:.2f}%")
print(f"  Within 1sigma: {err < delta_CKM_err_deg}")
print(f"  Interpretation: delta_CP = 5 * theta_Cabibbo")
print(f"  The factor 5 = a_1 (intersection number of 600-cell)")

print(f"\n  --- Special analysis: arctan(phi^2) ---")
val = np.degrees(np.arctan(PHI**2))
err = abs(val - delta_CKM_deg)
err_pct = err / delta_CKM_deg * 100
print(f"  arctan(phi^2) = {val:.4f} deg")
print(f"  Error: {err:.2f} deg = {err_pct:.2f}%")
print(f"  phi^2 = phi+1 = {PHI**2:.4f}")

print(f"\n  --- Special analysis: 2*arctan(phi) ---")
val = 2*np.degrees(np.arctan(PHI))
err = abs(val - delta_CKM_deg)
err_pct = err / delta_CKM_deg * 100
print(f"  2*arctan(phi) = {val:.4f} deg")
print(f"  Error: {err:.2f} deg = {err_pct:.2f}%")

# ============================================================
# Step 5: Connection to geometry
# ============================================================
print("\n" + "=" * 70)
print("Step 5: Geometric Interpretation")
print("=" * 70)

# Key angles in the 600-cell
print("\n  Key angles in 600-cell geometry:")
# Dihedral angle of tetrahedron
dihedral_tet = np.degrees(np.arccos(1/3))
print(f"  Tetrahedral dihedral: arccos(1/3) = {dihedral_tet:.2f} deg")

# Dihedral angle of 600-cell
# Between adjacent tetrahedral cells
dihedral_600 = np.degrees(np.pi - np.arccos(1/(2*PHI**2)))
print(f"  600-cell dihedral (approx): {dihedral_600:.2f} deg")

# Angle subtended by an edge at center
edge_angle = np.degrees(np.arccos(PHI/2))
print(f"  Edge angle (arccos(phi/2)): {edge_angle:.2f} deg")
print(f"  Compare to delta_CKM: {delta_CKM_deg:.1f} deg")
err_pct = abs(edge_angle - delta_CKM_deg) / delta_CKM_deg * 100
print(f"  Error: {err_pct:.1f}%")

# Vertex figure angles
# The vertex figure is an icosahedron
# Edge-midpoint angle in icosahedron
ico_edge = np.degrees(np.arctan(PHI))  # ~58.3 deg
print(f"\n  Icosahedron edge-to-center: arctan(phi) = {ico_edge:.2f} deg")
err_pct = abs(ico_edge - delta_CKM_deg) / delta_CKM_deg * 100
print(f"  Error vs delta_CKM: {err_pct:.1f}%")

# ============================================================
# CONCLUSIONS
# ============================================================
print("\n" + "=" * 70)
print("CONCLUSIONS")
print("=" * 70)

best = candidates_CKM[0]
print(f"\nCKM delta_CP = {delta_CKM_deg:.1f} deg:")
print(f"  Best candidate: {best[2]} = {best[1]:.2f} deg ({best[0]:.2f}% error)")
print()

# Highlight the most physically meaningful
print(f"  Most meaningful candidates:")
for e, v, f, c in candidates_CKM[:10]:
    if c in ('cabibbo-mult', 'phi-expr', 'arctan', 'special', 'SM-param'):
        meaningful = " ***" if e < 2 else ""
        print(f"    {e:>5.2f}%  {v:>8.2f}  {f}{meaningful}")

print(f"\nPMNS delta_CP = {delta_PMNS_deg:.0f} deg (poorly measured):")
best_pmns = candidates_PMNS[0]
print(f"  Best candidate: {best_pmns[2]} = {best_pmns[1]:.2f} deg ({best_pmns[0]:.2f}% error)")
print(f"  (Note: PMNS delta has ~20% uncertainty, so many candidates match)")

print(f"""
CATEGORY ASSESSMENT:
- If best CKM < 1%: potential PATTERN
- If best CKM 1-3%: SUGGESTIVE
- If best CKM > 3%: likely NOISE (too many free expressions)
""")
