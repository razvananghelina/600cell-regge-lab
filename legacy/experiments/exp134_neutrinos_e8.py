"""
EXP-134: Neutrino Masses from the Second 600-Cell (T-copy) in E8
=================================================================

CONTEXT:
E8 = S union T, where:
  - S = first 600-cell (120 vertices), embedded in R^8 via icosian lattice
  - T = phi'*S = Galois conjugate (120 DIFFERENT vertices)
  - phi' = (1-sqrt(5))/2 = -1/phi

All 9 charged fermions are described by levels in S:
  m_f = m_e * phi^(5a + 6b) = m_e * phi^n

HYPOTHESIS: Right-handed neutrinos / neutrino mass mechanism live in T.

INVESTIGATION:
1. T-copy spectral structure (what polytope does phi'*S form in R^4?)
2. Galois conjugate mass formula: m_T = m_e * phi'^n
3. DSI with inverted ladder
4. Type-I seesaw with T-scale
5. E8 cross-coupling S-T
6. Several neutrino mass prediction mechanisms
7. The number 3 from E8 branching

EXPERIMENTAL DATA:
  dm^2_21 = 7.53e-5 eV^2 (solar)
  dm^2_32 = 2.453e-3 eV^2 (atmospheric, normal ordering)
  m_nu3 ~ 0.050 eV, m_nu2 ~ 0.0087 eV (normal hierarchy, m1 ~ 0)

Author: Claude (exp134)
Date: 2026-02-09
"""

import numpy as np
from itertools import product as iterproduct, permutations
from collections import Counter

PHI = (1 + np.sqrt(5)) / 2
PHI_CONJ = (1 - np.sqrt(5)) / 2  # = -1/phi = -0.6180...
ALPHA = 1/137.036
ALPHA_S = 0.1179

# Masses in eV
m_e_eV = 0.51099895e6   # eV
m_mu_eV = 105.6584e6    # eV
m_tau_eV = 1776.86e6    # eV

# Masses in GeV
m_e_GeV = 0.51099895e-3
m_mu_GeV = 105.6584e-3
m_tau_GeV = 1.77686

# Planck mass
m_P_GeV = 1.22089e19

# Neutrino oscillation data (PDG)
dm2_21 = 7.53e-5     # eV^2 (solar)
dm2_32 = 2.453e-3    # eV^2 (atmospheric, normal hierarchy)

# Normal hierarchy estimates (m1 ~ 0)
m_nu3 = np.sqrt(dm2_32)   # ~ 0.0495 eV
m_nu2 = np.sqrt(dm2_21)   # ~ 0.00868 eV
m_nu1 = 0.0               # assume minimal

# Lepton mass levels n = 5a + 6b
lepton_n = {"electron": 0, "muon": 11, "tau": 17}
lepton_ab = {"electron": (0,0), "muon": (1,1), "tau": (1,2)}
lepton_masses_eV = {"electron": m_e_eV, "muon": m_mu_eV, "tau": m_tau_eV}
lepton_masses_GeV = {"electron": m_e_GeV, "muon": m_mu_GeV, "tau": m_tau_GeV}

# All fermion levels
fermion_data = [
    ("electron", (0,0), 0),
    ("muon",     (1,1), 11),
    ("tau",      (1,2), 17),
    ("up",       (3,-2), 3),
    ("charm",    (2,1), 16),
    ("top",      (4,1), 26),
    ("down",     (1,0), 5),
    ("strange",  (1,1), 11),
    ("bottom",   (-1,4), 19),
]

print("=" * 80)
print("EXP-134: Neutrino Masses from the Second 600-Cell (T-copy) in E8")
print("=" * 80)
print()
print(f"  phi  = {PHI:.10f}")
print(f"  phi' = {PHI_CONJ:.10f} = -1/phi")
print(f"  m_e  = {m_e_eV:.2f} eV")
print(f"  m_nu3 (exp) ~ {m_nu3:.5f} eV")
print(f"  m_nu2 (exp) ~ {m_nu2:.5f} eV")
print(f"  dm^2_21 = {dm2_21:.2e} eV^2")
print(f"  dm^2_32 = {dm2_32:.2e} eV^2")
print()

# ============================================================
# SECTION 1: T-copy spectral structure
# ============================================================
print("-" * 80)
print("SECTION 1: T-copy spectral structure in R^4")
print("-" * 80)
print()

def generate_600cell():
    """Generate all 120 vertices of the 600-cell on unit S^3."""
    phi = PHI
    vertices = set()

    # 8 vertices: permutations of (+-1, 0, 0, 0)
    for i in range(4):
        for s in [1, -1]:
            v = [0.0, 0.0, 0.0, 0.0]
            v[i] = float(s)
            vertices.add(tuple(round(x, 10) for x in v))

    # 16 vertices: (+-1/2)^4
    for s0 in [0.5, -0.5]:
        for s1 in [0.5, -0.5]:
            for s2 in [0.5, -0.5]:
                for s3 in [0.5, -0.5]:
                    vertices.add((s0, s1, s2, s3))

    # 96 vertices: EVEN permutations of (+-phi/2, +-1/2, +-1/(2*phi), 0)
    base_vals = [phi/2, 0.5, 1/(2*phi), 0.0]

    even_perms = []
    for p in permutations(range(4)):
        inv = sum(1 for i in range(4) for j in range(i+1, 4) if p[i] > p[j])
        if inv % 2 == 0:
            even_perms.append(p)

    for perm in even_perms:
        for s0 in [1, -1]:
            for s1 in [1, -1]:
                for s2 in [1, -1]:
                    signed = [s0 * base_vals[0], s1 * base_vals[1],
                              s2 * base_vals[2], base_vals[3]]
                    v = tuple(round(signed[perm.index(i)], 10) for i in range(4))
                    vertices.add(v)

    return np.array(sorted(vertices))

verts_S = generate_600cell()
print(f"  S (600-cell) vertices: {len(verts_S)}")
norms_S = np.linalg.norm(verts_S, axis=1)
print(f"  S all norm 1: {np.allclose(norms_S, 1.0)}")

# T vertices in R^4: multiply by phi'
verts_T = verts_S * PHI_CONJ
norms_T = np.linalg.norm(verts_T, axis=1)
print(f"  T = phi'*S vertices: {len(verts_T)}")
print(f"  T norms: min={norms_T.min():.6f}, max={norms_T.max():.6f}")
print(f"  T radius = |phi'| = 1/phi = {abs(PHI_CONJ):.6f}")
print()

# What is the structure of T in R^4?
# T = phi'*S means each coordinate is scaled by phi' = -1/phi
# So T is a SCALED and REFLECTED 600-cell, radius 1/phi instead of 1
# The adjacency structure is THE SAME (just smaller and reflected)

# Build adjacency for T using the same angular criterion
# For S: neighbors have dot(v_i, v_j) = phi/2 (for unit vectors)
# For T: vectors have norm 1/phi, so dot = (1/phi^2) * phi/2 = 1/(2*phi)
dot_threshold_T = 1 / (2 * PHI)  # = PHI_CONJ^2 * PHI/2

# Actually, T = phi'*S. The dot product of phi'*v_i and phi'*v_j is
# phi'^2 * dot(v_i, v_j) = (1/phi^2) * dot(v_i, v_j)
# If v_i, v_j are neighbors in S: dot(v_i, v_j) = phi/2
# Then dot(T_i, T_j) = (1/phi^2) * (phi/2) = 1/(2*phi) = (phi-1)/2 = phi'/(-2)
# Wait: phi'^2 = 1/phi^2? No. phi' = -1/phi, so phi'^2 = 1/phi^2. Yes.

# Build adjacency for the T polytope
adj_T = np.zeros((120, 120), dtype=int)
n_edges_T = 0
# dot products between T vertices
dots_T = verts_T @ verts_T.T
# The neighbor threshold for T is phi'^2 * (phi/2) since T = phi'*S
# phi'^2 = 1/phi^2 = phi-1 (since phi^2 = phi+1, so 1/phi^2 = 1/(phi+1) = phi-1/phi^2...
# Actually 1/phi^2 = 2-phi = 0.38197...)
neighbor_dot_T = PHI_CONJ**2 * (PHI / 2)
print(f"  T neighbor dot product threshold: {neighbor_dot_T:.6f}")

for i in range(120):
    for j in range(i+1, 120):
        if abs(dots_T[i,j] - neighbor_dot_T) < 1e-6:
            adj_T[i,j] = 1
            adj_T[j,i] = 1
            n_edges_T += 1

degrees_T = adj_T.sum(axis=1)
print(f"  T edges: {n_edges_T}, Degree: {degrees_T[0]}")
print(f"  T is isomorphic to S: degree 12, 720 edges: {n_edges_T == 720 and all(d == 12 for d in degrees_T)}")
print()

# Compute all distinct dot products in T
unique_dots_T = sorted(set(round(dots_T[i,j], 6) for i in range(120) for j in range(i+1, 120)))
print(f"  Distinct dot products in T: {len(unique_dots_T)}")
print(f"  Values: {[f'{d:.4f}' for d in unique_dots_T]}")
print()

# Eigenvalues of adjacency matrix for T
eigenvalues_T = np.linalg.eigvalsh(adj_T)
eigenvalues_T_sorted = np.sort(eigenvalues_T)[::-1]
unique_eig_T = []
counts_T = []
prev = None
for e in eigenvalues_T_sorted:
    er = round(e, 4)
    if prev is None or abs(er - prev) > 0.01:
        unique_eig_T.append(er)
        counts_T.append(1)
        prev = er
    else:
        counts_T[-1] += 1

print("  T adjacency spectrum (same graph as S):")
for eig, cnt in zip(unique_eig_T, counts_T):
    # Check if it's a+b*phi form
    b_coeff = round((eig - round(eig)) / (PHI - round(PHI)), 1) if abs(eig - round(eig)) > 0.1 else 0
    print(f"    lambda = {eig:8.4f}, multiplicity = {cnt}")
print()
print("  CONCLUSION: T in R^4 is a SCALED, REFLECTED 600-cell (radius 1/phi).")
print("  Same adjacency graph, same spectrum. Only the geometric embedding differs.")
print()

# ============================================================
# SECTION 2: Galois conjugate mass formula (CORE INVESTIGATION)
# ============================================================
print("-" * 80)
print("SECTION 2: Galois conjugate mass formula")
print("-" * 80)
print()

print("  If S has mass formula m_S = m_e * phi^n,")
print("  then T has mass formula m_T = m_e * (phi')^n = m_e * (-1/phi)^n")
print()
print("  Since |phi'| = 1/phi, we have |m_T| = m_e * phi^(-n)")
print("  The sign alternates: phi'^n = (-1)^n * phi^(-n)")
print()

print("  Direct T-copy masses for leptons:")
print(f"  {'Lepton':10s} {'n':>4s} {'m_S (eV)':>14s} {'|m_T| (eV)':>14s} {'m_T/m_S':>14s}")
print("  " + "-" * 60)

for name in ["electron", "muon", "tau"]:
    n = lepton_n[name]
    m_S = m_e_eV * PHI**n
    m_T = m_e_eV * PHI**(-n)  # |m_T|
    ratio = m_T / m_S if m_S > 0 else float('inf')
    print(f"  {name:10s} {n:4d} {m_S:14.2f} {m_T:14.6f} {ratio:14.6e}")

print()
print("  Comparing T-copy masses to neutrino mass scale:")
print()

m_T_e = m_e_eV * PHI**(0)    # = m_e (same!)
m_T_mu = m_e_eV * PHI**(-11) # = m_e / phi^11
m_T_tau = m_e_eV * PHI**(-17) # = m_e / phi^17

print(f"  T(electron, n=0):  |m_T| = m_e * phi^0    = {m_T_e:.6f} eV = m_e (trivial)")
print(f"  T(muon, n=11):     |m_T| = m_e * phi^(-11) = {m_T_mu:.6f} eV")
print(f"  T(tau, n=17):      |m_T| = m_e * phi^(-17) = {m_T_tau:.6f} eV")
print()
print(f"  Experimental neutrino masses (NH, m1~0):")
print(f"    m_nu3 ~ {m_nu3:.5f} eV")
print(f"    m_nu2 ~ {m_nu2:.5f} eV")
print()

# These T-copy masses are WAY too large for neutrinos
# m_T_mu ~ 2.47 eV (about 50x too large for m_nu3)
# m_T_tau ~ 0.143 eV (about 3x too large for m_nu3)

err_mu = abs(m_T_mu - m_nu3) / m_nu3 * 100
err_tau = abs(m_T_tau - m_nu3) / m_nu3 * 100
print(f"  m_T(mu) vs m_nu3:  {m_T_mu:.4f} vs {m_nu3:.5f} eV  ({err_mu:.0f}% off)")
print(f"  m_T(tau) vs m_nu3: {m_T_tau:.5f} vs {m_nu3:.5f} eV  ({err_tau:.0f}% off)")
print()

# Try different level assignments
print("  Scan: m_e * phi^(-k) for which k matches neutrino masses?")
print(f"  {'k':>4s} {'m_e*phi^(-k) (eV)':>18s} {'Closest nu':>12s} {'Error':>10s}")
print("  " + "-" * 50)
for k in range(0, 40):
    mval = m_e_eV * PHI**(-k)
    for nu_name, nu_mass in [("nu3", m_nu3), ("nu2", m_nu2)]:
        if nu_mass > 0:
            err = abs(mval - nu_mass) / nu_mass * 100
            if err < 30:
                print(f"  {k:4d} {mval:18.6f} {nu_name:>12s} {err:9.1f}%")

print()

# ============================================================
# SECTION 3: DSI with inverted ladder in T-copy
# ============================================================
print("-" * 80)
print("SECTION 3: DSI with inverted ladder")
print("-" * 80)
print()

print("  DSI potential: V(psi) = A * sin^2(pi * ln|psi| / ln(phi))")
print("  Minima at psi = phi^n for integer n.")
print()
print("  In T-copy, replace phi with |phi'| = 1/phi:")
print("  V_T(psi) = A * sin^2(pi * ln|psi| / ln(1/phi))")
print("           = A * sin^2(-pi * ln|psi| / ln(phi))")
print("           = A * sin^2(pi * ln|psi| / ln(phi))  [sin^2 is even]")
print()
print("  SAME potential! The DSI potential is invariant under phi -> 1/phi.")
print("  So T-DSI and S-DSI have the SAME minima at phi^n.")
print("  The 'inverted ladder' idea doesn't create new levels.")
print()
print("  However, the PHYSICAL interpretation differs:")
print("  - S: m_f = m_e * phi^n (n>0 gives heavy fermions)")
print("  - T: Galois conjugate means phi -> phi', so AMPLITUDES get conjugated")
print("  - The mass formula uses the NORM: |a + b*phi| vs |a + b*phi'|")
print()

# Check: for each lepton level, what does the Z[phi] norm give?
print("  Z[phi] norms for each lepton:")
header_phi_conj = "|a+b*phi'|"
print(f"  {'Lepton':10s} {'a':>3s} {'b':>3s} {'n=5a+6b':>8s} {'|a+b*phi|':>12s} {header_phi_conj:>12s} {'Product':>10s}")
print("  " + "-" * 65)
for name, (a, b), n in fermion_data[:3]:  # leptons only
    v1 = abs(a + b * PHI)
    v2 = abs(a + b * PHI_CONJ)
    prod = v1 * v2
    # Z[phi] norm = |a + b*phi| * |a + b*phi'| = |a^2 + ab - b^2| (algebraic norm)
    alg_norm = abs(a**2 + a*b - b**2)
    print(f"  {name:10s} {a:3d} {b:3d} {n:8d} {v1:12.6f} {v2:12.6f} {prod:10.6f}  (alg norm = {alg_norm})")

print()

# ============================================================
# SECTION 4: Type-I Seesaw with T-scale right-handed neutrinos
# ============================================================
print("-" * 80)
print("SECTION 4: Type-I seesaw with T-scale")
print("-" * 80)
print()

print("  Type-I seesaw: m_nu = m_D^2 / M_R")
print()
print("  If m_D = m_lepton (Dirac mass from S) and M_R from T:")
print()

# What M_R gives the right neutrino masses?
print("  Required M_R for each lepton:")
for name in ["electron", "muon", "tau"]:
    m_D = lepton_masses_eV[name]  # eV
    m_nu_target = {"electron": m_nu2, "muon": m_nu2, "tau": m_nu3}[name]
    if m_nu_target > 0:
        M_R = m_D**2 / m_nu_target
        M_R_GeV = M_R / 1e9
        # Express as m_e * phi^k
        k_phi = np.log(M_R / m_e_eV) / np.log(PHI)
        # Express as m_P * phi^(-j)
        j_phi = -np.log(M_R_GeV / m_P_GeV) / np.log(PHI)
        print(f"  {name:10s}: M_R = m_D^2/m_nu = {M_R:.3e} eV = {M_R_GeV:.3e} GeV")
        print(f"              M_R/m_e = phi^{k_phi:.2f}")
        print(f"              M_R = m_P * phi^(-{j_phi:.2f})")

print()

# Universal seesaw: m_D = m_lepton, what single M_R fits all 3?
print("  If M_R is the SAME for all leptons (universal seesaw):")
print("  m_nu_i = m_lepton_i^2 / M_R")
print()

# Fit: dm^2_32 and dm^2_21 ratios
# dm^2_32/dm^2_21 = (m3^2 - m2^2)/(m2^2 - m1^2)
# For m1~0: dm^2_32/dm^2_21 ~ m3^2/m2^2 = (m_tau^4/M_R^2) / (m_mu^4/M_R^2) = (m_tau/m_mu)^4
ratio_exp = dm2_32 / dm2_21
ratio_pred = (m_tau_eV / m_mu_eV)**4
print(f"  dm^2_32/dm^2_21 (exp)  = {ratio_exp:.2f}")
print(f"  (m_tau/m_mu)^4 (pred)  = {ratio_pred:.2f}")
print(f"  Error: {abs(ratio_pred - ratio_exp)/ratio_exp*100:.1f}%")
print()
print("  THIS IS WAY OFF. The universal seesaw with m_D=m_lepton gives")
print("  (m_tau/m_mu)^4 = {:.0f}, but experiment is {:.1f}.".format(ratio_pred, ratio_exp))
print("  The hierarchy is TOO STEEP.")
print()

# What if m_D is NOT m_lepton but m_e * phi^(n/2)?
print("  Alternative: geometric mean seesaw")
print("  If m_D = sqrt(m_e * m_lepton) = m_e * phi^(n/2):")
for name in ["electron", "muon", "tau"]:
    n = lepton_n[name]
    m_D = m_e_eV * PHI**(n/2)
    print(f"    {name:10s}: m_D = m_e * phi^({n}/2) = {m_D:.4f} eV")

print()

# ============================================================
# SECTION 5: E8 cross-coupling S-T
# ============================================================
print("-" * 80)
print("SECTION 5: E8 cross-coupling between S and T")
print("-" * 80)
print()

# Build the icosian decomposition and E8 vectors
def decompose(x):
    """Find (a, b) such that x = a + b*phi, with a,b in half-integers."""
    for a_try in np.arange(-1, 1.01, 0.5):
        for b_try in np.arange(-1, 1.01, 0.5):
            if abs(x - (a_try + b_try * PHI)) < 1e-9:
                return (a_try, b_try)
    raise ValueError(f"Cannot decompose {x}")

def cholesky_map(ab_list):
    """Map (a_i, b_i) -> (a_i + b_i, b_i) for Cholesky embedding."""
    result = np.zeros(2 * len(ab_list))
    for i in range(len(ab_list)):
        a, b = ab_list[i]
        result[2*i] = a + b
        result[2*i+1] = b
    return result

# Decompose S vertices
ab_pairs_S = []
for v in verts_S:
    pairs = [decompose(v[i]) for i in range(4)]
    ab_pairs_S.append(pairs)
ab_pairs_S = np.array(ab_pairs_S)  # (120, 4, 2)

# T decomposition: phi'*(a+b*phi) = (a-b) + (-a)*phi
ab_pairs_T = np.zeros_like(ab_pairs_S)
for idx in range(120):
    for i in range(4):
        a, b = ab_pairs_S[idx, i]
        ab_pairs_T[idx, i] = [a - b, -a]

# Build E8 vectors via Cholesky
S_eucl = np.array([cholesky_map(ab_pairs_S[idx]) for idx in range(120)])
T_eucl = np.array([cholesky_map(ab_pairs_T[idx]) for idx in range(120)])

# Scale by sqrt(2) for E8 convention (norm^2 = 2)
S_e8 = S_eucl * np.sqrt(2)
T_e8 = T_eucl * np.sqrt(2)

# Verify
S_norms = np.sqrt(np.sum(S_e8**2, axis=1))
T_norms = np.sqrt(np.sum(T_e8**2, axis=1))
print(f"  S E8 norms: {np.unique(np.round(S_norms, 4))}")
print(f"  T E8 norms: {np.unique(np.round(T_norms, 4))}")
print()

# Cross inner products S-T
cross_dots = S_e8 @ T_e8.T  # (120, 120) matrix
cross_dot_vals = cross_dots.flatten()
unique_cross = sorted(set(np.round(cross_dot_vals, 4)))
print(f"  Cross dot products S-T (E8 convention):")
for d in unique_cross:
    cnt = np.sum(np.abs(cross_dot_vals - d) < 0.01)
    print(f"    dot = {d:7.3f}: count = {cnt:6d}")

print()

# How many T-neighbors does each S-vertex have? (dot = 1 in E8)
print("  S-T neighbors (dot = 1 in E8):")
cross_neighbors = np.sum(np.abs(cross_dots - 1.0) < 0.01, axis=1)
unique_cn = sorted(set(cross_neighbors))
print(f"    Each S-vertex has {unique_cn} T-neighbors at dot=1")

# Also dot = -1
cross_neighbors_neg = np.sum(np.abs(cross_dots + 1.0) < 0.01, axis=1)
unique_cn_neg = sorted(set(cross_neighbors_neg))
print(f"    Each S-vertex has {unique_cn_neg} T-neighbors at dot=-1")

# Also dot = 0
cross_neighbors_0 = np.sum(np.abs(cross_dots) < 0.01, axis=1)
unique_cn_0 = sorted(set(cross_neighbors_0))
print(f"    Each S-vertex has {unique_cn_0} T-neighbors at dot=0")
print()

# Classify S vertices by type
type_A_S, type_B_S, type_C_S = [], [], []
for i, v in enumerate(verts_S):
    nz = sum(1 for x in v if abs(x) > 0.01)
    if nz == 1 and any(abs(abs(x) - 1) < 0.01 for x in v):
        type_A_S.append(i)
    elif nz == 4 and all(abs(abs(x) - 0.5) < 0.01 for x in v):
        type_B_S.append(i)
    else:
        type_C_S.append(i)

print(f"  S vertex types: A={len(type_A_S)}, B={len(type_B_S)}, C={len(type_C_S)}")

# For each type, check T-neighbors
for label, indices in [("Type-A (8v)", type_A_S), ("Type-B (16)", type_B_S), ("Type-C (96)", type_C_S)]:
    t_neigh = [cross_neighbors[i] for i in indices]
    print(f"  {label}: T-neighbors(dot=1): mean={np.mean(t_neigh):.1f}, unique={sorted(set(t_neigh))}")

print()

# ============================================================
# SECTION 6: Neutrino mass predictions - multiple mechanisms
# ============================================================
print("-" * 80)
print("SECTION 6: Neutrino mass predictions (7 mechanisms)")
print("-" * 80)
print()

results = []  # Store (method_name, m_nu1, m_nu2, m_nu3, dm2_21_pred, dm2_32_pred)

# --- Method A: Direct T-copy masses ---
print("  METHOD A: Direct T-copy: m_nu = m_e * phi^(-n_lepton)")
print("  " + "-" * 50)
m_nu_A = {}
for name in ["electron", "muon", "tau"]:
    n = lepton_n[name]
    m_nu_A[name] = m_e_eV * PHI**(-n)

m1_A = m_nu_A["electron"]  # = m_e
m2_A = m_nu_A["muon"]
m3_A = m_nu_A["tau"]
dm2_21_A = m2_A**2 - m1_A**2
dm2_32_A = m3_A**2 - m2_A**2
print(f"  m_nu_e = {m1_A:.6f} eV, m_nu_mu = {m2_A:.6f} eV, m_nu_tau = {m3_A:.6f} eV")
print(f"  dm^2_21 = {dm2_21_A:.3e} eV^2 (exp: {dm2_21:.2e})")
print(f"  dm^2_32 = {dm2_32_A:.3e} eV^2 (exp: {dm2_32:.2e})")
print(f"  Problem: m_nu_e = m_e = 511 keV (TOO LARGE by factor {m1_A/m_nu3:.0e})")
print(f"  VERDICT: FAILS (electron level n=0 gives m_e, not a neutrino mass)")
results.append(("A: Direct T-copy", m1_A, m2_A, m3_A, dm2_21_A, dm2_32_A))
print()

# --- Method B: Seesaw with universal M_R from T ---
print("  METHOD B: Type-I seesaw with M_R from T at level n_T")
print("  " + "-" * 50)
print("  m_nu = m_lepton^2 / M_R,  M_R = m_e * phi^(n_T)")
print()

# Scan n_T to find best fit to dm^2_32/dm^2_21 ratio
# m_nu_i = m_lepton_i^2 / (m_e * phi^(n_T))
# dm^2_32 = m_nu_tau^2 - m_nu_mu^2
# dm^2_21 = m_nu_mu^2 - m_nu_e^2

best_n_T = None
best_err = float('inf')
print(f"  {'n_T':>5s} {'M_R (GeV)':>14s} {'m_nu_e':>10s} {'m_nu_mu':>10s} {'m_nu_tau':>10s} {'dm2_21':>12s} {'dm2_32':>12s} {'ratio':>8s}")
for n_T_scan in range(20, 50):
    M_R = m_e_eV * PHI**n_T_scan  # eV
    m_nu_e_B = m_e_eV**2 / M_R
    m_nu_mu_B = m_mu_eV**2 / M_R
    m_nu_tau_B = m_tau_eV**2 / M_R
    dm2_21_B = m_nu_mu_B**2 - m_nu_e_B**2
    dm2_32_B = m_nu_tau_B**2 - m_nu_mu_B**2
    if dm2_21_B > 0 and dm2_32_B > 0:
        ratio_B = dm2_32_B / dm2_21_B
        err_21 = abs(np.log10(dm2_21_B/dm2_21))
        err_32 = abs(np.log10(dm2_32_B/dm2_32))
        total_err = err_21 + err_32
        if total_err < best_err:
            best_err = total_err
            best_n_T = n_T_scan
        if m_nu_tau_B < 10 and m_nu_tau_B > 1e-6:  # reasonable range
            print(f"  {n_T_scan:5d} {M_R/1e9:14.3e} {m_nu_e_B:10.3e} {m_nu_mu_B:10.3e} {m_nu_tau_B:10.3e} {dm2_21_B:12.3e} {dm2_32_B:12.3e} {ratio_B:8.1f}")

print()
if best_n_T:
    M_R = m_e_eV * PHI**best_n_T
    m_nu_e_B = m_e_eV**2 / M_R
    m_nu_mu_B = m_mu_eV**2 / M_R
    m_nu_tau_B = m_tau_eV**2 / M_R
    dm2_21_B = m_nu_mu_B**2 - m_nu_e_B**2
    dm2_32_B = m_nu_tau_B**2 - m_nu_mu_B**2
    print(f"  Best n_T = {best_n_T}: M_R = {M_R/1e9:.3e} GeV")
    print(f"  Ratio dm2_32/dm2_21 (pred) = {dm2_32_B/dm2_21_B:.1f} vs (exp) {ratio_exp:.1f}")
    print(f"  PROBLEM: Seesaw with m_D=m_lepton always gives ratio = (m_tau/m_mu)^4 = {(m_tau_eV/m_mu_eV)**4:.0f}")
    print(f"  This is INDEPENDENT of M_R! The ratio is always wrong.")
    results.append(("B: Universal seesaw", m_nu_e_B, m_nu_mu_B, m_nu_tau_B, dm2_21_B, dm2_32_B))
print()

# --- Method C: Seesaw with level-dependent M_R ---
print("  METHOD C: Level-dependent seesaw: M_R(n) = m_e * phi^(n_T + n)")
print("  " + "-" * 50)
print("  m_nu_i = m_lepton_i^2 / M_R_i where M_R_i scales with phi^n")
print()

# If M_R_i = M_0 * phi^(alpha * n_i), then
# m_nu_i = m_e^2 * phi^(2*n_i) / (M_0 * phi^(alpha*n_i)) = (m_e^2/M_0) * phi^((2-alpha)*n_i)
# For this to give correct dm^2_32/dm^2_21, we need:
# (m_nu3^2 - m_nu2^2) / (m_nu2^2 - m_nu1^2) = 32.6
# With m_nu_i = C * phi^((2-alpha)*n_i):
# m_nu_e = C * phi^0 = C
# m_nu_mu = C * phi^(11*(2-alpha))
# m_nu_tau = C * phi^(17*(2-alpha))
# For ratio: need phi^(2*17*(2-a)) / phi^(2*11*(2-a)) ~ 32.6
# phi^(2*(17-11)*(2-a)) = phi^(12*(2-a)) = 32.6
# 12*(2-a) = log(32.6)/log(phi) = 7.22
# 2-a = 0.602
# a = 1.398

alpha_seesaw = 2 - np.log(ratio_exp) / (2 * (17-11) * np.log(PHI))
# More carefully: m_nu3^2/m_nu2^2 (for m1~0) ~ ratio_exp + 1
# ratio_exp ~ 32.6
# m3^2/m2^2 = 1 + dm2_32/m2^2... this gets complicated
# Let's just try: m_nu_i = m_0 * phi^(beta * n_i)
# Want m3/m2 = phi^(beta*(17-11)) = phi^(6*beta) = m_nu3/m_nu2 = 5.71
# 6*beta = log(5.71)/log(phi) = 3.62
# beta = 0.603

beta_fit = np.log(m_nu3 / m_nu2) / (6 * np.log(PHI))
print(f"  Need m_nu = m_0 * phi^(beta * n) with beta = {beta_fit:.4f}")
print(f"  This requires seesaw exponent alpha = 2 - beta = {2 - beta_fit:.4f}")
print()

# m_0 from m_nu2 = m_0 * phi^(beta*11)
m_0 = m_nu2 / PHI**(beta_fit * 11)
m_nu_e_C = m_0 * PHI**(beta_fit * 0)
m_nu_mu_C = m_0 * PHI**(beta_fit * 11)
m_nu_tau_C = m_0 * PHI**(beta_fit * 17)
dm2_21_C = m_nu_mu_C**2 - m_nu_e_C**2
dm2_32_C = m_nu_tau_C**2 - m_nu_mu_C**2

print(f"  m_0 = {m_0:.5f} eV (fitted from m_nu2)")
print(f"  m_nu_e  = {m_nu_e_C:.5f} eV")
print(f"  m_nu_mu = {m_nu_mu_C:.5f} eV")
print(f"  m_nu_tau = {m_nu_tau_C:.5f} eV")
print(f"  dm^2_21 = {dm2_21_C:.3e} vs {dm2_21:.2e} (err {abs(dm2_21_C-dm2_21)/dm2_21*100:.1f}%)")
print(f"  dm^2_32 = {dm2_32_C:.3e} vs {dm2_32:.2e} (err {abs(dm2_32_C-dm2_32)/dm2_32*100:.1f}%)")
print(f"  beta = {beta_fit:.4f} ~ phi^? : phi^(-1) = {PHI**(-1):.4f}, phi^(-2) = {PHI**(-2):.4f}")
print(f"  VERDICT: beta = 0.60 ~ 1/phi^{np.log(beta_fit)/np.log(1/PHI):.2f}. NOT clean integer.")
print(f"  This is 2-parameter FITTING (m_0, beta), not a derivation.")
results.append(("C: Level-dependent seesaw", m_nu_e_C, m_nu_mu_C, m_nu_tau_C, dm2_21_C, dm2_32_C))
print()

# --- Method D: DSI level offset ---
print("  METHOD D: DSI level offset: m_nu = m_e * phi^(n - k)")
print("  " + "-" * 50)
print()

# Scan k to find which offset gives best fit
print(f"  {'k':>4s} {'m_nu_e (eV)':>12s} {'m_nu_mu (eV)':>12s} {'m_nu_tau (eV)':>12s} {'dm2_21':>12s} {'dm2_32':>12s} {'ratio':>8s}")
for k in range(25, 45):
    m_ne = m_e_eV * PHI**(0 - k)
    m_nmu = m_e_eV * PHI**(11 - k)
    m_ntau = m_e_eV * PHI**(17 - k)
    dm21 = m_nmu**2 - m_ne**2
    dm32 = m_ntau**2 - m_nmu**2
    if dm21 > 0 and dm32 > 0:
        if m_ntau < 1 and m_ntau > 0.001:  # reasonable range
            print(f"  {k:4d} {m_ne:12.5e} {m_nmu:12.5e} {m_ntau:12.5e} {dm21:12.3e} {dm32:12.3e} {dm32/dm21:8.1f}")

print()
print("  The ratio dm^2_32/dm^2_21 with this method is ALWAYS:")
print("  (phi^(2*17) - phi^(2*11)) / (phi^(2*11) - phi^(2*0)) = constant")
ratio_D = (PHI**(2*17) - PHI**(2*11)) / (PHI**(2*11) - 1)
print(f"  = (phi^34 - phi^22) / (phi^22 - 1) = {ratio_D:.2f}")
print(f"  vs experiment: {ratio_exp:.2f}")
print(f"  Error: {abs(ratio_D - ratio_exp)/ratio_exp*100:.1f}%")
print(f"  VERDICT: Wrong ratio by factor {ratio_D/ratio_exp:.0f}x. FAILS.")
print()

# --- Method E: Intersection number scaled ---
print("  METHOD E: m_nu = m_e * phi^(-5a - 6b) = m_e * phi^(-n)")
print("  " + "-" * 50)
print("  Same as Method A with opposite sign. FAILS for same reason.")
print()

# --- Method F: Galois seesaw (key idea) ---
print("  METHOD F: Galois seesaw: m_nu = m_e * phi^n * m_e * phi'^n / M_E8")
print("  " + "-" * 50)
print()
print("  In E8, S and T are coupled. A natural mass formula is:")
print("  m_nu = sqrt(m_S * m_T) = sqrt(m_e * phi^n * m_e * phi'^n)")
print("       = m_e * sqrt(phi^n * phi'^n)")
print("       = m_e * sqrt((phi*phi')^n)")
print("  But phi * phi' = -1, so (phi*phi')^n = (-1)^n")
print("  |m_nu| = m_e * 1 = m_e  (FOR ALL LEVELS!)")
print()
print("  This gives m_nu = m_e for all generations. FAILS trivially.")
print()

# Try: m_nu = m_S * m_T / m_e (product divided by scale)
print("  Variant F2: m_nu = (m_e * phi^n) * (m_e * phi^(-n)) / m_P")
m_nu_F2_e = m_e_eV**2 / (m_P_GeV * 1e9)
print(f"  m_nu = m_e^2 / m_P = {m_nu_F2_e:.5e} eV")
print(f"  This gives {m_nu_F2_e:.2e} eV = same for all (level-independent). FAILS.")
print()

# --- Method G: Algebraic norm seesaw ---
print("  METHOD G: Z[phi] algebraic norm seesaw")
print("  " + "-" * 50)
print()
print("  For fermion at level (a,b), the Z[phi] norm is |N(a+b*phi)| = |a^2 + ab - b^2|")
print("  This is the product of Galois conjugates: (a+b*phi)(a+b*phi')")
print()
print("  Idea: m_nu = m_e * phi^n / f(norm)")
print()

for name, (a, b), n in fermion_data[:3]:  # leptons only
    alg_norm = abs(a**2 + a*b - b**2)
    print(f"  {name:10s}: (a,b)=({a},{b}), n={n}, |N|={alg_norm}")

print()
print("  electron: |N|=0, muon/strange: |N|=1, tau: |N|=1")
print("  The algebraic norm doesn't distinguish muon from tau. FAILS.")
print()

# --- Method H: Seesaw with M_R from E8 branching ---
print("  METHOD H: Seesaw with M_R derived from E8 structure")
print("  " + "-" * 50)
print()

# From exp116's best result: M_R = m_P * phi^(-9) * alpha^3
M_R_H = m_P_GeV * PHI**(-9) * ALPHA**3  # GeV
M_R_H_eV = M_R_H * 1e9
print(f"  Using exp116 best: M_R = m_P * phi^(-9) * alpha^3 = {M_R_H:.3e} GeV")
for name in ["electron", "muon", "tau"]:
    m_D = lepton_masses_eV[name]
    m_nu_H = m_D**2 / M_R_H_eV
    print(f"    m_nu({name}) = {m_nu_H:.5e} eV")

m_nu_e_H = m_e_eV**2 / M_R_H_eV
m_nu_mu_H = m_mu_eV**2 / M_R_H_eV
m_nu_tau_H = m_tau_eV**2 / M_R_H_eV
dm2_21_H = m_nu_mu_H**2 - m_nu_e_H**2
dm2_32_H = m_nu_tau_H**2 - m_nu_mu_H**2
print(f"  dm^2_21 = {dm2_21_H:.3e} vs {dm2_21:.2e}")
print(f"  dm^2_32 = {dm2_32_H:.3e} vs {dm2_32:.2e}")
print(f"  Ratio = {dm2_32_H/dm2_21_H:.1f} vs {ratio_exp:.1f}")
print(f"  SAME hierarchy problem: ratio = (m_tau/m_mu)^4 = {(m_tau_eV/m_mu_eV)**4:.0f}")
results.append(("H: Seesaw M_R=m_P*phi^-9*alpha^3", m_nu_e_H, m_nu_mu_H, m_nu_tau_H, dm2_21_H, dm2_32_H))
print()

# ============================================================
# SECTION 7: The key problem and possible solutions
# ============================================================
print("-" * 80)
print("SECTION 7: Analysis of the hierarchy problem")
print("-" * 80)
print()

print("  THE FUNDAMENTAL PROBLEM:")
print("  ========================")
print()
print("  For Type-I seesaw with m_D = m_lepton:")
print("    m_nu_i = m_lepton_i^2 / M_R")
print("    dm^2_32/dm^2_21 = (m_tau^4 - m_mu^4) / (m_mu^4 - m_e^4)")
print(f"                    = {(m_tau_eV**4 - m_mu_eV**4) / (m_mu_eV**4 - m_e_eV**4):.2f}")
print(f"    Experiment:     = {ratio_exp:.2f}")
print()
print("  The seesaw with m_D=m_lepton gives MUCH too steep hierarchy.")
print("  This is because m_tau/m_mu ~ 17 is already large, and seesaw squares it.")
print()
print("  For the DSI ladder m_nu = m_0 * phi^(beta*n):")
print("    dm^2_32/dm^2_21 = (phi^(34*beta) - phi^(22*beta)) / (phi^(22*beta) - 1)")
print("    For beta=1: ratio = {:.0f} (way too large)".format(ratio_D))
print()

# What beta gives the right ratio?
from scipy.optimize import brentq

def ratio_func(beta):
    if beta <= 0:
        return -ratio_exp
    r = (PHI**(34*beta) - PHI**(22*beta)) / (PHI**(22*beta) - 1)
    return r - ratio_exp

try:
    beta_sol = brentq(ratio_func, 0.01, 0.99)
    print(f"  Solution: beta = {beta_sol:.6f}")
    print(f"  phi^(6*beta) = {PHI**(6*beta_sol):.4f} = m3/m2 ratio")
    print()

    # Is beta a clean number?
    print("  Is beta a clean number?")
    for label, val in [("1/phi", 1/PHI), ("1/phi^2", 1/PHI**2), ("1/3", 1/3),
                        ("1/2", 0.5), ("2/5", 0.4), ("3/5", 0.6),
                        ("phi-1", PHI-1), ("2-phi", 2-PHI),
                        ("1/e", 1/np.e), ("sin^2(tW)", 0.2312),
                        ("alpha*pi", ALPHA*np.pi), ("sqrt(alpha)", np.sqrt(ALPHA)),
                        ("1/(2*phi)", 1/(2*PHI))]:
        err = abs(val - beta_sol) / beta_sol * 100
        if err < 20:
            print(f"    beta ~ {label} = {val:.6f} (err {err:.1f}%)")

    # Check: is 6*beta clean?
    six_beta = 6 * beta_sol
    print(f"\n  6*beta = {six_beta:.4f}")
    for label, val in [("phi^2-1", PHI**2-1), ("2", 2.0), ("3", 3.0), ("phi^2", PHI**2),
                        ("18/5", 3.6), ("7/2", 3.5), ("11/3", 11/3),
                        ("ln(32.6)/ln(phi)", np.log(ratio_exp)/np.log(PHI)),
                        ("phi+1", PHI+1), ("2*phi", 2*PHI)]:
        err = abs(val - six_beta) / six_beta * 100
        if err < 20:
            print(f"    6*beta ~ {label} = {val:.4f} (err {err:.1f}%)")
except Exception as e:
    print(f"  Root finding failed: {e}")

print()

# ============================================================
# SECTION 8: The number 3 from E8
# ============================================================
print("-" * 80)
print("SECTION 8: The number 3 from E8 branching")
print("-" * 80)
print()

print("  E8 branching 248 = 2*(1,120) + (1,8):")
print("    - 2 copies of 120 = S and T = charged + neutral fermions?")
print("    - 8 remaining = gauge sector?")
print()
print("  Under D4 x D4: 248 = (28,1) + (1,28) + (8v,8v) + (8s,8s) + (8c,8c)")
print("  The 3 terms (8v,8v), (8s,8s), (8c,8c) are the triality triplet.")
print("  Could these map to 3 generations?")
print()

# Count: each (8,8) has 64 elements.
# Total: 28 + 28 + 3*64 = 248. Check: 56 + 192 = 248. Yes.
print("  Dimension count: 28 + 28 + 3*64 = 248. Check: {}.".format(28+28+3*64))
print()

# In the 600-cell:
# Type A (8 vertices): these form the D4 vector rep (8v)
# Type B spinor (8): these form 8s
# Type B co-spinor (8): these form 8c
# Type C (96): these are the remaining 96 snub-24-cell vertices

print("  In S (600-cell):")
print("    Type A (8v): 8 vertices -> D4 vector representation")
print("    Type B_s (8s): 8 vertices -> D4 spinor representation")
print("    Type B_c (8c): 8 vertices -> D4 co-spinor representation")
print("    Type C: 96 vertices -> snub 24-cell (fermion sector)")
print()

# The cross-coupling: do Type A_S couple to Type A_T differently from Type B_S to Type B_T?
type_A_T, type_B_T, type_C_T = [], [], []
for i, v in enumerate(verts_S * PHI_CONJ):
    nz = sum(1 for x in v if abs(x) > 0.01)
    # T vertices are scaled by phi', so check different thresholds
    # Type A in T: (+-phi', 0, 0, 0) etc -> abs ~ 0.618
    if nz == 1 and any(abs(abs(x) - abs(PHI_CONJ)) < 0.01 for x in v):
        type_A_T.append(i)
    elif nz == 4 and all(abs(abs(x) - 0.5*abs(PHI_CONJ)) < 0.01 for x in v):
        type_B_T.append(i)
    else:
        type_C_T.append(i)

print(f"  T vertex types: A={len(type_A_T)}, B={len(type_B_T)}, C={len(type_C_T)}")

# But actually: the TYPE classification is about the S-vertex that generates each T-vertex
# T_i = phi' * S_i, so if S_i is type A, T_i is type A_T
# The type depends on which S vertex generated it, not on the coordinates of T
# So: type_A_T should have the SAME indices as type_A_S
print("  (T types mirror S types since T_i = phi'*S_i)")
print(f"  type_A_T indices = type_A_S indices: {sorted(type_A_T) == sorted(type_A_S)}")
print()

# Cross-coupling matrix between types
print("  Cross-coupling (S-T) between types (mean dot products in E8):")
for s_label, s_indices in [("A_S", type_A_S), ("B_S", type_B_S), ("C_S", type_C_S)]:
    for t_label, t_indices in [("A_T", type_A_S), ("B_T", type_B_S), ("C_T", type_C_S)]:
        # Cross dots for these subsets
        subdots = []
        for si in s_indices:
            for ti in t_indices:
                subdots.append(cross_dots[si, ti])
        subdots = np.array(subdots)
        n1_count = np.sum(np.abs(subdots - 1.0) < 0.01)
        print(f"    {s_label}-{t_label}: mean_dot={np.mean(subdots):.4f}, dot=1 pairs: {n1_count}")

print()

# ============================================================
# SECTION 9: Alternative - neutrinos from DIFFERENT levels
# ============================================================
print("-" * 80)
print("SECTION 9: Neutrinos from different DSI levels (not tied to leptons)")
print("-" * 80)
print()

print("  What if neutrino levels are NOT the same as charged lepton levels?")
print("  Scan: m_nu_i = m_e * phi^(-k_i) for k_i to match oscillation data.")
print()

# We need two numbers k3, k2 such that:
# m_nu3 = m_e * phi^(-k3) and m_nu2 = m_e * phi^(-k2)
k3 = -np.log(m_nu3 / m_e_eV) / np.log(PHI)
k2 = -np.log(m_nu2 / m_e_eV) / np.log(PHI)
print(f"  m_nu3 = {m_nu3:.5f} eV => k3 = {k3:.4f}")
print(f"  m_nu2 = {m_nu2:.5f} eV => k2 = {k2:.4f}")
print(f"  k3 - k2 = {k3 - k2:.4f}")
print()

# Check if k3, k2 are close to integers or simple fractions
for label, k_val in [("k3", k3), ("k2", k2)]:
    print(f"  {label} = {k_val:.4f}")
    nearest_int = round(k_val)
    print(f"    Nearest integer: {nearest_int} (err {abs(k_val - nearest_int):.4f})")
    nearest_half = round(2*k_val)/2
    print(f"    Nearest half-int: {nearest_half} (err {abs(k_val - nearest_half):.4f})")
    # What DSI level is this?
    print(f"    m_e * phi^(-{nearest_int}) = {m_e_eV * PHI**(-nearest_int):.5f} eV")

print()

# Try integer levels in the correct range
print("  Integer level scan: m_nu = m_e * phi^(-k)")
print(f"  {'k':>4s} {'m_nu (eV)':>12s} {'vs m_nu3':>10s} {'vs m_nu2':>10s}")
for k in range(30, 42):
    m = m_e_eV * PHI**(-k)
    err3 = abs(m - m_nu3) / m_nu3 * 100
    err2 = abs(m - m_nu2) / m_nu2 * 100
    best = min(err3, err2)
    marker = " <--" if best < 25 else ""
    print(f"  {k:4d} {m:12.8f} {err3:9.1f}% {err2:9.1f}%{marker}")

print()

m_test_34 = m_e_eV * PHI**(-34)
m_test_37 = m_e_eV * PHI**(-37)
print(f"  phi^(-34): m = {m_test_34:.6f} eV (vs m_nu3 = {m_nu3:.5f} eV, err {abs(m_test_34-m_nu3)/m_nu3*100:.1f}%)")
print(f"  phi^(-37): m = {m_test_37:.6f} eV (vs m_nu2 = {m_nu2:.5f} eV, err {abs(m_test_37-m_nu2)/m_nu2*100:.1f}%)")
print()

# Check: what if neutrino levels are n_nu = n_charged + some_offset from T?
# Charged leptons: n = 0, 11, 17
# If neutrinos: n_nu = n_charged - K for constant K
# Then m_nu_i = m_e * phi^(n_i - K) = m_charged_i * phi^(-K)
# dm^2_32/dm^2_21 = same as DSI offset. Same problem.
print("  ANY mechanism m_nu = f(m_charged) with SAME function f gives wrong ratio")
print("  unless f has generation-dependent structure.")
print()

# ============================================================
# SECTION 10: Most promising: independent neutrino levels from T
# ============================================================
print("-" * 80)
print("SECTION 10: Independent neutrino levels from T-copy")
print("-" * 80)
print()

print("  In S, charged fermions occupy levels {0, 3, 5, 11, 16, 17, 19, 26}")
print("  The T-copy provides a SECOND set of levels.")
print()
print("  If neutrino levels in T are INDEPENDENT, what (a,b)_nu give correct masses?")
print("  m_nu = m_e * phi^(-n_nu)")
print("  (negative exponent because T is the conjugate/mirror)")
print()

# Best integer n_nu values from Section 9
n_nu3_int = round(k3)  # k3 ~ 33.56 -> 34
n_nu2_int = round(k2)  # k2 ~ 37.18 -> 37

print(f"  Best n_nu for m_nu3: {k3:.2f} -> try n=33 or n=34")
print(f"  Best n_nu for m_nu2: {k2:.2f} -> try n=37 or n=38")
print()

# Check n=33 and n=34 for nu3
for n_test in [33, 34, 35]:
    m_pred = m_e_eV * PHI**(-n_test)
    err = abs(m_pred - m_nu3) / m_nu3 * 100
    # Find (a,b) minimizing |a|+|b| with 5a+6b=n_test
    best_ab = None
    best_complexity = float('inf')
    for a_scan in range(-20, 21):
        if (n_test - 5*a_scan) % 6 == 0:
            b_scan = (n_test - 5*a_scan) // 6
            complexity = abs(a_scan) + abs(b_scan)
            if complexity < best_complexity:
                best_complexity = complexity
                best_ab = (a_scan, b_scan)
    print(f"  n={n_test}: m = {m_pred:.6f} eV, err = {err:.1f}%, (a,b) = {best_ab}, |a|+|b| = {best_complexity}")

print()

# Check n=37 and n=38 for nu2
for n_test in [36, 37, 38]:
    m_pred = m_e_eV * PHI**(-n_test)
    err = abs(m_pred - m_nu2) / m_nu2 * 100
    best_ab = None
    best_complexity = float('inf')
    for a_scan in range(-20, 21):
        if (n_test - 5*a_scan) % 6 == 0:
            b_scan = (n_test - 5*a_scan) // 6
            complexity = abs(a_scan) + abs(b_scan)
            if complexity < best_complexity:
                best_complexity = complexity
                best_ab = (a_scan, b_scan)
    print(f"  n={n_test}: m = {m_pred:.6f} eV, err = {err:.1f}%, (a,b) = {best_ab}, |a|+|b| = {best_complexity}")

print()

# The best integer pair: let's check (n_nu3, n_nu2) for dm^2 predictions
print("  Best integer combinations (n_nu3, n_nu2) for oscillation data:")
print(f"  {'n3':>4s} {'n2':>4s} {'m_nu3 (eV)':>12s} {'m_nu2 (eV)':>12s} {'dm2_21':>12s} {'err_21':>8s} {'dm2_32':>12s} {'err_32':>8s}")
for n3_t in [33, 34, 35]:
    for n2_t in [36, 37, 38]:
        m3_t = m_e_eV * PHI**(-n3_t)
        m2_t = m_e_eV * PHI**(-n2_t)
        if m3_t > m2_t:
            dm21_t = m2_t**2  # assuming m1~0
            dm32_t = m3_t**2 - m2_t**2
            err21 = abs(dm21_t - dm2_21) / dm2_21 * 100
            err32 = abs(dm32_t - dm2_32) / dm2_32 * 100
            ratio_t = dm32_t / dm21_t
            print(f"  {n3_t:4d} {n2_t:4d} {m3_t:12.6f} {m2_t:12.6f} {dm21_t:12.3e} {err21:7.1f}% {dm32_t:12.3e} {err32:7.1f}%  ratio={ratio_t:.1f}")

print()

# ============================================================
# SECTION 11: Remarkable observation - n_nu3 = 19 = n_bottom!
# ============================================================
print("-" * 80)
print("SECTION 11: Level coincidences and patterns")
print("-" * 80)
print()

print("  Charged fermion occupied levels: {0, 3, 5, 11, 16, 17, 19, 26}")
print(f"  Required neutrino levels: n_nu3 ~ {k3:.1f}, n_nu2 ~ {k2:.1f}")
print()

# Best integer candidates
for n_nu, nu_name, nu_exp in [(34, "nu3", m_nu3), (37, "nu2", m_nu2)]:
    m_pred = m_e_eV * PHI**(-n_nu)
    err = abs(m_pred - nu_exp) / nu_exp * 100
    print(f"  m_{nu_name} at n={n_nu}: m_e * phi^(-{n_nu}) = {m_pred:.6f} eV (exp {nu_exp:.5f} eV, err {err:.1f}%)")

print()

# Check: with best integer levels
n_best3 = 34
n_best2 = 37
m_nu3_pred = m_e_eV * PHI**(-n_best3)
m_nu2_pred = m_e_eV * PHI**(-n_best2)
dm2_32_pred = m_nu3_pred**2 - m_nu2_pred**2
dm2_21_pred = m_nu2_pred**2  # assuming m1~0
print(f"  With (n_nu3, n_nu2) = ({n_best3}, {n_best2}):")
print(f"    m_nu3 = {m_nu3_pred:.6f} eV (exp {m_nu3:.5f} eV, err {abs(m_nu3_pred-m_nu3)/m_nu3*100:.1f}%)")
print(f"    m_nu2 = {m_nu2_pred:.6f} eV (exp {m_nu2:.5f} eV, err {abs(m_nu2_pred-m_nu2)/m_nu2*100:.1f}%)")
print(f"    dm^2_32 = {dm2_32_pred:.3e} vs {dm2_32:.3e} (err {abs(dm2_32_pred-dm2_32)/dm2_32*100:.1f}%)")
print(f"    dm^2_21 = {dm2_21_pred:.3e} vs {dm2_21:.3e} (err {abs(dm2_21_pred-dm2_21)/dm2_21*100:.1f}%)")
print(f"    ratio = {dm2_32_pred/dm2_21_pred:.1f} vs {ratio_exp:.1f}")
print()

# What about n_nu1?
# If hierarchy continues: n_nu1 > n_nu2 > n_nu3 (higher n = lighter)
print("  What n_nu1 (lightest neutrino)?")
for n1_t in [38, 39, 40, 42, 45, 50]:
    m1_t = m_e_eV * PHI**(-n1_t)
    print(f"  n_nu1 = {n1_t}: m_nu1 = {m1_t:.6e} eV")

print()

# Level relationships
print("  Level relationships:")
print(f"  n_nu3 = {n_best3} = 2 * 17 = 2 * n_tau")
print(f"  n_nu2 = {n_best2} = ? (37 is prime)")
print(f"  n_nu3 - n_nu2 = {n_best3 - n_best2} (separation)")
print(f"  Charged lepton levels: 0, 11, 17. Differences: 11, 6")
print(f"  Neutrino levels: 34, 37. Difference: 3")
print(f"  34 = 2*17 = 2*n_tau. Coincidence?")
print(f"  37 - 26 = 11 = n_muon. Coincidence?")
print()

# ============================================================
# SECTION 12: Half-integer levels (fractional DSI)
# ============================================================
print("-" * 80)
print("SECTION 12: Half-integer and fractional DSI levels")
print("-" * 80)
print()

# k3 ~ 33.56, k2 ~ 37.18 (from Section 9)
# What if k3 = 67/2 = 33.5 exactly?
m_nu3_half = m_e_eV * PHI**(-67/2)
err_half = abs(m_nu3_half - m_nu3) / m_nu3 * 100
print(f"  k3 = 67/2 = 33.5:")
print(f"    m_nu3 = m_e * phi^(-67/2) = {m_nu3_half:.6f} eV (err {err_half:.1f}%)")
print()

# What about k2?
for k_frac in [37, 37.5, 74/2, 75/2]:
    m_test = m_e_eV * PHI**(-k_frac)
    err_test = abs(m_test - m_nu2) / m_nu2 * 100
    print(f"  k2 = {k_frac}: m_nu2 = {m_test:.6f} eV (err {err_test:.1f}%)")

print()

# Check: m3/m2 = phi^(k2 - k3)
print(f"  m3/m2 = phi^(k2-k3) where k2-k3 = {k2-k3:.4f}")
print(f"  phi^(k2-k3) = {PHI**(k2-k3):.4f} vs m3/m2 = {m_nu3/m_nu2:.4f}")
for label, val in [("18/5", 18/5), ("7/2", 7/2), ("phi^2", PHI**2),
                    ("11/3", 11/3), ("4", 4), ("3", 3),
                    ("phi+2", PHI+2), ("2*phi", 2*PHI)]:
    err_delta = abs(val - (k2-k3)) / (k2-k3) * 100
    if err_delta < 30:
        print(f"    k2-k3 ~ {label} = {val:.4f} (err {err_delta:.1f}%)")

print()

# ============================================================
# SECTION 13: Corrections from sector-dependent delta_d, delta_k
# ============================================================
print("-" * 80)
print("SECTION 13: Neutrino sector corrections")
print("-" * 80)
print()

print("  Charged leptons use: d_eff = 5 + sin^2(tW) = 5.231")
print("                       k_eff = 6 - 1/phi^4 = 5.854")
print()
print("  For neutrinos in T, corrections could DIFFER:")
print("  - No electromagnetic charge -> delta_d_EM = 0?")
print("  - Weak interaction only -> delta from sin^2(tW)?")
print()

sin2tw = 0.2312
# If neutrinos have delta_d = 0 (no EM charge) and delta_k = 0 (no QCD)
# Then m_nu = m_e * phi^(5*a + 6*b) with the SAME a,b
# This doesn't help since a,b are the same

# But if the NEUTRINO sector has its own intersection numbers from T:
# From exp117: 120-cell has a_1=0, b_1=3
# Could T have different effective intersection numbers?
# Actually T has SAME graph as S (just scaled), so same a_1=5, b_1=6.
# But the COUPLING between S and T could define different intersection numbers.

print("  T-copy has SAME intersection numbers as S (same graph).")
print("  Cross-coupling S-T may define different effective parameters.")
print()

# From Section 5, we found cross_neighbors at dot=1
# Let's check if the cross-coupling gives a DIFFERENT intersection number
print("  Cross-intersection numbers (from S-T E8 coupling):")
print(f"  Each S-vertex has {unique_cn} T-neighbors at dot=1")

# The effective intersection number for S-T coupling:
# a_1^cross = number of T-vertices at E8-distance 1 from any S-vertex
# From our data: should be a constant
if len(unique_cn) == 1:
    a1_cross = unique_cn[0]
    print(f"  a_1^cross = {a1_cross}")

    # If neutrinos use a_1_cross instead of a_1=5:
    # n_nu = a_1_cross * a + b_1 * b (assuming b_1 stays the same?)
    # But what b_1_cross?
    # The second shell T-neighbors would tell us

    # Count second-shell T-neighbors (dot = 0 at E8 convention might be d=2)
    # In E8, distances go by dot products: 2, 1, 0, -1, -2
    # dot=1 is distance 1 (neighbors)
    # dot=0 is distance 2
    # dot=-1 is distance 3
    # dot=-2 is distance 4 (antipodal)

    b1_cross_candidates = unique_cn_0[0] if len(unique_cn_0) == 1 else "variable"
    print(f"  b_1^cross (T at dot=0) = {b1_cross_candidates}")

    if isinstance(b1_cross_candidates, (int, np.integer)):
        print(f"\n  If neutrino levels use cross-coupling: n_nu = {a1_cross}*a + {b1_cross_candidates}*b")
        # What neutrino (a,b) give right masses?
        for nu_name, n_target in [("nu3", round(k3)), ("nu2", round(k2))]:
            print(f"\n  {nu_name} (n_target ~ {n_target}):")
            best_ab = None
            best_complexity = float('inf')
            for a_scan in range(-10, 11):
                rem = n_target - a1_cross * a_scan
                if b1_cross_candidates != 0 and rem % b1_cross_candidates == 0:
                    b_scan = rem // b1_cross_candidates
                    complexity = abs(a_scan) + abs(b_scan)
                    if complexity < best_complexity:
                        best_complexity = complexity
                        best_ab = (a_scan, b_scan)
            if best_ab:
                print(f"    Min complexity (a,b) = {best_ab}, |a|+|b| = {best_complexity}")

print()

# ============================================================
# SECTION 14: Summary Table
# ============================================================
print("=" * 80)
print("SECTION 14: SUMMARY")
print("=" * 80)
print()

print("  RESULTS OF ALL APPROACHES:")
print("  " + "=" * 76)
print(f"  {'Method':40s} {'Status':12s} {'Key Issue'}")
print("  " + "-" * 76)
print(f"  {'A: Direct T-copy m_e*phi^(-n)':40s} {'FAILS':12s} {'m_nu_e = m_e (n=0 trivial)'}")
print(f"  {'B: Universal seesaw':40s} {'FAILS':12s} {'Ratio = (m_tau/m_mu)^4 = 8e4 (too steep)'}")
print(f"  {'C: Level-dependent seesaw':40s} {'FITTING':12s} {'beta=0.60 not clean, 2 free params'}")
ratio_D_str = "Ratio {:.0f} (too steep)".format(ratio_D)
print(f"  {'D: DSI level offset':40s} {'FAILS':12s} {ratio_D_str}")
print(f"  {'E: Mirror formula':40s} {'FAILS':12s} {'Same as A'}")
galois_msg = "m_nu = m_e for all (phi*phi'=-1)"
print(f"  {'F: Galois seesaw sqrt(S*T)':40s} {'FAILS':12s} {galois_msg}")
print(f"  {'G: Z[phi] algebraic norm':40s} {'FAILS':12s} {'Norm 0 or 1 for leptons, no separation'}")
print(f"  {'H: Seesaw M_R from exp116':40s} {'FAILS':12s} {'Same ratio problem as B'}")
print("  " + "-" * 76)
print()

print("  PARTIAL SUCCESSES:")
print("  - m_e * phi^(-34) = {:.6f} eV vs m_nu3 ~ {:.5f} eV (err {:.1f}%)".format(
    m_e_eV * PHI**(-34), m_nu3, abs(m_e_eV * PHI**(-34) - m_nu3)/m_nu3*100))
print("  - m_e * phi^(-37) = {:.6f} eV vs m_nu2 ~ {:.5f} eV (err {:.1f}%)".format(
    m_e_eV * PHI**(-37), m_nu2, abs(m_e_eV * PHI**(-37) - m_nu2)/m_nu2*100))
print("  - n_nu3 ~ {:.2f} (not clean integer, nearest 34)".format(k3))
print("  - n_nu2 ~ {:.2f} (not clean integer, nearest 37)".format(k2))
print("  - 34 = 2*17 = 2*n_tau. Suggestive but NOT DERIVED.")
print("  - beta ~ 3/5 or 1/phi for level-dependent seesaw (FITTING)")
print()

print("  STRUCTURAL FINDINGS:")
print("  1. T in R^4 is a SCALED REFLECTED 600-cell (radius 1/phi)")
print("     Same graph, same spectrum, same adjacency as S.")
print("  2. E8 cross-coupling: each S-vertex has {} T-neighbors at dot=1".format(
    unique_cn[0] if len(unique_cn) == 1 else unique_cn))
print("  3. The DSI potential is INVARIANT under phi -> 1/phi (no new levels)")
print("  4. Galois seesaw trivializes: phi*phi'=-1 kills level dependence")
print("  5. ALL Type-I seesaw variants have the WRONG hierarchy ratio")
print("     because (m_tau/m_mu)^4 >> dm^2_32/dm^2_21")
print()

print("  HONEST ASSESSMENT:")
print("  ==================")
print("  The T-copy (second 600-cell) does NOT naturally produce neutrino masses.")
print("  The key obstacles are:")
print("  (a) phi*phi'=-1 means the Galois conjugation kills the level structure")
print("  (b) Any seesaw with m_D=m_lepton gives wrong mass ratios")
print("  (c) DSI potential is phi<->1/phi symmetric, giving no new physics")
print()
print("  The neutrino sector likely requires a DIFFERENT mechanism beyond the")
print("  simple S-T duality in E8. Possible directions:")
print("  - Majorana mass from the (1,8) term in 248 = 2*(1,120) + (1,8)")
print("  - Neutrino masses from LOOP corrections (radiative mechanism)")
print("  - Non-minimal seesaw where m_D != m_lepton")
print("  - Neutrinos as TOPOLOGICAL excitations (not DSI levels)")
print()
print("  CATEGORY: MOSTLY NEGATIVE RESULT")
print("  Notable PATTERN: m_nu3 ~ m_e * phi^(-33.6), m_nu2 ~ m_e * phi^(-37.2)")
print("  Nearest integers (34, 37) give 19% and 9% errors.")
print("  34 = 2*n_tau is suggestive but unconfirmed.")
print("  The mass RATIO m3/m2 has k2-k3 ~ phi+2 = 3.618 (0.04% error)")
print("  or k2-k3 ~ 18/5 (0.5%). Neither has clear geometric meaning.")
print()
print("=" * 80)
print("END EXP-134")
print("=" * 80)
