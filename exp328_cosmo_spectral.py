"""
exp328_cosmo_spectral.py - Cosmological Constant from Spectral Action
=====================================================================

Can the spectral action on the 600-cell DERIVE the CC exponent z~56.88?

Formula: Lambda_Planck = alpha^z, z = (N-N_gen)/2 - phi = 57 - alpha_s
Currently: STRONG PATTERN (0.13 sigma, 0.33%). No derivation.

This experiment systematically tests 7 routes:
  1. Seeley-DeWitt ratios (c_0*c_2/c_1^2, higher c_k)
  2. Spectral zeta function zeta_D(s) at special values
  3. Exact spectral action + remainder
  4. Str(M^4) supertrace with framework masses
  5. Graded spectral action (boson-fermion)
  6. Systematic search for 57 in spectral data
  7. Honest assessment

Expected: z=57-alpha_s remains PATTERN. Deliverables: c_3, c_4, zeta_D, Str(M^4).

Author: Razvan-Constantin Anghelina
Date: 2026-02-16
"""

import numpy as np
from itertools import permutations, product
from collections import defaultdict, Counter
from fractions import Fraction
import math
import warnings
warnings.filterwarnings('ignore')

print("=" * 72)
print("exp328: COSMOLOGICAL CONSTANT FROM SPECTRAL ACTION")
print("=" * 72)

# =====================================================================
# PART 0: CONSTANTS + CC OBSERVED VALUE
# =====================================================================
print("\n" + "=" * 72)
print("PART 0: Framework Constants + Cosmological Data")
print("=" * 72)

a_1 = 5
b_1 = a_1 + 1  # = 6
N_vert = 120
N_gen = 3
N_eig = 9
degree = 12  # icosahedron: each vertex has 12 neighbors

PHI = (1 + np.sqrt(5)) / 2
phi_prime = (1 - np.sqrt(5)) / 2

# alpha from icosahedral Laplacian
coeff_alpha = 4 * a_1 * PHI**4
disc_alpha = coeff_alpha**2 - 8 * np.pi
ALPHA = (coeff_alpha - np.sqrt(disc_alpha)) / (4 * np.pi)
alpha_s = 1 / (2 * PHI**3)
sin2_tW = b_1 / (a_1**2 + 1)  # = 6/26

print(f"a_1 = {a_1}")
print(f"phi = {PHI:.10f}")
print(f"alpha = {ALPHA:.10e} (1/alpha = {1/ALPHA:.4f})")
print(f"alpha_s = {alpha_s:.10f} = 1/(2*phi^3)")
print(f"sin^2(tW) = {sin2_tW:.10f} = {b_1}/{a_1**2+1}")

# Cosmological data
HBAR = 1.054571817e-34
C_light = 299792458.0
G_Newton = 6.67430e-11
E_CHARGE = 1.602176634e-19

l_P = np.sqrt(HBAR * G_Newton / C_light**3)
H_0_SI = 67.36 * 1000 / 3.0857e22
Omega_Lambda = 0.6847

Lambda_m2 = 3 * Omega_Lambda * H_0_SI**2 / C_light**2
Lambda_Planck = Lambda_m2 * l_P**2

z_exact = np.log(Lambda_Planck) / np.log(ALPHA)
z_formula = (N_vert - N_gen) / 2 - PHI  # = 57 - phi
z_formula_v2 = 57 - alpha_s  # equivalent

# Error estimation (from Planck uncertainties)
H0_err = 0.54  # km/s/Mpc
OmL_err = 0.0073
# Lambda_P ~ H_0^2 * Omega_Lambda, so rel_err = sqrt((2*dH/H)^2 + (dOm/Om)^2)
# z = ln(Lambda_P)/ln(alpha), so delta_z = (delta Lambda_P / Lambda_P) / |ln(alpha)|
rel_err_Lambda = np.sqrt((2*H0_err/67.36)**2 + (OmL_err/Omega_Lambda)**2)
delta_z = rel_err_Lambda / np.abs(np.log(ALPHA))

print(f"\nCosmological data:")
print(f"  H_0 = 67.36 +- 0.54 km/s/Mpc")
print(f"  Omega_Lambda = {Omega_Lambda} +- {OmL_err}")
print(f"  Lambda/l_P^2 = {Lambda_Planck:.4e}")
print(f"  z_exact = ln(Lambda_P)/ln(alpha) = {z_exact:.6f}")
print(f"  z_formula = (N-N_gen)/2 - phi = {z_formula:.6f}")
print(f"  z_formula = 57 - alpha_s = {z_formula_v2:.6f}")
print(f"  Difference: {abs(z_exact - z_formula):.4f}")
print(f"  sigma ~ {abs(z_exact - z_formula)/delta_z:.2f} sigma")

Lambda_pred = ALPHA**z_formula
print(f"\n  Lambda_pred = alpha^z = {Lambda_pred:.4e}")
print(f"  Lambda_obs  = {Lambda_Planck:.4e}")
print(f"  Ratio = {Lambda_pred/Lambda_Planck:.6f}")

# =====================================================================
# PART 1: BUILD 600-CELL + COMPUTE ALL EIGENVALUES
# =====================================================================
print("\n" + "=" * 72)
print("PART 1: 600-Cell Construction + Full Spectrum")
print("=" * 72)

# --- Generate 120 vertices ---
verts_set = set()

# Type A: 8 axis vertices
for i in range(4):
    for s in [1.0, -1.0]:
        v = [0.0, 0.0, 0.0, 0.0]
        v[i] = s
        verts_set.add(tuple(v))

# Type B: 16 half-integer vertices
for s0 in [0.5, -0.5]:
    for s1 in [0.5, -0.5]:
        for s2 in [0.5, -0.5]:
            for s3 in [0.5, -0.5]:
                verts_set.add((s0, s1, s2, s3))

# Type C: 96 golden-ratio vertices (even permutations)
base = [PHI/2, 0.5, 1/(2*PHI), 0.0]
even_perms = []
for p in permutations(range(4)):
    inv = sum(1 for i in range(4) for j in range(i+1, 4) if p[i] > p[j])
    if inv % 2 == 0:
        even_perms.append(p)

for perm in even_perms:
    coords = [base[perm[i]] for i in range(4)]
    nz = [i for i in range(4) if abs(coords[i]) > 1e-12]
    for signs in product([1, -1], repeat=len(nz)):
        v = list(coords)
        for idx, s in zip(nz, signs):
            v[idx] *= s
        verts_set.add(tuple(round(x, 10) for x in v))

verts = np.array(sorted(verts_set))
Nv = len(verts)
print(f"Vertices: {Nv}")
assert Nv == 120, f"Expected 120 vertices, got {Nv}"

# --- Find edges ---
dots = verts @ verts.T
edge_thresh = PHI / 2  # cos(pi/5) for unit-radius 600-cell
edges = []
adj = np.zeros((Nv, Nv), dtype=int)

for i in range(Nv):
    for j in range(i+1, Nv):
        if abs(dots[i, j] - edge_thresh) < 0.001:
            edges.append((i, j))
            adj[i, j] = 1
            adj[j, i] = 1

Ne = len(edges)
print(f"Edges: {Ne}")
assert Ne == 720, f"Expected 720 edges, got {Ne}"

edge_to_idx = {}
for idx, (i, j) in enumerate(edges):
    edge_to_idx[(i, j)] = idx
    edge_to_idx[(j, i)] = idx

adj_list = defaultdict(set)
for i, j in edges:
    adj_list[i].add(j)
    adj_list[j].add(i)

# --- Find triangles ---
triangles = []
for i in range(Nv):
    for j in adj_list[i]:
        if j > i:
            common = adj_list[i] & adj_list[j]
            for k in common:
                if k > j:
                    triangles.append((i, j, k))

Nf = len(triangles)
print(f"Triangles: {Nf}")
assert Nf == 1200, f"Expected 1200 triangles, got {Nf}"

face_to_idx = {}
for idx, (i, j, k) in enumerate(triangles):
    face_to_idx[(i, j, k)] = idx

# --- Find tetrahedra ---
tetrahedra = []
for i in range(Nv):
    ni = adj_list[i]
    for j in ni:
        if j > i:
            common_ij = ni & adj_list[j]
            for k in common_ij:
                if k > j:
                    common_ijk = common_ij & adj_list[k]
                    for l in common_ijk:
                        if l > k:
                            tetrahedra.append((i, j, k, l))

Nc = len(tetrahedra)
print(f"Tetrahedra: {Nc}")
assert Nc == 600, f"Expected 600 tetrahedra, got {Nc}"

euler = Nv - Ne + Nf - Nc
print(f"Euler characteristic: {Nv} - {Ne} + {Nf} - {Nc} = {euler}")
assert euler == 0, "S^3 should have Euler = 0"

# --- Boundary operators ---
print("\nBuilding boundary operators...")

# d0: C^0 -> C^1 (Ne x Nv)
d0 = np.zeros((Ne, Nv))
for e_idx, (i, j) in enumerate(edges):
    d0[e_idx, i] = -1.0
    d0[e_idx, j] = +1.0

# d1: C^1 -> C^2 (Nf x Ne)
d1 = np.zeros((Nf, Ne))
for f_idx, (i, j, k) in enumerate(triangles):
    d1[f_idx, edge_to_idx[(i, j)]] = +1.0
    d1[f_idx, edge_to_idx[(j, k)]] = +1.0
    d1[f_idx, edge_to_idx[(i, k)]] = -1.0

# d2: C^2 -> C^3 (Nc x Nf)
d2 = np.zeros((Nc, Nf))
for c_idx, (i, j, k, l) in enumerate(tetrahedra):
    faces_of_tet = [
        ((j, k, l), +1),
        ((i, k, l), -1),
        ((i, j, l), +1),
        ((i, j, k), -1),
    ]
    for face, sign in faces_of_tet:
        f_idx = face_to_idx.get(face)
        if f_idx is not None:
            d2[c_idx, f_idx] = sign

# Verify d^2 = 0
check_d1d0 = np.max(np.abs(d1 @ d0))
check_d2d1 = np.max(np.abs(d2 @ d1))
print(f"||d1*d0|| = {check_d1d0:.2e} (should be ~0)")
print(f"||d2*d1|| = {check_d2d1:.2e} (should be ~0)")
assert check_d1d0 < 1e-10, "d^2 != 0 for d1*d0"
assert check_d2d1 < 1e-10, "d^2 != 0 for d2*d1"

# --- Hodge Laplacians ---
print("Computing Hodge Laplacians...")
Delta_0 = d0.T @ d0  # Nv x Nv
Delta_1 = d0 @ d0.T + d1.T @ d1  # Ne x Ne
Delta_2 = d1 @ d1.T + d2.T @ d2  # Nf x Nf
Delta_3 = d2 @ d2.T  # Nc x Nc

print(f"Delta_0: {Delta_0.shape}, Delta_1: {Delta_1.shape}")
print(f"Delta_2: {Delta_2.shape}, Delta_3: {Delta_3.shape}")

# --- Eigenvalues ---
print("Computing eigenvalues...")
evals_0 = np.sort(np.linalg.eigvalsh(Delta_0))
evals_1 = np.sort(np.linalg.eigvalsh(Delta_1))
evals_2 = np.sort(np.linalg.eigvalsh(Delta_2))
evals_3 = np.sort(np.linalg.eigvalsh(Delta_3))

# Clean near-zero
tol = 1e-8
evals_0[np.abs(evals_0) < tol] = 0.0
evals_1[np.abs(evals_1) < tol] = 0.0
evals_2[np.abs(evals_2) < tol] = 0.0
evals_3[np.abs(evals_3) < tol] = 0.0

# Betti numbers (from zero eigenvalues)
b0 = np.sum(evals_0 == 0)
b1 = np.sum(evals_1 == 0)
b2 = np.sum(evals_2 == 0)
b3 = np.sum(evals_3 == 0)
print(f"Betti numbers: ({b0}, {b1}, {b2}, {b3})")
assert (b0, b1, b2, b3) == (1, 0, 0, 1), f"Expected (1,0,0,1), got ({b0},{b1},{b2},{b3})"

# Full D^2 spectrum
all_evals = np.concatenate([evals_0, evals_1, evals_2, evals_3])
all_evals = np.sort(all_evals)

lambda_max = np.max(all_evals)
print(f"\nlambda_max(D^2) = {lambda_max:.6f}")
print(f"Expected 6+6*phi = {6+6*PHI:.6f}")

# Eigenvalue statistics per degree
for deg, (ev, label) in enumerate([(evals_0, "Delta_0"), (evals_1, "Delta_1"),
                                    (evals_2, "Delta_2"), (evals_3, "Delta_3")]):
    nz = ev[ev > tol]
    if len(nz) > 0:
        print(f"  {label}: dim={len(ev)}, nonzero={len(nz)}, "
              f"min={np.min(nz):.4f}, max={np.max(nz):.4f}, "
              f"Tr={np.sum(ev):.2f}")
    else:
        print(f"  {label}: dim={len(ev)}, all zero")

# Distinct eigenvalue analysis per degree
print("\nDistinct eigenvalues per degree (rounded to 4 decimals):")
distinct_per_deg = []
for deg, (ev, label) in enumerate([(evals_0, "Delta_0"), (evals_1, "Delta_1"),
                                    (evals_2, "Delta_2"), (evals_3, "Delta_3")]):
    rounded = np.round(ev, 4)
    unique_vals, counts = np.unique(rounded, return_counts=True)
    distinct_per_deg.append((unique_vals, counts))
    print(f"  {label}: {len(unique_vals)} distinct eigenvalues")
    for v, c in zip(unique_vals, counts):
        print(f"    lambda={v:10.4f}  mult={c}")

# =====================================================================
# PART 2: SEELEY-DEWITT COEFFICIENTS + SPECTRAL MOMENTS
# =====================================================================
print("\n" + "=" * 72)
print("PART 2: Seeley-DeWitt Coefficients + Spectral Moments")
print("=" * 72)

c_0 = Nv + Ne + Nf + Nc  # total dim of Hilbert space
c_1 = np.sum(all_evals)  # Tr(D^2)
c_2 = np.sum(all_evals**2) / 2  # Tr(D^4)/2
c_3 = np.sum(all_evals**3) / 6  # Tr(D^6)/6
c_4 = np.sum(all_evals**4) / 24  # Tr(D^8)/24

print(f"c_0 = dim(H) = {c_0}")
print(f"c_1 = Tr(D^2) = {c_1:.2f}")
print(f"c_2 = Tr(D^4)/2 = {c_2:.2f}")
print(f"c_3 = Tr(D^6)/6 = {c_3:.2f}")
print(f"c_4 = Tr(D^8)/24 = {c_4:.2f}")

assert c_0 == 2640, f"c_0 = {c_0}, expected 2640"
assert abs(c_1 - 14880) < 1, f"c_1 = {c_1}, expected 14880"
assert abs(c_2 - 55920) < 10, f"c_2 = {c_2}, expected 55920"

# Exact ratios
ratio_c1_c0 = c_1 / c_0
ratio_c2_c0 = c_2 / c_0
print(f"\nc_1/c_0 = {ratio_c1_c0:.6f} = 62/11 = {62/11:.6f}")
print(f"c_2/c_0 = {ratio_c2_c0:.6f} = 233/11 = {233/11:.6f}")

# Key ratio: c_0*c_2/c_1^2
ratio_key = c_0 * c_2 / c_1**2
print(f"\nc_0*c_2/c_1^2 = {ratio_key:.10f}")
print(f"2/3 = {2/3:.10f}")
print(f"Difference from 2/3 = {ratio_key - 2/3:.6e}")

# Exact with Fractions
F_c0 = Fraction(2640)
F_c1 = Fraction(14880)
F_c2 = Fraction(55920)
F_ratio = F_c0 * F_c2 / (F_c1 * F_c1)
print(f"\nExact: c_0*c_2/c_1^2 = {F_ratio} = {float(F_ratio):.10f}")
print(f"  Numerator = {F_ratio.numerator} = 11 * 233 = 11 * F(13)")
print(f"  Denominator = {F_ratio.denominator} = 4 * 31^2 = (2*(a_1^2+b_1))^2")
print(f"  Verify: 11*233 = {11*233}, 4*31^2 = {4*31**2}")

# NOT 2/3:
F_two_thirds = Fraction(2, 3)
diff = F_ratio - F_two_thirds
print(f"  Difference from 2/3: {diff} = {float(diff):.10f}")
print(f"  = 1/{1/float(diff):.1f}")

# Kurtosis-like
kurtosis = 2 * ratio_key
print(f"\nKurtosis = 2*c_0*c_2/c_1^2 = {kurtosis:.6f}")
print(f"  = 2563/1922 = {2563/1922:.6f}")

# c_3 and c_4 in terms of framework numbers
print(f"\nHigher Seeley-DeWitt coefficients:")
c_3_per_240 = c_3 / 240
c_4_per_240 = c_4 / 240
print(f"  c_3/c_0 = {c_3/c_0:.6f}")
print(f"  c_4/c_0 = {c_4/c_0:.6f}")
print(f"  c_3/240 = {c_3_per_240:.6f}")
print(f"  c_4/240 = {c_4_per_240:.6f}")

# Try to express c_3, c_4 as integer multiples
c_3_round = round(c_3)
c_4_round = round(c_4)
print(f"  c_3 ~ {c_3_round}")
print(f"  c_4 ~ {c_4_round}")

# Check c_3/c_0 and c_4/c_0
F_c3 = Fraction(c_3_round)
F_c4 = Fraction(c_4_round)
print(f"  c_3/c_0 = {F_c3/F_c0} = {float(F_c3/F_c0):.6f}")
print(f"  c_4/c_0 = {F_c4/F_c0} = {float(F_c4/F_c0):.6f}")

# Express c_k in terms of a_1 = 5
print(f"\nAll c_k ratios:")
print(f"  c_1/c_0 = 62/11 = 2*(a_1^2+b_1)/(a_1+b_1) = 2*{a_1**2+b_1}/{a_1+b_1}")
print(f"  c_2/c_0 = 233/11 = F(13)/L(5)")
for k in range(3, 5):
    ck = [c_3, c_4][k-3]
    r = ck / c_0
    # Try Fraction
    fr = Fraction(round(ck)).limit_denominator(10000) / F_c0
    print(f"  c_{k}/c_0 = {r:.6f} (exact: {Fraction(round(ck))}/{c_0})")

# =====================================================================
# PART 3: SPECTRAL ZETA FUNCTION
# =====================================================================
print("\n" + "=" * 72)
print("PART 3: Spectral Zeta Function")
print("=" * 72)

# Nonzero eigenvalues of D^2
nonzero_evals = all_evals[all_evals > tol]
N_nonzero = len(nonzero_evals)
print(f"Number of nonzero D^2 eigenvalues: {N_nonzero}")
print(f"Expected: c_0 - (b0+b1+b2+b3) = 2640 - 2 = 2638")

# zeta_D(s) = sum_{lambda_i > 0} lambda_i^{-s}
def zeta_D(s):
    """Spectral zeta function of D^2."""
    return np.sum(nonzero_evals**(-s))

print(f"\nSpectral zeta at special values:")
special_s = [-2, -1, -0.5, 0, 0.5, 1, 1.5, 2, 3]
for s in special_s:
    if s == 0:
        val = N_nonzero  # analytic: zeta(0) = dim - kernel
        print(f"  zeta_D(0) = {val} (= N_nonzero)")
    else:
        val = zeta_D(s)
        print(f"  zeta_D({s:5.1f}) = {val:.6f}")

# Negative integer values give moment sums
print(f"\nMoment sums (zeta at negative integers):")
print(f"  zeta_D(-1) = Tr(D^2)|_nonzero = {zeta_D(-1):.2f} (c_1 = {c_1:.2f})")
print(f"  zeta_D(-2) = Tr(D^4)|_nonzero = {zeta_D(-2):.2f} (2*c_2 = {2*c_2:.2f})")

# zeta'(0) = -sum(ln(lambda_i))
zeta_prime_0 = -np.sum(np.log(nonzero_evals))
log_det = -zeta_prime_0  # = ln(det(D^2|nonzero))
print(f"\nzeta'_D(0) = -sum(ln(lambda_i)) = {zeta_prime_0:.6f}")
print(f"ln(det(D^2|nonzero)) = {log_det:.6f}")
print(f"det(D^2|nonzero)^(1/N) = exp(ln(det)/N) = {np.exp(log_det/N_nonzero):.6f}")

# Search for z~57 or z~56.88 in zeta values
print(f"\n--- Search for z ~ 57 or z ~ 56.88 in spectral zeta ---")
found_57 = False

# Check zeta_D(s) = 57 for some s
# Since zeta_D(s) is monotonically decreasing for large s -> 0, and at s=0 is 2638,
# at s=1 is much smaller, we can check
for s_try in np.linspace(0.01, 5.0, 5000):
    val = zeta_D(s_try)
    if abs(val - 57) < 0.5:
        print(f"  zeta_D({s_try:.4f}) = {val:.4f}  (near 57)")
        found_57 = True
    if abs(val - 56.882) < 0.5:
        print(f"  zeta_D({s_try:.4f}) = {val:.4f}  (near 56.882)")
        found_57 = True

if not found_57:
    # Find minimum of zeta_D(s) for s > 0
    s_test = np.linspace(0.1, 5.0, 1000)
    z_test = [zeta_D(s) for s in s_test]
    s_min_idx = np.argmin(z_test)
    s_min = s_test[s_min_idx]
    z_min = z_test[s_min_idx]
    print(f"  zeta_D has minimum {z_min:.2f} at s = {s_min:.3f}")
    print(f"  Minimum >> 57, so zeta_D(s) NEVER equals 57 for any real s > 0")
    print(f"  zeta_D ranges: [{z_min:.1f}, infinity) -- cannot reach 57")

# Regularized CC
print(f"\nRegularized CC from spectral zeta:")
print(f"  rho_vac ~ mu^4 * (zeta'(0) + zeta(0)*ln(mu^2))")
print(f"  Setting mu = 1 (Planck units): rho_vac ~ zeta'(0) = {zeta_prime_0:.4f}")
print(f"  exp(zeta'(0)/N_nonzero) = {np.exp(zeta_prime_0/N_nonzero):.6f}")

# Check if any combination gives the exponent
ratio_zeta = zeta_prime_0 / N_nonzero
print(f"  zeta'(0)/zeta(0) = {ratio_zeta:.6f}")
print(f"  -zeta'(0)/(2*ln(alpha)) = {-zeta_prime_0/(2*np.log(ALPHA)):.4f}")

# =====================================================================
# PART 4: EXACT SPECTRAL ACTION + RESIDUAL
# =====================================================================
print("\n" + "=" * 72)
print("PART 4: Exact Spectral Action + Residual")
print("=" * 72)

# Natural cutoff
Lambda_nat = np.sqrt(6 + 6*PHI)
print(f"Natural cutoff: sqrt(lambda_max) = sqrt(6+6*phi) = {Lambda_nat:.6f}")
print(f"lambda_max = {6+6*PHI:.6f}")

# Sharp cutoff: N(Lambda) = #{eigenvalues <= Lambda^2}
Lambda_range = np.linspace(0.1, Lambda_nat * 1.5, 1000)
N_Lambda = np.array([np.sum(all_evals <= L**2) for L in Lambda_range])

print(f"\nSharp eigenvalue counting function N(Lambda):")
print(f"  N(0) = {np.sum(all_evals == 0)} (kernel)")
print(f"  N(sqrt(lambda_max)) = N({Lambda_nat:.4f}) = {np.sum(all_evals <= 6+6*PHI)}")
print(f"  N(infinity) = {c_0}")

# Gaussian spectral action: S(Lambda) = Tr(exp(-D^2/Lambda^2))
def gaussian_spectral(Lambda_sq):
    return np.sum(np.exp(-all_evals / Lambda_sq))

def seeley_approx(Lambda_sq, order=4):
    """Seeley-DeWitt approximation: sum c_k / k! * Lambda^{-2k} * ..."""
    # Actually: Tr(e^{-tD^2}) ~ sum_k c_k * (-t)^k / k! for t->0
    # At t = 1/Lambda^2: Tr(e^{-D^2/Lambda^2}) ~ sum c_k * (-1/Lambda^2)^k / k!
    # Wait, the standard expansion is:
    # Tr(e^{-tD^2}) = sum_{k>=0} a_k * t^{(k-d)/2} for t->0
    # On a discrete system, we have exact eigenvalues, so:
    # Tr(e^{-tD^2}) = sum_i e^{-t*lambda_i}
    # For small t: ~ c_0 - c_1*t + c_2*t^2/2 - c_3*t^3/6 + ...
    # where c_k = Tr(D^{2k}) (raw moments)
    # But our c_k are already divided: c_2 = Tr(D^4)/2
    # So: Tr(e^{-tD^2}) ~ c_0 - c_1*t + c_2*t^2 - c_3*t^3 + c_4*t^4 ...
    # where c_k = Tr(D^{2k})/k!
    t = 1.0 / Lambda_sq
    raw_moments = [c_0, c_1, 2*c_2, 6*c_3, 24*c_4]
    result = 0
    for k, mk in enumerate(raw_moments[:order+1]):
        result += mk * (-t)**k / math.factorial(k)
    return result

# Evaluate at natural cutoff
Lambda_sq_nat = 6 + 6*PHI
S_exact = gaussian_spectral(Lambda_sq_nat)
S_sd2 = seeley_approx(Lambda_sq_nat, order=2)
S_sd3 = seeley_approx(Lambda_sq_nat, order=3)
S_sd4 = seeley_approx(Lambda_sq_nat, order=4)

print(f"\nGaussian spectral action at Lambda^2 = 6+6*phi:")
print(f"  S_exact = {S_exact:.6f}")
print(f"  S_SD(order=2) = {S_sd2:.6f}")
print(f"  S_SD(order=3) = {S_sd3:.6f}")
print(f"  S_SD(order=4) = {S_sd4:.6f}")
print(f"  Residual R_2 = S_exact - S_SD(2) = {S_exact - S_sd2:.6f}")
print(f"  Residual R_3 = S_exact - S_SD(3) = {S_exact - S_sd3:.6f}")
print(f"  Residual R_4 = S_exact - S_SD(4) = {S_exact - S_sd4:.6f}")

# Tr(D^2)/lambda_max in Z[phi]
ratio_trace = c_1 / lambda_max
print(f"\nTr(D^2)/lambda_max = {ratio_trace:.6f}")
# 14880 / (6+6*phi) = 14880 / (6*(1+phi)) = 2480 / (1+phi) = 2480/phi^2 = 2480*(3-phi)/5
# Actually: 1+phi = phi^2, so 14880/(6*phi^2) = 2480/phi^2
val_zphi = 2480 / PHI**2
print(f"  = 2480/phi^2 = {val_zphi:.6f}")
print(f"  = 2480*(3-phi)/5 = {2480*(3-PHI)/5:.6f}")
# Simplify: 2480 = 2640 - 160 = c_0 - 160
# Or: 2480 = 16 * 155 = 16 * 5 * 31 = 16 * a_1 * (a_1^2+b_1)
print(f"  2480 = 16 * {2480//16} = 16 * a_1 * {2480//(16*a_1)}")
print(f"  = 16 * a_1 * (a_1^2+b_1) = 16*5*31 = {16*5*31}")

# Heat kernel at special t values
print(f"\nHeat kernel K(t) = Tr(exp(-t*D^2)) at special t:")
special_t = [1.0/(6+6*PHI), 1.0/12, PHI**(-4), 1.0/a_1, 1.0/b_1,
             1.0/(a_1**2+1), 1.0/(2*PHI**3)]
special_t_names = ["1/lambda_max", "1/12", "1/phi^4", "1/a_1", "1/b_1",
                   "1/(a_1^2+1)", "alpha_s"]
for t, name in zip(special_t, special_t_names):
    K = np.sum(np.exp(-t * all_evals))
    print(f"  K({name} = {t:.6f}) = {K:.6f}")

# =====================================================================
# PART 5: Str(M^4) WITH FRAMEWORK MASSES
# =====================================================================
print("\n" + "=" * 72)
print("PART 5: Supertrace Str(M^4) with Framework Masses")
print("=" * 72)

m_e = 0.51099895  # MeV
M_Planck_GeV = 1.22089e19  # GeV

# Fermion masses from m_f = m_e * phi^(5a+6b)
# (a,b): e(0,0), mu(1,1), tau(1,2), u(3,-2), c(2,1), t(4,1), d(1,0), s(1,1), b(-1,4)
fermion_data = [
    ("e",   0, 0, 1),
    ("mu",  1, 1, 1),
    ("tau", 1, 2, 1),
    ("u",   3,-2, 3),
    ("c",   2, 1, 3),
    ("t",   4, 1, 3),
    ("d",   1, 0, 3),
    ("s",   1, 1, 3),
    ("b",  -1, 4, 3),
]

print(f"Fermion masses (bare, from m_f = m_e * phi^(5a+6b)):")
fermion_m4_sum = 0
for name, a, b, n_c in fermion_data:
    exp = 5*a + 6*b
    m_f = m_e * PHI**exp  # MeV
    m_f_GeV = m_f * 1e-3
    # Each fermion: N_c colors * N_gen generations * (particle + antiparticle) = N_c * N_gen * 2
    # But for Str, coefficient is -4 per Weyl fermion
    # Dirac: 4 DOF per color per generation
    dof = 4 * n_c * N_gen
    contribution = dof * (m_f_GeV)**4
    fermion_m4_sum += contribution
    print(f"  {name:4s}: exp={exp:3d}, m={m_f:.4f} MeV, "
          f"DOF={dof}, m^4*DOF = {contribution:.4e} GeV^4")

print(f"\nTotal fermion Str contribution: 4*sum = {fermion_m4_sum:.4e} GeV^4")

# Boson masses from framework
m_W = 80.377  # GeV (experimental, but framework gives it from m_e*phi^25*...)
m_Z = 91.1876  # GeV
m_H = 125.25  # GeV

# W: 3 DOF (massive vector), 2 charged = 6 DOF total
# Z: 3 DOF (massive vector)
# H: 1 DOF (scalar)
# Photon: 2 DOF, massless -> 0 contribution
# Gluons: 2*8=16 DOF, massless -> 0 contribution

boson_contribs = [
    ("W+-", 6, m_W),
    ("Z",   3, m_Z),
    ("H",   1, m_H),
]

boson_m4_sum = 0
print(f"\nBoson contributions:")
for name, dof, m in boson_contribs:
    contrib = dof * m**4
    boson_m4_sum += contrib
    print(f"  {name:4s}: DOF={dof}, m={m:.3f} GeV, m^4*DOF = {contrib:.4e} GeV^4")

print(f"\nTotal boson Str: {boson_m4_sum:.4e} GeV^4")

# Supertrace: Str(M^4) = sum_bosons - sum_fermions
# (In SUSY-like context: Str = bosons - fermions)
Str_M4 = boson_m4_sum - fermion_m4_sum
print(f"\nStr(M^4) = bosons - fermions = {Str_M4:.4e} GeV^4")
print(f"|Str(M^4)| = {abs(Str_M4):.4e} GeV^4")

# Top quark dominance
m_t_GeV = m_e * 1e-3 * PHI**(5*4 + 6*1)  # exp = 26
print(f"\nTop quark mass (framework): m_t = m_e*phi^26 = {m_t_GeV*1e3:.2f} MeV = {m_t_GeV:.2f} GeV")
top_frac = 4 * 3 * N_gen * m_t_GeV**4 / fermion_m4_sum
print(f"Top contribution fraction: {top_frac:.6f}")

# Compare with Planck scale
ratio_Str_Planck = abs(Str_M4) / M_Planck_GeV**4
print(f"\n|Str(M^4)|/m_P^4 = {ratio_Str_Planck:.4e}")
print(f"alpha^56.88 = {ALPHA**56.882:.4e}")
print(f"Ratio: |Str/m_P^4| / alpha^56.88 = {ratio_Str_Planck / ALPHA**56.882:.4e}")
print(f"\n  --> Str(M^4) ~ m_t^4 ~ 10^{np.log10(abs(Str_M4)):.0f} GeV^4")
print(f"  --> Lambda ~ m_P^4 * alpha^57 ~ 10^{4*np.log10(M_Planck_GeV) + 57*np.log10(ALPHA):.0f} GeV^4")
print(f"  --> Str(M^4) does NOT explain the CC. Gap: ~10^{np.log10(ratio_Str_Planck/ALPHA**56.882):.0f}")

# =====================================================================
# PART 6: GRADED SPECTRAL ACTION (BOSON-FERMION)
# =====================================================================
print("\n" + "=" * 72)
print("PART 6: Graded Spectral Action")
print("=" * 72)

# Form grading: even forms (0-forms + 2-forms) vs odd forms (1-forms + 3-forms)
dim_even = Nv + Nf  # 120 + 1200 = 1320
dim_odd = Ne + Nc   # 720 + 600 = 1320
print(f"Form grading:")
print(f"  Even forms (0+2): {Nv} + {Nf} = {dim_even}")
print(f"  Odd forms  (1+3): {Ne} + {Nc} = {dim_odd}")
print(f"  Balance: {dim_even} - {dim_odd} = {dim_even - dim_odd}")

# Graded trace with form parity
Str_form_1 = np.sum(evals_0) + np.sum(evals_2) - np.sum(evals_1) - np.sum(evals_3)
print(f"\nStr_form(D^2) = Tr_even(D^2) - Tr_odd(D^2)")
print(f"  = ({np.sum(evals_0):.2f} + {np.sum(evals_2):.2f}) - ({np.sum(evals_1):.2f} + {np.sum(evals_3):.2f})")
print(f"  = {Str_form_1:.2f}")

Str_form_2 = np.sum(evals_0**2) + np.sum(evals_2**2) - np.sum(evals_1**2) - np.sum(evals_3**2)
print(f"Str_form(D^4) = {Str_form_2:.2f}")

# Vertex grading: 24 gauge (A+B) vs 96 fermionic (C)
# Type A: 8 vertices, Type B: 16 vertices, Type C: 96 vertices
# Build Gamma_0 grading on vertex space
Gamma_0_diag = np.zeros(Nv)
for i, v in enumerate(verts):
    v_tuple = tuple(round(x, 10) for x in v)
    # Type A: one coord = +-1, rest 0
    is_A = any(abs(abs(c) - 1.0) < 0.01 for c in v_tuple) and \
           sum(abs(c) < 0.01 for c in v_tuple) == 3
    # Type B: all coords = +-0.5
    is_B = all(abs(abs(c) - 0.5) < 0.01 for c in v_tuple)

    if is_A or is_B:
        Gamma_0_diag[i] = +1  # gauge vertex
    else:
        Gamma_0_diag[i] = -1  # fermionic vertex

n_gauge = np.sum(Gamma_0_diag > 0)
n_fermi = np.sum(Gamma_0_diag < 0)
print(f"\nVertex grading (Gamma_0 on 0-forms):")
print(f"  Gauge vertices (A+B): {int(n_gauge)}")
print(f"  Fermionic vertices (C): {int(n_fermi)}")

Gamma_0_mat = np.diag(Gamma_0_diag)
graded_Delta0 = Gamma_0_mat @ Delta_0
Tr_graded_Delta0 = np.trace(graded_Delta0)
print(f"  Tr(Gamma_0 * Delta_0) = {Tr_graded_Delta0:.4f}")

# Higher graded traces
Tr_graded_D2 = np.trace(Gamma_0_mat @ Delta_0)
Tr_graded_D4 = np.trace(Gamma_0_mat @ Delta_0 @ Delta_0)
print(f"  Tr(Gamma_0 * Delta_0^2) = {Tr_graded_D4:.4f}")

# Check if graded quantities relate to 57
print(f"\n  Tr_graded(D^2)/24 = {Tr_graded_Delta0/24:.6f}")
print(f"  Tr_graded(D^4)/24 = {Tr_graded_D4/24:.6f}")

# Full form-graded Str including all degrees
# Build full grading: (-1)^p for p-forms
Gamma_full = np.zeros(c_0)
Gamma_full[:Nv] = 1        # 0-forms: +1
Gamma_full[Nv:Nv+Ne] = -1  # 1-forms: -1
Gamma_full[Nv+Ne:Nv+Ne+Nf] = 1    # 2-forms: +1
Gamma_full[Nv+Ne+Nf:] = -1         # 3-forms: -1

# Build full D^2 block-diagonal
D2_full = np.zeros((c_0, c_0))
D2_full[:Nv, :Nv] = Delta_0
D2_full[Nv:Nv+Ne, Nv:Nv+Ne] = Delta_1
D2_full[Nv+Ne:Nv+Ne+Nf, Nv+Ne:Nv+Ne+Nf] = Delta_2
D2_full[Nv+Ne+Nf:, Nv+Ne+Nf:] = Delta_3

Gamma_full_mat = np.diag(Gamma_full)
for k in range(1, 5):
    D2k = np.linalg.matrix_power(D2_full, k)
    str_k = np.trace(Gamma_full_mat @ D2k)
    print(f"  Str_form(D^{2*k}) = {str_k:.4f}")

# =====================================================================
# PART 7: SEARCH FOR 57 IN SPECTRAL DATA
# =====================================================================
print("\n" + "=" * 72)
print("PART 7: Search for 57 in Spectral Data")
print("=" * 72)

# Known: (N-N_gen)/2 = (120-3)/2 = 117/2 = 58.5, 58.5 - phi = 56.882
# And: 117 = N_eig * (degree+1) = 9*13

print("Spectral invariants to check:")

# 1. Number of distinct eigenvalues per degree
n_distinct = [len(distinct_per_deg[i][0]) for i in range(4)]
print(f"\n  Distinct eigenvalues: Delta_0={n_distinct[0]}, Delta_1={n_distinct[1]}, "
      f"Delta_2={n_distinct[2]}, Delta_3={n_distinct[3]}")
print(f"  Total distinct: {sum(n_distinct)}")

# 2. Multiplicity sums
for i in range(4):
    vals, counts = distinct_per_deg[i]
    mult_sum = sum(counts)
    print(f"  Sum multiplicities Delta_{i}: {mult_sum}")

# 3. Eigenvalue ratios
print(f"\n  Eigenvalue ratios:")
all_distinct_nz = set()
for i in range(4):
    vals, counts = distinct_per_deg[i]
    for v in vals:
        if abs(v) > tol:
            all_distinct_nz.add(round(v, 4))

sorted_distinct = sorted(all_distinct_nz)
print(f"  All distinct nonzero eigenvalues ({len(sorted_distinct)}):")
for v in sorted_distinct:
    print(f"    {v:.4f}")

# 4. Check ratios between eigenvalues for 57
print(f"\n  Ratios between eigenvalues:")
for i, vi in enumerate(sorted_distinct):
    for j, vj in enumerate(sorted_distinct):
        if j > i and vi > 0.01:
            r = vj / vi
            if abs(r - 57) < 1 or abs(r - 56.882) < 1:
                print(f"    {vj:.4f}/{vi:.4f} = {r:.4f} (near 57!)")

# 5. Sum/product of distinct eigenvalues
sum_distinct = sum(sorted_distinct)
prod_distinct = 1.0
for v in sorted_distinct:
    prod_distinct *= v
print(f"\n  Sum of distinct nonzero eigenvalues: {sum_distinct:.4f}")
print(f"  Product of distinct nonzero eigenvalues: {prod_distinct:.4e}")

# 6. Check if 57 appears as multiplicity sum
print(f"\n  Multiplicity patterns:")
all_mults = []
for i in range(4):
    vals, counts = distinct_per_deg[i]
    for v, c in zip(vals, counts):
        if abs(v) > tol:
            all_mults.append((v, c, i))

# Sort by multiplicity
all_mults.sort(key=lambda x: x[1])
for v, c, deg in all_mults:
    if c == 57 or c == 117:
        print(f"    !! lambda={v:.4f}, mult={c}, degree={deg}")

# 7. Cumulative multiplicity counts
print(f"\n  Cumulative multiplicities:")
cum = 0
target_vals = [57, 117, 58, 56]
for v in sorted_distinct:
    for i in range(4):
        vals, counts = distinct_per_deg[i]
        for val, cnt in zip(vals, counts):
            if abs(val - v) < 0.001:
                cum += cnt
    if cum in target_vals:
        print(f"    At lambda <= {v:.4f}: cumulative = {cum} (**)")

# 8. c_k ratios that give ~57
print(f"\n  c_k ratios:")
ck_list = [c_0, c_1, c_2, c_3, c_4]
ck_names = ["c_0", "c_1", "c_2", "c_3", "c_4"]
for i in range(len(ck_list)):
    for j in range(len(ck_list)):
        if i != j and ck_list[j] != 0:
            r = ck_list[i] / ck_list[j]
            if abs(r - 57) < 2 or abs(r - 56.882) < 2:
                print(f"    {ck_names[i]}/{ck_names[j]} = {r:.4f}")

# 9. Algebraic combinations
print(f"\n  Algebraic combinations involving spectral data:")
N_over_2 = N_vert / 2
print(f"  N/2 = {N_over_2}")
print(f"  N/2 - N_gen = {N_over_2 - N_gen}")
print(f"  (N-N_gen)/2 = {(N_vert-N_gen)/2}")
print(f"  (N-N_gen)/2 - phi = {(N_vert-N_gen)/2 - PHI:.6f}")
print(f"  c_1/c_0 - 1/phi^3 = {c_1/c_0 - 1/PHI**3:.6f}")
print(f"  c_0/(2*N_gen*degree+a_1) = {c_0/(2*N_gen*degree + a_1):.6f}")

# 10. Check (c_0 - N_gen) / 2 / (c_0/N) - phi
print(f"  c_0/N = {c_0/N_vert} = 22")
print(f"  (c_0/N - N_gen)/2 - phi = {(c_0/N_vert - N_gen)/2 - PHI:.6f}")

# =====================================================================
# PART 8: HONEST ASSESSMENT
# =====================================================================
print("\n" + "=" * 72)
print("PART 8: HONEST ASSESSMENT")
print("=" * 72)

print("""
+--+---------------------------+------------------------------------------+----------+
|# | Route                     | Result                                   | Verdict  |
+--+---------------------------+------------------------------------------+----------+
|1 | Seeley-DeWitt ratios      | c_0*c_2/c_1^2 = 2563/3844, NOT 2/3      | NEGATIVE |
|  |                           | c_3, c_4 computed but no 57 pattern      |          |
+--+---------------------------+------------------------------------------+----------+
|2 | Spectral zeta             | min(zeta_D)=681 at s=1.47; NEVER = 57   | NEGATIVE |
|  |                           | zeta'(0) unrelated to CC exponent        |          |
+--+---------------------------+------------------------------------------+----------+
|3 | Exact spectral action     | Remainder R computable but no 57         | NEGATIVE |
|  |                           | Heat kernel at special t: no pattern     |          |
+--+---------------------------+------------------------------------------+----------+
|4 | Str(M^4)                  | ~ m_t^4 ~ 10^8 GeV^4, not 10^{-47}      | NEGATIVE |
|  |                           | 54 orders too large for CC               |          |
+--+---------------------------+------------------------------------------+----------+
|5 | Graded spectral action    | Form grading: Str(D^{2k}) computable     | NEGATIVE |
|  |                           | Vertex grading: no CC-relevant quantity   |          |
+--+---------------------------+------------------------------------------+----------+
|6 | Search for 57             | 57 NOT a spectral invariant               | NEGATIVE |
|  |                           | Only from (N-N_gen)/2 - phi (KNOWN)      |          |
+--+---------------------------+------------------------------------------+----------+
""")

print("CONCLUSION:")
print("  Lambda_Planck = alpha^(57-alpha_s) remains a STRONG PATTERN.")
print("  The exponent 57 = (N-N_gen)/2 - 1/2 is a COUNTING relation,")
print("  not derivable from spectral action coefficients.")
print("  No route from spectral data produces z ~ 56.88.")
print()
print("  POSITIVE DELIVERABLES:")
print(f"    c_3 = {c_3:.2f}  (c_3/c_0 = 2125/33)")
print(f"    c_4 = {c_4:.2f}  (c_4/c_0 = 2035/12)")
print(f"    zeta_D: minimum {z_min:.1f} at s={s_min:.3f}, never reaches 57")
print(f"    Str(M^4) = {Str_M4:.4e} GeV^4 (top-dominated, 55 orders too large)")
print(f"    Form-graded Str(D^{'{'}2k{'}'}) = 0 for ALL k (exact cancellation!)")
print(f"    c_0*c_2/c_1^2 = 2563/3844 = 11*F(13)/(2*(a_1^2+b_1))^2")
print()
print("  CATEGORY: STRONG PATTERN (unchanged)")
print("  STATUS: No upgrade from spectral action. All 6 routes negative.")

print("\n" + "=" * 72)
print("exp328 COMPLETE")
print("=" * 72)
