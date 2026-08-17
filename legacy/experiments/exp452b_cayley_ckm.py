"""
exp452b: Cayley graph of 2I (120 nodes) — CKM exponents from spectral structure
================================================================================
The McKay graph (9 nodes, diam=7) is too small for n_13=12.
The Cayley graph has 120 nodes and phi in its spectrum — test here.
"""
import numpy as np
from math import sqrt, log, sin, cos, pi

phi = (1 + sqrt(5)) / 2
phi_p = (1 - sqrt(5)) / 2  # = -1/phi
ln_phi = log(phi)
a1 = 5; b1 = 6; N = 120

print("=" * 70)
print("exp452b: Cayley Graph of 2I — CKM from Spectral Structure")
print("=" * 70)

# ============================================================
# 1. Character table of 2I (9 irreps x 9 classes)
# ============================================================
# Classes: 1a, 2a, 3a, 4a, 5a, 5b, 6a, 10a, 10b
# Half-angles theta/2 for SU(2) character formula
half_angles = {
    '1a': 0,       '2a': pi,      '3a': pi/3,
    '4a': pi/2,    '5a': 4*pi/5,  '5b': 3*pi/5,
    '6a': 2*pi/3,  '10a': pi/5,   '10b': 2*pi/5
}
class_sizes = [1, 1, 20, 30, 12, 12, 20, 12, 12]
class_names = ['1a', '2a', '3a', '4a', '5a', '5b', '6a', '10a', '10b']

def chi_su2(j, alpha):
    """SU(2) character for spin j at half-angle alpha."""
    if abs(alpha) < 1e-12:
        return 2*j + 1
    if abs(alpha - pi) < 1e-12:
        return (-1)**(int(2*j)) * (2*j + 1)
    return sin((2*j+1)*alpha) / sin(alpha)

# Compute SU(2) characters for j = 0, 1/2, 1, ..., 4
su2_chars = {}
for j_half in range(9):  # j = 0, 0.5, 1, ..., 4
    j = j_half / 2.0
    su2_chars[j] = []
    for cn in class_names:
        su2_chars[j].append(chi_su2(j, half_angles[cn]))

# Irrep characters of 2I (dims [1,2,3,4,5,6,4,2,3])
# rho_0..rho_5 = V_0..V_{5/2} (direct restriction)
# rho_6 = V_4 - rho_4 = chi_4^SU2 - chi_2^SU2
# rho_7 = V_{7/2} - rho_5 = chi_{7/2}^SU2 - chi_{5/2}^SU2
# rho_8 = V_3 - rho_6 = chi_3^SU2 - chi_6
dims = [1, 2, 3, 4, 5, 6, 4, 2, 3]
chi = np.zeros((9, 9))  # chi[irrep][class]
for k in range(6):
    for c in range(9):
        chi[k][c] = su2_chars[k/2.0][c]

# rho_6 = V_4 - V_2
for c in range(9):
    chi[6][c] = su2_chars[4.0][c] - su2_chars[2.0][c]

# rho_7 = V_{7/2} - V_{5/2}
for c in range(9):
    chi[7][c] = su2_chars[3.5][c] - su2_chars[2.5][c]

# rho_8 = V_3 - rho_6
for c in range(9):
    chi[8][c] = su2_chars[3.0][c] - chi[6][c]

print("\nCharacter table of 2I:")
print(f"{'Irrep':>6s} {'dim':>3s}", end="")
for cn in class_names:
    print(f" {cn:>7s}", end="")
print()
for k in range(9):
    print(f"rho_{k:d}  {dims[k]:3d}", end="")
    for c in range(9):
        print(f" {chi[k][c]:7.3f}", end="")
    print()

# Verify dimensions
for k in range(9):
    assert abs(chi[k][0] - dims[k]) < 1e-10, f"dim mismatch for rho_{k}"

# Verify orthogonality
print("\nOrthogonality check (should be 120*delta_ij):")
for i in range(9):
    for j in range(i, min(i+2, 9)):
        ip = sum(class_sizes[c] * chi[i][c] * chi[j][c] for c in range(9))
        expected = N if i == j else 0
        status = "OK" if abs(ip - expected) < 1e-8 else "FAIL"
        if i == j or abs(ip) > 0.01:
            print(f"  <rho_{i}|rho_{j}> = {ip:.4f} (expect {expected}) {status}")

# ============================================================
# 2. Cayley eigenvalues (generators = class 10a, degree 12)
# ============================================================
gen_class = 7  # index of 10a in class_names
degree = 12

cayley_evals = []
for k in range(9):
    lam = degree * (1 - chi[k][gen_class] / dims[k])
    cayley_evals.append(lam)

print("\nCayley graph eigenvalues (class 10a generators):")
print(f"{'Irrep':>6s} {'dim':>3s} {'chi(10a)':>9s} {'lambda':>10s} {'mult':>5s}")
for k in range(9):
    mult = dims[k]**2
    print(f"rho_{k:d}  {dims[k]:3d} {chi[k][gen_class]:9.4f} {cayley_evals[k]:10.4f} {mult:5d}")

# ============================================================
# 3. Key eigenvalue ratios (log base phi)
# ============================================================
print("\n" + "=" * 70)
print("EIGENVALUE RATIOS (log base phi)")
print("=" * 70)

quark_nodes = {'u': 1, 'd': 2, 's': 3, 'c': 5, 't': 7, 'b': 8}
quark_labels = list(quark_nodes.keys())

print(f"\n{'Pair':>6s} {'lambda_i':>10s} {'lambda_j':>10s} {'ratio':>10s} {'log_phi':>10s}")
for qi, ni in quark_nodes.items():
    for qj, nj in quark_nodes.items():
        if ni < nj:
            li = cayley_evals[ni]
            lj = cayley_evals[nj]
            ratio = lj / li
            if ratio > 0:
                lp = log(ratio) / ln_phi
            else:
                lp = float('nan')
            print(f"{qi+qj:>6s} {li:10.4f} {lj:10.4f} {ratio:10.5f} {lp:10.4f}")

# Check exact ratios
print("\n--- Exact ratio checks ---")
lam_u = 12 - 6*phi; lam_d = 12 - 4*phi; lam_s = 9
lam_c = 14; lam_t = 6 + 6*phi; lam_b = 8 + 4*phi

# lambda_t / lambda_u
r_tu = lam_t / lam_u
print(f"lambda_t/lambda_u = (6+6phi)/(12-6phi) = {r_tu:.10f}")
print(f"  phi^4 = {phi**4:.10f}")
print(f"  EXACT? {abs(r_tu - phi**4) < 1e-10}")
# Proof: 6(1+phi) / (6(2-phi)) = phi^2 / (2-phi) = phi^2 * phi^2 = phi^4
# since 2-phi = 1/phi^2

# lambda_b / lambda_d
r_bd = lam_b / lam_d
print(f"\nlambda_b/lambda_d = (8+4phi)/(12-4phi) = {r_bd:.10f}")
print(f"  phi^2 = {phi**2:.10f}")
print(f"  EXACT? {abs(r_bd - phi**2) < 1e-10}")
# Proof: (10+2sqrt5)/(10-2sqrt5) = (5+sqrt5)^2/20 = (30+10sqrt5)/20 = phi^2

# ============================================================
# 4. Resistance distance R(e, g) via representation theory
# ============================================================
print("\n" + "=" * 70)
print("RESISTANCE DISTANCES on Cayley graph R(e, g)")
print("=" * 70)

# R(e,g) = (2/|G|) sum_{rho != rho_0} (d_rho / lambda_rho) * [d_rho - chi_rho(g)]
R_class = []
for c in range(9):
    R = 0
    for k in range(1, 9):  # skip rho_0
        R += (dims[k] / cayley_evals[k]) * (dims[k] - chi[k][c])
    R *= 2.0 / N
    R_class.append(R)

print(f"\n{'Class':>6s} {'Size':>5s} {'R(e,g)':>12s} {'R/ln(phi)':>10s}")
for c in range(9):
    rlp = R_class[c] / ln_phi
    print(f"{class_names[c]:>6s} {class_sizes[c]:5d} {R_class[c]:12.6f} {rlp:10.4f}")

# ============================================================
# 5. "Spectral propagator" between irreps
# ============================================================
print("\n" + "=" * 70)
print("SPECTRAL PROPAGATOR between irreps via Cayley Green's function")
print("=" * 70)

# G_rep(rho_i, rho_j) = (1/|G|) sum_C |C| * chi_i(C) * L+(e,C) * chi_j(C)
# where L+(e,C) = (1/|G|) sum_{rho!=0} (d_rho/lambda_rho) * chi_rho(C)

Lplus_class = []
for c in range(9):
    lp = 0
    for k in range(1, 9):
        lp += (dims[k] / cayley_evals[k]) * chi[k][c]
    lp /= N
    Lplus_class.append(lp)

G_rep = np.zeros((9, 9))
for i in range(9):
    for j in range(9):
        val = 0
        for c in range(9):
            val += class_sizes[c] * chi[i][c] * Lplus_class[c] * chi[j][c]
        G_rep[i][j] = val / N

print("\nSpectral propagator G_rep(rho_i, rho_j):")
print(f"{'':>8s}", end="")
for k in range(9):
    print(f" {'rho_'+str(k):>8s}", end="")
print()
for i in range(9):
    print(f"{'rho_'+str(i):>8s}", end="")
    for j in range(9):
        print(f" {G_rep[i][j]:8.5f}", end="")
    print()

# Extract quark-quark propagators
print("\nQuark-quark propagators:")
for qi, ni in quark_nodes.items():
    for qj, nj in quark_nodes.items():
        if ni <= nj:
            print(f"  G({qi},{qj}) = G(rho_{ni},rho_{nj}) = {G_rep[ni][nj]:.6f}")

# "Propagator distance" = -ln(G_ij / sqrt(G_ii*G_jj)) / ln(phi)
print("\nPropagator distance -ln(G_ij/sqrt(G_ii*G_jj))/ln(phi):")
gen_pairs = [('u','c',3), ('c','t',7), ('u','t',12),
             ('d','s',3), ('s','b',7), ('d','b',12)]
for qi, qj, n_target in gen_pairs:
    ni, nj = quark_nodes[qi], quark_nodes[qj]
    Gij = G_rep[ni][nj]
    Gii = G_rep[ni][ni]
    Gjj = G_rep[nj][nj]
    if Gii > 0 and Gjj > 0 and Gij != 0:
        normalized = abs(Gij) / sqrt(abs(Gii * Gjj))
        if normalized > 0 and normalized < 1:
            d = -log(normalized) / ln_phi
        else:
            d = float('nan')
    else:
        d = float('nan')
    print(f"  d({qi},{qj}) = {d:.4f}  (need n={n_target})")

# ============================================================
# 6. Key test: spectral decay n_ij = ln(lambda_j/lambda_i) / ln(phi)
# ============================================================
print("\n" + "=" * 70)
print("SPECTRAL DECAY: n_ij = ln(lambda_j/lambda_i) / ln(phi)")
print("=" * 70)

print("\nUp-type (u, c, t):")
print(f"  n(u->c) = ln({lam_c:.3f}/{lam_u:.3f})/ln(phi) = {log(lam_c/lam_u)/ln_phi:.4f}  (need 3)")
print(f"  n(c->t) = ln({lam_t:.3f}/{lam_c:.3f})/ln(phi) = {log(lam_t/lam_c)/ln_phi:.4f}  (need 7)")
print(f"  n(u->t) = ln({lam_t:.3f}/{lam_u:.3f})/ln(phi) = {log(lam_t/lam_u)/ln_phi:.4f}  (need 12)")
print(f"  SUM = {log(lam_c/lam_u)/ln_phi + log(lam_t/lam_c)/ln_phi:.4f}  (= n(u->t))")

print("\nDown-type (d, s, b):")
print(f"  n(d->s) = ln({lam_s:.3f}/{lam_d:.3f})/ln(phi) = {log(lam_s/lam_d)/ln_phi:.4f}  (need 3)")
print(f"  n(s->b) = ln({lam_b:.3f}/{lam_s:.3f})/ln(phi) = {log(lam_b/lam_s)/ln_phi:.4f}  (need 7)")
print(f"  n(d->b) = ln({lam_b:.3f}/{lam_d:.3f})/ln(phi) = {log(lam_b/lam_d)/ln_phi:.4f}  (need 12)")

# ============================================================
# 7. Products and sums of log-ratios
# ============================================================
print("\n" + "=" * 70)
print("COMBINED (up x down) spectral exponents")
print("=" * 70)
n_u_c = log(lam_c/lam_u)/ln_phi
n_c_t = log(lam_t/lam_c)/ln_phi
n_u_t = log(lam_t/lam_u)/ln_phi  # = 4 exactly
n_d_s = log(lam_s/lam_d)/ln_phi
n_s_b = log(lam_b/lam_s)/ln_phi
n_d_b = log(lam_b/lam_d)/ln_phi  # = 2 exactly

print(f"Up: {n_u_c:.4f}, {n_c_t:.4f}, {n_u_t:.4f}")
print(f"Down: {n_d_s:.4f}, {n_s_b:.4f}, {n_d_b:.4f}")

print(f"\nSums (up+down):")
print(f"  12: {n_u_c + n_d_s:.4f}  (need 3)")
print(f"  23: {n_c_t + n_s_b:.4f}  (need 7)")
print(f"  13: {n_u_t + n_d_b:.4f}  = 4+2 = 6  (need 12)")

print(f"\nProducts (up*down):")
print(f"  12: {n_u_c * n_d_s:.4f}  (need 3)")
print(f"  23: {n_c_t * n_s_b:.4f}  (need 7)")
print(f"  13: {n_u_t * n_d_b:.4f}  = 4*2 = 8  (need 12)")

print(f"\nScaled sums a1*(up+down)/2:")
n12_sc = a1 * (n_u_c + n_d_s) / 2
n23_sc = a1 * (n_c_t + n_s_b) / 2
n13_sc = a1 * (n_u_t + n_d_b) / 2
print(f"  12: {n12_sc:.4f}  (need 3)")
print(f"  23: {n23_sc:.4f}  (need 7)")
print(f"  13: {n13_sc:.1f}  = 5*3 = 15  (need 12)")

# ============================================================
# 8. Test: CKM from spectral DIFFERENCES lambda_j - lambda_i
# ============================================================
print("\n" + "=" * 70)
print("SPECTRAL DIFFERENCES lambda_j - lambda_i")
print("=" * 70)
print(f"Up: lam_c - lam_u = {lam_c - lam_u:.4f}")
print(f"Up: lam_t - lam_c = {lam_t - lam_c:.4f}")
print(f"Up: lam_t - lam_u = {lam_t - lam_u:.4f}")
print(f"Down: lam_s - lam_d = {lam_s - lam_d:.4f}")
print(f"Down: lam_b - lam_s = {lam_b - lam_s:.4f}")
print(f"Down: lam_b - lam_d = {lam_b - lam_d:.4f}")

# Check if differences have phi structure
print(f"\n(lam_t - lam_u) = {lam_t - lam_u:.6f} = 12*phi - 6 = 6*(2phi-1) = 6*sqrt(5)")
print(f"  6*sqrt(5) = {6*sqrt(5):.6f}")
print(f"  EXACT? {abs(lam_t - lam_u - 6*sqrt(5)) < 1e-10}")

print(f"\n(lam_b - lam_d) = {lam_b - lam_d:.6f} = 8*phi - 4 = 4*(2phi-1) = 4*sqrt(5)")
print(f"  4*sqrt(5) = {4*sqrt(5):.6f}")
print(f"  EXACT? {abs(lam_b - lam_d - 4*sqrt(5)) < 1e-10}")

print(f"\nRatio of differences: (lam_t-lam_u)/(lam_b-lam_d) = 6sqrt5/4sqrt5 = 3/2")
print(f"  Computed: {(lam_t-lam_u)/(lam_b-lam_d):.6f}")

# ============================================================
# 9. COMPREHENSIVE PATTERN SEARCH
# ============================================================
print("\n" + "=" * 70)
print("PATTERN SEARCH: combinations giving {3, 7, 12}")
print("=" * 70)

# Collect all spectral quantities
spectral = {
    'lam_u': lam_u, 'lam_d': lam_d, 'lam_s': lam_s,
    'lam_c': lam_c, 'lam_t': lam_t, 'lam_b': lam_b,
    'n_ut': 4, 'n_db': 2,  # exact log-ratios
    'phi': phi, 'a1': a1, 'b1': b1, 'N_gen': 3,
    'degree': 12, 'N_eig': 9
}

# Test: n_ij = a * n_up + b * n_down for CKM exponents
# n_up = {n_uc, n_ct, n_ut}, n_down = {n_ds, n_sb, n_db}
# We know n_ut = 4 (exact), n_db = 2 (exact)
# n_uc = 3.762, n_ct = 0.238, n_ds = 1.017, n_sb = 0.983

print(f"\nExact spectral exponents:")
print(f"  n_ut = log_phi(lam_t/lam_u) = 4  [EXACT: phi^4]")
print(f"  n_db = log_phi(lam_b/lam_d) = 2  [EXACT: phi^2]")
print(f"  n_uc = log_phi(lam_c/lam_u) = {n_u_c:.6f}")
print(f"  n_ct = log_phi(lam_t/lam_c) = {n_c_t:.6f}")
print(f"  n_ds = log_phi(lam_s/lam_d) = {n_d_s:.6f}")
print(f"  n_sb = log_phi(lam_b/lam_s) = {n_s_b:.6f}")

print(f"\nNote: n_uc + n_ct = {n_u_c + n_c_t:.6f} = n_ut = 4 [EXACT]")
print(f"Note: n_ds + n_sb = {n_d_s + n_s_b:.6f} = n_db = 2 [EXACT]")

# The CKM exponents are n_12=3, n_23=7, n_13=12
# n_ut=4, n_db=2 -> n_ut*n_db = 8 (not 12)
# n_ut + n_db = 6 (not 12)
# n_ut * n_db * N_gen = 24 (not 12)
# 2*n_ut + n_db = 10 = 2*a1 = n_12 + n_23 (I1!)
# n_ut * n_db + n_ut = 12 = n_13!
# n_ut * n_db - n_ut + n_db + 1 = 8-4+2+1 = 7 = n_23!

print("\n--- Testing algebraic combinations of n_ut=4, n_db=2 ---")
n_ut_ex = 4; n_db_ex = 2
tests = [
    ("n_ut + n_db", n_ut_ex + n_db_ex, "6"),
    ("2*n_ut + n_db", 2*n_ut_ex + n_db_ex, "10 = 2a1 = I1!"),
    ("n_ut * n_db", n_ut_ex * n_db_ex, "8"),
    ("n_ut * n_db + n_ut", n_ut_ex * n_db_ex + n_ut_ex, "12 = n_13!"),
    ("n_ut * n_db - 1", n_ut_ex * n_db_ex - 1, "7 = n_23!"),
    ("n_ut - 1", n_ut_ex - 1, "3 = n_12!"),
    ("n_ut^2 - n_ut", n_ut_ex**2 - n_ut_ex, "12 = n_13!"),
]
for name, val, comment in tests:
    match = ""
    if val in [3, 7, 12]:
        match = " *** MATCH ***"
    print(f"  {name:30s} = {val:5d}   ({comment}){match}")

# ============================================================
# 10. ASSESSMENT
# ============================================================
print("\n" + "=" * 70)
print("ASSESSMENT")
print("=" * 70)
print(f"""
TWO EXACT RESULTS on the Cayley graph of 2I:

  lambda_t / lambda_u = phi^4   [EXACT, algebraic proof]
  lambda_b / lambda_d = phi^2   [EXACT, algebraic proof]

These are the spectral ratios of the GEN-1 <-> GEN-3 quark eigenvalues.

From n_ut = 4 and n_db = 2 alone:
  n_12 = n_ut - 1     = 4 - 1 = 3   [MATCH]
  n_23 = n_ut*n_db - 1 = 8 - 1 = 7   [MATCH]
  n_13 = n_ut*n_db + n_ut = 8+4 = 12 [MATCH]

But also:
  n_12 = n_ut - 1 = 3
  n_23 = n_ut^2 - n_ut - 1 = 16-4-1 = 11 [NO]

Simplest consistent system:
  n_12 = n_ut - 1 = 3
  n_23 = n_ut * n_db - 1 = 7
  n_13 = n_ut * n_db + n_ut = n_ut(n_db + 1) = 4*3 = 12

Verification against concordance identities:
  I1: n_12 + n_23 = 3+7 = 10 = 2*a1   [OK]
  I2: n_12 * n_13 = 3*12 = 36 = b1^2   [OK]
  I3: n_23 + n_13 = 7+12 = 19 = a1^2-a1-1 [OK]

STATUS: The Cayley graph eigenvalue ratios phi^4 (up-type) and phi^2
(down-type) encode n_ut=4 and n_db=2. The CKM exponents {{3,7,12}}
can be expressed as:
  n_12 = n_ut - 1
  n_23 = n_ut*n_db - 1
  n_13 = n_ut*(n_db + 1)

However, the map (4,2) -> (3,7,12) involves the "-1" offset and
specific products. This is a PATTERN: the algebra is correct but
the physical mechanism selecting these particular combinations is
NOT derived from the spectral action or any Lagrangian principle.

CLASSIFICATION: PATTERN (promising structural connection, not derivation)
""")

# Verify: n_ut = 4 is itself derived
print("DERIVATION STATUS of n_ut = 4:")
print(f"  lambda_t/lambda_u = phi^4 follows from:")
print(f"  lambda_t = 6+6phi = 6*phi^2 (since 1+phi=phi^2)")
print(f"  lambda_u = 12-6phi = 6*(2-phi) = 6/phi^2 (since 2-phi=1/phi^2)")
print(f"  Ratio = phi^4. QED.")
print(f"  The exponent 4 = a1-1 = d_ST (spacetime dimension).")
print()
print(f"  lambda_b/lambda_d = phi^2 follows from:")
print(f"  lambda_b = 8+4phi = 2(4+2phi) = 2(5+sqrt5)")
print(f"  lambda_d = 12-4phi = 2(6-2phi) = 2(5-sqrt5)")
print(f"  Ratio = (5+sqrt5)/(5-sqrt5) = phi^2. QED.")
print(f"  The exponent 2 = n_db = d_ST/2.")
print()
print(f"  So: n_ut = d_ST, n_db = d_ST/2.")
print(f"  n_12 = d_ST - 1 = N_gen")
print(f"  n_23 = d_ST*(d_ST/2) - 1 = d_ST^2/2 - 1 = 7")
print(f"  n_13 = d_ST*(d_ST/2 + 1) = d_ST*(d_ST+2)/2 = 12")
print(f"  = degree of 600-cell (each vertex has 12 neighbors)")
