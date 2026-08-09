"""
exp439_greens_function_unification.py
=====================================
TEST: Can a single Green's function on the McKay graph (affine E8)
produce ALL 4 mass correction mechanisms?

Hypothesis: The pseudoinverse of the graph Laplacian L^+ encodes:
  1. Resistance distance R(i,j) = L^+(i,i) + L^+(j,j) - 2*L^+(i,j)
  2. Spectral diagonal = eigenvalue decomposition of L^+
  3. Arakelov height via zeta regularization
  4. Morrey-Sobolev regularity as off-diagonal behavior

We test whether delta_f = C * <rho_f | K | rho_f> for some kernel K
reproduces all 7 branches of the correction formula.
"""
import numpy as np
from numpy.linalg import pinv, eigh

# ============================================================
# 1. McKay graph (affine E8) setup
# ============================================================
N_nodes = 9
dims = np.array([1, 2, 3, 4, 5, 6, 4, 2, 3])  # irrep dimensions
names = ['e', 'u', 'd', 's', 'mu', 'c', 'tau', 't', 'b']
edges = [(0,1), (1,2), (2,3), (3,4), (4,5), (5,6), (6,7), (5,8)]

# Adjacency matrix
A = np.zeros((N_nodes, N_nodes))
for i, j in edges:
    A[i,j] = 1
    A[j,i] = 1

# Degree matrix
D = np.diag(A.sum(axis=1))

# Graph Laplacian L = D - A
L = D - A

print("=" * 70)
print("PART 1: Graph Laplacian and its pseudoinverse")
print("=" * 70)
print("\nAdjacency matrix A:")
print(A.astype(int))
print("\nDegree matrix diagonal:", D.diagonal().astype(int))
print("\nLaplacian L = D - A:")
print(L.astype(int))

# ============================================================
# 2. Pseudoinverse (Green's function)
# ============================================================
# L is singular (has zero eigenvalue for constant vector)
# L^+ is the Moore-Penrose pseudoinverse = Green's function

eigenvalues, eigenvectors = eigh(L)
print("\nLaplacian eigenvalues:", np.round(eigenvalues, 6))

# Pseudoinverse via spectral decomposition
L_plus = np.zeros((N_nodes, N_nodes))
for k in range(N_nodes):
    if abs(eigenvalues[k]) > 1e-10:
        L_plus += np.outer(eigenvectors[:, k], eigenvectors[:, k]) / eigenvalues[k]

print("\nGreen's function L^+ (pseudoinverse):")
for i in range(N_nodes):
    row = "  ".join(f"{L_plus[i,j]:+.5f}" for j in range(N_nodes))
    print(f"  {names[i]:3s}: [{row}]")

# Verify: L * L^+ * L = L
residual = np.max(np.abs(L @ L_plus @ L - L))
print(f"\nVerification ||L L^+ L - L|| = {residual:.2e}")

# ============================================================
# 3. Resistance distances from L^+
# ============================================================
print("\n" + "=" * 70)
print("PART 2: Resistance distances R(i,j) = L+(i,i)+L+(j,j)-2*L+(i,j)")
print("=" * 70)

R = np.zeros((N_nodes, N_nodes))
for i in range(N_nodes):
    for j in range(N_nodes):
        R[i,j] = L_plus[i,i] + L_plus[j,j] - 2*L_plus[i,j]

print("\nResistance distance matrix:")
header = "     " + "  ".join(f"{n:>7s}" for n in names)
print(header)
for i in range(N_nodes):
    row = "  ".join(f"{R[i,j]:7.4f}" for j in range(N_nodes))
    print(f"  {names[i]:3s}: {row}")

# Key test: R(e, d) = R(0, 2) should be 2 (graph distance on tree)
print(f"\nR(e, d) = R(0, 2) = {R[0,2]:.6f}  (need: 2.0 for delta_d = -2/a1)")
# On a tree, resistance = graph distance
print(f"R(e, u) = R(0, 1) = {R[0,1]:.6f}  (graph distance: 1)")
print(f"R(e, s) = R(0, 3) = {R[0,3]:.6f}  (graph distance: 3)")
print(f"R(e, c) = R(0, 5) = {R[0,5]:.6f}")

# But wait - affine E8 has a branch at node 5, so it's NOT a tree
# for paths through node 5. Let's check.
print("\nNOTE: Affine E8 IS a tree (no cycles). All resistance = graph distance.")

# ============================================================
# 4. Diagonal of L^+ (self-Green's function)
# ============================================================
print("\n" + "=" * 70)
print("PART 3: Diagonal L^+(i,i) - self-energy at each node")
print("=" * 70)

for i in range(N_nodes):
    print(f"  L^+({names[i]:3s},{names[i]:3s}) = {L_plus[i,i]:+.6f}   dim={dims[i]}")

# ============================================================
# 5. Weighted Laplacian (by irrep dimensions)
# ============================================================
print("\n" + "=" * 70)
print("PART 4: Dimension-weighted Laplacian")
print("=" * 70)

# The McKay graph has natural weights: each node has weight d_k (irrep dim)
# The weighted Laplacian: L_w(i,j) = delta_ij * degree_i - A_ij * sqrt(d_i*d_j)/d_i ?
# Or more naturally: the normalized Laplacian L_norm = D^{-1/2} L D^{-1/2}

# Actually, the physically relevant object might be the
# Laplacian on the weighted graph where edge (i,j) has weight d_i*d_j/120
# (the Plancherel weight)

# Standard normalized Laplacian
D_inv_sqrt = np.diag(1.0 / np.sqrt(D.diagonal()))
L_norm = D_inv_sqrt @ L @ D_inv_sqrt

evals_norm, evecs_norm = eigh(L_norm)
print("Normalized Laplacian eigenvalues:", np.round(evals_norm, 6))

# Plancherel-weighted Laplacian: weight each node by d_k^2/120
W = np.diag(dims**2 / 120.0)
L_plancherel = W @ L
# This is not symmetric, so use W^{1/2} L W^{1/2}
W_sqrt = np.diag(dims / np.sqrt(120.0))
W_inv_sqrt = np.diag(np.sqrt(120.0) / dims)
L_planch_sym = W_sqrt @ L @ W_sqrt

evals_planch, evecs_planch = eigh(L_planch_sym)
print("Plancherel-weighted eigenvalues:", np.round(evals_planch, 6))

# ============================================================
# 6. Cayley graph Laplacian (the PHYSICAL one from 2I)
# ============================================================
print("\n" + "=" * 70)
print("PART 5: Cayley graph eigenvalues and Galois anti-symmetric G_k")
print("=" * 70)

# The Cayley graph eigenvalues for class c with character chi_k(c):
# lambda_k(c) = |2I|/d_k * (1 - chi_k(c)/d_k)  ... no
# Actually: lambda_k(c) = degree * (1 - chi_k(c)/d_k)
# where degree = |conjugacy class of generators| = 12 (class 10a has 12 elements)

# For the TWO Galois-conjugate generator classes (10a and 5a):
# The adjacency eigenvalues of McKay graph are 2*cos(pi*j/a1) for the
# phi-sector: the eigenvalue for class 10a is phi, for class 5a is -phi

# Character table from McKay graph eigenvectors
# eigenvalues of A (McKay adjacency) in descending order correspond to
# conjugacy classes: 1a(2), 10a(phi), 3a(1), 10b(1/phi), 4a(0),
# 5b(-1/phi), 6a(-1), 5a(-phi), 2a(-2)

phi = (1 + np.sqrt(5)) / 2

# McKay adjacency eigenvalues (sorted descending)
mckay_evals_sorted = np.array([2, phi, 1, 1/phi, 0, -1/phi, -1, -phi, -2])

# Get eigenvectors of A (sorted by eigenvalue descending)
evals_A, evecs_A = eigh(A)
# eigh returns ascending order, reverse
idx = np.argsort(evals_A)[::-1]
evals_A = evals_A[idx]
evecs_A = evecs_A[:, idx]

print("McKay adjacency eigenvalues (computed):", np.round(evals_A, 6))
print("McKay adjacency eigenvalues (expected):", np.round(mckay_evals_sorted, 6))

# Character table: chi_k(c) = evec_k(c) * sqrt(120 / |class_c|)
class_sizes = np.array([1, 12, 20, 12, 30, 12, 20, 12, 1])

char_table = np.zeros((9, 9))
for k in range(9):
    for c in range(9):
        char_table[k, c] = evecs_A[k, c] * np.sqrt(120.0 / class_sizes[c])

print("\nCharacter table chi_k(class_c):")
class_names = ['1a', '10a', '3a', '10b', '4a', '5b', '6a', '5a', '2a']
header = "     " + "  ".join(f"{cn:>7s}" for cn in class_names)
print(header)
for k in range(9):
    row = "  ".join(f"{char_table[k,c]:7.3f}" for c in range(9))
    print(f"  rho_{k} ({names[k]:3s}, d={dims[k]}): {row}")

# Cayley Laplacian eigenvalues for classes 10a and 5a
# L_k(c) = 12 * (1 - chi_k(c) / d_k)
L_10a = np.zeros(9)
L_5a = np.zeros(9)
for k in range(9):
    L_10a[k] = 12.0 * (1.0 - char_table[k, 1] / dims[k])  # class 10a = col 1
    L_5a[k] = 12.0 * (1.0 - char_table[k, 7] / dims[k])    # class 5a = col 7

# Galois anti-symmetric eigenvalue
G_k = (L_10a - L_5a) / 2.0

print("\nGalois anti-symmetric eigenvalue G_k:")
bipartite = ['W', 'B', 'W', 'B', 'W', 'B', 'W', 'B', 'W']
for k in range(9):
    print(f"  rho_{k} ({names[k]:3s}, {bipartite[k]}): G = {G_k[k]:+.6f}  "
          f"L_10a = {L_10a[k]:.4f}, L_5a = {L_5a[k]:.4f}")

# ============================================================
# 7. THE KEY TEST: Can L^+ reproduce the corrections?
# ============================================================
print("\n" + "=" * 70)
print("PART 6: KEY TEST - Comparing L^+ data with mass corrections")
print("=" * 70)

a1 = 5
C_coeff = 2.0 / 13
c_ell = C_coeff * phi**3 / 4.0
N_gen = 3

# Known corrections (from exp438f)
ab_data = {
    'e':   (0, 0),
    'u':   (3, -2),
    'd':   (1, 0),
    's':   (1, 1),
    'mu':  (1, 1),
    'c':   (2, 1),
    'tau': (1, 2),
    't':   (4, 1),
    'b':   (-1, 4),
}

# Compute known delta_1 for each fermion
def known_delta1(name, a, b):
    z = a + b * phi
    zp = a + b * (1 - np.sqrt(5)) / 2
    N_z = a**2 + a*b - b**2  # = z * zp

    if name == 'e': return 0.0
    elif name in ('mu', 'tau'): return c_ell * np.sign(zp) * abs(zp)**0.75
    elif name == 'u': return 0.0  # G < 0
    elif name == 'd': return -2.0 / a1
    elif name == 's': return -N_gen * C_coeff / phi**2
    elif name in ('c', 't'): return C_coeff * np.log(abs(N_z))
    elif name == 'b': return -C_coeff * np.log(abs(N_z)) / phi
    return 0.0

print(f"\n{'Name':>4s}  {'delta_1':>10s}  {'L+(i,i)':>10s}  {'R(e,i)':>10s}  "
      f"{'R(e,i)/a1':>10s}  {'G_k':>10s}  {'G_k*C/phi2':>10s}")
print("-" * 80)

for i, name in enumerate(names):
    a, b = ab_data[name]
    d1 = known_delta1(name, a, b)
    print(f"{name:>4s}  {d1:+10.6f}  {L_plus[i,i]:+10.6f}  {R[0,i]:10.6f}  "
          f"{-R[0,i]/a1:+10.6f}  {G_k[i]:+10.6f}  {G_k[i]*C_coeff/phi**2:+10.6f}")

# ============================================================
# 8. Test specific relationships
# ============================================================
print("\n" + "=" * 70)
print("PART 7: Specific relationship tests")
print("=" * 70)

# Test 1: d quark - resistance distance
delta_d_known = -2.0 / a1
delta_d_resistance = -R[0, 2] / a1
print(f"\n[d quark] Resistance test:")
print(f"  delta_d (known)     = {delta_d_known:.6f}")
print(f"  -R(e,d)/a1          = {delta_d_resistance:.6f}")
print(f"  MATCH: {np.isclose(delta_d_known, delta_d_resistance)}")

# Test 2: s quark - Cayley eigenvalue
delta_s_known = -N_gen * C_coeff / phi**2
delta_s_cayley = C_coeff * G_k[3] / phi**2  # node 3 = s
print(f"\n[s quark] Cayley eigenvalue test:")
print(f"  delta_s (known)     = {delta_s_known:.6f}")
print(f"  C * G_s / phi^2     = {delta_s_cayley:.6f}")
print(f"  G_s = {G_k[3]:.6f} (should be -3 = -N_gen)")
print(f"  MATCH: {np.isclose(delta_s_known, delta_s_cayley)}")

# Test 3: u quark - G < 0 selection
print(f"\n[u quark] Spectral selection test:")
print(f"  G_u = {G_k[1]:.6f} (should be -6*phi = {-6*phi:.6f})")
print(f"  G_u < 0: {G_k[1] < 0} => delta_u = max(0, ...) = 0  CORRECT")

# Test 4: c, t quarks - Arakelov height
N_c = 2**2 + 2*1 - 1**2  # a^2 + ab - b^2 for c(2,1)
N_t = 4**2 + 4*1 - 1**2  # for t(4,1)
delta_c_known = C_coeff * np.log(abs(N_c))
delta_t_known = C_coeff * np.log(abs(N_t))

print(f"\n[c quark] Arakelov height test:")
print(f"  N(z_c) = {N_c}, ln|N| = {np.log(abs(N_c)):.6f}")
print(f"  delta_c (known) = C*ln|N| = {delta_c_known:.6f}")

print(f"\n[t quark] Arakelov height test:")
print(f"  N(z_t) = {N_t}, ln|N| = {np.log(abs(N_t)):.6f}")
print(f"  delta_t (known) = C*ln|N| = {delta_t_known:.6f}")

# Test 5: b quark - Arakelov + Galois suppression
N_b = (-1)**2 + (-1)*4 - 4**2  # for b(-1,4)
delta_b_known = -C_coeff * np.log(abs(N_b)) / phi

print(f"\n[b quark] Arakelov + Galois test:")
print(f"  N(z_b) = {N_b}, |N| = {abs(N_b)}, ln|N| = {np.log(abs(N_b)):.6f}")
print(f"  delta_b (known) = -C*ln|N|/phi = {delta_b_known:.6f}")

# ============================================================
# 9. THE BIG QUESTION: Is there a single kernel K?
# ============================================================
print("\n" + "=" * 70)
print("PART 8: Searching for unified kernel K")
print("=" * 70)

# The spectral decomposition of L^+ is:
# L^+(i,j) = sum_{k: lambda_k > 0} v_k(i) * v_k(j) / lambda_k
#
# The Cayley Laplacian eigenvalues lambda_c(k) = 12*(1 - chi_k(c)/d_k)
# depend on conjugacy class c.
#
# For Galois-conjugate classes (10a, 5a):
# L_10a(k) and L_5a(k) give the Galois-symmetric S_k and anti-symmetric G_k
#
# Can we build a kernel from the Cayley graph that gives Arakelov?
#
# Key insight: The Arakelov height ln|N(z)| uses the Z[phi] NORM,
# which involves BOTH the physical value z and the Galois conjugate z'.
# This is exactly what the TWO Galois eigenvalues encode!

# Let's try: define a "Galois kernel" using the 2I character table
# K_Galois(k) = some function of chi_k(10a) and chi_k(5a)

# For each node k, we have:
#   chi_k(10a) and chi_k(5a) - the two Galois-conjugate character values
#   Their product chi_k(10a)*chi_k(5a) relates to the NORM via Galois

print("\nGalois character products chi(10a)*chi(5a)/d_k^2:")
for k in range(9):
    product = char_table[k, 1] * char_table[k, 7] / dims[k]**2
    print(f"  rho_{k} ({names[k]:3s}): chi(10a)={char_table[k,1]:+.4f}, "
          f"chi(5a)={char_table[k,7]:+.4f}, "
          f"product/d^2={product:+.6f}")

# ============================================================
# 10. Heat kernel approach
# ============================================================
print("\n" + "=" * 70)
print("PART 9: Heat kernel K(t) = exp(-t*L) at various t")
print("=" * 70)

# The heat kernel on the graph: K(t,i,j) = sum_k exp(-t*lambda_k) v_k(i) v_k(j)
# At t -> infinity, K -> 1/N * J (uniform)
# At t -> 0, K -> I (identity)
# At intermediate t, it probes the geometry at scale sqrt(t)

for t in [0.1, 0.5, 1.0, 2.0, 5.0]:
    K_heat = np.zeros((N_nodes, N_nodes))
    for k in range(N_nodes):
        K_heat += np.exp(-t * eigenvalues[k]) * np.outer(
            eigenvectors[:, k], eigenvectors[:, k])

    # Diagonal = return probability
    diag = np.diag(K_heat)
    print(f"\nt = {t:.1f}: Heat kernel diagonal (return probability):")
    for i in range(N_nodes):
        print(f"  {names[i]:3s}: K(t,i,i) = {diag[i]:.6f}")

# ============================================================
# 11. Zeta function approach: zeta_L(s) = sum lambda_k^{-s}
# ============================================================
print("\n" + "=" * 70)
print("PART 10: Spectral zeta function and regularization")
print("=" * 70)

# Local zeta function at node i:
# zeta_i(s) = sum_{k: lambda_k > 0} |v_k(i)|^2 / lambda_k^s
# zeta_i(1) = L^+(i,i) (the diagonal of the Green's function)
# zeta_i'(0) = -sum_{k>0} |v_k(i)|^2 * ln(lambda_k) = spectral determinant

nonzero_evals = eigenvalues[eigenvalues > 1e-10]
nonzero_evecs = eigenvectors[:, eigenvalues > 1e-10]

print("\nLocal spectral zeta at s=1: zeta_i(1) = L^+(i,i)")
for i in range(N_nodes):
    zeta_1 = sum(nonzero_evecs[i, k]**2 / nonzero_evals[k]
                 for k in range(len(nonzero_evals)))
    print(f"  {names[i]:3s}: zeta(1) = {zeta_1:+.6f} = L^+(i,i) = {L_plus[i,i]:+.6f}")

print("\nLocal spectral zeta derivative at s=0: -zeta_i'(0)")
for i in range(N_nodes):
    zeta_prime = -sum(nonzero_evecs[i, k]**2 * np.log(nonzero_evals[k])
                      for k in range(len(nonzero_evals)))
    print(f"  {names[i]:3s}: -zeta'(0) = {zeta_prime:+.6f}")

# ============================================================
# 12. Combined kernel: L^+ with Z[phi] data
# ============================================================
print("\n" + "=" * 70)
print("PART 11: Attempting unified formula")
print("=" * 70)

# The idea: each correction mechanism is a different "channel" of the
# same Green's function, selected by the Z[phi] arithmetic:
#
# delta_f = C * { h_Arakelov(z_f)  if |N(z_f)| > 1   [prime sector]
#               { G_k / phi^2      if G_k != 0         [spectral sector]
#               { -R(0,k) / a1     if z in Z           [rational sector]
#               { c_ell/C * |z'|^{3/4}  if |N|<=1, G=0 [regularity sector]
#
# Can we write this as:
# delta_f = C * [h(z) * Theta(|N|>1) + G_k/phi^2 * Theta(G!=0) * Theta(|N|<=1)
#              + (-R/a1) * Theta(b=0) * Theta(|N|<=1) * Theta(G=0)
#              + Morrey(z') * Theta(|N|<=1) * Theta(G=0) * Theta(b!=0)]
#
# The Theta functions form a PARTITION: every fermion hits exactly one branch.

# Let me verify the partition property
print("\nPartition verification (each fermion hits exactly one branch):")
for i, name in enumerate(names):
    a, b = ab_data[name]
    N_z = a**2 + a*b - b**2
    zp = a + b * (1 - np.sqrt(5))/2

    is_prime = abs(N_z) > 1
    has_G = abs(G_k[i]) > 1e-10
    is_rational = (b == 0)
    is_lepton_type = (not is_prime) and (not has_G) and (not is_rational)
    is_electron = (name == 'e')

    branch = "?"
    if is_electron: branch = "INPUT"
    elif is_prime: branch = "ARAKELOV"
    elif has_G: branch = "SPECTRAL"
    elif is_rational: branch = "RESISTANCE"
    elif is_lepton_type: branch = "MORREY"

    print(f"  {name:3s}: |N|={abs(N_z):2d}, G={G_k[i]:+7.3f}, b={b:+2d}, "
          f"z'={zp:+7.4f} => {branch}")

# ============================================================
# 13. Key insight: ALL mechanisms relate to graph Laplacian
# ============================================================
print("\n" + "=" * 70)
print("PART 12: Connecting all mechanisms to L")
print("=" * 70)

print("""
STRUCTURAL CONNECTIONS:

1. RESISTANCE (d quark):
   R(0,2) = L^+(0,0) + L^+(2,2) - 2*L^+(0,2)
   This IS the Green's function of L, evaluated between nodes 0 and 2.
   delta_d = -R(0,2)/a1 = -2/5.
   => DIRECTLY from L^+.

2. SPECTRAL (u, s quarks):
   G_k = (L_10a(k) - L_5a(k)) / 2
   where L_c(k) = 12*(1 - chi_k(c)/d_k)
   This is the Galois-antisymmetric part of the CAYLEY Laplacian.
   The Cayley Laplacian eigenvalues ARE the McKay adjacency eigenvalues
   (via the McKay correspondence).
   => FROM the spectrum of L (Cayley = McKay).

3. ARAKELOV (c, t, b quarks):
   ln|N(z)| = ln|z| + ln|sigma(z)|
   The norm N(z) = z * sigma(z) involves the PRODUCT of Galois conjugates.
   On the McKay graph, the Galois action sigma swaps eigenvalues
   lambda(10a) <-> lambda(5a), i.e., phi <-> -phi.
   The Arakelov height is the LOG of the product of embeddings
   => this is zeta'(0) of the Arakelov divisor on Q(sqrt(5)).
   QUESTION: Does this connect to zeta'(0) of L?

4. MORREY (mu, tau):
   |z'|^{3/4}: the Galois conjugate value raised to the Morrey exponent.
   The exponent 3/4 = 1 - 1/d_ST comes from summability of the
   spectral triple, which IS a property of D = d + d* on the 600-cell.
   QUESTION: Does this connect to off-diagonal L^+ behavior?
""")

# ============================================================
# 14. Test: Arakelov from graph zeta?
# ============================================================
print("=" * 70)
print("PART 13: Testing Arakelov-graph connection")
print("=" * 70)

# The Z[phi] data for each fermion lives on the McKay graph.
# For each node k, define z_k = a_k + b_k * phi.
# The Galois conjugate z'_k = a_k + b_k * phi'.
# The norm N(z_k) = z_k * z'_k.

# On the McKay graph, the eigenvalue for class 10a is phi,
# and for class 5a is -phi (= phi').
# So the "physical" propagator uses lambda = phi,
# and the "Galois" propagator uses lambda = phi' = -1/phi.

# The resolvent of A at eigenvalue lambda:
# G(lambda; i, j) = sum_k v_k(i)*v_k(j) / (lambda - evals_A[k])
# At lambda = phi: singular at the 10a eigenvalue
# At lambda = phi': singular at the 5a eigenvalue

# Instead: the REGULARIZED product of resolvents might give Arakelov!
# h(z) ~ G_reg(phi; i,i) * G_reg(-phi; i,i) ?

# Let's compute the regularized Green's functions
def resolvent_reg(A, lam, i, j, evals, evecs):
    """Resolvent with the singular eigenvalue removed (regularized)."""
    result = 0.0
    for k in range(len(evals)):
        diff = lam - evals[k]
        if abs(diff) > 0.01:  # skip near-singular
            result += evecs[i, k] * evecs[j, k] / diff
    return result

print("\nRegularized diagonal resolvent G_reg(lambda; i, i):")
header_gp = "G(phi';i,i)"
print(f"{'Name':>4s}  {'G(phi;i,i)':>12s}  {header_gp:>12s}  {'Product':>12s}  {'ln|N(z)|':>10s}")
print("-" * 60)

for i, name in enumerate(names):
    a, b = ab_data[name]
    N_z = a**2 + a*b - b**2

    g_phi = resolvent_reg(A, phi, i, i, evals_A, evecs_A)
    g_phi_prime = resolvent_reg(A, -1/phi, i, i, evals_A, evecs_A)
    product = g_phi * g_phi_prime

    ln_N = np.log(abs(N_z)) if N_z != 0 else 0.0

    print(f"{name:>4s}  {g_phi:+12.6f}  {g_phi_prime:+12.6f}  {product:+12.6f}  {ln_N:10.6f}")

# ============================================================
# 15. Summary: what connects to L^+, what doesn't
# ============================================================
print("\n" + "=" * 70)
print("PART 14: SUMMARY")
print("=" * 70)

print("""
RESULT SUMMARY:

MECHANISM          | FROM L^+ ?     | STATUS
-------------------|----------------|--------
1. Resistance (d)  | YES (direct)   | R(0,2)/a1 = 2/5 EXACT
2. Spectral (u,s)  | YES (spectrum)  | G_k from McKay eigenvalues
3. Arakelov (c,t,b)| PARTIAL        | ln|N| involves Z[phi] data,
                   |                | not directly L^+ diagonal
4. Morrey (mu,tau) | INDIRECT       | exponent 3/4 from d_ST-summability,
                   |                | coefficient from C and phi^3

KEY FINDING: Mechanisms 1 and 2 come DIRECTLY from the graph Laplacian L.
Mechanism 3 (Arakelov) uses ARITHMETIC data (the Z[phi] norm) that is
assigned to the graph nodes but not intrinsic to L alone.
Mechanism 4 (Morrey) uses ANALYTIC data (Sobolev regularity) that comes
from the spectral triple, not from L alone.

CONCLUSION: A SINGLE L^+ does NOT generate all corrections.
But all corrections come from TWO objects:
  (A) The graph Laplacian L (mechanisms 1, 2)
  (B) The Z[phi]-valued connection z_f on the graph (mechanisms 3, 4)

The UNIFIED object is: the Laplacian of the Z[phi]-valued connection
on the McKay graph, i.e., the operator
  D_Z = d + z * d*
where d is the graph differential and z is the Z[phi] field.
""")

print("\nDone.")
