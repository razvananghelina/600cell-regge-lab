"""
exp474: Can 4*a1*phi^4 (the alpha equation linear coefficient) be derived from spectral data?

Goal: Find an exact spectral-action route to the tree-level coupling 1/alpha_0 = 4*a1*phi^4.
"""
import numpy as np

phi = (1 + np.sqrt(5)) / 2
a1 = 5
b1 = 6
N = 120

# Target
target = 4 * a1 * phi**4  # = 137.0820...

print("="*70)
print("  exp474: SPECTRAL DERIVATION OF 4*a1*phi^4")
print("="*70)
print(f"\n  Target: 4*a1*phi^4 = {target:.10f}")
print(f"  = 1/alpha_0 (tree-level inverse coupling)\n")

# ============================================================
# SECTION 1: Spectral action coefficients (simplicial D = d+d*)
# ============================================================
c0 = 2640   # Tr(I) = V+E+F+C
c1 = 14880  # Tr(D^2)
c2 = 55920  # (1/2)*Tr(D^4)
A0, A1, A2 = c0//240, c1//240, c2//240  # 11, 62, 233

print("--- Simplicial spectral coefficients ---")
print(f"  c0={c0}, c1={c1}, c2={c2}")
print(f"  A0={A0}, A1={A1}, A2={A2}")
print(f"  Diophantine: 2*A1^2+1 = {2*A1**2+1}, 3*A0*A2 = {3*A0*A2}")

# ============================================================
# SECTION 2: Icosahedral Laplacian eigenvalues
# ============================================================
L3  = a1 - np.sqrt(a1)   # 5 - sqrt(5) = 2.7639...
L5  = b1                  # 6
L3p = a1 + np.sqrt(a1)   # 5 + sqrt(5) = 7.2360...

print("\n--- Icosahedral Laplacian eigenvalues ---")
print(f"  L(3) = {L3:.10f}")
print(f"  L(5) = {L5:.10f}")
print(f"  L(3')= {L3p:.10f}")
print(f"  L(3)*L(3') = {L3*L3p:.10f} = 4*a1 = {4*a1}")
print(f"  L(3')/L(3) = {L3p/L3:.10f} = phi^2 = {phi**2:.10f}")
print(f"  L(3')^3/L(3) = {L3p**3/L3:.10f} = target = {target:.10f}")

# ============================================================
# SECTION 3: Cayley graph eigenvalues (adjacency matrix of 600-cell)
# ============================================================
# 2I has 9 irreps. Adjacency eigenvalues = 12*chi_std/dim
# chi_std for each irrep at the generating class:
cayley_eigs = {
    'rho0': (1, 12),      # trivial
    'rho1': (2, 6*phi),    # std (dim 2)
    'rho2': (3, 4*phi),    # Sym^2 (dim 3)
    'rho3': (2, 6*phi - 6),# dim 2'
    'rho4': (3, 8-4*phi),  # dim 3' (BRANCH)
    'rho5': (4, 3),        # dim 4
    'rho6': (5, 0),        # dim 5
    'rho7': (4, -3),       # dim 4
    'rho8': (6, -2),       # dim 6
}

# Actually let me use the known adjacency eigenvalues of 600-cell
# They are: 12, 6phi, 6phi', 4phi, 4phi', 3, 0, -2, -3
# with multiplicities 1, 4, 4, 9, 9, 16, 25, 25, 16
adj_eigs = [
    (12, 1),
    (6*phi, 4),
    (6*(1-phi), 4),     # 6*phi' = 6-6phi
    (4*phi, 9),
    (4*(1-phi), 9),     # 4*phi' = 4-4phi
    (3, 16),
    (0, 25),
    (-2, 25),
    (-3, 16),
]

print("\n--- Cayley graph adjacency eigenvalues ---")
total_mult = 0
for val, mult in adj_eigs:
    total_mult += mult
    print(f"  lambda = {val:10.6f}, mult = {mult}")
print(f"  Total multiplicity: {total_mult} (should be {N})")

# ============================================================
# SECTION 4: Systematic search for exact expressions
# ============================================================
print("\n" + "="*70)
print("  SYSTEMATIC SEARCH FOR 4*a1*phi^4")
print("="*70)

results = []

def check(name, value, category="?"):
    err = abs(value - target) / target * 100
    match = "EXACT" if err < 1e-8 else f"{err:.6f}%"
    if err < 0.01:
        results.append((name, value, match, category))
        print(f"  {match:>12s}  {name} = {value:.10f}  [{category}]")

# From icosahedral eigenvalues
check("L(3)*L(3')*phi^4", L3*L3p*phi**4, "STRUCTURAL (current)")
check("L(3')^3/L(3)", L3p**3/L3, "ALGEBRAIC identity")
check("L(3')*L(5)*phi^2", L3p*L5*phi**2, "?")  # = (5+sqrt5)*6*phi^2

# From spectral action coefficients
check("c1/c0 * phi^4 * a1", (c1/c0)*phi**4*a1, "?")
check("c2/c1 * phi^4 * a1", (c2/c1)*phi**4*a1, "?")
check("A2/A0 * phi^2", (A2/A0)*phi**2, "?")
check("A1 * phi^2", A1*phi**2, "?")
check("(A0+A1)/A0 * phi^4", (A0+A1)/A0 * phi**4, "?")
check("c2/(2*N*phi^2)", c2/(2*N*phi**2), "?")

# Product geometry: c_k(M x F) = c_k(M)*Tr(I_F) + ...
# Tr(D_F^2) = 8 = rank(E8)
TrDF2 = 8
check("c1/(2*N) * TrDF2 * phi^4 / a1", c1/(2*N) * TrDF2 * phi**4 / a1, "?")
check("A1 * TrDF2 / A0", A1*TrDF2/A0, "?")

# Kaluza-Klein on Hopf fiber
# Vol(S^1) = 2*pi, coupling = 1/Vol(fiber)
# alpha*alpha' = 1/(2*pi)
# So alpha + alpha' = 4*a1*phi^4 / (2*pi)
# The sum of roots should relate to fiber geometry

# From the coexact spectral gap
lam1_coex = 7 - 4*phi
check("N / lam1_coex", N/lam1_coex, "?")
check("(c1/c0) * N / TrDF2", (c1/c0)*N/TrDF2, "?")

# Deeper: relate to the Seeley-DeWitt identity
# 2*A1^2 + 1 = 3*A0*A2
# Can we extract phi^4 from this?
# Note: A1 = 62, A0 = 11, A2 = 233
# A1^2 = 3844, A0*A2 = 2563, 2*3844+1 = 7689 = 3*2563
check("(2*A1^2+1)/(3*A0) * phi^2/a1", (2*A1**2+1)/(3*A0)*phi**2/a1, "?")
check("A2 * phi^4 / TrDF2", A2*phi**4/TrDF2, "?")

# ============================================================
# SECTION 5: The key question - product geometry coupling
# ============================================================
print("\n" + "="*70)
print("  PRODUCT GEOMETRY APPROACH")
print("="*70)

# In NCG, the gauge coupling on M x F satisfies:
# 1/g^2 = f_0 * Tr(Q^2) where Q is the gauge generator on H_F
# For U(1): Tr(Y^2) over all fermion reps
# For SU(2): Tr(T^2)
# For SU(3): Tr(T_a^2)

# On the 600-cell, the "manifold part" gives c_2 = 55920
# The "finite part" (McKay graph) gives Tr(D_F^2) = 8

# The total coupling should involve c_2 * something_from_F

# In standard NCG: S_YM = f_0 * a(F) * int |F|^2
# where a(F) = Tr(D_F^2) for the appropriate projection

# So: 1/g^2 ∝ f_0 * Tr(proj_G(D_F^2))
# The overall normalization involves c_2 from the manifold part

# Let's try: 1/alpha_0 = c_2 / (something)
check("c_2 / (N * phi^2)", c2/(N*phi**2), "?")
check("c_2 / (2*N * phi^2)", c2/(2*N*phi**2), "?")

# What about c_1^2 / (c_0 * something)?
# c_1^2/c_0 = 14880^2/2640 = 83869.09...
check("c1^2/(c0*N*phi^2)", c1**2/(c0*N*phi**2), "?")

# Let me try: the coupling should be c_2/c_1 times something geometric
# c_2/c_1 = 233/62 = 3.7580...
check("(c2/c1)*a1*phi^4/(A2/A1)", (c2/c1)*a1*phi**4/(A2/A1), "?")

# ============================================================
# SECTION 6: Hopf fiber decomposition approach
# ============================================================
print("\n" + "="*70)
print("  HOPF FIBER DECOMPOSITION")
print("="*70)

# On the 600-cell, the Hopf fibration S^3 -> S^2 decomposes:
# 720 edges = 120 fiber + 600 cross-fiber
# Fiber/cross = 120/600 = 1/a1
# Each vertex: 2 fiber neighbors + 2*a1 = 10 cross-fiber

# The spectral action on the fiber C_10 cycle:
# Tr(L_fiber) = 10*2 = 20 (each vertex has degree 2)
# Tr(L_fiber^2) = 10*2*3 = 60

# Cross-fiber Laplacian:
# Tr(L_cross) = 120*10 = 1200
# DEG_cross = 10 = 2*a1

# The coupling 1/alpha involves the fiber holonomy (2*pi)
# and the cross-fiber eigenvalue product (4*a1*phi^4)

# KEY INSIGHT: On the Hopf fiber, the eigenvalues of L_fiber are
# lambda_k = 2(1 - cos(k*pi/a1)) for k=0,...,a1-1
# Product of nonzero eigenvalues = a1 * tau(C_10) / 10
# For C_10: tau = 10, so product = a1*10/10 = a1 = 5... not right

# Actually for a cycle C_n: product of nonzero eigenvalues = n (by matrix-tree)
# Wait, Kirchhoff: tau(C_n) = n, so (1/n)*prod(lambda_i) = tau = n
# => prod(lambda_i) = n^2
# For C_10: prod = 100

fiber_eigs = [2*(1 - np.cos(k*np.pi/a1)) for k in range(1, a1)]
fiber_prod = np.prod(fiber_eigs)
fiber_sum = np.sum(fiber_eigs)
print(f"  Hopf fiber C_10 nonzero eigenvalues: {[f'{e:.6f}' for e in fiber_eigs]}")
print(f"  Product: {fiber_prod:.6f}")
print(f"  Sum: {fiber_sum:.6f} (should be 2*a1 = {2*a1})")

# The fiber eigenvalues are: 2(1-cos(k*pi/5)) for k=1,2,3,4
# = 1/phi^2, 2+phi, 2-phi+2, phi^2+1... let me compute exactly
# k=1: 2(1-cos(pi/5)) = 2(1-phi/2) = 2-phi = 1/phi^2  ✓
# k=2: 2(1-cos(2pi/5)) = 2(1-(phi-1)/2) wait...
# cos(pi/5) = phi/2, cos(2pi/5) = (phi-1)/2 = 1/(2phi)
# k=2: 2(1-1/(2phi)) = 2 - 1/phi = 2 - (phi-1) = 3-phi = phi^2+1-phi... hmm
# Actually 1/phi = phi-1, so 2-1/phi = 2-(phi-1) = 3-phi

print(f"\n  Fiber eigenvalues (exact):")
print(f"    k=1: 1/phi^2 = {1/phi**2:.10f}, computed: {fiber_eigs[0]:.10f}")
print(f"    k=2: 3-phi  = {3-phi:.10f}, computed: {fiber_eigs[1]:.10f}")
print(f"    k=3: phi+1  = {phi+1:.10f}, computed: {fiber_eigs[2]:.10f}")
print(f"    k=4: phi^2+1= {phi**2+1:.10f}, computed: {fiber_eigs[3]:.10f}")
# Wait: k=3: 2(1-cos(3pi/5)) = 2(1-cos(pi-2pi/5)) = 2(1+cos(2pi/5)) = 2+1/phi = 2+phi-1 = 1+phi = phi^2
# k=4: 2(1-cos(4pi/5)) = 2(1-cos(pi-pi/5)) = 2(1+cos(pi/5)) = 2+phi = phi^2+1...
# Hmm: 2+phi = phi^2 + 1 since phi^2 = phi+1

# Fiber product = (1/phi^2)*(3-phi)*(1+phi)*(2+phi)
# = (1/phi^2)*(3-phi)*phi^2*(2+phi)
# = (3-phi)*(2+phi) = 6+3phi-2phi-phi^2 = 6+phi-phi^2 = 6+phi-(phi+1) = 5
# Wait: (3-phi)(2+phi) = 6+3phi-2phi-phi^2 = 6+phi-(phi+1) = 5
# So fiber product = (1/phi^2) * 5 * phi^2... no

# Let me just compute: product = (1/phi^2)*(3-phi)*phi^2*(2+phi)
# = (3-phi)*(2+phi) = 6+3phi-2phi-phi^2 = 6+phi-phi-1 = 5
# That's 5. But we computed 100 above? Let me check...
# Oh wait, the C_10 cycle has 10 vertices, so 9 nonzero eigenvalues, not 4!
# I only computed k=1,...,4 (the a1-1=4 distinct values)
# The full set for C_10 is k=1,...,9 with eigenvalues 2(1-cos(2*pi*k/10))

fiber_eigs_full = [2*(1 - np.cos(2*np.pi*k/10)) for k in range(1, 10)]
fiber_prod_full = np.prod(fiber_eigs_full)
print(f"\n  Full C_10 cycle eigenvalues (k=1..9): product = {fiber_prod_full:.2f}")
# Kirchhoff: tau(C_10) = 10, prod = 10 * 10 = 100

# So the Hopf fiber contributes product = 10^2 = (2*a1)^2
# Hmm, not directly useful.

# ============================================================
# SECTION 7: The definitive test - eigenvalue-ratio derivation
# ============================================================
print("\n" + "="*70)
print("  EIGENVALUE-RATIO APPROACH")
print("="*70)

# We KNOW that 4*a1*phi^4 = L(3)*L(3') * (L(3')/L(3))^2
# = L(3')^3 / L(3)
# The question is: WHY this specific combination?

# In the spectral action, the Yang-Mills term involves:
# S_YM = f_0 * Tr(F^2) where F is the curvature
# For U(1) on the Hopf fiber: F = d(connection)
# The connection 1-form has components in the L(3), L(5), L(3') eigenspaces
# The tree-level coupling should be the "stiffness" of the gauge connection

# Physical argument:
# The U(1) coupling on a fiber bundle is 1/g^2 = Vol(base)*stiffness(fiber)
# On S^3 -> S^2 with fiber S^1:
# - Stiffness of the fiber direction = L(3)*L(3') = 4*a1 (Galois product)
# - Vol ratio = (L(3')/L(3))^2 = phi^4 (ratio of broken/unbroken sectors SQUARED)

# WHY squared? In the spectral action:
# Tr(F^2) involves |nabla A|^2 where A is in the 3' direction
# The kinetic energy has TWO powers of the eigenvalue ratio

# Alternative: in KK, the coupling is:
# 1/alpha = Vol(total) / Vol(fiber) * (curvature correction)
# Vol(total) ∝ N = 120
# Vol(fiber) ∝ 2*a1 = 10 edges per fiber
# N/(2*a1) = 120/10 = 12 = vertex degree

# But 1/alpha_0 = 4*a1*phi^4 ≈ 137, not 12.
# So KK alone doesn't give it directly.

# Let me try: 1/alpha_0 = (N/(2*a1)) * (L(3')/L(3))^2 * ...
check("12 * phi^4", 12*phi**4, "KK: degree * phi^4")
# 12*phi^4 = 82.25, not right

# What about: from the PRODUCT geometry
# The spectral action on M x F at the c_2 level:
# S_c2 = c_2(M) * f_0 * (gauge kinetic terms with D_F)
# c_2(M) = 55920 for the 600-cell simplicial complex
# But f_0 is the test function moment (free parameter!)

# In the NCG framework, the coupling is:
# 1/g^2 = f_0 * a(F) / pi^2
# where a(F) depends on the finite geometry

# For our framework, a(F) = Tr(D_F^2)/2 = 4? Or something related.

# KEY REALIZATION: f_0 is the ONLY free parameter in the spectral action.
# If we SET f_0 such that 1/alpha = 4*a1*phi^4, then:
# f_0 = pi^2 * 4*a1*phi^4 / a(F)

# But this is FITTING f_0, not deriving alpha!
# Unless we can derive f_0 from the geometry...

# The paper mentions f_4 = N^4/(2*pi) as a pattern.
# Is there a similar expression for f_0?

# Let's check: if f_0 = N^0/(2*pi) = 1/(2*pi), then:
# 1/alpha = (1/(2*pi)) * pi^2 * Tr(D_F^2)/2 = pi/2 * 4 = 2*pi
# That gives 1/alpha = 2*pi ≈ 6.28. Not right.

# What if f_0 involves the icosahedral eigenvalues?
# f_0 = L(3')^2/(2*pi*L(3)) ?
# Then 1/alpha = L(3')^2/(2*pi*L(3)) * pi^2 * 4 = 2*pi*L(3')^2/L(3) * 4
# Hmm, getting complicated.

# ============================================================
# SECTION 8: The Cayley product formula
# ============================================================
print("\n" + "="*70)
print("  CAYLEY PRODUCT FORMULA")
print("="*70)

# From exp452b: Cayley spectral ratios lambda_t/lambda_u = phi^4
# The Cayley eigenvalue for rho_1 (the standard rep, dim 2) is 6*phi
# The Cayley eigenvalue for rho_8 (the branch 3', dim 3) is 6*phi' = 6-6*phi

# rho_2 (dim 3, SU(3) part): eigenvalue = 4*phi  (= 12 - 4*phi... wait)
# Actually adjacency eigenvalue for rho_k is lambda_k = DEG * chi_k(gen)/dim_k
# where DEG = 12

# Let me be precise. For the 600-cell Cayley graph:
# rho_0 (dim 1): lambda = 12
# rho_1 (dim 2): lambda = 6*phi ≈ 9.708
# rho_2 (dim 3): lambda = 12*phi/3 = 4*phi ≈ 6.472
# Wait, this depends on the character values.

# Character of the STANDARD (2-dim) rep at the generating class:
# The generating class of 2I consists of the 12 nearest neighbors
# For rho_1 (standard, dim 2): chi = phi (golden ratio character)
# adjacency eigenvalue = 12 * phi / 2... no.

# Actually for Cayley graphs: A*v_k = lambda_k * v_k
# where lambda_k = sum_{s in generators} chi_k(s) / dim(rho_k) * dim(rho_k)
# = sum_{s} chi_k(s) ... no.

# The adjacency eigenvalue of a Cayley graph G = Cay(Gamma, S) for irrep rho_k is:
# lambda_k = (1/dim(rho_k)) * Tr(sum_{s in S} rho_k(s))
# = |S| * chi_k(class of s) / dim(rho_k) when all generators are in one class
# Wait, that's not right either. Let me look at this from known values.

# Known adjacency eigenvalues of 600-cell:
# 12, 6+6/phi, -(6+6/phi), 2+4/phi, -(2+4/phi), 3, -3, 0, -2
# But in Z[phi]: 6+6/phi = 6+6(phi-1) = 6phi? No: 6/phi = 6(phi-1) = 6phi-6
# So 6+6/phi = 6+6phi-6 = 6phi. OK.

# Let me just use the known spectrum:
# {12, 6phi, 6phi', 4phi, 4phi', 3, 0, -2, -3}
# = {12, 9.708, -3.708, 6.472, -0.472, 3, 0, -2, -3}

# The Galois-conjugate pairs:
# (6phi, 6phi') with phi' = 1-phi = -0.618
# (4phi, 4phi')

# Tree-level coupling involves L(3)*L(3') from the ICOSAHEDRAL Laplacian
# which is the VERTEX FIGURE, not the Cayley graph itself.

# But the icosahedral eigenvalues and Cayley eigenvalues are related:
# Laplacian eigenvalue = DEG - adjacency eigenvalue
# For icosahedron (DEG = a1 = 5):
# L(3) = 5 - (adjacency eigenvalue of 3-dim irrep of A_5)
# The icosahedron is the Cayley graph of A_5 with generators = rotations by 2pi/5

# OK this is getting tangled. Let me just state what we know:

print("\n  KNOWN IDENTITIES:")
print(f"  4*a1*phi^4 = L(3)*L(3')*phi^4 = L(3')^3/L(3)")
print(f"  = {target:.10f}")
print()

# The cleanest spectral expression:
# From the Seeley-DeWitt coefficients and the Laplacian:
#
# c_1 = Tr(D^2) = sum over all p-simplices of local trace
# Per vertex: c_1/V = 14880/120 = 124 = dim(E8)/2
#
# c_1/V = 124. And 1/alpha_0 = 4*a1*phi^4 = 137.08
# Ratio: 137.08/124 = 1.1054
# = phi^4/phi^2*something? phi^2 = 2.618, phi^4/phi^2 = phi^2...

# IMPORTANT NEW IDEA:
# c_1 per vertex = 124 = dim(E8)/2
# 1/alpha_0 per "gauge channel" = ?
#
# In KK: 1/alpha = 1/alpha_0 - 2*pi*alpha (self-consistency)
# So 1/alpha_0 is the BARE coupling.
# The bare coupling should be: c_2 / (c_1 * normalization)?

# Let me try: c_2/c_1 * (something)
# c_2/c_1 = 55920/14880 = 233/62 = 3.7580...
# 3.758 * phi^4 / (something)?

# c_2*phi^4/(c_1*phi^2) = c_2*phi^2/c_1 = 55920*2.618/14880 = 9.839
# Not useful.

# What about from the reduced coefficients?
# A1 = 62, and 62*phi^2 = 162.3... not target
# A2/phi^2 = 233/2.618 = 89.00! = F(11)!
check("A2/phi^2", A2/phi**2, "Fibonacci?")
# That's 89, a Fibonacci number! But not our target.

# A0*A2/A1 = 11*233/62 = 41.338...
# phi^4 * A0*A2/A1 = 6.854 * 41.338 = 283.3... not useful

# What about cube of phi^2 * something simple?
# phi^6 = phi^4*phi^2 = 6.854*2.618 = 17.944 = 4*phi^3+4 = 4*4.236+4 = 20.944... wait
# phi^6 = phi^4*phi^2 = (phi+1)^2*(phi+1) = ... let me just compute
print(f"\n  phi^6 = {phi**6:.10f}")
print(f"  4*a1 = 20")
print(f"  target/20 = phi^4 = {phi**4:.10f}")

# So 4*a1*phi^4 = 20*phi^4 = 20*(phi+1)^2 = 20*(phi^2+2phi+1) = 20*(3+2phi+1)
# Hmm wait, phi+1 = phi^2, so phi^4 = (phi^2)^2 = (phi+1)^2 = phi^2+2phi+1 = (phi+1)+2phi+1 = 3phi+2
# Actually phi^4 = 3phi+2 is the Z[phi] expression.
# So 4*a1*phi^4 = 20*(3phi+2) = 60phi+40 in Z[phi].

# In terms of framework constants:
# = 60*phi + 40 = 12*a1*phi + 8*a1 = 4*a1*(3*phi+2) = 4*a1*phi^4
# That's circular.

# ============================================================
# FINAL ASSESSMENT
# ============================================================
print("\n" + "="*70)
print("  ASSESSMENT")
print("="*70)

print("""
  The identity 4*a1*phi^4 = L(3)*L(3') * (L(3')/L(3))^2 = L(3')^3/L(3)
  is an ALGEBRAIC FACT about the icosahedral Laplacian eigenvalues.

  The PHYSICAL question is: why does 1/alpha_0 equal this specific combination?

  Approaches tested:
  1. Spectral action coefficients (c0,c1,c2): These are INTEGERS, while
     4*a1*phi^4 is IRRATIONAL. No rational combination of c_k gives the target.
     STATUS: IMPOSSIBLE (by number theory)

  2. Product geometry (M x F): Requires the test function moment f_0,
     which is genuinely free in any NCG model. Alpha depends on f_0.
     STATUS: NEGATIVE (f_0 is free)

  3. Kaluza-Klein on Hopf fiber: The Vieta product alpha*alpha' = 1/(2*pi)
     IS derived. But the SUM (hence individual roots) requires 4*a1*phi^4.
     STATUS: PARTIAL (product derived, sum not)

  4. Eigenvalue ratio: L(3')/L(3) = phi^2 is EXACT algebra for a1=5.
     The SQUARE phi^4 = (phi^2)^2 could reflect the quadratic nature of
     the gauge kinetic term |F|^2 (two powers of the field strength).
     STATUS: STRUCTURAL MOTIVATION (not a theorem)

  CONCLUSION: 4*a1*phi^4 cannot be derived from the spectral action
  coefficients alone because it's irrational while c_k are integers.
  It CAN be derived from the icosahedral Laplacian eigenvalues, but
  the identification of this specific combination as 1/alpha_0 requires
  either:
  (a) A variational principle selecting the eigenvalue product*ratio^2, or
  (b) A Kaluza-Klein argument on the Hopf bundle (which gives the product
      alpha*alpha'=1/(2pi) but not the individual roots).

  The honest status remains: STRUCTURAL.
  The alpha equation coefficients are all DETERMINED by a1=5, but the
  equation itself is a self-consistency relation, not derivable from
  a single action principle without the test function f_0.
""")

# Quick bonus: is there a UNIQUE self-consistency condition?
print("  BONUS: Self-consistency uniqueness")
print(f"  The equation 2*pi*x^2 - 4*a1*phi^4*x + 1 = 0 is the UNIQUE")
print(f"  quadratic in alpha whose coefficients are:")
print(f"    - quadratic: Vol(S^1) = 2*pi (Hopf holonomy)")
print(f"    - linear: L(3)*L(3')*phi^4 = N*phi^4/b1 (Cayley product)")
print(f"    - constant: 1 (normalization/Vieta)")
print(f"  Each coefficient has a distinct geometric origin.")
print(f"  The COMBINATION is forced if all three constraints apply simultaneously.")
print(f"  This is the strongest honest claim: structural self-consistency.")
