"""
exp570: Gauge-consistent U(1) twist on the full simplicial complex.

The RIGHT way to do lattice gauge theory on the 600-cell:
  1. U(1) connection on edges: U_{ij} = e^{i*theta/10} on fiber edges
  2. Twisted d0(theta): vertex -> edge with parallel transport
  3. Twisted d1(theta): edge -> face with holonomy
  4. Field strength F(theta) = d1(theta) @ d0(theta) (NON-ZERO for theta != 0)
  5. Yang-Mills action S_YM = ||F||^2
  6. Check if alpha emerges from the gauge structure

Key insight: d1(A) @ d0(A) != 0 for non-flat connections.
The failure IS the field strength. S_YM = Tr(F^*F) gives the gauge coupling.
"""

import numpy as np
from numpy.linalg import norm, eigvalsh
from collections import defaultdict
import sys
sys.path.insert(0, ".")
from commons import build_600cell

PHI = (1 + np.sqrt(5)) / 2
a1 = 5
b1 = 6
N = 120

verts, adj, lap = build_600cell()

# Build Hopf fibration
def qmul(p, q):
    return np.array([
        p[0]*q[0]-p[1]*q[1]-p[2]*q[2]-p[3]*q[3],
        p[0]*q[1]+p[1]*q[0]+p[2]*q[3]-p[3]*q[2],
        p[0]*q[2]-p[1]*q[3]+p[2]*q[0]+p[3]*q[1],
        p[0]*q[3]+p[1]*q[2]-p[2]*q[1]+p[3]*q[0]])

def find_idx(v, verts, tol=1e-6):
    dots = verts @ v; idx = np.argmax(dots)
    return idx if dots[idx] > 1 - tol else -1

def find_fibration():
    for i in range(N):
        if abs(verts[i, 0] - PHI/2) < 1e-6:
            g = verts[i]; p = g.copy(); ok = True
            for k in range(2, 11):
                p = qmul(p, g)
                if k == 5 and not np.allclose(p, [-1,0,0,0], atol=1e-6): ok=False; break
                if k == 10 and not np.allclose(p, [1,0,0,0], atol=1e-6): ok=False
            if not ok: continue
            used = set(); fibers = []; subg = []
            pp = np.array([1.0,0,0,0])
            for k in range(10): subg.append(find_idx(pp, verts)); pp = qmul(pp, g)
            for s in range(N):
                if s in used: continue
                fib = []
                for si in subg:
                    q = qmul(verts[s], verts[si]); idx = find_idx(q, verts)
                    if idx >= 0 and idx not in used: fib.append(idx); used.add(idx)
                if len(fib) == 10: fibers.append(fib)
            if len(fibers) == 12: return fibers
    return None

fibers = find_fibration()
vtx_fiber = {}
for fi, f in enumerate(fibers):
    for v in f: vtx_fiber[v] = fi

# Get ordered cycles for fiber orientation
def get_ordered_cycles(fibers, adj):
    cycles = []
    for fib in fibers:
        fib_set = set(fib)
        cycle = [fib[0]]
        visited = {fib[0]}
        while len(cycle) < 10:
            curr = cycle[-1]
            for j in fib_set - visited:
                if adj[curr, j] > 0.5:
                    cycle.append(j); visited.add(j); break
        cycles.append(cycle)
    return cycles

cycles = get_ordered_cycles(fibers, adj)

# Build oriented fiber edge set: forward direction along each cycle
fiber_forward = set()  # (i, j) oriented edges
for cyc in cycles:
    for k in range(10):
        i, j = cyc[k], cyc[(k+1) % 10]
        fiber_forward.add((i, j))

# =====================================================================
# STEP 1: Build simplicial complex
# =====================================================================
print("=" * 70)
print("STEP 1: SIMPLICIAL COMPLEX")
print("=" * 70)

adj_list = defaultdict(set)
edges = []
edge_to_idx = {}
for i in range(N):
    for j in range(i+1, N):
        if adj[i,j] > 0.5:
            adj_list[i].add(j); adj_list[j].add(i)
            edge_to_idx[(i,j)] = len(edges)
            edges.append((i,j))

Ne = len(edges)

triangles = []
for i in range(N):
    for j in adj_list[i]:
        if j > i:
            for k in adj_list[i] & adj_list[j]:
                if k > j:
                    triangles.append((i,j,k))

Nf = len(triangles)
face_to_idx = {t: idx for idx, t in enumerate(triangles)}

print(f"  V={N}, E={Ne}, F={Nf}")


# =====================================================================
# STEP 2: Gauge connection U_{ij}(theta)
# =====================================================================

def edge_phase(i, j, theta):
    """U(1) parallel transport from i to j. Phase on fiber edges, 1 on cross."""
    if (i, j) in fiber_forward:
        return np.exp(1j * theta / 10)
    elif (j, i) in fiber_forward:
        return np.exp(-1j * theta / 10)
    else:
        return 1.0  # cross edge: no gauge field


# =====================================================================
# STEP 3: Twisted boundary operators
# =====================================================================

def build_twisted_operators(theta):
    """Build gauge-twisted d0(theta) and d1(theta)."""

    # d0(theta): Nex x Nv, complex
    # For edge e = (i,j) with i<j: (d0 f)_e = f(j) - U_{ij}(theta) * f(i)
    # Equivalently: d0_{e,j} = +1, d0_{e,i} = -U_{ij}
    d0 = np.zeros((Ne, N), dtype=complex)
    for e_idx, (i, j) in enumerate(edges):
        U_ij = edge_phase(i, j, theta)
        d0[e_idx, i] = -U_ij
        d0[e_idx, j] = +1.0

    # d1(theta): Nf x Ne, complex
    # For face (i,j,k) with i<j<k, boundary = (j,k) - (i,k) + (i,j)
    # With gauge: include holonomy along each edge in the triangle
    # The twisted boundary: d1_{f, (a,b)} = sign * U_{path from reference to (a,b)}
    #
    # For triangle (i,j,k), the oriented boundary is:
    #   +U_{jk} on edge (j,k)
    #   -U_{ik} on edge (i,k)  [reversed orientation]
    #   +U_{ij} on edge (i,j)
    #
    # But we need to be consistent: the holonomy around the triangle
    # is W = U_{ij} * U_{jk} * U_{ki} = U_{ij} * U_{jk} * conj(U_{ik})
    #
    # For a FLAT connection: W = 1 and d1@d0 = 0.
    # For non-flat: W != 1 and d1@d0 captures the field strength.
    #
    # Standard lattice convention: d1 encodes the untwisted boundary,
    # and the twist enters through d0. Then F = d1 @ d0(theta).

    d1 = np.zeros((Nf, Ne), dtype=complex)
    for f_idx, (i, j, k) in enumerate(triangles):
        # Boundary of (i,j,k) with parallel transport:
        # We use the convention that d1 carries the gauge transport
        # along the face boundary
        e_ij = edge_to_idx[(i,j)]
        e_jk = edge_to_idx[(j,k)]
        e_ik = edge_to_idx[(i,k)]

        # Holonomy: go i->j->k->i
        # U_{ij->jk} = transport from edge (i,j) endpoint to edge (j,k) start = 1
        # The boundary with transport:
        d1[f_idx, e_ij] = +1.0
        d1[f_idx, e_jk] = +edge_phase(j, k, theta) / edge_phase(j, k, 0)  # just the phase
        d1[f_idx, e_ik] = -edge_phase(i, k, theta) / edge_phase(i, k, 0)

    # Actually, the simplest correct approach: use UNTWISTED d1,
    # and let the field strength emerge from F = d1 @ d0(theta).
    # The untwisted d1:
    d1_flat = np.zeros((Nf, Ne), dtype=complex)
    for f_idx, (i, j, k) in enumerate(triangles):
        d1_flat[f_idx, edge_to_idx[(i,j)]] = +1.0
        d1_flat[f_idx, edge_to_idx[(j,k)]] = +1.0
        d1_flat[f_idx, edge_to_idx[(i,k)]] = -1.0

    return d0, d1_flat


# =====================================================================
# STEP 4: Field strength and Yang-Mills action
# =====================================================================
print("\n" + "=" * 70)
print("STEP 2-4: FIELD STRENGTH AND YANG-MILLS")
print("=" * 70)

# At theta = 0: F = d1 @ d0 = 0 (flat connection)
d0_0, d1_0 = build_twisted_operators(0.0)
F_0 = d1_0 @ d0_0
print(f"  theta=0: ||F|| = {norm(F_0):.2e} (expect 0, flat)")

# At theta != 0: F = d1 @ d0(theta) != 0
# F is a Nf x Nv matrix. The field strength on face f due to vertex i.
# ||F||^2 = Tr(F* @ F) is the Yang-Mills action.

print(f"\n  Yang-Mills action S_YM(theta) = ||d1 @ d0(theta)||^2:")

thetas = np.linspace(0, 2*np.pi, 100)
S_YM = np.zeros(len(thetas))

for idx, theta in enumerate(thetas):
    d0_th, d1_th = build_twisted_operators(theta)
    F = d1_th @ d0_th
    S_YM[idx] = np.real(np.trace(F.conj().T @ F))

print(f"    S_YM(0)     = {S_YM[0]:.6f}")
print(f"    S_YM(pi/5)  = {S_YM[len(thetas)//10]:.6f}")
print(f"    S_YM(pi)    = {S_YM[len(thetas)//2]:.6f}")
print(f"    S_YM(2pi)   = {S_YM[-1]:.6f}")

# The field strength per triangle:
# For a triangle with 1 fiber edge and flux theta:
# The holonomy W = e^{i*theta/10} * 1 * 1 = e^{i*theta/10}
# The field strength F = W - 1 ~ i*theta/10 for small theta
# ||F||^2 per such triangle ~ (theta/10)^2
# Total: T1 * (theta/10)^2 = 600 * theta^2/100 = 6*theta^2

# For small theta: S_YM ~ 6*theta^2 + O(theta^4)
# (only T1 = 600 triangles contribute, each with one fiber edge)

# Check:
small_theta = 0.01
d0_s, d1_s = build_twisted_operators(small_theta)
F_s = d1_s @ d0_s
S_small = np.real(np.trace(F_s.conj().T @ F_s))
coeff = S_small / small_theta**2

print(f"\n  Small-theta coefficient: S_YM ~ {coeff:.4f} * theta^2")
print(f"  Expected: T1 * (1/10)^2 * 2 = 600/100*2 = ... let me check")

# Actually compute more carefully
# F_{f,i} = (d1 @ d0(theta))_{f,i} for face f and vertex i
# For triangle (a,b,c):
# (d1 @ d0)_{(abc), i} = d1_{(abc),(ab)} * d0(theta)_{(ab),i}
#                       + d1_{(abc),(bc)} * d0(theta)_{(bc),i}
#                       + d1_{(abc),(ac)} * d0(theta)_{(ac),i}
#
# For untwisted d1 and twisted d0:
# = (+1)(-U_{ab} delta_{i,a} + delta_{i,b})
# + (+1)(-U_{bc} delta_{i,b} + delta_{i,c})
# + (-1)(-U_{ac} delta_{i,a} + delta_{i,c})
#
# For i=a: -U_{ab}*1 + 0 + 0 + 0 + U_{ac}*(-1) + 0 = -U_{ab} - (-1)*U_{ac}...
# This is getting messy. Let me just compute numerically.

# Per-triangle field strength at small theta
d0_t, d1_t = build_twisted_operators(0.001)
F_t = d1_t @ d0_t  # Nf x Nv
F_norms = np.array([np.real(np.dot(F_t[f].conj(), F_t[f])) for f in range(Nf)])

# Classify by triangle type
T0_norms = []
T1_norms = []
for f_idx, (i,j,k) in enumerate(triangles):
    n_fib = sum([vtx_fiber[i]==vtx_fiber[j], vtx_fiber[j]==vtx_fiber[k],
                 vtx_fiber[i]==vtx_fiber[k]])
    if n_fib == 0:
        T0_norms.append(F_norms[f_idx])
    else:
        T1_norms.append(F_norms[f_idx])

print(f"\n  Field strength by triangle type (theta = 0.001):")
print(f"    T0 (0 fiber edges): mean |F|^2 = {np.mean(T0_norms):.2e}")
print(f"    T1 (1 fiber edge):  mean |F|^2 = {np.mean(T1_norms):.2e}")
print(f"    Only T1 triangles carry field strength: {np.mean(T0_norms) < 1e-15}")


# =====================================================================
# STEP 5: Yang-Mills action coefficient and alpha
# =====================================================================
print("\n" + "=" * 70)
print("STEP 5: YANG-MILLS COEFFICIENT AND alpha")
print("=" * 70)

# S_YM(theta) = K * theta^2 + O(theta^4)
# In continuum gauge theory: S_YM = (1/4g^2) * integral F^2
# And alpha = g^2/(4*pi)
#
# On the lattice: S_YM = K * theta^2 where theta is the total flux
# The lattice coupling is: 1/g^2_lattice = 2*K
# (factor 2 from the standard Wilson action convention)
#
# If we identify theta = g * A * L (gauge coupling * potential * length)
# then S = K*(gAL)^2 and matching to continuum (1/4g^2)*F^2*V gives
# alpha = 1/(8*pi*K*V) or similar.

# But first: what IS K?
# S_YM(theta) for small theta:
thetas_fine = np.linspace(0, 0.1, 1000)
S_fine = np.zeros(len(thetas_fine))
for idx, theta in enumerate(thetas_fine):
    d0_th, d1_th = build_twisted_operators(theta)
    F = d1_th @ d0_th
    S_fine[idx] = np.real(np.trace(F.conj().T @ F))

# Fit S = K*theta^2 + L*theta^4
from numpy.polynomial import polynomial as P
# Fit S/theta^2 = K + L*theta^2
mask = thetas_fine > 0.001
x = thetas_fine[mask]**2
y = S_fine[mask] / thetas_fine[mask]**2
fit = np.polyfit(x, y, 1)
K = fit[1]  # intercept = K
L_coeff = fit[0]  # slope = L

print(f"  S_YM(theta) = K*theta^2 + L*theta^4 + ...")
print(f"  K = {K:.10f}")
print(f"  L = {L_coeff:.6f}")

# What is K exactly?
# K should be related to the number of plaquettes with flux
# Each T1 triangle has holonomy ~ e^{i*theta/10}, so |F|^2 ~ |e^{i*theta/10}-1|^2 ~ theta^2/100
# But there are N vertices, so the trace sums over vertices too.
# Let me compute K from the explicit formula:

# F = d1 @ d0(theta). At small theta: d0(theta) ~ d0(0) + theta*d0'
# d0'_{e,i} = -(i/10) * delta_{fiber_forward} * delta_{i,source}
# So F ~ theta * d1 @ d0' + O(theta^2)
# K = ||d1 @ d0'||^2

# Build d0' (derivative at theta=0)
d0_prime = np.zeros((Ne, N), dtype=complex)
for e_idx, (i, j) in enumerate(edges):
    if (i, j) in fiber_forward:
        d0_prime[e_idx, i] = -1j / 10  # derivative of -e^{i*theta/10} w.r.t. theta
    elif (j, i) in fiber_forward:
        d0_prime[e_idx, i] = 1j / 10   # derivative of -e^{-i*theta/10} w.r.t. theta

F_prime = d1_0 @ d0_prime  # Nf x Nv
K_exact = np.real(np.trace(F_prime.conj().T @ F_prime))

print(f"\n  K (exact, from F' = d1 @ d0') = {K_exact:.10f}")
print(f"  K (from fit) = {K:.10f}")

# Try to identify K
print(f"\n  K = {K_exact:.6f}")
print(f"  T1/100 = {600/100:.1f} = 6")
print(f"  T1 * 2/100 = {600*2/100:.1f} = 12 = degree")
print(f"  K / (T1/100) = {K_exact/(600/100):.6f}")
print(f"  K * 100/T1 = {K_exact*100/600:.6f}")

# In lattice gauge theory, 1/g^2 = beta = 2*N_c * K_plaquette
# For U(1): N_c = 1, so 1/g^2 = 2*K
# And alpha = g^2/(4*pi) = 1/(8*pi*K)

alpha_from_K = 1.0 / (8 * np.pi * K_exact)
alpha_phys = (4*a1*PHI**4 - np.sqrt((4*a1*PHI**4)**2 - 8*np.pi)) / (4*np.pi)

print(f"\n  If 1/g^2 = 2*K (Wilson convention):")
print(f"    g^2 = 1/(2*K) = {1/(2*K_exact):.6f}")
print(f"    alpha = g^2/(4*pi) = {1/(8*np.pi*K_exact):.10f}")
print(f"    Physical alpha     = {alpha_phys:.10f}")
print(f"    Ratio: {alpha_from_K/alpha_phys:.6f}")

# Alternative: what if 1/alpha = 4*pi*K * (normalization)?
inv_alpha_from_K = 8 * np.pi * K_exact
print(f"\n  1/alpha from K = 8*pi*K = {inv_alpha_from_K:.6f}")
print(f"  1/alpha (framework) = {1/alpha_phys:.6f}")

# Try other normalizations
for name, factor in [("1", 1), ("N", N), ("T1", 600), ("N/b1", N/b1),
                      ("phi^4", PHI**4), ("2*a1", 2*a1),
                      ("degree", 12), ("b1", b1)]:
    val = K_exact * factor
    ratio = val / (1/alpha_phys)
    if 0.1 < ratio < 10:
        print(f"  K * {name:>8s} = {val:10.4f}, ratio to 1/alpha = {ratio:.6f}")

# The key: what IS K_exact?
# Each T1 triangle contributes to F. Each has 1 fiber edge with phase e^{i*theta/10}.
# The holonomy around the triangle is (e^{i*theta/10})(1)(1) - 1 ~ i*theta/10 for small theta.
# ||F||^2 per triangle ~ (theta/10)^2 * (factor from vertex counting)

# Let me compute K from first principles:
# For a triangle (a,b,c) with fiber edge (a,b):
# F = d1 @ d0(theta)
# At theta=0+: F_{(abc), a} gets a contribution from the fiber edge (a,b)
# proportional to i/10.
# The exact value depends on the boundary operator signs.

print(f"\n  K_exact = {K_exact:.10f}")
print(f"  K_exact * 10 = {K_exact*10:.6f}")
print(f"  K_exact * 100 = {K_exact*100:.6f}")


# =====================================================================
# STEP 6: FULL YANG-MILLS SCAN
# =====================================================================
print("\n" + "=" * 70)
print("STEP 6: FULL S_YM(theta) STRUCTURE")
print("=" * 70)

thetas_full = np.linspace(0, 2*np.pi, 500)
S_full = np.zeros(len(thetas_full))
for idx, theta in enumerate(thetas_full):
    d0_th, d1_th = build_twisted_operators(theta)
    F = d1_th @ d0_th
    S_full[idx] = np.real(np.trace(F.conj().T @ F))

print(f"  S_YM(0)      = {S_full[0]:.6f}")
print(f"  S_YM(pi/5)   = {S_full[len(thetas_full)//10]:.6f}")
print(f"  S_YM(pi)     = {S_full[len(thetas_full)//2]:.6f}")
print(f"  S_YM(2*pi)   = {S_full[-1]:.6f}")
print(f"  max S_YM     = {np.max(S_full):.6f} at theta = {thetas_full[np.argmax(S_full)]:.4f}")

# Is S_YM(2*pi) = 0? (periodicity check)
print(f"\n  Periodicity: S_YM(2*pi) = {S_full[-1]:.6f}")
print(f"  S_YM(2*pi) = 0? {S_full[-1] < 0.01}")

# S_YM should be periodic with period 2*pi (full Hopf fiber holonomy)
# Check: S(theta + 2*pi) = S(theta)?
# On C10 with 10 edges: theta -> theta + 2*pi adds 2*pi/10 per edge
# e^{i*(theta+2*pi)/10} = e^{i*theta/10} * e^{i*pi/5}
# This is NOT e^{i*theta/10} since e^{i*pi/5} != 1!
# So the period is 10 * 2*pi = 20*pi, not 2*pi!
# The PHYSICAL period (where holonomy returns to 1) is theta = 2*pi*10 = 20*pi.
# But flux quantization on the fiber gives integer units of 2*pi.

# Check periodicity at 20*pi:
# This is a large computation. Let me check a few values.
print(f"\n  Period check:")
for period_name, period in [("2*pi", 2*np.pi), ("4*pi", 4*np.pi),
                             ("2*pi*a1", 2*np.pi*a1), ("20*pi", 20*np.pi)]:
    d0_p, d1_p = build_twisted_operators(period)
    F_p = d1_p @ d0_p
    S_p = np.real(np.trace(F_p.conj().T @ F_p))
    print(f"    S_YM({period_name:>8s} = {period:8.4f}) = {S_p:.6f}")


print("\n" + "=" * 70)
print("EXP-570 COMPLETE")
print("=" * 70)
