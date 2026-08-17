"""
exp560: Can the alpha equation emerge from the spectral action on Box?

IDEA: Introduce a U(1) gauge flux Phi on the Hopf fiber.
The modified fiber adjacency gets Peierls phases: A_fiber(Phi).
The spectral action S(Phi) = Tr f(Box(Phi)/Lambda) is a function of Phi.
If dS/dPhi = 0 at a special value related to alpha, the KK identification
becomes a consequence, not a postulate.

APPROACH:
  1. Build Box(Phi) = b1*A_fiber(Phi) - A with phases on fiber edges
  2. Compute eigenvalues for each Phi
  3. Compute S(Phi) for various test functions f
  4. Find extrema and check for alpha

BLIND: compute first, compare after.
"""

import numpy as np
from numpy.linalg import eigvalsh
import sys
sys.path.insert(0, ".")
from commons import build_600cell

PHI_GR = (1 + np.sqrt(5)) / 2  # golden ratio (renamed to avoid confusion with flux)
a1 = 5
b1 = 6
N = 120

verts, adj, lap = build_600cell()

# === Quaternion tools and Hopf fibration ===
def qmul(p, q):
    w = p[0]*q[0] - p[1]*q[1] - p[2]*q[2] - p[3]*q[3]
    x = p[0]*q[1] + p[1]*q[0] + p[2]*q[3] - p[3]*q[2]
    y = p[0]*q[2] - p[1]*q[3] + p[2]*q[0] + p[3]*q[1]
    z = p[0]*q[3] + p[1]*q[2] - p[2]*q[1] + p[3]*q[0]
    return np.array([w, x, y, z])

def find_idx(v, verts, tol=1e-6):
    dots = verts @ v
    idx = np.argmax(dots)
    return idx if dots[idx] > 1 - tol else -1

# Find Hopf fibration with ORDERED fibers (need orientation for gauge field)
def find_ordered_fibration():
    """Find fibration and return fibers as ORDERED cycles."""
    target_w = PHI_GR / 2.0
    for i in range(N):
        if abs(verts[i, 0] - target_w) < 1e-6:
            g = verts[i]
            power = g.copy()
            ok = True
            for k in range(2, 11):
                power = qmul(power, g)
                if k == 5 and not np.allclose(power, [-1,0,0,0], atol=1e-6):
                    ok = False; break
                if k == 10 and not np.allclose(power, [1,0,0,0], atol=1e-6):
                    ok = False
            if not ok:
                continue

            # Build ORDERED fibers using right multiplication by g
            subg_indices = []
            p = np.array([1.0, 0, 0, 0])
            for k in range(10):
                idx = find_idx(p, verts)
                subg_indices.append(idx)
                p = qmul(p, g)

            used = set()
            fibers_ordered = []
            for s in range(N):
                if s in used:
                    continue
                fiber = []
                q = verts[s].copy()
                for k in range(10):
                    idx = find_idx(q, verts)
                    if idx < 0 or idx in used:
                        break
                    fiber.append(idx)
                    used.add(idx)
                    q = qmul(q, g)
                if len(fiber) == 10:
                    fibers_ordered.append(fiber)

            if len(fibers_ordered) == 12 and len(used) == N:
                return fibers_ordered
    return None

print("Building 600-cell and Hopf fibration...")
fibers = find_ordered_fibration()
assert fibers is not None, "Could not find ordered fibration"

# Build UNMODIFIED A_fiber for reference
A_fiber = np.zeros((N, N))
for fib in fibers:
    for i in fib:
        for j in fib:
            if i != j and adj[i,j] > 0.5:
                A_fiber[i,j] = 1.0


# =====================================================================
# Build Box(Phi) with U(1) flux on fiber edges
# =====================================================================

def build_box_with_flux(fibers, adj, flux):
    """
    Build Box(Phi) = b1*A_fiber(Phi) - A where A_fiber has Peierls phases.

    Each fiber is an ordered cycle [v0, v1, ..., v9].
    Edge v_k -> v_{k+1} gets phase e^{i*flux/10} (flux per edge = total flux / 10).
    Edge v_{k+1} -> v_k gets phase e^{-i*flux/10} (Hermitian conjugate).

    Returns: Hermitian matrix Box(Phi), shape (N, N), complex.
    """
    phase_per_edge = flux / 10  # total flux divided by edges per fiber
    A_fiber_phi = np.zeros((N, N), dtype=complex)

    for fib in fibers:
        for k in range(10):
            i = fib[k]
            j = fib[(k + 1) % 10]
            # Check this is actually an edge
            if adj[i, j] > 0.5:
                A_fiber_phi[i, j] = np.exp(1j * phase_per_edge)
                A_fiber_phi[j, i] = np.exp(-1j * phase_per_edge)
            else:
                # Fiber ordering might not match graph edges; try reverse
                # The decagonal fiber has edges between consecutive vertices
                # but our ordering might skip vertices. Use brute force.
                pass

    # Check: should have same nonzero pattern as A_fiber
    nonzero_match = np.sum(np.abs(A_fiber_phi) > 0.5)
    nonzero_orig = np.sum(A_fiber > 0.5)

    if nonzero_match < nonzero_orig:
        # Fallback: find fiber edges directly and assign phases
        A_fiber_phi = np.zeros((N, N), dtype=complex)
        for fib in fibers:
            # Find the cycle ordering from adjacency
            fib_set = set(fib)
            cycle = [fib[0]]
            visited = {fib[0]}
            while len(cycle) < 10:
                curr = cycle[-1]
                found_next = False
                for j in fib_set - visited:
                    if adj[curr, j] > 0.5:
                        cycle.append(j)
                        visited.add(j)
                        found_next = True
                        break
                if not found_next:
                    break

            if len(cycle) == 10:
                for k in range(10):
                    i = cycle[k]
                    j = cycle[(k + 1) % 10]
                    A_fiber_phi[i, j] = np.exp(1j * phase_per_edge)
                    A_fiber_phi[j, i] = np.exp(-1j * phase_per_edge)

    Box_phi = b1 * A_fiber_phi - adj.astype(complex)
    return Box_phi


# =====================================================================
# PART 1: Verify flux=0 reproduces original Box
# =====================================================================
print("\n" + "=" * 70)
print("PART 1: VERIFICATION (flux = 0)")
print("=" * 70)

Box0 = build_box_with_flux(fibers, adj, 0.0)
Box_real = b1 * A_fiber - adj
diff = np.max(np.abs(Box0 - Box_real))
print(f"  ||Box(0) - Box_original|| = {diff:.2e}")

evals0 = np.sort(np.linalg.eigvalsh(Box0))
ker0 = np.sum(np.abs(evals0) < 1e-6)
print(f"  ker(Box(0)) = {ker0} (expect 9)")


# =====================================================================
# PART 2: Spectral action S(Phi) scan
# =====================================================================
print("\n" + "=" * 70)
print("PART 2: SPECTRAL ACTION S(Phi)")
print("=" * 70)

# Scan flux from 0 to 2*pi
n_flux = 200
flux_vals = np.linspace(0, 2*np.pi, n_flux)

# Test functions for spectral action
# f1: heat kernel at t=1: f(x) = e^{-x}  => S = Tr(exp(-Box))
# f2: zeta at s=2: f(x) = |x|^{-2}  => S = sum 1/|mu_i|^2 (exclude zeros)
# f3: quadratic: f(x) = x^2  => S = Tr(Box^2)  -- trivially constant?
# f4: log determinant: f(x) = log|x|  => S = log|det'(Box)|

S_tr2 = np.zeros(n_flux)      # Tr(Box^2)
S_det = np.zeros(n_flux)      # log|det'|
S_eta = np.zeros(n_flux)      # eta invariant: sum sign(mu)
S_ker = np.zeros(n_flux)      # kernel dimension
S_zeta2 = np.zeros(n_flux)    # sum 1/|mu|^2

for i, phi in enumerate(flux_vals):
    Box_phi = build_box_with_flux(fibers, adj, phi)
    evals = np.sort(np.linalg.eigvalsh(Box_phi))

    S_tr2[i] = np.sum(evals**2)
    S_ker[i] = np.sum(np.abs(evals) < 1e-4)

    nonzero = evals[np.abs(evals) > 1e-6]
    if len(nonzero) > 0:
        S_det[i] = np.sum(np.log(np.abs(nonzero)))
        S_eta[i] = np.sum(np.sign(nonzero))
        S_zeta2[i] = np.sum(np.abs(nonzero)**(-2))

# Find extrema of each spectral function
print(f"\n  Scanning {n_flux} flux values in [0, 2*pi]...")
print(f"\n  Tr(Box^2):")
print(f"    min = {np.min(S_tr2):.4f} at Phi = {flux_vals[np.argmin(S_tr2)]:.6f}")
print(f"    max = {np.max(S_tr2):.4f} at Phi = {flux_vals[np.argmax(S_tr2)]:.6f}")
print(f"    range = {np.max(S_tr2) - np.min(S_tr2):.4f}")

print(f"\n  log|det'(Box)|:")
idx_min_det = np.argmin(S_det[S_det != 0])
idx_max_det = np.argmax(S_det[S_det != 0])
print(f"    min = {np.min(S_det[S_det != 0]):.4f}")
print(f"    max = {np.max(S_det[S_det != 0]):.4f}")

# Find local extrema of log|det'|
from scipy.signal import argrelextrema
if len(S_det) > 3:
    local_max = argrelextrema(S_det, np.greater, order=3)[0]
    local_min = argrelextrema(S_det, np.less, order=3)[0]

    print(f"\n  Local maxima of log|det'(Box)|:")
    for idx in local_max[:10]:
        print(f"    Phi = {flux_vals[idx]:.6f} ({flux_vals[idx]/np.pi:.4f}*pi), "
              f"S = {S_det[idx]:.4f}")

    print(f"\n  Local minima of log|det'(Box)|:")
    for idx in local_min[:10]:
        print(f"    Phi = {flux_vals[idx]:.6f} ({flux_vals[idx]/np.pi:.4f}*pi), "
              f"S = {S_det[idx]:.4f}")

# Check kernel dimension
print(f"\n  Kernel dimension vs flux:")
unique_kers = set()
for i, phi in enumerate(flux_vals):
    k = int(S_ker[i])
    if k > 0 or i == 0:
        unique_kers.add((flux_vals[i], k))

for phi, k in sorted(unique_kers)[:20]:
    if k > 0:
        print(f"    Phi = {phi:.6f}: ker = {k}")


# =====================================================================
# PART 3: Fine scan around extrema
# =====================================================================
print("\n" + "=" * 70)
print("PART 3: FINE SCAN -- LOOKING FOR ALPHA")
print("=" * 70)

# The physical alpha: Phi_phys = 2*pi*alpha? Or Phi_phys = alpha * (pi/a1)?
alpha_phys = (4*a1*PHI_GR**4 - np.sqrt((4*a1*PHI_GR**4)**2 - 8*np.pi)) / (4*np.pi)
print(f"\n  Physical alpha = {alpha_phys:.10f} = 1/{1/alpha_phys:.4f}")
print(f"  2*pi*alpha = {2*np.pi*alpha_phys:.10f}")
print(f"  alpha*pi/a1 = {alpha_phys*np.pi/a1:.10f}")
print(f"  alpha*2*pi*a1 = {alpha_phys*2*np.pi*a1:.10f}")

# Fine scan near Phi = 0 (small flux regime)
n_fine = 500
flux_fine = np.linspace(0, 0.5, n_fine)
S_det_fine = np.zeros(n_fine)
S_tr2_fine = np.zeros(n_fine)

for i, phi in enumerate(flux_fine):
    Box_phi = build_box_with_flux(fibers, adj, phi)
    evals = np.sort(np.linalg.eigvalsh(Box_phi))
    nonzero = evals[np.abs(evals) > 1e-6]
    S_det_fine[i] = np.sum(np.log(np.abs(nonzero))) if len(nonzero) > 0 else 0
    S_tr2_fine[i] = np.sum(evals**2)

# Numerical derivative of log|det'|
dS_det = np.gradient(S_det_fine, flux_fine)
# Find zeros of derivative (extrema)
zero_crossings = []
for i in range(len(dS_det) - 1):
    if dS_det[i] * dS_det[i+1] < 0:
        # Linear interpolation
        phi_zero = flux_fine[i] - dS_det[i] * (flux_fine[i+1] - flux_fine[i]) / (dS_det[i+1] - dS_det[i])
        zero_crossings.append(phi_zero)

print(f"\n  Extrema of log|det'(Box)| in [0, 0.5]:")
for phi_z in zero_crossings:
    print(f"    Phi = {phi_z:.10f}")
    print(f"    Phi/(2*pi) = {phi_z/(2*np.pi):.10f}")
    print(f"    Phi * a1/pi = {phi_z * a1/np.pi:.10f}")

# Also check: what is dS/dPhi at Phi = 2*pi*alpha?
phi_alpha = 2*np.pi*alpha_phys
if phi_alpha < 0.5:
    # Find nearest point
    idx_alpha = np.argmin(np.abs(flux_fine - phi_alpha))
    print(f"\n  At Phi = 2*pi*alpha = {phi_alpha:.6f}:")
    print(f"    dS/dPhi = {dS_det[idx_alpha]:.6f}")
    print(f"    S = {S_det_fine[idx_alpha]:.6f}")


# =====================================================================
# PART 4: ANALYTIC APPROACH -- Perturbation theory
# =====================================================================
print("\n" + "=" * 70)
print("PART 4: PERTURBATIVE ANALYSIS (small flux)")
print("=" * 70)

# For small flux Phi, Box(Phi) = Box(0) + Phi*B1 + Phi^2*B2/2 + ...
# where B1 = d(Box)/dPhi|_0 and B2 = d^2(Box)/dPhi^2|_0
#
# B1 = b1 * d(A_fiber)/dPhi|_0
# For A_fiber(Phi): (A_f)_ij = exp(i*Phi/10) on forward edges, exp(-i*Phi/10) on backward
# d(A_f)/dPhi|_0 = (i/10)*A_forward - (i/10)*A_backward = (i/10)*(A_forward - A_backward)
# This is i/10 times the "current operator" J on the fiber

# Build the current operator J = A_forward - A_backward
J = np.zeros((N, N), dtype=complex)
for fib in fibers:
    # Get cycle ordering
    fib_set = set(fib)
    cycle = [fib[0]]
    visited = {fib[0]}
    while len(cycle) < 10:
        curr = cycle[-1]
        for j in fib_set - visited:
            if adj[curr, j] > 0.5:
                cycle.append(j)
                visited.add(j)
                break

    if len(cycle) == 10:
        for k in range(10):
            i = cycle[k]
            j = cycle[(k + 1) % 10]
            J[i, j] = 1j   # forward
            J[j, i] = -1j  # backward

# B1 = b1 * (i/10) * J... but actually d(Box)/dPhi = b1 * d(A_fiber)/dPhi
# d(A_fiber(Phi))_ij/dPhi = (i/10)*exp(i*Phi/10) for forward edge
#                         = (-i/10)*exp(-i*Phi/10) for backward edge
# At Phi=0: d(A_fiber)/dPhi = (1/10)*J (where J has i and -i)
B1 = b1 * J / 10  # This is the first-order perturbation (anti-Hermitian)

# B2 = b1 * d^2(A_fiber)/dPhi^2|_0
# d^2(A_f)/dPhi^2|_0 = (i/10)^2 * (exp(i*0) on forward + exp(-i*0) on backward)
#                     = -1/100 * A_fiber
B2 = -b1 * A_fiber / 100

print(f"  First-order perturbation B1 = b1*J/10")
print(f"  ||B1||_F = {np.linalg.norm(B1):.6f}")
print(f"  B1 is anti-Hermitian: {np.allclose(B1, -B1.conj().T)}")

print(f"\n  Second-order perturbation B2 = -b1*A_fiber/100")
print(f"  ||B2||_F = {np.linalg.norm(B2):.6f}")

# Perturbative expansion of Tr(Box^2):
# Tr(Box(Phi)^2) = Tr(Box^2) + 2*Phi*Tr(Box*B1) + Phi^2*(Tr(B1^2) + Tr(Box*B2)) + ...
# Since B1 is anti-Hermitian and Box is real symmetric:
# Tr(Box*B1) = sum_ij Box_ij * B1_ji (complex in general)

Box_real_mat = b1 * A_fiber - adj
tr_box_B1 = np.trace(Box_real_mat @ B1)
tr_B1_sq = np.trace(B1 @ B1)
tr_box_B2 = np.trace(Box_real_mat @ B2)

print(f"\n  Tr(Box*B1) = {tr_box_B1:.6f} (expect 0 by symmetry)")
print(f"  Tr(B1^2) = {tr_B1_sq:.6f}")
print(f"  Tr(Box*B2) = {tr_box_B2:.6f}")
print(f"  Tr(B1^2) + Tr(Box*B2) = {tr_B1_sq + tr_box_B2:.6f}")

# So Tr(Box(Phi)^2) = 7200 + Phi^2*(Tr(B1^2)+Tr(Box*B2)) + O(Phi^4)
coeff_phi2 = (tr_B1_sq + tr_box_B2).real
print(f"\n  Tr(Box(Phi)^2) = 7200 + {coeff_phi2:.4f}*Phi^2 + O(Phi^4)")
print(f"  This is a quadratic in Phi^2: {'minimum at 0' if coeff_phi2 > 0 else 'maximum at 0'}")

# For the spectral determinant perturbation:
# d/dPhi log|det'(Box)| = Tr(Box^{-1} * dBox/dPhi) restricted to nonzero eigenspace
# = Tr'(Box^{-1} * B1) * Phi (to first order)

# Compute Box^{-1} on nonzero eigenspace
evals_box0, evecs_box0 = np.linalg.eigh(Box_real_mat)
nonzero_mask = np.abs(evals_box0) > 1e-6
evals_inv = np.zeros_like(evals_box0)
evals_inv[nonzero_mask] = 1.0 / evals_box0[nonzero_mask]
Box_inv = evecs_box0 @ np.diag(evals_inv) @ evecs_box0.T

tr_Boxinv_B1 = np.trace(Box_inv @ B1)
print(f"\n  Tr'(Box^{{-1}} * B1) = {tr_Boxinv_B1:.10f}")
print(f"  (Should be 0 by time-reversal symmetry of Box)")


# =====================================================================
# PART 5: SPECTRAL ASYMMETRY AND GAUGE COUPLING
# =====================================================================
print("\n" + "=" * 70)
print("PART 5: ETA INVARIANT AND GAUGE COUPLING")
print("=" * 70)

# The eta invariant eta(s) = sum sign(mu_i)|mu_i|^{-s} is related to
# the spectral asymmetry. For the Dirac operator on a circle with flux Phi:
# eta(0, Phi) = 1 - 2*{Phi/(2*pi)} (fractional part)
# This is the Atiyah-Patodi-Singer eta invariant.

# For our Box operator, eta(0) at flux Phi:
eta_vals = np.zeros(n_flux)
for i, phi in enumerate(flux_vals):
    Box_phi = build_box_with_flux(fibers, adj, phi)
    evals = np.sort(np.linalg.eigvalsh(Box_phi))
    nonzero = evals[np.abs(evals) > 1e-4]
    eta_vals[i] = np.sum(np.sign(nonzero))

print(f"  eta(0) at Phi=0: {eta_vals[0]:.0f}")
print(f"  eta(0) range: [{np.min(eta_vals):.0f}, {np.max(eta_vals):.0f}]")

# Find where eta jumps (spectral flow = eigenvalue crossing zero)
jumps = []
for i in range(len(eta_vals)-1):
    if abs(eta_vals[i+1] - eta_vals[i]) > 0.5:
        jumps.append((flux_vals[i], flux_vals[i+1], eta_vals[i], eta_vals[i+1]))

print(f"\n  Spectral flow (eta jumps):")
for phi1, phi2, e1, e2 in jumps:
    phi_mid = (phi1+phi2)/2
    print(f"    Phi ~ {phi_mid:.4f} ({phi_mid/np.pi:.4f}*pi): "
          f"eta {e1:.0f} -> {e2:.0f}")

# Total spectral flow from 0 to 2*pi
total_flow = int(eta_vals[-1] - eta_vals[0])
print(f"\n  Total spectral flow [0, 2*pi]: {total_flow}")
print(f"  Expected from index theorem: related to Chern number c1 of Hopf bundle")


# =====================================================================
# PART 6: KEY RATIO SEARCH
# =====================================================================
print("\n" + "=" * 70)
print("PART 6: SEARCHING FOR alpha IN SPECTRAL ACTION")
print("=" * 70)

# Compute Tr(Box(Phi)^2) accurately vs Phi
S2_accurate = np.zeros(n_flux)
for i, phi in enumerate(flux_vals):
    Box_phi = build_box_with_flux(fibers, adj, phi)
    evals = np.linalg.eigvalsh(Box_phi)
    S2_accurate[i] = np.sum(evals**2)

# Tr(Box(Phi)^2) should be a smooth function of Phi
# Find its structure
print(f"  Tr(Box(Phi)^2):")
print(f"    Phi=0:    {S2_accurate[0]:.4f}")
print(f"    Phi=pi/5: {S2_accurate[n_flux//10]:.4f}")
print(f"    Phi=pi:   {S2_accurate[n_flux//2]:.4f}")
print(f"    Phi=2pi:  {S2_accurate[-1]:.4f}")

# The ratio S2(Phi)/S2(0) at special values
print(f"\n  Ratios S2(Phi)/S2(0):")
for label, idx in [("0", 0), ("pi/5", n_flux//10), ("pi/2", n_flux//4),
                    ("pi", n_flux//2), ("3pi/2", 3*n_flux//4), ("2pi", -1)]:
    ratio = S2_accurate[idx] / S2_accurate[0]
    print(f"    Phi={label:>5s}: S2/S2(0) = {ratio:.10f}")

# Check: is Tr(Box(Phi)^2) = 7200 - K*cos(Phi) for some K?
# By Fourier analysis of S2 vs Phi
from numpy.fft import fft
S2_fft = fft(S2_accurate - np.mean(S2_accurate))
print(f"\n  Fourier analysis of Tr(Box(Phi)^2):")
for k in range(6):
    amp = np.abs(S2_fft[k]) / n_flux
    phase = np.angle(S2_fft[k])
    if amp > 0.01:
        print(f"    k={k}: amplitude = {amp:.6f}, phase = {phase:.4f}")

# Most important: check if there's a connection to alpha at any extremum
# The spectral action S = Tr f(Box/Lambda) extremized over Lambda gives
# conditions on the spectrum. But we can also look for:
# S(Phi) has a SADDLE POINT at Phi_* where S(Phi_*) encodes alpha

# Final check: compute det'(Box(Phi)) / det'(Box(0)) at special Phi values
print(f"\n  log|det'(Box(Phi))| - log|det'(Box(0))| at special Phi:")
S_det_0 = S_det[0]
for label, phi_val in [("pi/a1=pi/5", np.pi/a1),
                        ("2*pi*alpha", 2*np.pi*alpha_phys),
                        ("alpha", alpha_phys),
                        ("pi", np.pi),
                        ("2*pi", 2*np.pi)]:
    Box_phi = build_box_with_flux(fibers, adj, phi_val)
    evals = np.linalg.eigvalsh(Box_phi)
    nonzero = evals[np.abs(evals) > 1e-6]
    det_val = np.sum(np.log(np.abs(nonzero)))
    delta = det_val - S_det_0
    print(f"    Phi={label:>15s} ({phi_val:.6f}): delta = {delta:.6f}")

print("\n" + "=" * 70)
print("EXP-560 COMPLETE")
print("=" * 70)
