"""
exp561: Fine-scan -- is alpha the critical flux where Box loses a kernel mode?

From exp560: the first spectral flow (eta jump) occurs at Phi ~ 0.047,
suspiciously close to 2*pi*alpha ~ 0.0459.

QUESTION: Is the first zero-crossing of a Box eigenvalue at EXACTLY Phi = 2*pi*alpha?

BLIND: find the crossing point first, then compare with 2*pi*alpha.
"""

import numpy as np
from numpy.linalg import eigvalsh
import sys
sys.path.insert(0, ".")
from commons import build_600cell

PHI_GR = (1 + np.sqrt(5)) / 2
a1 = 5
b1 = 6
N = 120

verts, adj, lap = build_600cell()

# === Build Hopf fibration with ordered cycles ===
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

# Find ordered fibration
def find_ordered_fibration():
    target_w = PHI_GR / 2.0
    for i in range(N):
        if abs(verts[i, 0] - target_w) < 1e-6:
            g = verts[i]
            p = g.copy()
            ok = True
            for k in range(2, 11):
                p = qmul(p, g)
                if k == 5 and not np.allclose(p, [-1,0,0,0], atol=1e-6):
                    ok = False; break
                if k == 10 and not np.allclose(p, [1,0,0,0], atol=1e-6):
                    ok = False
            if not ok:
                continue
            used = set()
            fibers = []
            for s in range(N):
                if s in used: continue
                fib = []
                q = verts[s].copy()
                for k in range(10):
                    idx = find_idx(q, verts)
                    if idx < 0 or idx in used: break
                    fib.append(idx)
                    used.add(idx)
                    q = qmul(q, g)
                if len(fib) == 10:
                    fibers.append(fib)
            if len(fibers) == 12 and len(used) == N:
                return fibers
    return None

fibers = find_ordered_fibration()
assert fibers is not None

# Build ordered cycles from fibers using adjacency
def get_ordered_cycles(fibers, adj):
    """Convert fiber vertex lists to properly ordered cycles using adjacency."""
    cycles = []
    for fib in fibers:
        fib_set = set(fib)
        cycle = [fib[0]]
        visited = {fib[0]}
        while len(cycle) < 10:
            curr = cycle[-1]
            found = False
            for j in fib_set - visited:
                if adj[curr, j] > 0.5:
                    cycle.append(j)
                    visited.add(j)
                    found = True
                    break
            if not found:
                break
        if len(cycle) == 10:
            cycles.append(cycle)
    return cycles

cycles = get_ordered_cycles(fibers, adj)
assert len(cycles) == 12


def build_box_with_flux(cycles, adj, flux):
    """Build Box(Phi) = b1*A_fiber(Phi) - A with Peierls phases on fiber edges."""
    phase = flux / 10
    A_f = np.zeros((N, N), dtype=complex)
    for cyc in cycles:
        for k in range(10):
            i = cyc[k]
            j = cyc[(k + 1) % 10]
            A_f[i, j] = np.exp(1j * phase)
            A_f[j, i] = np.exp(-1j * phase)
    return b1 * A_f - adj.astype(complex)


# =====================================================================
# STEP 1: Identify which eigenvalue crosses zero
# =====================================================================
print("=" * 70)
print("STEP 1: IDENTIFY THE CROSSING EIGENVALUE")
print("=" * 70)

# At Phi=0, the 9 kernel eigenvalues are exactly 0.
# Track the 9 near-zero eigenvalues as Phi increases slightly.
Box0 = build_box_with_flux(cycles, adj, 0.0)
evals0 = np.sort(eigvalsh(Box0))
kernel_indices = np.where(np.abs(evals0) < 1e-6)[0]
print(f"\n  Kernel at Phi=0: indices {kernel_indices}, dim = {len(kernel_indices)}")

# Track eigenvalues near zero for small Phi
print(f"\n  Eigenvalue tracking near zero (9 kernel modes):")
print(f"  {'Phi':>10s}", end="")
for k in range(9):
    print(f"  {'mu_'+str(k):>10s}", end="")
print()

for phi_test in [0, 0.001, 0.005, 0.01, 0.02, 0.03, 0.04, 0.045, 0.05, 0.06]:
    Box_phi = build_box_with_flux(cycles, adj, phi_test)
    evals = np.sort(eigvalsh(Box_phi))
    # Find the 9 eigenvalues closest to zero
    near_zero = evals[kernel_indices[0]:kernel_indices[-1]+1]
    print(f"  {phi_test:10.4f}", end="")
    for val in near_zero:
        print(f"  {val:10.6f}", end="")
    print()


# =====================================================================
# STEP 2: Fine scan for the FIRST zero crossing
# =====================================================================
print("\n" + "=" * 70)
print("STEP 2: FINE SCAN [0.04, 0.05] -- 10000 points")
print("=" * 70)

n_scan = 10000
phi_lo, phi_hi = 0.04, 0.05
flux_scan = np.linspace(phi_lo, phi_hi, n_scan)

# Track the eigenvalue that's about to cross zero
# At Phi=0.04, which eigenvalue is closest to zero?
Box_test = build_box_with_flux(cycles, adj, 0.04)
evals_test = np.sort(eigvalsh(Box_test))
# Find the most negative eigenvalue among the "kernel" group
near_zero_test = evals_test[kernel_indices[0]:kernel_indices[-1]+1]
crossing_idx_in_kernel = np.argmin(near_zero_test)
crossing_global_idx = kernel_indices[0] + crossing_idx_in_kernel
print(f"  Eigenvalue closest to crossing: index {crossing_global_idx}")
print(f"  Value at Phi=0.04: {evals_test[crossing_global_idx]:.10f}")

# Track this eigenvalue through the scan
crossing_evals = np.zeros(n_scan)
for i, phi in enumerate(flux_scan):
    Box_phi = build_box_with_flux(cycles, adj, phi)
    evals = np.sort(eigvalsh(Box_phi))
    crossing_evals[i] = evals[crossing_global_idx]

# Find zero crossing by linear interpolation
crossing_phi = None
for i in range(n_scan - 1):
    if crossing_evals[i] * crossing_evals[i+1] < 0:
        # Linear interpolation
        phi_cross = flux_scan[i] - crossing_evals[i] * (flux_scan[i+1] - flux_scan[i]) / (crossing_evals[i+1] - crossing_evals[i])
        crossing_phi = phi_cross
        print(f"\n  ZERO CROSSING FOUND!")
        print(f"    Phi_cross = {phi_cross:.12f}")
        print(f"    mu({flux_scan[i]:.6f}) = {crossing_evals[i]:.10f}")
        print(f"    mu({flux_scan[i+1]:.6f}) = {crossing_evals[i+1]:.10f}")
        break

if crossing_phi is None:
    # Check if eigenvalue is still positive or already negative throughout
    print(f"\n  No crossing in [{phi_lo}, {phi_hi}]")
    print(f"  Eigenvalue range: [{np.min(crossing_evals):.8f}, {np.max(crossing_evals):.8f}]")
    # Extend search
    print(f"\n  Extending search to [0.0, 0.1]...")
    n_ext = 10000
    flux_ext = np.linspace(0.0, 0.1, n_ext)
    cross_ext = np.zeros(n_ext)
    for i, phi in enumerate(flux_ext):
        Box_phi = build_box_with_flux(cycles, adj, phi)
        evals = np.sort(eigvalsh(Box_phi))
        cross_ext[i] = evals[crossing_global_idx]

    for i in range(n_ext - 1):
        if cross_ext[i] * cross_ext[i+1] < 0:
            phi_cross = flux_ext[i] - cross_ext[i] * (flux_ext[i+1] - flux_ext[i]) / (cross_ext[i+1] - cross_ext[i])
            crossing_phi = phi_cross
            print(f"  CROSSING at Phi = {phi_cross:.12f}")
            break

    if crossing_phi is None:
        # Track all 9 kernel eigenvalues to find which one crosses first
        print(f"\n  Tracking ALL 9 kernel eigenvalues...")
        for ki in range(9):
            gi = kernel_indices[ki]
            vals = np.zeros(1000)
            fluxes = np.linspace(0, 0.1, 1000)
            for j, phi in enumerate(fluxes):
                Box_phi = build_box_with_flux(cycles, adj, phi)
                evals = np.sort(eigvalsh(Box_phi))
                vals[j] = evals[gi]
            # Find crossing
            for j in range(999):
                if vals[j] * vals[j+1] < 0:
                    pc = fluxes[j] - vals[j]*(fluxes[j+1]-fluxes[j])/(vals[j+1]-vals[j])
                    print(f"    Kernel eigenvalue {ki} (idx {gi}): crosses at Phi = {pc:.10f}")
                    if crossing_phi is None:
                        crossing_phi = pc
                    break
            else:
                print(f"    Kernel eigenvalue {ki} (idx {gi}): no crossing in [0, 0.1], "
                      f"range [{np.min(vals):.6f}, {np.max(vals):.6f}]")


# =====================================================================
# STEP 3: COMPARE WITH 2*pi*alpha (BLIND)
# =====================================================================
print("\n" + "=" * 70)
print("STEP 3: BLIND COMPARISON")
print("=" * 70)

# Alpha equation
B = 4 * a1 * PHI_GR**4
alpha_phys = (B - np.sqrt(B**2 - 8*np.pi)) / (4*np.pi)
two_pi_alpha = 2 * np.pi * alpha_phys

print(f"\n  MEASURED crossing point:  Phi_cross = {crossing_phi:.12f}" if crossing_phi else "\n  No crossing found!")

if crossing_phi is not None:
    print(f"  PREDICTED 2*pi*alpha:     Phi_pred  = {two_pi_alpha:.12f}")
    print(f"  Difference:               delta     = {crossing_phi - two_pi_alpha:.6e}")
    print(f"  Relative:                 delta/Phi = {(crossing_phi - two_pi_alpha)/two_pi_alpha:.6e}")
    ratio = crossing_phi / two_pi_alpha
    print(f"  Ratio Phi_cross/(2*pi*alpha) = {ratio:.10f}")

    if abs(ratio - 1) < 0.001:
        print(f"\n  *** MATCH TO 0.1% -- alpha IS the critical flux! ***")
    elif abs(ratio - 1) < 0.01:
        print(f"\n  ** Close match to 1% -- suggestive but not exact **")
    elif abs(ratio - 1) < 0.05:
        print(f"\n  * Within 5% -- near-miss, investigate further *")
    else:
        print(f"\n  NEGATIVE: crossing at {crossing_phi:.6f} != 2*pi*alpha = {two_pi_alpha:.6f}")
        print(f"  Ratio = {ratio:.6f}, discrepancy = {abs(ratio-1)*100:.1f}%")

    # Also check other interpretations
    print(f"\n  Other possible identifications:")
    for label, val in [
        ("alpha", alpha_phys),
        ("2*pi*alpha", two_pi_alpha),
        ("alpha*pi/a1", alpha_phys*np.pi/a1),
        ("pi/(a1*1/alpha)", np.pi*alpha_phys/a1),
        ("1/(2*a1*1/alpha)", 1/(2*a1/alpha_phys)),
        ("sqrt(alpha)*pi/a1", np.sqrt(alpha_phys)*np.pi/a1),
    ]:
        ratio_check = crossing_phi / val
        if 0.5 < ratio_check < 2.0:
            match = "***" if abs(ratio_check - 1) < 0.01 else ""
            print(f"    Phi_cross / ({label}) = {ratio_check:.10f} {match}")

    # Check if crossing is at a simple fraction of pi/a1
    ratio_pi5 = crossing_phi / (np.pi / a1)
    print(f"\n  Phi_cross / (pi/a1) = {ratio_pi5:.10f}")
    print(f"  Phi_cross / (pi/a1^2) = {crossing_phi / (np.pi / a1**2):.10f}")

    # Check algebraic expressions
    print(f"\n  Algebraic checks:")
    print(f"    Phi_cross * a1 / pi = {crossing_phi * a1 / np.pi:.10f}")
    print(f"    Phi_cross * a1^2 / pi = {crossing_phi * a1**2 / np.pi:.10f}")
    print(f"    Phi_cross * 137 = {crossing_phi * 137:.10f}")
    print(f"    1/Phi_cross = {1/crossing_phi:.6f}")


# =====================================================================
# STEP 4: Refine with bisection (if crossing found)
# =====================================================================
if crossing_phi is not None:
    print("\n" + "=" * 70)
    print("STEP 4: BISECTION REFINEMENT (20 iterations)")
    print("=" * 70)

    lo, hi = crossing_phi - 0.001, crossing_phi + 0.001
    # Verify bracket
    Box_lo = build_box_with_flux(cycles, adj, lo)
    Box_hi = build_box_with_flux(cycles, adj, hi)
    val_lo = np.sort(eigvalsh(Box_lo))[crossing_global_idx]
    val_hi = np.sort(eigvalsh(Box_hi))[crossing_global_idx]

    if val_lo * val_hi > 0:
        # Adjust bracket
        lo, hi = 0.0, 0.1
        val_lo = 0  # we know eigenvalue is 0 at Phi=0
        Box_hi2 = build_box_with_flux(cycles, adj, hi)
        val_hi = np.sort(eigvalsh(Box_hi2))[crossing_global_idx]

    for iteration in range(30):
        mid = (lo + hi) / 2
        Box_mid = build_box_with_flux(cycles, adj, mid)
        val_mid = np.sort(eigvalsh(Box_mid))[crossing_global_idx]

        if abs(val_mid) < 1e-14:
            break
        if val_lo * val_mid < 0:
            hi = mid
            val_hi = val_mid
        else:
            lo = mid
            val_lo = val_mid

    crossing_refined = (lo + hi) / 2
    print(f"  Refined crossing: Phi = {crossing_refined:.15f}")
    print(f"  2*pi*alpha:             {two_pi_alpha:.15f}")
    print(f"  Difference:             {crossing_refined - two_pi_alpha:.6e}")
    print(f"  Ratio:                  {crossing_refined/two_pi_alpha:.12f}")

    match_pct = abs(crossing_refined/two_pi_alpha - 1) * 100
    print(f"\n  Discrepancy: {match_pct:.4f}%")

    if match_pct < 0.01:
        print(f"  *** EXACT MATCH (< 0.01%) ***")
    elif match_pct < 0.1:
        print(f"  ** Very close (< 0.1%) **")
    elif match_pct < 1:
        print(f"  * Close (< 1%) *")
    else:
        print(f"  NEGATIVE: {match_pct:.2f}% discrepancy")

print("\n" + "=" * 70)
print("EXP-561 COMPLETE")
print("=" * 70)
