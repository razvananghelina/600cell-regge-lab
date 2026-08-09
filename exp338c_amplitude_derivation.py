"""
exp338c: Derive mode amplitudes from physics
=============================================
From exp338b, Model F gives RMS=0.033 with 4 mode amplitudes.
The assignment: e->0, mu->3, tau->5, u->4, c->2, t->8, d->7, s->1, b->6

Now: can we derive the 4 amplitudes from a 1-2 parameter function?
If c_m = f(lambda_m, d_m) with 1-2 params, this becomes a genuine prediction.

Also try: one-loop with gauge + fermion contributions (competing signs).
"""
import numpy as np
from itertools import product, permutations

PHI = (1 + np.sqrt(5)) / 2
ALPHA_S = 1 / (2 * PHI**3)

# ================================================================
# Setup from exp338b
# ================================================================
dims = [1, 2, 3, 4, 5, 6, 4, 2, 3]

CG_self = np.array([
    [1, 0, 0, 0, 0, 0, 0, 0, 0],
    [1, 0, 1, 0, 0, 0, 0, 0, 0],
    [1, 0, 1, 0, 1, 0, 0, 0, 0],
    [1, 0, 1, 0, 1, 0, 1, 0, 1],
    [1, 0, 1, 0, 2, 0, 2, 0, 1],
    [1, 0, 2, 0, 3, 0, 2, 0, 2],
    [1, 0, 1, 0, 1, 0, 1, 0, 1],
    [1, 0, 0, 0, 0, 0, 0, 0, 1],
    [1, 0, 0, 0, 1, 0, 0, 0, 1],
])

char_C1 = [1, PHI, PHI, 1, 0, -1, -1, -1/PHI, -1/PHI]
cayley_eigs = np.array([12 * (1 - char_C1[k] / dims[k]) for k in range(9)])

even_modes = [2, 4, 6, 8]
lam = np.array([cayley_eigs[m] for m in even_modes])  # [5.528, 12.0, 15.0, 14.472]
d_m = np.array([dims[m] for m in even_modes])          # [3, 5, 4, 3]

fermion_names = ['e', 'mu', 'tau', 'u', 'c', 't', 'd', 's', 'b']
fermion_T3 = np.array([-0.5, -0.5, -0.5, 0.5, 0.5, 0.5, -0.5, -0.5, -0.5])
fermion_Nc = np.array([0, 0, 0, 1, 1, 1, 1, 1, 1])
delta_exp = np.array([0.000, 0.080, -0.055, -0.004, 0.247, 0.456, -0.402, -0.177, -0.278])

# Best T3-signed assignment from Model F
assign_F = [0, 3, 5, 4, 2, 8, 7, 1, 6]
# e->0, mu->3, tau->5, u->4, c->2, t->8, d->7, s->1, b->6

print("="*72)
print("exp338c: Derive Mode Amplitudes from Physics")
print("="*72)

# Best-fit amplitudes from Model F
target_c = np.array([0.1189, 0.0259, -0.1774, 0.2635])
print(f"\nTarget amplitudes (from Model F fit):")
print(f"  c_2 = {target_c[0]:.4f} (lambda={lam[0]:.3f}, d=3)")
print(f"  c_4 = {target_c[1]:.4f} (lambda={lam[1]:.3f}, d=5)")
print(f"  c_6 = {target_c[2]:.4f} (lambda={lam[2]:.3f}, d=4)")
print(f"  c_8 = {target_c[3]:.4f} (lambda={lam[3]:.3f}, d=3)")

# Express in framework units
print(f"\nIn units of alpha_s = {ALPHA_S:.6f}:")
for i, m in enumerate(even_modes):
    print(f"  c_{m}/alpha_s = {target_c[i]/ALPHA_S:.4f}")


# ================================================================
# Build the prediction function for assignment F
# ================================================================
def predict_delta(c_modes, assignment=assign_F):
    """Given 4 mode amplitudes, predict deltas using T3*Nc*CG structure."""
    pred = np.zeros(9)
    for f in range(9):
        k = assignment[f]
        for idx, m in enumerate(even_modes):
            pred[f] += fermion_T3[f] * fermion_Nc[f] * CG_self[k, m] * dims[m] * c_modes[idx]
    return pred


# ================================================================
# SCAN 1: Single-parameter functional forms c_m = f(lambda_m)
# ================================================================
print("\n" + "="*72)
print("SCAN 1: One-parameter functions c_m = f(lambda_m)")
print("="*72)

results = []

# 1a: c_m = A/lambda_m
for A in np.linspace(-5, 5, 10000):
    c = A / lam
    pred = predict_delta(c)
    rms = np.sqrt(np.mean((delta_exp - pred)**2))
    results.append(('A/lam', 1, A, rms, c.copy()))

# 1b: c_m = A * lambda_m
for A in np.linspace(-0.5, 0.5, 10000):
    c = A * lam
    pred = predict_delta(c)
    rms = np.sqrt(np.mean((delta_exp - pred)**2))
    results.append(('A*lam', 1, A, rms, c.copy()))

# 1c: c_m = A / d_m
for A in np.linspace(-3, 3, 10000):
    c = A / d_m
    pred = predict_delta(c)
    rms = np.sqrt(np.mean((delta_exp - pred)**2))
    results.append(('A/d', 1, A, rms, c.copy()))

# 1d: c_m = A * d_m
for A in np.linspace(-0.5, 0.5, 10000):
    c = A * d_m
    pred = predict_delta(c)
    rms = np.sqrt(np.mean((delta_exp - pred)**2))
    results.append(('A*d', 1, A, rms, c.copy()))

# 1e: c_m = A * (lambda_m - 12) (centered at rho_4)
for A in np.linspace(-0.5, 0.5, 10000):
    c = A * (lam - 12)
    pred = predict_delta(c)
    rms = np.sqrt(np.mean((delta_exp - pred)**2))
    results.append(('A*(lam-12)', 1, A, rms, c.copy()))

# 1f: c_m = A * chi_m(C1) / d_m (character-weighted)
chi_even = np.array([char_C1[m] for m in even_modes])
for A in np.linspace(-2, 2, 10000):
    c = A * chi_even / d_m
    pred = predict_delta(c)
    rms = np.sqrt(np.mean((delta_exp - pred)**2))
    results.append(('A*chi/d', 1, A, rms, c.copy()))

# Sort by RMS
results.sort(key=lambda x: x[3])
print(f"\nTop 10 one-parameter fits:")
print(f"{'Form':>15} {'Params':>7} {'Param':>8} {'RMS':>8}")
print("-"*45)
seen = set()
for form, npar, param, rms, c in results[:30]:
    if form not in seen:
        seen.add(form)
        print(f"{form:>15} {npar:>7} {param:>8.4f} {rms:>8.4f}")
    if len(seen) >= 10:
        break


# ================================================================
# SCAN 2: Two-parameter functional forms
# ================================================================
print("\n" + "="*72)
print("SCAN 2: Two-parameter functions")
print("="*72)

best_2p = {}

# 2a: c_m = A / (r * lambda_m + 1) (propagator)
best = (999, 0, 0)
for A in np.linspace(-3, 3, 300):
    for r in np.logspace(-3, 2, 300):
        c = A / (r * lam + 1)
        pred = predict_delta(c)
        rms = np.sqrt(np.mean((delta_exp - pred)**2))
        if rms < best[0]:
            best = (rms, A, r)
best_2p['A/(r*lam+1)'] = best

# 2b: c_m = A + B * lambda_m (affine)
best = (999, 0, 0)
for A in np.linspace(-1, 1, 500):
    for B in np.linspace(-0.1, 0.1, 500):
        c = A + B * lam
        pred = predict_delta(c)
        rms = np.sqrt(np.mean((delta_exp - pred)**2))
        if rms < best[0]:
            best = (rms, A, B)
best_2p['A+B*lam'] = best

# 2c: c_m = A/d_m + B*d_m
best = (999, 0, 0)
for A in np.linspace(-2, 2, 500):
    for B in np.linspace(-0.5, 0.5, 500):
        c = A / d_m + B * d_m
        pred = predict_delta(c)
        rms = np.sqrt(np.mean((delta_exp - pred)**2))
        if rms < best[0]:
            best = (rms, A, B)
best_2p['A/d+B*d'] = best

# 2d: c_m = A * (lambda_m - B) (shifted)
best = (999, 0, 0)
for A in np.linspace(-0.3, 0.3, 500):
    for B in np.linspace(0, 20, 500):
        c = A * (lam - B)
        pred = predict_delta(c)
        rms = np.sqrt(np.mean((delta_exp - pred)**2))
        if rms < best[0]:
            best = (rms, A, B)
best_2p['A*(lam-B)'] = best

# 2e: c_m = A * chi_m/d_m + B
best = (999, 0, 0)
for A in np.linspace(-3, 3, 500):
    for B in np.linspace(-1, 1, 500):
        c = A * chi_even / d_m + B
        pred = predict_delta(c)
        rms = np.sqrt(np.mean((delta_exp - pred)**2))
        if rms < best[0]:
            best = (rms, A, B)
best_2p['A*chi/d+B'] = best

# 2f: c_m = A / (lambda_m + B)
best = (999, 0, 0)
for A in np.linspace(-50, 50, 500):
    for B in np.linspace(-5, 20, 500):
        if all(lam[i] + B > 0.01 for i in range(4)):
            c = A / (lam + B)
            pred = predict_delta(c)
            rms = np.sqrt(np.mean((delta_exp - pred)**2))
            if rms < best[0]:
                best = (rms, A, B)
best_2p['A/(lam+B)'] = best

# 2g: c_m = A * lambda_m^p (power law)
best = (999, 0, 0)
for A in np.linspace(-2, 2, 500):
    for p in np.linspace(-3, 3, 500):
        c = A * lam**p
        pred = predict_delta(c)
        rms = np.sqrt(np.mean((delta_exp - pred)**2))
        if rms < best[0]:
            best = (rms, A, p)
best_2p['A*lam^p'] = best

# 2h: c_m = A/lam + B/d_m
best = (999, 0, 0)
for A in np.linspace(-5, 5, 500):
    for B in np.linspace(-3, 3, 500):
        c = A / lam + B / d_m
        pred = predict_delta(c)
        rms = np.sqrt(np.mean((delta_exp - pred)**2))
        if rms < best[0]:
            best = (rms, A, B)
best_2p['A/lam+B/d'] = best

# 2i: c_m = (A + B*d_m) / lambda_m
best = (999, 0, 0)
for A in np.linspace(-5, 5, 500):
    for B in np.linspace(-3, 3, 500):
        c = (A + B * d_m) / lam
        pred = predict_delta(c)
        rms = np.sqrt(np.mean((delta_exp - pred)**2))
        if rms < best[0]:
            best = (rms, A, B)
best_2p['(A+B*d)/lam'] = best

# 2j: c_m = A*(-1)^(m/2) + B (alternating sign!)
# m = 2,4,6,8 -> m/2 = 1,2,3,4 -> (-1)^(m/2) = -1,1,-1,1
alt_sign = np.array([-1, 1, -1, 1])
best = (999, 0, 0)
for A in np.linspace(-0.5, 0.5, 500):
    for B in np.linspace(-0.5, 0.5, 500):
        c = A * alt_sign + B
        pred = predict_delta(c)
        rms = np.sqrt(np.mean((delta_exp - pred)**2))
        if rms < best[0]:
            best = (rms, A, B)
best_2p['A*(-1)^k+B'] = best

# Sort 2-param results
sorted_2p = sorted(best_2p.items(), key=lambda x: x[1][0])
print(f"\n{'Form':>20} {'RMS':>8} {'param1':>10} {'param2':>10}")
print("-"*55)
for form, (rms, p1, p2) in sorted_2p:
    print(f"{form:>20} {rms:>8.4f} {p1:>10.4f} {p2:>10.4f}")


# ================================================================
# SCAN 3: Galois-motivated amplitudes
# ================================================================
print("\n" + "="*72)
print("SCAN 3: Framework-motivated functional forms")
print("="*72)

# The Cayley eigenvalues have algebraic forms:
# lam_2 = 12 - 4*phi = 12 - 4*1.618 = 5.528
# lam_4 = 12
# lam_6 = 15
# lam_8 = 12 + 4/phi = 12 + 2.472 = 14.472
#
# Note: lam_2 and lam_8 are GALOIS CONJUGATES!
# sigma(lam_2) = sigma(12-4*phi) = 12-4*phi' = 12+4/phi = lam_8
#
# And lam_4 = 12 (rational), lam_6 = 15 (rational)

# Galois conjugation maps: lam_2 <-> lam_8, lam_4 <-> lam_4, lam_6 <-> lam_6
# If the amplitude has a Galois structure:
# c_2 and c_8 should be Galois conjugates!
# And c_4, c_6 should be rational.

print("\nGalois structure of Cayley eigenvalues:")
print(f"  lam_2 = 12 - 4*phi = {12 - 4*PHI:.4f}")
print(f"  lam_4 = 12")
print(f"  lam_6 = 15")
print(f"  lam_8 = 12 + 4/phi = {12 + 4/PHI:.4f}")
print(f"  sigma(lam_2) = lam_8 (Galois conjugate pair)")

# If c_m = f(lam_m) and f is defined over Q(phi):
# Then c_2 = a + b*phi, c_8 = a + b*phi' = a - b/phi (Galois conjugate)
# And c_4, c_6 are rational.
# This gives 4 params reduced to: a, b (for the pair) + c_4, c_6 = 4 params.
# No reduction yet.

# But if f is a RATIONAL function of lambda:
# c_m = (alpha + beta*lam) / (gamma + delta*lam)
# Then c_2 and c_8 are automatically Galois conjugates!
# And we have 3 effective params (ratio removes one).

# Let's check: do the target amplitudes satisfy c_2 <-> c_8 Galois?
# c_2 = 0.1189, c_8 = 0.2635
# If c = a + b*phi: 0.1189 = a + 1.618*b, 0.2635 = a - 0.618*b
# Subtracting: 0.1189 - 0.2635 = (1.618+0.618)*b = 2.236*b
# b = -0.1446/2.236 = -0.0647
# a = 0.1189 - 1.618*(-0.0647) = 0.1189 + 0.1046 = 0.2235
# Check: a - 0.618*b = 0.2235 + 0.0400 = 0.2635 ✓

a_pair = 0.2235
b_pair = -0.0647
print(f"\nGalois decomposition of target amplitudes:")
print(f"  c_2 = {a_pair:.4f} + ({b_pair:.4f})*phi = {a_pair + b_pair*PHI:.4f} (target: 0.1189)")
print(f"  c_8 = {a_pair:.4f} + ({b_pair:.4f})*phi' = {a_pair + b_pair*(1-PHI):.4f} (target: 0.2635)")
print(f"  (c_2 + c_8)/2 = {(0.1189+0.2635)/2:.4f} (= a = Galois-invariant part)")
print(f"  (c_2 - c_8)/(phi-phi') = {(0.1189-0.2635)/(PHI-(1-PHI)):.4f} (= b)")
print(f"  c_4 = {target_c[1]:.4f} (rational)")
print(f"  c_6 = {target_c[2]:.4f} (rational)")

# Now try: c_m = (A + B*lam_m) / (C*lam_m + 1) with 3 params
# This automatically makes (c_2, c_8) a Galois pair
print("\n--- Rational function c_m = (A + B*lam) / (C*lam + 1) ---")
best_3p = (999, 0, 0, 0)
for A in np.linspace(-2, 2, 100):
    for B in np.linspace(-0.3, 0.3, 100):
        for C in np.linspace(-0.1, 0.5, 100):
            denom = C * lam + 1
            if np.all(np.abs(denom) > 0.01):
                c = (A + B * lam) / denom
                pred = predict_delta(c)
                rms = np.sqrt(np.mean((delta_exp - pred)**2))
                if rms < best_3p[0]:
                    best_3p = (rms, A, B, C)

rms, A, B, C = best_3p
c_best3 = (A + B * lam) / (C * lam + 1)
pred_best3 = predict_delta(c_best3)
print(f"  Best: A={A:.4f}, B={B:.4f}, C={C:.4f}, RMS={rms:.4f}")
print(f"  c_m = {np.round(c_best3, 4)}")
for f in range(9):
    print(f"  {fermion_names[f]:>4}: pred={pred_best3[f]:>+.4f}, exp={delta_exp[f]:>+.4f}")

# Try: c_m = A / lam_m + B * lam_m + C
print("\n--- Polynomial c_m = A/lam + B*lam + C ---")
best_3p2 = (999, 0, 0, 0)
for A in np.linspace(-5, 5, 150):
    for B in np.linspace(-0.1, 0.1, 150):
        for C in np.linspace(-2, 2, 150):
            c = A / lam + B * lam + C
            pred = predict_delta(c)
            rms = np.sqrt(np.mean((delta_exp - pred)**2))
            if rms < best_3p2[0]:
                best_3p2 = (rms, A, B, C)

rms, A, B, C = best_3p2
c_best3b = A / lam + B * lam + C
pred_best3b = predict_delta(c_best3b)
print(f"  Best: A={A:.4f}, B={B:.4f}, C={C:.4f}, RMS={rms:.4f}")
print(f"  c_m = {np.round(c_best3b, 4)}")
for f in range(9):
    print(f"  {fermion_names[f]:>4}: pred={pred_best3b[f]:>+.4f}, exp={delta_exp[f]:>+.4f}")


# ================================================================
# SCAN 4: Physical one-loop (gauge + fermion, competing signs)
# ================================================================
print("\n" + "="*72)
print("SCAN 4: Gauge vs Fermion one-loop (competing signs)")
print("="*72)

# c_m = A_gauge * d_m / (r*lam + 1) + A_ferm / (r*lam + s)
# Gauge boson loop: coefficient proportional to d_m (representation dimension)
# Fermion loop: coefficient from CG coupling
# These have OPPOSITE signs (gauge negative, fermion positive)

# Simplified: c_m = A * d_m + B (2 params, no propagator dependence)
best_df = (999, 0, 0)
for A in np.linspace(-0.3, 0.3, 1000):
    for B in np.linspace(-1, 1, 1000):
        c = A * d_m + B
        pred = predict_delta(c)
        rms = np.sqrt(np.mean((delta_exp - pred)**2))
        if rms < best_df[0]:
            best_df = (rms, A, B)

print(f"\nc_m = A*d_m + B: A={best_df[1]:.4f}, B={best_df[2]:.4f}, RMS={best_df[0]:.4f}")
c_df = best_df[1] * d_m + best_df[2]
pred_df = predict_delta(c_df)
print(f"  c_m = {np.round(c_df, 4)}")
for f in range(9):
    print(f"  {fermion_names[f]:>4}: pred={pred_df[f]:>+.4f}, exp={delta_exp[f]:>+.4f}")

# Now: c_m = A * d_m / (r*lam + 1) + B / (r*lam + 1)  [same propagator]
# = (A*d_m + B) / (r*lam + 1)
print("\n--- c_m = (A*d_m + B) / (r*lam + 1) ---")
best_3 = (999, 0, 0, 0)
for A in np.linspace(-0.5, 0.5, 200):
    for B in np.linspace(-2, 2, 200):
        for r in np.logspace(-3, 1, 200):
            c = (A * d_m + B) / (r * lam + 1)
            pred = predict_delta(c)
            rms = np.sqrt(np.mean((delta_exp - pred)**2))
            if rms < best_3[0]:
                best_3 = (rms, A, B, r)

rms, A, B, r = best_3
c_3 = (A * d_m + B) / (r * lam + 1)
pred_3 = predict_delta(c_3)
print(f"  A={A:.4f}, B={B:.4f}, r={r:.4f}, RMS={rms:.4f}")
print(f"  c_m = {np.round(c_3, 4)}")
for f in range(9):
    print(f"  {fermion_names[f]:>4}: pred={pred_3[f]:>+.4f}, exp={delta_exp[f]:>+.4f}")


# ================================================================
# SCAN 5: ALL 40320 assignments with OPTIMAL amplitudes (least squares)
# ================================================================
print("\n" + "="*72)
print("SCAN 5: Optimal amplitudes for ALL 8! assignments")
print("="*72)

# KEY INSIGHT: For a given assignment, pred = M @ c where M is 9x4.
# We can solve for optimal c by least squares! No grid search needed.
# This makes 40320 assignments trivially fast.

remaining_irreps = list(range(1, 9))
dims_even = np.array([dims[m] for m in even_modes])  # [3, 5, 4, 3]

# 2-param linear forms: c = F @ [p1, p2] where F is 4x2
# Form: c = p1*v1 + p2*v2
forms_linear = {
    'A*alt+B': (alt_sign, np.ones(4)),
    'A*lam+B': (lam, np.ones(4)),
    'A/d+B*d': (1.0/d_m, d_m),
    'A*chi+B': (chi_even/d_m, np.ones(4)),
    'A/lam+B/d': (1.0/lam, 1.0/d_m),
    'A*lam+B*d': (lam, d_m.astype(float)),
}

# Also scan with FREE 4 amplitudes (no functional constraint)
print("\n--- FREE optimal amplitudes (4 params, lstsq) ---")
best_free = (999, None, None)
count = 0
for perm in permutations(remaining_irreps):
    assignment = [0] + list(perm)
    # Build design matrix M[f, idx] = T3*Nc * CG_self[k_f, m] * d_m
    M = np.zeros((9, 4))
    for f in range(9):
        k = assignment[f]
        for idx, m_idx in enumerate(even_modes):
            M[f, idx] = fermion_T3[f] * fermion_Nc[f] * CG_self[k, m_idx] * dims[m_idx]
    # Solve M @ c = delta_exp by least squares
    result = np.linalg.lstsq(M, delta_exp, rcond=None)
    c_opt = result[0]
    pred = M @ c_opt
    rms = np.sqrt(np.mean((delta_exp - pred)**2))
    if rms < best_free[0]:
        best_free = (rms, list(assignment), c_opt.copy())
    count += 1

rms, assignment, c_opt = best_free
M = np.zeros((9, 4))
for f in range(9):
    k = assignment[f]
    for idx, m_idx in enumerate(even_modes):
        M[f, idx] = fermion_T3[f] * fermion_Nc[f] * CG_self[k, m_idx] * dims[m_idx]
pred = M @ c_opt
print(f"  Scanned {count} assignments")
print(f"  Best RMS = {rms:.6f}")
print(f"  Assignment: ", end="")
for f in range(9):
    print(f"{fermion_names[f]}->rho_{assignment[f]} ", end="")
print(f"\n  c_m = {np.round(c_opt, 6)}")
for f in range(9):
    print(f"  {fermion_names[f]:>4}: pred={pred[f]:>+.4f}, exp={delta_exp[f]:>+.4f}, err={pred[f]-delta_exp[f]:>+.4f}")

# Show amplitudes in framework units
print(f"\n  Amplitudes in units of alpha_s ({ALPHA_S:.6f}):")
for idx, m in enumerate(even_modes):
    print(f"    c_{m} = {c_opt[idx]:.6f} = {c_opt[idx]/ALPHA_S:.4f} * alpha_s")

# Now scan 2-param linear forms over all assignments
for form_name, (v1, v2) in forms_linear.items():
    print(f"\n--- Form: {form_name} (2 params, lstsq over 8! assignments) ---")
    F = np.column_stack([v1, v2])  # 4x2
    best_form = (999, None, None)

    for perm in permutations(remaining_irreps):
        assignment = [0] + list(perm)
        M = np.zeros((9, 4))
        for f in range(9):
            k = assignment[f]
            for idx, m_idx in enumerate(even_modes):
                M[f, idx] = fermion_T3[f] * fermion_Nc[f] * CG_self[k, m_idx] * dims[m_idx]
        # pred = M @ F @ params = (M@F) @ params
        MF = M @ F  # 9x2
        result = np.linalg.lstsq(MF, delta_exp, rcond=None)
        params = result[0]
        pred = MF @ params
        rms = np.sqrt(np.mean((delta_exp - pred)**2))
        if rms < best_form[0]:
            best_form = (rms, list(assignment), params.copy())

    rms, assignment, params = best_form
    c_form = F @ params
    M = np.zeros((9, 4))
    for f in range(9):
        k = assignment[f]
        for idx, m_idx in enumerate(even_modes):
            M[f, idx] = fermion_T3[f] * fermion_Nc[f] * CG_self[k, m_idx] * dims[m_idx]
    pred = M @ c_form
    print(f"  Best: RMS={rms:.4f}, params=[{params[0]:.4f}, {params[1]:.4f}]")
    print(f"  Assignment: ", end="")
    for f in range(9):
        print(f"{fermion_names[f]}->rho_{assignment[f]} ", end="")
    print(f"\n  c_m = {np.round(c_form, 4)}")
    for f in range(9):
        print(f"  {fermion_names[f]:>4}: pred={pred[f]:>+.4f}, exp={delta_exp[f]:>+.4f}, err={pred[f]-delta_exp[f]:>+.4f}")


# ================================================================
# SCAN 6: Zero-parameter framework formulas over ALL assignments
# ================================================================
print("\n" + "="*72)
print("SCAN 6: Zero-parameter framework formulas + best assignment")
print("="*72)

# Test many zero-param formulas, find best assignment for each
zero_param_formulas = {}

# Cayley propagator: c_m = 1/lam_m (normalized)
zero_param_formulas['1/lam'] = 1.0 / lam
zero_param_formulas['-1/lam'] = -1.0 / lam
# Character-weighted
zero_param_formulas['chi/d'] = chi_even / d_m
zero_param_formulas['-chi/d'] = -chi_even / d_m
# alpha_s scaled
zero_param_formulas['alpha_s*alt'] = ALPHA_S * alt_sign
zero_param_formulas['-alpha_s*alt'] = -ALPHA_S * alt_sign
# Dimension weighted
zero_param_formulas['1/(a1*d)'] = 1.0 / (5 * d_m)
zero_param_formulas['alt/(a1*d)'] = alt_sign / (5 * d_m)
# Heat kernel: exp(-t*lam)
for t_val in [0.1, 0.5, 1.0]:
    zero_param_formulas[f'exp(-{t_val}*lam)'] = np.exp(-t_val * lam)
    zero_param_formulas[f'-exp(-{t_val}*lam)'] = -np.exp(-t_val * lam)
# Spectral zeta: 1/lam^s
for s_val in [0.5, 1.5, 2.0]:
    zero_param_formulas[f'1/lam^{s_val}'] = 1.0 / lam**s_val
    zero_param_formulas[f'-1/lam^{s_val}'] = -1.0 / lam**s_val
# Green function: (lam - lam_mean)
lam_mean = np.mean(lam)
zero_param_formulas['lam-mean'] = lam - lam_mean
zero_param_formulas['mean-lam'] = lam_mean - lam
# Galois: phi-related
zero_param_formulas['phi/lam'] = PHI / lam
zero_param_formulas['1/(phi*lam)'] = 1.0 / (PHI * lam)
# alpha_s / lam
zero_param_formulas['alpha_s/lam'] = ALPHA_S / lam
zero_param_formulas['-alpha_s/lam'] = -ALPHA_S / lam
# 1/(N*lam) where N=120
zero_param_formulas['1/(N*lam)'] = 1.0 / (120 * lam)
# Physical: alpha_s * chi/d
zero_param_formulas['alpha_s*chi/d'] = ALPHA_S * chi_even / d_m

# For EACH formula, find the best overall scale AND best assignment
# Since pred is linear in c, the optimal scale is:
# scale = (delta_exp . (M@c)) / ||M@c||^2
print(f"\n{'Formula':>20} {'Best RMS':>10} {'Assignment':>60}")
print("-" * 95)

best_zero = (999, None, None, None, None)
for fname, c_raw in zero_param_formulas.items():
    best_this = (999, None, None)
    for perm in permutations(remaining_irreps):
        assignment = [0] + list(perm)
        M = np.zeros((9, 4))
        for f in range(9):
            k = assignment[f]
            for idx, m_idx in enumerate(even_modes):
                M[f, idx] = fermion_T3[f] * fermion_Nc[f] * CG_self[k, m_idx] * dims[m_idx]
        raw_pred = M @ c_raw
        # Optimal overall scale
        denom = np.dot(raw_pred, raw_pred)
        if denom > 1e-20:
            scale = np.dot(delta_exp, raw_pred) / denom
        else:
            scale = 0.0
        pred = scale * raw_pred
        rms = np.sqrt(np.mean((delta_exp - pred)**2))
        if rms < best_this[0]:
            best_this = (rms, list(assignment), scale)
    rms, assign, scale = best_this
    assign_str = " ".join(f"{fermion_names[f]}:{assign[f]}" for f in range(9))
    print(f"{fname:>20} {rms:>10.4f}  {assign_str}")
    if rms < best_zero[0]:
        best_zero = (rms, assign, scale, fname, c_raw)

# Show best zero-param in detail
rms, assignment, scale, fname, c_raw = best_zero
c_scaled = scale * c_raw
M = np.zeros((9, 4))
for f in range(9):
    k = assignment[f]
    for idx, m_idx in enumerate(even_modes):
        M[f, idx] = fermion_T3[f] * fermion_Nc[f] * CG_self[k, m_idx] * dims[m_idx]
pred = M @ c_scaled
print(f"\nBest zero-param formula: {fname} (with optimal scale {scale:.6f})")
print(f"  RMS = {rms:.4f}")
print(f"  Assignment: ", end="")
for f in range(9):
    print(f"{fermion_names[f]}->rho_{assignment[f]} ", end="")
print(f"\n  c_m = {np.round(c_scaled, 6)}")
for f in range(9):
    print(f"  {fermion_names[f]:>4}: pred={pred[f]:>+.4f}, exp={delta_exp[f]:>+.4f}, err={pred[f]-delta_exp[f]:>+.4f}")

# Best single-function result
print(f"\n{'='*72}")
print(f"SUMMARY: Best results by number of parameters")
print(f"{'='*72}")

# Collect all results
all_results = [
    ("T3*Nc*alpha_s/ln(phi)", 0, 0.169),
    ("T3*Nc*F(r=0) CG", 0, 0.142),
    ("T3*x quarks", 1, 0.127),
]

# Add best 1-param
best_1p = min(results, key=lambda x: x[3])
all_results.append((f"1p: {best_1p[0]}", 1, best_1p[3]))

# Add best 2-param
best_2p_item = sorted_2p[0]
all_results.append((f"2p: {best_2p_item[0]}", 2, best_2p_item[1][0]))

# 4 free modes (Model F from exp338b)
all_results.append(("Model F (4 modes)", 4, 0.033))

# Bare
all_results.append(("Bare formula", 0, 0.247))

all_results.sort(key=lambda x: (x[1], x[2]))
print(f"\n{'Model':>30} {'Params':>7} {'RMS':>8}")
print("-"*50)
for name, npar, rms in all_results:
    print(f"{name:>30} {npar:>7} {rms:>8.4f}")

print("\n\nDone.")
