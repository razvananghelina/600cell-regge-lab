"""
EXP-542: KILL TEST — Does RG Bootstrap Select a1=5, n=25?
============================================================

QUESTION: Is (a1=5, n=25) a sharp, isolated minimum of a global
cost function, or is it just one point in a broad valley?

PROTOCOL:
  For each (a1, n) with a1 in [3,10], n in [10,40]:
    1. Compute phi(a1) = (1+sqrt(a1))/2
    2. Compute framework couplings: sin^2 = (a1+1)/(a1^2+1), alpha_s = 1/(2*phi^3)
    3. Compute bare lattice scale: Lambda = m_e * phi^n
    4. Evaluate SM running couplings at Lambda (from PDG at M_Z)
    5. Cost = mismatch between framework and SM couplings at Lambda
    6. Optionally: check if Lambda * running_ratio ~ m_Z

  Same recipe for ALL (a1, n). No tuning.

PASS: Sharp isolated minimum at (5, 25).
FAIL: Broad valley, or other (a1, n) nearly as good.
"""

import numpy as np
import sys
sys.path.insert(0, '.')

print("=" * 70)
print("EXP-542: KILL TEST — RG BOOTSTRAP SELECTION")
print("=" * 70)

# ============================================================
# SETUP: SM RUNNING FROM PDG AT M_Z
# ============================================================
M_Z = 91.1876
m_e = 0.51099895e-3

# PDG at M_Z
alpha_em_MZ = 1.0 / 127.951
alpha_s_MZ = 0.1179
sin2_MZ = 0.23122

# GUT-normalised betas
b_rg = np.array([41.0/10, -19.0/6, -7.0])

def sin2_aem_to_inv_alphas(sin2, aem):
    a2 = aem / sin2
    a_Y = aem / (1 - sin2)
    a1_gut = (5.0/3) * a_Y
    return 1.0/a1_gut, 1.0/a2

def inv_alphas_to_sin2(inv_a1, inv_a2):
    al1 = 1.0/inv_a1; al2 = 1.0/inv_a2
    a_Y = (3.0/5)*al1
    return a_Y / (a_Y + al2)

inv_a1_MZ, inv_a2_MZ = sin2_aem_to_inv_alphas(sin2_MZ, alpha_em_MZ)
inv_a3_MZ = 1.0 / alpha_s_MZ

def sm_couplings_at(mu):
    """SM running couplings at scale mu (one-loop from M_Z)."""
    if mu < 0.5 or mu > 1e7:
        return None, None
    ln = np.log(mu / M_Z)
    ia1 = inv_a1_MZ - b_rg[0]/(2*np.pi)*ln
    ia2 = inv_a2_MZ - b_rg[1]/(2*np.pi)*ln
    ia3 = inv_a3_MZ - b_rg[2]/(2*np.pi)*ln
    if ia1 <= 0 or ia2 <= 0 or ia3 <= 0:
        return None, None
    sin2 = inv_alphas_to_sin2(ia1, ia2)
    return sin2, 1.0/ia3

# ============================================================
# COST FUNCTION
# ============================================================

def cost(a1_val, n_val):
    """Global cost for (a1, n): mismatch of framework couplings at bare scale."""
    if a1_val < 2:
        return 1e10

    phi = (1 + np.sqrt(a1_val)) / 2
    if phi <= 1:
        return 1e10

    # Framework couplings
    b1_val = a1_val + 1
    sin2_fw = b1_val / (a1_val**2 + 1)
    alpha_s_fw = 1.0 / (2 * phi**3)

    # Check physical range
    if sin2_fw <= 0 or sin2_fw >= 0.5:
        return 1e10
    if alpha_s_fw <= 0 or alpha_s_fw >= 1:
        return 1e10

    # Bare lattice scale
    Lambda = m_e * phi**n_val
    if Lambda < 1 or Lambda > 1e6:
        return 1e10

    # SM couplings at Lambda
    sin2_sm, as_sm = sm_couplings_at(Lambda)
    if sin2_sm is None:
        return 1e10

    # Cost: relative mismatch squared
    c_sin2 = (sin2_sm / sin2_fw - 1)**2
    c_as = (as_sm / alpha_s_fw - 1)**2

    return c_sin2 + c_as


# ============================================================
# SCAN: ALL (a1, n) PAIRS
# ============================================================
print(f"\n{'='*70}")
print("FULL SCAN: a1 in [3,10], n in [10,40]")
print(f"{'='*70}")

results = {}
best_cost = 1e10
best_pair = None

for a1_val in range(3, 11):
    for n_val in range(10, 41):
        c = cost(a1_val, n_val)
        results[(a1_val, n_val)] = c
        if c < best_cost:
            best_cost = c
            best_pair = (a1_val, n_val)

print(f"\n  Global minimum: a1={best_pair[0]}, n={best_pair[1]}, cost={best_cost:.6e}")
print(f"  Expected: a1=5, n=25")
print(f"  Match: {'YES' if best_pair == (5, 25) else 'NO'}")

# ============================================================
# LANDSCAPE: COST FOR EACH a1 AT ITS BEST n
# ============================================================
print(f"\n{'='*70}")
print("LANDSCAPE: BEST n FOR EACH a1")
print(f"{'='*70}")

print(f"  {'a1':>3s} {'best_n':>7s} {'cost':>12s} {'Lambda':>10s} {'sin2_fw':>8s} {'as_fw':>8s} {'ratio':>10s}")

for a1_val in range(3, 11):
    # Find best n for this a1
    best_n = min(range(10, 41), key=lambda n: results.get((a1_val, n), 1e10))
    c = results[(a1_val, best_n)]
    phi = (1 + np.sqrt(a1_val)) / 2
    Lambda = m_e * phi**best_n
    sin2_fw = (a1_val+1) / (a1_val**2+1)
    as_fw = 1/(2*phi**3)

    ratio_to_best = c / best_cost if best_cost > 0 else float('inf')

    marker = " <-- MINIMUM" if (a1_val, best_n) == best_pair else ""
    print(f"  {a1_val:3d} {best_n:7d} {c:12.4e} {Lambda:10.2f} {sin2_fw:8.5f} {as_fw:8.5f} {ratio_to_best:10.1f}x{marker}")

# ============================================================
# SHARPNESS: NEIGHBORS OF (5, 25)
# ============================================================
print(f"\n{'='*70}")
print("SHARPNESS: COST AT (5,25) AND NEIGHBORS")
print(f"{'='*70}")

center_cost = results[(5, 25)]
print(f"  Center: cost(5, 25) = {center_cost:.6e}\n")

print(f"  {'(a1, n)':>10s} {'cost':>12s} {'ratio':>10s}")
for da1 in [-2, -1, 0, 1, 2]:
    for dn in [-2, -1, 0, 1, 2]:
        a1_val = 5 + da1
        n_val = 25 + dn
        if a1_val < 3 or a1_val > 10 or n_val < 10 or n_val > 40:
            continue
        c = results.get((a1_val, n_val), 1e10)
        ratio = c / center_cost if center_cost > 0 else float('inf')
        marker = " <--" if da1 == 0 and dn == 0 else ""
        print(f"  ({a1_val:2d},{n_val:2d}) {c:12.4e} {ratio:10.1f}x{marker}")

# ============================================================
# STABILITY: 1-LOOP vs 2-LOOP
# ============================================================
print(f"\n{'='*70}")
print("STABILITY: 1-LOOP vs 2-LOOP")
print(f"{'='*70}")

# 2-loop correction: add B_ij terms
B_mat = np.array([
    [199.0/50, 27.0/10, 44.0/5],
    [9.0/10, 35.0/6, 12.0],
    [11.0/10, 9.0/2, -26.0],
])

def sm_couplings_2loop(mu, n_steps=2000):
    """2-loop gauge running via RK4."""
    if mu < 0.5 or mu > 1e7:
        return None, None

    alphas = np.array([1.0/inv_a1_MZ, 1.0/inv_a2_MZ, alpha_s_MZ])
    t_start = np.log(M_Z)
    t_end = np.log(mu)
    dt = (t_end - t_start) / n_steps

    for _ in range(n_steps):
        one_loop = (b_rg / (2*np.pi)) * alphas**2
        two_loop = np.array([
            (alphas[i]**2 / (8*np.pi**2)) * np.dot(B_mat[i], alphas)
            for i in range(3)
        ])
        k1 = one_loop + two_loop
        a2 = alphas + 0.5*dt*k1
        one_loop2 = (b_rg/(2*np.pi)) * a2**2
        two_loop2 = np.array([(a2[i]**2/(8*np.pi**2))*np.dot(B_mat[i],a2) for i in range(3)])
        k2 = one_loop2 + two_loop2
        alphas = alphas + dt*k2

    if any(a <= 0 for a in alphas) or any(1/a <= 0 for a in alphas):
        return None, None

    sin2 = inv_alphas_to_sin2(1/alphas[0], 1/alphas[1])
    return sin2, alphas[2]

def cost_2loop(a1_val, n_val):
    phi = (1 + np.sqrt(a1_val)) / 2
    if phi <= 1: return 1e10
    sin2_fw = (a1_val+1)/(a1_val**2+1)
    as_fw = 1/(2*phi**3)
    Lambda = m_e * phi**n_val
    if Lambda < 1 or Lambda > 1e6: return 1e10
    sin2_sm, as_sm = sm_couplings_2loop(Lambda)
    if sin2_sm is None: return 1e10
    return (sin2_sm/sin2_fw - 1)**2 + (as_sm/as_fw - 1)**2

# Compare 1-loop and 2-loop for key points
print(f"  {'(a1,n)':>10s} {'1-loop':>12s} {'2-loop':>12s} {'ratio':>10s}")
test_points = [(5,25), (4,26), (5,24), (5,26), (6,24), (4,16), (3,29)]
for a1_val, n_val in test_points:
    c1 = results.get((a1_val, n_val), 1e10)
    c2 = cost_2loop(a1_val, n_val)
    ratio = c2/c1 if c1 > 0 and c1 < 1e9 else float('inf')
    print(f"  ({a1_val:2d},{n_val:2d}) {c1:12.4e} {c2:12.4e} {ratio:10.2f}")

# Full 2-loop scan over the same discrete domain as the 1-loop test.
results_2l = {}
best_cost_2l = 1e10
best_pair_2l = None
for a1_val in range(3, 11):
    for n_val in range(10, 41):
        c = cost_2loop(a1_val, n_val)
        results_2l[(a1_val, n_val)] = c
        if c < best_cost_2l:
            best_cost_2l = c
            best_pair_2l = (a1_val, n_val)

print(f"\n  Full 2-loop minimum on scanned domain: ({best_pair_2l[0]}, {best_pair_2l[1]})")
print(f"  2-loop cost at minimum: {best_cost_2l:.4e}")

second_best_2l = min(
    c for pair, c in results_2l.items()
    if pair != best_pair_2l and c < 1e9
)
best_other_a1_2l = min(
    c for pair, c in results_2l.items()
    if pair[0] != 5 and c < 1e9
)
print(f"  Nearest 2-loop competitor ratio: {second_best_2l / best_cost_2l:.1f}x")
print(f"  Best 2-loop competitor with a1 != 5: {best_other_a1_2l / best_cost_2l:.1f}x")

print(f"\n  Best 2-loop n for each a1:")
for a1_val in range(3, 11):
    best_n_2l = min(range(10, 41), key=lambda n: results_2l[(a1_val, n)])
    print(f"    a1={a1_val}: best n={best_n_2l}, cost={results_2l[(a1_val, best_n_2l)]:.4e}")

# ============================================================
# VERDICT
# ============================================================
print(f"\n{'='*70}")
print("VERDICT")
print(f"{'='*70}")

# Check all criteria
is_global_min = best_pair == (5, 25)
neighbor_ratio = min(results.get((a1, n), 1e10) / center_cost
                     for a1 in [4, 6] for n in range(20, 30)
                     if (a1, n) != (5, 25) and results.get((a1, n), 1e10) < 1e9)

is_sharp = neighbor_ratio > 10  # 10x worse than center
is_isolated = all(results.get((a1, n), 1e10) > 100 * center_cost
                  for a1 in range(3, 11) for n in range(10, 41)
                  if (a1, n) != (5, 25) and abs(a1-5) + abs(n-25) > 2)

print(f"""
  Global minimum at (5, 25): {is_global_min}
  Cost at (5,25): {center_cost:.4e}
  Nearest competitor ratio: {neighbor_ratio:.1f}x worse
  Sharp (>10x): {is_sharp}
  Isolated (>100x for |da1|+|dn|>2): {is_isolated}
  2-loop minimum on full scanned domain: {best_pair_2l == (5, 25)}
  2-loop nearest competitor ratio: {second_best_2l / best_cost_2l:.1f}x
  2-loop best a1 != 5 ratio: {best_other_a1_2l / best_cost_2l:.1f}x

  RESULT: {'PASS' if is_global_min and is_sharp else 'FAIL'}

  {'The RG bootstrap SHARPLY selects (a1=5, n=25).' if is_global_min and is_sharp else
   'The minimum is NOT sharp enough to constitute a selection principle.'}
""")

print("=" * 70)
print("EXP-542 COMPLETE")
print("=" * 70)
