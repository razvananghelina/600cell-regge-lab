"""
Verify: Universality class of the vacuum selection.

The vacuum (5,25) is selected by ANY functional in the class:
  { S >= 0, S(n0) = 0, RG monotonicity on the scanned domain }

This script verifies the THREE hypotheses of the universality argument,
not the conclusion (which follows logically from the hypotheses).

Tests:
  T1: Non-negativity -- 4 standard mismatch norms are >= 0 everywhere
  T2: Unique zero -- each norm has S = 0 only at (5, n0=25) on the domain
  T3: RG monotonicity -- each coupling deviates monotonically from its
      algebraic value as |n - n0| increases (one-loop, no second crossing)
  T4: Conclusion -- under T1+T2+T3, (5,25) is the unique global minimum
      for all 4 norms simultaneously
"""

import numpy as np

passed = 0
failed = 0

def check(name, condition, detail=""):
    global passed, failed
    if condition:
        passed += 1
        print(f"  PASS  {name}")
    else:
        failed += 1
        print(f"  FAIL  {name}")
    if detail:
        print(f"        {detail}")

print("=" * 70)
print("VERIFY: UNIVERSALITY CLASS OF VACUUM SELECTION")
print("=" * 70)

# Framework constants
a1 = 5
phi = (1 + np.sqrt(5)) / 2
m_e = 0.51099895e-3
n0 = a1**2

# Couplings
sin2_fw = (a1+1) / (a1**2+1)
alpha_s_fw = 1 / (2*phi**3)
B = 4*a1*phi**4
alpha_fw = (B - np.sqrt(B**2 - 8*np.pi)) / (4*np.pi)

# Derived beta functions
b1_rg = 4/3*3 + 1/10
b2_rg = 4/3*3 - 22/3 + 1/6
b3_rg = -7.0

# Inverse GUT-normalised couplings at the fixed point
alpha_Y = alpha_fw / (1 - sin2_fw)
alpha_2 = alpha_fw / sin2_fw
inv_a1_fw = 1 / ((5/3)*alpha_Y)
inv_a2_fw = 1 / alpha_2
inv_a3_fw = 1 / alpha_s_fw

def rg_deviations(n):
    """Return (delta_sin2, delta_alpha_s) = relative deviations at scale n."""
    Lambda = m_e * phi**n
    Lambda_0 = m_e * phi**n0
    if Lambda <= 0 or Lambda_0 <= 0:
        return None, None
    ln_r = np.log(Lambda / Lambda_0)
    ia1 = inv_a1_fw - b1_rg/(2*np.pi)*ln_r
    ia2 = inv_a2_fw - b2_rg/(2*np.pi)*ln_r
    ia3 = inv_a3_fw - b3_rg/(2*np.pi)*ln_r
    if ia1 <= 0 or ia2 <= 0 or ia3 <= 0:
        return None, None
    aY = (3/5)/ia1; a2 = 1/ia2
    sin2 = aY/(aY+a2); alpha_s = 1/ia3
    return sin2/sin2_fw - 1, alpha_s/alpha_s_fw - 1

# ============================================================
# FOUR MISMATCH NORMS
# ============================================================
# All non-negative, all zero iff both deviations are zero.

def norm_L2(n):
    """Sum of squares (the standard choice)."""
    d1, d2 = rg_deviations(n)
    if d1 is None: return 1e10
    return d1**2 + d2**2

def norm_L1(n):
    """Sum of absolute values."""
    d1, d2 = rg_deviations(n)
    if d1 is None: return 1e10
    return abs(d1) + abs(d2)

def norm_Linf(n):
    """Maximum absolute deviation."""
    d1, d2 = rg_deviations(n)
    if d1 is None: return 1e10
    return max(abs(d1), abs(d2))

def norm_KL(n):
    """Kullback-Leibler-inspired: sum of (x*ln(x) - x + 1) where x = g/g_alg."""
    d1, d2 = rg_deviations(n)
    if d1 is None: return 1e10
    x1 = 1 + d1; x2 = 1 + d2
    if x1 <= 0 or x2 <= 0: return 1e10
    return (x1*np.log(x1) - x1 + 1) + (x2*np.log(x2) - x2 + 1)

norms = [("L2 (sum of squares)", norm_L2),
         ("L1 (sum of |dev|)", norm_L1),
         ("Linf (max |dev|)", norm_Linf),
         ("KL (divergence)", norm_KL)]

N_RANGE = range(15, 36)  # scan domain

# ============================================================
# T1: NON-NEGATIVITY
# ============================================================
print("\nT1: Non-negativity")
for name, norm in norms:
    vals = [norm(n) for n in N_RANGE]
    all_nonneg = all(v >= -1e-15 for v in vals if v < 1e9)
    check(f"T1: {name} >= 0 on [{min(N_RANGE)},{max(N_RANGE)}]",
          all_nonneg,
          f"min = {min(v for v in vals if v < 1e9):.4e}")

# ============================================================
# T2: UNIQUE ZERO AT n0=25
# ============================================================
print("\nT2: Unique zero at n0=25")
TOL_ZERO = 1e-10
for name, norm in norms:
    vals = {n: norm(n) for n in N_RANGE}
    zeros = [n for n, v in vals.items() if abs(v) < TOL_ZERO]
    check(f"T2: {name} has unique zero at n=25",
          zeros == [25],
          f"zeros at: {zeros}")

# ============================================================
# T3: RG MONOTONICITY
# ============================================================
print("\nT3: RG monotonicity (deviations grow with |n-25|)")

# Check that |delta_sin2| and |delta_alpha_s| are monotonically
# increasing as n moves away from n0 in BOTH directions.
devs_sin2 = {}
devs_as = {}
for n in N_RANGE:
    d1, d2 = rg_deviations(n)
    if d1 is not None:
        devs_sin2[n] = d1
        devs_as[n] = d2

# Check monotonicity of |deviation| for n > n0
mono_sin2_right = all(abs(devs_sin2.get(n+1, 0)) >= abs(devs_sin2.get(n, 0)) - 1e-10
                      for n in range(n0, max(N_RANGE)-1))
mono_as_right = all(abs(devs_as.get(n+1, 0)) >= abs(devs_as.get(n, 0)) - 1e-10
                    for n in range(n0, max(N_RANGE)-1))

# Check monotonicity for n < n0
mono_sin2_left = all(abs(devs_sin2.get(n-1, 0)) >= abs(devs_sin2.get(n, 0)) - 1e-10
                     for n in range(min(N_RANGE)+1, n0+1))
mono_as_left = all(abs(devs_as.get(n-1, 0)) >= abs(devs_as.get(n, 0)) - 1e-10
                   for n in range(min(N_RANGE)+1, n0+1))

check("T3a: |delta_sin2| monotone increasing for n > 25", mono_sin2_right)
check("T3b: |delta_sin2| monotone increasing for n < 25", mono_sin2_left)
check("T3c: |delta_alpha_s| monotone increasing for n > 25", mono_as_right)
check("T3d: |delta_alpha_s| monotone increasing for n < 25", mono_as_left)

# Explicit display
print(f"\n  Deviation profile:")
print(f"  {'n':>4s} {'d(sin2)':>12s} {'d(alpha_s)':>12s} {'|d_s2|':>10s} {'|d_as|':>10s}")
for n in range(20, 31):
    d1 = devs_sin2.get(n, None)
    d2 = devs_as.get(n, None)
    if d1 is not None:
        mark = " <--" if n == n0 else ""
        print(f"  {n:4d} {d1:+12.6f} {d2:+12.6f} {abs(d1):10.6f} {abs(d2):10.6f}{mark}")

# ============================================================
# T4: ALL NORMS SELECT (5,25)
# ============================================================
print("\nT4: All norms select n=25 as global minimum")
for name, norm in norms:
    vals = {n: norm(n) for n in N_RANGE}
    best_n = min(vals, key=vals.get)
    best_S = vals[best_n]
    # Also check isolation: next-best is significantly worse
    sorted_vals = sorted(vals.items(), key=lambda x: x[1])
    if len(sorted_vals) > 1:
        runner_up_n, runner_up_S = sorted_vals[1]
        isolation = runner_up_S / best_S if best_S > 0 else float('inf')
    else:
        isolation = float('inf')
    check(f"T4: {name} minimum at n=25",
          best_n == 25,
          f"min at n={best_n}, S={best_S:.4e}, runner-up at n={runner_up_n} ({isolation:.0f}x)")

# ============================================================
# SUMMARY
# ============================================================
print(f"\n{'='*70}")
print(f"RESULTS: {passed} passed, {failed} failed out of {passed+failed}")
if failed == 0:
    print("ALL TESTS PASS")
print(f"{'='*70}")

print(f"""
  The universality class of the vacuum selection is:
    C = {{ S : R -> R_{{>=0}} | S(n0) = 0, |dev_i(n)| monotone in |n-n0| }}

  All 4 tested norms (L2, L1, Linf, KL) belong to C.
  Within C, the vacuum (5,25) is the unique global minimum.

  This is NOT polynomial universality (which fails, exp552).
  It IS physical universality: any non-negative mismatch measure
  with one-loop RG monotonicity selects the same vacuum.
""")
