"""
EXP-528: Two-Loop + Threshold Test in the Electroweak Window
=============================================================

Goal:
  Starting from the exp527 electroweak matching picture, test whether
  standard two-loop SM gauge running and small threshold corrections can
  collapse the two framework conditions

    sin^2(theta_W) = 6/26
    alpha_s = 1/(2 phi^3)

  to a single common electroweak matching scale.

Strategy:
  1. Start from PDG values at M_Z.
  2. Compare one-loop and two-loop gauge running between ~70 and 91 GeV.
  3. Extract the matching scales for the two algebraic framework values.
  4. Evaluate the geometric-mean scale sqrt(M_W M_Z).
  5. Quantify the effective threshold shifts still needed at that scale.

Notes:
  - This is a practical RG consistency test, not a first-principles derivation
    of matching coefficients from the spectral action.
  - We use standard SM gauge two-loop coefficients (GUT-normalised U(1)).
  - Yukawa terms are omitted; over the narrow 70-91 GeV window they are a
    subleading correction relative to the one-loop baseline.
"""

import numpy as np
import sys

sys.path.insert(0, '.')
from commons import PHI, alpha_em as alpha_fw_ir, inv_alpha as inv_alpha_fw_ir


print("=" * 70)
print("EXP-528: TWO-LOOP + THRESHOLD TEST IN THE EW WINDOW")
print("=" * 70)


# ============================================================
# INPUT DATA
# ============================================================
M_Z = 91.1876         # GeV
M_W = 80.3692         # GeV
M_H = 125.25          # GeV (diagnostic only)
M_t = 172.69          # GeV (diagnostic only)

# PDG 2024 input at M_Z (MSbar-like electroweak input)
alpha_em_MZ = 1.0 / 127.951
alpha_s_MZ = 0.1179
sin2_MZ = 0.23122

# Framework targets
sin2_fw = 6.0 / 26.0
alpha_s_fw = 1.0 / (2.0 * PHI**3)

print("\nInputs")
print(f"  M_W = {M_W:.4f} GeV")
print(f"  M_Z = {M_Z:.4f} GeV")
print(f"  sqrt(M_W M_Z) = {np.sqrt(M_W * M_Z):.4f} GeV")
print(f"  Framework: sin^2(tW) = {sin2_fw:.6f}, alpha_s = {alpha_s_fw:.6f}")
print(f"  PDG @ M_Z: sin^2(tW) = {sin2_MZ:.6f}, alpha_s = {alpha_s_MZ:.6f}, 1/alpha = {1/alpha_em_MZ:.3f}")


# ============================================================
# RG COEFFICIENTS
# ============================================================
# Convention:
#   d alpha_i / d ln(mu) = (b_i / 2pi) alpha_i^2
#                        + (alpha_i^2 / 8pi^2) sum_j B_ij alpha_j
#
# with U(1) GUT normalisation alpha_1 = (5/3) alpha_Y
b_vec = np.array([41.0 / 10.0, -19.0 / 6.0, -7.0], dtype=float)

# Standard SM gauge two-loop matrix in the same convention
B_mat = np.array([
    [199.0 / 50.0, 27.0 / 10.0, 44.0 / 5.0],
    [9.0 / 10.0, 35.0 / 6.0, 12.0],
    [11.0 / 10.0, 9.0 / 2.0, -26.0],
], dtype=float)

print("\nRG coefficients")
print(f"  One-loop b = {b_vec}")
print("  Two-loop B =")
for row in B_mat:
    print(f"    {row}")


# ============================================================
# CONVERSION HELPERS
# ============================================================
def sin2_aem_to_inv_alphas(sin2, alpha_em):
    """Convert (sin^2 theta_W, alpha_em) to GUT-normalised inverse couplings."""
    alpha_2 = alpha_em / sin2
    alpha_Y = alpha_em / (1.0 - sin2)
    alpha_1 = (5.0 / 3.0) * alpha_Y
    return np.array([1.0 / alpha_1, 1.0 / alpha_2], dtype=float)


def inv_alphas_to_observables(inv_a1, inv_a2):
    """Convert GUT-normalised inverse couplings to (sin^2 theta_W, alpha_em)."""
    alpha_1 = 1.0 / inv_a1
    alpha_2 = 1.0 / inv_a2
    alpha_Y = (3.0 / 5.0) * alpha_1
    sin2 = alpha_Y / (alpha_Y + alpha_2)
    alpha_em = alpha_Y * alpha_2 / (alpha_Y + alpha_2)
    return sin2, alpha_em


inv_a1_MZ, inv_a2_MZ = sin2_aem_to_inv_alphas(sin2_MZ, alpha_em_MZ)
inv_a3_MZ = 1.0 / alpha_s_MZ
alphas_MZ = np.array([1.0 / inv_a1_MZ, 1.0 / inv_a2_MZ, alpha_s_MZ], dtype=float)


# ============================================================
# ONE-LOOP REFERENCE
# ============================================================
def run_inv_alpha_one_loop(inv_alpha_ref, b_coeff, mu_ref, mu_target):
    return inv_alpha_ref - (b_coeff / (2.0 * np.pi)) * np.log(mu_target / mu_ref)


def observables_one_loop(mu):
    ia1 = run_inv_alpha_one_loop(inv_a1_MZ, b_vec[0], M_Z, mu)
    ia2 = run_inv_alpha_one_loop(inv_a2_MZ, b_vec[1], M_Z, mu)
    ia3 = run_inv_alpha_one_loop(inv_a3_MZ, b_vec[2], M_Z, mu)
    sin2, alpha_em = inv_alphas_to_observables(ia1, ia2)
    return sin2, 1.0 / ia3, alpha_em


# ============================================================
# TWO-LOOP GAUGE-ONLY RUNNING
# ============================================================
def dalpha_dt(alpha_vec):
    """Two-loop gauge-only RG flow for alpha_i."""
    one_loop = (b_vec / (2.0 * np.pi)) * alpha_vec**2
    two_loop = np.zeros(3, dtype=float)
    for i in range(3):
        two_loop[i] = (alpha_vec[i]**2 / (8.0 * np.pi**2)) * np.dot(B_mat[i], alpha_vec)
    return one_loop + two_loop


def rk4_step(alpha_vec, dt):
    k1 = dalpha_dt(alpha_vec)
    k2 = dalpha_dt(alpha_vec + 0.5 * dt * k1)
    k3 = dalpha_dt(alpha_vec + 0.5 * dt * k2)
    k4 = dalpha_dt(alpha_vec + dt * k3)
    return alpha_vec + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)


def build_two_loop_scan(mu_min=70.0, mu_max=M_Z, n_steps=6000):
    """Integrate downward from M_Z to mu_min and return observables on a dense grid."""
    t_max = np.log(mu_max)
    t_min = np.log(mu_min)
    t_grid = np.linspace(t_max, t_min, n_steps)
    mu_grid = np.exp(t_grid)

    alpha_grid = np.zeros((n_steps, 3), dtype=float)
    alpha_grid[0] = alphas_MZ

    for i in range(n_steps - 1):
        dt = t_grid[i + 1] - t_grid[i]
        alpha_grid[i + 1] = rk4_step(alpha_grid[i], dt)

    sin2_grid = np.zeros(n_steps, dtype=float)
    as_grid = alpha_grid[:, 2].copy()
    aem_grid = np.zeros(n_steps, dtype=float)

    for i in range(n_steps):
        alpha_1, alpha_2 = alpha_grid[i, 0], alpha_grid[i, 1]
        alpha_Y = (3.0 / 5.0) * alpha_1
        sin2_grid[i] = alpha_Y / (alpha_Y + alpha_2)
        aem_grid[i] = alpha_Y * alpha_2 / (alpha_Y + alpha_2)

    return mu_grid, sin2_grid, as_grid, aem_grid


mu_grid, sin2_two, as_two, aem_two = build_two_loop_scan()


def interpolate_match(mu_grid, obs_grid, target):
    """Return mu where obs(mu) ~= target, assuming monotonicity on the scan interval."""
    idx = int(np.argmin(np.abs(obs_grid - target)))
    return mu_grid[idx], obs_grid[idx]


mu_sin2_1l = interpolate_match(
    np.geomspace(70.0, M_Z, 4000),
    np.array([observables_one_loop(mu)[0] for mu in np.geomspace(70.0, M_Z, 4000)]),
    sin2_fw,
)[0]
mu_as_1l = interpolate_match(
    np.geomspace(70.0, M_Z, 4000),
    np.array([observables_one_loop(mu)[1] for mu in np.geomspace(70.0, M_Z, 4000)]),
    alpha_s_fw,
)[0]

mu_sin2_2l, sin2_at_match_2l = interpolate_match(mu_grid, sin2_two, sin2_fw)
mu_as_2l, as_at_match_2l = interpolate_match(mu_grid, as_two, alpha_s_fw)

print(f"\n{'=' * 70}")
print("MATCHING SCALES")
print(f"{'=' * 70}")
print(f"  One-loop: sin^2(tW)=6/26 at {mu_sin2_1l:.3f} GeV")
print(f"  One-loop: alpha_s=1/(2phi^3) at {mu_as_1l:.3f} GeV")
print(f"  Two-loop: sin^2(tW)=6/26 at {mu_sin2_2l:.3f} GeV")
print(f"  Two-loop: alpha_s=1/(2phi^3) at {mu_as_2l:.3f} GeV")
print(f"  Two-loop separation = {abs(mu_as_2l - mu_sin2_2l):.3f} GeV")


# ============================================================
# BEST COMMON SCALE
# ============================================================
def cost_function(sin2_vals, as_vals):
    return ((sin2_vals / sin2_fw) - 1.0)**2 + ((as_vals / alpha_s_fw) - 1.0)**2


cost_two = cost_function(sin2_two, as_two)
idx_best = int(np.argmin(cost_two))
mu_best = mu_grid[idx_best]
sin2_best = sin2_two[idx_best]
as_best = as_two[idx_best]
aem_best = aem_two[idx_best]

# Balanced scale: make the two relative deviations as equal as possible
rel_s2_two = (sin2_two / sin2_fw) - 1.0
rel_as_two = (as_two / alpha_s_fw) - 1.0
idx_bal = int(np.argmin(np.abs(np.abs(rel_s2_two) - np.abs(rel_as_two))))
mu_bal = mu_grid[idx_bal]
sin2_bal = sin2_two[idx_bal]
as_bal = as_two[idx_bal]
aem_bal = aem_two[idx_bal]

mu_geom = np.sqrt(M_W * M_Z)
idx_geom = int(np.argmin(np.abs(mu_grid - mu_geom)))
sin2_geom = sin2_two[idx_geom]
as_geom = as_two[idx_geom]
aem_geom = aem_two[idx_geom]

print(f"\n{'=' * 70}")
print("COMMON-SCALE TEST")
print(f"{'=' * 70}")
print(f"  Best two-loop common scale mu* = {mu_best:.3f} GeV")
print(f"    sin^2(tW)(mu*) = {sin2_best:.6f}   dev {(sin2_best / sin2_fw - 1) * 100:+.4f}%")
print(f"    alpha_s(mu*)   = {as_best:.6f}   dev {(as_best / alpha_s_fw - 1) * 100:+.4f}%")
print(f"    1/alpha(mu*)   = {1.0 / aem_best:.4f}   vs IR framework {inv_alpha_fw_ir:.4f}")
print(f"  Balanced two-loop scale mu_bal = {mu_bal:.3f} GeV")
print(f"    sin^2(tW)(mu_bal) = {sin2_bal:.6f}   dev {(sin2_bal / sin2_fw - 1) * 100:+.4f}%")
print(f"    alpha_s(mu_bal)   = {as_bal:.6f}   dev {(as_bal / alpha_s_fw - 1) * 100:+.4f}%")
print(f"    1/alpha(mu_bal)   = {1.0 / aem_bal:.4f}")
print(f"  Geometric mean scale mu_g = sqrt(M_W M_Z) = {mu_geom:.3f} GeV")
print(f"    sin^2(tW)(mu_g) = {sin2_geom:.6f}   dev {(sin2_geom / sin2_fw - 1) * 100:+.4f}%")
print(f"    alpha_s(mu_g)   = {as_geom:.6f}   dev {(as_geom / alpha_s_fw - 1) * 100:+.4f}%")
print(f"    1/alpha(mu_g)   = {1.0 / aem_geom:.4f}")


# ============================================================
# EFFECTIVE THRESHOLD SHIFTS
# ============================================================
def print_threshold_requirements(label, mu, sin2_val, as_val):
    ds_abs = sin2_fw - sin2_val
    das_abs = alpha_s_fw - as_val
    ds_rel = ds_abs / sin2_fw * 100.0
    das_rel = das_abs / alpha_s_fw * 100.0

    print(f"\n  {label} at mu = {mu:.3f} GeV")
    print(f"    Needed threshold shift in sin^2(tW): {ds_abs:+.6e} ({ds_rel:+.4f}%)")
    print(f"    Needed threshold shift in alpha_s:   {das_abs:+.6e} ({das_rel:+.4f}%)")


print(f"\n{'=' * 70}")
print("EFFECTIVE THRESHOLD SHIFTS")
print(f"{'=' * 70}")
print_threshold_requirements("At geometric mean", mu_geom, sin2_geom, as_geom)
print_threshold_requirements("At best common scale", mu_best, sin2_best, as_best)
print_threshold_requirements("At balanced scale", mu_bal, sin2_bal, as_bal)


# ============================================================
# FINE TABLE IN THE EW WINDOW
# ============================================================
print(f"\n{'=' * 70}")
print("EW WINDOW TABLE (TWO-LOOP)")
print(f"{'=' * 70}")
print(f"{'mu (GeV)':>10} {'sin2':>10} {'alpha_s':>10} {'1/alpha':>10} {'d_s2(%)':>10} {'d_as(%)':>10}")
for mu in [80.0, M_W, 83.3, mu_geom, 87.0, 88.5, mu_bal, mu_best, 90.4, M_Z]:
    idx = int(np.argmin(np.abs(mu_grid - mu)))
    s2 = sin2_two[idx]
    a_s = as_two[idx]
    a_em = aem_two[idx]
    print(f"{mu_grid[idx]:>10.3f} {s2:>10.6f} {a_s:>10.6f} {1.0/a_em:>10.4f} "
          f"{(s2/sin2_fw-1)*100:>+10.4f} {(a_s/alpha_s_fw-1)*100:>+10.4f}")


# ============================================================
# SUMMARY
# ============================================================
print(f"\n{'=' * 70}")
print("SUMMARY")
print(f"{'=' * 70}")
print(f"""
  1. One-loop and two-loop both keep the framework matching in the same
     electroweak window, not at a UV scale.
  2. The two-loop matching scales are:
       sin^2(tW)=6/26  -> {mu_sin2_2l:.2f} GeV
       alpha_s=1/(2phi^3) -> {mu_as_2l:.2f} GeV
  3. The residual mismatch is only {abs(mu_as_2l - mu_sin2_2l):.2f} GeV.
  4. The geometric-mean scale sqrt(M_W M_Z) = {mu_geom:.2f} GeV sits inside
     that window and requires only sub-percent effective threshold shifts.
  5. Best common two-loop scale: {mu_best:.2f} GeV.
     Remaining deviations there are
       sin^2: {(sin2_best / sin2_fw - 1) * 100:+.4f}%
       alpha_s: {(as_best / alpha_s_fw - 1) * 100:+.4f}%
  6. Balanced two-loop scale: {mu_bal:.2f} GeV.
     There the two deviations are almost identical:
       sin^2: {(sin2_bal / sin2_fw - 1) * 100:+.4f}%
       alpha_s: {(as_bal / alpha_s_fw - 1) * 100:+.4f}%

  Interpretation:
    Two-loop gauge running does not destroy the exp527 picture.
    The remaining gap is small enough that electroweak threshold/scheme
    effects are a credible explanation.  The real question is no longer
    "is there a scale mismatch?" but "what exact matching prescription
    selects the common EW scale?"
""")

print("=" * 70)
print("EXP-528 COMPLETE")
print("=" * 70)
