"""
EXP-103: STG - TOATE Predictiile vs GR
========================================

Metrica STG: ds^2 = -dt^2/C^2 + C^2*(dr^2 + r^2*dOmega^2)
  C(r) = 1 + M/r = (r+M)/r

Metrica Schwarzschild: ds^2 = -(1-2M/r)dt^2 + dr^2/(1-2M/r) + r^2*dOmega^2

Comparam FIECARE predictie: STG vs GR vs Newton.
TOATE formulele STG sunt DERIVATE din metrica, zero fitting.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq

M = 1.0  # Masa in unitati geometrizate (G=c=1)

print("=" * 70)
print("EXP-103: STG vs GR - TOATE PREDICTIILE")
print("=" * 70)
print(f"Metrica STG: ds^2 = -dt^2/C^2 + C^2*dr^2, C = 1 + M/r")
print(f"Schwarzschild: ds^2 = -(1-2M/r)dt^2 + dr^2/(1-2M/r) + r^2*dOmega^2")
print(f"M = {M} (unitati geometrizate)")


# ==============================================================
# 1. ORIZONT DE EVENIMENTE
# ==============================================================

print("\n" + "=" * 70)
print("1. ORIZONT DE EVENIMENTE")
print("=" * 70)

print(f"""
  GR (Schwarzschild):
    g_tt = -(1 - 2M/r) = 0  la  r_S = 2M = {2*M}
    => ORIZONT la r = 2M. Sub orizont, spatiul si timpul se inverseaza.

  STG:
    g_tt = -1/C^2 = -r^2/(r+M)^2
    g_tt = 0  doar la  r = 0  (singularitatea punctuala)
    Pentru r > 0: g_tt < 0 INTOTDEAUNA
    => FARA ORIZONT. Lumina poate scapa de la orice r > 0.

  DERIVAT: direct din forma metricii. Nu exista solutie r > 0 pt g_tt = 0.

  Implicatie: in STG, "gaurile negre" sunt obiecte ultra-compacte
  dar fara orizont - lumina e foarte redshiftata dar nu captata complet.
""")


# ==============================================================
# 2. REDSHIFT GRAVITATIONAL
# ==============================================================

print("=" * 70)
print("2. REDSHIFT GRAVITATIONAL")
print("=" * 70)

def redshift_stg(r_emit, r_obs=np.inf):
    """STG: 1+z = C_emit/C_obs. DERIVAT din g_tt."""
    C_e = 1 + M / r_emit
    C_o = 1 + M / r_obs if r_obs < np.inf else 1.0
    return C_e / C_o - 1

def redshift_gr(r_emit, r_obs=np.inf):
    """GR: 1+z = sqrt(g_tt_obs / g_tt_emit)."""
    if r_emit <= 2 * M:
        return np.inf
    f_e = 1 - 2 * M / r_emit
    f_o = 1 - 2 * M / r_obs if r_obs < np.inf else 1.0
    return np.sqrt(f_o / f_e) - 1

print(f"\n  {'r/M':>8} {'z_STG':>10} {'z_GR':>10} {'z_Newton':>10} {'STG/GR':>10}")
print(f"  {'-'*8} {'-'*10} {'-'*10} {'-'*10} {'-'*10}")
for r in [100, 50, 20, 10, 5, 3, 2.5, 2.1, 2.0, 1.5, 1.0, 0.5]:
    z_s = redshift_stg(r * M)
    z_g = redshift_gr(r * M) if r > 2 else np.inf
    z_n = M / (r * M)  # Newton: z = GM/(rc^2) = M/r
    ratio = z_s / z_g if z_g != np.inf and z_g > 0 else np.nan
    print(f"  {r:8.1f} {z_s:10.4f} {z_g:10.4f} {z_n:10.4f} {ratio:10.4f}")

print(f"""
  La camp slab (r >> M): z_STG = z_GR = M/r (identice!)
  La r = 2M: z_STG = 0.50, z_GR = infinit (orizont GR)
  La r = M:  z_STG = 1.00, z_GR = N/A (sub orizont)

  DERIVAT: z_STG = C(r_emit) - 1 = M/r_emit
  Verificare observationala: Pound-Rebka (camp slab) - ambele OK.
  Diferenta la camp puternic (stele neutronice, quasi-BH).
""")


# ==============================================================
# 3. POTENTIALUL EFECTIV SI ISCO
# ==============================================================

print("=" * 70)
print("3. ISCO (Innermost Stable Circular Orbit)")
print("=" * 70)

def V_eff_stg(r, L):
    """STG effective potential for massive particle."""
    C = 1 + M / r
    return 1.0 / (2 * C**2) + L**2 / (2 * C**4 * r**2)

def V_eff_gr(r, L):
    """Schwarzschild effective potential."""
    if r <= 2 * M:
        return np.inf
    return (1 - 2*M/r) * (1 + L**2/r**2) / 2

def L_circular_stg(r):
    """Angular momentum for circular orbit in STG. DERIVAT."""
    # From dV_eff/dr = 0: L^2 = M*(r+M)^2 / (r - M)
    if r <= M:
        return np.inf
    return np.sqrt(M * (r + M)**2 / (r - M))

def L_circular_gr(r):
    """Angular momentum for circular orbit in Schwarzschild."""
    if r <= 3 * M:
        return np.inf
    return np.sqrt(M * r**2 / (r - 3*M))

# ISCO: find where d^2V_eff/dr^2 = 0 for circular orbit
# For STG: L^2 = M(r+M)^2/(r-M). ISCO when this is minimum of V_eff's V_eff.
# Numerically: find r where circular orbit becomes unstable

def find_isco_stg():
    """Find ISCO by checking stability of circular orbits."""
    # Circular orbit exists for r > M (photon sphere at r=M)
    # Check d²V_eff/dr² > 0 for stability
    dr = 0.001
    for r_test in np.arange(M + 0.01, 10*M, 0.001):
        L = L_circular_stg(r_test)
        if np.isinf(L):
            continue
        V_m = V_eff_stg(r_test - dr, L)
        V_0 = V_eff_stg(r_test, L)
        V_p = V_eff_stg(r_test + dr, L)
        d2V = (V_p - 2*V_0 + V_m) / dr**2
        if d2V > 0:
            return r_test
    return None

def find_isco_gr():
    """GR ISCO = 6M (analytic). Verify numerically."""
    dr = 0.001
    for r_test in np.arange(3*M + 0.01, 10*M, 0.001):
        L = L_circular_gr(r_test)
        if np.isinf(L):
            continue
        V_m = V_eff_gr(r_test - dr, L)
        V_0 = V_eff_gr(r_test, L)
        V_p = V_eff_gr(r_test + dr, L)
        d2V = (V_p - 2*V_0 + V_m) / dr**2
        if d2V > 0:
            return r_test
    return None

isco_stg = find_isco_stg()
isco_gr = find_isco_gr()

print(f"""
  STG - Orbita circulara: L^2 = M*(r+M)^2 / (r - M)
    DERIVAT: din dV_eff/dr = 0
    Sfera fotonilor (r minim): r = M = {M:.1f} (GR: r = 3M = {3*M:.1f})
    ISCO (numeric): r = {isco_stg:.3f}M  (GR: r = {isco_gr:.3f}M, analitic: 6M)

  GR - Orbita circulara: L^2 = M*r^2 / (r - 3M)
    Sfera fotonilor: r = 3M
    ISCO: r = 6M (analitic)
""")


# ==============================================================
# 4. DEFLEXA LUMINII
# ==============================================================

print("=" * 70)
print("4. DEFLEXA LUMINII")
print("=" * 70)

def light_deflection_ode(b, M, theory='stg'):
    """Compute light deflection by solving the Binet ODE as an IVP.

    Instead of integrating 1/sqrt (with singularity at u_max), solve:
      d²u/dφ² + u = RHS(u)
    as an IVP with u(0)=0, u'(0)=1/b (incoming from infinity at φ=0).
    Integrate until u returns to 0 (outgoing). Total angle = φ_exit.
    Deflection = φ_exit - π.

    STG: RHS = 2M(1+Mu)³/b²   (from orbit eq. (du/dφ)² = (1+Mu)⁴/b² - u²)
    GR:  RHS = 3Mu²              (from orbit eq. (du/dphi)^2 = 1/b^2 - u^2 + 2Mu^3)
    NOTE: For light (null geodesic), GR has NO M/b^2 term (that's for massive particles).
    """
    from scipy.integrate import solve_ivp

    if theory == 'stg':
        def rhs(phi, y):
            u, up = y
            duup = 2*M*(1 + M*u)**3 / b**2 - u
            return [up, duup]
    else:
        def rhs(phi, y):
            u, up = y
            duup = 3*M*u**2 - u  # GR Binet equation for NULL geodesic
            return [up, duup]

    # IVP: u(0)=0, u'(0)=1/b (incoming ray)
    y0 = [0.0, 1.0/b]

    # Event: u crosses zero again (outgoing ray)
    def u_zero(phi, y):
        return y[0]
    u_zero.terminal = True
    u_zero.direction = -1  # Detect u going from + to 0

    # Integrate from phi=0 to phi ~ 2*pi (generous upper bound)
    sol = solve_ivp(rhs, [0, 2*np.pi], y0, events=u_zero,
                    rtol=1e-12, atol=1e-14, max_step=0.01)

    if sol.t_events[0].size > 0:
        phi_exit = sol.t_events[0][0]
        return phi_exit - np.pi
    else:
        return np.inf  # Light captured


def light_deflection_stg(b, M=1.0):
    """Light deflection in STG metric."""
    # STG photon sphere: b_crit where circular photon orbit exists
    # At r = M: b_stg = (1+M/M)² * M = 4M (from L/E with C²r²)
    return light_deflection_ode(b, M, 'stg')


def light_deflection_gr(b, M=1.0):
    """Light deflection in Schwarzschild."""
    b_crit = 3 * np.sqrt(3) * M
    if b <= b_crit * 1.01:
        return np.inf
    return light_deflection_ode(b, M, 'gr')


print(f"\n  {'b/M':>8} {'d_STG(rad)':>14} {'d_GR(rad)':>14} {'STG/GR':>10} {'4M/b':>10}")
print(f"  {'-'*8} {'-'*14} {'-'*14} {'-'*10} {'-'*10}")

for b_ratio in [5, 10, 20, 50, 100, 500, 1000, 10000]:
    b = b_ratio * M
    d_stg = light_deflection_stg(b, M)
    d_gr = light_deflection_gr(b, M)
    d_weak = 4 * M / b
    ratio = d_stg / d_gr if 0 < d_gr < 10 else np.nan
    print(f"  {b_ratio:8d} {d_stg:14.6e} {d_gr:14.6e} {ratio:10.4f} {d_weak:10.6e}")

# Verify weak-field limit explicitly
print(f"\n  Verificare limita camp slab (b/M -> infinit):")
for b_ratio in [1000, 5000, 10000, 50000]:
    b = b_ratio * M
    d_stg = light_deflection_stg(b, M)
    ratio_to_4Mb = d_stg / (4*M/b)
    print(f"    b = {b_ratio:>6d}M:  d_STG/(4M/b) = {ratio_to_4Mb:.6f}")

print(f"""
  La camp slab (b >> M): d_STG = d_GR = 4M/b (identice!)
  Aceasta e predictia celebra a Einstein (1.75 arcsec pt Soare).
  STG reproduce EXACT rezultatul GR in limita camp slab.
  La camp puternic: difera (STG < GR in general).

  DERIVAT: din ecuatia geodezicei nule pe metrica STG.
  Ecuatia orbitei: d^2u/dphi^2 + u = 2M(1+Mu)^3/b^2
""")


# ==============================================================
# 5. PRECESIA ORBITALA (analog Mercur)
# ==============================================================

print("=" * 70)
print("5. PRECESIA ORBITALA (tip Mercur)")
print("=" * 70)

def orbit_precession(r0, e, theory='stg', n_orbits=5):
    """Simulate orbit and measure precession per orbit.
    r0: semi-latus rectum (average radius)
    e: eccentricity
    """
    # Start at perihelion
    r_peri = r0 * (1 - e)
    r_apo = r0 * (1 + e)

    if theory == 'stg':
        L = L_circular_stg(r0)
        def acc_r(r, L):
            C = 1 + M/r
            # Radial acceleration from effective potential
            # a_r = -dV_eff/dr = M*r/(r+M)^3 - L^2*(r-M)*r^2/(r+M)^5
            # Actually, compute numerically
            dr = 1e-6
            return -(V_eff_stg(r + dr, L) - V_eff_stg(r - dr, L)) / (2 * dr)
    else:
        L = L_circular_gr(r0)
        def acc_r(r, L):
            dr = 1e-6
            return -(V_eff_gr(r + dr, L) - V_eff_gr(r - dr, L)) / (2 * dr)

    # Use 2D orbit integration
    # Start at perihelion: (r_peri, 0), velocity purely tangential
    x, y = r_peri, 0.0
    vx = 0.0
    vy = L / (r_peri * (1 + M/r_peri)**2) if theory == 'stg' else L / r_peri

    # For eccentric orbit, adjust speed
    # v_peri for given L: L = C^2 * r * v_phi (STG) or L = r * v_phi (GR)
    if theory == 'stg':
        C_peri = 1 + M/r_peri
        vy = L / (C_peri**2 * r_peri)
    else:
        vy = L / r_peri

    # Now modify L to get the desired eccentricity
    # For Keplerian: L = sqrt(M * r0 * (1-e^2))
    # We compute energy at perihelion and apo, adjust
    # Actually, just use the circular L and give a radial kick for eccentricity
    # Simpler: set L for circular, then start at r_peri = r0*(1-e) with that L
    # The orbit will naturally be eccentric

    # Start at perihelion (r_peri, 0) with tangential velocity
    if theory == 'stg':
        v_circ_r0 = r0 * np.sqrt(M) / (r0 + M)**1.5
    else:
        v_circ_r0 = np.sqrt(M / r0)
    v_peri = v_circ_r0 * np.sqrt((1 + e) / (1 - e)) if e > 0 else v_circ_r0

    x, y = r_peri, 0.0
    vx, vy = 0.0, v_peri

    # Angular momentum (conserved in both theories)
    L = x * vy - y * vx  # = r_peri * v_peri

    # Estimate period
    if theory == 'stg':
        T_est = 2 * np.pi * (r0 + M)**1.5 / np.sqrt(M)
    else:
        T_est = 2 * np.pi * r0**1.5 / np.sqrt(M)
    dt = T_est * 5e-5  # Finer dt for precession accuracy
    n_steps = int(n_orbits * T_est / dt) + 1

    # Cumulative angle tracking + perihelion detection
    # This correctly counts full revolutions (unlike arctan2 alone)
    cumul_angle = 0.0
    prev_angle = np.arctan2(y, x)  # initial angle
    prev_r = np.sqrt(x**2 + y**2)
    pprev_r = prev_r
    perihelion_cumul_angles = []  # cumulative angle at each perihelion

    # Initial acceleration
    r = np.sqrt(x**2 + y**2)
    if theory == 'stg':
        F = M * r / (r + M)**3
        ax, ay = -F * x / r, -F * y / r
    else:
        F_newton = M / r**2
        F_gr_corr = 3 * M * L**2 / r**4
        F_total = F_newton + F_gr_corr
        ax, ay = -F_total * x / r, -F_total * y / r

    for step in range(n_steps):
        x += vx * dt + 0.5 * ax * dt**2
        y += vy * dt + 0.5 * ay * dt**2
        r = np.sqrt(x**2 + y**2)

        if theory == 'stg':
            F = M * r / (r + M)**3
            ax_new = -F * x / r
            ay_new = -F * y / r
        else:
            if r <= 2 * M:
                ax_new, ay_new = 0, 0
            else:
                F_newton = M / r**2
                F_gr_corr = 3 * M * L**2 / r**4
                F_total = F_newton + F_gr_corr
                ax_new = -F_total * x / r
                ay_new = -F_total * y / r

        vx += 0.5 * (ax + ax_new) * dt
        vy += 0.5 * (ay + ay_new) * dt
        ax, ay = ax_new, ay_new

        # Track cumulative angle
        curr_angle = np.arctan2(y, x)
        d_angle = curr_angle - prev_angle
        if d_angle > np.pi: d_angle -= 2 * np.pi
        elif d_angle < -np.pi: d_angle += 2 * np.pi
        cumul_angle += d_angle
        prev_angle = curr_angle

        # Detect perihelion (local minimum in r)
        if prev_r < pprev_r and prev_r < r:
            perihelion_cumul_angles.append(cumul_angle)

        pprev_r = prev_r
        prev_r = r

    # Compute precession from consecutive perihelion cumulative angles
    if len(perihelion_cumul_angles) >= 3:
        angles = np.array(perihelion_cumul_angles)
        d_angles = np.diff(angles)
        # Skip first interval (might be inaccurate)
        if len(d_angles) > 2:
            d_angles = d_angles[1:]
        # Each d_angle should be ~2pi (one revolution). Precession = excess over 2pi.
        precession_per_orbit = np.mean(d_angles) - 2 * np.pi
        return precession_per_orbit, len(perihelion_cumul_angles)
    return np.nan, 0


# ---- ANALYTIC STG PRECESSION ----
# STG orbit equation: d²u/dφ² + u = M / (L²(1+Mu)³)
# where u = 1/r, L = r²(dφ/dt) = conserved angular momentum
#
# For circular orbit at r0: L² = M*r0⁴/(r0+M)³
# Linearizing around u0 = 1/r0:
#   d²(δu)/dφ² + ω²δu = 0
#   where ω² = 1 + 3M/(r0+M)
#
# Since ω > 1, perihelion returns BEFORE full revolution.
# Precession per orbit: δφ = 2π/ω - 2π = 2π(1/ω - 1)
# At leading order (r0 >> M): δφ_STG ≈ -3πM/r0
# Compare GR: δφ_GR = +6πM/(r0(1-e²))
#
# CRITICAL: STG predicts RETROGRADE precession (opposite to GR)!

print(f"\n  --- Derivare analitica precesie STG ---")
print(f"  Ecuatia orbitei: d^2u/dphi^2 + u = M / (L^2*(1+Mu)^3)")
print(f"  Frecventa radiala: omega^2 = 1 + 3M/(r0+M) > 1")
print(f"  => Precesie RETROGRADA: dphi = 2*pi*(1/omega - 1) < 0")
print(f"  La ordinul 1: dphi_STG ~ -3*pi*M/r0  (GR: +6*pi*M/r0)")
print()

# Test precession at different distances with e=0.2
e_test = 0.2
print(f"  Eccentricitate e = {e_test}")
print(f"  {'r0/M':>8} {'prec_STG':>12} {'STG_analyt':>12} {'prec_GR':>12} {'GR_analyt':>12}")
print(f"  {'-'*8} {'-'*12} {'-'*12} {'-'*12} {'-'*12}")

for r0_ratio in [20, 50, 100, 200, 500]:
    r0 = r0_ratio * M
    p_stg, n_stg = orbit_precession(r0, e_test, 'stg', n_orbits=10)
    p_gr, n_gr = orbit_precession(r0, e_test, 'gr', n_orbits=10)
    # GR analytic: delta_phi = 6*pi*M / (r0*(1-e^2))
    p_gr_th = 6 * np.pi * M / (r0 * (1 - e_test**2))
    # STG analytic: delta_phi = 2*pi*(1/omega - 1)
    omega_stg = np.sqrt(1 + 3*M/(r0 + M))
    p_stg_th = 2 * np.pi * (1.0/omega_stg - 1)
    print(f"  {r0_ratio:8d} {p_stg:12.6f} {p_stg_th:12.6f} {p_gr:12.6f} {p_gr_th:12.6f}")

print(f"""
  REZULTAT CRITIC:
  ================
  GR:  precesie PROGRADA  dphi = +6*pi*M/(a*(1-e^2))  [Mercury: +43"/century = OK]
  STG: precesie RETROGRADA dphi ~ -3*pi*M/r0           [ar prezice ~-21.5"/century]

  Cauza: forta STG F = Mr/(r+M)^3 scade MAI REPEDE decat 1/r^2 la distante mici.
  Fortele care scad mai repede => precesie retrograda (Bertrand's theorem).
  GR are un termen EXTRA (-3ML^2/r^4) care face forta EFECTIVA mai puternica
  la distante mici => precesie prograda.

  CATEGORIE: DERIVAT (din ecuatia orbitei STG)
  STATUS: PROBLEMA DESCHISA - STG da semnul gresit al precesiei.
""")


# ==============================================================
# 6. VITEZA ORBITALA CIRCULARA
# ==============================================================

print("=" * 70)
print("6. VITEZA ORBITALA (curbe de rotatie)")
print("=" * 70)

r = np.linspace(3*M, 100*M, 200)
v_stg = r * np.sqrt(M) / (r + M)**1.5
v_gr = np.sqrt(M / r)  # Newtonian (GR circular orbit in Schwarzschild coords is more complex)

# GR exact: v^2 = M/(r-2M) * r ... actually v = sqrt(M/r) / sqrt(1-3M/r) for Schwarzschild
v_gr_exact = np.where(r > 3*M, np.sqrt(M/r) / np.sqrt(1 - 3*M/r), np.nan)

print(f"\n  {'r/M':>8} {'v_STG':>10} {'v_Newton':>10} {'v_GR':>10}")
for rr in [3, 5, 10, 20, 50, 100]:
    vs = rr*M * np.sqrt(M) / (rr*M + M)**1.5
    vn = np.sqrt(M / (rr*M))
    vg = np.sqrt(M/(rr*M)) / np.sqrt(1 - 3*M/(rr*M)) if rr > 3 else np.nan
    print(f"  {rr:8d} {vs:10.6f} {vn:10.6f} {vg:10.6f}")

print(f"""
  La r >> M: v_STG = v_Newton = sqrt(M/r)
  STG: v = r*sqrt(M)/(r+M)^(3/2) -> 0 la r->0 (viteza SCADE spre centru!)
  GR: v -> c la r -> 3M (diverge la sfera fotonilor)
  Newton: v -> infinit la r -> 0

  Curba de rotatie STG: SCADE la distante mici (efect unic!)
""")


# ==============================================================
# 7. SUMAR COMPARATIV
# ==============================================================

print("=" * 70)
print("7. SUMAR: STG vs GR vs NEWTON")
print("=" * 70)

print(f"""
  {'Proprietate':<30} {'Newton':>15} {'GR (Schwarzschild)':>20} {'STG':>20}
  {'-'*30} {'-'*15} {'-'*20} {'-'*20}
  {'Forta (forma)':<30} {'M/r^2':>15} {'Complex':>20} {'Mr/(r+M)^3':>20}
  {'Forta la r->0':<30} {'Infinit':>15} {'N/A (orizont)':>20} {'0 (finita!)':>20}
  {'Forta maxima':<30} {'N/A':>15} {'N/A':>20} {'4/(27M) la M/2':>20}
  {'Orizont':<30} {'N/A':>15} {'r = 2M':>20} {'NU EXISTA':>20}
  {'Singularitate':<30} {'Da (r=0)':>15} {'Da (r=0)':>20} {'Nu (F->0)':>20}
  {'Sfera fotonilor':<30} {'N/A':>15} {'r = 3M':>20} {'r = M':>20}
  {'ISCO':<30} {'N/A':>15} {'r = 6M':>20} {'r ~ {isco_stg:.1f}M':>20}
  {'Kepler III':<30} {'T^2 ~ r^3':>15} {'T^2 ~ r^3 (+corr)':>20} {'T^2 ~ (r+M)^3':>20}
  {'Redshift (r=2M)':<30} {'0.5':>15} {'Infinit':>20} {'0.5':>20}
  {'Deflexie lumina':<30} {'2M/b':>15} {'4M/b':>20} {'4M/b':>20}
  {'Potential minim':<30} {'-Infinit':>15} {'N/A':>20} {'-1/2':>20}
  {'Limita camp slab':<30} {'=':>15} {'= Newton':>20} {'= Newton':>20}
""")

print(f"""
  DIFERENTE CHEIE STG vs GR:
  ==========================
  1. STG nu are orizont de evenimente => "gauri negre" fara orizont
  2. STG: forta -> 0 la r -> 0 => fara singularitate fizica
  3. STG: sfera fotonilor la r = M (nu 3M)
  4. STG: Kepler modificat r -> r+M (GR are o corectie diferita)
  5. STG: potential finit V_min = -1/2 (GR diverge)

  IDENTITATI la camp slab (r >> M):
  ==================================
  - Forta: M/r^2 (toate trei)
  - Redshift: M/r (toate trei)
  - Deflexie lumina: 4M/b (STG = GR, Newton da 2M/b)
  - Kepler: T^2 ~ r^3 (toate trei)
""")


# ==============================================================
# PLOTS
# ==============================================================

fig, axes = plt.subplots(2, 3, figsize=(18, 10))
fig.suptitle('EXP-103: STG vs GR - All Predictions\n'
             r'STG: $ds^2 = -dt^2/C^2 + C^2 dr^2$, $C = 1+M/r$ | '
             r'GR: $ds^2 = -(1-2M/r)dt^2 + \frac{dr^2}{1-2M/r}$',
             fontsize=11)

# 1. Force comparison
ax = axes[0, 0]
r = np.linspace(0.01, 15, 1000)
F_stg = M * r / (r + M)**3
F_newton = M / r**2
F_gr_approx = M / r**2 * np.where(r > 2*M, 1, np.nan)  # approximate
ax.plot(r/M, F_stg, 'r-', linewidth=2, label='STG: Mr/(r+M)^3')
ax.plot(r/M, F_newton, 'b--', linewidth=1.5, label='Newton: M/r^2')
ax.axvline(x=2, color='gray', linestyle=':', alpha=0.5, label='r=2M (GR horizon)')
ax.axvline(x=0.5, color='green', linestyle=':', alpha=0.5, label='r=M/2 (STG F_max)')
ax.set_xlabel('r / M')
ax.set_ylabel('Force')
ax.set_title('Gravitational Force')
ax.legend(fontsize=7)
ax.grid(True, alpha=0.3)
ax.set_ylim(0, 0.5)
ax.set_xlim(0, 15)

# 2. Effective potential
ax = axes[0, 1]
r = np.linspace(0.5, 30, 1000)
for L_val in [3.0, 4.0, 5.0]:
    V_s = np.array([V_eff_stg(ri, L_val) for ri in r])
    V_g = np.array([V_eff_gr(ri, L_val) if ri > 2*M else np.nan for ri in r])
    ax.plot(r/M, V_s, 'r-', linewidth=1.5, alpha=0.7)
    ax.plot(r/M, V_g, 'b--', linewidth=1.5, alpha=0.7)
ax.plot([], [], 'r-', label='STG')
ax.plot([], [], 'b--', label='GR')
ax.set_xlabel('r / M')
ax.set_ylabel('V_eff')
ax.set_title(f'Effective Potential (L = 3, 4, 5)')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)
ax.set_ylim(0.3, 0.7)

# 3. Redshift
ax = axes[0, 2]
r = np.linspace(0.5, 20, 500)
z_stg = np.array([redshift_stg(ri) for ri in r])
z_gr = np.array([redshift_gr(ri) if ri > 2*M else np.nan for ri in r])
z_newton = M / r
ax.semilogy(r/M, 1+z_stg, 'r-', linewidth=2, label='STG: 1+z = 1+M/r')
ax.semilogy(r/M, 1+z_gr, 'b--', linewidth=2, label='GR: 1+z = 1/sqrt(1-2M/r)')
ax.semilogy(r/M, 1+z_newton, 'g:', linewidth=1.5, label='Newton: 1+z = 1+M/r')
ax.axvline(x=2, color='blue', linestyle=':', alpha=0.3)
ax.set_xlabel('r / M')
ax.set_ylabel('1 + z')
ax.set_title('Gravitational Redshift')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

# 4. Light deflection
ax = axes[1, 0]
b_vals = np.logspace(0.7, 4, 50)
d_stg_arr = []
d_gr_arr = []
for b in b_vals:
    try:
        d_s = light_deflection_stg(b*M, M)
        d_g = light_deflection_gr(b*M, M)
    except:
        d_s = np.nan
        d_g = np.nan
    d_stg_arr.append(d_s)
    d_gr_arr.append(d_g)
d_stg_arr = np.array(d_stg_arr)
d_gr_arr = np.array(d_gr_arr)
d_weak = 4*M / (b_vals * M)

ax.loglog(b_vals, np.abs(d_stg_arr), 'r-', linewidth=2, label='STG')
ax.loglog(b_vals, np.abs(d_gr_arr), 'b--', linewidth=2, label='GR')
ax.loglog(b_vals, d_weak, 'g:', linewidth=1.5, label='4M/b (weak field)')
ax.set_xlabel('b / M')
ax.set_ylabel('Deflection angle (rad)')
ax.set_title('Light Deflection')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

# 5. Orbital velocity
ax = axes[1, 1]
r = np.linspace(1.5, 50, 500)
v_stg_arr = r * np.sqrt(M) / (r + M)**1.5
v_newton_arr = np.sqrt(M / r)
v_gr_arr = np.where(r > 3*M, np.sqrt(M/r) / np.sqrt(1 - 3*M/r), np.nan)
ax.plot(r/M, v_stg_arr, 'r-', linewidth=2, label='STG')
ax.plot(r/M, v_newton_arr, 'g:', linewidth=1.5, label='Newton')
ax.plot(r/M, v_gr_arr, 'b--', linewidth=2, label='GR')
ax.axvline(x=1, color='red', linestyle=':', alpha=0.3, label='STG photon sphere')
ax.axvline(x=3, color='blue', linestyle=':', alpha=0.3, label='GR photon sphere')
ax.set_xlabel('r / M')
ax.set_ylabel('v_circular')
ax.set_title('Circular Orbital Velocity')
ax.legend(fontsize=7)
ax.grid(True, alpha=0.3)
ax.set_ylim(0, 1.0)

# 6. No horizon diagram
ax = axes[1, 2]
r = np.linspace(0.1, 10, 500)
gtt_stg = -r**2 / (r + M)**2
gtt_gr = -(1 - 2*M/r)
ax.plot(r/M, gtt_stg, 'r-', linewidth=2, label=r'STG: $g_{tt} = -r^2/(r+M)^2$')
ax.plot(r/M, gtt_gr, 'b--', linewidth=2, label=r'GR: $g_{tt} = -(1-2M/r)$')
ax.axhline(y=0, color='black', linewidth=0.5)
ax.fill_between(r/M, gtt_gr, 0, where=(gtt_gr > 0), color='blue', alpha=0.1, label='GR: g_tt > 0 (inside horizon)')
ax.set_xlabel('r / M')
ax.set_ylabel(r'$g_{tt}$')
ax.set_title('Metric component g_tt (horizon = where g_tt = 0)')
ax.legend(fontsize=7)
ax.grid(True, alpha=0.3)
ax.set_xlim(0, 10)
ax.set_ylim(-1.2, 0.5)
ax.annotate('GR horizon\nr = 2M', xy=(2, 0), fontsize=8, color='blue',
            ha='center', va='bottom')
ax.annotate('STG: NO horizon\ng_tt < 0 always', xy=(1, -0.25), fontsize=8, color='red',
            ha='center')

plt.tight_layout()
plt.savefig('exp103_stg_predictions.png', dpi=150, bbox_inches='tight')
print(f"\nSalvat: exp103_stg_predictions.png")

plt.show()
