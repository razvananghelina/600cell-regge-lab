"""
EXP-104: Precesia orbitala in 3 metrici STG vs GR
==================================================

Intrebare: STG defineste C(r) = 1 + M/r, dar metrica COMPLETA poate fi:
  A) Conforma:  ds^2 = -dt^2/C^2 + C^2(dr^2 + r^2 dOmega^2)
  B) Mixta:     ds^2 = -dt^2/C^2 + C^2 dr^2 + r^2 dOmega^2
  GR: Schwarzschild ds^2 = -(1-2M/r)dt^2 + dr^2/(1-2M/r) + r^2 dOmega^2

Ecuatiile orbitei (Binet) pentru particula masiva:
  A) d^2u/dphi^2 + u = M / (L_A^2 * (1+Mu)^3)
     L_A^2 = M*(r0+M)^2 / (r0 - M)  [circular orbit]
     omega^2 = 1 + 3M/(r0+M) > 1  =>  RETROGRAD

  B) d^2u/dphi^2 + u = u + (M - L_B^2 * u) / (L_B^2 * (1+Mu)^3)
     L_B^2 = M*r0  [circular orbit, = Newtonian!]
     omega^2 = 1/(1+M/r0)^3 < 1  =>  PROGRAD

  GR) d^2u/dphi^2 + u = M/L_GR^2 + 3*M*u^2
     L_GR^2 = M*r0^2 / (r0 - 3M)  [circular orbit]
     omega^2 = 1 - 6M/r0 < 1  =>  PROGRAD

Predictii analitice (camp slab, r0 >> M):
  A) dphi = -3*pi*M/r0  (retrograd)
  B) dphi = +3*pi*M/r0  (prograd)
  GR) dphi = +6*pi*M/r0  (prograd)

DERIVAT: ecuatia orbitei din fiecare metrica prin calculul geodezicelor.
"""

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

M = 1.0

# ==============================================================
# DERIVATIONS (documented)
# ==============================================================

print("=" * 70)
print("EXP-104: PRECESIA ORBITALA - 3 METRICI")
print("=" * 70)

print("""
  METRICA A (conforma): ds^2 = -dt^2/C^2 + C^2(dr^2 + r^2 dOmega^2)
  METRICA B (mixta):    ds^2 = -dt^2/C^2 + C^2 dr^2 + r^2 dOmega^2
  GR (Schwarzschild):   ds^2 = -(1-2M/r)dt^2 + dr^2/(1-2M/r) + r^2 dOmega^2

  Diferenta: partea ANGULARA -- C^2*r^2 dOmega^2 (A) vs r^2 dOmega^2 (B)
  Aceasta schimba definitia momentului cinetic si ecuatia orbitei.
""")

# ==============================================================
# METRIC A: Conformal STG
# ==============================================================
# L_A = C^2 * r^2 * dphi/dtau
# Circular orbit: L_A^2 = M*(r0+M)^2/(r0-M), exists for r0 > M
# Orbit eq: d^2u/dphi^2 + u = M/(L_A^2*(1+Mu)^3)

def L2_circular_A(r0):
    """Squared angular momentum for circular orbit, metric A."""
    if r0 <= M:
        return None
    return M * (r0 + M)**2 / (r0 - M)

def binet_rhs_A(u, L2):
    """RHS of d^2u/dphi^2 + u = RHS for metric A."""
    return M / (L2 * (1 + M*u)**3)

def omega2_A(r0):
    """Radial oscillation frequency squared, metric A."""
    return 1 + 3*M/(r0 + M)

# ==============================================================
# METRIC B: Mixed STG
# ==============================================================
# L_B = r^2 * dphi/dtau (standard)
# Circular orbit: L_B^2 = M*r0 (Newtonian!)
# Orbit eq: d^2u/dphi^2 = (M - L^2*u) / (L^2*(1+Mu)^3)
# Rewritten: d^2u/dphi^2 + u = u + (M - L^2*u)/(L^2*(1+Mu)^3) = F_B(u)

def L2_circular_B(r0):
    """Squared angular momentum for circular orbit, metric B."""
    return M * r0

def binet_rhs_B(u, L2):
    """RHS of d^2u/dphi^2 + u = RHS for metric B."""
    return u + (M - L2*u) / (L2 * (1 + M*u)**3)

def omega2_B(r0):
    """Radial oscillation frequency squared, metric B."""
    return 1.0 / (1 + M/r0)**3

# ==============================================================
# GR: Schwarzschild
# ==============================================================
# L_GR = r^2 * dphi/dtau
# Circular orbit: L_GR^2 = M*r0^2/(r0-3M), exists for r0 > 3M
# Orbit eq: d^2u/dphi^2 + u = M/L^2 + 3*M*u^2

def L2_circular_GR(r0):
    """Squared angular momentum for circular orbit, GR."""
    if r0 <= 3*M:
        return None
    return M * r0**2 / (r0 - 3*M)

def binet_rhs_GR(u, L2):
    """RHS of d^2u/dphi^2 + u = RHS for GR."""
    return M / L2 + 3*M*u**2

def omega2_GR(r0):
    """Radial oscillation frequency squared, GR."""
    return 1 - 6*M/r0


# ==============================================================
# NUMERICAL ORBIT SOLVER
# ==============================================================

def compute_precession_ode(r0, e, metric='B', n_orbits=5):
    """Compute precession by solving the Binet ODE.

    Solve: d^2u/dphi^2 + u = RHS(u, L^2)
    as IVP: u(0) = u_peri, u'(0) = 0.
    Track perihelion positions (u' crosses zero from below).
    Precession = mean angle between perihelions - 2*pi.
    """
    # Circular orbit angular momentum
    if metric == 'A':
        L2 = L2_circular_A(r0)
        rhs_func = lambda u: binet_rhs_A(u, L2)
    elif metric == 'B':
        L2 = L2_circular_B(r0)
        rhs_func = lambda u: binet_rhs_B(u, L2)
    else:  # GR
        L2 = L2_circular_GR(r0)
        rhs_func = lambda u: binet_rhs_GR(u, L2)

    if L2 is None:
        return np.nan, 0

    u0 = 1.0 / r0

    # For eccentric orbit starting at perihelion:
    # u_peri > u0 (closer to mass). For Keplerian: u_peri = u0*(1+e)
    # This is approximate but good enough for small e
    u_peri = u0 * (1 + e)

    # ODE: y = [u, du/dphi]
    def ode(phi, y):
        u, up = y
        if u < 0:
            u = 1e-10
        duup = rhs_func(u) - u  # d^2u/dphi^2 = RHS - u
        return [up, duup]

    # Detect perihelion: du/dphi crosses zero from + to -
    # (u reaches maximum = closest approach)
    def event_perihelion(phi, y):
        return y[1]  # du/dphi = 0
    event_perihelion.direction = -1  # + to -

    # Integrate for n_orbits * 2*pi
    phi_max = n_orbits * 2 * np.pi * 1.5  # generous upper bound
    sol = solve_ivp(ode, [0, phi_max], [u_peri, 0.0],
                    events=event_perihelion,
                    rtol=1e-12, atol=1e-14, max_step=0.05)

    perihelion_phis = sol.t_events[0]
    n_peri = len(perihelion_phis)

    if n_peri >= 3:
        # Angle between consecutive perihelions
        d_phis = np.diff(perihelion_phis)
        # Skip first interval (might be affected by initial conditions)
        if len(d_phis) > 2:
            d_phis = d_phis[1:]
        # Precession = excess over 2*pi per orbit
        precession = np.mean(d_phis) - 2 * np.pi
        return precession, n_peri
    return np.nan, n_peri


# ==============================================================
# ANALYTIC PRECESSION FORMULAS
# ==============================================================

def precession_analytic_A(r0):
    """Analytic precession for metric A."""
    omega = np.sqrt(omega2_A(r0))
    return 2 * np.pi * (1.0/omega - 1)

def precession_analytic_B(r0):
    """Analytic precession for metric B."""
    omega = np.sqrt(omega2_B(r0))
    return 2 * np.pi * (1.0/omega - 1)

def precession_analytic_GR(r0, e=0):
    """Analytic precession for GR (exact at e=0)."""
    omega = np.sqrt(omega2_GR(r0))
    if omega <= 0:
        return np.inf
    return 2 * np.pi * (1.0/omega - 1)

def precession_GR_standard(r0, e=0):
    """Standard GR formula: 6*pi*M / (r0*(1-e^2))."""
    return 6 * np.pi * M / (r0 * (1 - e**2))


# ==============================================================
# DERIVATIONS (print)
# ==============================================================

print("=" * 70)
print("1. ECUATIILE ORBITEI (Binet) - DERIVARE")
print("=" * 70)

print("""
  METRICA A (conforma):
  ---------------------
  Metrica: g_tt = -1/C^2, g_rr = C^2, g_theta = C^2*r^2
  Moment cinetic: L_A = C^2 * r^2 * dphi/dtau
  Normalizare: -E^2*C^2 + C^2*rdot^2 + L_A^2/(C^2*r^2) = -1

  Ecuatia orbitei:
    d^2u/dphi^2 + u = M / (L_A^2 * (1+Mu)^3)

  Orbita circulara: L_A^2 = M*(r0+M)^2 / (r0 - M)
  Frecventa: omega^2 = 1 + 3M/(r0+M) > 1  =>  RETROGRAD


  METRICA B (mixta):
  ------------------
  Metrica: g_tt = -1/C^2, g_rr = C^2, g_theta = r^2
  Moment cinetic: L_B = r^2 * dphi/dtau  (standard)
  Normalizare: -E^2*C^2 + C^2*rdot^2 + L_B^2/r^2 = -1

  Ecuatia orbitei:
    d^2u/dphi^2 = (M - L_B^2 * u) / (L_B^2 * (1+Mu)^3)

  Orbita circulara: L_B^2 = M*r0  (identic Newton!)
  Frecventa: omega^2 = 1/(1+M/r0)^3 < 1  =>  PROGRAD


  GR (Schwarzschild):
  -------------------
  Ecuatia orbitei: d^2u/dphi^2 + u = M/L^2 + 3Mu^2
  Orbita circulara: L^2 = M*r0^2/(r0-3M)
  Frecventa: omega^2 = 1 - 6M/r0 < 1  =>  PROGRAD
""")


# ==============================================================
# 2. NUMERIC PRECESSION TEST
# ==============================================================

print("=" * 70)
print("2. PRECESIA NUMERICA vs ANALITICA")
print("=" * 70)

e_test = 0.1  # Small eccentricity for clean comparison

print(f"\n  Eccentricitate e = {e_test}")
print(f"  {'r0/M':>6} | {'A_num':>10} {'A_anal':>10} | {'B_num':>10} {'B_anal':>10} | {'GR_num':>10} {'GR_anal':>10} {'GR_std':>10}")
print(f"  {'-'*6}-+-{'-'*10}-{'-'*10}-+-{'-'*10}-{'-'*10}-+-{'-'*10}-{'-'*10}-{'-'*10}")

results = []
for r0_ratio in [20, 50, 100, 200, 500, 1000]:
    r0 = r0_ratio * M

    # Metric A
    pA, nA = compute_precession_ode(r0, e_test, 'A', n_orbits=8)
    pA_th = precession_analytic_A(r0)

    # Metric B
    pB, nB = compute_precession_ode(r0, e_test, 'B', n_orbits=8)
    pB_th = precession_analytic_B(r0)

    # GR
    pGR, nGR = compute_precession_ode(r0, e_test, 'GR', n_orbits=8)
    pGR_th = precession_analytic_GR(r0)
    pGR_std = precession_GR_standard(r0, e_test)

    results.append((r0_ratio, pA, pA_th, pB, pB_th, pGR, pGR_th, pGR_std))

    print(f"  {r0_ratio:6d} | {pA:10.6f} {pA_th:10.6f} | {pB:10.6f} {pB_th:10.6f} | {pGR:10.6f} {pGR_th:10.6f} {pGR_std:10.6f}")


# ==============================================================
# 3. RATIOS
# ==============================================================

print(f"\n\n  Rapoarte precesie (camp slab, r >> M):")
print(f"  {'r0/M':>6} | {'A/GR':>10} {'B/GR':>10} | {'A_th/GR_th':>12} {'B_th/GR_th':>12}")
print(f"  {'-'*6}-+-{'-'*10}-{'-'*10}-+-{'-'*12}-{'-'*12}")

for r0_ratio, pA, pA_th, pB, pB_th, pGR, pGR_th, pGR_std in results:
    rAG = pA / pGR if abs(pGR) > 1e-15 else np.nan
    rBG = pB / pGR if abs(pGR) > 1e-15 else np.nan
    rAG_th = pA_th / pGR_th if abs(pGR_th) > 1e-10 else np.nan
    rBG_th = pB_th / pGR_th if abs(pGR_th) > 1e-10 else np.nan
    print(f"  {r0_ratio:6d} | {rAG:10.4f} {rBG:10.4f} | {rAG_th:12.6f} {rBG_th:12.6f}")


# ==============================================================
# 4. LIGHT DEFLECTION FOR ALL 3 METRICS
# ==============================================================

print(f"\n\n{'='*70}")
print("3. DEFLEXA LUMINII - TOATE METRICILE")
print("=" * 70)

print("""
  Ecuatiile orbitei pentru LUMINA (geodezice nule):

  Metrica A: d^2u/dphi^2 + u = 2M(1+Mu)^3/b^2
  Metrica B: d^2u/dphi^2 + u/(1+Mu)^3 = 0
             (rewritten: d^2u/dphi^2 + u = u*(1 - 1/(1+Mu)^3))
  GR:        d^2u/dphi^2 + u = 3Mu^2

  La camp slab (1PN): toate dau 4M/b (identic!)
  La camp puternic: difera.
""")

def light_deflection(b, metric='B'):
    """Compute light deflection by solving Binet ODE."""
    from scipy.integrate import solve_ivp

    if metric == 'A':
        def rhs(phi, y):
            u, up = y
            return [up, 2*M*(1 + M*u)**3 / b**2 - u]
    elif metric == 'B':
        def rhs(phi, y):
            u, up = y
            # d^2u/dphi^2 = -u/(1+Mu)^3
            # rewritten: d^2u/dphi^2 = u*(1 - 1/(1+Mu)^3) - u = -u/(1+Mu)^3
            return [up, -u / (1 + M*u)**3]
    else:  # GR
        def rhs(phi, y):
            u, up = y
            return [up, 3*M*u**2 - u]

    y0 = [0.0, 1.0/b]

    def u_zero(phi, y):
        return y[0]
    u_zero.terminal = True
    u_zero.direction = -1

    sol = solve_ivp(rhs, [0, 2*np.pi], y0, events=u_zero,
                    rtol=1e-12, atol=1e-14, max_step=0.01)

    if sol.t_events[0].size > 0:
        return sol.t_events[0][0] - np.pi
    return np.inf


print(f"  {'b/M':>8} | {'d_A(rad)':>12} {'d_B(rad)':>12} {'d_GR(rad)':>12} | {'A/GR':>8} {'B/GR':>8} {'4M/b':>12}")
print(f"  {'-'*8}-+-{'-'*12}-{'-'*12}-{'-'*12}-+-{'-'*8}-{'-'*8}-{'-'*12}")

for b_ratio in [10, 20, 50, 100, 500, 1000, 10000]:
    b = b_ratio * M
    dA = light_deflection(b, 'A')
    dB = light_deflection(b, 'B')
    dGR = light_deflection(b, 'GR')
    d_weak = 4*M/b
    rAG = dA/dGR if 0 < dGR < 10 else np.nan
    rBG = dB/dGR if 0 < dGR < 10 else np.nan
    print(f"  {b_ratio:8d} | {dA:12.6e} {dB:12.6e} {dGR:12.6e} | {rAG:8.4f} {rBG:8.4f} {d_weak:12.6e}")


# ==============================================================
# 5. SUMMARY
# ==============================================================

print(f"\n\n{'='*70}")
print("4. SUMAR COMPARATIV")
print("=" * 70)

print("""
  Proprietate            Metrica A (conf)    Metrica B (mixta)     GR (Schwarzschild)
  -------------------    ----------------    -----------------     ------------------
  Moment cinetic L       C^2*r^2*dphi/dtau   r^2*dphi/dtau        r^2*dphi/dtau
  L_circ^2               M(r0+M)^2/(r0-M)   M*r0 (Newtonian!)    M*r0^2/(r0-3M)
  Precesie (semn)        RETROGRAD (-)       PROGRAD (+)           PROGRAD (+)
  Precesie (ord.1)       -3*pi*M/r0          +3*pi*M/r0           +6*pi*M/r0
  Precesie/GR            -1/2                +1/2                  1
  Deflexie lumina        4M/b                4M/b                  4M/b
  Mercury (43"/cent)     -21.5"/cent !!      +21.5"/cent           +43"/cent

  CONCLUZIE:
  ==========
  - Metrica A: EXCLUSA (precesie retrograd, contrazice Mercury)
  - Metrica B: semn CORECT dar magnitudine 1/2 din GR
  - GR: corect (confirmat observational)

  IMPLICATIE PENTRU TEORIA STG:
  =============================
  Daca STG e corecta, metrica trebuie sa fie de tip B (mixta),
  dar cu o CORECTIE suplimentara care sa dubleze precesia.

  Posibilitati:
  1. Termenul de ordin 2 in C poate modifica precesia
  2. Topologia discreta (600-cell) poate adauga un termen
  3. O metrica intermediara (C^alpha * r^2 dOmega^2, alpha != 0,1)
""")


# ==============================================================
# 6. SEARCH FOR OPTIMAL ALPHA
# ==============================================================

print("=" * 70)
print("5. CAUTARE: CE ALPHA DA PRECESIA GR?")
print("=" * 70)

print("""
  Metrica generalizata: ds^2 = -dt^2/C^2 + C^2 dr^2 + C^{2*alpha} r^2 dOmega^2
  alpha=0 => Metrica B (precesie +3*pi*M/r0)
  alpha=1 => Metrica A (precesie -3*pi*M/r0)

  Cautam alpha care da precesie = +6*pi*M/r0 (GR)
""")

# For general alpha, the angular momentum is:
# L = C^alpha * r^2 * dphi/dtau
# The orbit equation changes. Let me derive it for general alpha.
#
# Normalization: -E^2*C^2 + C^2*rdot^2 + L^2/(C^{2alpha}*r^2) = -1
# rdot^2 = E^2 - 1/C^2 - L^2/(C^{2+2alpha}*r^2)
#
# Using u = 1/r, dphi/dtau = L*u^2/C^{2alpha}:
# L^2(du/dphi)^2 = [E^2 - (1 + L^2*u^2/C^{2alpha})/C^2] * C^{4alpha}
#
# This gets complicated. Let me just do it numerically:
# solve the orbit ODE for different alpha values and measure precession.

def compute_precession_alpha(r0, e, alpha, n_orbits=8):
    """Compute precession for metric with C^{2*alpha} angular part."""
    # Circular orbit condition derived from dV_eff/dr = 0
    # For general alpha, this needs to be solved numerically

    # V_eff(r) = [1 + L^2*u^2/C^{2alpha}] / (2C^2)
    # where C = 1+M/r and u = 1/r

    # For circular orbit at r0, compute L^2 numerically:
    u0 = 1.0/r0
    C0 = 1 + M*u0

    # Condition: dV_eff/dr = 0 at r0
    # Compute numerically
    def V_eff(r, L2):
        u = 1.0/r
        C = 1 + M*u
        return (1 + L2*u**2/C**(2*alpha)) / (2*C**2)

    # dV_eff/dr = 0 => solve for L^2
    dr = 1e-6
    # dV/dr = (V(r+dr) - V(r-dr))/(2dr) = 0 for circular
    # dV/dr = dV/dL2 = 0 is wrong. We need:
    # dV/dr(L2) = 0 at r=r0

    # V_eff = [C^{-2} + L^2*u^2*C^{-2-2alpha}] / 2
    # dV/dr = d/dr[C^{-2}]/2 + L^2/2 * d/dr[u^2*C^{-2-2alpha}]
    # = Mr/(r+M)^3 + L^2/2 * d/dr[u^2*C^{-2-2alpha}]

    # Let h(r) = u^2 * C^{-2-2alpha} = r^{-2} * (1+M/r)^{-2-2alpha}
    #          = r^{-2} * (r/(r+M))^{2+2alpha} = r^{2alpha} / (r+M)^{2+2alpha}
    # dh/dr = [2alpha*r^{2alpha-1}*(r+M)^{2+2alpha} - r^{2alpha}*(2+2alpha)*(r+M)^{1+2alpha}]
    #         / (r+M)^{2*(2+2alpha)}
    # = r^{2alpha-1} * [2alpha*(r+M) - (2+2alpha)*r] / (r+M)^{3+2alpha}
    # = r^{2alpha-1} * [2alpha*r + 2alpha*M - 2r - 2alpha*r] / (r+M)^{3+2alpha}
    # = r^{2alpha-1} * [2alpha*M - 2r] / (r+M)^{3+2alpha}
    # = 2*r^{2alpha-1} * (alpha*M - r) / (r+M)^{3+2alpha}

    # dV/dr = Mr/(r+M)^3 + L^2 * r^{2alpha-1} * (alpha*M - r) / (r+M)^{3+2alpha} = 0
    # L^2 = -Mr * (r+M)^{3+2alpha} / ((r+M)^3 * r^{2alpha-1} * (alpha*M - r))
    # L^2 = -M * r^{2-2alpha} * (r+M)^{2alpha} / (alpha*M - r)
    # L^2 = M * r^{2-2alpha} * (r+M)^{2alpha} / (r - alpha*M)

    if r0 <= alpha * M:
        return np.nan, 0

    L2 = M * r0**(2-2*alpha) * (r0+M)**(2*alpha) / (r0 - alpha*M)

    # Orbit ODE: derive d^2u/dphi^2 from the energy equation
    # This is complex for general alpha. Use numerical differentiation.

    # Energy at circular orbit:
    E2 = (1 + L2*u0**2/C0**(2*alpha)) / C0**2

    def energy_eq_rhs(u):
        """RHS of L^2*(du/dphi)^2 = ..."""
        C = 1 + M*u
        return E2 - (1 + L2*u**2/C**(2*alpha)) / C**2

    # d^2u/dphi^2 = (1/2) * d/du[energy_eq_rhs] / L^2 ... no
    # From L^2(du/dphi)^2 = F(u), differentiating:
    # 2L^2(du/dphi)(d^2u/dphi^2) = F'(u)(du/dphi)
    # d^2u/dphi^2 = F'(u) / (2L^2)

    def d2u_dphi2(u):
        du = 1e-8 * max(abs(u), 1e-6)
        F_plus = energy_eq_rhs(u + du)
        F_minus = energy_eq_rhs(u - du)
        dF_du = (F_plus - F_minus) / (2*du)
        return dF_du / (2*L2)

    # Perihelion: u_peri = u0*(1+e)
    u_peri = u0 * (1 + e)

    def ode(phi, y):
        u, up = y
        if u < 1e-12:
            u = 1e-12
        return [up, d2u_dphi2(u)]

    def event_peri(phi, y):
        return y[1]
    event_peri.direction = -1

    phi_max = n_orbits * 2 * np.pi * 1.5
    sol = solve_ivp(ode, [0, phi_max], [u_peri, 0.0],
                    events=event_peri,
                    rtol=1e-11, atol=1e-13, max_step=0.05)

    perihelion_phis = sol.t_events[0]
    n_peri = len(perihelion_phis)

    if n_peri >= 3:
        d_phis = np.diff(perihelion_phis)
        if len(d_phis) > 2:
            d_phis = d_phis[1:]
        precession = np.mean(d_phis) - 2*np.pi
        return precession, n_peri
    return np.nan, n_peri


# Test for different alpha values at r0 = 100M
r0_test = 100 * M
e_test_alpha = 0.1
GR_target = precession_GR_standard(r0_test, e_test_alpha)

print(f"\n  Test la r0 = {r0_test:.0f}M, e = {e_test_alpha}")
print(f"  GR target: {GR_target:.6f} rad/orbit")
print(f"\n  {'alpha':>8} {'precesie':>12} {'prec/GR':>10}")
print(f"  {'-'*8} {'-'*12} {'-'*10}")

alpha_values = np.linspace(-0.5, 1.5, 21)
prec_vs_alpha = []

for alpha in alpha_values:
    p, n = compute_precession_alpha(r0_test, e_test_alpha, alpha, n_orbits=8)
    ratio = p / GR_target if abs(GR_target) > 1e-15 and not np.isnan(p) else np.nan
    prec_vs_alpha.append((alpha, p, ratio))
    print(f"  {alpha:8.3f} {p:12.6f} {ratio:10.4f}")

# Find alpha where precession matches GR
print(f"\n  Cautam alpha unde precesie/GR = 1.0...")
from scipy.optimize import brentq

def prec_ratio_minus_1(alpha):
    p, n = compute_precession_alpha(r0_test, e_test_alpha, alpha, n_orbits=8)
    if np.isnan(p):
        return 10.0  # far from target
    return p / GR_target - 1.0

# Search in region where ratio crosses 1
# From the scan, find approximate location
found_alpha = None
for i in range(len(prec_vs_alpha) - 1):
    a1, p1, r1 = prec_vs_alpha[i]
    a2, p2, r2 = prec_vs_alpha[i+1]
    if not np.isnan(r1) and not np.isnan(r2) and (r1 - 1.0) * (r2 - 1.0) < 0:
        try:
            alpha_opt = brentq(prec_ratio_minus_1, a1, a2, xtol=1e-4)
            found_alpha = alpha_opt
        except:
            pass

if found_alpha is not None:
    p_opt, n_opt = compute_precession_alpha(r0_test, e_test_alpha, found_alpha, n_orbits=8)
    ratio_opt = p_opt / GR_target
    print(f"\n  GASIT: alpha = {found_alpha:.4f}")
    print(f"  Precesie = {p_opt:.6f} rad/orbit (GR target: {GR_target:.6f})")
    print(f"  Ratio = {ratio_opt:.6f}")

    # Verify at other r0 values
    print(f"\n  Verificare la alte r0:")
    print(f"  {'r0/M':>6} {'prec_alpha':>12} {'prec_GR':>12} {'ratio':>10}")
    for r0r in [50, 200, 500]:
        r0v = r0r * M
        pv, nv = compute_precession_alpha(r0v, e_test_alpha, found_alpha, n_orbits=8)
        GRv = precession_GR_standard(r0v, e_test_alpha)
        rv = pv / GRv if abs(GRv) > 1e-15 and not np.isnan(pv) else np.nan
        print(f"  {r0r:6d} {pv:12.6f} {GRv:12.6f} {rv:10.4f}")
else:
    print(f"\n  Nu s-a gasit alpha care da precesia GR exact.")
    print(f"  Cel mai aproape:")
    closest = min(prec_vs_alpha, key=lambda x: abs(x[2]-1.0) if not np.isnan(x[2]) else 999)
    print(f"  alpha = {closest[0]:.3f}, ratio = {closest[2]:.4f}")


# ==============================================================
# PLOTS
# ==============================================================

fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle('EXP-104: Precession in 3 STG Metrics vs GR', fontsize=12)

# Plot 1: Precession vs r0 for all metrics
ax = axes[0, 0]
r0_range = np.linspace(20, 500, 50)
for label, func, color, ls in [
    ('Metric A (retrograd)', precession_analytic_A, 'red', '-'),
    ('Metric B (prograd)', precession_analytic_B, 'blue', '-'),
    ('GR (Schwarzschild)', lambda r: precession_GR_standard(r, 0), 'black', '--'),
]:
    prec = [func(r*M) for r in r0_range]
    ax.plot(r0_range, prec, color=color, linestyle=ls, linewidth=2, label=label)
ax.axhline(y=0, color='gray', linestyle=':', alpha=0.5)
ax.set_xlabel('r0 / M')
ax.set_ylabel('Precession (rad/orbit)')
ax.set_title('Precession vs distance')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

# Plot 2: Precession ratio vs r0
ax = axes[0, 1]
for label, func, color in [
    ('A/GR', precession_analytic_A, 'red'),
    ('B/GR', precession_analytic_B, 'blue'),
]:
    ratio = [func(r*M) / precession_GR_standard(r*M, 0) for r in r0_range]
    ax.plot(r0_range, ratio, color=color, linewidth=2, label=label)
ax.axhline(y=-0.5, color='red', linestyle=':', alpha=0.5, label='A/GR -> -1/2')
ax.axhline(y=0.5, color='blue', linestyle=':', alpha=0.5, label='B/GR -> +1/2')
ax.axhline(y=1.0, color='black', linestyle='--', alpha=0.5, label='GR/GR = 1')
ax.set_xlabel('r0 / M')
ax.set_ylabel('Precession / GR precession')
ax.set_title('Ratio to GR')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

# Plot 3: Precession vs alpha
ax = axes[1, 0]
alphas = [x[0] for x in prec_vs_alpha]
ratios = [x[2] for x in prec_vs_alpha]
ax.plot(alphas, ratios, 'ko-', linewidth=2, markersize=4)
ax.axhline(y=1.0, color='green', linestyle='--', linewidth=2, label='GR target')
ax.axhline(y=0.5, color='blue', linestyle=':', label='Metric B (alpha=0)')
ax.axhline(y=-0.5, color='red', linestyle=':', label='Metric A (alpha=1)')
ax.axhline(y=0, color='gray', linestyle=':')
if found_alpha is not None:
    ax.axvline(x=found_alpha, color='green', linestyle='-', alpha=0.7,
               label=f'alpha={found_alpha:.3f}')
ax.set_xlabel('alpha (exponent on angular metric)')
ax.set_ylabel('Precession / GR precession')
ax.set_title(f'Precession ratio vs alpha (r0={r0_test:.0f}M)')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

# Plot 4: Light deflection comparison
ax = axes[1, 1]
b_range = np.logspace(1, 4, 30)
dA_arr, dB_arr, dGR_arr = [], [], []
for b_ratio in b_range:
    b = b_ratio * M
    dA_arr.append(light_deflection(b, 'A'))
    dB_arr.append(light_deflection(b, 'B'))
    dGR_arr.append(light_deflection(b, 'GR'))
dA_arr = np.array(dA_arr)
dB_arr = np.array(dB_arr)
dGR_arr = np.array(dGR_arr)
d_weak = 4*M/(b_range*M)

ax.loglog(b_range, dA_arr, 'r-', linewidth=2, label='Metric A')
ax.loglog(b_range, dB_arr, 'b-', linewidth=2, label='Metric B')
ax.loglog(b_range, dGR_arr, 'k--', linewidth=2, label='GR')
ax.loglog(b_range, d_weak, 'g:', linewidth=1, label='4M/b (weak field)')
ax.set_xlabel('b / M (impact parameter)')
ax.set_ylabel('Deflection (rad)')
ax.set_title('Light Deflection')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('exp104_precession_3metrics.png', dpi=150, bbox_inches='tight')
print(f"\nSalvat: exp104_precession_3metrics.png")
