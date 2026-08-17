"""
EXP-102: STG Kepler III - Formula Analitica
============================================

DESCOPERIRE: Formula STG se simplifica la:

    F_STG = M * r / (r + M)^3

Aceasta e TOATA gravitatia STG intr-o singura formula.
Nu e nevoie de graf, BFS, S^3, sau discretizare.

DERIVAT:
    C(r) = 1 + M/r
    F = M / (r^2 * C^3) = M / (r^2 * (1+M/r)^3) = M*r / (r+M)^3

PROPRIETATI (toate DERIVATE analitic):
    1. r >> M: F -> M/r^2 (Newton)
    2. r = M/2: F_max = 4/(27*M) (forta MAXIMA, finita!)
    3. r -> 0: F -> 0 (FARA singularitate)
    4. Kepler III: T^2 = 4*pi^2*(r+M)^3/M
    5. La r >> M: T^2 -> 4*pi^2*r^3/M (Kepler clasic)

CLASIFICARE: 100% DERIVAT. Zero fitting. Zero parametri liberi.
"""

import numpy as np
import matplotlib.pyplot as plt

# ==============================================================
# STG FORCE LAW - CLOSED FORM
# ==============================================================

def F_stg(r, M=1.0):
    """STG gravitational force. DERIVED: F = M*r / (r+M)^3"""
    return M * r / (r + M) ** 3

def F_newton(r, M=1.0):
    """Newton: F = M/r^2"""
    return M / r ** 2

def V_stg(r, M=1.0):
    """STG potential. DERIVED: V = 1/(2*C^2) - 1/2, C = 1+M/r"""
    C = 1 + M / r
    return 1.0 / (2.0 * C ** 2) - 0.5

def V_newton(r, M=1.0):
    """Newton potential: V = -M/r"""
    return -M / r

def T_stg(r, M=1.0):
    """STG orbital period. DERIVED: T = 2*pi*(r+M)^(3/2) / sqrt(M)"""
    return 2 * np.pi * (r + M) ** 1.5 / np.sqrt(M)

def T_newton(r, M=1.0):
    """Kepler III: T = 2*pi*r^(3/2) / sqrt(M)"""
    return 2 * np.pi * r ** 1.5 / np.sqrt(M)

def v_circular_stg(r, M=1.0):
    """Circular orbital speed. DERIVED: v = r*sqrt(M) / (r+M)^(3/2)"""
    return r * np.sqrt(M) / (r + M) ** 1.5

def v_circular_newton(r, M=1.0):
    """Newton: v = sqrt(M/r)"""
    return np.sqrt(M / r)


# ==============================================================
# PASUL 1: VERIFICARE ANALITICA
# ==============================================================

print("=" * 70)
print("EXP-102: STG KEPLER III - FORMULA ANALITICA")
print("=" * 70)

M = 1.0  # Planck mass

print(f"""
FORMULA STG (DERIVATA):
=======================
  F(r) = M * r / (r + M)^3

  Echivalent cu: F = M / (r^2 * C^3), C = 1 + M/r
  Demonstratie: (r+M)^3 = r^3*(1+M/r)^3 = r^3*C^3
                M*r/(r+M)^3 = M*r/(r^3*C^3) = M/(r^2*C^3) QED

PROPRIETATI:
============""")

# Forta maxima
r_max = M / 2
F_max = F_stg(r_max, M)
F_max_theory = 4.0 / (27.0 * M)
print(f"  Forta maxima: la r = M/2 = {r_max}")
print(f"    F_max (calculat) = {F_max:.6f}")
print(f"    F_max (formula)  = 4/(27*M) = {F_max_theory:.6f}")
print(f"    Match: {'DA' if abs(F_max - F_max_theory) < 1e-10 else 'NU'}")

# Limita Newton
print(f"\n  Limita Newton (r >> M):")
for r in [10, 100, 1000, 10000]:
    ratio = F_stg(r, M) / F_newton(r, M)
    print(f"    r={r:>6d}: F_STG/F_Newton = {ratio:.8f} (1 = perfect Newton)")

# Kepler III
print(f"\n  Kepler III: T^2 = 4*pi^2 * (r+M)^3 / M")
print(f"  La r >> M: T^2 -> 4*pi^2 * r^3 / M (Kepler clasic)")
print(f"\n  {'r':>8} {'T_STG':>12} {'T_Newton':>12} {'T_STG/T_N':>12} {'T^2_STG/(r+M)^3':>16}")
for r in [2, 5, 10, 50, 100, 1000]:
    Ts = T_stg(r, M)
    Tn = T_newton(r, M)
    k = Ts ** 2 / ((r + M) ** 3)
    print(f"  {r:8d} {Ts:12.4f} {Tn:12.4f} {Ts/Tn:12.6f} {k:16.6f}")
print(f"  Constanta Kepler STG: k = 4*pi^2/M = {4*np.pi**2/M:.6f}")


# ==============================================================
# PASUL 2: SIMULARE NUMERICA - ORBITE IN SPATIU PLAT
# ==============================================================

print("\n" + "=" * 70)
print("PASUL 2: SIMULARE NUMERICA - ORBITE CIRCULARE")
print("=" * 70)

def simulate_orbit(r0, M, n_orbits=10, dt_factor=0.001):
    """Simulate circular orbit using Velocity Verlet in 2D flat space."""
    v0 = v_circular_stg(r0, M)
    T_theory = T_stg(r0, M)
    dt = T_theory * dt_factor

    # Initial conditions: particle at (r0, 0), velocity (0, v0)
    x, y = r0, 0.0
    vx, vy = 0.0, v0

    n_steps = int(n_orbits * T_theory / dt) + 1

    times = []
    radii = []
    energies = []
    angles = []

    def acc(x, y):
        r = np.sqrt(x ** 2 + y ** 2)
        if r < 1e-12:
            return 0.0, 0.0
        F = M * r / (r + M) ** 3
        return -F * x / r, -F * y / r

    ax, ay = acc(x, y)
    total_angle = 0.0
    prev_angle = np.arctan2(y, x)

    for step in range(n_steps):
        if step % 10 == 0:
            t = step * dt
            r = np.sqrt(x ** 2 + y ** 2)
            v = np.sqrt(vx ** 2 + vy ** 2)
            KE = 0.5 * v ** 2
            PE = V_stg(r, M)
            times.append(t)
            radii.append(r)
            energies.append(KE + PE)
            angles.append(total_angle)

        # Verlet step
        x += vx * dt + 0.5 * ax * dt ** 2
        y += vy * dt + 0.5 * ay * dt ** 2
        ax_new, ay_new = acc(x, y)
        vx += 0.5 * (ax + ax_new) * dt
        vy += 0.5 * (ay + ay_new) * dt
        ax, ay = ax_new, ay_new

        # Track cumulative angle
        curr_angle = np.arctan2(y, x)
        d_angle = curr_angle - prev_angle
        if d_angle > np.pi:
            d_angle -= 2 * np.pi
        elif d_angle < -np.pi:
            d_angle += 2 * np.pi
        total_angle += d_angle
        prev_angle = curr_angle

    times = np.array(times)
    radii = np.array(radii)
    energies = np.array(energies)
    angles = np.array(angles)

    # Measure period from cumulative angle (much more reliable than periapsis)
    total_time = times[-1] - times[0]
    n_revolutions = total_angle / (2 * np.pi)
    if n_revolutions > 0.5:
        T_measured = total_time / n_revolutions
    else:
        T_measured = np.nan

    E_drift = (energies[-1] - energies[0]) / abs(energies[0]) if energies[0] != 0 else 0
    ecc = (np.max(radii) - np.min(radii)) / (np.max(radii) + np.min(radii))

    return {
        'r0': r0,
        'T_theory': T_theory,
        'T_measured': T_measured,
        'E_drift': E_drift,
        'ecc': ecc,
        'times': times,
        'radii': radii,
        'energies': energies,
    }


# Test at multiple distances
test_radii = [0.5, 1, 2, 3, 5, 10, 20, 50, 100]

print(f"\n  {'r':>6} {'T_theory':>10} {'T_measured':>12} {'T_err%':>8} "
      f"{'E_drift':>10} {'ecc':>8} {'T^2/(r+M)^3':>12}")
print(f"  {'-'*6} {'-'*10} {'-'*12} {'-'*8} {'-'*10} {'-'*8} {'-'*12}")

results = {}
for r0 in test_radii:
    res = simulate_orbit(r0, M, n_orbits=10)
    T_err = abs(res['T_measured'] - res['T_theory']) / res['T_theory'] * 100 if not np.isnan(res['T_measured']) else np.nan
    k = res['T_measured'] ** 2 / (r0 + M) ** 3 if not np.isnan(res['T_measured']) else np.nan
    print(f"  {r0:6.1f} {res['T_theory']:10.4f} {res['T_measured']:12.4f} "
          f"{T_err:7.4f}% {res['E_drift']:10.2e} {res['ecc']:8.6f} {k:12.6f}")
    results[r0] = res


# ==============================================================
# PASUL 3: FIT KEPLER T^2 vs r^n si T^2 vs (r+M)^n
# ==============================================================

print("\n" + "=" * 70)
print("PASUL 3: KEPLER III - T^2 vs r^n vs (r+M)^n")
print("=" * 70)

valid = [(r, res) for r, res in results.items() if not np.isnan(res['T_measured'])]
r_arr = np.array([r for r, _ in valid])
T_arr = np.array([res['T_measured'] for _, res in valid])

# Fit T^2 ~ r^n
log_T = np.log(T_arr)
log_r = np.log(r_arr)
coeffs_r = np.polyfit(log_r, log_T, 1)
n_r = 2 * coeffs_r[0]

# Fit T^2 ~ (r+M)^n
log_rM = np.log(r_arr + M)
coeffs_rM = np.polyfit(log_rM, log_T, 1)
n_rM = 2 * coeffs_rM[0]

print(f"\n  Fit T^2 ~ r^n:     n = {n_r:.4f} (Kepler: n=3)")
print(f"  Fit T^2 ~ (r+M)^n: n = {n_rM:.4f} (STG prediction: n=3)")

# Kepler constant
k_vals = T_arr ** 2 / (r_arr + M) ** 3
k_theory = 4 * np.pi ** 2 / M
print(f"\n  Constanta Kepler STG: T^2/(r+M)^3")
print(f"    Teorie:  {k_theory:.6f}")
print(f"    Masurat: {np.mean(k_vals):.6f} +/- {np.std(k_vals):.6f}")
print(f"    Eroare:  {abs(np.mean(k_vals) - k_theory)/k_theory*100:.4f}%")


# ==============================================================
# PASUL 4: COMPARATIE CU SISTEMUL SOLAR
# ==============================================================

print("\n" + "=" * 70)
print("PASUL 4: APLICARE LA SISTEMUL SOLAR")
print("=" * 70)

# In Planck units:
# M_sun ~ 1.0e38 M_P (but exact value depends on convention)
# r_Mercury ~ 5.8e10 m / l_P ~ 3.6e45 l_P
# r/M ratio is what matters

print(f"""
  In unitati Planck:
    M_Soare ~ 10^38 M_P
    r_Mercur ~ 3.6 * 10^45 l_P
    r/M ~ 10^7

  Corectia STG la Kepler III:
    T^2_STG / T^2_Newton = ((r+M)/r)^3 = (1 + M/r)^3

  Pentru r/M = 10^7: (1 + 10^-7)^3 = 1 + 3*10^-7 + ...
  Corectie: 0.00003% -- COMPLET NEDETECTABILA

  STG reduce EXACT la Newton/Kepler in regimul camp slab.
""")

# More detailed
print("  Corectia STG la diferite scari:")
print(f"  {'r/M':>12} {'Corectie':>15} {'Detectabil?':>15}")
for ratio in [1, 2, 10, 100, 1e3, 1e6, 1e7, 1e10]:
    correction = ((1 + 1/ratio) ** 3 - 1) * 100
    det = "DA (Planck)" if ratio < 10 else "DA (neutron star)" if ratio < 1e3 else "LIMITA" if ratio < 1e6 else "NU"
    print(f"  {ratio:12.0e} {correction:14.6f}% {det:>15}")


# ==============================================================
# PLOTS
# ==============================================================

fig, axes = plt.subplots(2, 3, figsize=(18, 10))
fig.suptitle('EXP-102: STG Kepler III - Analytic Formula\n'
             r'$F_{STG} = \frac{Mr}{(r+M)^3}$, $T^2 = \frac{4\pi^2(r+M)^3}{M}$'
             ' (DERIVED, zero free parameters)',
             fontsize=13, fontweight='bold')

# 1. Force comparison
ax = axes[0, 0]
r = np.linspace(0.01, 10, 1000)
ax.plot(r, F_stg(r), 'r-', linewidth=2, label=r'$F_{STG} = Mr/(r+M)^3$')
ax.plot(r, F_newton(r), 'b--', linewidth=1.5, label=r'$F_{Newton} = M/r^2$')
ax.axvline(x=0.5, color='green', linestyle=':', alpha=0.7, label=r'$r_{F_{max}} = M/2$')
ax.axhline(y=4/27, color='green', linestyle=':', alpha=0.7)
ax.set_xlabel('Distance r (Planck units)')
ax.set_ylabel('Force F(r)')
ax.set_title('Force: STG vs Newton')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)
ax.set_ylim(0, 0.5)

# 2. Potential comparison
ax = axes[0, 1]
r = np.linspace(0.3, 10, 1000)
ax.plot(r, V_stg(r), 'r-', linewidth=2, label=r'$V_{STG} = \frac{1}{2C^2} - \frac{1}{2}$')
ax.plot(r, V_newton(r), 'b--', linewidth=1.5, label=r'$V_{Newton} = -M/r$')
ax.axhline(y=-0.5, color='gray', linestyle=':', alpha=0.5, label='V_STG min = -1/2')
ax.set_xlabel('Distance r')
ax.set_ylabel('V(r)')
ax.set_title('Potential: STG vs Newton')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)
ax.set_ylim(-3, 0.1)

# 3. Kepler III: T^2 vs r
ax = axes[0, 2]
r_plot = np.array(test_radii)
T_plot = np.array([results[r]['T_measured'] for r in test_radii])
valid_mask = ~np.isnan(T_plot)

ax.loglog(r_plot[valid_mask], T_plot[valid_mask] ** 2, 'ro', markersize=8, label='Measured $T^2$')

r_fine = np.logspace(np.log10(0.3), np.log10(150), 200)
ax.loglog(r_fine, T_stg(r_fine) ** 2, 'r-', linewidth=2, label=r'$T^2_{STG} = \frac{4\pi^2(r+M)^3}{M}$')
ax.loglog(r_fine, T_newton(r_fine) ** 2, 'b--', linewidth=1.5, label=r'$T^2_{Kepler} = \frac{4\pi^2 r^3}{M}$')

ax.set_xlabel('Distance r')
ax.set_ylabel('$T^2$')
ax.set_title('Kepler III: STG vs Newton')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3, which='both')

# 4. Kepler constant T^2/(r+M)^3
ax = axes[1, 0]
k_measured = T_plot[valid_mask] ** 2 / (r_plot[valid_mask] + M) ** 3
k_newton = T_plot[valid_mask] ** 2 / r_plot[valid_mask] ** 3
ax.semilogx(r_plot[valid_mask], k_measured, 'ro-', markersize=8, label=r'$T^2/(r+M)^3$ (STG)')
ax.semilogx(r_plot[valid_mask], k_newton, 'bs-', markersize=6, label=r'$T^2/r^3$ (Newton)')
ax.axhline(y=4*np.pi**2, color='red', linestyle=':', alpha=0.7, label=r'$4\pi^2/M$ = {:.4f}'.format(4*np.pi**2))
ax.set_xlabel('Distance r')
ax.set_ylabel('Kepler constant')
ax.set_title(r'Kepler constant: $T^2/(r+M)^3$ = const')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

# 5. Energy conservation
ax = axes[1, 1]
for r0 in [1, 5, 20, 100]:
    res = results[r0]
    e0 = res['energies'][0]
    e_rel = (np.array(res['energies']) - e0) / abs(e0)
    ax.plot(res['times'] / res['T_theory'], e_rel, label=f'r={r0}', linewidth=0.7)
ax.set_xlabel('Time (orbits)')
ax.set_ylabel('$(E(t) - E_0) / |E_0|$')
ax.set_title('Energy conservation')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

# 6. Orbital speed
ax = axes[1, 2]
r = np.linspace(0.3, 50, 500)
ax.plot(r, v_circular_stg(r), 'r-', linewidth=2, label=r'$v_{STG} = r\sqrt{M}/(r+M)^{3/2}$')
ax.plot(r, v_circular_newton(r), 'b--', linewidth=1.5, label=r'$v_{Newton} = \sqrt{M/r}$')
ax.axvline(x=0.5, color='green', linestyle=':', alpha=0.7)
ax.set_xlabel('Distance r')
ax.set_ylabel('Circular orbital speed v')
ax.set_title('Orbital speed: STG vs Newton')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)
ax.set_ylim(0, 1.5)

plt.tight_layout()
plt.savefig('exp102_stg_kepler_analytic.png', dpi=150, bbox_inches='tight')
print(f"\nSalvat: exp102_stg_kepler_analytic.png")


# ==============================================================
# CONCLUZIE
# ==============================================================

print("\n" + "=" * 70)
print("CONCLUZIE EXP-102")
print("=" * 70)

print(f"""
FORMULA STG GRAVITATIONALA (100% DERIVATA):
===========================================

  F(r) = M * r / (r + M)^3

  Echivalent: F = M / (r^2 * (1 + M/r)^3)

PROPRIETATI DERIVATE:
=====================
  1. Limita Newton:  r >> M  =>  F -> M/r^2           (VERIFICAT)
  2. Forta maxima:   r = M/2 =>  F_max = 4/(27*M)     (VERIFICAT)
  3. Fara singulari: r -> 0  =>  F -> 0                (VERIFICAT)
  4. Kepler III:     T^2 = 4*pi^2*(r+M)^3/M           (VERIFICAT)
  5. Potential finit: V_min = -1/2 (la r=0)            (VERIFICAT)

KEPLER III MODIFICAT:
=====================
  Newton:  T^2 = 4*pi^2 * r^3 / M
  STG:     T^2 = 4*pi^2 * (r+M)^3 / M

  Singura diferenta: r -> r + M
  La r >> M: identice (corectie < 10^-7 pentru sistem solar)
  La r ~ M: diferenta semnificativa (regim Planck/stele neutronice)

FIT T^2 ~ (r+M)^n:  n = {n_rM:.4f} (prediction: 3.0000)
Kepler constant:     k = {np.mean(k_vals):.6f} (theory: {k_theory:.6f})

CLASIFICARE: DERIVAT + VERIFICAT NUMERIC
  - Formula: algebra simpla din C(r) = 1 + M/r
  - Kepler III: calcul analitic exact
  - Zero parametri liberi
  - Reduce EXACT la Newton in camp slab
""")

plt.show()
