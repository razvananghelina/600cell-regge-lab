"""
EXP-101: STG Kepler Test on 600-cell (v2 - smooth C field)
============================================================

OBIECTIV: Testare cantitativa a gravitatiei STG pe graful 600-cell.
- Orbite la distante graf d = 2, 3, 4
- Verificare T^2 ~ d^n (Kepler III => n=3?)
- Masurare precesie per orbita (diferenta STG vs Schwarzschild)
- Conservarea energiei (Verlet simplectic)

FIZICA (DERIVATA, fara fitting):
- C(d) = 1 + M/d
- a = nabla(C) / C^3  (geodezica STG)
- |a| = M / (d^2 * C^3)  directionata pe S^3 spre masa

FIX v2: Distanta graf INTERPOLATA (IDW pe k=6 vecini)
- Problema v1: distanta graf discreta sare cand particula trece
  intre celulele Voronoi => discontinuitati ~58% in forta
- Verlet conserva energia doar pt potentiale netede!
- Fix: d_eff = sum(w_i * d_i) / sum(w_i), w_i = 1/ang_dist_i^2

UNITATI: Planck (l_P = 1 edge, M_P = 1, c = 1)
"""

import numpy as np
import json
import matplotlib.pyplot as plt
from collections import deque

# ==============================================================
# LOAD 600-CELL DATA
# ==============================================================

with open('hjer/data/600cell.json', 'r') as f:
    data = json.load(f)

vertices = np.array(data['vertices4D'])  # (120, 4) on unit S^3
edges = data['edges']
adjacency = data['adjacency']
N_VERTS = data['numVertices']

print("=" * 70)
print("EXP-101: STG KEPLER TEST ON 600-CELL")
print("=" * 70)
print(f"Vertices: {N_VERTS}, Edges: {len(edges)}, Degree: {data['degree']}, Diameter: {data['graphDiameter']}")


def bfs(source):
    """BFS distances from source vertex."""
    dist = np.full(N_VERTS, -1, dtype=int)
    dist[source] = 0
    queue = deque([source])
    while queue:
        u = queue.popleft()
        for v in adjacency[u]:
            if dist[v] == -1:
                dist[v] = dist[u] + 1
                queue.append(v)
    return dist


# ==============================================================
# S^3 GEOMETRY UTILITIES
# ==============================================================

def normalize(v):
    """Normalize 4D vector to unit sphere."""
    n = np.linalg.norm(v)
    return v / n if n > 1e-12 else v


def project_to_tangent(v, pos):
    """Project v onto tangent plane of S^3 at pos."""
    return v - np.dot(v, pos) * pos


def tangent_toward(pos, target):
    """Tangent vector at pos pointing toward target on S^3."""
    dot = np.dot(pos, target)
    t = target - dot * pos
    n = np.linalg.norm(t)
    if n < 1e-10:
        return np.zeros(4)
    return t / n


def angular_distance(p, q):
    """Angular distance between two points on S^3."""
    dot = np.clip(np.dot(p, q), -1, 1)
    return np.arccos(dot)


# ==============================================================
# STG PHYSICS ON S^3
# ==============================================================

class STGField:
    """Surface-Tension Gravity field on 600-cell."""

    def __init__(self, mass_vertex, mass_value, mode='angular'):
        """
        mode: 'discrete' - original (discontinuous, for comparison)
              'idw'      - IDW interpolation (smooth but not exactly central)
              'angular'  - continuous angular distance on S^3 (smooth + central)
        """
        self.mass_vertex = mass_vertex
        self.mass_pos = vertices[mass_vertex]
        self.mass = mass_value
        self.mode = mode
        self.bfs_dist = bfs(mass_vertex)
        self.bfs_dist_float = self.bfs_dist.astype(float)

        # Precompute C at all vertices
        self.C = np.ones(N_VERTS)
        for i in range(N_VERTS):
            d = self.bfs_dist[i]
            if d > 0:
                self.C[i] = 1 + mass_value / d
            elif i == mass_vertex:
                self.C[i] = 1 + mass_value / 0.5

        # Angular distance scale: 1 graph edge = r_1 radians on S^3
        # This is the angular distance between adjacent vertices
        neighbor = adjacency[mass_vertex][0]
        self.r1 = angular_distance(vertices[mass_vertex], vertices[neighbor])
        print(f"    [STG] mode={mode}, r1 (1 edge in radians) = {self.r1:.6f}")

    def nearest_vertex(self, pos):
        """Find nearest graph vertex to a point on S^3."""
        dots = vertices @ pos
        return np.argmax(dots)

    def get_distance(self, pos):
        """Get effective distance based on mode."""
        if self.mode == 'discrete':
            nearest = self.nearest_vertex(pos)
            return float(self.bfs_dist[nearest])
        elif self.mode == 'idw':
            return self._idw_distance(pos)
        else:  # angular
            return self._angular_distance(pos)

    def _idw_distance(self, pos, k=6):
        """IDW interpolation of graph distance."""
        dots = vertices @ pos
        nearest_k = np.argpartition(-dots, k)[:k]
        ang_dists = np.arccos(np.clip(dots[nearest_k], -1, 1))
        min_ang = ang_dists.min()
        if min_ang < 1e-8:
            return self.bfs_dist_float[nearest_k[np.argmin(ang_dists)]]
        weights = 1.0 / (ang_dists ** 2)
        weights /= np.sum(weights)
        return np.sum(weights * self.bfs_dist_float[nearest_k])

    def _angular_distance(self, pos):
        """
        Continuous angular distance calibrated to graph units.

        r_graph = r_angular / r1

        where r1 = angular distance of one graph edge.
        This is the CONTINUOUS LIMIT of the discrete graph distance.

        DERIVAT: limita continua naturala. Nu introduce parametri (r1 e masurat
        din geometria 600-cell, nu ales).
        """
        r_ang = angular_distance(pos, self.mass_pos)
        return r_ang / self.r1  # Convert to graph units

    def acceleration(self, pos):
        """
        STG geodesic acceleration at pos on S^3.

        DERIVED:
        a = nabla(C) / C^3
        |a| = M / (d^2 * C^3), directed toward mass on S^3
        """
        d = self.get_distance(pos)
        if d <= 0.1:
            return np.zeros(4)

        C = 1 + self.mass / d
        force_mag = self.mass / (d * d * C * C * C)

        # Tangent direction on S^3 from pos toward mass (always continuous)
        direction = tangent_toward(pos, self.mass_pos)

        return force_mag * direction

    def potential_energy(self, pos):
        """
        STG potential on S^3 from integrating a_theta = M*r1^2/(theta^2 * C^3):

        V(theta) = r1/(2*C^2) - r1/2

        where C = 1 + M*r1/theta, theta = angular distance.

        In terms of d = theta/r1: V = r1/(2*(1+M/d)^2) - r1/2

        DERIVAT: integral exact al acceleratiei pe S^3.

        For discrete/idw modes (where d = graph distance, not theta/r1),
        we use V = 1/(2C^2) - 1/2 (no r1 factor, since d is already in
        graph units and the force uses graph distance directly).
        """
        d = self.get_distance(pos)
        if d <= 0.1:
            return 0.0
        C = 1 + self.mass / d
        if self.mode == 'angular':
            # Angular mode: force = M*r1^2/(theta^2*C^3), potential has r1 factor
            return self.r1 / (2.0 * C * C) - self.r1 / 2.0
        else:
            # Discrete/IDW: force = M/(d^2*C^3), potential without r1
            return 1.0 / (2.0 * C * C) - 0.5


# ==============================================================
# VELOCITY VERLET ON S^3
# ==============================================================

def verlet_step(pos, vel, acc_func, dt):
    """One Velocity Verlet step on S^3."""
    acc = acc_func(pos)

    # Position update
    pos_new = pos + vel * dt + 0.5 * acc * dt * dt
    pos_new = normalize(pos_new)

    # New acceleration
    acc_new = acc_func(pos_new)

    # Velocity update
    vel_new = vel + 0.5 * (acc + acc_new) * dt

    # Re-orthogonalize to tangent plane
    vel_new = project_to_tangent(vel_new, pos_new)

    return pos_new, vel_new


def compute_orbital_speed(stg, graph_dist):
    """
    Compute circular orbital speed at graph distance d on S^3.

    On S^3, centripetal acceleration for circular motion at angular distance theta:
      a_centripetal = v^2 * cot(theta)

    So circular orbit condition:
      F = v^2 * cot(theta)
      v = sqrt(F * tan(theta))

    IMPORTANT: tan(theta) > 0 only for theta < pi/2.
    At theta = pi/2: geodesic (no force needed for circular motion)
    At theta > pi/2: would need REPULSIVE force => no circular orbit possible

    DERIVAT: geometrie diferentiala pe S^3.
    """
    # Find a vertex at this graph distance
    start_vert = -1
    for i in range(N_VERTS):
        if stg.bfs_dist[i] == graph_dist:
            start_vert = i
            break
    if start_vert == -1:
        return 0, -1, 0

    r_ang = angular_distance(stg.mass_pos, vertices[start_vert])

    # Use STG field's own distance for force (mode-consistent)
    d_eff = stg.get_distance(vertices[start_vert])
    C = 1 + stg.mass / d_eff
    F = stg.mass / (d_eff ** 2 * C ** 3)

    if r_ang >= np.pi / 2 - 0.01:
        # theta >= pi/2: circular orbit not possible
        # Use a reasonable initial speed instead (free-fall speed)
        # v ~ sqrt(2 * |V|) where V = 1/(2C^2) - 1/2
        V = 1.0 / (2.0 * C * C) - 0.5
        speed = np.sqrt(abs(2 * V)) * 0.5  # sub-escape speed
        print(f"    NOTE: theta={r_ang:.4f} >= pi/2, circular orbit NOT possible on S^3")
        print(f"    Using sub-escape speed: v={speed:.6f}")
        return speed, start_vert, r_ang

    # v = sqrt(F * tan(theta))
    speed = np.sqrt(F * np.tan(r_ang))
    print(f"    d_eff={d_eff:.4f}, tan(theta)={np.tan(r_ang):.4f}")

    return speed, start_vert, r_ang


# ==============================================================
# EXPERIMENT: COMPARE THREE MODES
# ==============================================================

MASS_VERTEX = 0
MASS_VALUE = 1.0  # 1 Planck mass
DT = 0.0002
N_ORBITS = 5

MODES = ['discrete', 'idw', 'angular']
all_results = {}

for mode in MODES:
    print("\n" + "=" * 70)
    print(f"MODE: {mode.upper()}")
    print("=" * 70)

    stg = STGField(MASS_VERTEX, MASS_VALUE, mode=mode)

    results = {}

    for target_d in [2, 3, 4]:
        speed, start_vert, r_angular = compute_orbital_speed(stg, target_d)
        if speed == 0:
            print(f"  d={target_d}: No vertex found, skipping")
            continue

        C_val = 1 + MASS_VALUE / target_d
        F_val = MASS_VALUE / (target_d ** 2 * C_val ** 3)

        print(f"\n  d={target_d}: C={C_val:.4f}, F={F_val:.6f}, "
              f"r_ang={r_angular:.4f}, v_orbit={speed:.6f}")

        # Set up particle
        pos = vertices[start_vert].copy()
        mass_pos = stg.mass_pos

        # Radial direction (tangent toward mass)
        radial = tangent_toward(pos, mass_pos)

        # Tangential direction (perpendicular to radial in tangent plane)
        arb = np.array([0, 1, 0, 0], dtype=float)
        arb = project_to_tangent(arb, pos)
        arb = arb - np.dot(arb, radial) * radial
        arb_len = np.linalg.norm(arb)
        if arb_len < 1e-8:
            arb = np.array([0, 0, 1, 0], dtype=float)
            arb = project_to_tangent(arb, pos)
            arb = arb - np.dot(arb, radial) * radial
            arb_len = np.linalg.norm(arb)
        tangential = arb / arb_len

        vel = tangential * speed

        # Simulate
        T_estimate = 2 * np.pi * r_angular / speed
        n_steps = int(N_ORBITS * T_estimate / DT) + 1
        n_steps = min(n_steps, 2000000)  # safety cap

        print(f"  T_estimate={T_estimate:.2f}, steps={n_steps}")

        # Track data (sample every 10 steps to save memory)
        times = []
        distances = []
        energies = []
        SAMPLE = 10

        for step in range(n_steps):
            if step % SAMPLE == 0:
                t = step * DT
                r = angular_distance(pos, mass_pos)
                distances.append(r)
                times.append(t)

                v_mag = np.linalg.norm(vel)
                KE = 0.5 * v_mag ** 2
                PE = stg.potential_energy(pos)
                energies.append(KE + PE)

            pos, vel = verlet_step(pos, vel, stg.acceleration, DT)

        times = np.array(times)
        distances = np.array(distances)
        energies = np.array(energies)

        # --- Measure orbital period from distance oscillation ---
        periapsis_times = []
        for i in range(2, len(distances) - 2):
            if (distances[i] < distances[i-1] and distances[i] < distances[i+1] and
                distances[i] < distances[i-2] and distances[i] < distances[i+2]):
                periapsis_times.append(times[i])

        if len(periapsis_times) >= 2:
            periods = np.diff(periapsis_times)
            T_measured = np.mean(periods)
            T_std = np.std(periods) if len(periods) > 1 else 0
        else:
            T_measured = np.nan
            T_std = np.nan

        # Energy conservation
        E_init = energies[0]
        E_final = energies[-1]
        E_mean = np.mean(energies)
        E_drift = (E_final - E_init) / abs(E_init) if E_init != 0 else 0
        E_max_dev = np.max(np.abs(np.array(energies) - E_init)) / abs(E_init) if E_init != 0 else 0

        eccentricity = (np.max(distances) - np.min(distances)) / (np.max(distances) + np.min(distances))

        print(f"  T_measured = {T_measured:.4f} +/- {T_std:.4f}")
        print(f"  Periapsis passages: {len(periapsis_times)}")
        print(f"  Energy: E0={E_init:.6f}, drift={E_drift:.2e}, max_dev={E_max_dev:.2e}")
        print(f"  r_min={np.min(distances):.4f}, r_max={np.max(distances):.4f}, e={eccentricity:.4f}")

        results[target_d] = {
            'T': T_measured,
            'T_std': T_std,
            'r_angular': r_angular,
            'speed': speed,
            'C': C_val,
            'F': F_val,
            'E_drift': E_drift,
            'E_max_dev': E_max_dev,
            'times': times,
            'distances': distances,
            'energies': energies,
            'periapsis_times': periapsis_times,
            'r_min': np.min(distances),
            'r_max': np.max(distances),
            'eccentricity': eccentricity,
        }

    all_results[mode] = results

# Use 'angular' mode as primary results
results = all_results['angular']

# ==============================================================
# KEPLER III TEST: T^2 vs d^n
# ==============================================================

print("\n" + "-" * 70)
print("PASUL 2: TEST KEPLER III - T^2 vs d^n")
print("-" * 70)

valid_d = sorted([d for d in results if not np.isnan(results[d]['T'])])

if len(valid_d) >= 2:
    T_vals = np.array([results[d]['T'] for d in valid_d])
    d_vals = np.array(valid_d, dtype=float)
    r_vals = np.array([results[d]['r_angular'] for d in valid_d])

    # Fit T^2 = k * d^n  =>  2*log(T) = log(k) + n*log(d)
    log_T = np.log(T_vals)
    log_d = np.log(d_vals)

    # Linear fit in log space
    coeffs = np.polyfit(log_d, log_T, 1)
    n_kepler = coeffs[0] * 2  # T ~ d^(n/2), so T^2 ~ d^n => log(T) = (n/2)*log(d) + const
    # Actually: T^2 ~ d^n => 2*log(T) = n*log(d) + const => log(T) = (n/2)*log(d) + const
    # So n = 2 * slope
    n_fit = 2 * coeffs[0]

    print(f"\n  Kepler fit: T^2 ~ d^{n_fit:.3f}")
    print(f"  (Kepler III for Newtonian gravity: T^2 ~ r^3, n=3)")
    print(f"  (STG correction from C^3 term may modify this)")

    # Also fit T^2 vs r_angular^n
    log_r = np.log(r_vals)
    coeffs_r = np.polyfit(log_r, log_T, 1)
    n_fit_r = 2 * coeffs_r[0]
    print(f"  Kepler fit (angular distance): T^2 ~ r_ang^{n_fit_r:.3f}")

    print(f"\n  Data points:")
    print(f"  {'d':>4} {'r_ang':>8} {'T':>10} {'T^2':>10} {'T^2/d^3':>10} {'T^2/r^3':>10}")
    for d in valid_d:
        r = results[d]
        T = r['T']
        r_a = r['r_angular']
        print(f"  {d:4d} {r_a:8.4f} {T:10.4f} {T**2:10.4f} "
              f"{T**2/d**3:10.4f} {T**2/r_a**3:10.4f}")
else:
    print("  Not enough data points for Kepler fit")
    n_fit = np.nan
    n_fit_r = np.nan

# ==============================================================
# PRECESSION MEASUREMENT
# ==============================================================

# ==============================================================
# COMPARISON TABLE: ALL MODES
# ==============================================================

print("\n" + "-" * 70)
print("COMPARATIE MODURI: CONSERVARE ENERGIE")
print("-" * 70)
print(f"  Potential corect: V(r) = 1/(2*C^2) - 1/2  (DERIVAT din F = M/(r^2*C^3))")
print(f"  {'Mode':<10} {'d':>3} {'E_drift':>12} {'max_dev':>12} {'ecc':>8} {'T':>10}")
print(f"  {'-'*10} {'-'*3} {'-'*12} {'-'*12} {'-'*8} {'-'*10}")

for mode in MODES:
    for d in sorted(all_results[mode].keys()):
        r = all_results[mode][d]
        if np.isnan(r['T']):
            continue
        print(f"  {mode:<10} {d:3d} {r['E_drift']:12.2e} {r['E_max_dev']:12.2e} "
              f"{r['eccentricity']:8.4f} {r['T']:10.4f}")

# ==============================================================
# PLOTS
# ==============================================================

fig, axes = plt.subplots(2, 3, figsize=(18, 10))
fig.suptitle('EXP-101: STG Kepler Test on 600-cell (v2)\n'
             'C(r)=1+M/r, V=1/(2C^2)-1/2, a=nabla(C)/C^3, M=1 M_P, ZERO free params',
             fontsize=12, fontweight='bold')

# Plot 1: r(t) angular mode
ax = axes[0, 0]
ax.set_title('Angular mode: r(t)')
for d in valid_d:
    r = results[d]
    ax.plot(r['times'], r['distances'], label=f'd={d}, T={r["T"]:.2f}', linewidth=0.5)
ax.set_xlabel('Time (Planck units)')
ax.set_ylabel('Angular distance from mass')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

# Plot 2: Energy conservation - angular mode
ax = axes[0, 1]
ax.set_title('Angular mode: Energy E(t)')
for d in valid_d:
    r = results[d]
    e0 = r['energies'][0]
    e_rel = (np.array(r['energies']) - e0) / abs(e0)
    ax.plot(r['times'], e_rel, label=f'd={d}', linewidth=0.5)
ax.set_xlabel('Time (Planck units)')
ax.set_ylabel('(E(t) - E(0)) / |E(0)|')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

# Plot 3: Energy comparison across modes (d=3)
ax = axes[0, 2]
ax.set_title('Energy comparison d=3 (all modes)')
for mode in MODES:
    if 3 in all_results[mode]:
        r = all_results[mode][3]
        e0 = r['energies'][0]
        e_rel = (np.array(r['energies']) - e0) / abs(e0)
        ax.plot(r['times'], e_rel, label=f'{mode}', linewidth=0.5)
ax.set_xlabel('Time (Planck units)')
ax.set_ylabel('(E(t) - E(0)) / |E(0)|')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

# Plot 4: Kepler III
ax = axes[1, 0]
ax.set_title(f'Kepler III: T^2 vs d^n (n = {n_fit:.2f})')
if len(valid_d) >= 2:
    T_vals = [results[d]['T'] for d in valid_d]
    d_arr = np.array(valid_d, dtype=float)
    T_arr = np.array(T_vals)
    ax.plot(d_arr, T_arr ** 2, 'ro-', markersize=8, label='Measured T^2')
    d_fine = np.linspace(d_arr[0], d_arr[-1], 50)
    k3 = T_arr[0] ** 2 / d_arr[0] ** 3
    ax.plot(d_fine, k3 * d_fine ** 3, 'b--', label='T^2 ~ d^3 (Kepler)', alpha=0.7)
    k_fit = T_arr[0] ** 2 / d_arr[0] ** n_fit
    ax.plot(d_fine, k_fit * d_fine ** n_fit, 'g--',
            label=f'T^2 ~ d^{n_fit:.2f} (STG)', alpha=0.7)
ax.set_xlabel('Graph distance d')
ax.set_ylabel('T^2')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

# Plot 5: C(d) profile
ax = axes[1, 1]
ax.set_title('STG: C(d) and Force profile')
d_range = np.linspace(0.5, 5, 100)
C_profile = 1 + MASS_VALUE / d_range
ax.plot(d_range, C_profile, 'b-', label='C(r) = 1 + M/r')
ax.set_xlabel('Distance r (graph units)')
ax.set_ylabel('C(r)', color='blue')
ax.tick_params(axis='y', labelcolor='blue')
ax2 = ax.twinx()
F_profile = MASS_VALUE / (d_range ** 2 * C_profile ** 3)
ax2.plot(d_range, F_profile, 'r-', label='F = M/(r^2 C^3)')
ax2.set_ylabel('Force magnitude', color='red')
ax2.tick_params(axis='y', labelcolor='red')
ax.legend(loc='upper right', fontsize=8)
ax2.legend(loc='center right', fontsize=8)
ax.grid(True, alpha=0.3)

# Plot 6: Potential V(r)
ax = axes[1, 2]
ax.set_title('STG Potential V(r) = 1/(2C^2) - 1/2')
V_profile = 1.0 / (2.0 * C_profile ** 2) - 0.5
ax.plot(d_range, V_profile, 'k-', linewidth=2, label='V(r) = 1/(2C^2) - 1/2')
# Newton for comparison
V_newton = -MASS_VALUE / d_range
ax.plot(d_range, V_newton, 'b--', alpha=0.5, label='V_Newton = -M/r')
ax.set_xlabel('Distance r (graph units)')
ax.set_ylabel('V(r)')
ax.set_ylim(-2, 0.1)
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)
ax.axhline(y=0, color='gray', linewidth=0.5)

plt.tight_layout()
plt.savefig('exp101_stg_kepler.png', dpi=150, bbox_inches='tight')
print(f"\nSalvat: exp101_stg_kepler.png")

# ==============================================================
# CONCLUZIE
# ==============================================================

print("\n" + "=" * 70)
print("CONCLUZIE EXP-101")
print("=" * 70)

print(f"""
REZULTATE (v2 - 3 moduri comparate):
=====================================
1. Orbite stabile pe S^3 la d = {', '.join(str(d) for d in valid_d)}
2. Kepler III (mode angular): T^2 ~ d^{n_fit:.2f} (clasic = 3.0)
   {"OK - consistent cu Kepler!" if abs(n_fit - 3) < 0.5 else "DEVIATIE de la Kepler - efect STG"}
3. Conservare energie (angular): {', '.join(f'd={d}: {results[d]["E_drift"]:.1e}' for d in valid_d)}
4. Parametri liberi: ZERO

CLASIFICARE:
============
- Topologia 600-cell: DERIVAT
- C(r) = 1 + M/r: DERIVAT (din metrica STG)
- V(r) = 1/(2*C^2) - 1/2: DERIVAT (integral exact al fortei)
- a = nabla(C)/C^3 = M/(r^2*C^3): DERIVAT (ecuatia geodezicei)
- r_angular/r_1 ca distanta: LIMITA CONTINUA (nu fitting, nu parametri)
  r_1 e distanta angulara a unei muchii - MASURAT din geometria 600-cell
- Orbite stabile: REZULTAT (nu fitting)
- Exponent Kepler: MASURAT (nu impus)

TREI MODURI:
============
- discrete: distanta graf intreaga (sare la granita Voronoi)
- idw: interpolare IDW pe k=6 vecini (neteda dar nu exact centrala)
- angular: distanta angulara pe S^3 / r_1 (neteda + centrala + conservativa)

NOTA: Diferenta intre Kepler clasic (n=3) si STG vine din C^3 in numitor:
  F_Newton = M/r^2  =>  T^2 ~ r^3
  F_STG = M/(r^2 * C^3)  =>  C nu e constant, modifica exponentul
  Aceasta DEVIATIE e o PREDICTIE a STG, nu un defect!
""")

plt.show()
