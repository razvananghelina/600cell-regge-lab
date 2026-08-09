"""
EXP-105: Origine geometrica a exponentului metric alpha
========================================================

Din exp104: metrica STG generalizata
  ds^2 = -dt^2/C^2 + C^2 dr^2 + C^{2*alpha} r^2 dOmega^2
cu alpha = 1/2 reproduce precesia GR la ordinul 1 in M/r.

INTREBARI:
1. Are alpha=1/2 origine geometrica din 600-cell?
2. Apare alpha_EM (constanta structurii fine) in corectii?
3. La ordinul 2 (2PN), care e diferenta STG vs GR?

ABORDARE:
- Pe 600-cell, pt fiecare vertex la distanta d de masa,
  numara vecini la d-1, d, d+1 (radiali vs angulari)
- Raportul angular/radial ar putea determina alpha
- Verifica conexiunea cu alpha_EM = 1/137.036
"""

import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def pr(s=''):
    sys.stdout.write(str(s) + '\n')
    sys.stdout.flush()

PHI = (1 + np.sqrt(5)) / 2
ALPHA_EM = 7.2973525693e-3  # Fine structure constant

# ==============================================================
# 1. GENERARE 600-CELL SI ANALIZA VECINILOR
# ==============================================================

pr("=" * 70)
pr("EXP-105: ORIGINE GEOMETRICA A EXPONENTULUI METRIC")
pr("=" * 70)

# Generate 600-cell vertices (from hjer/data/generate_data.py - verified correct)
def generate_600cell():
    """Generate 120 vertices of the 600-cell on unit S^3."""
    from itertools import permutations
    vertices = []
    # 16 vertices: (+-0.5, +-0.5, +-0.5, +-0.5)
    for i in range(16):
        signs = np.array([1 if (i >> j) & 1 else -1 for j in range(4)])
        vertices.append(signs * 0.5)
    # 8 vertices: permutations of (+-1, 0, 0, 0)
    for i in range(4):
        v = np.zeros(4); v[i] = 1.0
        vertices.append(v); vertices.append(-v)
    # 96 vertices: EVEN permutations of (phi/2, 0.5, 1/(2*phi), 0) with ALL signs
    v_base = np.array([PHI, 1.0, 1.0/PHI, 0.0]) * 0.5
    even_perms = []
    for p in permutations([0,1,2,3]):
        inv_count = sum(1 for i in range(4) for j in range(i+1,4) if p[i] > p[j])
        if inv_count % 2 == 0:
            even_perms.append(p)
    for p in even_perms:
        perm_vals = v_base[list(p)]
        for i in range(16):
            signs = np.array([1 if (i >> j) & 1 else -1 for j in range(4)])
            vertices.append(perm_vals * signs)
    vertices = np.array(vertices)
    _, idx = np.unique(np.round(vertices, 8), axis=0, return_index=True)
    return vertices[idx]

vertices = generate_600cell()
n_verts = len(vertices)
pr(f"\n  600-cell: {n_verts} vertecsi")

# Build adjacency using Euclidean distance (min dist method from generate_data.py)
dists_sq = np.zeros((n_verts, n_verts))
for i in range(n_verts):
    for j in range(i+1, n_verts):
        d = np.sum((vertices[i] - vertices[j])**2)
        dists_sq[i,j] = d
        dists_sq[j,i] = d
upper = dists_sq[np.triu_indices(n_verts, k=1)]
min_dist_sq = np.min(upper)
adj = (dists_sq < min_dist_sq * 1.1) & (dists_sq > 0)
degree = adj.sum(axis=1)
pr(f"  Grad: min={degree.min()}, max={degree.max()}, mean={degree.mean():.1f}")
pr(f"  Edge dist^2: {min_dist_sq:.6f}, angle: {np.arccos(1-min_dist_sq/2)*180/np.pi:.2f} deg")

# BFS from a vertex
def bfs(start, adj_matrix):
    n = adj_matrix.shape[0]
    dist = np.full(n, -1, dtype=int)
    dist[start] = 0
    queue = [start]
    head = 0
    while head < len(queue):
        v = queue[head]
        head += 1
        for u in range(n):
            if adj_matrix[v, u] and dist[u] == -1:
                dist[u] = dist[v] + 1
                queue.append(u)
    return dist


# ==============================================================
# 2. FRACTIE ANGULARA VS RADIALA PER DISTANTA
# ==============================================================

pr(f"\n{'='*70}")
pr("1. VECINI ANGULARI vs RADIALI pe 600-cell")
pr("=" * 70)

# For each vertex as "mass", compute neighbor statistics at each distance
mass_vertex = 0
dist = bfs(mass_vertex, adj)
max_dist = dist.max()
pr(f"\n  Masa la vertex {mass_vertex}, diametru graf = {max_dist}")

pr(f"\n  {'d':>4} {'N(d)':>6} {'n_in':>6} {'n_out':>6} {'n_ang':>6} {'f_ang':>8} {'f_rad':>8}")
pr(f"  {'-'*4} {'-'*6} {'-'*6} {'-'*6} {'-'*6} {'-'*8} {'-'*8}")

# For vertices at distance d, count neighbors at d-1, d, d+1
angular_fractions = []
for d in range(1, max_dist + 1):
    verts_at_d = np.where(dist == d)[0]
    n_d = len(verts_at_d)
    if n_d == 0:
        continue

    total_in = 0    # neighbors at d-1 (radially inward)
    total_out = 0   # neighbors at d+1 (radially outward)
    total_ang = 0   # neighbors at d (angular, same shell)
    total_neighbors = 0

    for v in verts_at_d:
        neighbors = np.where(adj[v])[0]
        for u in neighbors:
            total_neighbors += 1
            du = dist[u]
            if du == d - 1:
                total_in += 1
            elif du == d + 1:
                total_out += 1
            elif du == d:
                total_ang += 1
            # du could also be d-2 or d+2 etc (shouldn't happen for BFS tree)

    avg_in = total_in / n_d
    avg_out = total_out / n_d
    avg_ang = total_ang / n_d
    f_ang = avg_ang / 12.0  # fraction of 12 neighbors that are angular
    f_rad = (avg_in + avg_out) / 12.0

    angular_fractions.append((d, n_d, avg_in, avg_out, avg_ang, f_ang))
    pr(f"  {d:4d} {n_d:6d} {avg_in:6.2f} {avg_out:6.2f} {avg_ang:6.2f} {f_ang:8.4f} {f_rad:8.4f}")


# Average over all distances (weighted by number of vertices)
total_verts = sum(nv for _, nv, _, _, _, _ in angular_fractions)
avg_f_ang = sum(nv * fa for _, nv, _, _, _, fa in angular_fractions) / total_verts
pr(f"\n  Fractie angulara medie (ponderata): {avg_f_ang:.6f}")
pr(f"  Fractie angulara medie = {avg_f_ang:.4f}")
pr(f"  1/2 = 0.5000")
pr(f"  Diferenta: {avg_f_ang - 0.5:.6f}")


# ==============================================================
# 3. AVERAGE OVER ALL MASS POSITIONS
# ==============================================================

pr(f"\n{'='*70}")
pr("2. MEDIA PESTE TOATE POZITIILE MASEI")
pr("=" * 70)

# Precompute adjacency list (faster than np.where each time)
adj_list = [np.where(adj[v])[0] for v in range(n_verts)]

all_f_ang = []
pr(f"  Computing BFS from all {n_verts} vertices...")
for mass_v in range(n_verts):
    if mass_v % 30 == 0:
        pr(f"    vertex {mass_v}/{n_verts}...")
    d = bfs(mass_v, adj)
    total_n = 0
    total_ang_n = 0
    for v in range(n_verts):
        dv = d[v]
        if dv <= 0:
            continue
        for u in adj_list[v]:
            total_n += 1
            if d[u] == dv:
                total_ang_n += 1
    if total_n > 0:
        all_f_ang.append(total_ang_n / total_n)

mean_f_ang = np.mean(all_f_ang)
std_f_ang = np.std(all_f_ang)
pr(f"\n  Media fractie angulara (toate masele): {mean_f_ang:.6f} +/- {std_f_ang:.6f}")
pr(f"  1/2 = 0.500000")
pr(f"  Diferenta de la 1/2: {mean_f_ang - 0.5:.6f}")


# ==============================================================
# 4. CONEXIUNE CU ALPHA_EM
# ==============================================================

pr(f"\n{'='*70}")
pr("3. CONEXIUNE CU ALPHA_EM")
pr("=" * 70)

pr(f"\n  alpha_EM = {ALPHA_EM:.10f} = 1/{1/ALPHA_EM:.3f}")
pr(f"  20*phi^4 = {20*PHI**4:.6f}")
pr(f"  Formula: 2*pi*alpha^2 - 20*phi^4*alpha + 1 = 0")

# Check various combinations
pr(f"\n  Posibile forme alpha_metric:")
pr(f"  1/2 exact:                          0.500000")
pr(f"  1/2 - alpha_EM:                     {0.5 - ALPHA_EM:.6f}")
pr(f"  1/2 - alpha_EM/pi:                  {0.5 - ALPHA_EM/np.pi:.6f}")
pr(f"  1/2 - alpha_EM/(2*pi):              {0.5 - ALPHA_EM/(2*np.pi):.6f}")
pr(f"  1/2 * (1 - alpha_EM):               {0.5 * (1-ALPHA_EM):.6f}")
pr(f"  1/2 * (1 - 2*alpha_EM):             {0.5 * (1-2*ALPHA_EM):.6f}")
pr(f"  Fractie angulara 600-cell:           {mean_f_ang:.6f}")

# Check if angular fraction relates to alpha_EM
diff = mean_f_ang - 0.5
pr(f"\n  Diferenta f_ang - 1/2 = {diff:.6f}")
pr(f"  alpha_EM = {ALPHA_EM:.6f}")
pr(f"  diff / alpha_EM = {diff/ALPHA_EM:.4f}")
pr(f"  diff * 137 = {diff * 137:.4f}")


# ==============================================================
# 5. ANALIZA 2PN: STG (alpha=1/2) vs GR
# ==============================================================

pr(f"\n{'='*70}")
pr("4. ANALIZA 2PN: STG (alpha=1/2) vs GR")
pr("=" * 70)

pr("""
  La 1PN (ordinul M/r), STG cu alpha=1/2 = GR (EXACT).
  La 2PN (ordinul M^2/r^2), difera.

  STG metric (alpha=1/2):
    g_tt = -r^2/(r+M)^2 = -(1 - 2M/r + 3M^2/r^2 - 4M^3/r^3 + ...)
    g_rr = (1+M/r)^2    = 1 + 2M/r + M^2/r^2 + ...
    g_th = (1+M/r)*r^2  = r^2 + Mr + ...

  GR Schwarzschild:
    g_tt = -(1 - 2M/r)    [exact, no M^2 term!]
    g_rr = 1/(1-2M/r)     = 1 + 2M/r + 4M^2/r^2 + 8M^3/r^3 + ...
    g_th = r^2             [exact]

  Diferente la 2PN:
    g_tt: STG has +3M^2/r^2, GR has 0  => coeff diff = 3
    g_rr: STG has +M^2/r^2, GR has +4M^2/r^2  => coeff diff = -3
    g_th: STG has +M*r, GR has 0  => coeff diff = M*r

  NOTA: Coordonatele STG si GR sunt DIFERITE!
  Pentru comparatie corecta, trebuie transformat la acelasi sistem.
""")

# Numerical: compute precession ratio at different r0 for alpha=1/2
from scipy.integrate import solve_ivp

M = 1.0

def compute_precession_exact(r0, e, alpha_met, n_orbits=10):
    """Exact numerical precession for general alpha metric."""
    u0 = 1.0/r0
    C0 = 1 + M*u0

    # L^2 for circular orbit
    if r0 <= alpha_met * M:
        return np.nan
    L2 = M * r0**(2-2*alpha_met) * (r0+M)**(2*alpha_met) / (r0 - alpha_met*M)

    # Energy at circular orbit
    E2 = (1 + L2*u0**2/C0**(2*alpha_met)) / C0**2

    def energy_rhs(u):
        C = 1 + M*u
        return E2 - (1 + L2*u**2/C**(2*alpha_met)) / C**2

    def d2u_dphi2(u):
        du = 1e-8 * max(abs(u), 1e-6)
        return (energy_rhs(u + du) - energy_rhs(u - du)) / (2*du) / (2*L2)

    u_peri = u0 * (1 + e)

    def ode(phi, y):
        u, up = y
        if u < 1e-12: u = 1e-12
        return [up, d2u_dphi2(u)]

    def event_peri(phi, y):
        return y[1]
    event_peri.direction = -1

    mstep = min(0.1, 2*np.pi/50)  # at least 50 steps per orbit
    sol = solve_ivp(ode, [0, n_orbits*2*np.pi*1.5], [u_peri, 0.0],
                    events=event_peri, rtol=1e-11, atol=1e-13, max_step=mstep)

    peri_phis = sol.t_events[0]
    if len(peri_phis) >= 3:
        d_phis = np.diff(peri_phis)
        if len(d_phis) > 2: d_phis = d_phis[1:]
        return np.mean(d_phis) - 2*np.pi
    return np.nan


# GR precession (exact via ODE)
def precession_GR_ode(r0, e, n_orbits=10):
    """Exact GR precession from Schwarzschild orbit equation."""
    if r0 <= 3*M:
        return np.nan
    L2 = M * r0**2 / (r0 - 3*M)

    def rhs_gr(u):
        return M/L2 + 3*M*u**2

    u0 = 1.0/r0
    u_peri = u0 * (1 + e)

    def ode(phi, y):
        u, up = y
        return [up, rhs_gr(u) - u]

    def event_peri(phi, y):
        return y[1]
    event_peri.direction = -1

    mstep = min(0.1, 2*np.pi/50)
    sol = solve_ivp(ode, [0, n_orbits*2*np.pi*1.5], [u_peri, 0.0],
                    events=event_peri, rtol=1e-11, atol=1e-13, max_step=mstep)

    peri_phis = sol.t_events[0]
    if len(peri_phis) >= 3:
        d_phis = np.diff(peri_phis)
        if len(d_phis) > 2: d_phis = d_phis[1:]
        return np.mean(d_phis) - 2*np.pi
    return np.nan


e_test = 0.05
pr(f"\n  Eccentricitate e = {e_test}")
pr(f"  {'r0/M':>8} {'STG(a=0.5)':>12} {'GR':>12} {'STG/GR':>10} {'diff':>12} {'diff*r0':>10}")
pr(f"  {'-'*8} {'-'*12} {'-'*12} {'-'*10} {'-'*12} {'-'*10}")

for r0_ratio in [20, 50, 100, 200]:
    pr(f"    r0={r0_ratio}M...")
    r0 = r0_ratio * M
    p_stg = compute_precession_exact(r0, e_test, 0.5, n_orbits=5)
    p_gr = precession_GR_ode(r0, e_test, n_orbits=5)
    if not np.isnan(p_stg) and not np.isnan(p_gr) and abs(p_gr) > 1e-15:
        ratio = p_stg / p_gr
        diff = p_stg - p_gr
        diff_scaled = diff * r0  # should be ~constant if 2PN
        pr(f"  {r0_ratio:8d} {p_stg:12.8f} {p_gr:12.8f} {ratio:10.6f} {diff:12.2e} {diff_scaled:10.4f}")


# ==============================================================
# 6. WHAT VALUE OF ALPHA MAKES STG MATCH GR AT EACH r0?
# ==============================================================

pr(f"\n{'='*70}")
pr("5. ALPHA OPTIM PER r0 (2PN dependence)")
pr("=" * 70)

from scipy.optimize import brentq

pr(f"  {'r0/M':>8} {'alpha_opt':>12} {'alpha-0.5':>12} {'(a-0.5)*r0/M':>14}")
pr(f"  {'-'*8} {'-'*12} {'-'*12} {'-'*14}")

alpha_opts = []
for r0_ratio in [30, 50, 100, 200]:
    pr(f"    alpha_opt r0={r0_ratio}M...")
    r0 = r0_ratio * M
    p_gr = precession_GR_ode(r0, e_test, n_orbits=5)
    if np.isnan(p_gr) or abs(p_gr) < 1e-15:
        continue

    def target(alpha):
        p = compute_precession_exact(r0, e_test, alpha, n_orbits=5)
        if np.isnan(p):
            return 10.0
        return p / p_gr - 1.0

    try:
        a_opt = brentq(target, 0.3, 0.7, xtol=1e-5)
        diff = a_opt - 0.5
        scaled = diff * r0_ratio
        alpha_opts.append((r0_ratio, a_opt, diff, scaled))
        pr(f"  {r0_ratio:8d} {a_opt:12.6f} {diff:12.6f} {scaled:14.4f}")
    except:
        pr(f"  {r0_ratio:8d} {'FAILED':>12}")

if len(alpha_opts) >= 2:
    # Check if (alpha-0.5)*r0 is constant (2PN behavior)
    scaled_vals = [s for _, _, _, s in alpha_opts]
    pr(f"\n  (alpha-0.5)*r0/M values: {[f'{s:.4f}' for s in scaled_vals]}")
    mean_scaled = np.mean(scaled_vals)
    std_scaled = np.std(scaled_vals)
    pr(f"  Media: {mean_scaled:.4f} +/- {std_scaled:.4f}")

    # Is this related to alpha_EM?
    pr(f"\n  Comparatie cu alpha_EM:")
    pr(f"  (alpha-0.5)*r0/M ~ {mean_scaled:.4f}")
    pr(f"  -pi/2 = {-np.pi/2:.4f}")
    pr(f"  -3/2 = {-1.5:.4f}")
    pr(f"  -alpha_EM * 137 = {-ALPHA_EM * 137:.4f}")
    pr(f"  alpha_EM * 20*phi^4 = {ALPHA_EM * 20*PHI**4:.4f}")


# ==============================================================
# SUMMARY
# ==============================================================

pr(f"\n{'='*70}")
pr("6. SUMAR")
pr("=" * 70)

pr(f"""
  REZULTATE:
  ==========
  1. Fractie angulara medie pe 600-cell: {mean_f_ang:.6f}
     (1/2 = 0.500000, diferenta = {mean_f_ang - 0.5:.6f})

  2. La 1PN: alpha_metric = 1/2 reproduce GR EXACT
     (precesie = 6*pi*M/r0, deflexie = 4M/b)

  3. La 2PN: STG (alpha=1/2) difera de GR
     (corectia depinde de M/r0)

  4. Constanta structurii fine: alpha_EM = {ALPHA_EM:.10f} = 1/{1/ALPHA_EM:.3f}
     Din 600-cell: 2*pi*alpha^2 - 20*phi^4*alpha + 1 = 0
""")

# Plot
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
fig.suptitle('EXP-105: Geometric Origin of Metric Exponent', fontsize=12)

# Angular fraction per distance
ax = axes[0]
dists = [af[0] for af in angular_fractions]
f_angs = [af[5] for af in angular_fractions]
n_verts_per_d = [af[1] for af in angular_fractions]
ax.bar(dists, f_angs, color='steelblue', alpha=0.7, label='Angular fraction')
ax.axhline(y=0.5, color='red', linestyle='--', label='1/2')
ax.axhline(y=mean_f_ang, color='green', linestyle=':', label=f'Mean = {mean_f_ang:.4f}')
ax.set_xlabel('Graph distance d from mass')
ax.set_ylabel('Fraction of angular neighbors')
ax.set_title('Angular vs Radial Neighbors (600-cell)')
ax.legend()
ax.grid(True, alpha=0.3)

# STG/GR ratio vs r0
ax = axes[1]
if len(alpha_opts) >= 2:
    r0s = [ao[0] for ao in alpha_opts]
    aopt = [ao[1] for ao in alpha_opts]
    ax.plot(r0s, aopt, 'ko-', linewidth=2, markersize=6, label='alpha_opt(r0)')
    ax.axhline(y=0.5, color='red', linestyle='--', label='alpha = 1/2')
    ax.set_xlabel('r0 / M')
    ax.set_ylabel('Optimal alpha')
    ax.set_title('2PN: alpha needed to match GR')
    ax.legend()
    ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('exp105_alpha_geometric.png', dpi=150, bbox_inches='tight')
pr(f"Salvat: exp105_alpha_geometric.png")
