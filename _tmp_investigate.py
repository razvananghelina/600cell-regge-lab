import numpy as np

PHI = (1 + np.sqrt(5)) / 2

# ============================================
# PART 1: Exact angular fraction
# ============================================
print('='*70)
print('PART 1: EXACT ANGULAR FRACTION')
print('='*70)

# From exp105 data (for mass at vertex 0):
# d=1: 12 verts, avg_ang=5.00 => total_ang = 60
# d=2: 32 verts, avg_ang=5.625 => total_ang = 180
# d=3: 42 verts, avg_ang=5.714 => total_ang = 240
# d=4: 32 verts, avg_ang=7.50 => total_ang = 240  
# d=5: 1 vert, avg_ang=0 => total_ang = 0

total_ang = 60 + 180 + 240 + 240 + 0
print(f'Total angular neighbor entries (d>=1): {total_ang}')

# Total neighbor entries from d>=1 vertices: 119 * 12 = 1428
total_all = 119 * 12
print(f'Total neighbor entries (d>=1): {total_all}')

f_ang = total_ang / total_all
print(f'f_ang = {total_ang}/{total_all} = {f_ang:.10f}')

# Simplify fraction
from math import gcd
g = gcd(total_ang, total_all)
print(f'f_ang = {total_ang//g}/{total_all//g} = {total_ang/total_all:.10f}')

# Now INCLUDING the mass vertex:
total_ang_with_mass = total_ang + 0  # mass has 12 neighbors, all at d=1 (0 angular)
total_all_with_mass = 120 * 12  # all 120 vertices, each with 12 neighbors
print(f'\nWith mass vertex included:')
print(f'f_ang = {total_ang_with_mass}/{total_all_with_mass} = {total_ang_with_mass/total_all_with_mass:.10f}')
print(f'= {total_ang_with_mass//gcd(total_ang_with_mass,total_all_with_mass)}/{total_all_with_mass//gcd(total_ang_with_mass,total_all_with_mass)}')
print(f'= 1/2 exactly? {total_ang_with_mass/total_all_with_mass == 0.5}')

# INTERPRETATION:
# Each edge connects vertices at distances (d_i, d_j) from mass.
# If d_i == d_j: angular edge. If |d_i-d_j|==1: radial edge.
# Total edges = 720. Angular edges = total_ang/2 = 360.
# 360/720 = 1/2 exactly!
print(f'\nAngular edges = {total_ang//2}')
print(f'Total edges = 720')
print(f'Angular/Total = {total_ang//2}/720 = {total_ang/2/720}')
print(f'EXACTLY HALF of all 600-cell edges are angular!')

# Check what 60/119 is in terms of phi
print(f'\n--- Expressions for 60/119 ---')
print(f'60/119 = {60/119:.10f}')
print(f'1/2 = {0.5:.10f}')
print(f'60/119 - 1/2 = {60/119 - 0.5:.10f}')
print(f'1/(2*119) = {1/(2*119):.10f}')
print(f'So 60/119 = 1/2 + 1/(2*119) = 1/2 + 1/238')

# Check phi expressions
print(f'\nPhi expressions:')
print(f'phi/(2*phi+1) = {PHI/(2*PHI+1):.10f}')
print(f'5/phi^4 = {5/PHI**4:.10f}')  
print(f'1/(2-1/120) = {1/(2-1/120):.10f}')
print(f'60/119 = {60/119:.10f}')

# ============================================
# PART 2: 2PN coefficient analysis
# ============================================
print(f'\n{"="*70}')
print('PART 2: ANALYTIC 2PN DERIVATION')
print('='*70)

print("""
For alpha=1/2 STG metric:
  A(r) = 1/C^2,  B(r) = C^2,  F(r) = C * r^2
  where C = 1 + M/r

Geodesic conserved quantities:
  E = A * dt/ds = (dt/ds)/C^2
  L = F * dph/ds = C*r^2 * dph/ds

Normalization: -1 = -E^2/A + B*(dr/ds)^2 + L^2/F

=> (dr/ds)^2 = [E^2*C^2 - 1 - L^2*u^2/C] / C^2
   where u = 1/r

Using u substitution with dphi:
  dr/dphi = (dr/ds)/(dphi/ds) = (dr/ds)*F/L = (dr/ds)*C/(L*u^2)
  du/dphi = -u^2 * dr/dphi... 

Let me derive (du/dphi)^2 directly.

From normalization:
  C^2*(dr/ds)^2 = E^2*C^2 - 1 - L^2*u^2/C

dphi/ds = L/F = L*u^2/C

So dr/dphi = (dr/ds)/(dphi/ds) = (dr/ds)*C/(L*u^2)

(dr/dphi)^2 = C^2*(dr/ds)^2 * C^2 / (L^2*u^4)
            = [E^2*C^2 - 1 - L^2*u^2/C] * C^2 / (L^2*u^4)

du/dphi = -u^2/r^2 * dr/dphi = -u^2 * dr/dphi  (since u = 1/r)

Wait: du = -dr/r^2 = -u^2 dr
So du/dphi = -u^2 * dr/dphi

(du/dphi)^2 = u^4 * (dr/dphi)^2
= u^4 * [E^2*C^2 - 1 - L^2*u^2/C] * C^2 / (L^2*u^4)
= C^2/L^2 * [E^2*C^2 - 1 - L^2*u^2/C]

Let h(u) = (du/dphi)^2 = C^2/L^2 * [E^2*C^2 - 1 - L^2*u^2/C]

where C = 1 + Mu.

d^2u/dphi^2 = h'(u)/2

h(u) = E^2*C^4/L^2 - C^2/L^2 - u^2*C

dC/du = M

h'(u) = E^2*4C^3*M/L^2 - 2C*M/L^2 - 2u*C - u^2*M

d^2u/dphi^2 = (1/2)*[4E^2*M*C^3/L^2 - 2MC/L^2 - 2uC - Mu^2]
""")

# Let's expand perturbatively around circular orbit
# Circular orbit at r=r0: u0 = 1/r0, C0 = 1 + M*u0
# Conditions: h(u0) = 0 and h'(u0) = 0

from scipy.integrate import solve_ivp

M = 1.0

def precession_stg_half(r0, e=0.05, n_orbits=5):
    u0 = 1.0/r0
    C0 = 1 + M*u0
    if r0 <= M/2:
        return float('nan')
    L2 = M * r0 * (r0 + M) / (r0 - M/2)
    E2 = (1 + L2*u0**2/C0) / C0**2
    
    def h(u):
        C = 1 + M*u
        return C**2/L2 * (E2*C**2 - 1 - L2*u**2/C)
    
    def d2u_dphi2(u):
        du = 1e-9 * max(abs(u), 1e-8)
        return (h(u+du) - h(u-du)) / (4*du)  # h'/2
    
    u_peri = u0 * (1 + e)
    
    def ode(phi, y):
        u, up = y
        if u < 1e-12: u = 1e-12
        return [up, d2u_dphi2(u)]
    
    def event(phi, y):
        return y[1]
    event.direction = -1
    
    sol = solve_ivp(ode, [0, n_orbits*2*np.pi*1.5], [u_peri, 0.0],
                    events=event, rtol=1e-12, atol=1e-14, max_step=0.1)
    
    phis = sol.t_events[0]
    if len(phis) >= 3:
        dphis = np.diff(phis)
        if len(dphis) > 2: dphis = dphis[1:]
        return np.mean(dphis) - 2*np.pi
    return float('nan')

def precession_gr(r0, e=0.05, n_orbits=5):
    if r0 <= 3*M:
        return float('nan')
    L2 = M * r0**2 / (r0 - 3*M)
    
    def rhs(u):
        return M/L2 + 3*M*u**2
    
    u0 = 1.0/r0
    u_peri = u0 * (1 + e)
    
    def ode(phi, y):
        u, up = y
        return [up, rhs(u) - u]
    
    def event(phi, y):
        return y[1]
    event.direction = -1
    
    sol = solve_ivp(ode, [0, n_orbits*2*np.pi*1.5], [u_peri, 0.0],
                    events=event, rtol=1e-12, atol=1e-14, max_step=0.1)
    
    phis = sol.t_events[0]
    if len(phis) >= 3:
        dphis = np.diff(phis)
        if len(dphis) > 2: dphis = dphis[1:]
        return np.mean(dphis) - 2*np.pi
    return float('nan')

print(f'\n--- 2PN Numerical Analysis ---')
print(f'{"r0/M":>8} {"STG":>14} {"GR":>14} {"ratio":>10} {"diff*r0^2":>12}')
print(f'{"-"*8} {"-"*14} {"-"*14} {"-"*10} {"-"*12}')

results = []
for r0_val in [30, 50, 100, 200, 500, 1000]:
    r0 = r0_val * M
    ps = precession_stg_half(r0, e=0.01, n_orbits=8)
    pg = precession_gr(r0, e=0.01, n_orbits=8)
    if not np.isnan(ps) and not np.isnan(pg) and abs(pg) > 1e-15:
        ratio = ps/pg
        diff = ps - pg
        # 1PN: diff ~ 0, 2PN: diff ~ c*(M/r0)^2, so diff*r0^2/M^2 should be const
        diff_scaled = diff * r0_val**2
        results.append((r0_val, ps, pg, ratio, diff_scaled))
        print(f'{r0_val:8d} {ps:14.10f} {pg:14.10f} {ratio:10.6f} {diff_scaled:12.6f}')

if len(results) >= 2:
    scaled = [r[4] for r in results]
    print(f'\ndiff*r0^2/M^2 values: {[f"{s:.4f}" for s in scaled]}')
    # Check if this converges
    print(f'Last value: {scaled[-1]:.6f}')
    print(f'6*pi = {6*np.pi:.6f}')
    print(f'-6*pi = {-6*np.pi:.6f}')
    print(f'-3*pi = {-3*np.pi:.6f}')
    print(f'-2*pi = {-2*np.pi:.6f}')
    print(f'-pi = {-np.pi:.6f}')
    print(f'-9*pi = {-9*np.pi:.6f}')
    print(f'-18*pi = {-18*np.pi:.6f}')
    print(f'-6*pi^2 = {-6*np.pi**2:.6f}')
    print(f'-3*pi^2 = {-3*np.pi**2:.6f}')
    
    # GR 1PN: 6*pi*M/r0
    # If STG = GR*(1 + c2*M/r0), then diff = 6*pi*M/r0 * c2*M/r0 = 6*pi*c2*M^2/r0^2
    # diff*r0^2 = 6*pi*c2
    c2_estimate = scaled[-1] / (6*np.pi)
    print(f'\nc2 = diff*r0^2/(6*pi) = {c2_estimate:.6f}')
    print(f'GR 2PN coefficient for e~0: 51/4 = {51/4}')
    print(f'STG 2PN coefficient: {51/4 + c2_estimate:.6f}')
    print(f'Difference from GR: {c2_estimate:.6f}')
    print(f'= -3 exactly? {abs(c2_estimate + 3) < 0.1}')
    print(f'= -pi? {abs(c2_estimate + np.pi) < 0.1}')
    print(f'= -10/3? {abs(c2_estimate + 10/3) < 0.1}')
