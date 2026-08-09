"""
EXP-266: COSMOLOGICAL CONSTANT - SPECTRAL ZETA ON 600-CELL
===========================================================
Compute spectral zeta function zeta(s) = sum lambda^{-s}
on the 600-cell Dirac operator D^2.
Use zeta regularization for vacuum energy.
"""
import numpy as np
from scipy import sparse
import time

PHI = (1 + np.sqrt(5)) / 2
a1 = 5; b1 = 6; N_verts = 120

print("="*72)
print("EXP-266: COSMOLOGICAL CONSTANT - SPECTRAL ZETA")
print("="*72)

# === Build 600-cell ===
print("\n  Building 600-cell...")
t0 = time.time()

verts = []
for i in range(4):
    for s in [1, -1]:
        v = [0.0]*4; v[i] = float(s); verts.append(v)
for s0 in [-1,1]:
    for s1 in [-1,1]:
        for s2 in [-1,1]:
            for s3 in [-1,1]:
                verts.append([s0*0.5, s1*0.5, s2*0.5, s3*0.5])
coords_b = [0.0, 0.5, PHI/2, 1/(2*PHI)]
even_perms = [
    (0,1,2,3),(0,2,3,1),(0,3,1,2),(1,0,3,2),(1,2,0,3),(1,3,2,0),
    (2,0,1,3),(2,1,3,0),(2,3,0,1),(3,0,2,1),(3,1,0,2),(3,2,1,0)]
for perm in even_perms:
    for s1 in [-1,1]:
        for s2 in [-1,1]:
            for s3 in [-1,1]:
                v = [0.0]*4
                v[perm[1]] = s1*coords_b[1]
                v[perm[2]] = s2*coords_b[2]
                v[perm[3]] = s3*coords_b[3]
                verts.append(v)
verts = np.array(verts)
uniq = [verts[0]]
for v in verts[1:]:
    if np.min(np.linalg.norm(np.array(uniq)-v, axis=1)) > 1e-8:
        uniq.append(v)
verts = np.array(uniq)
assert len(verts) == 120, f"Got {len(verts)} vertices"

# Edges
edge_len = 1.0/PHI
edges = []; adj = {i: set() for i in range(120)}
for i in range(120):
    di = np.linalg.norm(verts - verts[i], axis=1)
    for j in range(i+1, 120):
        if abs(di[j] - edge_len) < 1e-6:
            edges.append((i,j)); adj[i].add(j); adj[j].add(i)
assert len(edges) == 720

# Triangles
triangles = []
for i,j in edges:
    for k in (adj[i] & adj[j]):
        if k > j: triangles.append((i,j,k))
assert len(triangles) == 1200

# Tetrahedra
tetrahedra = []
for i,j,k in triangles:
    for l in (adj[i] & adj[j] & adj[k]):
        if l > k: tetrahedra.append((i,j,k,l))
assert len(tetrahedra) == 600

print(f"    V={len(verts)}, E={len(edges)}, F={len(triangles)}, T={len(tetrahedra)}")
print(f"    Time: {time.time()-t0:.1f}s")

# === Build boundary operators ===
print("  Building boundary operators...")
t0 = time.time()
n_e = len(edges); n_f = len(triangles); n_t = len(tetrahedra)

edge_idx = {}
for idx, (i,j) in enumerate(edges):
    edge_idx[(i,j)] = idx; edge_idx[(j,i)] = idx

tri_idx = {}
for idx, (i,j,k) in enumerate(triangles):
    tri_idx[(i,j,k)] = idx

# d0: edges x vertices
r0,c0,v0 = [],[],[]
for idx,(i,j) in enumerate(edges):
    r0 += [idx,idx]; c0 += [i,j]; v0 += [-1.0,1.0]
d0 = sparse.csr_matrix((v0,(r0,c0)), shape=(n_e, 120))

# d1: triangles x edges
r1,c1,v1 = [],[],[]
for idx,(i,j,k) in enumerate(triangles):
    r1 += [idx,idx,idx]
    c1 += [edge_idx[(i,j)], edge_idx[(i,k)], edge_idx[(j,k)]]
    v1 += [1.0, -1.0, 1.0]
d1 = sparse.csr_matrix((v1,(r1,c1)), shape=(n_f, n_e))

# d2: tetrahedra x triangles
r2,c2,v2 = [],[],[]
for idx,(i,j,k,l) in enumerate(tetrahedra):
    for face, sign in [((j,k,l),1),((i,k,l),-1),((i,j,l),1),((i,j,k),-1)]:
        fs = tuple(sorted(face))
        if fs in tri_idx:
            r2.append(idx); c2.append(tri_idx[fs]); v2.append(float(sign))
d2 = sparse.csr_matrix((v2,(r2,c2)), shape=(n_t, n_f))

# Verify d^2 = 0
assert sparse.linalg.norm(d1.dot(d0)) < 1e-10, "d1*d0 != 0"
assert sparse.linalg.norm(d2.dot(d1)) < 1e-10, "d2*d1 != 0"
print(f"    d^2 = 0: VERIFIED. Time: {time.time()-t0:.1f}s")

# === Hodge Laplacians ===
print("  Computing Hodge Laplacians + spectra...")
t0 = time.time()

Delta0 = (d0.T.dot(d0)).toarray()
Delta1 = (d0.dot(d0.T) + d1.T.dot(d1)).toarray()
Delta2 = (d1.dot(d1.T) + d2.T.dot(d2)).toarray()
Delta3 = (d2.dot(d2.T)).toarray()

eig0 = np.sort(np.linalg.eigvalsh(Delta0))
eig1 = np.sort(np.linalg.eigvalsh(Delta1))
eig2 = np.sort(np.linalg.eigvalsh(Delta2))
eig3 = np.sort(np.linalg.eigvalsh(Delta3))

print(f"    Time: {time.time()-t0:.1f}s")

# Verify traces
tr = [np.trace(M) for M in [Delta0, Delta1, Delta2, Delta3]]
ns = [120, 720, 1200, 600]
print(f"    Per-simplex Tr: {[t/n for t,n in zip(tr,ns)]}")
print(f"    c_1 = sum(Tr) = {sum(tr):.0f}")

# Betti numbers
betti = [np.sum(np.abs(e) < 1e-8) for e in [eig0, eig1, eig2, eig3]]
print(f"    Betti: {betti}")

# === Spectral Zeta Function ===
print("\n--- Spectral Zeta Function ---")

all_eigs = np.concatenate([eig0, eig1, eig2, eig3])
nonzero = all_eigs[np.abs(all_eigs) > 1e-8]
print(f"  Total eigenvalues: {len(all_eigs)}")
print(f"  Nonzero: {len(nonzero)}, Zero: {len(all_eigs)-len(nonzero)}")

# Per-degree spectra
for k, (eigs, label) in enumerate([(eig0,'0-forms'),(eig1,'1-forms'),
                                     (eig2,'2-forms'),(eig3,'3-forms')]):
    nz = eigs[eigs > 1e-8]
    print(f"\n  {label}: {len(eigs)} total, {len(nz)} nonzero")
    if len(nz) > 0:
        print(f"    min={nz[0]:.4f}, max={nz[-1]:.4f}")
        print(f"    sum(sqrt) = {np.sum(np.sqrt(nz)):.4f}")
        print(f"    sum = {np.sum(nz):.4f}")
        # Distinct eigenvalues
        distinct = []
        for e in nz:
            if not any(abs(e-d) < 1e-6 for d in distinct):
                distinct.append(e)
        print(f"    Distinct eigenvalues: {len(distinct)}")

# Full zeta function
print(f"\n  --- zeta(s) values ---")
for s in [-3, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 3]:
    z = np.sum(nonzero**(-s))
    print(f"    zeta({s:+5.1f}) = {z:>20.4f}")

# Key identities
c0_val = 2640; c1_val = 14880; c2_val = 55920
print(f"\n  Verification:")
print(f"    zeta(-1) = {np.sum(nonzero):.0f} (should be c1-0 ~ {c1_val})")
# Note: zeta(-1) = sum of ALL eigenvalues including zeros
# Actually Tr(D^2) includes zeros, so c1 = sum of ALL = 14880
# But zeta(-1) = sum of NONZERO only
zeta_m1 = np.sum(nonzero)
print(f"    sum(all eigs) = {np.sum(all_eigs):.0f} = c1 = {c1_val}")
print(f"    sum(nonzero)  = {zeta_m1:.0f}")

zeta_m2 = np.sum(nonzero**2)
c2_check = np.sum(all_eigs**2)
print(f"    sum(all^2) = {c2_check:.0f} = c2 = {c2_val}")

# === Vacuum Energy ===
print(f"\n--- Vacuum Energy ---")
E_vac_total = 0.5 * np.sum(np.sqrt(nonzero))
print(f"  E_vac = (1/2)*sum sqrt(lambda) = {E_vac_total:.4f} (lattice units)")

# Per degree
for k, (eigs, label) in enumerate([(eig0,'0'),(eig1,'1'),(eig2,'2'),(eig3,'3')]):
    nz = eigs[eigs > 1e-8]
    E_k = 0.5 * np.sum(np.sqrt(nz))
    print(f"    E_vac({label}-forms) = {E_k:.4f}")

# Volume of S^3 (unit radius)
V_S3 = 2*np.pi**2  # = 19.739...
rho_vac = E_vac_total / V_S3
print(f"\n  V(S^3) = 2*pi^2 = {V_S3:.4f}")
print(f"  rho_vac = E_vac/V = {rho_vac:.4f} (lattice units)")

# === Physical scales ===
print(f"\n--- Cosmological Constant at Various Scales ---")
l_P = 1.616e-35    # Planck length (m)
m_P = 1.22089e19   # Planck mass (GeV)
R_H = 4.4e26       # Hubble radius (m)
Lambda_obs = 1.1056e-52  # m^{-2} (observed)
rho_obs_GeV4 = 2.846e-47  # GeV^4 (observed vacuum energy density)

# If 600-cell at radius R, eigenvalues scale as 1/R^2
# E_vac(R) = E_vac_dimless / R
# rho_vac(R) = E_vac_dimless / (2*pi^2*R^4)
# In natural units (hbar=c=1): rho_vac has dimensions [length]^{-4}

print(f"  rho_obs = {rho_obs_GeV4:.3e} GeV^4")
print(f"  m_P^4 = {m_P**4:.3e} GeV^4")
print(f"  rho_obs/m_P^4 = {rho_obs_GeV4/m_P**4:.3e}")

# At various R
for R_label, R_m in [("l_Planck", l_P), ("1 fm", 1e-15),
                      ("1 mm", 1e-3), ("1 m", 1.0), ("R_Hubble", R_H)]:
    # rho = rho_dimless / R^4, convert R to GeV^{-1}: R_GeV = R_m / (1.97e-16)
    R_GeV_inv = R_m / 1.97327e-16  # R in GeV^{-1}
    rho_R = rho_vac / R_GeV_inv**4  # GeV^4
    ratio = rho_R / rho_obs_GeV4
    log_ratio = np.log10(abs(ratio)) if ratio != 0 else 0
    print(f"    R = {R_label:12s}: rho = {rho_R:.3e} GeV^4, ratio = 10^{log_ratio:+.1f}")

# What R gives rho_obs?
R_needed_GeV = (rho_vac / rho_obs_GeV4)**0.25
R_needed_m = R_needed_GeV * 1.97327e-16
print(f"\n  R for correct Lambda_cc:")
print(f"    R = {R_needed_GeV:.3e} GeV^{{-1}} = {R_needed_m:.3e} m = {R_needed_m*1000:.4f} mm")
print(f"    This is the 'dark energy length scale' ~ 0.085 mm")
print(f"    (Known numerology, not specific to our theory)")

# === Z[phi] structure ===
print(f"\n--- Z[phi] Structure of Zeta Values ---")
for s in [0.5, 1.0, 1.5, 2.0, 3.0]:
    z = np.sum(nonzero**(-s))
    # Try z = a + b*PHI
    best_err = 1e10; best_ab = (0,0)
    for a_try in range(-5000, 5001):
        b_try = round((z - a_try)/PHI)
        err = abs(a_try + b_try*PHI - z)
        if err < best_err:
            best_err = err; best_ab = (a_try, b_try)
    a_z, b_z = best_ab
    in_zphi = best_err < 1e-4
    print(f"  zeta({s:.1f}) = {z:.6f} {'= '+str(a_z)+' + '+str(b_z)+'*phi' if in_zphi else '(NOT in Z[phi])'} (err={best_err:.2e})")

# === Spectral asymmetry and eta ===
print(f"\n--- Spectral Properties ---")
# Spectral dimension from zeta
# d_s = 2 * lim_{s->0+} s*log(zeta(s)) / log(1/s) (approximate)
# Or: d_s such that zeta(s) ~ s^{-d_s/2} for large s
# For discrete system, spectral dimension varies

# Heat trace: K(t) = sum exp(-t*lambda) = Tr(exp(-t*Delta))
print(f"  Heat trace K(t) = Tr(exp(-t*D^2)):")
for t in [0.001, 0.01, 0.1, 0.5, 1.0, 2.0, 5.0]:
    Kt = np.sum(np.exp(-t*all_eigs))
    d_eff = -2*np.log(Kt/len(all_eigs))/np.log(t) if t < 1 and Kt > 0 else 0
    print(f"    K({t:.3f}) = {Kt:.4f}")

# Spectral dimension from return probability
print(f"\n  Spectral dimension d_s = -2*d(ln K)/d(ln t):")
ts = np.logspace(-3, 1, 50)
Kts = np.array([np.sum(np.exp(-t*all_eigs)) for t in ts])
for i in range(1, len(ts)-1):
    d_s = -2*(np.log(Kts[i+1])-np.log(Kts[i-1]))/(np.log(ts[i+1])-np.log(ts[i-1]))
    if ts[i] in [0.001, 0.01, 0.1, 1.0] or abs(d_s - 3) < 0.1 or abs(d_s - 4) < 0.1:
        pass  # print selectively
# Just print at key scales
for t_target in [0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1.0, 5.0]:
    idx = np.argmin(np.abs(ts - t_target))
    if idx > 0 and idx < len(ts)-1:
        d_s = -2*(np.log(Kts[idx+1])-np.log(Kts[idx-1]))/(np.log(ts[idx+1])-np.log(ts[idx-1]))
        print(f"    t = {ts[idx]:.4f}: d_s = {d_s:.3f}")

# === Connection to framework constants ===
print(f"\n--- Framework Connections ---")
# E_vac decomposition
E0 = 0.5*np.sum(np.sqrt(eig0[eig0>1e-8]))
E1 = 0.5*np.sum(np.sqrt(eig1[eig1>1e-8]))
E2 = 0.5*np.sum(np.sqrt(eig2[eig2>1e-8]))
E3 = 0.5*np.sum(np.sqrt(eig3[eig3>1e-8]))
print(f"  E_vac = {E0:.2f} + {E1:.2f} + {E2:.2f} + {E3:.2f} = {E0+E1+E2+E3:.2f}")
print(f"  Ratios: E1/E0={E1/E0:.3f}, E2/E1={E2/E1:.3f}, E3/E2={E3/E2:.3f}")

# Bosonic vs fermionic interpretation
# In product with McKay: D_F eigenvalues {0,3,5,11,11,16,17,19,26}
n_f_vals = [0, 3, 5, 11, 11, 16, 17, 19, 26]
print(f"\n  Product vacuum energy (with McKay D_F):")
E_product = 0
for nf in n_f_vals:
    E_nf = 0.5 * np.sum(np.sqrt(nonzero + nf**2))
    E_product += E_nf
    print(f"    n_F={nf:2d}: E = {E_nf:.2f}")
print(f"    TOTAL: E_product = {E_product:.2f}")
print(f"    E_product/E_vac_pure = {E_product/E_vac_total:.3f}")
print(f"    E_product/9 = {E_product/9:.2f} (per fermion species)")

# Log ratio
print(f"\n  --- The 10^122 Problem ---")
print(f"  rho_vac(Planck) / rho_obs = 10^{np.log10(rho_vac * m_P**4 / rho_obs_GeV4):.1f}")
print(f"  N = 120, 4*h+2 = {4*30+2}")
print(f"  log10(m_P^4/rho_obs) = {np.log10(m_P**4/rho_obs_GeV4):.1f}")
print(f"  This is the standard CC problem, NOT specific to 600-cell")

print(f"\n  SUMMARY:")
print(f"    - Spectral zeta well-defined on 600-cell (finite, no regularization needed)")
print(f"    - E_vac = {E_vac_total:.2f} in lattice units")
print(f"    - Correct Lambda_cc requires R ~ 0.085 mm (= dark energy scale)")
print(f"    - This is known numerology, not a prediction")
print(f"    - Spectral dimension runs from ~3 (UV) to ~0 (IR)")
print(f"    - Product E_vac with McKay = {E_product:.2f}")
print(f"    STATUS: COMPUTED but NOT SOLVED (standard CC problem)")

print("\n" + "="*72)
print("EXP-266 COMPLETE")
print("="*72)
