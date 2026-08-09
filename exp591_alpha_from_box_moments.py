#!/usr/bin/env python3
"""
EXP-591: ECUATIA ALPHA DIN MOMENTELE SPECTRULUI BOX
=====================================================

IDEEA: Box = L_cross - a1*L_fiber a fost sursa pentru c^2=a1, foton, gravitatie.
Coeficientii ecuatiei alpha sunt toti din spectrul acelorasi Laplacieni.

INTREBARE: Apar coeficientii ecuatiei alpha ca RATII de Tr(Box^n)?
Daca da -> alpha iese din spectral action pe Box, nu din KK postulat.

Ecuatia: 2*pi*alpha^2 - 4*a1*phi^4*alpha + 1 = 0
  A = 2*pi          (topologic, c1=1 Hopf)
  B = 4*a1*phi^4     (spectral)

Tr(Box^n) pentru n=0,1,2,3,4,5 + ratii + cautare.
"""
import numpy as np
from numpy.linalg import eigh
import sys
sys.path.insert(0, ".")
from commons import build_600cell
from commons.constants import PHI, SQRT5, a1, b1, N, N_gen

ALPHA_EXP = 1/137.035999084

print("=" * 80)
print("EXP-591: ECUATIA ALPHA DIN MOMENTELE LUI BOX")
print("=" * 80)

# Build everything
verts, adj, lap = build_600cell()

# Build Hopf fibration
def qmul(p, q):
    return np.array([
        p[0]*q[0]-p[1]*q[1]-p[2]*q[2]-p[3]*q[3],
        p[0]*q[1]+p[1]*q[0]+p[2]*q[3]-p[3]*q[2],
        p[0]*q[2]-p[1]*q[3]+p[2]*q[0]+p[3]*q[1],
        p[0]*q[3]+p[1]*q[2]-p[2]*q[1]+p[3]*q[0]])

def find_idx(v, vs, tol=1e-6):
    dots = vs @ v; idx = np.argmax(dots)
    return idx if dots[idx] > 1 - tol else -1

fibers = None
for i in range(N):
    if abs(verts[i, 0] - PHI/2) < 1e-6:
        g = verts[i]; p = g.copy(); ok = True
        for k in range(2, 11):
            p = qmul(p, g)
            if k == 5 and not np.allclose(p, [-1,0,0,0], atol=1e-6): ok=False; break
            if k == 10 and not np.allclose(p, [1,0,0,0], atol=1e-6): ok=False
        if not ok: continue
        used = set(); fibers_tmp = []; subg = []
        pp = np.array([1.0,0,0,0])
        for k in range(10): subg.append(find_idx(pp, verts)); pp = qmul(pp, g)
        for s in range(N):
            if s in used: continue
            fib = []
            for si in subg:
                q = qmul(verts[s], verts[si]); idx = find_idx(q, verts)
                if idx >= 0 and idx not in used: fib.append(idx); used.add(idx)
            if len(fib) == 10: fibers_tmp.append(fib)
        if len(fibers_tmp) == 12:
            fibers = fibers_tmp
            break

A_fiber = np.zeros((N, N))
for fib in fibers:
    for i in fib:
        for j in fib:
            if i != j and adj[i,j] > 0.5: A_fiber[i,j] = 1.0

# Box = b1 * A_fiber - adj (wave operator)
Box = b1 * A_fiber - adj

# Also build L_fiber and L_cross
L_fiber = 2*np.eye(N) - A_fiber  # fiber Laplacian (each vertex has 2 fiber neighbors)
A_cross = adj - A_fiber
# Cross degree: each vertex has 10 cross neighbors
L_cross = 10*np.eye(N) - A_cross

# Verify Box = L_cross - a1*L_fiber
Box_check = L_cross - a1 * L_fiber
print(f"\n  Box = L_cross - a1*L_fiber? {np.allclose(Box, Box_check)}")

# =====================================================================
# PART 1: Trace moments of Box
# =====================================================================
print("\n" + "=" * 80)
print("PART 1: Tr(Box^n) MOMENTE")
print("=" * 80)

# Compute eigenvalues of Box
evals_box = np.sort(np.linalg.eigvalsh(Box))

# Compute trace moments
moments = {}
for n in range(8):
    tr = np.sum(evals_box**n)
    moments[n] = tr
    # Express as known quantities
    print(f"  Tr(Box^{n}) = {tr:>18.4f}", end="")
    if n == 0: print(f"  = N = {N}")
    elif n == 1: print(f"  (traceless)")
    elif n == 2: print(f"  = {tr/N:.1f} * N = {int(round(tr/N))} * {N}")
    elif n == 3: print(f"  = N^2 = {N**2}? {'YES' if abs(tr - N**2) < 1 else 'NO: ' + str(tr)}")
    else: print()

# Key moments
T0 = moments[0]  # = N = 120
T2 = moments[2]  # should be 7200
T3 = moments[3]  # should be N^2 = 14400 (trace-bootstrap!)
T4 = moments[4]
T5 = moments[5]

print(f"\n  Ratii cheie:")
print(f"    Tr(Box^2)/N = {T2/N:.4f}")
print(f"    Tr(Box^3)/N^2 = {T3/N**2:.4f}")
print(f"    Tr(Box^2)/N^2 = {T2/N**2:.4f}")
print(f"    Tr(Box^4)/Tr(Box^2) = {T4/T2:.4f}")
print(f"    Tr(Box^4)/N^2 = {T4/N**2:.4f}")
print(f"    Tr(Box^3)/Tr(Box^2) = {T3/T2:.4f}")

# =====================================================================
# PART 2: Cautam B = 4*a1*phi^4 in ratii de momente
# =====================================================================
print("\n" + "=" * 80)
print("PART 2: B = 4*a1*phi^4 IN RATII DE MOMENTE?")
print("=" * 80)

B = 4*a1*PHI**4  # = 137.082
A = 2*np.pi       # = 6.283

print(f"\n  Target B = 4*a1*phi^4 = {B:.6f}")
print(f"  Target A = 2*pi = {A:.6f}")
print(f"  B^2 = {B**2:.4f}")
print(f"  B^2 - 4*A = {B**2 - 4*A:.4f} (discriminant)")

# Try ratii de momente
print(f"\n  Ratii Tr(Box^m)/Tr(Box^n):")
for m in range(2, 8):
    for n in range(0, m):
        if abs(moments[n]) < 1e-10: continue
        ratio = moments[m] / moments[n]
        err_B = abs(ratio - B) / B * 100 if abs(ratio) > 1 else 999
        err_B2 = abs(ratio - B**2) / B**2 * 100 if abs(ratio) > 100 else 999
        if err_B < 5:
            print(f"    Tr(Box^{m})/Tr(Box^{n}) = {ratio:.4f}  ~B?  err = {err_B:.2f}%")
        if err_B2 < 5:
            print(f"    Tr(Box^{m})/Tr(Box^{n}) = {ratio:.4f}  ~B^2?  err = {err_B2:.2f}%")

# Try B = Tr(Box^m) / (N^k * something)
print(f"\n  Cautam B = Tr(Box^m) / (N^k * phi^p):")
for m in range(2, 7):
    for k in range(0, 4):
        for p in range(-6, 7):
            val = moments[m] / (N**k * PHI**p) if N**k * PHI**p != 0 else 0
            if abs(val) > 100 and abs(val) < 200:
                err = abs(val - B) / B * 100
                if err < 0.1:
                    print(f"    Tr(Box^{m})/(N^{k}*phi^{p}) = {val:.6f}  err = {err:.4f}%")

# =====================================================================
# PART 3: B^2 in momente (discriminant)
# =====================================================================
print("\n" + "=" * 80)
print("PART 3: B^2 SI DISCRIMINANTUL IN MOMENTE")
print("=" * 80)

B2 = B**2  # = 18791.5
disc = B2 - 4*A  # = 18766.3

print(f"  B^2 = {B2:.4f}")
print(f"  B^2 - 4A = {disc:.4f}")
print(f"  Tr(Box^2) = {T2:.4f}")
print(f"  Tr(Box^4) = {T4:.4f}")

# Key: is B^2 related to Tr(Box^n)/Tr(Box^m)?
print(f"\n  B^2 = {B2:.2f}. Candidati din momente:")
print(f"    Tr(Box^4)/Tr(Box^2) = {T4/T2:.4f}")
print(f"    Tr(Box^4)/N = {T4/N:.4f}")
print(f"    Tr(Box^4)/N^2 = {T4/N**2:.4f}")
print(f"    Tr(Box^2)^2/N^2 = {T2**2/N**2:.4f}")
print(f"    Tr(Box^2)*phi^4/N = {T2*PHI**4/N:.4f}")

# Check: B^2 = Tr(Box^2)*phi^4/N?
val = T2 * PHI**4 / N
print(f"\n  TEST: B^2 = Tr(Box^2)*phi^4/N = {T2}*{PHI**4:.4f}/{N} = {val:.4f}")
print(f"  B^2 = {B2:.4f}")
print(f"  Match: {abs(val - B2)/B2*100:.4f}%")

# Check: B = sqrt(Tr(Box^2)*phi^4/N)?
if val > 0:
    B_from_moments = np.sqrt(val)
    print(f"  B = sqrt(above) = {B_from_moments:.6f}")
    print(f"  B target = {B:.6f}")
    print(f"  Match: {abs(B_from_moments - B)/B*100:.4f}%")

# =====================================================================
# PART 4: Momente pe sector (fiber zero mode vs KK)
# =====================================================================
print("\n" + "=" * 80)
print("PART 4: MOMENTE PE SECTOARE (gauge vs KK)")
print("=" * 80)

# Project to fiber zero mode
P0 = np.zeros((N, N))
for fib in fibers:
    for i in fib:
        for j in fib:
            P0[i, j] = 1.0 / len(fib)

# Gauge sector Box
Box_gauge = P0 @ Box @ P0
evals_gauge = np.sort(np.linalg.eigvalsh(Box_gauge))
evals_gauge_nonzero = [e for e in evals_gauge if abs(e) > 0.01]

# KK sector Box
P_KK = np.eye(N) - P0
Box_KK = P_KK @ Box @ P_KK
evals_KK = np.sort(np.linalg.eigvalsh(Box_KK))
evals_KK_nonzero = [e for e in evals_KK if abs(e) > 0.01]

print(f"\n  Gauge sector (fiber zero mode): {len(evals_gauge_nonzero)} nonzero eigenvalues")
print(f"  KK sector (fiber nonzero modes): {len(evals_KK_nonzero)} nonzero eigenvalues")

# Moments per sector
for name, evals_sector in [("Gauge", evals_gauge), ("KK", evals_KK), ("Full", evals_box)]:
    print(f"\n  {name} sector moments:")
    for n in range(5):
        tr = np.sum(np.array(evals_sector)**n)
        print(f"    Tr(Box_{name}^{n}) = {tr:.4f}")

# =====================================================================
# PART 5: Spectral action ratii
# =====================================================================
print("\n" + "=" * 80)
print("PART 5: SPECTRAL ACTION RATII")
print("=" * 80)

# The spectral action coefficients (from paper)
c0 = 2640
c1 = 14880
c2 = 55920

print(f"\n  Spectral action coefficients: c0={c0}, c1={c1}, c2={c2}")
print(f"  c1/c0 = {c1/c0:.4f}")
print(f"  c2/c1 = {c2/c1:.4f}")
print(f"  c2/c0 = {c2/c0:.4f}")

# Known: G_N = N/c1 = 1/124
print(f"  G_N = N/c1 = {N}/{c1} = 1/{c1//N}")

# Is alpha related to c0, c1, c2?
print(f"\n  Cautam 1/alpha in ratii de c0, c1, c2:")
for expr, val in [
    ("c1/(c0/N)", c1/(c0/N)),
    ("c2/(c1/N)", c2/(c1/N)),
    ("c1/N", c1/N),
    ("c0/N", c0/N),
    ("c2/N", c2/N),
    ("c2/c1 * N/b1", c2/c1 * N/b1),
    ("c1^2/(c0*c2)", c1**2/(c0*c2)),
    ("c1*phi^4/c0", c1*PHI**4/c0),
    ("N*c1/(c0*b1)", N*c1/(c0*b1)),
]:
    err = abs(val - 1/ALPHA_EXP) / (1/ALPHA_EXP) * 100
    flag = "  <-- CLOSE!" if err < 1 else ""
    print(f"    {expr:>25s} = {val:>12.4f}  (err from 1/alpha: {err:.2f}%){flag}")

# =====================================================================
# PART 6: The gauge sector spectral action
# =====================================================================
print("\n" + "=" * 80)
print("PART 6: SPECTRAL ACTION PE SECTORUL GAUGE")
print("=" * 80)

# Tr(Box_gauge^2) should give the Yang-Mills term for U(1)
Tr_gauge_2 = sum(e**2 for e in evals_gauge)
Tr_gauge_4 = sum(e**4 for e in evals_gauge)
Tr_KK_2 = sum(e**2 for e in evals_KK)
Tr_full_2 = sum(e**2 for e in evals_box)

print(f"\n  Tr(Box_gauge^2) = {Tr_gauge_2:.4f}")
print(f"  Tr(Box_KK^2) = {Tr_KK_2:.4f}")
print(f"  Tr(Box_full^2) = {Tr_full_2:.4f}")
print(f"  Tr_gauge/Tr_full = {Tr_gauge_2/Tr_full_2:.6f} = 1/{Tr_full_2/Tr_gauge_2:.4f}")

# Ratio Tr_full/Tr_gauge = a1?
print(f"\n  Tr(Box_full^2)/Tr(Box_gauge^2) = {Tr_full_2/Tr_gauge_2:.4f}")
print(f"  a1 = {a1}")
print(f"  N/n_base = {N/len(fibers)}")

# Can we get B from gauge sector?
print(f"\n  B din gauge sector:")
print(f"    Tr_gauge_2 * phi^4 / n_base = {Tr_gauge_2 * PHI**4 / len(fibers):.4f}")
print(f"    B = {B:.4f}")
print(f"    Tr_gauge_2 / n_base = {Tr_gauge_2/len(fibers):.4f}")
print(f"    Tr_gauge_2 / N = {Tr_gauge_2/N:.4f}")
print(f"    Tr_gauge_4 / Tr_gauge_2 = {Tr_gauge_4/Tr_gauge_2:.4f}")

# =====================================================================
# PART 7: THE KEY TEST — alpha from Tr(Box^n)?
# =====================================================================
print("\n" + "=" * 80)
print("PART 7: TEST CHEIE — alpha din Tr(Box^n)")
print("=" * 80)

# The alpha equation: A*alpha^2 - B*alpha + 1 = 0
# Can we express A and B as spectral data of Box?

# We know:
# B = (N/b1)*phi^4 = 20*phi^4
# A = 2*pi

# From the traces:
# Tr(Box^2) = 7200
# N = 120, b1 = 6

# B = (N/b1) * phi^4 = Tr(Box^0)/b1 * phi^4 = (120/6)*phi^4 = 20*phi^4
# This is: B = [Tr(Box^0)/b1] * [1/lambda_1(fiber)]^2
# where b1 and lambda_1 are from the FIBER part of Box

# What about A = 2*pi?
# 2*pi = holonomy of Hopf fiber = 2*a1 * (pi/a1) = 10 * pi/5
# Is there a spectral way to get 2*pi from Box?

# Heat kernel: Tr(exp(-t*Box)) at specific t?
print(f"\n  Heat kernel Tr(exp(-t*Box)):")
for t in [0.01, 0.1, 0.5, 1.0, 2.0, 5.0]:
    hk = sum(np.exp(-t*e) for e in evals_box)
    print(f"    t={t:.2f}: Tr(exp(-t*Box)) = {hk:.4f}")

# Spectral zeta: sum(lambda^{-s})
print(f"\n  Spectral zeta zeta_Box(s) = sum(|lambda|^(-s)) [nonzero only]:")
evals_nonzero = [e for e in evals_box if abs(e) > 0.01]
for s in [0.5, 1.0, 1.5, 2.0]:
    zeta = sum(abs(e)**(-s) for e in evals_nonzero)
    print(f"    zeta_Box({s}) = {zeta:.6f}")

# Zeta determinant: exp(-zeta'(0))
# This is related to the regularized determinant

# =====================================================================
# PART 8: Tr(Box^2) decomposition
# =====================================================================
print("\n" + "=" * 80)
print("PART 8: DECOMPOZITIEA Tr(Box^2)")
print("=" * 80)

# Box = L_cross - a1*L_fiber
# Box^2 = L_cross^2 - 2*a1*L_cross*L_fiber + a1^2*L_fiber^2
# Tr(Box^2) = Tr(L_cross^2) - 2*a1*Tr(L_cross*L_fiber) + a1^2*Tr(L_fiber^2)

# Since [L_fiber, L_cross] = 0, the cross term is real

Tr_Lc2 = np.trace(L_cross @ L_cross)
Tr_Lf2 = np.trace(L_fiber @ L_fiber)
Tr_LcLf = np.trace(L_cross @ L_fiber)

print(f"  Tr(L_cross^2) = {Tr_Lc2:.4f}")
print(f"  Tr(L_fiber^2) = {Tr_Lf2:.4f}")
print(f"  Tr(L_cross*L_fiber) = {Tr_LcLf:.4f}")
print(f"  Tr(Box^2) = {Tr_Lc2} - 2*{a1}*{Tr_LcLf} + {a1}^2*{Tr_Lf2}")
print(f"           = {Tr_Lc2} - {2*a1*Tr_LcLf} + {a1**2*Tr_Lf2}")
print(f"           = {Tr_Lc2 - 2*a1*Tr_LcLf + a1**2*Tr_Lf2:.4f}")
print(f"  Direct:    {T2:.4f}")

# Eigenvalues of L_fiber and L_cross
evals_Lf = np.sort(np.linalg.eigvalsh(L_fiber))
evals_Lc = np.sort(np.linalg.eigvalsh(L_cross))
print(f"\n  L_fiber: min={evals_Lf[0]:.4f}, max={evals_Lf[-1]:.4f}, Tr={sum(evals_Lf):.4f}")
print(f"  L_cross: min={evals_Lc[0]:.4f}, max={evals_Lc[-1]:.4f}, Tr={sum(evals_Lc):.4f}")

# The spectral gap of L_fiber = 1/phi^2
gap_fiber = sorted(set(np.round(evals_Lf, 6)))[1]  # smallest nonzero
print(f"\n  Fiber spectral gap: {gap_fiber:.6f} = 1/phi^2 = {1/PHI**2:.6f}")

# =====================================================================
# PART 9: FORMULA CANDIDAT
# =====================================================================
print("\n" + "=" * 80)
print("PART 9: FORMULA CANDIDAT PENTRU ALPHA DIN BOX")
print("=" * 80)

# The observation from Part 3:
# B^2 = Tr(Box^2) * phi^4 / N?
val_B2 = T2 * PHI**4 / N
print(f"\n  Test: B^2 = Tr(Box^2)*phi^4/N")
print(f"    LHS = B^2 = {B**2:.6f}")
print(f"    RHS = {T2}*{PHI**4:.6f}/{N} = {val_B2:.6f}")
print(f"    Ratio: {B**2/val_B2:.6f}")

# Actually: B = (N/b1)*phi^4 so B^2 = N^2*phi^8/b1^2
# And Tr(Box^2) = 7200 = N*60 = N*|A5| = N*(N/2)
# So: Tr(Box^2)*phi^4/N = (N^2/2)*phi^4/N = N*phi^4/2
# And: B^2 = N^2*phi^8/b1^2
# Ratio: B^2/(Tr*phi^4/N) = N*phi^4/b1^2 / (1/2) = 2*N*phi^4/b1^2

# Let me try: B = Tr(Box^2)/(N*something)
for denom_name, denom in [
    ("b1*phi^2", b1*PHI**2),
    ("N/phi^2", N/PHI**2),
    ("Tr(Box^3)/Tr(Box^2)", T3/T2),
    ("2*a1*phi^2", 2*a1*PHI**2),
    ("degree", 12.0),
    ("a1*(a1-1)", a1*(a1-1)),
    ("a1*phi^2", a1*PHI**2),
    ("b1*sqrt(a1)", b1*np.sqrt(a1)),
]:
    val = T2 / (N * denom) if N*denom != 0 else 0
    if 100 < val < 200:
        err = abs(val - B)/B*100
        print(f"    Tr(Box^2)/(N*{denom_name}) = {T2}/(120*{denom:.4f}) = {val:.4f}  err={err:.3f}%")

# =====================================================================
# PART 10: VERDICT
# =====================================================================
print("\n" + "=" * 80)
print("V E R D I C T")
print("=" * 80)

# Let me check the most promising leads
print(f"""
  MOMENTE CALCULATE:
  Tr(Box^0) = {moments[0]:.0f} = N
  Tr(Box^1) = {moments[1]:.4f} (traceless)
  Tr(Box^2) = {moments[2]:.0f} = N*(N/2) = N*|A5|
  Tr(Box^3) = {moments[3]:.0f} = N^2 (trace-bootstrap!)
  Tr(Box^4) = {moments[4]:.0f}
  Tr(Box^5) = {moments[5]:.0f}

  RATII RELEVANTE:
  Tr(Box^2)/N = {T2/N:.0f} = N/2 = |A5| = 60
  Tr(Box^3)/N^2 = {T3/N**2:.4f}
  Tr(Box^4)/Tr(Box^2) = {T4/T2:.4f}
  Tr(Box^2)/b1^2 = {T2/b1**2:.4f} = N*(N/2)/b1^2

  TREE-LEVEL B:
  B = (N/b1)*phi^4 = {B:.6f}
  B in termeni spectral: N/(b1*lambda_1^2) unde lambda_1 = 1/phi^2
  Echivalent: Tr(Box^0) / (b1 * [spectral gap fiber]^2)

  OBSERVATIE CHEIE:
  Tr(Box^2) = {T2:.0f} = {int(T2/N)}*N
  {int(T2/N)} = N/2 = a1! / 2 = |A5|

  B = (N/b1)*phi^4 = {B:.4f}
  B = [Tr(Box^0)/b1] * phi^4

  Coeficientul B NU apare ca raport de momente Tr(Box^m)/Tr(Box^n).
  Apare ca: [moment de ordin 0] / [parametru structural b1] * [spectral gap]^{{-2}}

  Asta inseamna ca alpha ecuatia combina:
  - Date TOPOLOGICE (2*pi din c1=1)
  - Date SPECTRALE de ordin 0 (N = Tr(Box^0))
  - Date STRUCTURALE (b1 = fiber stiffness)
  - Date SPECTRALE ale FIBREI (gap = 1/phi^2)

  NU iese dintr-un singur Tr(Box^n).
  Dar TOATE ingredientele sunt din structura spectrala a 600-cell-ului.
""")
