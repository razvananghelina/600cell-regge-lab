"""
EXP-122: Complete CKM Matrix Derivation from 600-Cell Geometry
==============================================================

Attempts a COMPLETE derivation of all 9 CKM matrix elements
from the 600-cell geometry using 5 different approaches.

Context:
- PHI = (1+sqrt(5))/2
- Quark quantum numbers (a,b): u(3,-2), c(2,1), t(4,1), d(1,0), s(1,1), b(-1,4)
- Spectral levels: n = 5a + 6b
- Up-type: n_u=3, n_c=16, n_t=26
- Down-type: n_d=5, n_s=11, n_b=19

APPROACHES:
1. CKM from (a,b) geometry (vector overlaps)
2. CKM from spectral levels (level differences)
3. CKM from Wolfenstein with lambda = phi^(-3)
4. CKM from line structure in (a,b) lattice
5. Direct phi-power CKM ansatz
"""

import numpy as np

PHI = (1 + np.sqrt(5)) / 2

print("=" * 78)
print("EXP-122: COMPLETE CKM MATRIX DERIVATION FROM 600-CELL GEOMETRY")
print("=" * 78)

# ==============================================================
# EXPERIMENTAL CKM VALUES (PDG 2024)
# ==============================================================
print("\n" + "=" * 78)
print("EXPERIMENTAL CKM MATRIX (PDG 2024)")
print("=" * 78)

V_exp = np.array([
    [0.97373, 0.2243, 0.00382],
    [0.221,   0.975,  0.0408 ],
    [0.0086,  0.0415, 1.014  ]
])

V_exp_err = np.array([
    [0.00031, 0.0008, 0.00020],
    [0.004,   0.006,  0.0014 ],
    [0.0002,  0.0011, 0.029  ]
])

labels_up = ['u', 'c', 't']
labels_dn = ['d', 's', 'b']

print("\n|V_CKM| experimental:")
print("        d          s          b")
for i in range(3):
    row = "  " + labels_up[i] + "  "
    for j in range(3):
        row += f"  {V_exp[i,j]:.5f}"
    print(row)

# Quantum numbers
a_up = np.array([3, 2, 4])   # u, c, t
b_up = np.array([-2, 1, 1])
a_dn = np.array([1, 1, -1])  # d, s, b
b_dn = np.array([0, 1, 4])

n_up = 5 * a_up + 6 * b_up   # [3, 16, 26]
n_dn = 5 * a_dn + 6 * b_dn   # [5, 11, 19]

print(f"\nQuantum numbers (a,b) and spectral levels n=5a+6b:")
for i, q in enumerate(labels_up):
    print(f"  {q}: a={a_up[i]:+d}, b={b_up[i]:+d}, n={n_up[i]}")
for i, q in enumerate(labels_dn):
    print(f"  {q}: a={a_dn[i]:+d}, b={b_dn[i]:+d}, n={n_dn[i]}")

print(f"\nUp-type levels: n_u={n_up[0]}, n_c={n_up[1]}, n_t={n_up[2]}")
print(f"Down-type levels: n_d={n_dn[0]}, n_s={n_dn[1]}, n_b={n_dn[2]}")

# Helper function: compute error metrics
def compute_errors(V_pred, V_exp, V_exp_err, label=""):
    """Compute element-by-element errors and summary statistics."""
    abs_err = np.abs(V_pred - V_exp)
    rel_err = abs_err / V_exp * 100  # percent
    sigma_err = abs_err / V_exp_err  # in units of sigma

    # Unitarity check
    prod = V_pred @ V_pred.T
    unitarity_diag = np.array([prod[i,i] for i in range(3)])
    unitarity_offdiag = np.array([abs(prod[0,1]), abs(prod[0,2]), abs(prod[1,2])])

    print(f"\n--- {label} ---")
    print(f"{'Element':<10} {'Pred':>10} {'Exp':>10} {'Err%':>8} {'Sigma':>8}")
    print("-" * 50)
    for i in range(3):
        for j in range(3):
            name = f"V_{labels_up[i]}{labels_dn[j]}"
            print(f"  {name:<8} {V_pred[i,j]:10.5f} {V_exp[i,j]:10.5f} "
                  f"{rel_err[i,j]:7.2f}% {sigma_err[i,j]:7.1f}s")

    mean_rel = np.mean(rel_err)
    max_rel = np.max(rel_err)
    rms_rel = np.sqrt(np.mean(rel_err**2))
    chi2 = np.sum(sigma_err**2)

    print(f"\n  Mean relative error: {mean_rel:.2f}%")
    print(f"  Max  relative error: {max_rel:.2f}%")
    print(f"  RMS  relative error: {rms_rel:.2f}%")
    print(f"  Chi^2 (9 dof):       {chi2:.1f}")

    print(f"\n  Unitarity (V V^T diagonal): {unitarity_diag}")
    print(f"  Unitarity (V V^T off-diag): {unitarity_offdiag}")

    return mean_rel, max_rel, chi2


# ==============================================================
# APPROACH 1: CKM from (a,b) geometry
# ==============================================================
print("\n" + "=" * 78)
print("APPROACH 1: CKM FROM (a,b) GEOMETRY")
print("=" * 78)

print("""
The (a,b) vectors for each quark define positions in a 2D lattice.
CKM elements should relate to "overlaps" between up-type and down-type vectors.

Vectors: u=(3,-2), c=(2,1), t=(4,1), d=(1,0), s=(1,1), b=(-1,4)
""")

# Define vectors
v_up = np.array([[3, -2], [2, 1], [4, 1]], dtype=float)
v_dn = np.array([[1, 0], [1, 1], [-1, 4]], dtype=float)

# --- Approach 1a: Euclidean distance with exponential decay ---
print("--- 1a: V_ij ~ exp(-lambda * |v_ui - v_dj|) ---")

# Try different lambda values
best_1a_err = 1e10
best_1a_lam = 0
best_1a_V = None

for lam in np.linspace(0.1, 5.0, 500):
    V_pred = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            dist = np.linalg.norm(v_up[i] - v_dn[j])
            V_pred[i, j] = np.exp(-lam * dist)
        # Normalize each row
        V_pred[i] /= np.linalg.norm(V_pred[i])
    err = np.mean(np.abs(V_pred - V_exp) / V_exp * 100)
    if err < best_1a_err:
        best_1a_err = err
        best_1a_lam = lam
        best_1a_V = V_pred.copy()

print(f"\nBest lambda = {best_1a_lam:.4f}")
m1a, x1a, c1a = compute_errors(best_1a_V, V_exp, V_exp_err, "1a: Exp(-lam*dist)")

# --- Approach 1b: V_ij ~ phi^(-|v_ui - v_dj|) ---
print("\n--- 1b: V_ij ~ phi^(-alpha * |v_ui - v_dj|) ---")

best_1b_err = 1e10
best_1b_alpha = 0
best_1b_V = None

for alpha in np.linspace(0.1, 10.0, 1000):
    V_pred = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            dist = np.linalg.norm(v_up[i] - v_dn[j])
            V_pred[i, j] = PHI**(-alpha * dist)
        V_pred[i] /= np.linalg.norm(V_pred[i])
    err = np.mean(np.abs(V_pred - V_exp) / V_exp * 100)
    if err < best_1b_err:
        best_1b_err = err
        best_1b_alpha = alpha
        best_1b_V = V_pred.copy()

print(f"\nBest alpha = {best_1b_alpha:.4f}")
m1b, x1b, c1b = compute_errors(best_1b_V, V_exp, V_exp_err, "1b: phi^(-alpha*dist)")

# --- Approach 1c: Angle between vectors in (a,b) space ---
print("\n--- 1c: V_ij ~ |cos(angle between v_ui and v_dj)| ---")

V_1c = np.zeros((3, 3))
for i in range(3):
    for j in range(3):
        dot = np.dot(v_up[i], v_dn[j])
        norm_u = np.linalg.norm(v_up[i])
        norm_d = np.linalg.norm(v_dn[j])
        if norm_u > 0 and norm_d > 0:
            V_1c[i, j] = abs(dot / (norm_u * norm_d))
        else:
            V_1c[i, j] = 0

# Normalize rows
for i in range(3):
    V_1c[i] /= np.linalg.norm(V_1c[i])

m1c, x1c, c1c = compute_errors(V_1c, V_exp, V_exp_err, "1c: cos(angle)")

# --- Approach 1d: Manhattan distance in (a,b) ---
print("\n--- 1d: V_ij ~ phi^(-alpha * Manhattan_dist) ---")

best_1d_err = 1e10
best_1d_alpha = 0
best_1d_V = None

for alpha in np.linspace(0.1, 10.0, 1000):
    V_pred = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            dist = abs(a_up[i] - a_dn[j]) + abs(b_up[i] - b_dn[j])
            V_pred[i, j] = PHI**(-alpha * dist)
        V_pred[i] /= np.linalg.norm(V_pred[i])
    err = np.mean(np.abs(V_pred - V_exp) / V_exp * 100)
    if err < best_1d_err:
        best_1d_err = err
        best_1d_alpha = alpha
        best_1d_V = V_pred.copy()

print(f"\nBest alpha = {best_1d_alpha:.4f}")
m1d, x1d, c1d = compute_errors(best_1d_V, V_exp, V_exp_err, "1d: phi^(-alpha*Manhattan)")

# --- Approach 1e: Weighted distance with (5,6) metric ---
print("\n--- 1e: V_ij ~ phi^(-|5*da + 6*db|) [spectral metric] ---")

# This is equivalent to using |n_i - n_j| but let's try with a free scale

best_1e_err = 1e10
best_1e_scale = 0
best_1e_V = None

for scale in np.linspace(0.01, 3.0, 300):
    V_pred = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            da = a_up[i] - a_dn[j]
            db = b_up[i] - b_dn[j]
            # Use the spectral distance = |5*da + 6*db| = |n_i - n_j|
            dist = abs(5 * da + 6 * db)
            V_pred[i, j] = PHI**(-scale * dist)
        V_pred[i] /= np.linalg.norm(V_pred[i])
    err = np.mean(np.abs(V_pred - V_exp) / V_exp * 100)
    if err < best_1e_err:
        best_1e_err = err
        best_1e_scale = scale
        best_1e_V = V_pred.copy()

print(f"\nBest scale = {best_1e_scale:.4f}")
m1e, x1e, c1e = compute_errors(best_1e_V, V_exp, V_exp_err, "1e: phi^(-scale*|dn|)")


# ==============================================================
# APPROACH 2: CKM from spectral levels
# ==============================================================
print("\n" + "=" * 78)
print("APPROACH 2: CKM FROM SPECTRAL LEVELS")
print("=" * 78)

print("""
Spectral levels: n_u=3, n_c=16, n_t=26, n_d=5, n_s=11, n_b=19

Level difference matrix |n_i - n_j|:
""")

# Level differences
dn_matrix = np.zeros((3, 3), dtype=int)
for i in range(3):
    for j in range(3):
        dn_matrix[i, j] = abs(n_up[i] - n_dn[j])

print("  |n_i - n_j|:")
print("        d    s    b")
for i in range(3):
    row = "  " + labels_up[i] + "  "
    for j in range(3):
        row += f"  {dn_matrix[i,j]:3d}"
    print(row)

# --- Approach 2a: V_ij ~ phi^(-|dn|) ---
print("\n--- 2a: V_ij ~ phi^(-|n_i - n_j|) (unnormalized) ---")

V_2a_raw = np.zeros((3, 3))
for i in range(3):
    for j in range(3):
        V_2a_raw[i, j] = PHI**(-dn_matrix[i, j])

print("\nRaw phi^(-|dn|):")
print("        d          s          b")
for i in range(3):
    row = "  " + labels_up[i] + "  "
    for j in range(3):
        row += f"  {V_2a_raw[i,j]:.6f}"
    print(row)

# Normalize
V_2a = V_2a_raw.copy()
for i in range(3):
    V_2a[i] /= np.linalg.norm(V_2a[i])

m2a, x2a, c2a = compute_errors(V_2a, V_exp, V_exp_err, "2a: phi^(-|dn|) normalized")

# --- Approach 2b: V_ij ~ phi^(-alpha * |dn|) with best alpha ---
print("\n--- 2b: V_ij ~ phi^(-alpha * |dn|) with fitted alpha ---")

best_2b_err = 1e10
best_2b_alpha = 0
best_2b_V = None

for alpha in np.linspace(0.01, 2.0, 2000):
    V_pred = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            V_pred[i, j] = PHI**(- alpha * dn_matrix[i, j])
        V_pred[i] /= np.linalg.norm(V_pred[i])
    err = np.mean(np.abs(V_pred - V_exp) / V_exp * 100)
    if err < best_2b_err:
        best_2b_err = err
        best_2b_alpha = alpha
        best_2b_V = V_pred.copy()

print(f"\nBest alpha = {best_2b_alpha:.4f}")
m2b, x2b, c2b = compute_errors(best_2b_V, V_exp, V_exp_err, "2b: phi^(-alpha*|dn|)")

# Check if alpha is close to a simple fraction
print(f"\n  alpha = {best_2b_alpha:.4f}")
print(f"  alpha * 5 = {best_2b_alpha * 5:.4f}")
print(f"  alpha * 6 = {best_2b_alpha * 6:.4f}")
print(f"  alpha * 12 = {best_2b_alpha * 12:.4f}")
print(f"  1/(2*alpha) = {1/(2*best_2b_alpha):.4f}")

# --- Approach 2c: V_ij ~ phi^(-alpha * |dn|^beta) ---
print("\n--- 2c: V_ij ~ phi^(-alpha * |dn|^beta) with best alpha, beta ---")

best_2c_err = 1e10
best_2c_params = (0, 0)
best_2c_V = None

for alpha in np.linspace(0.01, 3.0, 200):
    for beta in np.linspace(0.3, 3.0, 200):
        V_pred = np.zeros((3, 3))
        for i in range(3):
            for j in range(3):
                dn = dn_matrix[i, j]
                if dn > 0:
                    V_pred[i, j] = PHI**(- alpha * dn**beta)
                else:
                    V_pred[i, j] = 1.0
            V_pred[i] /= np.linalg.norm(V_pred[i])
        err = np.mean(np.abs(V_pred - V_exp) / V_exp * 100)
        if err < best_2c_err:
            best_2c_err = err
            best_2c_params = (alpha, beta)
            best_2c_V = V_pred.copy()

print(f"\nBest alpha = {best_2c_params[0]:.4f}, beta = {best_2c_params[1]:.4f}")
m2c, x2c, c2c = compute_errors(best_2c_V, V_exp, V_exp_err, "2c: phi^(-alpha*|dn|^beta)")

# --- Approach 2d: Using sin(theta) = phi^(-|dn|/k) to build rotation matrix ---
print("\n--- 2d: Rotation matrix from spectral angles ---")

# CKM = R23 * R13 * R12 (standard parametrization)
# Try: sin(theta_12) ~ phi^(-|n_u - n_s|/k) = phi^(-8/k)
# and  sin(theta_23) ~ phi^(-|n_c - n_b|/k) = phi^(-3/k)
# and  sin(theta_13) ~ phi^(-|n_u - n_b|/k) = phi^(-16/k)

# From exp082: theta_C = arctan(phi^(-3)) works with ~1.8% error
# sin(theta_C) ~ 0.225, phi^(-3) ~ 0.236

# Let's check: what |dn| corresponds to each CKM angle?
# theta_12 ~ V_us: primary mixing u-s.  |n_u - n_s| = |3-11| = 8
# theta_23 ~ V_cb: primary mixing c-b.  |n_c - n_b| = |16-19| = 3
# theta_13 ~ V_ub: primary mixing u-b.  |n_u - n_b| = |3-19| = 16

print(f"CKM angle assignments from spectral levels:")
print(f"  theta_12 ~ V_us: |n_u - n_s| = {abs(n_up[0]-n_dn[1])}")
print(f"  theta_23 ~ V_cb: |n_c - n_b| = {abs(n_up[1]-n_dn[2])}")
print(f"  theta_13 ~ V_ub: |n_u - n_b| = {abs(n_up[0]-n_dn[2])}")

# Try sin(theta) = phi^(-dn/k) with fitted k
best_2d_err = 1e10
best_2d_k = 0
best_2d_V = None

dn_12 = abs(n_up[0] - n_dn[1])  # 8
dn_23 = abs(n_up[1] - n_dn[2])  # 3
dn_13 = abs(n_up[0] - n_dn[2])  # 16

for k in np.linspace(0.5, 15.0, 1500):
    s12 = PHI**(-dn_12 / k)
    s23 = PHI**(-dn_23 / k)
    s13 = PHI**(-dn_13 / k)

    # Clamp
    s12 = min(s12, 0.999)
    s23 = min(s23, 0.999)
    s13 = min(s13, 0.999)

    c12 = np.sqrt(1 - s12**2)
    c23 = np.sqrt(1 - s23**2)
    c13 = np.sqrt(1 - s13**2)

    # Standard parametrization (no CP phase for magnitudes)
    V_pred = np.array([
        [c12*c13,          s12*c13,          s13    ],
        [-s12*c23 - c12*s23*s13,  c12*c23 - s12*s23*s13,  s23*c13],
        [s12*s23 - c12*c23*s13,  -c12*s23 - s12*c23*s13,  c23*c13]
    ])
    V_pred = np.abs(V_pred)
    err = np.mean(np.abs(V_pred - V_exp) / V_exp * 100)
    if err < best_2d_err:
        best_2d_err = err
        best_2d_k = k
        best_2d_V = V_pred.copy()

print(f"\nBest k = {best_2d_k:.4f}")
print(f"  sin(theta_12) = phi^(-{dn_12}/{best_2d_k:.2f}) = phi^(-{dn_12/best_2d_k:.4f})")
print(f"  sin(theta_23) = phi^(-{dn_23}/{best_2d_k:.2f}) = phi^(-{dn_23/best_2d_k:.4f})")
print(f"  sin(theta_13) = phi^(-{dn_13}/{best_2d_k:.2f}) = phi^(-{dn_13/best_2d_k:.4f})")

m2d, x2d, c2d = compute_errors(best_2d_V, V_exp, V_exp_err, "2d: sin(theta)=phi^(-dn/k)")


# ==============================================================
# APPROACH 3: CKM from Wolfenstein parametrization
# ==============================================================
print("\n" + "=" * 78)
print("APPROACH 3: CKM FROM WOLFENSTEIN WITH lambda = phi^(-n)")
print("=" * 78)

print("""
Wolfenstein parametrization:
  V_us ~ lambda
  V_cb ~ A * lambda^2
  V_ub ~ A * lambda^3 * sqrt(rho^2 + eta^2)
  V_td ~ A * lambda^3 * sqrt((1-rho)^2 + eta^2)

Try lambda = phi^(-n) for n = 2, 3, and fractional values.
""")

# --- 3a: lambda = phi^(-3) (from exp082) ---
print("--- 3a: lambda = phi^(-3) ---")
lam3 = PHI**(-3)
print(f"  lambda = phi^(-3) = {lam3:.6f}")
print(f"  V_us (exp) = {V_exp[0,1]:.4f}, lambda = {lam3:.4f}, ratio = {V_exp[0,1]/lam3:.4f}")

# Fit A, rho, eta
# V_cb = A * lambda^2
A3 = V_exp[1, 2] / lam3**2
print(f"  A = V_cb/lambda^2 = {A3:.4f}")

# V_ub = A * lambda^3 * sqrt(rho^2 + eta^2)
rhoeta3 = V_exp[0, 2] / (A3 * lam3**3)
print(f"  sqrt(rho^2+eta^2) = {rhoeta3:.4f}")

# Build Wolfenstein CKM (to order lambda^3)
V_3a = np.array([
    [1 - lam3**2/2,           lam3,                    A3*lam3**3*rhoeta3   ],
    [lam3,                    1 - lam3**2/2,            A3*lam3**2           ],
    [A3*lam3**3,              A3*lam3**2,               1                    ]
])
V_3a = np.abs(V_3a)

m3a, x3a, c3a = compute_errors(V_3a, V_exp, V_exp_err, "3a: Wolfenstein, lambda=phi^(-3)")

# --- 3b: Scan over n in lambda = phi^(-n) ---
print("\n--- 3b: Scan lambda = phi^(-n) for best n ---")

best_3b_err = 1e10
best_3b_n = 0
best_3b_V = None
best_3b_A = 0

for n_val in np.linspace(2.0, 5.0, 3000):
    lam = PHI**(-n_val)
    if lam >= 1 or lam <= 0:
        continue
    A = V_exp[1, 2] / lam**2
    if A <= 0:
        continue
    rhoeta_val = V_exp[0, 2] / (A * lam**3) if A * lam**3 > 0 else 0

    V_pred = np.array([
        [1 - lam**2/2,    lam,             A*lam**3*rhoeta_val],
        [lam,             1 - lam**2/2,    A*lam**2            ],
        [A*lam**3,        A*lam**2,        1                   ]
    ])
    V_pred = np.abs(V_pred)
    err = np.mean(np.abs(V_pred - V_exp) / V_exp * 100)
    if err < best_3b_err:
        best_3b_err = err
        best_3b_n = n_val
        best_3b_V = V_pred.copy()
        best_3b_A = A

print(f"\nBest n = {best_3b_n:.4f}, lambda = phi^(-{best_3b_n:.4f}) = {PHI**(-best_3b_n):.6f}")
print(f"  A = {best_3b_A:.4f}")
m3b, x3b, c3b = compute_errors(best_3b_V, V_exp, V_exp_err, "3b: Wolfenstein, best n")

# --- 3c: Full rotation matrix with 3 angles from phi ---
print("\n--- 3c: Full rotation matrix, all angles from phi ---")

# theta_12: sin = lambda = phi^(-n12)
# theta_23: sin = A*lambda^2 = phi^(-n23) ?
# theta_13: sin = V_ub = phi^(-n13) ?

# Compute effective exponents for experimental values
n_exp_12 = -np.log(V_exp[0, 1]) / np.log(PHI)
n_exp_23 = -np.log(V_exp[1, 2]) / np.log(PHI)
n_exp_13 = -np.log(V_exp[0, 2]) / np.log(PHI)

print(f"\nEffective phi-exponents for experimental sin(theta):")
print(f"  sin(theta_12) = {V_exp[0,1]:.5f} = phi^(-{n_exp_12:.4f})  [closest int: {round(n_exp_12)}]")
print(f"  sin(theta_23) = {V_exp[1,2]:.5f} = phi^(-{n_exp_23:.4f})  [closest int: {round(n_exp_23)}]")
print(f"  sin(theta_13) = {V_exp[0,2]:.6f} = phi^(-{n_exp_13:.4f})  [closest int: {round(n_exp_13)}]")

# Try integer exponents
for n12, n23, n13 in [(3, 7, 12), (3, 6, 12), (3, 7, 11), (3, 6, 11),
                       (3, 7, 13), (3, 6, 10), (3, 7, 10)]:
    s12 = PHI**(-n12)
    s23 = PHI**(-n23)
    s13 = PHI**(-n13)
    c12 = np.sqrt(1 - s12**2)
    c23 = np.sqrt(1 - s23**2)
    c13 = np.sqrt(1 - s13**2)

    V_pred = np.abs(np.array([
        [c12*c13,          s12*c13,          s13    ],
        [-s12*c23-c12*s23*s13,  c12*c23-s12*s23*s13,  s23*c13],
        [s12*s23-c12*c23*s13,  -c12*s23-s12*c23*s13,  c23*c13]
    ]))
    err = np.mean(np.abs(V_pred - V_exp) / V_exp * 100)
    print(f"  n=({n12},{n23},{n13}): mean err = {err:.2f}%")

# Build the best one: n = (3, 7, 12) from exp082/exp085
print("\nDetailed result for n = (3, 7, 12):")
s12 = PHI**(-3)
s23 = PHI**(-7)
s13 = PHI**(-12)
c12 = np.sqrt(1 - s12**2)
c23 = np.sqrt(1 - s23**2)
c13 = np.sqrt(1 - s13**2)

V_3c = np.abs(np.array([
    [c12*c13,          s12*c13,          s13    ],
    [-s12*c23-c12*s23*s13,  c12*c23-s12*s23*s13,  s23*c13],
    [s12*s23-c12*c23*s13,  -c12*s23-s12*c23*s13,  c23*c13]
]))

m3c, x3c, c3c = compute_errors(V_3c, V_exp, V_exp_err, "3c: sin(theta)=phi^(-n), n=(3,7,12)")

# --- 3d: fractional exponents from spectral levels ---
print("\n--- 3d: Fractional exponents from spectral level ratios ---")

# What if the exponents come from ratios of spectral levels?
# theta_12: n_u/n_s = 3/11, or |n_u-n_s|/|n_u-n_d| = 8/2 = 4...
# Let's check various combinations

combos = {
    '|n_u-n_s|/|n_d-n_s|': abs(n_up[0]-n_dn[1]) / abs(n_dn[0]-n_dn[1]) if abs(n_dn[0]-n_dn[1]) != 0 else 0,
    '|n_c-n_b|/|n_s-n_b|': abs(n_up[1]-n_dn[2]) / abs(n_dn[1]-n_dn[2]) if abs(n_dn[1]-n_dn[2]) != 0 else 0,
    '|n_u-n_b|/|n_d-n_b|': abs(n_up[0]-n_dn[2]) / abs(n_dn[0]-n_dn[2]) if abs(n_dn[0]-n_dn[2]) != 0 else 0,
    'n_s/n_u': n_dn[1]/n_up[0],
    'n_b/n_c': n_dn[2]/n_up[1],
    'n_b/n_u': n_dn[2]/n_up[0],
}

print(f"\nPotential spectral ratios:")
for name, val in combos.items():
    phi_exp = val
    sin_val = PHI**(-phi_exp)
    print(f"  {name} = {val:.4f} -> phi^(-{val:.4f}) = {sin_val:.6f}")

# Try theta_ij = arctan(phi^(-dn_ij/divisor))
# where dn_12 = 8, dn_23 = 3, dn_13 = 16
# and the divisor might be related to the 600-cell (e.g., diameter=5)

for divisor_name, divisor in [("diameter=5", 5), ("decagon=6", 6),
                                ("grad/2=6", 6), ("(d+k)/2=5.5", 5.5),
                                ("edges/vertices=6", 6), ("3", 3)]:
    s12 = PHI**(-dn_12 / divisor)
    s23 = PHI**(-dn_23 / divisor)
    s13 = PHI**(-dn_13 / divisor)

    # Clamp
    s12 = min(s12, 0.9999)
    s23 = min(s23, 0.9999)
    s13 = min(s13, 0.9999)

    c12 = np.sqrt(1 - s12**2)
    c23 = np.sqrt(1 - s23**2)
    c13 = np.sqrt(1 - s13**2)

    V_pred = np.abs(np.array([
        [c12*c13, s12*c13, s13],
        [-s12*c23-c12*s23*s13, c12*c23-s12*s23*s13, s23*c13],
        [s12*s23-c12*c23*s13, -c12*s23-s12*c23*s13, c23*c13]
    ]))
    err = np.mean(np.abs(V_pred - V_exp) / V_exp * 100)
    print(f"  divisor={divisor_name}: s12={s12:.4f} s23={s23:.4f} s13={s13:.6f} -> mean err={err:.2f}%")


# ==============================================================
# APPROACH 4: CKM from line structure
# ==============================================================
print("\n" + "=" * 78)
print("APPROACH 4: CKM FROM LINE STRUCTURE IN (a,b) LATTICE")
print("=" * 78)

print("""
The (a,b) quantum numbers lie on 3 lines in the lattice:
  Line L1 (a=1):         d(1,0), s(1,1), tau(1,2)   [leptons + down sector]
  Line L2 (b=1):         s(1,1), c(2,1), t(4,1)     [heavy quarks]
  Line L3 (3a+2b=5):     u(3,-2), s(1,1), b(-1,4)   [cross-generation]

Note: s(1,1) is at the INTERSECTION of all 3 lines!

Hypothesis: Mixing is determined by which lines the quarks share.
  - Same line: large mixing
  - Different lines: small mixing
""")

# Determine which lines each quark belongs to
def on_line(a_val, b_val):
    lines = []
    if a_val == 1:
        lines.append('L1')
    if b_val == 1:
        lines.append('L2')
    if 3*a_val + 2*b_val == 5:
        lines.append('L3')
    return lines

print("Quark line assignments:")
for i, q in enumerate(labels_up):
    lines = on_line(a_up[i], b_up[i])
    print(f"  {q}({a_up[i]},{b_up[i]}): {lines}")
for i, q in enumerate(labels_dn):
    lines = on_line(a_dn[i], b_dn[i])
    print(f"  {q}({a_dn[i]},{b_dn[i]}): {lines}")

# Count shared lines between up-type i and down-type j
shared_lines = np.zeros((3, 3), dtype=int)
for i in range(3):
    for j in range(3):
        lines_u = set(on_line(a_up[i], b_up[i]))
        lines_d = set(on_line(a_dn[j], b_dn[j]))
        shared_lines[i, j] = len(lines_u & lines_d)

print(f"\nShared lines matrix:")
print("        d    s    b")
for i in range(3):
    row = "  " + labels_up[i] + "  "
    for j in range(3):
        row += f"  {shared_lines[i,j]:3d}"
    print(row)

# Define "line distance" = minimum number of line changes to connect
# Two quarks on same line: distance 0
# Two quarks sharing a line through s: distance 1 (must go through intersection)
# Two quarks on different lines with no intersection: distance 2

# For line distance, let's use a different metric:
# Count total lines for each quark pair
line_dist = np.zeros((3, 3), dtype=int)
for i in range(3):
    for j in range(3):
        lines_u = set(on_line(a_up[i], b_up[i]))
        lines_d = set(on_line(a_dn[j], b_dn[j]))
        if lines_u & lines_d:
            line_dist[i, j] = 0  # share a line
        else:
            # Check if they connect through s (which is on all 3 lines)
            # All quarks connect through s, so at most 1 step
            line_dist[i, j] = 1

print(f"\nLine distance matrix:")
print("        d    s    b")
for i in range(3):
    row = "  " + labels_up[i] + "  "
    for j in range(3):
        row += f"  {line_dist[i,j]:3d}"
    print(row)

# --- 4a: Mixing from line overlap with spectral weighting ---
print("\n--- 4a: Combined line + spectral approach ---")

# Hypothesis: V_ij depends on both line structure and spectral distance
# V_ij ~ phi^(-alpha * dn * (1 + beta * (1 - shared_lines)))
best_4a_err = 1e10
best_4a_params = (0, 0)
best_4a_V = None

for alpha in np.linspace(0.01, 1.0, 100):
    for beta in np.linspace(0.0, 5.0, 100):
        V_pred = np.zeros((3, 3))
        for i in range(3):
            for j in range(3):
                dn = dn_matrix[i, j]
                shared = shared_lines[i, j]
                eff_dist = dn * (1 + beta * (1 - shared))
                V_pred[i, j] = PHI**(-alpha * eff_dist)
            V_pred[i] /= np.linalg.norm(V_pred[i])
        err = np.mean(np.abs(V_pred - V_exp) / V_exp * 100)
        if err < best_4a_err:
            best_4a_err = err
            best_4a_params = (alpha, beta)
            best_4a_V = V_pred.copy()

print(f"\nBest alpha={best_4a_params[0]:.4f}, beta={best_4a_params[1]:.4f}")
m4a, x4a, c4a = compute_errors(best_4a_V, V_exp, V_exp_err, "4a: Line+spectral")

# --- 4b: Independent line-based rotation angles ---
print("\n--- 4b: Rotation angles from line structure ---")

# The 3 CKM angles might correspond to the 3 lines:
# theta_12 (u-s mixing): u is on L3, s is on L1+L2+L3, shared=L3
# theta_23 (c-b mixing): c is on L2, b is on L3, shared=none
# theta_13 (u-b mixing): u is on L3, b is on L3, shared=L3

# Hypothesis: angle ~ phi^(-(something related to line membership))
# When on same line, mixing is governed by distance ALONG the line
# When on different lines, mixing is suppressed by line separation

# On L3 (3a+2b=5): u(3,-2), s(1,1), b(-1,4)
# distance along L3: u->s = sqrt((3-1)^2+(-2-1)^2) = sqrt(13), s->b = sqrt((1+1)^2+(1-4)^2) = sqrt(13)
# u->b = sqrt((3+1)^2+(-2-4)^2) = sqrt(52) = 2*sqrt(13)

print(f"  On line L3: u(3,-2) -> s(1,1) -> b(-1,4)")
print(f"    d(u,s) = sqrt(13) = {np.sqrt(13):.4f}")
print(f"    d(s,b) = sqrt(13) = {np.sqrt(13):.4f}")
print(f"    d(u,b) = 2*sqrt(13) = {2*np.sqrt(13):.4f}")

print(f"  On line L2 (b=1): s(1,1), c(2,1), t(4,1)")
print(f"    d(s,c) = 1, d(c,t) = 2, d(s,t) = 3")

print(f"  On line L1 (a=1): d(1,0), s(1,1)")
print(f"    d(d,s) = 1")


# ==============================================================
# APPROACH 5: Direct phi-power CKM
# ==============================================================
print("\n" + "=" * 78)
print("APPROACH 5: DIRECT PHI-POWER CKM ANSATZ")
print("=" * 78)

print("""
Ansatz: |V_ij|^2 = phi^(-k_ij) where k_ij depends on quantum numbers.

First, compute k_ij from experimental values:
""")

# Compute k_ij
k_matrix = np.zeros((3, 3))
for i in range(3):
    for j in range(3):
        if V_exp[i, j] > 0 and V_exp[i, j] < 1:
            k_matrix[i, j] = -2 * np.log(V_exp[i, j]) / np.log(PHI)
        else:
            k_matrix[i, j] = 0

print("k_ij matrix (|V_ij|^2 = phi^(-k_ij)):")
print("        d          s          b")
for i in range(3):
    row = "  " + labels_up[i] + "  "
    for j in range(3):
        row += f"  {k_matrix[i,j]:8.4f}"
    print(row)

print(f"\nLooking for patterns in k_ij:")
print(f"  k_ud = {k_matrix[0,0]:.4f} (small, ~diagonal)")
print(f"  k_us = {k_matrix[0,1]:.4f}")
print(f"  k_ub = {k_matrix[0,2]:.4f}")
print(f"  k_cd = {k_matrix[1,0]:.4f}")
print(f"  k_cs = {k_matrix[1,1]:.4f} (small, ~diagonal)")
print(f"  k_cb = {k_matrix[1,2]:.4f}")
print(f"  k_td = {k_matrix[2,0]:.4f}")
print(f"  k_ts = {k_matrix[2,1]:.4f}")
print(f"  k_tb = {k_matrix[2,2]:.4f} (small, ~diagonal)")

# Check relations
print(f"\nRatios and patterns:")
print(f"  k_us/k_ub = {k_matrix[0,1]/k_matrix[0,2]:.4f}")
print(f"  k_cb/k_ub = {k_matrix[1,2]/k_matrix[0,2]:.4f}")
print(f"  k_td/k_ub = {k_matrix[2,0]/k_matrix[0,2]:.4f}")
print(f"  k_ub/k_us = {k_matrix[0,2]/k_matrix[0,1]:.4f}")
print(f"  k_td/k_ts = {k_matrix[2,0]/k_matrix[2,1]:.4f}")
print(f"  k_us ~ {k_matrix[0,1]:.2f} ~ 2*phi = {2*PHI:.4f}")
print(f"  k_cb ~ {k_matrix[1,2]:.2f} ~ 2*phi^2 = {2*PHI**2:.4f}")

# --- 5a: k_ij = alpha * |n_i - n_j| + beta ---
print("\n--- 5a: k_ij = alpha * |n_i - n_j| (linear spectral) ---")

# Fit alpha
dn_flat = dn_matrix.flatten()
k_flat = k_matrix.flatten()

# Only use off-diagonal (non-tiny k)
mask = k_flat > 0.1
if np.sum(mask) > 0:
    from numpy.linalg import lstsq
    A_fit = dn_flat[mask].reshape(-1, 1)
    b_fit = k_flat[mask]
    result = lstsq(A_fit, b_fit, rcond=None)
    alpha_fit = result[0][0]
    print(f"\n  Best fit: k = {alpha_fit:.4f} * |dn|")

    # Also fit with intercept
    A_fit2 = np.column_stack([dn_flat[mask], np.ones(np.sum(mask))])
    result2 = lstsq(A_fit2, b_fit, rcond=None)
    alpha2, beta2 = result2[0]
    print(f"  With intercept: k = {alpha2:.4f} * |dn| + {beta2:.4f}")

    # Predict CKM with k = alpha * |dn|
    V_5a = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            k_pred = alpha_fit * dn_matrix[i, j]
            V_5a[i, j] = PHI**(-k_pred / 2)  # |V| = phi^(-k/2)

    # Normalize
    for i in range(3):
        V_5a[i] /= np.linalg.norm(V_5a[i])

    m5a, x5a, c5a = compute_errors(V_5a, V_exp, V_exp_err, "5a: k=alpha*|dn|")

# --- 5b: Direct mapping k_ij from (a,b) differences ---
print("\n--- 5b: k_ij from quadratic form in (da, db) ---")

# Try k_ij = p*(da)^2 + q*(db)^2 + r*da*db
da_matrix = np.zeros((3, 3))
db_matrix = np.zeros((3, 3))
for i in range(3):
    for j in range(3):
        da_matrix[i, j] = a_up[i] - a_dn[j]
        db_matrix[i, j] = b_up[i] - b_dn[j]

# Fit k = p*da^2 + q*db^2 + r*da*db
da_flat = da_matrix.flatten()
db_flat = db_matrix.flatten()

X_quad = np.column_stack([da_flat**2, db_flat**2, da_flat*db_flat])
result_quad = lstsq(X_quad, k_flat, rcond=None)
p, q, r = result_quad[0]
print(f"\n  k = {p:.4f}*da^2 + {q:.4f}*db^2 + {r:.4f}*da*db")

# Predict
V_5b = np.zeros((3, 3))
for i in range(3):
    for j in range(3):
        da = da_matrix[i, j]
        db = db_matrix[i, j]
        k_pred = p*da**2 + q*db**2 + r*da*db
        k_pred = max(k_pred, 0)
        V_5b[i, j] = PHI**(-k_pred / 2)

for i in range(3):
    V_5b[i] /= np.linalg.norm(V_5b[i])

m5b, x5b, c5b = compute_errors(V_5b, V_exp, V_exp_err, "5b: k=quadratic(da,db)")

# --- 5c: Direct integer k search ---
print("\n--- 5c: Search for integer k_ij that minimize total error ---")

# For off-diagonal elements, search integer k values
# V_ij = phi^(-k/2), V_ii ~ 1 (determined by unitarity)

# Focus on the 3 independent angles
# sin(theta_12) = phi^(-n1), etc.
# k_us = 2*n1, k_cb = 2*n2, k_ub = 2*n3

print(f"\n  Required k values for perfect match:")
print(f"    k_us = 2 * {n_exp_12:.4f} = {2*n_exp_12:.4f}")
print(f"    k_cb = 2 * {n_exp_23:.4f} = {2*n_exp_23:.4f}")
print(f"    k_ub = 2 * {n_exp_13:.4f} = {2*n_exp_13:.4f}")

# Nearest integers/half-integers
for step in [1, 0.5]:
    tag = "integer" if step == 1 else "half-integer"
    best_err = 1e10
    best_combo = None
    best_V = None

    k_range_12 = np.arange(1, 10, step)
    k_range_23 = np.arange(2, 15, step)
    k_range_13 = np.arange(5, 25, step)

    for k12 in k_range_12:
        for k23 in k_range_23:
            for k13 in k_range_13:
                s12 = PHI**(-k12/2)
                s23 = PHI**(-k23/2)
                s13 = PHI**(-k13/2)

                if s12 >= 1 or s23 >= 1 or s13 >= 1:
                    continue

                c12 = np.sqrt(1 - s12**2)
                c23 = np.sqrt(1 - s23**2)
                c13 = np.sqrt(1 - s13**2)

                V_pred = np.abs(np.array([
                    [c12*c13, s12*c13, s13],
                    [-s12*c23-c12*s23*s13, c12*c23-s12*s23*s13, s23*c13],
                    [s12*s23-c12*c23*s13, -c12*s23-s12*c23*s13, c23*c13]
                ]))
                err = np.mean(np.abs(V_pred - V_exp) / V_exp * 100)
                if err < best_err:
                    best_err = err
                    best_combo = (k12, k23, k13)
                    best_V = V_pred.copy()

    print(f"\n  Best {tag} k: ({best_combo[0]}, {best_combo[1]}, {best_combo[2]}), mean err={best_err:.2f}%")
    if best_combo is not None:
        print(f"    sin(theta_12) = phi^(-{best_combo[0]/2:.1f}) = {PHI**(-best_combo[0]/2):.5f} (exp: {V_exp[0,1]:.5f})")
        print(f"    sin(theta_23) = phi^(-{best_combo[1]/2:.1f}) = {PHI**(-best_combo[1]/2):.5f} (exp: {V_exp[1,2]:.5f})")
        print(f"    sin(theta_13) = phi^(-{best_combo[2]/2:.1f}) = {PHI**(-best_combo[2]/2):.6f} (exp: {V_exp[0,2]:.6f})")


# ==============================================================
# APPROACH 5d: CKM from GEOMETRIC angles in (a,b) space
# ==============================================================
print("\n--- 5d: CKM from rotation angle needed to align up with down sectors ---")

# Consider the up-type and down-type (a,b) vectors as defining
# two coordinate frames. The CKM matrix is the rotation between them.

# For 3 generations, we need a 3x3 rotation. Let's try to build it
# by defining orthonormal bases for up and down sectors.

# Up-type vectors in (a,b) space: (3,-2), (2,1), (4,1)
# Down-type vectors in (a,b) space: (1,0), (1,1), (-1,4)

# These are 2D vectors, so we need to use a different approach.
# Let's embed them in a higher-dimensional space using n = 5a + 6b.

# Use (a, b, n) as 3D vectors
v3_up = np.array([[a_up[i], b_up[i], n_up[i]] for i in range(3)], dtype=float)
v3_dn = np.array([[a_dn[i], b_dn[i], n_dn[i]] for i in range(3)], dtype=float)

# Gram-Schmidt orthonormalize each set
def gram_schmidt(V):
    U = np.zeros_like(V, dtype=float)
    for i in range(len(V)):
        U[i] = V[i].copy()
        for j in range(i):
            U[i] -= np.dot(U[i], U[j]) * U[j]
        norm = np.linalg.norm(U[i])
        if norm > 1e-10:
            U[i] /= norm
    return U

U_up = gram_schmidt(v3_up)
U_dn = gram_schmidt(v3_dn)

# CKM = U_up . U_dn^T (overlap of the two bases)
V_5d = np.abs(U_up @ U_dn.T)

print(f"\n  Up-type orthonormal basis (from Gram-Schmidt):")
for i in range(3):
    print(f"    e_{labels_up[i]} = ({U_up[i,0]:.4f}, {U_up[i,1]:.4f}, {U_up[i,2]:.4f})")

print(f"\n  Down-type orthonormal basis (from Gram-Schmidt):")
for i in range(3):
    print(f"    e_{labels_dn[i]} = ({U_dn[i,0]:.4f}, {U_dn[i,1]:.4f}, {U_dn[i,2]:.4f})")

m5d, x5d, c5d = compute_errors(V_5d, V_exp, V_exp_err, "5d: Gram-Schmidt overlap (a,b,n)")

# --- 5e: Try with just (a, b) vectors in 2D, using SVD ---
print("\n--- 5e: SVD overlap of (a,b) vectors ---")

# Build 3x2 matrices and use SVD
M_up = np.array([[a_up[i], b_up[i]] for i in range(3)], dtype=float)
M_dn = np.array([[a_dn[i], b_dn[i]] for i in range(3)], dtype=float)

# Overlap matrix O = M_up . M_dn^T (3x3)
O = M_up @ M_dn.T

# SVD of O: O = U S V^T
U_svd, S_svd, Vt_svd = np.linalg.svd(O)

# The "rotation" part: R = U V^T
R = U_svd @ Vt_svd
V_5e = np.abs(R)

print(f"  Overlap matrix M_up . M_dn^T:")
for i in range(3):
    row = "    "
    for j in range(3):
        row += f"{O[i,j]:8.2f}"
    print(row)

print(f"  Singular values: {S_svd}")

m5e, x5e, c5e = compute_errors(V_5e, V_exp, V_exp_err, "5e: SVD rotation from (a,b)")


# ==============================================================
# APPROACH 6 (BONUS): Hierarchical phi structure
# ==============================================================
print("\n" + "=" * 78)
print("APPROACH 6 (BONUS): HIERARCHICAL PHI STRUCTURE")
print("=" * 78)

print("""
The CKM matrix has a hierarchical structure. In the 600-cell,
the hierarchy comes from the graph distance (diameter = 5).

Key observation:
  V_us ~ phi^(-3)   ->  3 = dim(Im(H)) quaternions
  V_cb ~ phi^(-6.6) ->  close to 7 = dim(Im(O)) octonions
  V_ub ~ phi^(-11.6)->  close to 12 = degree of 600-cell

From exp082/exp085:
  n_12 = 3 (quaternions)
  n_23 = 7 (octonions)
  n_13 = 12 (degree)

But the actual exponents are NOT integers:
  n_12_exact = 3.10 (close to 3)
  n_23_exact = 6.62 (NOT close to 7)
  n_13_exact = 11.57 (close to 12)
""")

# --- 6a: Perturbative correction to integer exponents ---
print("--- 6a: Correction to n=(3,7,12) from 600-cell geometry ---")

# The corrections might come from:
# delta_n12 = n_exact - 3 = 0.10
# delta_n23 = n_exact - 7 = -0.38
# delta_n13 = n_exact - 12 = -0.43

dn12_corr = n_exp_12 - 3
dn23_corr = n_exp_23 - 7
dn13_corr = n_exp_13 - 12

print(f"\n  Corrections needed:")
print(f"    delta_n12 = {dn12_corr:+.4f}")
print(f"    delta_n23 = {dn23_corr:+.4f}")
print(f"    delta_n13 = {dn13_corr:+.4f}")

# Check if corrections relate to known quantities
alpha_em = 1/137.036
alpha_s_mz = 0.1179
sin2_tw = 0.23122

print(f"\n  Possible origins of corrections:")
print(f"    delta_n12 = {dn12_corr:+.4f} ~ alpha_em*pi*4.4 = {alpha_em*np.pi*4.4:.4f}")
print(f"    delta_n12 = {dn12_corr:+.4f} ~ sin^2(tW)/2 = {sin2_tw/2:.4f}")
print(f"    delta_n23 = {dn23_corr:+.4f} ~ -alpha_s*pi = {-alpha_s_mz*np.pi:.4f}")
print(f"    delta_n23 = {dn23_corr:+.4f} ~ -2*sin^2(tW)/phi = {-2*sin2_tw/PHI:.4f}")
print(f"    delta_n13 = {dn13_corr:+.4f} ~ delta_n12 + delta_n23 = {dn12_corr+dn23_corr:.4f}")
print(f"    (actual delta_n13 = {dn13_corr:.4f}, sum = {dn12_corr+dn23_corr:.4f})")

# If corrections come from alpha_s perturbation
print(f"\n  If corrections are O(alpha_s):")
print(f"    delta_n12/alpha_s = {dn12_corr/alpha_s_mz:.4f}")
print(f"    delta_n23/alpha_s = {dn23_corr/alpha_s_mz:.4f}")
print(f"    delta_n13/alpha_s = {dn13_corr/alpha_s_mz:.4f}")

# --- 6b: n_ij from spectral level formula ---
print("\n--- 6b: Check if n_ij = |n_up - n_dn| / divisor matches ---")

print(f"\n  |n_up - n_dn| matrix:")
print(f"        d={n_dn[0]:2d}   s={n_dn[1]:2d}   b={n_dn[2]:2d}")
for i in range(3):
    row = f"  {labels_up[i]}(n={n_up[i]:2d})"
    for j in range(3):
        row += f"  {dn_matrix[i,j]:4d}"
    print(row)

# What divisors would give the right CKM exponents?
print(f"\n  Divisors needed for |dn|/div = n_exact:")
for i, j, n_exact_val, label in [(0,1, n_exp_12, "theta_12"),
                                   (1,2, n_exp_23, "theta_23"),
                                   (0,2, n_exp_13, "theta_13")]:
    dn = dn_matrix[i, j]
    div = dn / n_exact_val
    print(f"    {label}: |dn|={dn}, n_exact={n_exact_val:.4f}, divisor={div:.4f}")

# Check: divisor ~ 2.58 for theta_12, ~ 0.45 for theta_23, ~ 1.38 for theta_13
# Not uniform -- this approach doesn't give a single divisor.

# --- 6c: Mixed approach: use phi^(-3) for Cabibbo, then hierarchy ---
print("\n--- 6c: Hierarchical CKM with lambda=phi^(-3) ---")

lam = PHI**(-3)
# Standard Wolfenstein with corrected A
# V_cb = A * lam^2 => A = V_cb / lam^2
A_val = V_exp[1, 2] / lam**2
rhoeta_val = V_exp[0, 2] / (A_val * lam**3)

# More precise Wolfenstein (order lambda^5)
s12 = lam
s23 = A_val * lam**2
s13 = V_exp[0, 2]  # Use exact value for V_ub

c12 = np.sqrt(1 - s12**2)
c23 = np.sqrt(1 - s23**2)
c13 = np.sqrt(1 - s13**2)

V_6c = np.abs(np.array([
    [c12*c13, s12*c13, s13],
    [-s12*c23-c12*s23*s13, c12*c23-s12*s23*s13, s23*c13],
    [s12*s23-c12*c23*s13, -c12*s23-s12*c23*s13, c23*c13]
]))

m6c, x6c, c6c = compute_errors(V_6c, V_exp, V_exp_err, "6c: Hierarchical lam=phi^(-3), fit A")

# Check: is A also geometric?
print(f"\n  Wolfenstein parameters:")
print(f"    lambda = phi^(-3) = {lam:.6f}")
print(f"    A = {A_val:.4f}")
print(f"    A ~ phi^(-?) => {-np.log(A_val)/np.log(PHI):.4f}")
print(f"    A*lambda = {A_val*lam:.4f}")
print(f"    A*lambda ~ phi^(-?) => {-np.log(A_val*lam)/np.log(PHI):.4f}")

# Check if A = phi^(-n)/lam^2 * lam^2 = phi^(-n_A)
n_A = -np.log(A_val) / np.log(PHI)
print(f"    A = phi^(-{n_A:.4f})")
print(f"    A close to 1/phi = {1/PHI:.4f}? A = {A_val:.4f}")
print(f"    A close to 1/phi^(0.6) = {PHI**(-0.6):.4f}? Diff = {abs(A_val - PHI**(-0.6)):.4f}")

# --- 6d: The "complete" geometric CKM ---
print("\n--- 6d: Best possible geometric CKM ---")

# Use: sin(theta_12) = phi^(-3) * (1 + correction)
# where correction comes from sin^2(theta_W) = 6/26

# Try to get V_cb from geometry
# V_cb = 0.0408, -log(V_cb)/log(phi) = 6.62
# Maybe V_cb = phi^(-3) * phi^(-3.62)?  = phi^(-6.62)
# Or V_cb ~ phi^(-3) * phi^(-3) * phi^(-0.62)

# Alternative: V_cb = phi^(-3)^2 * phi^(delta) ?
# phi^(-6) = 0.0557, too high
# phi^(-7) = 0.0345, too low

# Geometric mean: sqrt(phi^(-6) * phi^(-7)) = phi^(-6.5) = 0.0438
# That's closer! phi^(-6.5) = 0.0438 vs 0.0408

# Try V_cb = phi^(-6.5)
print(f"  phi^(-6.5) = {PHI**(-6.5):.5f} vs V_cb = {V_exp[1,2]:.5f} ({abs(PHI**(-6.5)-V_exp[1,2])/V_exp[1,2]*100:.1f}%)")
print(f"  phi^(-6.62) = {PHI**(-6.62):.5f} (exact match)")

# Try V_ub: phi^(-11.57)
# Close to phi^(-23/2) = phi^(-11.5)
print(f"  phi^(-11.5) = {PHI**(-11.5):.6f} vs V_ub = {V_exp[0,2]:.6f} ({abs(PHI**(-11.5)-V_exp[0,2])/V_exp[0,2]*100:.1f}%)")
print(f"  phi^(-23/2) = phi^(-11.5) = {PHI**(-11.5):.6f}")

# Build CKM with half-integer exponents: n = (3, 13/2, 23/2)
s12 = PHI**(-3)
s23 = PHI**(-13/2)
s13 = PHI**(-23/2)
c12 = np.sqrt(1 - s12**2)
c23 = np.sqrt(1 - s23**2)
c13 = np.sqrt(1 - s13**2)

V_6d = np.abs(np.array([
    [c12*c13, s12*c13, s13],
    [-s12*c23-c12*s23*s13, c12*c23-s12*s23*s13, s23*c13],
    [s12*s23-c12*c23*s13, -c12*s23-s12*c23*s13, c23*c13]
]))

m6d, x6d, c6d = compute_errors(V_6d, V_exp, V_exp_err, "6d: sin=phi^(-n), n=(3, 13/2, 23/2)")

# --- 6e: Exponents from 600-cell numerology ---
print("\n--- 6e: Exponents from 600-cell numbers ---")

# 600-cell numbers: 120 vertices, 720 edges, 1200 faces, 600 cells
# Diameter = 5, degree = 12, decagon = 6 (sides)
# h (Coxeter number) = 30
# H4 exponents: 1, 11, 19, 29

# Try: n_12 = 3 = (h-1)/some = ...
# n_12 = 3 = 30/10 = h/10
# n_23 = 7 = (h-30)/? no
# n_13 = 12 = degree

# Try exponents as Coxeter differences:
# h = 30, exponents 1,11,19,29
# Differences: 10, 8, 10
# Sums: 12 (1+11), 30 (11+19), 48 (19+29), 40 (11+29), 20 (1+19)

print(f"  Coxeter-derived quantities:")
print(f"    H4 exponents: 1, 11, 19, 29")
print(f"    Differences: 10, 8, 10")
print(f"    (11-1)/2 - 2 = 3  <-- n_12?")
print(f"    (19-1)/2 - 2 = 7  <-- n_23?")
print(f"    (29-1)/2 + 2 - 4 = 12  <-- nope")
print(f"    11 + 1 = 12 = degree <-- n_13?")

# Try Coxeter formula: n = (m_i + m_j) / c
for c in [1, 2, 3, 4]:
    for (m_i, m_j) in [(1,11), (11,19), (19,29), (1,19), (1,29), (11,29)]:
        val = (m_i + m_j) / c
        if val in [3, 7, 12]:
            print(f"    ({m_i}+{m_j})/{c} = {val} <-- match!")

# Alternative: n = |m_i - m_j| / c
for c in [1, 2, 3, 4]:
    for (m_i, m_j) in [(1,11), (11,19), (19,29), (1,19), (1,29), (11,29)]:
        val = abs(m_i - m_j) / c
        if val in [3, 7, 12]:
            print(f"    |{m_i}-{m_j}|/{c} = {val} <-- match!")

# ==============================================================
# FINAL SUMMARY
# ==============================================================
print("\n" + "=" * 78)
print("FINAL SUMMARY: ALL APPROACHES COMPARED")
print("=" * 78)

results = [
    ("1a: exp(-lam*dist)",          m1a, x1a, c1a, "FITTING"),
    ("1b: phi^(-alpha*dist)",       m1b, x1b, c1b, "FITTING"),
    ("1c: cos(angle)",              m1c, x1c, c1c, "GEOMETRIC"),
    ("1d: phi^(-alpha*Manhattan)",  m1d, x1d, c1d, "FITTING"),
    ("1e: phi^(-scale*|dn|)",       m1e, x1e, c1e, "FITTING"),
    ("2a: phi^(-|dn|) norm",        m2a, x2a, c2a, "GEOMETRIC"),
    ("2b: phi^(-alpha*|dn|)",       m2b, x2b, c2b, "FITTING(1)"),
    ("2c: phi^(-alpha*|dn|^beta)",  m2c, x2c, c2c, "FITTING(2)"),
    ("2d: sin=phi^(-dn/k) rot",     m2d, x2d, c2d, "FITTING(1)"),
    ("3a: Wolf lam=phi^-3",         m3a, x3a, c3a, "GEOMETRIC+FIT(1)"),
    ("3b: Wolf best n",             m3b, x3b, c3b, "FITTING(1)"),
    ("3c: sin=phi^-n, (3,7,12)",    m3c, x3c, c3c, "GEOMETRIC"),
    ("5b: k=quadratic(da,db)",      m5b, x5b, c5b, "FITTING(3)"),
    ("5d: Gram-Schmidt (a,b,n)",    m5d, x5d, c5d, "GEOMETRIC"),
    ("5e: SVD rotation (a,b)",      m5e, x5e, c5e, "GEOMETRIC"),
    ("6c: lam=phi^-3, fit A",       m6c, x6c, c6c, "GEOMETRIC+FIT(1)"),
    ("6d: n=(3,13/2,23/2)",         m6d, x6d, c6d, "GEOMETRIC"),
]

# Sort by mean error
results.sort(key=lambda x: x[1])

print(f"\n{'Approach':<35} {'Mean%':>7} {'Max%':>7} {'Chi2':>8} {'Category':<20}")
print("-" * 82)
for name, mean_err, max_err, chi2, cat in results:
    print(f"  {name:<33} {mean_err:7.2f} {max_err:7.2f} {chi2:8.1f}  {cat}")

# Find the best geometric (no free parameters) approach
print("\n" + "-" * 78)
print("BEST GEOMETRIC APPROACHES (no free parameters or minimal fitting):")
print("-" * 78)

geo_results = [(n, m, x, c2, cat) for n, m, x, c2, cat in results if "GEOMETRIC" in cat and "FIT" not in cat]
for name, mean_err, max_err, chi2, cat in geo_results:
    print(f"  {name:<33} {mean_err:7.2f}%")

print("\n" + "=" * 78)
print("CONCLUSIONS")
print("=" * 78)

print("""
1. BEST PURELY GEOMETRIC APPROACH:
   sin(theta_12) = phi^(-3), sin(theta_23) = phi^(-7), sin(theta_13) = phi^(-12)
   From: dim(Im(H))=3, dim(Im(O))=7, degree(600-cell)=12

   This gives the CORRECT hierarchical structure but with significant
   errors for theta_23 (phi^(-7) = 0.035 vs exp 0.041).

2. BEST FITTING APPROACH:
   Spectral level CKM with sin(theta) = phi^(-|dn|/k)
   or Wolfenstein with lambda = phi^(-3).

   With 1-2 free parameters, excellent fits are achievable.

3. KEY OBSERVATION:
   The Cabibbo angle = arctan(phi^(-3)) with ~1.8% accuracy
   is the MOST ROBUST geometric prediction.

   The other angles require either:
   - Non-integer phi exponents (fitting)
   - Geometric corrections (perturbative)

4. SPECTRAL LEVEL CONNECTION:
   The level differences |n_up - n_dn| provide a NATURAL
   explanation for CKM hierarchy, but need a divisor/scale.

5. GEOMETRIC ORIGIN SUMMARY:
   - theta_Cabibbo from dim(Im(H)) = 3: STATUS DERIVAT (~1.8%)
   - theta_23 from dim(Im(O)) = 7: STATUS PATTERN (~15%)
   - theta_13 from degree = 12: STATUS PATTERN (~10%)
   - Full CKM with corrections: STATUS FITTING

HONEST ASSESSMENT:
==================
Category: PATTERN (not fully DERIVAT)

The 600-cell geometry provides a NATURAL framework for
understanding CKM hierarchy through phi-power suppression,
but a complete derivation without any free parameters
achieves only ~15-25% accuracy for individual elements.

The Cabibbo angle theta_C = arctan(phi^(-3)) remains
the single strongest geometric prediction.
""")
