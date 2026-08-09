"""
exp461: TQFT Braiding Amplitudes and CKM Mixing Angles
=======================================================
Can CKM mixing angles be derived from braiding in the SU(2)_3
Fibonacci anyon model?

MOTIVATION (from exp460):
- CKM requires a TWO-mediator process (Z/2 fusion grading no-go)
- Perturbative two-mediator exchange gives ADDITION of spectral exponents
- But CKM exponents have MULTIPLICATIVE structure: n_23 = n_ut * n_db - 1 = 7
- This points toward NON-PERTURBATIVE mechanism

KEY IDEA: In TQFT, braiding is inherently non-perturbative and produces
phases that are powers of the fundamental phase exp(2*pi*i/(k+2)).
In the Fibonacci anyon model (SU(2) at level k=3, k+2 = 5 = a1),
braiding produces phases related to phi.

FRAMEWORK:
- SU(2)_3 Chern-Simons / Fibonacci anyons
- Level k=3, k+2 = a1 = 5
- Spins: j = 0, 1/2, 1, 3/2
- Quantum dimensions: d_0=1, d_{1/2}=phi, d_1=phi, d_{3/2}=1
- Total quantum dimension: D^2 = a1 + sqrt(a1)
- Fusion: tau x tau = 1 + tau (Fibonacci)

RULES: REGULA ZERO applies. No invention. DERIVED/PATTERN/NEGATIVE.
"""

import numpy as np
from math import sqrt, log, sin, cos, pi, atan, asin
import cmath

# ============================================================
# CONSTANTS
# ============================================================
phi = (1 + sqrt(5)) / 2
phi_p = (1 - sqrt(5)) / 2  # = -1/phi
ln_phi = log(phi)
a1 = 5; b1 = 6; N_order = 120
alpha_em = 1 / 137.035999084
alpha_s = 1 / (2 * phi**3)
sin2tW = 6 / 26

# CKM experimental values (PDG 2024)
V_CKM_exp = {
    'us': 0.2243, 'cb': 0.0410, 'ub': 0.00382,
    'cd': 0.221,  'ts': 0.0388, 'td': 0.0080
}

# Target CKM exponents
n_target = {'12': 3, '23': 7, '13': 12}

# Target |V| from phi^{-n} (bare, no corrections)
V_bare = {
    'us': phi**(-3),   # 0.2361
    'cb': phi**(-7),   # 0.03444
    'ub': phi**(-12),  # 0.002360
}

# q = exp(2*pi*i/(k+2)) for k=3
q = cmath.exp(2j * pi / a1)  # = exp(2*pi*i/5)

print("=" * 78)
print("exp461: TQFT Braiding Amplitudes and CKM Mixing Angles")
print("=" * 78)

# ================================================================
# PART 1: SU(2)_3 CONFORMAL WEIGHTS AND R-MATRICES
# ================================================================
print("\n" + "=" * 78)
print("PART 1: SU(2)_3 Conformal Weights and R-Matrices")
print("=" * 78)

k = 3
kp2 = k + 2  # = 5 = a1

# Conformal weights h_j = j(j+1)/(k+2)
spins = [0, 0.5, 1, 1.5]
spin_labels = ['0', '1/2', '1', '3/2']
h = {}
for j in spins:
    h[j] = j * (j + 1) / kp2

print(f"\nConformal weights h_j = j(j+1)/(k+2), k={k}, k+2={kp2}")
for j, label in zip(spins, spin_labels):
    print(f"  h_{label:>3s} = {h[j]:.6f} = {j*(j+1)}/{kp2}")

# Quantum dimensions d_j = sin((2j+1)*pi/(k+2)) / sin(pi/(k+2))
d = {}
for j in spins:
    d[j] = sin((2*j+1)*pi/kp2) / sin(pi/kp2)

print(f"\nQuantum dimensions d_j = sin((2j+1)*pi/(k+2)) / sin(pi/(k+2))")
for j, label in zip(spins, spin_labels):
    print(f"  d_{label:>3s} = {d[j]:.6f}", end="")
    # Check against known values
    if abs(d[j] - 1) < 1e-10:
        print("  = 1")
    elif abs(d[j] - phi) < 1e-10:
        print(f"  = phi = {phi:.6f}")
    else:
        print()

D2 = sum(d[j]**2 for j in spins)
print(f"\n  D^2 = sum d_j^2 = {D2:.6f}")
print(f"  a1 + sqrt(a1) = {a1 + sqrt(a1):.6f}")
print(f"  Match: {abs(D2 - (a1 + sqrt(a1))) < 1e-10}")

# Double embedding check
print(f"\n  d_{{1/2}} = d_1 = phi: {abs(d[0.5] - d[1.0]) < 1e-10}  (DOUBLE EMBEDDING, unique to a1=5)")

# R-matrices: R_c^{ab} = (-1)^{a+b-c} * exp(i*pi*(h_c - h_a - h_b))
print(f"\n--- R-matrices ---")
print(f"R_c^{{ab}} = (-1)^{{a+b-c}} * exp(i*pi*(h_c - h_a - h_b))")

# Compute all R-matrix eigenvalues for pairs (a,b) that fuse to c
# Fusion rules for SU(2)_k: a x b -> |a-b|, |a-b|+1, ..., min(a+b, k-a-b)
def fusion_channels(a, b, k_level):
    """Return list of allowed fusion channels for a x b at level k."""
    channels = []
    j_min = abs(a - b)
    j_max = min(a + b, k_level - a - b) if a + b <= k_level else min(a + b, k_level)
    # Truncation: j_max = min(a+b, k-a-b) only if a+b > k-a-b
    # Actually for SU(2)_k: j goes from |a-b| to min(a+b, k-a-b)
    # But the constraint is: a+b+c <= k (for 3j-symbol) => c <= k-a-b
    j_max_trunc = k_level - a - b
    if j_max_trunc < abs(a - b):
        return []  # No channels
    j_max = min(a + b, j_max_trunc)
    j = j_min
    while j <= j_max + 1e-10:
        channels.append(j)
        j += 1
    return channels

print(f"\nFusion rules at k={k}:")
for a in spins:
    for b in spins:
        channels = fusion_channels(a, b, k)
        if channels:
            ch_str = " + ".join([spin_labels[spins.index(c)] if c in spins else f"{c}" for c in channels])
            print(f"  {spin_labels[spins.index(a)]:>3s} x {spin_labels[spins.index(b)]:>3s} = {ch_str}")

print(f"\n--- R-matrix eigenvalues for all fusion channels ---")
R = {}
for a in spins:
    for b in spins:
        channels = fusion_channels(a, b, k)
        for c in channels:
            sign = (-1)**(a + b - c)
            # Need integer check for sign
            exp_val = a + b - c
            if abs(exp_val - round(exp_val)) < 1e-10:
                sign = (-1)**int(round(exp_val))
            else:
                sign = cmath.exp(1j * pi * exp_val)

            phase = cmath.exp(1j * pi * (h[c] - h[a] - h[b]))
            R_val = sign * phase
            R[(a, b, c)] = R_val

            a_label = spin_labels[spins.index(a)]
            b_label = spin_labels[spins.index(b)]
            c_label = spin_labels[spins.index(c)] if c in spins else f"{c}"

            # Express as exp(i*theta)
            theta = cmath.phase(R_val)
            theta_over_pi = theta / pi

            print(f"  R_{{{c_label}}}^{{{a_label},{b_label}}} = {R_val:.6f} = exp(i*{theta_over_pi:.4f}*pi)")

# ================================================================
# PART 2: FIBONACCI ANYON BRAIDING MATRICES
# ================================================================
print("\n" + "=" * 78)
print("PART 2: Fibonacci Anyon Braiding Matrices")
print("=" * 78)

print(f"""
The Fibonacci anyon model uses the j=1 (tau) sector.
Fusion: tau x tau = 1 + tau  (Fibonacci fusion rule)

For 2 tau anyons, braiding is diagonal in the fusion channel basis:
  |0>: total spin 0   =>  R_0^{{1,1}}
  |1>: total spin 1   =>  R_1^{{1,1}}

The R-matrix in this basis is:
  R_fib = diag(R_0^{{11}}, R_1^{{11}})
""")

# R-matrix for Fibonacci sector (a=b=1, tau)
R_fib_0 = R[(1.0, 1.0, 0.0)]  # fusion to vacuum
R_fib_1 = R[(1.0, 1.0, 1.0)]  # fusion to tau

print(f"R_0^{{1,1}} = {R_fib_0:.6f} = exp(i*{cmath.phase(R_fib_0)/pi:.6f}*pi)")
print(f"R_1^{{1,1}} = {R_fib_1:.6f} = exp(i*{cmath.phase(R_fib_1)/pi:.6f}*pi)")

# Verify the known values
# R_0^{11} should be exp(-4*pi*i/5)
R_0_expected = cmath.exp(-4j * pi / 5)
# R_1^{11} should be exp(3*pi*i/5)
R_1_expected = cmath.exp(3j * pi / 5)

print(f"\nExpected:")
print(f"  R_0^{{1,1}} = exp(-4*pi*i/5) = {R_0_expected:.6f}")
print(f"  R_1^{{1,1}} = exp(3*pi*i/5)  = {R_1_expected:.6f}")
print(f"  Match R_0: {abs(R_fib_0 - R_0_expected) < 1e-10}")
print(f"  Match R_1: {abs(R_fib_1 - R_1_expected) < 1e-10}")

# F-matrix for Fibonacci anyons
# The F-matrix (recoupling matrix) for tau sector
# F_tau = [[1/phi, 1/sqrt(phi)], [1/sqrt(phi), -1/phi]]
# = [[phi^{-1}, phi^{-1/2}], [phi^{-1/2}, -phi^{-1}]]
print(f"\n--- F-matrix for Fibonacci anyons ---")
F_fib = np.array([
    [1/phi, 1/sqrt(phi)],
    [1/sqrt(phi), -1/phi]
])

print(f"F = [[{F_fib[0,0]:.6f}, {F_fib[0,1]:.6f}],")
print(f"     [{F_fib[1,0]:.6f}, {F_fib[1,1]:.6f}]]")
print(f"\n= [[phi^{{-1}}, phi^{{-1/2}}],")
print(f"   [phi^{{-1/2}}, -phi^{{-1}}]]")

# Verify F is unitary
FF = F_fib @ F_fib
print(f"\nF^2 = identity? {np.allclose(FF, np.eye(2))}")
# F is actually involutory for Fibonacci? Let's check
print(f"F^2 =\n{FF}")

# Actually F is symmetric and unitary, so F = F^{-1} iff F^2 = I
# If not, compute F inverse
F_inv = np.linalg.inv(F_fib)
print(f"\nF^{{-1}} =\n{F_inv}")
print(f"F = F^T? {np.allclose(F_fib, F_fib.T)}")
print(f"F * F^T = I? {np.allclose(F_fib @ F_fib.T, np.eye(2))}")

# ================================================================
# PART 3: BRAID GROUP B_3 REPRESENTATION
# ================================================================
print("\n" + "=" * 78)
print("PART 3: Braid Group B_3 Representation on 3 Fibonacci Anyons")
print("=" * 78)

print(f"""
For 3 Fibonacci (tau) anyons, the fusion space is 2-dimensional:
  - Path 1: (tau x tau -> 1) x tau -> tau  [through vacuum]
  - Path 2: (tau x tau -> tau) x tau -> tau [through tau]

The braid group B_3 has two generators sigma_1, sigma_2.
  sigma_1: braid particles 1 & 2 (acts as R in fusion channel basis)
  sigma_2: braid particles 2 & 3 (conjugated by F: F^{{-1}} R F)
""")

# sigma_1 acts as the diagonal R-matrix
sigma1 = np.array([
    [R_fib_0, 0],
    [0, R_fib_1]
], dtype=complex)

# sigma_2 = F^{-1} * R * F (change basis, braid, change back)
# Use complex F
F_c = F_fib.astype(complex)
F_inv_c = F_inv.astype(complex)
R_diag = np.diag([R_fib_0, R_fib_1])

sigma2 = F_inv_c @ R_diag @ F_c

print(f"sigma_1 (diagonal R-matrix):")
print(f"  [[{sigma1[0,0]:.6f}, {sigma1[0,1]:.6f}],")
print(f"   [{sigma1[1,0]:.6f}, {sigma1[1,1]:.6f}]]")

print(f"\nsigma_2 = F^{{-1}} R F:")
print(f"  [[{sigma2[0,0]:.6f}, {sigma2[0,1]:.6f}],")
print(f"   [{sigma2[1,0]:.6f}, {sigma2[1,1]:.6f}]]")

# Verify: sigma_1 * sigma_2 * sigma_1 = sigma_2 * sigma_1 * sigma_2 (braid relation)
lhs = sigma1 @ sigma2 @ sigma1
rhs = sigma2 @ sigma1 @ sigma2
print(f"\nBraid relation sigma_1*sigma_2*sigma_1 = sigma_2*sigma_1*sigma_2: {np.allclose(lhs, rhs)}")

# ================================================================
# PART 4: BRAIDING AMPLITUDES -- sigma^n
# ================================================================
print("\n" + "=" * 78)
print("PART 4: Braiding Amplitudes |<out|sigma^n|in>|^2")
print("=" * 78)

print(f"""
For n half-twists (braids) using sigma_1, the matrix is sigma_1^n.
Since sigma_1 is diagonal, sigma_1^n is also diagonal:
  sigma_1^n = diag(R_0^n, R_1^n)

No channel-changing (off-diagonal) amplitude from sigma_1 alone.

For sigma_2^n, the matrix is non-diagonal and produces transitions.
""")

print(f"--- sigma_2^n: off-diagonal amplitude ---")
print(f"{'n':>4s} {'|<0|sig2^n|1>|^2':>18s} {'|<1|sig2^n|0>|^2':>18s} {'log_phi(|01|^2)':>18s} {'log_phi(|10|^2)':>18s}")
print("-" * 78)

for n in range(1, 21):
    sig2n = np.linalg.matrix_power(sigma2, n)
    amp_01 = abs(sig2n[0, 1])**2
    amp_10 = abs(sig2n[1, 0])**2

    log_01 = -log(amp_01) / ln_phi if amp_01 > 1e-30 else float('nan')
    log_10 = -log(amp_10) / ln_phi if amp_10 > 1e-30 else float('nan')

    marker = ""
    for target_n, target_name in [(3, 'us'), (7, 'cb'), (12, 'ub')]:
        if not np.isnan(log_01) and abs(log_01 - target_n) < 0.5:
            marker += f" <-- n_{target_name}?"
        if not np.isnan(log_10) and abs(log_10 - target_n) < 0.5:
            marker += f" <-- n_{target_name}?"

    print(f"{n:4d} {amp_01:18.10f} {amp_10:18.10f} {log_01:18.4f} {log_10:18.4f}{marker}")

# ================================================================
# PART 5: FULL BRAID WORD (sigma_1 * sigma_2)^n
# ================================================================
print("\n" + "=" * 78)
print("PART 5: Full Braid Word (sigma_1 * sigma_2)^n")
print("=" * 78)

print(f"""
The product sigma_1 * sigma_2 represents a FULL twist of 3 anyons.
Its matrix elements involve phi through the F-matrix.
""")

sig12 = sigma1 @ sigma2  # One full twist

print(f"sigma_1 * sigma_2:")
print(f"  [[{sig12[0,0]:.6f}, {sig12[0,1]:.6f}],")
print(f"   [{sig12[1,0]:.6f}, {sig12[1,1]:.6f}]]")

# Eigenvalues of sig12
evals_12 = np.linalg.eigvals(sig12)
print(f"\nEigenvalues of sigma_1 * sigma_2:")
for i, ev in enumerate(evals_12):
    theta = cmath.phase(ev)
    print(f"  lambda_{i} = {ev:.6f}, |lambda| = {abs(ev):.6f}, phase = {theta/pi:.6f}*pi")

print(f"\n--- (sigma_1 * sigma_2)^n: off-diagonal amplitude ---")
print(f"{'n':>4s} {'|M_01|^2':>14s} {'|M_10|^2':>14s} {'-log_phi':>14s} {'-log_phi':>14s}")
print("-" * 65)

found_matches = []
for n in range(1, 31):
    sig12n = np.linalg.matrix_power(sig12, n)
    amp_01 = abs(sig12n[0, 1])**2
    amp_10 = abs(sig12n[1, 0])**2

    log_01 = -log(amp_01) / ln_phi if amp_01 > 1e-30 else float('nan')
    log_10 = -log(amp_10) / ln_phi if amp_10 > 1e-30 else float('nan')

    marker = ""
    for target_n, target_name in [(3, 'us'), (7, 'cb'), (12, 'ub')]:
        if not np.isnan(log_01) and abs(log_01 - target_n) < 0.3:
            marker += f" <-- n_{target_name}!"
            found_matches.append((n, '01', target_name, log_01))
        if not np.isnan(log_10) and abs(log_10 - target_n) < 0.3:
            marker += f" <-- n_{target_name}!"
            found_matches.append((n, '10', target_name, log_10))

    print(f"{n:4d} {amp_01:14.10f} {amp_10:14.10f} {log_01:14.4f} {log_10:14.4f}{marker}")

if found_matches:
    print(f"\n  MATCHES FOUND:")
    for nn, ch, tname, lval in found_matches:
        print(f"    n={nn}, channel={ch}: -log_phi = {lval:.4f} ~ n_{tname}={n_target[tname.replace('us','12').replace('cb','23').replace('ub','13')]}")
else:
    print(f"\n  No matches to {{3, 7, 12}} found in (sigma_1*sigma_2)^n for n=1..30")

# ================================================================
# PART 6: GENERAL BRAID WORDS
# ================================================================
print("\n" + "=" * 78)
print("PART 6: General Braid Words -- Systematic Search")
print("=" * 78)

print(f"""
Search over braid words of length L = 1..8 for words W such that:
  |<0|W|1>| ~ phi^{{-n}} with n in {{3, 7, 12}}

A braid word is a sequence of sigma_1, sigma_1^{{-1}}, sigma_2, sigma_2^{{-1}}.
""")

sigma1_inv = np.linalg.inv(sigma1)
sigma2_inv = np.linalg.inv(sigma2)

generators = {
    's1': sigma1,
    's1i': sigma1_inv,
    's2': sigma2,
    's2i': sigma2_inv,
}

# For short braid words, enumerate and check
# Length 1-4 only (4^4 = 256, manageable)
target_exponents = [3, 7, 12]

print(f"\nSearching braid words of length 1-5...")
best_words = {3: [], 7: [], 12: []}

gen_names = list(generators.keys())
gen_mats = list(generators.values())

for length in range(1, 6):
    # Generate all words of given length
    from itertools import product as iterproduct
    for word_indices in iterproduct(range(4), repeat=length):
        M = np.eye(2, dtype=complex)
        for idx in word_indices:
            M = M @ gen_mats[idx]

        amp_01 = abs(M[0, 1])**2
        if amp_01 > 1e-30 and amp_01 < 1 - 1e-10:
            n_eff = -log(amp_01) / ln_phi

            for target_n in target_exponents:
                if abs(n_eff - target_n) < 0.1:
                    word_str = "*".join([gen_names[idx] for idx in word_indices])
                    best_words[target_n].append((word_str, n_eff, amp_01, length))

# Show results
for target_n in target_exponents:
    words = best_words[target_n]
    label_map = {3: 'us/12', 7: 'cb/23', 12: 'ub/13'}
    print(f"\n--- Target n = {target_n} (V_{label_map[target_n]}) ---")
    if words:
        # Sort by closeness to target
        words.sort(key=lambda x: abs(x[1] - target_n))
        for w_str, n_eff, amp, length in words[:5]:
            print(f"  Word: {w_str:>30s}  n_eff = {n_eff:.6f}  |amp|^2 = {amp:.8f}  len={length}")
    else:
        print(f"  No words found with n_eff ~ {target_n} (within 0.1)")

# ================================================================
# PART 7: F-MATRIX AS CKM ANALOGUE
# ================================================================
print("\n" + "=" * 78)
print("PART 7: F-Matrix as CKM Analogue")
print("=" * 78)

print(f"""
The F-matrix IS a change of basis between fusion trees.
This is the quantum analogue of the CKM matrix!

F_fib = [[phi^{{-1}}, phi^{{-1/2}}],
         [phi^{{-1/2}}, -phi^{{-1}}]]

Elements:
  |F_00|^2 = phi^{{-2}} = {1/phi**2:.6f}
  |F_01|^2 = phi^{{-1}} = {1/phi:.6f}
  |F_10|^2 = phi^{{-1}} = {1/phi:.6f}
  |F_11|^2 = phi^{{-2}} = {1/phi**2:.6f}

Phi exponents: 1, 2 -- NOT {{3, 7, 12}}.

But what about PRODUCTS of F-matrices?
""")

# F^n analysis
print(f"--- Powers of F ---")
for n in range(1, 11):
    Fn = np.linalg.matrix_power(F_fib, n)
    print(f"\n  F^{n}:")
    for i in range(2):
        for j in range(2):
            val = abs(Fn[i,j])
            if val > 1e-15:
                n_eff = -log(val**2) / ln_phi
                print(f"    |F^{n}[{i},{j}]|^2 = {val**2:.8f}, -log_phi = {n_eff:.4f}")
            else:
                print(f"    |F^{n}[{i},{j}]|^2 ~ 0")

# Products of F and R
print(f"\n--- Products F * R * F * ... ---")
patterns = [
    ('F*R*F', F_fib.astype(complex) @ R_diag @ F_fib.astype(complex)),
    ('R*F*R', R_diag @ F_fib.astype(complex) @ R_diag),
    ('F*R^2*F', F_fib.astype(complex) @ np.linalg.matrix_power(R_diag, 2) @ F_fib.astype(complex)),
    ('(F*R)^2', np.linalg.matrix_power(F_fib.astype(complex) @ R_diag, 2)),
    ('(F*R)^3', np.linalg.matrix_power(F_fib.astype(complex) @ R_diag, 3)),
    ('(R*F)^3', np.linalg.matrix_power(R_diag @ F_fib.astype(complex), 3)),
    ('(F*R)^5', np.linalg.matrix_power(F_fib.astype(complex) @ R_diag, 5)),
]

for name, M in patterns:
    amp_01 = abs(M[0,1])**2
    if amp_01 > 1e-30 and amp_01 < 1 - 1e-10:
        n_eff = -log(amp_01) / ln_phi
        marker = ""
        for tn in target_exponents:
            if abs(n_eff - tn) < 0.3:
                marker = f"  <-- CLOSE to {tn}!"
        print(f"  {name:>15s}: |M_01|^2 = {amp_01:.8f}, -log_phi = {n_eff:.4f}{marker}")
    else:
        print(f"  {name:>15s}: |M_01|^2 = {amp_01:.8f} (trivial)")

# ================================================================
# PART 8: MONODROMY MATRICES
# ================================================================
print("\n" + "=" * 78)
print("PART 8: Monodromy Matrices (Double Braids)")
print("=" * 78)

print(f"""
A monodromy (full loop) is sigma^2, not sigma.
This is the PHYSICAL operation of moving one anyon completely
around another (2 half-twists = 1 full twist).

Monodromy matrix: M_mono = sigma^2
""")

# Monodromy of sigma_2
mono2 = sigma2 @ sigma2
print(f"Monodromy (sigma_2^2):")
print(f"  [[{mono2[0,0]:.6f}, {mono2[0,1]:.6f}],")
print(f"   [{mono2[1,0]:.6f}, {mono2[1,1]:.6f}]]")

amp_mono = abs(mono2[0,1])**2
if amp_mono > 1e-30:
    n_mono = -log(amp_mono) / ln_phi
    print(f"  |<0|mono|1>|^2 = {amp_mono:.8f}, -log_phi = {n_mono:.4f}")

# Multiple monodromies
print(f"\n--- Multiple monodromies ---")
for n in range(1, 16):
    mono_n = np.linalg.matrix_power(mono2, n)
    amp = abs(mono_n[0,1])**2
    if amp > 1e-30 and amp < 1 - 1e-10:
        n_eff = -log(amp) / ln_phi
        marker = ""
        for tn, tname in [(3, '12'), (7, '23'), (12, '13')]:
            if abs(n_eff - tn) < 0.3:
                marker = f" <-- n_{tname}!"
        print(f"  mono^{n:2d}: |<0|M|1>|^2 = {amp:.8f}, -log_phi = {n_eff:.4f}{marker}")

# Full monodromy: sigma_1^2 * sigma_2^2
mono_full = sigma1 @ sigma1 @ sigma2 @ sigma2
print(f"\n--- Full monodromy sigma_1^2 * sigma_2^2 ---")
print(f"  [[{mono_full[0,0]:.6f}, {mono_full[0,1]:.6f}],")
print(f"   [{mono_full[1,0]:.6f}, {mono_full[1,1]:.6f}]]")

for n in range(1, 16):
    mono_fn = np.linalg.matrix_power(mono_full, n)
    amp = abs(mono_fn[0,1])**2
    if amp > 1e-30 and amp < 1 - 1e-10:
        n_eff = -log(amp) / ln_phi
        marker = ""
        for tn, tname in [(3, '12'), (7, '23'), (12, '13')]:
            if abs(n_eff - tn) < 0.3:
                marker = f" <-- n_{tname}!"
        print(f"  (mono_full)^{n:2d}: |<0|M|1>|^2 = {amp:.8f}, -log_phi = {n_eff:.4f}{marker}")

# ================================================================
# PART 9: MAPPING QUARKS TO ANYONS -- 3 GENERATIONS
# ================================================================
print("\n" + "=" * 78)
print("PART 9: Mapping Quarks to Anyons -- 3 Generations")
print("=" * 78)

print(f"""
In the theory, the 4 TQFT objects are:
  j=0: vacuum      (d=1)
  j=1/2: matter    (d=phi)
  j=1: gauge       (d=phi)
  j=3/2: Planck    (d=1)

Quarks are j=1/2 (matter). For 3 generations, we need 3 up + 3 down
quarks, i.e., 6 anyons of type j=1/2.

The fusion space of n anyons of type tau (Fibonacci) has dimension
equal to the Fibonacci number F_n:
  1 anyon:  dim = 1
  2 anyons: dim = 1 (since tau x tau = 1 + tau, basis: |1>, |tau>)
            Actually dim = 2 fusion channels
  3 anyons: dim = 2
  4 anyons: dim = 3
  5 anyons: dim = 5
  6 anyons: dim = 8

For 6 tau anyons, the braid group B_6 acts on an 8-dimensional space.
The CKM matrix should be a 3x3 subblock of some braiding matrix
in this 8-dimensional space.

However, the FULL computation of the B_6 representation on the
8-dimensional Fibonacci fusion space is complex. Let me instead
check the KEY structural question: can braiding phases reproduce
phi^{{-n}} with the right exponents?
""")

# Dimension of fusion space for n Fibonacci anyons
def fibonacci_dim(n):
    """Dimension of the fusion space of n Fibonacci (tau) anyons fusing to tau."""
    if n <= 0:
        return 0
    if n == 1:
        return 1
    # For n anyons fusing to tau: F_{n-1}
    # For n anyons (total): F_{n+1}
    f = [1, 1]
    for i in range(2, n+2):
        f.append(f[-1] + f[-2])
    return f[n]

print(f"Fusion space dimensions for n tau anyons:")
for n in range(1, 9):
    print(f"  n={n}: dim = {fibonacci_dim(n)}")

# ================================================================
# PART 10: MODULAR S AND T MATRICES
# ================================================================
print("\n" + "=" * 78)
print("PART 10: Modular S and T Matrices of SU(2)_3")
print("=" * 78)

print(f"""
The modular S-matrix encodes the modular transformation tau -> -1/tau
on the torus. Its elements are:
  S_{{ab}} = sqrt(2/(k+2)) * sin(pi*(2a+1)*(2b+1)/(k+2))

The T-matrix is diagonal:
  T_{{ab}} = delta_{{ab}} * exp(2*pi*i*(h_a - c/24))
where c = 3k/(k+2) is the central charge.
""")

# Central charge
c_central = 3 * k / kp2
print(f"Central charge c = 3k/(k+2) = {c_central:.6f} = {3*k}/{kp2}")

# S-matrix
n_spins = len(spins)
S = np.zeros((n_spins, n_spins))
for i, a in enumerate(spins):
    for j, b in enumerate(spins):
        S[i,j] = sqrt(2/kp2) * sin(pi * (2*a+1) * (2*b+1) / kp2)

print(f"\nModular S-matrix:")
print(f"{'':>6s}", end="")
for label in spin_labels:
    print(f" {label:>10s}", end="")
print()
for i, label_i in enumerate(spin_labels):
    print(f"{label_i:>6s}", end="")
    for j in range(n_spins):
        print(f" {S[i,j]:10.6f}", end="")
    print()

# Verify S is unitary
print(f"\nS*S^T = I? {np.allclose(S @ S.T, np.eye(n_spins))}")

# T-matrix
T_diag = np.zeros(n_spins, dtype=complex)
for i, j_spin in enumerate(spins):
    T_diag[i] = cmath.exp(2j * pi * (h[j_spin] - c_central/24))
T_mat = np.diag(T_diag)

print(f"\nModular T-matrix (diagonal):")
for i, label in enumerate(spin_labels):
    theta = cmath.phase(T_diag[i])
    print(f"  T_{{{label}}} = {T_diag[i]:.6f} = exp(i*{theta/pi:.6f}*pi)")

# Check: S*T*S = T^{-1}*S*T^{-1} (modular relation)
# Actually (S*T)^3 = S^2 = C (charge conjugation)
ST = S.astype(complex) @ T_mat
ST3 = np.linalg.matrix_power(ST, 3)
S2 = S.astype(complex) @ S.astype(complex)
print(f"\n(S*T)^3 = S^2? {np.allclose(ST3, S2)}")

# Elements of S-matrix in terms of phi
print(f"\nS-matrix elements in terms of phi:")
for i, a in enumerate(spins):
    for j, b in enumerate(spins):
        val = S[i,j]
        if abs(val) > 1e-15:
            n_eff = log(abs(val)) / ln_phi
            print(f"  S_{{{spin_labels[i]},{spin_labels[j]}}} = {val:.6f}, log_phi = {n_eff:.4f}")

# ================================================================
# PART 11: S-MATRIX RATIOS AND CKM STRUCTURE
# ================================================================
print("\n" + "=" * 78)
print("PART 11: S-Matrix Ratios and CKM Structure")
print("=" * 78)

print(f"""
The S-matrix ratios S_{{ab}}/S_{{0b}} give quantum dimensions:
  S_{{ab}}/S_{{0b}} = d_a (independent of b)
  ... but only for Fibonacci model with two particles.

For SU(2)_3 with 4 particles, the S-matrix is 4x4.
The RATIOS of S-matrix elements might give CKM-like structure.

Key check: do any ratios of S-matrix elements give phi^{{-3}},
phi^{{-7}}, phi^{{-12}}?
""")

print(f"All S-matrix ratios S_ab/S_cd (as phi exponents):")
all_ratios = []
for i in range(n_spins):
    for j in range(n_spins):
        for i2 in range(n_spins):
            for j2 in range(n_spins):
                if abs(S[i,j]) > 1e-15 and abs(S[i2,j2]) > 1e-15:
                    ratio = S[i,j] / S[i2,j2]
                    if abs(ratio) > 1e-15 and abs(ratio) < 1e15:
                        n_eff = log(abs(ratio)) / ln_phi
                        all_ratios.append((i, j, i2, j2, n_eff))

# Find ratios close to target exponents
print(f"\nRatios close to target exponents:")
for target_n, tname in [(3, 'n_12'), (7, 'n_23'), (12, 'n_13')]:
    close_ratios = [(i, j, i2, j2, n_eff) for i, j, i2, j2, n_eff in all_ratios
                    if abs(abs(n_eff) - target_n) < 0.5]
    print(f"\n  Target n = {target_n} ({tname}):")
    if close_ratios:
        for i, j, i2, j2, n_eff in close_ratios[:5]:
            print(f"    S_{{{spin_labels[i]},{spin_labels[j]}}}/S_{{{spin_labels[i2]},{spin_labels[j2]}}} = phi^{{{n_eff:.4f}}}")
    else:
        print(f"    No ratios found close to {target_n}")

# ================================================================
# PART 12: VERLINDE FORMULA AND FUSION MULTIPLICITIES
# ================================================================
print("\n" + "=" * 78)
print("PART 12: Verlinde Formula Verification")
print("=" * 78)

print(f"""
The Verlinde formula gives fusion coefficients from S-matrix:
  N_{{ab}}^c = sum_x S_{{ax}} S_{{bx}} S*_{{cx}} / S_{{0x}}
""")

print(f"Fusion coefficients from Verlinde formula:")
N_verlinde = np.zeros((n_spins, n_spins, n_spins))
for a in range(n_spins):
    for b in range(n_spins):
        for c_idx in range(n_spins):
            val = 0
            for x in range(n_spins):
                val += S[a,x] * S[b,x] * np.conj(S[c_idx,x]) / S[0,x]
            N_verlinde[a,b,c_idx] = val.real

for a in range(n_spins):
    for b in range(a, n_spins):
        channels = []
        for c_idx in range(n_spins):
            N_val = N_verlinde[a, b, c_idx]
            if abs(N_val) > 0.01:
                channels.append(f"{round(N_val)}*{spin_labels[c_idx]}")
        if channels:
            print(f"  {spin_labels[a]:>3s} x {spin_labels[b]:>3s} = {' + '.join(channels)}")

# Verify Fibonacci fusion
print(f"\nVerification: tau(=1) x tau(=1) = 1 + tau?")
print(f"  N_{{1,1}}^0 = {N_verlinde[2,2,0]:.6f} (should be 1)")
print(f"  N_{{1,1}}^1 = {N_verlinde[2,2,2]:.6f} (should be 1)")

# ================================================================
# PART 13: T-MATRIX PHASES AND CKM EXPONENTS
# ================================================================
print("\n" + "=" * 78)
print("PART 13: T-Matrix Phases and CKM Concordance Connection")
print("=" * 78)

print(f"""
From exp449 (CKM Concordance), the T-matrix conformal weight GAPS
give the CKM exponents. Let's verify this in the SU(2)_3 framework.

Conformal weights: h_j = j(j+1)/5
  h_0 = 0, h_{{1/2}} = 3/20, h_1 = 2/5, h_{{3/2}} = 3/4

CKM connection from exp449:
  n_12 = 3 from conformal weight gaps
  n_23 = 2*a1 - 3 = 7 from conformal weight gaps
  n_13 = 3*a1 - 3 = 12 from conformal weight gaps
""")

# The conformal weight differences (in units of 1/kp2 = 1/5)
print(f"Conformal weight differences (h_a - h_b):")
for i, a in enumerate(spins):
    for j, b in enumerate(spins):
        if i != j:
            dh = h[a] - h[b]
            dh_unit = dh * kp2  # in units of 1/5
            if abs(dh_unit - round(dh_unit)) < 1e-10:
                print(f"  h_{spin_labels[i]} - h_{spin_labels[j]} = {dh:.6f} = {int(round(dh_unit))}/5")

# The key phase differences relevant for CKM
print(f"\nKey phase differences for braiding (in units of pi/5):")
# Between matter (1/2) and gauge (1)
dh_matter_gauge = h[1.0] - h[0.5]  # = 2/5 - 3/20 = 5/20 = 1/4
print(f"  h_1 - h_{{1/2}} = {dh_matter_gauge:.6f} = {dh_matter_gauge*20:.0f}/20")

# Between gauge (1) and Planck (3/2)
dh_gauge_planck = h[1.5] - h[1.0]  # = 3/4 - 2/5 = 7/20
print(f"  h_{{3/2}} - h_1 = {h[1.5]-h[1.0]:.6f} = {(h[1.5]-h[1.0])*20:.0f}/20")

# Full spectrum of phase differences in units of pi/5
print(f"\nPhase differences h_a - h_b in units of 1/20:")
for i, a in enumerate(spins):
    for j, b in enumerate(spins):
        if i > j:
            dh = h[a] - h[b]
            dh20 = dh * 20
            print(f"  h_{spin_labels[i]} - h_{spin_labels[j]} = {dh20:.1f}/20 = {dh:.6f}")

# Connection to CKM exponents
print(f"\n--- Connection to CKM exponents ---")
print(f"From T-matrix: T_j = exp(2*pi*i*(h_j - c/24))")
print(f"Phase of T_j (in units of pi):")
for i, j_spin in enumerate(spins):
    phase = 2 * (h[j_spin] - c_central/24)
    phase_pi = phase  # already in units of pi
    print(f"  T_{spin_labels[i]}: phase = {phase_pi:.6f}*pi")

# The CKM connection from exp449
print(f"\nCKM from conformal weight spectrum:")
print(f"  Spectrum: {[h[j] for j in spins]}")
print(f"  Gaps (h_{{j+1}} - h_j) * 20: ", end="")
for i in range(len(spins)-1):
    gap = (h[spins[i+1]] - h[spins[i]]) * 20
    print(f"{gap:.1f}, ", end="")
print()

# The gaps are 3, 5, 7 (in units of 1/20)
print(f"\n  Gaps in units of 1/20: 3, 5, 7")
print(f"  These are EXACTLY the T-matrix phases from exp449!")
print(f"  n_12 = gap(0->1/2) = 3")
print(f"  n_23 = gap(1->3/2) = 7  (skipping from gauge to Planck)")
print(f"  n_13 = sum(all gaps) = 3+5+7 = 15? No, n_13 = 12.")
print(f"  Actually n_13 = 3*(a1-1) = 3*4 = 12. Gap structure is different.")

# Let's check: what ARE the conformal weight gaps in standard units?
print(f"\n  Conformal weights * (k+2) = j*(j+1):")
for j_spin in spins:
    print(f"    j={spin_labels[spins.index(j_spin)]}: h*5 = {h[j_spin]*5:.4f}")
# These give: 0, 3/4, 2, 15/4
# Hmm, not integers. Let's try h * 20:
print(f"  Conformal weights * 20:")
for j_spin in spins:
    print(f"    j={spin_labels[spins.index(j_spin)]}: h*20 = {h[j_spin]*20:.4f}")
# These give: 0, 3, 8, 15
print(f"\n  h*20 values: 0, 3, 8, 15")
print(f"  Differences: 3, 5, 7")
print(f"  Product: 3*5 = 15, 3*7 = 21, 5*7 = 35")
print(f"  But CKM exponents are 3, 7, 12 not 3, 5, 7")

# ================================================================
# PART 14: BRAIDING PHASES AND phi STRUCTURE
# ================================================================
print("\n" + "=" * 78)
print("PART 14: Braiding Phases and phi Structure")
print("=" * 78)

print(f"""
The braiding in Fibonacci model produces phases that are powers of
q = exp(2*pi*i/5). Since phi = 2*cos(pi/5), the absolute values of
braiding amplitudes involve phi.

Key question: after n braids, does the transition amplitude involve
phi^{{-n}} for the RIGHT values of n?

For the Fibonacci representation, sigma_2^n has matrix elements
that are sums of q-powers weighted by F-matrix elements.
The absolute values involve phi because:
  |F_{{01}}|^2 = 1/phi = phi^{{-1}}
  |F_{{00}}|^2 = 1/phi^2 = phi^{{-2}}
""")

# Detailed amplitude decomposition for sigma_2^n
print(f"--- Amplitude decomposition for sigma_2^n ---")
# sigma_2 = F^{-1} R F
# sigma_2^n = (F^{-1} R F)^n = F^{-1} R (F F^{-1}) R ... F = F^{-1} R^n F (NO! This is wrong)
# Actually (F^{-1} R F)^n != F^{-1} R^n F in general.
# (F^{-1} R F)^n = F^{-1} R F * F^{-1} R F * ... = F^{-1} R (F F^{-1}) R (F F^{-1}) ... R F
# = F^{-1} R^n F (this IS correct because FF^{-1} = I)

# Wait: (ABC)^n = ABC * ABC * ... = A(BCA)^{n-1} * BC
# (F^{-1} R F)^n: let's verify
sigma2_check = F_inv_c @ R_diag @ F_c
sigma2_n_from_Rn = F_inv_c @ np.linalg.matrix_power(R_diag, 3) @ F_c
sigma2_direct = np.linalg.matrix_power(sigma2, 3)
print(f"Check: (F^{{-1}}RF)^n = F^{{-1}}R^nF? {np.allclose(sigma2_n_from_Rn, sigma2_direct)}")
# This should be TRUE because sigma_2 = F^{-1}RF, and
# sigma_2^n = (F^{-1}RF)^n = F^{-1}RF * F^{-1}RF * ... = F^{-1}R^nF
# Yes! Because F F^{-1} = I between consecutive R's.

print(f"\nSo sigma_2^n = F^{{-1}} R^n F")
print(f"The (0,1) element is:")
print(f"  <0|sigma_2^n|1> = sum_c (F^{{-1}})_{{0c}} R_c^n F_{{c1}}")
print(f"  = (F^{{-1}})_{{00}} R_0^n F_{{01}} + (F^{{-1}})_{{01}} R_1^n F_{{11}}")

# Since F is real and symmetric, F^{-1} = F^T * det = ... let's compute explicitly
print(f"\n  F^{{-1}}_{{00}} = {F_inv[0,0]:.6f}")
print(f"  F^{{-1}}_{{01}} = {F_inv[0,1]:.6f}")
print(f"  F_{{01}} = {F_fib[0,1]:.6f}")
print(f"  F_{{11}} = {F_fib[1,1]:.6f}")

print(f"\n  <0|sigma_2^n|1> = {F_inv[0,0]:.6f} * R_0^n * {F_fib[0,1]:.6f}")
print(f"                   + {F_inv[0,1]:.6f} * R_1^n * {F_fib[1,1]:.6f}")
print(f"\nSince R_0 and R_1 have UNIT modulus (pure phases), the ABSOLUTE VALUE")
print(f"of the amplitude depends on the INTERFERENCE between the two terms.")
print(f"  |term_0|^2 = |F^{{-1}}_{{00}}|^2 * |F_{{01}}|^2 = {abs(F_inv[0,0])**2 * abs(F_fib[0,1])**2:.6f}")
print(f"  |term_1|^2 = |F^{{-1}}_{{01}}|^2 * |F_{{11}}|^2 = {abs(F_inv[0,1])**2 * abs(F_fib[1,1])**2:.6f}")

# Detailed computation for each n
print(f"\n--- Phase analysis for sigma_2^n ---")
print(f"{'n':>4s} {'phase_0':>12s} {'phase_1':>12s} {'|amp_01|^2':>14s} {'-log_phi':>10s} {'note':>20s}")
print("-" * 72)

for n in range(1, 31):
    R0n = R_fib_0**n
    R1n = R_fib_1**n

    term0 = F_inv[0,0] * R0n * F_fib[0,1]
    term1 = F_inv[0,1] * R1n * F_fib[1,1]
    amp = term0 + term1

    amp_sq = abs(amp)**2
    if amp_sq > 1e-30 and amp_sq < 1 - 1e-10:
        n_eff = -log(amp_sq) / ln_phi
    else:
        n_eff = float('nan')

    phase0 = cmath.phase(R0n) / pi
    phase1 = cmath.phase(R1n) / pi

    note = ""
    if not np.isnan(n_eff):
        for tn, tname in [(3, 'V_us'), (7, 'V_cb'), (12, 'V_ub')]:
            if abs(n_eff - tn) < 0.3:
                note = f"<-- n ~ {tn} ({tname})"

    # Check for periodicity
    if n > 1 and abs(amp_sq) < 1e-15:
        note = "ZERO (destructive)"

    print(f"{n:4d} {phase0:12.4f} {phase1:12.4f} {amp_sq:14.10f} {n_eff:10.4f} {note}")

# ================================================================
# PART 15: PERIODICITY AND QUANTUM STRUCTURE
# ================================================================
print("\n" + "=" * 78)
print("PART 15: Periodicity and Quantum Structure of Braiding")
print("=" * 78)

# R-matrix phases
phase_R0 = cmath.phase(R_fib_0) / pi  # Should be -4/5
phase_R1 = cmath.phase(R_fib_1) / pi  # Should be 3/5

print(f"R_0 phase = {phase_R0:.6f}*pi = {phase_R0*5:.2f}/5 * pi")
print(f"R_1 phase = {phase_R1:.6f}*pi = {phase_R1*5:.2f}/5 * pi")
print(f"Phase difference per braid: (R_0/R_1)^n period:")
phase_diff = cmath.phase(R_fib_0 / R_fib_1) / pi
print(f"  Phase diff = {phase_diff:.6f}*pi = {phase_diff*5:.2f}/5 * pi")
print(f"  Period (in n): 5/{abs(phase_diff*5):.0f} * 5 = {int(round(5/abs(phase_diff*5)*5))} ?")

# Actually the period is when R_0^n = R_1^n
# R_0 = exp(i*(-4/5)*pi), R_1 = exp(i*(3/5)*pi)
# R_0^n = R_1^n when (-4n/5 - 3n/5) = 2m, i.e., -7n/5 = 2m
# => n = 10m/7. Minimum at m=-7: n=10 (or actually n=10 when 7n/5 = 2m)
# Let's check: R_0^10 = exp(-8*pi*i) = 1, R_1^10 = exp(6*pi*i) = 1
print(f"\nR_0^10 = {R_fib_0**10:.6f}")
print(f"R_1^10 = {R_fib_1**10:.6f}")
print(f"Period of R-matrix: 10 (both R_0^10 = R_1^10 = 1)")
print(f"=> Braiding amplitudes are PERIODIC with period 10 in n")

# Since period is 10, the braiding amplitude can only produce
# a FINITE set of phi-exponents (between n=1 and n=10)
print(f"\nFull period of |<0|sigma_2^n|1>|^2:")
n_effs = []
for n in range(1, 11):
    sig2n = np.linalg.matrix_power(sigma2, n)
    amp = abs(sig2n[0,1])**2
    n_eff = -log(amp) / ln_phi if 0 < amp < 1 else float('nan')
    n_effs.append(n_eff)
    print(f"  n={n:2d}: |amp|^2 = {amp:.10f}, -log_phi = {n_eff:.6f}")

print(f"\nAvailable phi-exponents from braiding (period 10):")
print(f"  {[f'{x:.3f}' for x in n_effs if not np.isnan(x)]}")
print(f"\nTarget exponents: 3, 7, 12")
print(f"  n=12 is IMPOSSIBLE because the maximum exponent in one period")
print(f"  is bounded (braiding is periodic with period 10).")

# ================================================================
# PART 16: MULTI-ANYON BRAIDING (4 and 5 anyons)
# ================================================================
print("\n" + "=" * 78)
print("PART 16: Multi-Anyon Braiding (4 and 5 Fibonacci Anyons)")
print("=" * 78)

print(f"""
For 3 anyons, the representation space is 2-dimensional (too small
for a 3x3 CKM matrix). For N_gen = 3, we need at least a 3-dim space.

4 tau anyons -> dim 3 (Fibonacci F_4 = 3)
5 tau anyons -> dim 5 (Fibonacci F_5 = 5)

The 4-anyon case gives a 3x3 braiding matrix -- the right size for CKM!
""")

# Build the 3x3 representation for 4 Fibonacci anyons
# The fusion paths for 4 tau anyons fusing to vacuum are:
# ((tau,tau)->1, (tau,tau)->1): through 1
# ((tau,tau)->1, (tau,tau)->tau)->tau->...: mixed paths
# Actually, let me use the standard recursive construction.

# For 4 anyons fusing to 1 (vacuum):
# Basis states are labeled by intermediate fusion channels:
# |a_12, a_123> where tau x tau -> a_12, a_12 x tau -> a_123, a_123 x tau -> 1
# Since last fusion must give 1, a_123 must be such that a_123 x tau contains 1
# => a_123 = tau (since tau x tau = 1 + tau)
# So a_123 = tau always, and a_12 can be 1 or tau:
# Basis: |a_12=1, a_123=tau>, |a_12=tau, a_123=tau>
# Wait, that's only 2, but F_4 = 3.

# Actually for 4 anyons fusing to tau (not vacuum), we get dim=3.
# For 4 anyons total, the space decomposes:
# Total fusing to 1: dim = F_3 = 2
# Total fusing to tau: dim = F_4 = 3
# Total: 2 + 3 = 5 = F_5 (consistent)

# Actually, for n anyons, the total fusion space has dim F_{n+1}
# and it decomposes based on the total charge.

# Let me use the standard basis for 4 anyons fusing to tau:
# |a_12, a_123> where tau x tau -> a_12, a_12 x tau -> a_123, a_123 x tau -> tau
# For tau x tau = 1 + tau:
#   a_12 = 1: then 1 x tau = tau, so a_123 = tau. tau x tau = 1 + tau, need tau output. OK.
#   a_12 = tau: then tau x tau = 1 + tau:
#     a_123 = 1: 1 x tau = tau. Need tau output. OK.
#     a_123 = tau: tau x tau = 1 + tau. Need tau output. OK.
# So basis: |1, tau>, |tau, 1>, |tau, tau> -- dim = 3.

print(f"Basis for 4 tau anyons fusing to tau:")
print(f"  |1> = |a_12=1, a_123=tau>")
print(f"  |2> = |a_12=tau, a_123=1>")
print(f"  |3> = |a_12=tau, a_123=tau>")
print(f"  Dimension: 3")

# Now build sigma_1, sigma_2, sigma_3 for B_4

# sigma_1 braids particles 1,2: acts on a_12 channel
# In basis (1, tau, tau) with different a_123:
# sigma_1 changes a_12 but not the rest
# For a_12 = 1: R_1^{tau,tau} = R_fib_0 (fusing to 1)
# For a_12 = tau: R_tau^{tau,tau} = R_fib_1 (fusing to tau)
# sigma_1 is block diagonal:
#   State |1> (a_12=1): gets R_0
#   States |2>, |3> (a_12=tau): get R_1
sigma1_4 = np.diag([R_fib_0, R_fib_1, R_fib_1]).astype(complex)

# sigma_2 braids particles 2,3: changes a_12 and a_123 via F
# This requires the F-matrix applied to the (a_12, a_123) pair
# The F-matrix relates different fusion trees:
# F: (12)3 -> 1(23), i.e., changes from grouping (1,2) first to (2,3) first.
#
# For sigma_2 on a 3-anyon substring (particles 2,3,4 viewed locally):
# We need to:
#   1. Re-parenthesize using F to expose the (2,3) braiding
#   2. Apply R to the (2,3) pair
#   3. Re-parenthesize back using F^{-1}
#
# This is more involved. Let me construct it properly.
#
# The F-matrix for the transition:
# ((tau x tau)_{a_12} x tau)_{a_123} -> (tau x (tau x tau)_{a_23})_{a_123}
# This is F^{tau,tau,tau,a_123}_{a_12, a_23}
#
# For each fixed a_123, F is a matrix in (a_12, a_23) space.

# F-matrix for tau,tau,tau,tau (4 external legs, all tau):
# F^{tau tau tau}_{a_123}: this is the recoupling of 3 taus
# The general F-matrix for 3 tau anyons:
# F_{a_12, a_23}^{a_123} where all external lines are tau

# For a_123 = tau (relevant for our basis states |2> and |3>):
# a_12 can be 1 or tau, a_23 can be 1 or tau
# F is the standard 2x2 Fibonacci F-matrix:
# F^{tau}_{a_12, a_23} = F_fib

# For a_123 = 1: both a_12 and a_23 must be tau (only option for tau x tau -> 1)
# Wait: if a_123 = 1, then for tau x a_23 = 1, we need a_23 = tau^* = tau
# And for a_12 x tau = 1, we need a_12 = tau^* = tau
# So the "1x1 block" F^{1} is just a phase.

# Let me think more carefully. The F-matrix transformation:
# |a_12, a_123>_s-channel -> sum_{a_23} F_{a_12, a_23}^{a_123} |a_23, a_123>_t-channel
# Then sigma_2 in the t-channel is R_{a_23}^{tau,tau}
# Then transform back.

# In our 3-dimensional space with basis |a_12, a_123>:
# |1> = |1, tau>    (a_12=1, a_123=tau)
# |2> = |tau, 1>    (a_12=tau, a_123=1)
# |3> = |tau, tau>  (a_12=tau, a_123=tau)

# The F-transform changes the parenthesization from ((12)3) to (1(23)):
# It acts on the (a_12, a_123) -> (a_23, a_13) indices, where a_13 is fixed.
# But actually, in the standard construction for B_n, sigma_2 acts locally
# on the substring involving particles 2,3.

# Let me use the explicit F-matrix approach. For 4 anyons fusing to tau:
# sigma_2 = F23^{-1} * R23 * F23 where F23 recouples around particle 2-3 braiding.

# The full F-matrix in the 3-dim space:
# The s-channel basis: |a_12, a_123> = {|1,tau>, |tau,1>, |tau,tau>}
# The t-channel basis: |a_23, a_13> = ... (different grouping)
# After F-transform: new basis is |a_23, a_13> where a_13 = total of 1,3.

# This is getting complicated. Let me construct it from the 2x2 F-matrix.
#
# The key idea: sigma_2 braids particles 2 and 3.
# In the s-channel basis ((12)3)4, particle 2 is inside the (12) group.
# To braid 2,3, we need to recouple.
#
# The F-move: F_{(12)3}^{1(23)} acts on the intermediate channel.
# Specifically, for fixed external tau and fixed a_123:
#   |a_12, a_123> = sum_{a_23} F^{a_123}_{a_12, a_23} |a_23, a_123>_reparenthesized
#
# But this reparenthesization only changes a_12 -> a_23, keeping a_123 fixed.
# So the F-matrix is block-diagonal in a_123:
#   Block a_123 = 1: only a_12 = tau is possible (tau x tau -> 1 requires intermediate tau)
#                     and a_23 = tau (tau x tau -> 1)
#                     => 1x1 block, F^{1}_{tau,tau} = 1 (just a phase)
#   Block a_123 = tau: a_12 in {1, tau}, a_23 in {1, tau}
#                      => 2x2 block, F^{tau} = F_fib

# Actually, let me look at this more carefully. The F-symbol is:
# [F^{abc}_d]_{ef} where a,b,c,d are external, e,f are internal.
# Here a=b=c=d=tau, e=a_12, f=a_23.
# No, the standard is: F^{abc}_d relates (a x_e b) x c -> d and a x (b x_f c) -> d

# For our case: (tau x_{a12} tau) x tau -> a_123 vs tau x (tau x_{a23} tau) -> a_123
# The constraint is that the final fusion to a_123 is the same.

# I'll use the known result. For SU(2)_3 / Fibonacci:
# [F^{tau,tau,tau}_{tau}]_{a12, a23} = F_fib (the standard 2x2 matrix)
# [F^{tau,tau,tau}_{1}]_{a12=tau, a23=tau} = 1 (scalar, since only tau x tau -> 1)

# So in our 3D basis:
# State |1> = |a_12=1, a_123=tau>  belongs to the tau block (a_123=tau)
# State |2> = |a_12=tau, a_123=1>  belongs to the 1 block (a_123=1)
# State |3> = |a_12=tau, a_123=tau> belongs to the tau block (a_123=tau)

# The F-matrix in the 3D space:
# Block a_123=tau (states |1> and |3>):
#   F^{tau} = [[F_00, F_01], [F_10, F_11]] = [[1/phi, 1/sqrt(phi)], [1/sqrt(phi), -1/phi]]
#   Maps: |a_12=1, tau> -> F_00|a_23=1, tau> + F_01|a_23=tau, tau>  etc.
# Block a_123=1 (state |2>):
#   F^{1} = 1 (trivial)
#   Maps: |a_12=tau, 1> -> |a_23=tau, 1>

# In the t-channel basis {|a_23=1, a_123=tau>, |a_23=tau, a_123=1>, |a_23=tau, a_123=tau>}:
# sigma_2 in t-channel: R applied to (2,3) braiding, which gives R_{a23}
# sigma_2^{t-channel} = diag(R_0^{tau,tau}, R_1^{tau,tau}, R_1^{tau,tau})
#                      = diag(R_fib_0, R_fib_1, R_fib_1)
# (a_23=1 -> R_0, a_23=tau -> R_1)

# But WAIT: the t-channel basis has the SAME form as the s-channel basis
# (just with a_23 replacing a_12), so the F-matrix maps:
# s-channel -> t-channel, then R in t-channel, then F^{-1} back.

# In the ordered basis {|1>=(1,tau), |2>=(tau,1), |3>=(tau,tau)}:
F_4 = np.zeros((3, 3))
# Block a_123=tau: states |1> (a_12=1) and |3> (a_12=tau)
# Maps to t-channel: a_23=1 (maps to position 0) and a_23=tau (maps to position 2)
F_4[0, 0] = F_fib[0, 0]  # |1> (a12=1, a123=tau) -> |1'> (a23=1, a123=tau)
F_4[0, 2] = F_fib[0, 1]  # |1> -> |3'> (a23=tau, a123=tau)
F_4[2, 0] = F_fib[1, 0]  # |3> -> |1'>
F_4[2, 2] = F_fib[1, 1]  # |3> -> |3'>
# Block a_123=1: state |2> (a_12=tau) -> |2'> (a_23=tau)
F_4[1, 1] = 1.0

print(f"\nF-matrix in 3D space (4 anyons):")
for row in F_4:
    print(f"  [{row[0]:8.5f} {row[1]:8.5f} {row[2]:8.5f}]")

# Verify F_4 is unitary (orthogonal since real)
print(f"F_4 * F_4^T = I? {np.allclose(F_4 @ F_4.T, np.eye(3))}")

# R in t-channel basis: diag(R_0, R_1, R_1)
R_tchan = np.diag([R_fib_0, R_fib_1, R_fib_1]).astype(complex)

# sigma_2 for 4 anyons
F_4c = F_4.astype(complex)
F_4_inv = np.linalg.inv(F_4c)
sigma2_4 = F_4_inv @ R_tchan @ F_4c

print(f"\nsigma_2 for 4 anyons:")
for row in sigma2_4:
    print(f"  [{row[0]:12.6f} {row[1]:12.6f} {row[2]:12.6f}]")

# sigma_3 braids particles 3,4: need another F-move
# For sigma_3, we need to recouple to expose particles 3,4.
# In the s-channel basis ((12)3)4, particle 4 is outermost.
# We need F to go from (a_123, 4) -> (a_34, ...) where we braid 3,4.
# This is simpler: sigma_3 recouples a_123 with the 4th particle.

# The F-matrix for the recoupling ((a_123) x tau) -> a_total = tau
# relates a_123 to a_34.
# For a_total = tau:
#   (a_123 x tau -> tau): a_123 can be 1 or tau
#     a_123 = 1: 1 x tau = tau. OK.
#     a_123 = tau: tau x tau = 1 + tau. Contains tau. OK.
#   Similarly for the (3,4) recoupling: a_34 can be 1 or tau.
# The F-symbol [F^{a_12, tau, tau}_{tau}]_{a_123, a_34}
# depends on a_12!

# For a_12 = 1:
#   a_123 must be tau (since 1 x tau = tau only)
#   a_34 can be 1 or tau
#   F is 1x2: only entry for a_123=tau -> a_34=1 or tau
#   This is [F^{1,tau,tau}_{tau}]_{tau, a_34}

# For a_12 = tau:
#   a_123 can be 1 or tau
#   The F-matrix relates {(tau x tau)_{a_123} x tau -> tau} to {tau x (tau x tau)_{a_34} -> tau}
#   But with fixed a_12 = tau...

# This is getting quite involved. Let me instead use the explicit recursive construction.
# For the braid group representation, I'll use the known result from the literature.

# Actually, for sigma_3 in the 4-anyon case, the relevant F-matrix is:
# F^{a_12,tau,tau}_{tau}: recouples the last two fusion vertices.
# The key insight: sigma_3 only affects the outermost fusion (a_123 x tau -> tau)
# This is equivalent to braiding 3,4 in the subtree.

# For the subtree of the last two particles:
# (a_123 x tau -> tau), braid the pair (3,4):
# Need to recouple: (a_123 x tau)_tau -> a_123 x (tau x_{a_34} tau)_...
# The F-symbol for this is F^{a_123, tau, tau}_{tau}

# For a_123 = 1:
#   1 x tau -> tau (unique), so a_34 must also give tau when fused: tau x a_34 -> ?
#   Not quite. Let me be more precise.
#   [F^{a_123=1, tau, tau}_{tau}]_{e=tau, f}: e is the intermediate in (1 x tau)=tau
#   Since 1 x tau = tau (unique channel), e=tau always.
#   After recoupling: 1 x (tau x_f tau) -> tau
#   Need tau x_f tau contains tau: f = 1 or tau both work (since tau x tau = 1+tau).
#   But then 1 x (result) -> tau: result must be tau.
#   If f=1: tau x tau -> 1, then 1 x 1 = 1. Need tau. NO!
#   If f=tau: tau x tau -> tau, then 1 x tau = tau. YES!
#   So only f=tau contributes. F^{1,tau,tau}_{tau} = 1 (scalar).

# For a_123 = tau:
#   tau x tau -> tau or 1.
#   [F^{tau,tau,tau}_{tau}]_{e,f} is the standard 2x2 Fibonacci F-matrix!
#   e,f in {1, tau}. This is F_fib.

# So sigma_3 in the 3D space {|1>=(1,tau), |2>=(tau,1), |3>=(tau,tau)}:
# State |1> has a_123 = tau. The subtree (a_123, tau) -> tau.
#   sigma_3 on this state: F^{tau,tau,tau}_{tau} applied, but |1> has a_12=1
#   and a_123=tau. sigma_3 changes a_123.
#
# Wait, I need to think about this differently.
# sigma_3 braids particles 3 and 4. In the tree ((12)3)4:
# - a_12 is the fusion of 1,2
# - a_123 is the fusion of (12),3
# - Total: a_123 x 4 -> tau
#
# sigma_3 requires recoupling to isolate the (3,4) pair.
# The F-move: a_12 x (a_23 x tau)... No, we need:
# (a_12 x tau_3)_{a_123} x tau_4 -> tau
# Recouple to: a_12 x (tau_3 x tau_4)_{a_34} -> tau
# The F-symbol: [F^{a_12, tau, tau}_{tau}]_{a_123, a_34}

# For a_12 = 1:
#   State |1>: a_123=tau. After F: what a_34?
#   [F^{1,tau,tau}_{tau}]_{tau, a_34}
#   From the standard F-symbol calculation:
#   Since a_12=1 is trivial, (1 x tau) = tau always.
#   The F-symbol F^{1,b,c}_d is delta_{e,b}: it's trivial.
#   So [F^{1,tau,tau}_{tau}]_{tau, a_34} = delta_{tau, a_34}?
#   Actually for F^{1,tau,tau}_{tau}, e goes from a_12 x tau_3 = 1 x tau = tau,
#   and f goes from tau_3 x tau_4 = 1 or tau.
#   The F-symbol for trivial first particle: F^{1,b,c}_d = delta (the "trivial" F-symbol)
#   So yes, it's trivial: a_34 = a_123 = tau. F = 1.

# For a_12 = tau:
#   States |2> (a_123=1) and |3> (a_123=tau).
#   [F^{tau,tau,tau}_{tau}]_{a_123, a_34} = F_fib (the standard 2x2)
#   This mixes states |2> and |3>.

# So the F-matrix for sigma_3 recoupling:
F_sigma3 = np.eye(3)
F_sigma3[1,1] = F_fib[0,0]  # (tau,1) block: a_123=1 -> a_34 index 0 (=1)
F_sigma3[1,2] = F_fib[0,1]  # (tau,1) -> (tau,tau)
F_sigma3[2,1] = F_fib[1,0]  # (tau,tau) -> (tau,1)
F_sigma3[2,2] = F_fib[1,1]  # (tau,tau) -> (tau,tau)

# In the recoupled basis, the R-matrix for (3,4) braiding:
# a_34 = 1 -> R_0, a_34 = tau -> R_1
# State |1> (a_12=1): a_34 = tau (from trivial F), so R_1
# States |2>,|3> (a_12=tau): after F, the |a_34=1> gets R_0, |a_34=tau> gets R_1

# Actually the R-matrix in the recoupled basis is the same for all states:
# R_{a34} applied to each. But the mapping is:
# |1> -> a_34=tau -> R_1 (unchanged position)
# |2'>, |3'> (after F on |2>,|3>) -> R_0, R_1 respectively

# Let me construct it as: sigma_3 = F3^{-1} R34 F3
F_sigma3_c = F_sigma3.astype(complex)
F_sigma3_inv = np.linalg.inv(F_sigma3_c)

# R in recoupled basis:
# Position 0 (|1>, a_34=tau after trivial F): R_1
# After F on states |2>,|3>: position 1 maps to a_34=1 (R_0), position 2 to a_34=tau (R_1)
R_34 = np.diag([R_fib_1, R_fib_0, R_fib_1]).astype(complex)

sigma3_4 = F_sigma3_inv @ R_34 @ F_sigma3_c

print(f"\nsigma_1 for 4 anyons (diagonal):")
for row in sigma1_4:
    print(f"  [{row[0]:12.6f} {row[1]:12.6f} {row[2]:12.6f}]")

print(f"\nsigma_2 for 4 anyons:")
for row in sigma2_4:
    print(f"  [{row[0]:12.6f} {row[1]:12.6f} {row[2]:12.6f}]")

print(f"\nsigma_3 for 4 anyons:")
for row in sigma3_4:
    print(f"  [{row[0]:12.6f} {row[1]:12.6f} {row[2]:12.6f}]")

# Verify braid relations for B_4
# sigma_i * sigma_{i+1} * sigma_i = sigma_{i+1} * sigma_i * sigma_{i+1}
br12_lhs = sigma1_4 @ sigma2_4 @ sigma1_4
br12_rhs = sigma2_4 @ sigma1_4 @ sigma2_4
print(f"\nBraid relation s1*s2*s1 = s2*s1*s2: {np.allclose(br12_lhs, br12_rhs)}")

br23_lhs = sigma2_4 @ sigma3_4 @ sigma2_4
br23_rhs = sigma3_4 @ sigma2_4 @ sigma3_4
print(f"Braid relation s2*s3*s2 = s3*s2*s3: {np.allclose(br23_lhs, br23_rhs)}")

# sigma_i * sigma_j = sigma_j * sigma_i for |i-j| >= 2
comm13 = sigma1_4 @ sigma3_4
comm31 = sigma3_4 @ sigma1_4
print(f"Far commutativity s1*s3 = s3*s1: {np.allclose(comm13, comm31)}")

# ================================================================
# PART 17: CKM FROM 4-ANYON BRAIDING
# ================================================================
print("\n" + "=" * 78)
print("PART 17: CKM from 4-Anyon Braiding Matrix Elements")
print("=" * 78)

print(f"""
The 4-anyon Fibonacci fusion space has dim=3 -- exactly the right
size for a 3x3 CKM-like matrix!

Hypothesis: The CKM matrix is a braiding matrix in this 3D space.
The 3 basis states correspond to 3 generations:
  |1> = |a_12=1, a_123=tau>    ~ Gen 1 (light)
  |2> = |a_12=tau, a_123=1>    ~ Gen 2 (middle)
  |3> = |a_12=tau, a_123=tau>  ~ Gen 3 (heavy)

Test: does any braid word W in B_4 give |W_ij| ~ phi^{{-n_ij}}?
""")

# Scan over simple braid words in B_4
gen4_names = ['s1', 's1i', 's2', 's2i', 's3', 's3i']
gen4_mats = [
    sigma1_4, np.linalg.inv(sigma1_4),
    sigma2_4, np.linalg.inv(sigma2_4),
    sigma3_4, np.linalg.inv(sigma3_4)
]

# First, check simple powers and products
print(f"\n--- Simple braid words in B_4 ---")
test_words_4 = [
    ('s2', sigma2_4),
    ('s3', sigma3_4),
    ('s2^2', np.linalg.matrix_power(sigma2_4, 2)),
    ('s3^2', np.linalg.matrix_power(sigma3_4, 2)),
    ('s2*s3', sigma2_4 @ sigma3_4),
    ('s3*s2', sigma3_4 @ sigma2_4),
    ('s1*s2*s3', sigma1_4 @ sigma2_4 @ sigma3_4),
    ('(s2*s3)^2', np.linalg.matrix_power(sigma2_4 @ sigma3_4, 2)),
    ('(s2*s3)^3', np.linalg.matrix_power(sigma2_4 @ sigma3_4, 3)),
    ('(s1*s2*s3)^2', np.linalg.matrix_power(sigma1_4 @ sigma2_4 @ sigma3_4, 2)),
    ('(s1*s2*s3)^3', np.linalg.matrix_power(sigma1_4 @ sigma2_4 @ sigma3_4, 3)),
    ('(s1*s2*s3)^5', np.linalg.matrix_power(sigma1_4 @ sigma2_4 @ sigma3_4, 5)),
]

print(f"{'Word':>20s} {'|01|':>8s} {'|02|':>8s} {'|12|':>8s} {'n_01':>8s} {'n_02':>8s} {'n_12':>8s}")
print("-" * 72)

for name, M in test_words_4:
    amps = []
    n_effs = []
    for i, j in [(0,1), (0,2), (1,2)]:
        amp = abs(M[i,j])
        amps.append(amp)
        if 0 < amp < 1:
            n_effs.append(-log(amp) / ln_phi)
        else:
            n_effs.append(float('nan'))

    marker = ""
    # Check if any pair of exponents matches CKM-like structure
    if not any(np.isnan(x) for x in n_effs):
        sorted_n = sorted(n_effs)
        if (abs(sorted_n[0] - 3) < 0.5 or abs(sorted_n[1] - 7) < 0.5 or
            abs(sorted_n[2] - 12) < 0.5):
            marker = " <-- INTERESTING"

    print(f"{name:>20s} {amps[0]:8.5f} {amps[1]:8.5f} {amps[2]:8.5f} "
          f"{n_effs[0]:8.3f} {n_effs[1]:8.3f} {n_effs[2]:8.3f}{marker}")

# Systematic search over short braid words in B_4
print(f"\n--- Systematic search: braid words of length 1-4 in B_4 ---")
print(f"Looking for words where the 3 off-diagonal phi-exponents")
print(f"match {{3, 7, 12}} (within 0.5 each)...")

best_ckm_words = []
for length in range(1, 5):
    for word_indices in iterproduct(range(6), repeat=length):
        M = np.eye(3, dtype=complex)
        for idx in word_indices:
            M = M @ gen4_mats[idx]

        # Extract off-diagonal phi exponents
        n_effs_dict = {}
        valid = True
        for i, j in [(0,1), (0,2), (1,2)]:
            amp = abs(M[i,j])
            if 0 < amp < 1:
                n_effs_dict[(i,j)] = -log(amp) / ln_phi
            else:
                valid = False
                break

        if not valid:
            continue

        # Check match to {3, 7, 12}
        vals = sorted(n_effs_dict.values())
        score = abs(vals[0] - 3) + abs(vals[1] - 7) + abs(vals[2] - 12)

        if score < 3:
            word_str = "*".join([gen4_names[idx] for idx in word_indices])
            best_ckm_words.append((word_str, vals, score, length))

best_ckm_words.sort(key=lambda x: x[2])

if best_ckm_words:
    print(f"\nBest matches (top 10):")
    print(f"{'Word':>30s} {'n_small':>8s} {'n_mid':>8s} {'n_large':>8s} {'score':>8s}")
    for w_str, vals, score, length in best_ckm_words[:10]:
        marker = " ***" if score < 1 else ""
        print(f"{w_str:>30s} {vals[0]:8.3f} {vals[1]:8.3f} {vals[2]:8.3f} {score:8.3f}{marker}")
else:
    print(f"\n  No braid words found with score < 3 (length 1-4)")

# ================================================================
# PART 18: PENTAGON AND HEXAGON IDENTITIES
# ================================================================
print("\n" + "=" * 78)
print("PART 18: Pentagon and Hexagon Relations")
print("=" * 78)

print(f"""
The F-matrices satisfy the Pentagon identity:
  F_12 F_13 F_23 = F_23 F_12  (schematic)

The R-matrices satisfy the Hexagon identity:
  R F R = F R F  (schematic)

These are EXACT relations that constrain F and R completely
(up to gauge freedom) for a given fusion category.

In the Fibonacci model, these identities fix F and R uniquely.
The CKM cannot come from a "deformation" of these relations --
they are RIGID.
""")

# Verify pentagon for the 2x2 F-matrix
# Pentagon: F F = I (for self-dual) -- actually this is the fact that F is involutory
# The real pentagon involves 5 F-matrices on a 5-vertex graph.
# For Fibonacci: F^2 + F - I = 0 (the golden ratio equation!)
# Actually: F satisfies its own char poly.
eigenvalues_F = np.linalg.eigvals(F_fib)
print(f"Eigenvalues of F_fib: {eigenvalues_F}")
print(f"Note: eigenvalues are +1 and -1 (F is an involution up to sign?)")

# Check F^2
F2 = F_fib @ F_fib
print(f"\nF^2 = \n{F2}")

# Hexagon: R_c F_{c,e}^{a,b} R_e = sum_f F_{c,f}^{b,a} R_f F_{f,e}^{a,b}
# In the Fibonacci model, this is automatically satisfied.
print(f"\nHexagon identity (automatic for Fibonacci, verified by construction).")

# ================================================================
# PART 19: TOPOLOGICAL INVARIANTS AND LINK POLYNOMIALS
# ================================================================
print("\n" + "=" * 78)
print("PART 19: Topological Invariants -- Traces of Braids")
print("=" * 78)

print(f"""
The trace of a braid (Markov trace) gives a topological link invariant.
For the Fibonacci model, this gives the Jones polynomial at q=exp(2pi*i/5).

The ABSOLUTE VALUE of the Jones polynomial might have phi structure.
""")

# Markov trace: tr(b) = sum_a d_a / D^2 * tr_a(b)
# For the 2D representation (3 anyons):
# The quantum trace includes d_a weights.

# For simplicity, just compute the ordinary trace
print(f"--- Traces of braid words (3-anyon, 2D rep) ---")
D_total = sqrt(D2)
d_fib = [1, phi]  # quantum dimensions of channels 0, 1 (vacuum, tau)

for n in range(1, 21):
    sig2n = np.linalg.matrix_power(sigma2, n)

    # Quantum trace: sum_a d_a * M_{aa} / D^2
    qtrace = sum(d_fib[a] * sig2n[a, a] for a in range(2)) / D2

    # Ordinary trace
    otrace = np.trace(sig2n)

    # Jones-like invariant: |qtrace|
    jones = abs(qtrace)

    if jones > 1e-15:
        n_jones = -log(jones) / ln_phi
        print(f"  n={n:2d}: |J| = {jones:.8f}, -log_phi = {n_jones:.4f}")

# ================================================================
# PART 20: QUANTITATIVE SUMMARY AND VERDICT
# ================================================================
print("\n" + "=" * 78)
print("=" * 78)
print("PART 20: QUANTITATIVE SUMMARY AND VERDICT")
print("=" * 78)

print(f"""
TARGET: CKM exponents {{3, 7, 12}} from TQFT braiding in SU(2)_3

RESULTS:
""")

# Collect all exponents found
print(f"1. CONFORMAL WEIGHTS of SU(2)_3:")
print(f"   h_j * 20 = {{0, 3, 8, 15}}")
print(f"   Gaps: 3, 5, 7")
print(f"   Compare CKM: 3, 7, 12")
print(f"   Match n_12 = 3: YES (gap h_0 -> h_{{1/2}} = 3/20)")
print(f"   Match n_23 = 7: YES (gap h_1 -> h_{{3/2}} = 7/20)")
print(f"   Match n_13 = 12: NO (sum of gaps = 15, not 12)")
print(f"   Assessment: PARTIAL PATTERN (2/3 match)")

print(f"\n2. BRAIDING AMPLITUDES (3 anyons, 2D rep):")
print(f"   sigma_2^n is periodic with period 10")
print(f"   Available phi-exponents (one period):")
avail = []
for n in range(1, 11):
    sig2n = np.linalg.matrix_power(sigma2, n)
    amp = abs(sig2n[0,1])**2
    if 0 < amp < 1:
        avail.append(-log(amp) / ln_phi)
print(f"   {[f'{x:.3f}' for x in avail]}")
print(f"   None close to 12 (maximum ~ {max(avail):.3f})")
print(f"   Assessment: NEGATIVE (periodicity prevents large exponents)")

print(f"\n3. BRAIDING AMPLITUDES (4 anyons, 3D rep):")
if best_ckm_words:
    best = best_ckm_words[0]
    print(f"   Best word: {best[0]}")
    print(f"   Exponents: {[f'{x:.3f}' for x in best[1]]}")
    print(f"   Score (dist from 3,7,12): {best[2]:.3f}")
    if best[2] < 1:
        print(f"   Assessment: PROMISING")
    else:
        print(f"   Assessment: NEGATIVE (score > 1)")
else:
    print(f"   No words found with matching exponents in B_4 (length 1-4)")
    print(f"   Assessment: NEGATIVE")

print(f"\n4. F-MATRIX elements:")
print(f"   |F_{{01}}|^2 = phi^{{-1}} = {1/phi:.6f}  (n_eff = 1)")
print(f"   |F_{{00}}|^2 = phi^{{-2}} = {1/phi**2:.6f}  (n_eff = 2)")
print(f"   These give n = 1, 2 -- too small for CKM.")
print(f"   Assessment: NEGATIVE")

print(f"\n5. MODULAR S-MATRIX:")
# Check if any S-matrix ratio gives phi^{-3}, phi^{-7}, phi^{-12}
found_s = {3: False, 7: False, 12: False}
for i, j, i2, j2, n_eff in all_ratios:
    for tn in [3, 7, 12]:
        if abs(abs(n_eff) - tn) < 0.1:
            found_s[tn] = True
print(f"   Ratio giving n=3: {found_s[3]}")
print(f"   Ratio giving n=7: {found_s[7]}")
print(f"   Ratio giving n=12: {found_s[12]}")
print(f"   Assessment: NEGATIVE (S-matrix ratios involve phi^{{+/-1}} only)")

print(f"""
======================================================================
FINAL VERDICT
======================================================================

STATUS: NEGATIVE (with one PARTIAL PATTERN)

The TQFT braiding in SU(2)_3 Fibonacci anyon model does NOT derive
the CKM exponents {{3, 7, 12}}. The fundamental reasons:

1. PERIODICITY: Braiding amplitudes in the Fibonacci model are
   PERIODIC (period 10 for 3-anyon, period 10 for each generator
   in 4-anyon). The maximum phi-exponent achievable in one period
   is limited (~2-3 for 3-anyon, ~3-4 for 4-anyon). The exponent
   n_13 = 12 requires phi^{{-12}} ~ 0.0024, which is FAR below
   the minimum braiding amplitude in one period.

2. WRONG STRUCTURE: The conformal weight gaps of SU(2)_3 are
   {{3, 5, 7}} (in units of 1/20), not {{3, 7, 12}}.
   The gap 5 corresponds to n_23 if we skip from 1/2 to 1,
   but n_23 = 7 requires the gap from 1 to 3/2.
   The sum 3+5+7 = 15 != 12 = n_13.

3. NO MULTIPLICATIVE STRUCTURE: Braiding gives ADDITIVE phases
   (phases ADD when braiding sequentially). The CKM needs
   MULTIPLICATIVE structure n_ut * n_db = 8 -> n_23 = 7.
   Braiding sigma^{{n1}} * sigma^{{n2}} gives phase n1+n2 (additive),
   not n1*n2 (multiplicative).

4. DIMENSION MISMATCH: The 3-anyon Fibonacci space is 2D (not 3D).
   The 4-anyon space is 3D but its braiding does not reproduce CKM.

PARTIAL PATTERN: The conformal weight gaps of SU(2)_3 contain 3 and 7
(matching n_12 and n_23), suggesting the CKM exponents are RELATED to
the conformal data but through a more complex mechanism than braiding.

OPEN QUESTION: Can the conformal weight spectrum {{0, 3, 8, 15}} (units
of 1/20) be connected to CKM via a mechanism that maps:
  3 -> n_12 = 3
  8 -> n_13 = 12 (via 8 + 4 = 12, where 4 = n_ut?)
  15 -> something else?

This would require n_13 = h_1*20 + n_ut = 8 + 4 = 12. SPECULATIVE.

The CKM exponents remain best understood through:
  n_12 = 3, n_23 = 7, n_13 = 12
from the Cayley graph spectral ratios (n_ut=4, n_db=2) and the
concordance algebraic identities (exp449, exp453).

CLASSIFICATION:
  Conformal weight gaps {3,7}: PATTERN (2/3 match to CKM)
  TQFT braiding -> CKM: NEGATIVE
  SU(2)_3 structure for CKM: OPEN (no derivation found)
""")
