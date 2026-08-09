"""
exp476: Is the alpha equation a fit or a derivation?

Exhaustive scan: how many 'natural' framework equations give alpha ~ 1/137?
"""
import numpy as np
from itertools import product as iprod

phi = (1 + np.sqrt(5)) / 2
a1 = 5
b1 = 6
N = 120
h = 30
rank_E8 = 8
dim_E8 = 248
Neig = 9
alpha_exp = 1 / 137.035999084

print("="*70)
print("  exp476: IS THE ALPHA EQUATION A FIT?")
print("="*70)

# ============================================================
# TASK 3 FIRST: How many expressions land near 137?
# ============================================================
print("\n" + "="*70)
print("  TASK 3: Natural expressions near 137")
print("="*70)

expressions_near_137 = {}

# Scan: c * a1^p * b1^q * phi^r
for c in [1, 2, 3, 4, 5, 6, 8, 10, 12, 15, 20, 24, 30, 60, 120]:
    for p_a1 in range(3):  # 0, 1, 2
        for p_b1 in range(3):
            for r in range(-2, 8):
                val = c * (a1**p_a1) * (b1**p_b1) * (phi**r)
                if 136 < val < 138:
                    expr = f"{c} * a1^{p_a1} * b1^{p_b1} * phi^{r}"
                    expressions_near_137[expr] = val

# Also with pi
for c in [1, 2, 3, 4, 5, 6, 8, 10, 12]:
    for p_a1 in range(3):
        for p_b1 in range(3):
            for r in range(-2, 8):
                val = c * (a1**p_a1) * (b1**p_b1) * (phi**r) * np.pi
                if 136 < val < 138:
                    expr = f"{c} * a1^{p_a1} * b1^{p_b1} * phi^{r} * pi"
                    expressions_near_137[expr] = val

# With special constants
for c in [1, 2, 3, 4]:
    for special, sname in [(N, 'N'), (h, 'h'), (rank_E8, 'rank'), (dim_E8, 'dim_E8'), (Neig, 'Neig')]:
        for r in range(-4, 8):
            val = c * special * (phi**r)
            if 136 < val < 138:
                expr = f"{c} * {sname} * phi^{r}"
                expressions_near_137[expr] = val

# Ratios of special constants * phi^r
for num in [N, h, rank_E8, dim_E8, Neig, a1, b1]:
    for den in [a1, b1, 2, 3, 4]:
        if den == 0: continue
        for r in range(-2, 8):
            val = (num/den) * (phi**r)
            if 136 < val < 138:
                expr = f"({num}/{den}) * phi^{r}"
                expressions_near_137[expr] = val

# Sort by distance to 1/alpha
target_B = 1/alpha_exp  # 137.036
sorted_exprs = sorted(expressions_near_137.items(), key=lambda x: abs(x[1] - target_B))

print(f"\n  Expressions giving values in [136, 138]:")
print(f"  Total found: {len(expressions_near_137)}")
print(f"\n  Top 20 closest to 1/alpha = {target_B:.6f}:")
for i, (expr, val) in enumerate(sorted_exprs[:20]):
    err = abs(val - target_B) / target_B * 100
    print(f"    {i+1:2d}. {expr:40s} = {val:.6f}  (err {err:.4f}%)")

# How many within various tolerances?
within_1pct = sum(1 for _, v in expressions_near_137.items() if abs(v - target_B)/target_B < 0.01)
within_01pct = sum(1 for _, v in expressions_near_137.items() if abs(v - target_B)/target_B < 0.001)
within_001pct = sum(1 for _, v in expressions_near_137.items() if abs(v - target_B)/target_B < 0.0001)

print(f"\n  Within 1%:    {within_1pct}")
print(f"  Within 0.1%:  {within_01pct}")
print(f"  Within 0.01%: {within_001pct}")

# Where does 4*a1*phi^4 rank?
val_claimed = 4 * a1 * phi**4
err_claimed = abs(val_claimed - target_B) / target_B * 100
rank_claimed = sum(1 for _, v in expressions_near_137.items()
                   if abs(v - target_B) < abs(val_claimed - target_B)) + 1
print(f"\n  4*a1*phi^4 = {val_claimed:.6f}, error = {err_claimed:.4f}%, rank = {rank_claimed}/{len(expressions_near_137)}")

# ============================================================
# TASK 1: Exhaustive equation scan
# ============================================================
print("\n" + "="*70)
print("  TASK 1: Full equation scan A*alpha^2 - B*alpha + C = 0")
print("="*70)

# Generate A values
A_values = {}
for n in [1, 2, 3, 4, 5, 6]:
    A_values[f"{n}*pi"] = n * np.pi
A_values["a1*pi"] = a1 * np.pi
A_values["b1*pi"] = b1 * np.pi
for n in [1, 2, 3]:
    A_values[f"pi/{n}"] = np.pi / n

# Generate B values (extensive)
B_values = {}
for c in [1, 2, 3, 4, 5, 6, 8, 10, 12, 15, 20, 24, 30]:
    for p in range(3):
        for q in range(3):
            for r in range(7):
                val = c * (a1**p) * (b1**q) * (phi**r)
                if 50 < val < 500:  # reasonable range for B
                    B_values[f"{c}*a1^{p}*b1^{q}*phi^{r}"] = val

# Add special combinations
for special, sname in [(N, 'N'), (h, 'h'), (rank_E8, 'rank'), (dim_E8, 'dim_E8')]:
    for r in range(7):
        val = special * (phi**r)
        if 50 < val < 500:
            B_values[f"{sname}*phi^{r}"] = val
    for den in [a1, b1, 2, 3]:
        for r in range(7):
            val = (special/den) * (phi**r)
            if 50 < val < 500:
                B_values[f"{sname}/{den}*phi^{r}"] = val

# Generate C values
C_values = {"1": 1, "1/2": 0.5, "2": 2, "a1": a1, "1/a1": 1/a1}

total_combos = len(A_values) * len(B_values) * len(C_values)
print(f"  A values: {len(A_values)}")
print(f"  B values: {len(B_values)}")
print(f"  C values: {len(C_values)}")
print(f"  Total combinations: {total_combos}")

hits = {0.01: [], 0.001: [], 0.0001: [], 0.00001: []}

for a_name, A in A_values.items():
    for b_name, B in B_values.items():
        for c_name, C in C_values.items():
            disc = B**2 - 4*A*C
            if disc < 0:
                continue
            root1 = (B - np.sqrt(disc)) / (2*A)
            root2 = (B + np.sqrt(disc)) / (2*A)

            for root in [root1, root2]:
                if root <= 0:
                    continue
                err = abs(root - alpha_exp) / alpha_exp
                entry = (a_name, b_name, c_name, root, 1/root, err)
                for tol in hits:
                    if err < tol:
                        hits[tol].append(entry)

print(f"\n  Results:")
for tol in sorted(hits.keys()):
    print(f"    Within {tol*100:.3f}%: {len(hits[tol])} equations")

# Show the best hits
if hits[0.001]:
    print(f"\n  Top equations within 0.1% of alpha:")
    best = sorted(hits[0.001], key=lambda x: x[5])[:15]
    for a_name, b_name, c_name, root, inv_root, err in best:
        print(f"    A={a_name:10s} B={b_name:25s} C={c_name:5s} -> 1/alpha={inv_root:.4f} (err={err*100:.5f}%)")

# ============================================================
# TASK 2: Count and analyze
# ============================================================
print("\n" + "="*70)
print("  TASK 2: Hit analysis")
print("="*70)

# Among hits at 0.01% level, which A values appear?
if hits[0.0001]:
    A_dist = {}
    C_dist = {}
    for a_name, b_name, c_name, root, inv_root, err in hits[0.0001]:
        A_dist[a_name] = A_dist.get(a_name, 0) + 1
        C_dist[c_name] = C_dist.get(c_name, 0) + 1

    print(f"\n  Among {len(hits[0.0001])} equations within 0.01%:")
    print(f"  A coefficient distribution:")
    for name, count in sorted(A_dist.items(), key=lambda x: -x[1]):
        print(f"    {name}: {count}")
    print(f"  C coefficient distribution:")
    for name, count in sorted(C_dist.items(), key=lambda x: -x[1]):
        print(f"    {name}: {count}")

# ============================================================
# TASK 4: Vieta product analysis
# ============================================================
print("\n" + "="*70)
print("  TASK 4: Vieta product alpha*alpha' = C/A")
print("="*70)

vieta_values = {}
for a_name, A in A_values.items():
    for c_name, C in C_values.items():
        vp = C / A
        vieta_values[f"C/A = {c_name}/{a_name}"] = vp

# The claimed value
print(f"  Claimed: alpha*alpha' = 1/(2*pi) = {1/(2*np.pi):.10f}")
print(f"\n  All Vieta products C/A:")
for name, vp in sorted(vieta_values.items(), key=lambda x: abs(x[1] - 1/(2*np.pi))):
    dist = abs(vp - 1/(2*np.pi))
    marker = " <-- CLAIMED" if dist < 1e-15 else ""
    if vp < 1 and vp > 0.001:  # only show reasonable range
        print(f"    {name:30s} = {vp:.10f}{marker}")

# Is 1/(2*pi) forced?
# In KK: gauge coupling^2 = 1/Vol(fiber)
# For S^1 fiber of circumference 2*pi: alpha*alpha' = 1/(2*pi)
# For S^1 fiber of circumference pi: alpha*alpha' = 1/pi
# For S^1 fiber of circumference 4*pi: alpha*alpha' = 1/(4*pi)

# The Hopf fiber IS a circle of circumference 2*pi (total holonomy = 2*pi)
# This is derived from: 2*a1 edges * pi/a1 per edge = 2*pi
# The angle pi/a1 per edge is derived from cos(pi/5) = phi/2 (specific to a1=5)

print(f"\n  Is 2*pi forced?")
print(f"  Hopf fiber: 2*a1 = {2*a1} edges, each subtending pi/a1 = {np.pi/a1:.6f} rad")
print(f"  Total holonomy: 2*a1 * pi/a1 = 2*pi (EXACT, independent of a1)")
print(f"  So the Vieta product 1/(2*pi) is TOPOLOGICAL (c1 = 1 of the Hopf bundle)")

# ============================================================
# TASK 5: The verdict
# ============================================================
print("\n" + "="*70)
print("  TASK 5: HONEST VERDICT")
print("="*70)

# Count B-values near 137 that are "simple"
simple_near_137 = []
for expr, val in sorted_exprs:
    # Count the "complexity" of the expression
    parts = expr.replace("*", " ").split()
    n_parts = len([p for p in parts if p not in ['pi']])
    if abs(val - target_B)/target_B < 0.001:  # within 0.1%
        simple_near_137.append((expr, val, n_parts))

print(f"\n  B-values within 0.1% of 1/alpha ({target_B:.3f}):")
for expr, val, complexity in simple_near_137:
    err = abs(val - target_B)/target_B * 100
    print(f"    {expr:40s} = {val:.6f} (err {err:.4f}%, complexity {complexity})")

n_equations_01 = len(hits[0.001]) if hits[0.001] else 0
n_equations_001 = len(hits[0.0001]) if hits[0.0001] else 0
n_B_near = len(simple_near_137)
n_B_total = len(expressions_near_137)

print(f"""
  SUMMARY:
  --------
  Total B-expressions in [136,138]: {n_B_total}
  B-expressions within 0.1% of 1/alpha: {n_B_near}
  Full equations within 0.1% of alpha: {n_equations_01}
  Full equations within 0.01% of alpha: {n_equations_001}

  1. Is 4*a1*phi^4 unique?
     There are {n_B_total} expressions landing in [136,138].
     Within 0.1%: {n_B_near} expressions.
     4*a1*phi^4 ranks #{rank_claimed} by proximity to 1/alpha.

  2. Vieta product 1/(2*pi)?
     FORCED: The Hopf fiber holonomy is 2*pi for ANY a1
     (total angle = 2*a1 * pi/a1 = 2*pi). So A = 2*pi, C = 1
     gives alpha*alpha' = 1/(2*pi) as a TOPOLOGICAL identity.

  3. Does the equation PREDICT alpha or FIT it?
     The equation 2*pi*x^2 - B*x + 1 = 0 with B near 137 gives
     x near 1/137 TAUTOLOGICALLY (since x ~ 1/B for large B).
     The CONTENT is whether B = 4*a1*phi^4 is the UNIQUE geometric
     expression, not whether the equation works.

  VERDICT:
""")

if n_B_near <= 3:
    print("  4*a1*phi^4 is HIGHLY SELECTIVE among natural framework expressions.")
    print("  The equation is NOT a generic fit.")
elif n_B_near <= 10:
    print("  4*a1*phi^4 is MODERATELY selective - a few competitors exist.")
    print("  The equation is STRUCTURAL but not uniquely forced.")
else:
    print("  MANY expressions land near 137. The B coefficient is NOT unique.")
    print("  The equation may be a FIT dressed in geometric language.")
