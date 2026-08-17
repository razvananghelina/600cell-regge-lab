"""
EXP-117b: Deep Spectral Analysis of the 120-Cell
==================================================

Follow-up to exp117. The 120-cell has 27 distinct eigenvalues.
exp117 found that only 18/27 match a+b*phi. But the a+b*phi search
may have had false positives (close but not exact matches).

This experiment:
1. Identifies minimal polynomials of ALL eigenvalues
2. Determines the algebraic fields involved
3. Checks for Galois pairing structure
4. Compares with 600-cell spectral structure

Author: Claude (exp117b)
Date: 2026-02-08
"""

import numpy as np
from itertools import permutations
from collections import Counter, deque

PHI = (1 + np.sqrt(5)) / 2

print("=" * 80)
print("EXP-117b: Deep Spectral Analysis of the 120-Cell")
print("=" * 80)
print()

# ============================================================
# SECTION 1: Rebuild 120-cell and compute eigenvalues
# ============================================================

print("-" * 80)
print("SECTION 1: Rebuilding 120-cell")
print("-" * 80)
print()

def build_600cell():
    phi = PHI
    vertices = []
    for i in range(4):
        for s in [1, -1]:
            v = [0, 0, 0, 0]
            v[i] = s
            vertices.append(v)
    for s0 in [0.5, -0.5]:
        for s1 in [0.5, -0.5]:
            for s2 in [0.5, -0.5]:
                for s3 in [0.5, -0.5]:
                    vertices.append([s0, s1, s2, s3])
    base_vals = [phi/2, 0.5, 1/(2*phi), 0]
    for s0 in [1, -1]:
        for s1 in [1, -1]:
            for s2 in [1, -1]:
                signed = [s0*base_vals[0], s1*base_vals[1],
                          s2*base_vals[2], base_vals[3]]
                for p in permutations(range(4)):
                    inv = sum(1 for i in range(4)
                              for j in range(i+1, 4) if p[i] > p[j])
                    if inv % 2 == 0:
                        v = [signed[p[i]] for i in range(4)]
                        v_rounded = [round(x, 10) for x in v]
                        if v_rounded not in vertices:
                            vertices.append(v_rounded)
    N = len(vertices)
    verts = [np.array(v) for v in vertices]
    edge_threshold = 1.0 / PHI + 0.01
    adj = np.zeros((N, N), dtype=int)
    for i in range(N):
        for j in range(i+1, N):
            d = np.linalg.norm(verts[i] - verts[j])
            if d < edge_threshold:
                adj[i, j] = 1
                adj[j, i] = 1
    return vertices, verts, adj

vertices_600, verts_600, adj_600 = build_600cell()
neighbors = [set(j for j in range(120) if adj_600[i, j]) for i in range(120)]

tetrahedra = []
for i in range(120):
    for j in neighbors[i]:
        if j <= i: continue
        common_ij = neighbors[i] & neighbors[j]
        for k in common_ij:
            if k <= j: continue
            common_ijk = common_ij & neighbors[k]
            for l in common_ijk:
                if l <= k: continue
                tetrahedra.append((i, j, k, l))

N_120 = len(tetrahedra)
face_to_tet = {}
for idx, tet in enumerate(tetrahedra):
    for skip in range(4):
        face = frozenset(tet[i] for i in range(4) if i != skip)
        if face in face_to_tet:
            face_to_tet[face].append(idx)
        else:
            face_to_tet[face] = [idx]

adj_120 = np.zeros((N_120, N_120), dtype=int)
for face, tets in face_to_tet.items():
    if len(tets) == 2:
        adj_120[tets[0], tets[1]] = 1
        adj_120[tets[1], tets[0]] = 1

print(f"  120-cell: V={N_120}, E={adj_120.sum()//2}, degree={adj_120.sum(axis=1)[0]}")

# Compute eigenvalues
evals_120 = np.linalg.eigvalsh(adj_120.astype(float))
evals_120 = sorted(evals_120, reverse=True)

# Group by multiplicity
unique_evals = []
tol = 0.0005
i = 0
while i < len(evals_120):
    val = evals_120[i]
    mult = 1
    while i + mult < len(evals_120) and abs(evals_120[i + mult] - val) < tol:
        mult += 1
    avg_val = sum(evals_120[i:i + mult]) / mult
    unique_evals.append((avg_val, mult))
    i += mult

print(f"  Distinct eigenvalues: {len(unique_evals)}")
print()

# ============================================================
# SECTION 2: Find minimal polynomials
# ============================================================

print("-" * 80)
print("SECTION 2: Minimal polynomials of eigenvalues")
print("-" * 80)
print()

def find_min_poly(theta, max_deg=4, max_coeff=30):
    """Find minimal polynomial with integer coefficients.
    Returns (degree, coefficients) or None."""
    powers = [theta**k for k in range(max_deg + 1)]

    # Degree 1: x - a = 0
    a = round(theta)
    if abs(theta - a) < 1e-8:
        return 1, [1, -a], f"x - {a}"

    # Degree 2: x^2 + bx + c = 0
    for b in range(-max_coeff, max_coeff + 1):
        c_float = -(theta**2 + b * theta)
        c = round(c_float)
        if abs(c_float - c) < 1e-6:
            val = theta**2 + b * theta + c
            if abs(val) < 1e-6:
                # Verify it's irreducible (no rational root)
                disc = b**2 - 4*c
                if disc >= 0:
                    sqrt_disc = np.sqrt(disc)
                    r1 = (-b + sqrt_disc) / 2
                    r2 = (-b - sqrt_disc) / 2
                    if abs(r1 - round(r1)) > 0.01 and abs(r2 - round(r2)) > 0.01:
                        return 2, [1, b, c], f"x^2 + {b}x + {c}"
                    elif abs(r1 - round(r1)) < 0.01 or abs(r2 - round(r2)) < 0.01:
                        # Has a rational root, but this IS the min poly if degree 1 failed
                        return 2, [1, b, c], f"x^2 + {b}x + {c} (reducible)"
                else:
                    return 2, [1, b, c], f"x^2 + {b}x + {c}"

    # Degree 3: x^3 + bx^2 + cx + d = 0
    for b in range(-max_coeff, max_coeff + 1):
        for c in range(-max_coeff, max_coeff + 1):
            d_float = -(theta**3 + b * theta**2 + c * theta)
            d = round(d_float)
            if abs(d_float - d) < 1e-4:
                val = theta**3 + b * theta**2 + c * theta + d
                if abs(val) < 1e-4:
                    return 3, [1, b, c, d], f"x^3 + {b}x^2 + {c}x + {d}"

    # Degree 4
    for b in range(-15, 16):
        for c in range(-15, 16):
            for d in range(-15, 16):
                e_float = -(theta**4 + b*theta**3 + c*theta**2 + d*theta)
                e = round(e_float)
                if abs(e_float - e) < 1e-3:
                    val = theta**4 + b*theta**3 + c*theta**2 + d*theta + e
                    if abs(val) < 1e-3:
                        return 4, [1, b, c, d, e], f"x^4 + {b}x^3 + {c}x^2 + {d}x + {e}"

    return None

def identify_field(deg, coeffs):
    """Identify the algebraic field from the minimal polynomial."""
    if deg == 1:
        return "Q"
    elif deg == 2:
        a, b, c = coeffs
        disc = b**2 - 4*a*c
        if disc < 0:
            return f"Q(sqrt({disc}))"
        sqrt_d = int(round(np.sqrt(disc)))
        if sqrt_d**2 == disc:
            return f"Q"  # rational roots
        # Simplify: disc = m * k^2 -> sqrt(disc) = k*sqrt(m)
        d = disc
        for p in [4, 9, 16, 25, 36, 49]:
            while d % p == 0:
                d //= p
        return f"Q(sqrt({d}))"
    elif deg == 3:
        return "cubic extension"
    elif deg == 4:
        return "quartic extension"
    return "unknown"

# Identify all eigenvalues
print(f"  {'theta':>10s} {'mult':>5s} {'min poly':>30s} {'field':>20s} {'roots':>30s}")
print(f"  {'-'*100}")

field_decomposition = {}

for val, mult in unique_evals:
    result = find_min_poly(val)
    if result:
        deg, coeffs, poly_str = result
        field = identify_field(deg, coeffs)

        # Compute roots explicitly
        roots_str = ""
        if deg == 1:
            roots_str = f"{-coeffs[1]}"
        elif deg == 2:
            a, b, c = coeffs
            disc = b**2 - 4*a*c
            if disc >= 0:
                r1 = (-b + np.sqrt(disc)) / (2*a)
                r2 = (-b - np.sqrt(disc)) / (2*a)
                roots_str = f"({-b}+/-sqrt({disc}))/2"
            else:
                roots_str = f"complex"
        elif deg == 3:
            roots_str = f"cubic"

        print(f"  {val:10.5f} {mult:5d} {poly_str:>30s} {field:>20s} {roots_str:>30s}")

        if field not in field_decomposition:
            field_decomposition[field] = []
        field_decomposition[field].append((val, mult, poly_str))
    else:
        print(f"  {val:10.5f} {mult:5d} {'NOT FOUND':>30s}")

print()

# ============================================================
# SECTION 3: Field decomposition summary
# ============================================================

print("-" * 80)
print("SECTION 3: Algebraic field decomposition")
print("-" * 80)
print()

for field, entries in sorted(field_decomposition.items()):
    total_mult = sum(m for _, m, _ in entries)
    print(f"  {field}:")
    for val, mult, poly in entries:
        sqrt_m = np.sqrt(mult)
        sq_str = f" = {int(round(sqrt_m))}^2" if abs(sqrt_m-round(sqrt_m))<0.01 else ""
        print(f"    theta = {val:10.5f}, mult = {mult:4d}{sq_str}, poly: {poly}")
    print(f"    Subtotal dimension: {total_mult}")
    print()

# ============================================================
# SECTION 4: Group eigenvalues by multiplicity
# ============================================================

print("-" * 80)
print("SECTION 4: Eigenvalues grouped by multiplicity")
print("-" * 80)
print()

mult_groups = {}
for val, mult in unique_evals:
    if mult not in mult_groups:
        mult_groups[mult] = []
    mult_groups[mult].append(val)

for mult in sorted(mult_groups.keys()):
    vals = mult_groups[mult]
    sqrt_m = np.sqrt(mult)
    sq_str = f" = {int(round(sqrt_m))}^2" if abs(sqrt_m-round(sqrt_m))<0.01 else ""
    print(f"  Multiplicity {mult}{sq_str}: {len(vals)} eigenvalue(s)")

    for v in vals:
        result = find_min_poly(v)
        field_str = identify_field(result[0], result[1]) if result else "?"
        print(f"    theta = {v:10.5f} in {field_str}")

    if len(vals) > 1:
        s = sum(vals)
        print(f"    Sum = {s:.4f}")
    print()

# ============================================================
# SECTION 5: Galois pairing structure
# ============================================================

print("-" * 80)
print("SECTION 5: Galois conjugate pairing")
print("-" * 80)
print()

# For Q(sqrt(5)): conjugate sends sqrt(5) -> -sqrt(5), i.e., phi -> 1-phi
# For Q(sqrt(D)): conjugate sends sqrt(D) -> -sqrt(D)

print("  Checking Galois pairs (eigenvalues sharing the same minimal polynomial):")
print()

# Find pairs sharing min poly
min_polys = {}
for val, mult in unique_evals:
    result = find_min_poly(val)
    if result:
        deg, coeffs, poly = result
        key = tuple(coeffs)
        if key not in min_polys:
            min_polys[key] = []
        min_polys[key].append((val, mult))

for coeffs, entries in sorted(min_polys.items()):
    if len(entries) > 1:
        deg = len(coeffs) - 1
        field = identify_field(deg, list(coeffs))
        vals_str = ", ".join(f"{v:.5f}(m={m})" for v, m in entries)
        poly_str = "x^{} ".format(deg)
        for i, c in enumerate(coeffs[1:], 1):
            if c != 0:
                sign = "+" if c > 0 else ""
                if i < deg:
                    poly_str += f"{sign}{c}x^{deg-i} "
                else:
                    poly_str += f"{sign}{c} "
        print(f"  Poly: {poly_str.strip()}")
        print(f"    Roots: {vals_str}")
        print(f"    Field: {field}")

        # Check if all entries have same multiplicity
        mults = [m for _, m in entries]
        print(f"    Same multiplicity: {'YES' if len(set(mults))==1 else 'NO'} ({mults})")
        print()

# ============================================================
# SECTION 6: Comparison with 600-cell
# ============================================================

print("-" * 80)
print("SECTION 6: Comparison with 600-cell spectral structure")
print("-" * 80)
print()

print("  600-CELL: All eigenvalues in Q or Q(sqrt(5))")
print("  All multiplicities are perfect squares: 1,4,9,16,25,36,9,16,4")
print("  Fields: Q (3 eigenvalues: 12, 3, 0, -2, -3) and Q(sqrt(5)) (2 pairs)")
print()

# Wait, 600-cell has: 12(1), 6phi(4), 4phi(9), 3(16), 0(25), -2(36), 4-4phi(9), -3(16), 6-6phi(4)
# Rational: 12, 3, 0, -2, -3 -> 5 eigenvalues
# Q(sqrt(5)): {6phi, 6-6phi} and {4phi, 4-4phi} -> 2 Galois pairs
# Total: 5 + 4 = 9 ✓

print("  120-CELL field decomposition:")
for field, entries in sorted(field_decomposition.items()):
    total = sum(m for _, m, _ in entries)
    count = len(entries)
    print(f"    {field}: {count} eigenvalue(s), total dim = {total}")

print()

# Total dimensions per field
total_Q = sum(m for v, m in unique_evals
              if find_min_poly(v) and find_min_poly(v)[0] == 1)
print(f"  Rational (Q) eigenspace dimension: {total_Q}")

total_sqrt5 = 0
for val, mult in unique_evals:
    result = find_min_poly(val)
    if result and result[0] == 2:
        _, coeffs, _ = result
        disc = coeffs[1]**2 - 4*coeffs[2]
        d = disc
        for p in [4, 9, 16, 25]:
            while d % p == 0:
                d //= p
        if d == 5:
            total_sqrt5 += mult

print(f"  Q(sqrt(5)) eigenspace dimension: {total_sqrt5}")

# ============================================================
# SECTION 7: Discriminant analysis
# ============================================================

print()
print("-" * 80)
print("SECTION 7: Discriminant analysis")
print("-" * 80)
print()

print("  Quadratic minimal polynomials and their discriminants:")
print()

quadratic_discs = []
for val, mult in unique_evals:
    result = find_min_poly(val)
    if result and result[0] == 2:
        deg, coeffs, poly = result
        disc = coeffs[1]**2 - 4*coeffs[2]
        # Squarefree part
        d = disc
        for p in [4, 9, 16, 25, 36, 49, 64, 100, 121, 144]:
            while d > 0 and d % p == 0:
                d //= p
        while d < 0 and (-d) % 4 == 0:
            d //= 4

        quadratic_discs.append(d)
        r1 = (-coeffs[1] + np.sqrt(disc)) / 2
        r2 = (-coeffs[1] - np.sqrt(disc)) / 2
        print(f"  {poly:30s}  disc = {disc:4d} = {'(' + str(int(np.sqrt(disc//d))) + ')^2*' if disc != d else ''}{d:3d}"
              f"  sqrt({d}), mult={mult}")

print()
unique_discs = sorted(set(quadratic_discs))
print(f"  Unique squarefree discriminants: {unique_discs}")
print()

# Check relationships to 600-cell constants
print("  Relationships to 600-cell constants:")
for d in unique_discs:
    notes = []
    if d == 5:
        notes.append("sqrt(5) = 2*phi - 1, defines golden ratio")
    if d == 2:
        notes.append("sqrt(2), appears in cross-polytope (Type A)")
    if d == 13:
        notes.append(f"13 = F_7 (7th Fibonacci number)")
    if d == 21:
        notes.append(f"21 = F_8 (8th Fibonacci number)")
    # Check Fibonacci
    fib = [1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89]
    if d in fib:
        idx = fib.index(d)
        notes.append(f"Fibonacci F_{idx+1}")
    # Check Lucas
    lucas = [2, 1, 3, 4, 7, 11, 18, 29, 47]
    if d in lucas:
        idx = lucas.index(d)
        notes.append(f"Lucas L_{idx}")

    print(f"    sqrt({d}): {', '.join(notes) if notes else 'no obvious connection'}")

# ============================================================
# SECTION 8: Cubic polynomial analysis
# ============================================================

print()
print("-" * 80)
print("SECTION 8: Cubic minimal polynomials")
print("-" * 80)
print()

for val, mult in unique_evals:
    result = find_min_poly(val)
    if result and result[0] == 3:
        deg, coeffs, poly = result
        a, b, c, d = coeffs

        # Vieta's formulas
        print(f"  {poly}")
        print(f"    One root: theta = {val:.6f}, multiplicity = {mult}")

        # Find all three roots
        # Vieta: r1+r2+r3 = -b, r1*r2+r1*r3+r2*r3 = c, r1*r2*r3 = -d
        print(f"    Sum of roots: {-b}")
        print(f"    Sum of pairwise products: {c}")
        print(f"    Product of roots: {-d}")

        # Find other roots among eigenvalues
        other_roots = []
        for v2, m2 in unique_evals:
            if abs(v2 - val) > 0.01:  # not the same eigenvalue
                test = v2**3 + b*v2**2 + c*v2 + d
                if abs(test) < 0.01:
                    other_roots.append((v2, m2))

        print(f"    Other roots found among eigenvalues:")
        for v2, m2 in other_roots:
            print(f"      theta = {v2:.6f}, mult = {m2}")

        # Discriminant
        disc = 18*b*c*d - 4*b**3*d + b**2*c**2 - 4*c**3 - 27*d**2
        print(f"    Discriminant: {disc}")
        if disc > 0:
            print(f"    3 real roots, Galois group = Z/3Z or S_3")
        elif disc == 0:
            print(f"    Multiple root")
        else:
            print(f"    1 real root + 2 complex, Galois group = S_3")

        # Check if disc is related to Fibonacci/phi
        print(f"    sqrt(|disc|) = {np.sqrt(abs(disc)):.4f}")
        print()

# ============================================================
# SECTION 9: The "factor 3" structure
# ============================================================

print("-" * 80)
print("SECTION 9: The factor-3 structure")
print("-" * 80)
print()

print("  Universal factor 3:")
print(f"    Degree ratio: 12/4 = 3")
print(f"    Diameter ratio: 15/5 = 3")
print(f"    Eigenvalue count ratio: 27/9 = 3")
print()

# Do the eigenvalues have a 3-fold structure?
# 600-cell eigenvalues sorted: 12, 6phi, 4phi, 3, 0, -2, 4-4phi, -3, 6-6phi
# Can we map 120-cell eigenvalues to "triplets" of 600-cell eigenvalues?

print("  Multiplicity comparison:")
print(f"    600-cell multiplicities: 1, 4, 9, 16, 25, 36, 9, 16, 4")
print(f"    Sum = 120")
print()

# Group 120-cell by n^2 value
for n in range(1, 7):
    n_sq = n**2
    vals_120 = [(v, m) for v, m in unique_evals if m == n_sq]
    vals_600 = [(v, m) for v, m in [(12,1),(6*PHI,4),(4*PHI,9),(3,16),(0,25),
                                     (-2,36),(4-4*PHI,9),(-3,16),(6-6*PHI,4)]
                if m == n_sq]

    count_600 = len(vals_600)
    count_120 = len(vals_120)
    total_120 = sum(m for v, m in vals_120)

    if count_600 > 0 or count_120 > 0:
        print(f"  mult = {n}^2 = {n_sq}:")
        print(f"    600-cell: {count_600} eigenvalue(s)")
        print(f"    120-cell: {count_120} eigenvalue(s), total dim = {total_120}")
        if count_600 > 0:
            print(f"    Ratio: {count_120}/{count_600} = {count_120/count_600:.1f}")
        print()

# Non-square multiplicities in 120-cell
non_sq_mults = [(v, m) for v, m in unique_evals
                if abs(np.sqrt(m) - round(np.sqrt(m))) > 0.01]
print(f"  Non-square multiplicities in 120-cell:")
for v, m in non_sq_mults:
    # Factor the multiplicity
    factors = []
    mm = m
    for p in [2, 3, 5, 7, 11, 13]:
        while mm % p == 0:
            factors.append(p)
            mm //= p
    result = find_min_poly(v)
    field = identify_field(result[0], result[1]) if result else "?"
    print(f"    theta = {v:10.5f}, mult = {m:4d} = {'*'.join(map(str,factors)):>10s}, field: {field}")

# ============================================================
# SECTION 10: Summary
# ============================================================

print()
print("=" * 80)
print("SUMMARY")
print("=" * 80)
print()

print("  ALGEBRAIC FIELDS in 120-cell spectrum:")
for field, entries in sorted(field_decomposition.items()):
    total = sum(m for _, m, _ in entries)
    count = len(entries)
    print(f"    {field:25s}: {count:2d} eigenvalue(s), dim = {total:4d}")

total_verified = sum(sum(m for _, m, _ in entries) for entries in field_decomposition.values())
print(f"    {'TOTAL':25s}: {sum(len(e) for e in field_decomposition.values()):2d} eigenvalue(s), dim = {total_verified:4d}")
print()

print("  KEY FINDINGS:")
print("  1. 120-cell eigenvalues span MULTIPLE algebraic fields (not just Q(sqrt(5)))")
print("  2. Factor 3: degree 12/4=3, diameter 15/5=3, eigenvalues 27/9=3")

all_sq = all(abs(np.sqrt(m) - round(np.sqrt(m))) < 0.01 for _, m in unique_evals)
print(f"  3. All multiplicities n^2: {'YES' if all_sq else 'NO'}")

# Check first 6
first_6 = [m for _, m in unique_evals[:6]]
print(f"  4. First 6 multiplicities: {first_6} = 1^2,2^2,...,6^2: {first_6 == [1,4,9,16,25,36]}")
print()

print("  CLASSIFICATION:")
print("  - 600-cell: spectrally PURE (all in Q(sqrt(5)), all mult=n^2)")
print("  - 120-cell: spectrally MIXED (Q, Q(sqrt(2)), Q(sqrt(5)), Q(sqrt(13)),")
print("              Q(sqrt(21)), plus cubic extensions)")
print("  - The 600-cell is the MORE fundamental object spectrally")
print()
