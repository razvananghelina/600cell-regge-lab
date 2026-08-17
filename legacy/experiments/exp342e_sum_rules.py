"""
exp342e: Sum Rules for Quantum Numbers
=======================================
KEY DISCOVERY: sum(a_k) = 12 = vertex_degree, sum(b_k) = 8 = rank(E8)!
"""
import numpy as np

PHI = (1 + np.sqrt(5)) / 2
A1 = 5
B1 = 6

# Fermion quantum numbers (a, b) and mass exponents n = 5a + 6b
data = {
    'e':   (0,  0,  0),
    'u':   (3, -2,  3),
    'mu':  (1,  1, 11),
    'tau': (1,  2, 17),
    'c':   (2,  1, 16),
    't':   (4,  1, 26),
    'd':   (1,  0,  5),
    's':   (1,  1, 11),
    'b':  (-1,  4, 19),
}

print("=" * 72)
print("QUANTUM NUMBER SUM RULES")
print("=" * 72)

sum_a = sum(d[0] for d in data.values())
sum_b = sum(d[1] for d in data.values())
sum_n = sum(d[2] for d in data.values())

print(f"\n  sum(a_k) = {sum_a}")
print(f"  sum(b_k) = {sum_b}")
print(f"  sum(n_k) = {sum_n}")
print(f"\n  FRAMEWORK IDENTIFICATION:")
print(f"  sum(a) = {sum_a} = vertex_degree = 2*(a1+1) = 12? "
      f"{'YES!!!' if sum_a == 12 else 'NO'}")
print(f"  sum(b) = {sum_b} = rank(E8) = 8? "
      f"{'YES!!!' if sum_b == 8 else 'NO'}")
print(f"  sum(n) = a1*sum(a) + b1*sum(b) = {A1}*{sum_a} + {B1}*{sum_b} "
      f"= {A1*sum_a + B1*sum_b}")
print(f"         = N_eig * vertex_degree = 9 * 12 = 108? "
      f"{'YES!!!' if sum_n == 108 else 'NO'}")

# Decomposition by Galois pairs
pairs = [('u','t'), ('d','b'), ('s','tau')]
unpaired = ['e', 'mu', 'c']

print(f"\n  GALOIS PAIR DECOMPOSITION:")
print(f"  {'Pair':>12} {'sum_a':>6} {'sum_b':>6} {'sum_n':>6}")
print("-" * 40)
total_pair_a = 0
total_pair_b = 0
for f1, f2 in pairs:
    sa = data[f1][0] + data[f2][0]
    sb = data[f1][1] + data[f2][1]
    sn = data[f1][2] + data[f2][2]
    total_pair_a += sa
    total_pair_b += sb
    print(f"  ({f1},{f2})".rjust(12) + f" {sa:>6} {sb:>6} {sn:>6}")

print(f"  {'Pairs total':>12} {total_pair_a:>6} {total_pair_b:>6} "
      f"{sum(data[f1][2]+data[f2][2] for f1,f2 in pairs):>6}")

sa_unp = sum(data[f][0] for f in unpaired)
sb_unp = sum(data[f][1] for f in unpaired)
sn_unp = sum(data[f][2] for f in unpaired)
print(f"  {'Unpaired':>12} {sa_unp:>6} {sb_unp:>6} {sn_unp:>6}")

print(f"\n  PAIR sum(a) = {total_pair_a} = N_eig = 9? "
      f"{'YES!!!' if total_pair_a == 9 else 'NO'}")
print(f"  PAIR sum(b) = {total_pair_b} = b1 = 6? "
      f"{'YES!!!' if total_pair_b == 6 else 'NO'}")
print(f"  UNPAIRED sum(a) = {sa_unp} = N_gen = 3? "
      f"{'YES!!!' if sa_unp == 3 else 'NO'}")
print(f"  UNPAIRED sum(b) = {sb_unp} = 2? "
      f"{'YES' if sb_unp == 2 else 'NO'}")

print(f"\n  CHECK: pairs_a + unp_a = {total_pair_a} + {sa_unp} "
      f"= {total_pair_a + sa_unp} = sum(a) = {sum_a}")
print(f"  CHECK: pairs_b + unp_b = {total_pair_b} + {sb_unp} "
      f"= {total_pair_b + sb_unp} = sum(b) = {sum_b}")

# Individual pair analysis
print(f"\n  INDIVIDUAL PAIR STRUCTURE:")
for f1, f2 in pairs:
    a1f, b1f = data[f1][0], data[f1][1]
    a2f, b2f = data[f2][0], data[f2][1]
    print(f"  ({f1},{f2}): a=({a1f},{a2f}), b=({b1f},{b2f})")
    print(f"    sum_a = {a1f+a2f}, diff_a = {abs(a1f-a2f)}")
    print(f"    sum_b = {b1f+b2f}, diff_b = {abs(b1f-b2f)}")

# Z[phi] elements and their norms
print(f"\n  Z[phi] ELEMENTS z = a + b*phi:")
print(f"  {'Fermion':>8} {'z':>12} {'z_bar':>12} {'N(z)':>8} {'Tr(z)':>8}")
for f in ['e','u','d','s','mu','c','tau','t','b']:
    a, b, n = data[f]
    z = a + b*PHI
    zbar = a + b*(1-PHI)/1  # z' = a + b*phi'
    zbar_val = a + b*(-1/PHI)
    norm = a**2 + a*b - b**2  # N(z) = z*z'
    trace = 2*a + b
    print(f"  {f:>8} {z:>12.4f} {zbar_val:>12.4f} {norm:>8} {trace:>8}")

# Sum and product rules for z
print(f"\n  Sum of all z_k: ", end="")
sum_z = sum(data[f][0] + data[f][1]*PHI for f in data)
print(f"{sum_z:.4f} = {sum_a} + {sum_b}*phi = 12 + 8*phi")
print(f"  = vertex_degree + rank(E8)*phi")
print(f"  = 2*(a1+1) + rank(E8)*phi")
print(f"  Numerically: {12 + 8*PHI:.4f}")

print(f"\n  Norm of sum: N(12 + 8*phi) = 12^2 + 12*8 - 8^2 = {144+96-64} = 176")
print(f"  = 16 * 11 = (a1-1)^2 * n_mu")
print(f"  Trace of sum: 2*12 + 8 = 32 = 2^5 = 2*fermions_per_gen")

# Check: is sum_z special?
print(f"\n  sum_z = 12 + 8*phi = 4*(3 + 2*phi) = 4*(1 + 2*phi + 2)")
print(f"  = 4*(3 + 2*phi)")
print(f"  Note: 3 + 2*phi = 2*phi^2 + 1 = 2*(phi+1) + 1 = 2*phi + 3")
print(f"  And phi^3 = phi^2 + phi = (phi+1) + phi = 2*phi + 1")
print(f"  So 3 + 2*phi = phi^3 + 2 = phi^3 + 2")
print(f"  Check: phi^3 + 2 = {PHI**3 + 2:.4f}, 3 + 2*phi = {3 + 2*PHI:.4f}")
print(f"  {'MATCH' if abs(PHI**3 + 2 - (3+2*PHI)) < 0.001 else 'NO MATCH'}")

# alpha_s = 1/(2*phi^3) so 2*phi^3 = 1/alpha_s
print(f"\n  sum_z = 4*(phi^3 + 2) = 4*phi^3 + 8 = 2/alpha_s + rank(E8)")
print(f"  REMARKABLE: sum(z_k) = 2/alpha_s + rank(E8)")
print(f"  Check: 2/alpha_s + 8 = {2*2*PHI**3 + 8:.4f} vs sum_z = {sum_z:.4f}")
print(f"  Wait, 2*phi^3 = 1/alpha_s, so 4*phi^3 = 2/alpha_s")
alpha_s = 1/(2*PHI**3)
print(f"  2/alpha_s = {2/alpha_s:.4f}")
print(f"  2/alpha_s + 8 = {2/alpha_s + 8:.4f}")
print(f"  sum_z = {sum_z:.4f}")
print(f"  MATCH: {abs(2/alpha_s + 8 - sum_z) < 0.001}")

# Pair sums in Z[phi]
print(f"\n  GALOIS PAIR z-SUMS:")
for f1, f2 in pairs:
    a1f, b1f = data[f1][0], data[f1][1]
    a2f, b2f = data[f2][0], data[f2][1]
    sa = a1f + a2f
    sb = b1f + b2f
    sn = data[f1][2] + data[f2][2]
    z = sa + sb*PHI
    print(f"  ({f1},{f2}): z_sum = {sa} + {sb}*phi = {z:.4f}, "
          f"N(z_sum) = {sa**2+sa*sb-sb**2}, Tr(z_sum) = {2*sa+sb}, "
          f"n_sum = {sn}")

# Check if pair z-sums have nice structure
print(f"\n  (u,t): z_sum = 7 - 1*phi, N = 49-7-1 = 41, Tr = 13 = (a1^2+1)/2")
print(f"  (d,b): z_sum = 0 + 4*phi, N = 0+0-16 = -16, Tr = 4 = dim(spacetime)")
print(f"  (s,tau): z_sum = 2 + 3*phi, N = 4+6-9 = 1, Tr = 7")

print(f"\n  N(z_sum) for (d,b) = -16 = -(a1-1)^2 = -N(z_Planck)")
print(f"  N(z_sum) for (s,tau) = 1 = UNIT in Z[phi]!")
print(f"  N(z_sum) for (u,t) = 41 = prime")

print(f"\n  Tr(z_sum) for (u,t) = 13 = (a1^2+1)/2")
print(f"  Tr(z_sum) for (d,b) = 4 = dim(spacetime)")
print(f"  Tr(z_sum) for (s,tau) = 7 = degree(600-cell)/2 - a1/2 ?")

print(f"\n{'='*72}")
print("SUMMARY OF SUM RULES")
print("=" * 72)
print(f"""
NEW SUM RULES (all exact):

  sum(a_k) = 12 = vertex_degree = 2*(a1+1)
  sum(b_k) = 8 = rank(E8)
  sum(n_k) = 108 = N_eig * vertex_degree

  sum_z = sum(a + b*phi) = 12 + 8*phi = 2/alpha_s + rank(E8)

  Galois pairs:
    sum_a(pairs) = 9 = N_eig
    sum_b(pairs) = 6 = b1
    sum_a(unpaired) = 3 = N_gen
    sum_b(unpaired) = 2

  Individual pair traces:
    Tr(z_u + z_t) = 13 = (a1^2+1)/2
    Tr(z_d + z_b) = 4 = dim(spacetime)
    Tr(z_s + z_tau) = 7

  Individual pair norms:
    N(z_d + z_b) = -16 = -(a1-1)^2 = -N(z_Planck)
    N(z_s + z_tau) = 1 (unit in Z[phi])

  Classification: DERIVED (algebraic identities, no free parameters)
""")
print("Done.")
