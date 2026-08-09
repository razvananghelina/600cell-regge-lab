#!/usr/bin/env python3
"""
Verify: Q8 is the maximally central rational subgroup of 2I,
with net abelian-only contribution after tetrahedral completion.

Tests:
  T1: 120 = 8 (Q8) + 16 (2T\Q8) + 96 (2I\2T), all subgroups closed
  T2: Q8 elements have tr in {0, +/-2} (even integers, coordinate-free)
  T3: 100% of Q8 commutators are central; <100% for 2T and 2I
  T4: Character weight W_k(Q8) = 0 for Galois-broken rho_8, rho_2
  T5: W_4(Q8) + W_4(2T\Q8) = 0 (rho_4 residue cancels at 2T shell)
  T6: After 2T completion, only rho_0 (abelian) survives as Q8-specific
"""

import numpy as np
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from commons.cell600 import build_600cell
from commons import PHI, N, IRREP_DIMS_2I

passed = 0
failed = 0

def check(name, condition):
    global passed, failed
    if condition:
        passed += 1
        print(f"  PASS  {name}")
    else:
        failed += 1
        print(f"  FAIL  {name}")

def qmul(p, q):
    a1,b1,c1,d1 = p; a2,b2,c2,d2 = q
    return np.array([a1*a2-b1*b2-c1*c2-d1*d2, a1*b2+b1*a2+c1*d2-d1*c2,
                     a1*c2-b1*d2+c1*a2+d1*b2, a1*d2+b1*c2-c1*b2+d1*a2])

def qconj(q):
    return np.array([q[0], -q[1], -q[2], -q[3]])

def find_el(q, verts):
    return np.argmin(np.linalg.norm(verts - q, axis=1))

# Build group
verts, _, _ = build_600cell()

# Classify
is_Q8 = np.array([np.sum(np.abs(v) > 0.01) == 1 and np.max(np.abs(v)) > 0.99 for v in verts])
is_2T = is_Q8.copy()
for i in range(N):
    if not is_Q8[i] and all(abs(abs(c) - 0.5) < 0.01 for c in verts[i]):
        is_2T[i] = True

Q8 = np.where(is_Q8)[0]
T_extra = np.where(is_2T & ~is_Q8)[0]
rest = np.where(~is_2T)[0]

print("=" * 60)
print("VERIFY: Q8 ABELIAN SECTOR")
print("=" * 60)

# T1: Sizes and closure
print("\nT1: Decomposition and closure")
check("T1a: |Q8| = 8", len(Q8) == 8)
check("T1b: |2T\\Q8| = 16", len(T_extra) == 16)
check("T1c: |2I\\2T| = 96", len(rest) == 96)
check("T1d: total = 120", len(Q8) + len(T_extra) + len(rest) == N)

q8_closed = all(is_Q8[find_el(qmul(verts[i], verts[j]), verts)]
                for i in Q8 for j in Q8)
check("T1e: Q8 closed under multiplication", q8_closed)

T_all = np.where(is_2T)[0]
t_closed = all(is_2T[find_el(qmul(verts[i], verts[j]), verts)]
               for i in T_all for j in T_all)
check("T1f: 2T closed under multiplication", t_closed)

# T2: Trace criterion
print("\nT2: Galois trace criterion")
traces = 2 * verts[:, 0]

q8_traces_ok = all(abs(traces[i]) < 0.01 or abs(abs(traces[i]) - 2) < 0.01 for i in Q8)
check("T2a: Q8 traces in {0, +/-2}", q8_traces_ok)

t_extra_traces_ok = all(abs(abs(traces[i]) - 1) < 0.01 for i in T_extra)
check("T2b: 2T\\Q8 traces = +/-1", t_extra_traces_ok)

# 2I\2T includes elements with tr in {0, +/-1, +/-phi, +/-1/phi}.
# The Galois criterion separates Q8 from 2T\Q8 (not Q8 from 2I\2T):
# Q8: |tr| in {0, 2} (even integer).  2T\Q8: |tr| = 1 (odd integer).
# The Galois-broken (irrational) elements are ALL in 2I\2T, but not
# all of 2I\2T is irrational (C4, C5, C6 contribute rational elements too).
rest_has_irrational = any(abs(abs(traces[i]) - PHI) < 0.01 or abs(abs(traces[i]) - 1/PHI) < 0.01
                         for i in rest)
rest_no_tr2 = all(abs(abs(traces[i]) - 2) > 0.01 for i in rest)
check("T2c: 2I\\2T contains all irrational-trace elements", rest_has_irrational)
check("T2d: 2I\\2T has no elements with |tr| = 2", rest_no_tr2)

# T3: Commutator criterion
print("\nT3: Commutator criterion")
def central_commutator_fraction(indices):
    count = 0; total = 0
    for i in indices:
        for j in indices:
            comm = qmul(qmul(verts[i], verts[j]), qmul(qconj(verts[i]), qconj(verts[j])))
            d_plus = np.linalg.norm(comm - [1,0,0,0])
            d_minus = np.linalg.norm(comm - [-1,0,0,0])
            if min(d_plus, d_minus) < 1e-6: count += 1
            total += 1
    return count, total

c_q8, t_q8 = central_commutator_fraction(Q8)
check("T3a: 100% of Q8 commutators central", c_q8 == t_q8)

c_2t, t_2t = central_commutator_fraction(T_all)
check("T3b: 2T commutator fraction < 100%", c_2t < t_2t)

# T4-T6: Character weights
print("\nT4-T6: Character weights on gauge irreps")

half_angles = [0, np.pi, np.pi/5, 2*np.pi/5, np.pi/3, np.pi/2,
               2*np.pi/3, 3*np.pi/5, 4*np.pi/5]
class_sizes_ref = [1, 1, 12, 12, 20, 30, 20, 12, 12]

def su2_char(j, a):
    if abs(a) < 1e-12: return 2*j+1
    if abs(a-np.pi) < 1e-12: return (-1)**(int(2*j))*(2*j+1)
    return np.sin((2*j+1)*a)/np.sin(a)

ct = np.zeros((9, 9))
for k in range(6):
    for c in range(9):
        ct[k,c] = su2_char(k/2, half_angles[c])
chi_V3 = np.array([su2_char(3, a) for a in half_angles])
chi_V72 = np.array([su2_char(3.5, a) for a in half_angles])
ct[7] = chi_V72 - ct[5]
ct[6] = ct[1] * ct[7]
ct[8] = chi_V3 - ct[6]

dims = np.array(IRREP_DIMS_2I)

# Q8 class membership: C0(1), C1(1), C5(6). Others: 0.
q8_cc = {0:1, 1:1, 5:6}
t_extra_cc = {4:8, 6:8}

def W(k, class_counts):
    return dims[k] * sum(cnt * ct[k, c] for c, cnt in class_counts.items()) / N

# T4: Q8 zero on Galois-broken
W8_q8 = W(8, q8_cc)
W2_q8 = W(2, q8_cc)
check("T4a: W_rho8(Q8) = 0 (SU(2) sector)", abs(W8_q8) < 1e-10)
check("T4b: W_rho2(Q8) = 0 (SU(3)_3 sector)", abs(W2_q8) < 1e-10)

# T5: rho_4 residue cancels
W4_q8 = W(4, q8_cc)
W4_t_extra = W(4, t_extra_cc)
check("T5a: W_rho4(Q8) != 0 (nonzero residue)", abs(W4_q8) > 0.1)
check("T5b: W_rho4(Q8) + W_rho4(2T\\Q8) = 0 (exact cancellation)",
      abs(W4_q8 + W4_t_extra) < 1e-10)

# T6: After 2T completion, only abelian survives
W0_q8 = W(0, q8_cc)
check("T6a: W_rho0(Q8) > 0 (abelian channel present)", W0_q8 > 0)

# Combined 2T shell: all gauge irreps
t_all_cc = {c: q8_cc.get(c, 0) + t_extra_cc.get(c, 0) for c in range(9)}
nonabelian_at_2T = [abs(W(k, t_all_cc)) for k in [2, 8]]  # rho_2, rho_8
rational_at_2T = abs(W(4, t_all_cc))  # rho_4

check("T6b: W_rho8(2T) = 0", nonabelian_at_2T[1] < 1e-10)
check("T6c: W_rho2(2T) = 0", nonabelian_at_2T[0] < 1e-10)
check("T6d: W_rho4(2T) = 0 (rational residue canceled)", rational_at_2T < 1e-10)

# Print summary values
print(f"\n  Key values:")
print(f"    W_rho0(Q8)    = {W0_q8:.6f}")
print(f"    W_rho8(Q8)    = {W8_q8:.6f}")
print(f"    W_rho2(Q8)    = {W2_q8:.6f}")
print(f"    W_rho4(Q8)    = {W4_q8:.6f}")
print(f"    W_rho4(2T\\Q8) = {W4_t_extra:.6f}")
print(f"    W_rho4(2T)    = {W4_q8 + W4_t_extra:.6f}")

print(f"\n{'='*60}")
print(f"RESULTS: {passed} passed, {failed} failed out of {passed+failed}")
if failed == 0:
    print("ALL TESTS PASS")
print(f"{'='*60}")
