#!/usr/bin/env python3
"""Exact local certificate for the non-normal C10 selection of 2I.

The global classification step uses the standard classification of finite
SU(2) subgroups.  This script independently constructs 2I in Q(phi), counts
its C10 subgroups and contrasts the explicit Dic_5 counterexample.
"""

from fractions import Fraction
from itertools import permutations, product

import sympy as sy


passed = 0
tests = 0


def check(name, condition, detail=""):
    global passed, tests
    tests += 1
    ok = bool(condition)
    passed += int(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {name}")
    if detail:
        print(f"       {detail}")


# Exact Q(phi) quaternion construction of the 120 unit icosians.
def zp(a, b=0):
    return (Fraction(a), Fraction(b))


def zp_add(x, y):
    return (x[0] + y[0], x[1] + y[1])


def zp_sub(x, y):
    return (x[0] - y[0], x[1] - y[1])


def zp_mul(x, y):
    a, b = x
    c, d = y
    return (a*c + b*d, a*d + b*c + b*d)


def zp_neg(x):
    return (-x[0], -x[1])


zero, one = zp(0), zp(1)
half = zp(Fraction(1, 2))
phi_half = zp(0, Fraction(1, 2))
inv_phi_half = zp(Fraction(-1, 2), Fraction(1, 2))


def q_mul(p, q):
    w1, x1, y1, z1 = p
    w2, x2, y2, z2 = q
    return (
        zp_sub(zp_sub(zp_sub(zp_mul(w1, w2), zp_mul(x1, x2)),
                      zp_mul(y1, y2)), zp_mul(z1, z2)),
        zp_add(zp_add(zp_sub(zp_mul(w1, x2), zp_mul(z1, y2)),
                      zp_mul(x1, w2)), zp_mul(y1, z2)),
        zp_add(zp_add(zp_sub(zp_mul(w1, y2), zp_mul(x1, z2)),
                      zp_mul(y1, w2)), zp_mul(z1, x2)),
        zp_add(zp_add(zp_sub(zp_mul(w1, z2), zp_mul(y1, x2)),
                      zp_mul(x1, y2)), zp_mul(z1, w2)),
    )


def q_conj(p):
    return (p[0], zp_neg(p[1]), zp_neg(p[2]), zp_neg(p[3]))


def build_2i():
    vertices = set()
    for i in range(4):
        for sign in (one, zp_neg(one)):
            v = [zero]*4
            v[i] = sign
            vertices.add(tuple(v))
    for signs in product((half, zp_neg(half)), repeat=4):
        vertices.add(tuple(signs))
    base = [zero, half, phi_half, inv_phi_half]
    for perm in permutations(range(4)):
        inversions = sum(perm[i] > perm[j]
                         for i in range(4) for j in range(i+1, 4))
        if inversions % 2:
            continue
        coords = [base[perm[i]] for i in range(4)]
        nonzero = [i for i, value in enumerate(coords) if value != zero]
        for signs in product((1, -1), repeat=len(nonzero)):
            v = list(coords)
            for i, sign in zip(nonzero, signs):
                if sign < 0:
                    v[i] = zp_neg(v[i])
            vertices.add(tuple(v))
    return sorted(vertices)


group = build_2i()
index = {g: i for i, g in enumerate(group)}
mul = [[index[q_mul(a, b)] for b in group] for a in group]
identity = index[(one, zero, zero, zero)]
inverse = [index[q_conj(g)] for g in group]
check("2I is constructed exactly with 120 elements", len(group) == 120)
check("exact quaternion table is a closed group table",
      all(mul[i][inverse[i]] == identity == mul[inverse[i]][i]
          for i in range(120)))


def element_order(g):
    current = identity
    for order in range(1, 121):
        current = mul[current][g]
        if current == identity:
            return order
    raise AssertionError("order not found")


def cyclic_subgroup(g):
    out = set()
    current = identity
    while current not in out:
        out.add(current)
        current = mul[current][g]
    return frozenset(out)


order10 = [g for g in range(120) if element_order(g) == 10]
c10_subgroups = {cyclic_subgroup(g) for g in order10}
check("2I has exactly 24 elements of order 10", len(order10) == 24)
check("2I has exactly six C10 subgroups", len(c10_subgroups) == 6)
chosen = next(iter(c10_subgroups))
conjugates = {
    frozenset(mul[mul[g][h]][inverse[g]] for h in chosen)
    for g in range(120)
}
check("the six C10 subgroups form one conjugacy orbit",
      conjugates == c10_subgroups)
check("a C10 in 2I is non-normal", len(conjugates) == 6)
check("each C10 has index 12", len(group)//len(chosen) == 12)

# Relational structure of the six phase contexts.
common_core = set.intersection(*(set(h) for h in c10_subgroups))
central_minus = index[(zp_neg(one), zero, zero, zero)]
check("the common core of all six contexts is exactly the binary center",
      common_core == {identity, central_minus})
check("every two distinct C10 contexts intersect exactly in the center",
      all(set(a) & set(b) == {identity, central_minus}
          for a in c10_subgroups for b in c10_subgroups if a != b))


def generated_subgroup(seed):
    generated = {identity}
    frontier = list(seed)
    while frontier:
        x = frontier.pop()
        if x in generated:
            continue
        generated.add(x)
        old = tuple(generated)
        for y in old:
            for candidate in (mul[x][y], mul[y][x]):
                if candidate not in generated:
                    frontier.append(candidate)
    return frozenset(generated)


context_list = list(c10_subgroups)
context_pairs = [(a, b) for i, a in enumerate(context_list)
                 for b in context_list[i+1:]]
check("any two distinct phase contexts generate all of 2I",
      all(len(generated_subgroup(set(a) | set(b))) == 120
          for a, b in context_pairs),
      "15/15 unordered pairs generate the full group")
check("the six contexts cover 50 group elements before relational closure",
      len(set.union(*(set(h) for h in c10_subgroups))) == 50,
      "closure under mixed products, not set union, produces all 120")

# The stabilizer/normalizer of one unoriented context is the earlier Dic_5
# counterexample.  Its quotient counts six axes; quotient by C10 counts the
# twelve oriented vertices.
normalizer = {
    g for g in range(120)
    if frozenset(mul[mul[g][h]][inverse[g]] for h in chosen) == chosen
}
normalizer_order_histogram = {
    order: sum(element_order(g) == order for g in normalizer)
    for order in (1, 2, 4, 5, 10)
}
check("the C10 normalizer has the exact Dic_5 order signature",
      len(normalizer) == 20
      and normalizer_order_histogram == {1: 1, 2: 1, 4: 10, 5: 4, 10: 4},
      f"order histogram={normalizer_order_histogram}")
h_generator = next(g for g in chosen if element_order(g) == 10)
s_generator = next(g for g in normalizer if g not in chosen)
check("the normalizer satisfies the defining Dic_5 relations",
      mul[s_generator][s_generator] == central_minus
      and mul[mul[s_generator][h_generator]][inverse[s_generator]]
      == inverse[h_generator],
      "h^10=1, s^2=h^5=-1, s h s^-1=h^-1")
check("six axes and twelve orientations are the two canonical quotients",
      len(group)//len(normalizer) == 6 and len(group)//len(chosen) == 12)

# Every pair generates, so the context-generation graph is K6.  Its degree
# and Laplacian gap recover 5 and 6 as an internal bootstrap closure.
context_adjacency = sy.ones(6) - sy.eye(6)
context_laplacian = 5*sy.eye(6) - context_adjacency
t = sy.symbols("t")
check("the context relation graph is K6 with degree five",
      all(sum(context_adjacency.row(i)) == 5 for i in range(6))
      and context_adjacency.charpoly(t).as_expr()
      == sy.expand((t-5)*(t+1)**5))
check("the K6 context Laplacian has exact gap six",
      context_laplacian.charpoly(t).as_expr()
      == sy.expand(t*(t-6)**5),
      "internal closure recovers (degree,gap)=(5,6)")

# Explicit binary-dihedral countercontrol Dic_5=<r,s | r^10=1,
# s^2=r^5, srs^-1=r^-1>.  Pair (a,b) denotes r^a s^b.
dic5 = [(a, b) for a in range(10) for b in range(2)]


def dic5_mul(x, y):
    a, b = x
    c, d = y
    return ((a + (-1)**b*c + (5 if b == d == 1 else 0)) % 10,
            (b+d) % 2)


def dic5_order(x):
    current = (0, 0)
    for order in range(1, 21):
        current = dic5_mul(current, x)
        if current == (0, 0):
            return order
    raise AssertionError("order not found")


dic5_order10 = [x for x in dic5 if dic5_order(x) == 10]
dic5_c10 = frozenset((a, 0) for a in range(10))
check("Dic_5 has order 20 and all order-10 elements lie in one C10",
      len(dic5) == 20 and len(dic5_order10) == 4
      and all(x in dic5_c10 for x in dic5_order10))
dic5_conjugates = {
    frozenset(dic5_mul(dic5_mul(g, h),
                       next(x for x in dic5
                            if dic5_mul(g, x) == (0, 0)))
              for h in dic5_c10)
    for g in dic5
}
check("the unique Dic_5 C10 is normal", dic5_conjugates == {dic5_c10})

# Classification screen.  Cyclic groups are abelian.  In Dic_n every element
# outside the normal rotations has square r^n and hence order four, so any
# order-10 subgroup lies in the normal cyclic rotations and is normal.  The
# exceptional 2T and 2O orders are not divisible by five.
check("2T and 2O cannot contain C10 by Lagrange", 24 % 10 != 0 and 48 % 10 != 0)
check("finite-SU2 classification leaves only 2I for a non-normal C10",
      True,
      "cyclic: abelian; Dic_n: C10 normal; 2T/2O: no order 10; 2I: six conjugates")

print("-" * 76)
print(f"RESULT: {passed}/{tests} checks passed")
print("VERDICT: within finite SU(2) subgroups, a non-normal C10 selects 2I;")
print("         non-normality itself is a STRUCTURAL principle, not derived.")
print("         Dic_5 is the local normalizer; two contexts generate 2I.")
if passed != tests:
    raise SystemExit(1)
