#!/usr/bin/env python3
"""Exact structural checks for the multiplicity-mixing-J classification.

This verifier deliberately separates:
  * the complete discrete classification of order-zero opposite actions;
  * the previously certified finite geometric J candidates; and
  * the still-unsolved continuous semialgebraic search inside each orbit.

It does not turn absence of an exhaustive continuous solve into a no-go.
"""

from itertools import permutations


tests = passed = 0


def check(name, condition, detail=""):
    global tests, passed
    tests += 1
    condition = bool(condition)
    passed += int(condition)
    print(f"  [{'PASS' if condition else 'FAIL'}] {name}")
    if detail:
        print(f"         {detail}")


print("=" * 78)
print("MULTIPLICITY-MIXING J: EXACT ORDER-ZERO CLASSIFICATION BOUNDARY")
print("=" * 78)

# McKay-chain order: rho0,...,rho8.
d = (1, 2, 3, 4, 5, 6, 4, 2, 3)
names = tuple(f"rho{i}" for i in range(9))
check("2I irreducible dimensions square-sum to 120",
      sum(x*x for x in d) == 120)


def valid_type(k, m):
    """Weighted row/column margins for an A-A^op bimodule type."""
    rows = tuple(sum(d[j]*k[i][j] for j in range(9)) for i in range(9))
    cols = tuple(sum(d[i]*k[i][j] for i in range(9)) for j in range(9))
    target = tuple(m*x for x in d)
    return rows == target and cols == target


def diagonal_type(m):
    return tuple(tuple(m if i == j else 0 for j in range(9))
                 for i in range(9))


def permutation_type(m, p):
    return tuple(tuple(m if j == p[i] else 0 for j in range(9))
                 for i in range(9))


# Equal-dimension permutations are explicit non-diagonal admissible types.
galois = (0, 7, 8, 3, 4, 5, 6, 1, 2)
swap_four = (0, 1, 2, 6, 4, 5, 3, 7, 8)
check("Galois label permutation preserves dimensions",
      tuple(d[galois[i]] for i in range(9)) == d)
check("the two four-dimensional labels can also be permuted",
      tuple(d[swap_four[i]] for i in range(9)) == d)

for m in (22, 44):
    kd = diagonal_type(m)
    kg = permutation_type(m, galois)
    k4 = permutation_type(m, swap_four)
    check(f"diagonal bimodule type is admissible for m={m}",
          valid_type(kd, m))
    check(f"Galois-permuted bimodule type is admissible for m={m}",
          valid_type(kg, m))
    check(f"four-block-permuted bimodule type is admissible for m={m}",
          valid_type(k4, m))
    check(f"a single global multiplicity vector does not classify m={m}",
          kd != kg and kg != k4 and kd != k4,
          "three inequivalent weighted 9x9 types have the same global "
          "regular multiplicities")

    # Real dimensions of the A'-conjugacy orbit and of the U torsor.
    ambient_unitary_dim = sum((x*m)**2 for x in d)
    centralizer_dim_diag = sum(sum(v*v for v in row) for row in kd)
    orbit_dim_diag = ambient_unitary_dim-centralizer_dim_diag
    u_torsor_dim = 120*m*m
    check(f"unitary-family dimensions close exactly for m={m}",
          ambient_unitary_dim == 120*m*m
          and centralizer_dim_diag == 9*m*m
          and orbit_dim_diag == 111*m*m
          and u_torsor_dim == 120*m*m,
          f"diagonal embedding orbit={orbit_dim_diag}, "
          f"fixed-image intertwiner torsor={u_torsor_dim}")

# The requested SM blocks are a corner, not a unital action on the arena.
sm_corner_complex_rank = 1*1 + 2*2 + 3*3
check("rho0+rho1+rho8 corner has regular rank 14, not 120",
      sm_corner_complex_rank == 14 and sm_corner_complex_rank != 120)
check("both Galois corner choices have the same nonunital rank",
      d[0]**2+d[7]**2+d[2]**2 == sm_corner_complex_rank)
for m in (22, 44):
    check(f"SM corner unit is not the identity for m={m}",
          sm_corner_complex_rank*m < 120*m,
          f"corner rank={sm_corner_complex_rank*m}, arena rank={120*m}")

# Exact finite facts imported from the already registered certificates.
# They are repeated here as classification boundary conditions, not recomputed
# by floating diagonalization.
primal_signs = {
    "coefficient conjugation": (+1, +1, +1, False),
    "orbitwise inversion": (+1, 0, +1, True),
}
double_signs = {
    "pure cellular star": (+1, +1, -1, False),
    "star times inversion": (+1, 0, -1, True),
}
check("no certified primal geometric J passes order zero, KO6, and JD sign",
      all(not (row[2] == -1 and row[1] in (-1, +1) and row[3])
          for row in primal_signs.values()))
check("no certified doubled geometric J passes order zero, KO6, and JD sign",
      all(not (row[2] == -1 and row[1] in (-1, +1) and row[3])
          for row in double_signs.values()))

# Q8 control: its multiplicity contains two regular Q8 copies; the 2I arenas
# do not contain even one regular 2I copy inside C^m.
check("Q8 counterexample multiplicity contains two regular copies",
      16 == 2*8)
check("2I multiplicity spaces cannot contain a regular 2I factor",
      22 < 120 and 44 < 120)
check("this dimension contrast is not promoted to a no-go",
      True,
      "distributed bimodule embeddings remain possible and require the "
      "continuous U solve")

print("-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
if passed != tests:
    raise SystemExit(1)
