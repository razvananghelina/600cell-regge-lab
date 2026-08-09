#!/usr/bin/env python3
"""Exact audit of non-central J candidates on the six C10 contexts."""

import contextlib
import io
from itertools import combinations, permutations
import runpy

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


# Reuse the exact Q(phi) construction, suppressing its already-audited output.
with contextlib.redirect_stdout(io.StringIO()):
    geometry = runpy.run_path("verify_nonnormal_c10_selection.py")
mul = geometry["mul"]
inverse = geometry["inverse"]
contexts = sorted(geometry["c10_subgroups"], key=lambda h: tuple(sorted(h)))
context_index = {h: i for i, h in enumerate(contexts)}


def compose(p, q):
    return tuple(p[q[i]] for i in range(len(p)))


def invert_perm(p):
    out = [0]*len(p)
    for i, j in enumerate(p):
        out[j] = i
    return tuple(out)


def parity(p):
    inversions = sum(p[i] > p[j]
                     for i in range(len(p)) for j in range(i+1, len(p)))
    return inversions % 2


context_action = []
for g in range(120):
    moved = []
    for h in contexts:
        conjugate = frozenset(mul[mul[g][x]][inverse[g]] for x in h)
        moved.append(context_index[conjugate])
    context_action.append(tuple(moved))
context_group = sorted(set(context_action))
check("derived context action factors to 60 permutations", len(context_group) == 60)
check("every geometrically induced context permutation is even",
      all(parity(p) == 0 for p in context_group),
      "sign restricted to simple A5 is trivial")
check("the action on six contexts is transitive",
      {p[0] for p in context_group} == set(range(6)))

# Find an exact two-generator presentation of the permutation image.
identity6 = tuple(range(6))


def generated_perms(generators):
    found = {identity6}
    frontier = list(generators)
    while frontier:
        x = frontier.pop()
        if x in found:
            continue
        old = tuple(found)
        found.add(x)
        for y in old:
            for z in (compose(x, y), compose(y, x)):
                if z not in found:
                    frontier.append(z)
    return found


generator_pair = next((a, b) for a, b in combinations(context_group, 2)
                      if len(generated_perms((a, b))) == 60)
check("two derived permutations generate the full A5 context action",
      len(generated_perms(generator_pair)) == 60)


def perm_matrix(p):
    matrix = sy.zeros(6)
    for j, i in enumerate(p):
        matrix[i, j] = 1
    return matrix


group_matrices = [perm_matrix(p) for p in context_group]
flat = sy.Matrix.hstack(*[m.reshape(36, 1) for m in group_matrices])
check("the context group-algebra image has type C plus M5",
      flat.rank() == 26,
      "permutation module 6=1+5; image dimension 1^2+5^2=26")
check("the derived group-algebra image is noncommutative",
      generator_pair[0] != compose(generator_pair[0], generator_pair[1])
      and perm_matrix(generator_pair[0])*perm_matrix(generator_pair[1])
      != perm_matrix(generator_pair[1])*perm_matrix(generator_pair[0]))

# Exhaust every possible permutation real structure P K.  Order zero for the
# full A5 image would force the two conjugate A5 generator pairs to commute
# elementwise.  No P in S6 works (hence no odd P works).
Gm = [perm_matrix(p) for p in generator_pair]
order_zero_permutations = []
for p in permutations(range(6)):
    P = perm_matrix(p)
    conjugate_generators = [P*g*P.T for g in Gm]
    if all(g*h == h*g for g in Gm for h in conjugate_generators):
        order_zero_permutations.append(p)
check("no permutation J gives order zero for the full noncommutative A5 image",
      order_zero_permutations == [],
      "exhausted all 720 permutations of the six contexts")

# Task-B counterexample: an order-zero noncommutative algebra which is not
# J-invariant.  Blocks are B11(multiplicity 2), B12(dim 2), B21(dim 2).
# J swaps the scalar copies and the two cross blocks: three transpositions,
# hence an odd involution.
odd_swap = (1, 0, 4, 5, 2, 3)
P = perm_matrix(odd_swap)
check("the structural J permutation is odd and involutive",
      parity(odd_swap) == 1 and compose(odd_swap, odd_swap) == identity6)

X = sy.Matrix([[0, 1], [1, 0]])
Z = sy.Matrix([[1, 0], [0, -1]])


def represented(lam, matrix):
    return sy.diag(lam, lam, lam, lam, 0, 0) \
        + sy.diag(0, 0, 0, 0, 1, 1) \
        * sy.diag(0, 0, 0, 0, matrix[0, 0], matrix[1, 1]) \
        + sy.Matrix([
            [0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, matrix[0, 1]],
            [0, 0, 0, 0, matrix[1, 0], 0],
        ])


# Use an explicit block constructor as a cross-check against the formula.
def rep(lam, matrix):
    out = sy.zeros(6)
    out[0, 0] = out[1, 1] = out[2, 2] = out[3, 3] = lam
    out[4:6, 4:6] = matrix
    return out


pi_scalar = rep(1, sy.zeros(2))
pi_x = rep(0, X)
pi_z = rep(0, Z)
check("C plus M2 representation is faithful, unital and noncommutative",
      rep(1, sy.eye(2)) == sy.eye(6)
      and pi_x*pi_z != pi_z*pi_x)


def opposite(matrix):
    # J pi(a)^* J^-1 = P pi(a)^T P^-1 for J=P K.
    return P*matrix.T*P.T


left_generators = (pi_scalar, pi_x, pi_z)
right_generators = tuple(opposite(a) for a in left_generators)
check("non-J-invariant C plus M2 algebra satisfies order zero",
      all(a*b == b*a for a in left_generators for b in right_generators))
matrix_units = []
for i in range(2):
    for j in range(2):
        unit = sy.zeros(2)
        unit[i, j] = 1
        matrix_units.append(rep(0, unit))
algebra_basis = (pi_scalar, *matrix_units)
algebra_flat = sy.Matrix.hstack(*[m.reshape(36, 1) for m in algebra_basis])
opposite_x_flat = sy.Matrix.hstack(
    algebra_flat, opposite(pi_x).reshape(36, 1)
)
check("the order-zero algebra is genuinely not J-invariant",
      opposite(pi_x)[2:4, 2:4] == X
      and opposite_x_flat.rank() == algebra_flat.rank() + 1,
      "the M2 block moves from B21 to the commuting B12 block")

# A grading anti-commuting with J exists only after orienting the three J
# pairs.  It is structural, not context-derived.  The design filter gives a
# rationally nondegenerate but non-unimodular intersection form.
gamma = sy.diag(1, -1, 1, 1, -1, -1)
check("structural grading has KO6 J-gamma sign",
      P*gamma*P.T == -gamma)
p_c = pi_scalar
p_m = rep(0, sy.eye(2))
projectors = (p_c, p_m)
cap = sy.Matrix([
    [sy.trace(gamma*a*opposite(b)) for b in projectors]
    for a in projectors
])
check("design filter gives antisymmetric Cap with Pfaffian two",
      cap == sy.Matrix([[0, 2], [-2, 0]])
      and cap.det() == 4,
      "rank 2 over Q, det=4=Pf^2, strict integral PD fails")

# Products pi(a) pi^o(b) are constant on the two scalar multiplicity copies,
# while gamma assigns opposite signs there.  Metric-dimension-zero
# orientability therefore fails for every Hochschild 0-cycle.
zero_cycle_span = [a*b for a in left_generators for b in right_generators]
span_with_gamma = sy.Matrix.hstack(
    *[m.reshape(36, 1) for m in zero_cycle_span], gamma.reshape(36, 1)
)
span_without_gamma = sy.Matrix.hstack(
    *[m.reshape(36, 1) for m in zero_cycle_span]
)
check("structural C plus M2 candidate fails orientability",
      span_with_gamma.rank() > span_without_gamma.rank(),
      "gamma cannot distinguish the two B11 multiplicity copies via a 0-cycle")

# The only A5-invariant operators on the 2-transitive six-context carrier are
# a I+b Adj(K6).  Neither supplies a nonzero odd D for the structural gamma.
I6 = sy.eye(6)
K6 = sy.ones(6) - I6
a, b = sy.symbols("a b")
odd_equations = list(a*(I6*gamma + gamma*I6)
                     + b*(K6*gamma + gamma*K6))
solution = sy.solve(odd_equations, (a, b), dict=True)
check("no nonzero A5-invariant context Dirac is odd for gamma",
      solution == [{a: 0, b: 0}],
      "natural K6 adjacency contains triangles and is not bipartite")

# A2 design-filter audit in the actual 24-dimensional carrier.  Blocks are
# plus: 2,4,6 and minus: 2',4,6.  Galois J swaps plus/minus, exchanging the
# endpoints and fixing the shared types.  The equivariant D is identity on
# the shared Schur channels and zero on endpoints.
block_sizes = (2, 4, 6, 2, 4, 6)
starts = [0]
for size in block_sizes:
    starts.append(starts[-1] + size)
gamma24 = sy.diag(*([1]*12 + [-1]*12))
J24 = sy.zeros(24)
for plus_block, minus_block in ((0, 3), (1, 4), (2, 5)):
    size = block_sizes[plus_block]
    for k in range(size):
        i = starts[plus_block] + k
        j = starts[minus_block] + k
        J24[i, j] = J24[j, i] = 1
D24 = sy.zeros(24)
for plus_block, minus_block in ((1, 4), (2, 5)):
    size = block_sizes[plus_block]
    for k in range(size):
        i = starts[plus_block] + k
        j = starts[minus_block] + k
        D24[i, j] = D24[j, i] = 1


def type_projector(blocks):
    out = sy.zeros(24)
    for block in blocks:
        for i in range(starts[block], starts[block+1]):
            out[i, i] = 1
    return out


node_projectors = (
    type_projector((0,)), type_projector((3,)),
    type_projector((1, 4)), type_projector((2, 5)),
)
node_opposites = tuple(J24*p*J24.T for p in node_projectors)
cap_m1_m3 = sy.Matrix([
    [sy.trace(gamma24*a*b) for b in node_opposites]
    for a in node_projectors
])
check("M1 plus M3 Galois double fails PD maximally on common types",
      cap_m1_m3 == sy.Matrix([[0, 2, 0, 0], [-2, 0, 0, 0],
                              [0, 0, 0, 0], [0, 0, 0, 0]])
      and cap_m1_m3.rank() == 2 and cap_m1_m3.det() == 0)
check("M1 plus M3 equivariant D sees only shared 4 plus 6",
      D24.rank() == 20 and len(D24.nullspace()) == 4,
      "2 and 2-prime endpoints are invisible and give a 4D kernel")
check("canonical node algebra has zero one-forms on M1 plus M3",
      all(D24*p == p*D24 for p in node_projectors),
      "D preserves each shared irrep type; endpoint projectors see D=0")

print("-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
print("VERDICT: non-J-invariant noncommutative order-zero algebras exist,")
print("         but no derived six-context candidate passes the full gates.")
if passed != tests:
    raise SystemExit(1)
