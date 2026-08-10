#!/usr/bin/env python3
"""Exact symmetry audit of the order-ten-subgroup Hopf fibrations.

The protocol and its falsifiers were frozen in commit 35b39c7 before this
implementation was run.  All load-bearing arithmetic is in Q(sqrt(5)); no
floating tolerance enters the selector or isotropy conclusions.
"""

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from itertools import permutations, product
import json
from pathlib import Path

import sympy as sp


OUTPUT = Path(__file__).with_name("hopf_symmetry_selector.json")
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    condition = bool(condition)
    passed += int(condition)
    print(f"  [{'PASS' if condition else 'FAIL'}] {label}")
    if detail:
        print(f"         {detail}")


@dataclass(frozen=True)
class Q5:
    """The exact number a+b*sqrt(5), with rational a and b."""

    a: Fraction = Fraction(0)
    b: Fraction = Fraction(0)

    def __init__(self, a=0, b=0):
        object.__setattr__(self, "a", Fraction(a))
        object.__setattr__(self, "b", Fraction(b))

    @staticmethod
    def coerce(value):
        return value if isinstance(value, Q5) else Q5(value)

    def __add__(self, other):
        other = self.coerce(other)
        return Q5(self.a+other.a, self.b+other.b)

    __radd__ = __add__

    def __neg__(self):
        return Q5(-self.a, -self.b)

    def __sub__(self, other):
        return self + (-self.coerce(other))

    def __rsub__(self, other):
        return self.coerce(other) - self

    def __mul__(self, other):
        other = self.coerce(other)
        return Q5(self.a*other.a+5*self.b*other.b,
                  self.a*other.b+self.b*other.a)

    __rmul__ = __mul__

    def inverse(self):
        norm = self.a*self.a-5*self.b*self.b
        if norm == 0:
            raise ZeroDivisionError("zero in Q(sqrt(5))")
        return Q5(self.a/norm, -self.b/norm)

    def __truediv__(self, other):
        return self*self.coerce(other).inverse()

    def __rtruediv__(self, other):
        return self.coerce(other)/self

    def __str__(self):
        if self.b == 0:
            return str(self.a)
        if self.a == 0:
            return f"{self.b}*sqrt(5)"
        sign = "+" if self.b > 0 else "-"
        return f"{self.a}{sign}{abs(self.b)}*sqrt(5)"


ZERO = Q5(0)
ONE = Q5(1)
HALF = Q5(Fraction(1, 2))
PHI_OVER_TWO = Q5(Fraction(1, 4), Fraction(1, 4))
INV_TWO_PHI = Q5(Fraction(-1, 4), Fraction(1, 4))
Quaternion = tuple[Q5, Q5, Q5, Q5]
Matrix = tuple[tuple[Q5, ...], ...]


def qkey(q):
    return tuple((coordinate.a, coordinate.b) for coordinate in q)


def qmul(left: Quaternion, right: Quaternion) -> Quaternion:
    w1, x1, y1, z1 = left
    w2, x2, y2, z2 = right
    return (
        w1*w2-x1*x2-y1*y2-z1*z2,
        w1*x2+x1*w2+y1*z2-z1*y2,
        w1*y2-x1*z2+y1*w2+z1*x2,
        w1*z2+x1*y2-y1*x2+z1*w2,
    )


def qconj(q: Quaternion) -> Quaternion:
    return (q[0], -q[1], -q[2], -q[3])


def qnorm2(q: Quaternion) -> Q5:
    return sum((coordinate*coordinate for coordinate in q), ZERO)


IDENTITY: Quaternion = (ONE, ZERO, ZERO, ZERO)


def build_exact_2i():
    vertices = set()

    for coordinate in range(4):
        for sign in (1, -1):
            value = [ZERO]*4
            value[coordinate] = Q5(sign)
            vertices.add(tuple(value))

    for signs in product((1, -1), repeat=4):
        vertices.add(tuple(Q5(sign)*HALF for sign in signs))

    base = (ZERO, HALF, PHI_OVER_TWO, INV_TWO_PHI)
    even_permutations = []
    for permutation in permutations(range(4)):
        inversions = sum(permutation[i] > permutation[j]
                         for i in range(4) for j in range(i+1, 4))
        if inversions % 2 == 0:
            even_permutations.append(permutation)
    for permutation in even_permutations:
        coordinates = [base[permutation[index]] for index in range(4)]
        nonzero = [index for index, value in enumerate(coordinates)
                   if value != ZERO]
        for signs in product((1, -1), repeat=len(nonzero)):
            value = list(coordinates)
            for index, sign in zip(nonzero, signs):
                value[index] = Q5(sign)*value[index]
            vertices.add(tuple(value))

    return tuple(sorted(vertices, key=qkey))


def element_order(element, limit=120):
    power = IDENTITY
    for order in range(1, limit+1):
        power = qmul(power, element)
        if power == IDENTITY:
            return order
    return None


def cyclic_subgroup(generator, order):
    result = []
    power = IDENTITY
    for _ in range(order):
        result.append(power)
        power = qmul(power, generator)
    return frozenset(result)


def subgroup_key(subgroup):
    return tuple(qkey(element) for element in sorted(subgroup, key=qkey))


def zero_matrix(dimension):
    return tuple(tuple(ZERO for _ in range(dimension))
                 for _ in range(dimension))


def identity_matrix(dimension):
    return tuple(tuple(ONE if row == col else ZERO
                       for col in range(dimension))
                 for row in range(dimension))


def madd(left: Matrix, right: Matrix) -> Matrix:
    return tuple(tuple(left[row][col]+right[row][col]
                       for col in range(len(left[0])))
                 for row in range(len(left)))


def msub(left: Matrix, right: Matrix) -> Matrix:
    return tuple(tuple(left[row][col]-right[row][col]
                       for col in range(len(left[0])))
                 for row in range(len(left)))


def mscale(scalar: Q5, matrix: Matrix) -> Matrix:
    return tuple(tuple(scalar*value for value in row) for row in matrix)


def mmul(left: Matrix, right: Matrix) -> Matrix:
    return tuple(tuple(sum((left[row][inner]*right[inner][col]
                            for inner in range(len(right))), ZERO)
                       for col in range(len(right[0])))
                 for row in range(len(left)))


def transpose(matrix: Matrix) -> Matrix:
    return tuple(tuple(matrix[row][col] for row in range(len(matrix)))
                 for col in range(len(matrix[0])))


def trace(matrix: Matrix) -> Q5:
    return sum((matrix[index][index] for index in range(len(matrix))), ZERO)


def outer(left, right) -> Matrix:
    return tuple(tuple(a*b for b in right) for a in left)


def line_projector(generator: Quaternion) -> Matrix:
    vector = generator[1:]
    norm_square = sum((value*value for value in vector), ZERO)
    return mscale(ONE/norm_square, outer(vector, vector))


def canonical_partition(blocks, index):
    return tuple(sorted(tuple(sorted(index[element] for element in block))
                        for block in blocks))


def coset_partition(group, subgroup, index, handedness):
    remaining = set(group)
    blocks = []
    while remaining:
        representative = min(remaining, key=qkey)
        if handedness == "qH":
            block = {qmul(representative, h) for h in subgroup}
        elif handedness == "Hq":
            block = {qmul(h, representative) for h in subgroup}
        else:
            raise ValueError(handedness)
        blocks.append(block)
        remaining -= block
    return canonical_partition(blocks, index)


def permute_partition(partition, vertex_permutation):
    return tuple(sorted(tuple(sorted(vertex_permutation[index]
                                    for index in block))
                        for block in partition))


def connected_components(number, permutations_):
    adjacency = [set() for _ in range(number)]
    for permutation in permutations_:
        for source, target in enumerate(permutation):
            adjacency[source].add(target)
            adjacency[target].add(source)
    components = []
    unseen = set(range(number))
    while unseen:
        seed = min(unseen)
        stack = [seed]
        component = set()
        while stack:
            vertex = stack.pop()
            if vertex in component:
                continue
            component.add(vertex)
            stack.extend(adjacency[vertex]-component)
        unseen -= component
        components.append(tuple(sorted(component)))
    return tuple(components)


def q5_json(value):
    return {"rational": str(value.a), "sqrt5": str(value.b)}


def matrix_json(matrix):
    return [[q5_json(value) for value in row] for row in matrix]


print("="*78)
print("EXACT HOPF SYMMETRY-SELECTOR AUDIT")
print("="*78)

group = build_exact_2i()
group_set = set(group)
group_index = {element: index for index, element in enumerate(group)}
check("exact coordinate construction has 120 distinct unit quaternions",
      len(group) == 120 and all(qnorm2(element) == ONE for element in group))

closure = all(qmul(left, right) in group_set
              for left in group for right in group)
inverse_identity = all(qmul(element, qconj(element)) == IDENTITY
                       and qmul(qconj(element), element) == IDENTITY
                       for element in group)
check("the exact coordinates close as a group and conjugation gives inverses",
      closure and inverse_identity)

orders = {element: element_order(element) for element in group}
order_distribution = Counter(orders.values())
order_ten_elements = [element for element in group if orders[element] == 10]
subgroups = sorted({cyclic_subgroup(element, 10)
                    for element in order_ten_elements}, key=subgroup_key)
check("all exact element orders resolve and order-ten elements form six C10s",
      None not in order_distribution and len(order_ten_elements) == 24
      and len(subgroups) == 6 and all(len(subgroup) == 10
                                     for subgroup in subgroups),
      f"orders={dict(sorted(order_distribution.items()))}; "
      f"order10={len(order_ten_elements)}, C10={len(subgroups)}")

projectors = []
generators_by_subgroup = []
for subgroup in subgroups:
    generators = sorted((element for element in subgroup
                         if orders[element] == 10), key=qkey)
    generator_projectors = {line_projector(generator)
                            for generator in generators}
    generators_by_subgroup.append(generators)
    projectors.append(next(iter(generator_projectors)))
    if len(generator_projectors) != 1:
        raise RuntimeError("one C10 produced more than one imaginary line")
check("each C10 fixes one generator-independent imaginary line",
      len(set(projectors)) == 6 and all(
          line_projector(generator) == projectors[index]
          for index, generators in enumerate(generators_by_subgroup)
          for generator in generators))

projector_axioms = all(
    projector == transpose(projector)
    and mmul(projector, projector) == projector
    and trace(projector) == ONE
    for projector in projectors
)
check("all six line tensors are exact symmetric rank-one projectors",
      projector_axioms)

qH_partitions = [coset_partition(group, subgroup, group_index, "qH")
                 for subgroup in subgroups]
Hq_partitions = [coset_partition(group, subgroup, group_index, "Hq")
                 for subgroup in subgroups]
partition_shape = all(len(partition) == 12
                      and all(len(block) == 10 for block in partition)
                      for partition in qH_partitions+Hq_partitions)
qH_family = set(qH_partitions)
Hq_family = set(Hq_partitions)
check("both handed coset constructions give six 12-by-10 partitions",
      partition_shape and len(qH_family) == 6 and len(Hq_family) == 6,
      f"qH={len(qH_family)}, Hq={len(Hq_family)}")
check("the two handed partition families are distinct",
      not (qH_family & Hq_family),
      f"overlap={len(qH_family & Hq_family)}, union={len(qH_family | Hq_family)}")

conjugation_permutation = tuple(
    group_index[qconj(element)] for element in group
)
mirrored_qH = {permute_partition(partition, conjugation_permutation)
               for partition in qH_partitions}
check("quaternion conjugation bijects qH fibrations with Hq fibrations",
      mirrored_qH == Hq_family)

subgroup_index = {subgroup: index for index, subgroup in enumerate(subgroups)}
action_permutations = set()
for group_element in group:
    inverse = qconj(group_element)
    permutation = []
    for subgroup in subgroups:
        conjugate_subgroup = frozenset(
            qmul(qmul(inverse, element), group_element)
            for element in subgroup
        )
        permutation.append(subgroup_index[conjugate_subgroup])
    action_permutations.add(tuple(permutation))

components = connected_components(len(subgroups), action_permutations)
constraint_rows = []
for permutation in action_permutations:
    for source, target in enumerate(permutation):
        row = [0]*len(subgroups)
        row[source] = 1
        row[target] -= 1
        constraint_rows.append(row)
constraint_rank = int(sp.Matrix(constraint_rows).rank())
fixed_dimension = len(subgroups)-constraint_rank
check("conjugation acts transitively on the six C10/projector choices",
      components == (tuple(range(6)),),
      f"distinct permutations={len(action_permutations)}, orbits={components}")
check("the invariant coefficient space is only the equal-weight line",
      fixed_dimension == 1 and constraint_rank == 5,
      f"constraint rank={constraint_rank}, fixed dimension={fixed_dimension}")

gram = tuple(tuple(trace(mmul(left, right)) for right in projectors)
             for left in projectors)
off_diagonal_gram = {gram[row][col] for row in range(6)
                     for col in range(6) if row != col}
check("the six axes are exact equiangular lines",
      off_diagonal_gram == {Q5(Fraction(1, 5))},
      "off-diagonal projector overlap="
      + ", ".join(map(str, sorted(off_diagonal_gram,
                                   key=lambda x: (x.a, x.b)))))

frame = zero_matrix(3)
for projector in projectors:
    frame = madd(frame, projector)
expected_frame = mscale(Q5(2), identity_matrix(3))
traceless_frame = msub(
    frame, mscale(trace(frame)/Q5(3), identity_matrix(3))
)
check("the exact six-axis frame operator is 2 I_3",
      frame == expected_frame and traceless_frame == zero_matrix(3),
      f"trace={trace(frame)}; exact traceless residual=0")

# Lift the frame identity independently at every exact 600-cell vertex.  A
# generator need not have unit imaginary part, hence each outer product is
# divided by its own exact squared norm.
right_lift_ok = True
left_lift_ok = True
I4 = identity_matrix(4)
for point in group:
    tangent = msub(I4, outer(point, point))
    right_sum = zero_matrix(4)
    left_sum = zero_matrix(4)
    for generators in generators_by_subgroup:
        generator = generators[0]
        pure_axis = (ZERO,)+generator[1:]
        axis_norm_square = qnorm2(pure_axis)
        right_vector = qmul(point, pure_axis)
        left_vector = qmul(pure_axis, point)
        right_sum = madd(right_sum,
                         mscale(ONE/axis_norm_square,
                                outer(right_vector, right_vector)))
        left_sum = madd(left_sum,
                        mscale(ONE/axis_norm_square,
                               outer(left_vector, left_vector)))
    right_lift_ok &= right_sum == mscale(Q5(2), tangent)
    left_lift_ok &= left_sum == mscale(Q5(2), tangent)
check("both chiral frame identities lift exactly at all 120 vertices",
      right_lift_ok and left_lift_ok)

payload = {
    "protocol_commit": "35b39c7",
    "arithmetic": "exact Q(sqrt(5)); entries are rational + sqrt5 coefficient",
    "group": {
        "size": len(group),
        "order_distribution": {str(key): value
                               for key, value in sorted(order_distribution.items())},
        "order_ten_elements": len(order_ten_elements),
        "cyclic_order_ten_subgroups": len(subgroups),
    },
    "fibration_partitions": {
        "qH": len(qH_family),
        "Hq": len(Hq_family),
        "overlap": len(qH_family & Hq_family),
        "union": len(qH_family | Hq_family),
        "mirror_bijection": mirrored_qH == Hq_family,
    },
    "conjugation_action": {
        "distinct_permutations": len(action_permutations),
        "orbits": [list(component) for component in components],
        "constraint_rank": constraint_rank,
        "fixed_coefficient_dimension": fixed_dimension,
    },
    "projector_gram": matrix_json(gram),
    "frame_operator": matrix_json(frame),
    "frame_trace": q5_json(trace(frame)),
    "traceless_frame": matrix_json(traceless_frame),
    "tangent_lift": {
        "q_times_u": right_lift_ok,
        "u_times_q": left_lift_ok,
        "vertices_checked": len(group),
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
check("exact audit JSON was written", OUTPUT.exists())

print("\n"+"-"*78)
print(f"RESULT: {passed}/{tests} checks passed")
print("DERIVED: there are six C10 choices and twelve handed coset partitions.")
print("DERIVED: each handed six-fibration orbit is an exact tight frame:")
print("         sum P_vertical = 2 P_tangent, average P_vertical = P_tangent/3.")
print("DERIVED NEGATIVE: unbroken full symmetry selects neither one fibration")
print("                  nor a quadratic anisotropy in their linear span.")
print("OPEN: a dynamical or boundary symmetry-breaking selector.")
raise SystemExit(0 if passed == tests else 1)
