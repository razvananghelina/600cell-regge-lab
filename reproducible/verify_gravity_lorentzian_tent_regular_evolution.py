#!/usr/bin/env python3
"""Regularity gate for the asymmetric Lorentzian tent pole equation.

Protocol commit: 1035c54.  The sign of the pole derivative was not computed
before registration.  This verifier proves or refutes only local implicit
solvability of the single internal equation; it does not construct the full
Lorentzian pre/post Legendre map.
"""

from fractions import Fraction
import itertools
import json
import math
from pathlib import Path

from flint import arb, ctx
import mpmath as mp
import networkx as nx


HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "gravity_lorentzian_tent_regular_evolution.json"
PROTOCOL_COMMIT = "1035c54"
WITNESS_COMMIT = "cc71574"
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    ok = bool(condition)
    passed += int(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {label}")
    if detail:
        print(f"       {detail}")


class Dual:
    """First-order dual number over Arb balls."""

    def __init__(self, value, derivative=0):
        self.v = value if isinstance(value, arb) else arb(value)
        self.d = derivative if isinstance(derivative, arb) else arb(derivative)

    @staticmethod
    def lift(other):
        return other if isinstance(other, Dual) else Dual(other)

    def __add__(self, other):
        other = self.lift(other)
        return Dual(self.v+other.v, self.d+other.d)

    __radd__ = __add__

    def __neg__(self):
        return Dual(-self.v, -self.d)

    def __sub__(self, other):
        return self + (-self.lift(other))

    def __rsub__(self, other):
        return self.lift(other) - self

    def __mul__(self, other):
        other = self.lift(other)
        return Dual(self.v*other.v, self.d*other.v+self.v*other.d)

    __rmul__ = __mul__

    def reciprocal(self):
        return Dual(1/self.v, -self.d/(self.v*self.v))

    def __truediv__(self, other):
        return self*self.lift(other).reciprocal()

    def __rtruediv__(self, other):
        return self.lift(other)/self

    def sqrt(self):
        root = self.v.sqrt()
        return Dual(root, self.d/(2*root))

    def acos(self):
        return Dual(self.v.acos(), -self.d/(1-self.v*self.v).sqrt())


def theta_dual(rho, q1, q2, q3):
    c1 = (1-rho-q1)/2
    c2 = (1-rho-q2)/2
    c3 = (1-rho-q3)/2
    denominator = 4*(c1*c1+rho)
    p11 = (4*c1*c1-4*c1*c2+4*c2*c2+3*rho)/denominator
    p22 = (4*c1*c1-4*c1*c3+4*c3*c3+3*rho)/denominator
    p12 = (
        2*c1*c1-2*c1*c2-2*c1*c3+4*c2*c3+rho
    )/denominator
    return (p12/(p11*p22).sqrt()).acos()


def weight_dual(rho, q):
    radicand = 4*rho+(1-q-rho)*(1-q-rho)
    return (1+q+rho)/(4*radicand.sqrt())


def global_E_dual(graph, triangles, rho, q_values):
    deficits = {node: Dual(2*arb.pi()) for node in graph}
    for triangle in triangles:
        for node in triangle:
            other = [vertex for vertex in triangle if vertex != node]
            deficits[node] -= theta_dual(
                rho, q_values[node], q_values[other[0]], q_values[other[1]]
            )
    total = Dual(0)
    for node in graph:
        total += deficits[node]*weight_dual(rho, q_values[node])
    return total


def arb_fraction(value):
    return arb(value.numerator)/value.denominator


def arb_ball_from_fractions(lower, upper):
    midpoint = (lower+upper)/2
    radius = (upper-lower)/2
    return arb(arb_fraction(midpoint), arb_fraction(radius))


def final_tetra_gram(values):
    return [
        [
            values[i] if i == j else (values[i]+values[j]-1)/2
            for j in range(3)
        ]
        for i in range(3)
    ]


def principal_minors(matrix):
    result = []
    for size in (1, 2, 3):
        for subset in itertools.combinations(range(3), size):
            if size == 1:
                result.append(matrix[subset[0]][subset[0]])
            elif size == 2:
                i, j = subset
                result.append(
                    matrix[i][i]*matrix[j][j]-matrix[i][j]*matrix[j][i]
                )
            else:
                a, b, c = matrix[0]
                _, d, e = matrix[1]
                _, _, f = matrix[2]
                result.append(a*d*f+2*b*c*e-a*e*e-d*c*c-f*b*b)
    return result


print("="*78)
print("LORENTZIAN TENT POLE REGULARITY / IMPLICIT-EVOLUTION GATE")
print("="*78)

ctx.prec = 220
graph = nx.icosahedral_graph()
triangles = sorted(
    tuple(sorted(clique))
    for clique in nx.enumerate_all_cliques(graph)
    if len(clique) == 3
)
u0 = 0
distances = nx.single_source_shortest_path_length(graph, u0)
shells = [
    sorted(node for node in graph if distances[node] == shell)
    for shell in range(4)
]

automorphisms = list(
    nx.algorithms.isomorphism.GraphMatcher(graph, graph).isomorphisms_iter()
)
stabilizer = [mapping for mapping in automorphisms if mapping[u0] == u0]
check(
    "the inherited local carrier and four stabilizer shells are rebuilt",
    graph.number_of_nodes() == 12
    and graph.number_of_edges() == 30
    and len(triangles) == 20
    and len(automorphisms) == 120
    and len(stabilizer) == 10
    and [len(shell) for shell in shells] == [1, 5, 5, 1],
)

rho0 = arb(1)/4
q_fixed = [None, arb(3)/2, arb(4)/5, arb(3)/2]
x_lo_fraction = Fraction(11, 25)
x_hi_fraction = Fraction(9, 20)
x_lo = arb_fraction(x_lo_fraction)
x_hi = arb_fraction(x_hi_fraction)
x_interval = arb_ball_from_fractions(x_lo_fraction, x_hi_fraction)


def shell_q_dual(x, derivative_shell=None):
    values = [x, Dual(q_fixed[1]), Dual(q_fixed[2]), Dual(q_fixed[3])]
    if derivative_shell is not None:
        values[derivative_shell] = Dual(values[derivative_shell].v, 1)
    return {node: values[distances[node]] for node in graph}


lo_equation = global_E_dual(
    graph, triangles, Dual(rho0), shell_q_dual(Dual(x_lo))
)
hi_equation = global_E_dual(
    graph, triangles, Dual(rho0), shell_q_dual(Dual(x_hi))
)
x_derivative_box = global_E_dual(
    graph, triangles, Dual(rho0), shell_q_dual(Dual(x_interval, 1))
)
check(
    "the prior unique root bracket is independently reproduced",
    lo_equation.v < 0 and hi_equation.v > 0 and x_derivative_box.d > 0,
    f"E(lo)={lo_equation.v}; E(hi)={hi_equation.v}; E_x={x_derivative_box.d}",
)

# The preregistered load-bearing calculation: differentiate with respect to
# rho over the complete x bracket before looking at any sign.
rho_derivative_whole = global_E_dual(
    graph, triangles, Dual(rho0, 1), shell_q_dual(Dual(x_interval))
).d


def strict_sign(value):
    if value > 0:
        return 1
    if value < 0:
        return -1
    return 0


whole_sign = strict_sign(rho_derivative_whole)
subdivision_derivatives = []
if whole_sign == 0:
    width = x_hi_fraction-x_lo_fraction
    for index in range(16):
        left = x_lo_fraction+width*index/16
        right = x_lo_fraction+width*(index+1)/16
        x_box = arb_ball_from_fractions(left, right)
        derivative = global_E_dual(
            graph, triangles, Dual(rho0, 1), shell_q_dual(Dual(x_box))
        ).d
        subdivision_derivatives.append(derivative)
    subdivision_signs = [strict_sign(value) for value in subdivision_derivatives]
    certified_sign = subdivision_signs[0] if (
        subdivision_signs[0] != 0
        and len(set(subdivision_signs)) == 1
    ) else 0
    method = "16 preregistered rational subintervals"
else:
    certified_sign = whole_sign
    subdivision_signs = []
    method = "one complete root-bracket interval"

check(
    "Arb separates the pole derivative from zero on the root bracket",
    certified_sign != 0,
    f"method={method}; E_rho={rho_derivative_whole}",
)

# Rigorous dyadic isolation of x*.  Every sign decision is made by Arb, not by
# the stored decimal from the earlier result.
root_lower = x_lo_fraction
root_upper = x_hi_fraction
for _ in range(180):
    midpoint = (root_lower+root_upper)/2
    midpoint_value = global_E_dual(
        graph,
        triangles,
        Dual(rho0),
        shell_q_dual(Dual(arb_fraction(midpoint))),
    ).v
    if midpoint_value < 0:
        root_lower = midpoint
    elif midpoint_value > 0:
        root_upper = midpoint
    else:
        raise RuntimeError("Arb could not decide a bisection sign")

root_ball = arb_ball_from_fractions(root_lower, root_upper)
q_root = shell_q_dual(Dual(root_ball))
root_equation = global_E_dual(graph, triangles, Dual(rho0), q_root)
root_rho_derivative = global_E_dual(
    graph, triangles, Dual(rho0, 1), q_root
).d
root_x_derivative = global_E_dual(
    graph, triangles, Dual(rho0), shell_q_dual(Dual(root_ball, 1))
).d
check(
    "the isolated root ball contains a zero and has nonzero E_rho and E_x",
    root_equation.v.contains(0)
    and strict_sign(root_rho_derivative) == certified_sign
    and root_x_derivative > 0,
    f"x*={root_ball}; E={root_equation.v}; E_rho={root_rho_derivative}",
)

# Twelve individual boundary derivatives.
individual_derivatives = {}
base_shell_values = [root_ball, q_fixed[1], q_fixed[2], q_fixed[3]]
for varied_node in graph:
    q_values = {
        node: Dual(
            base_shell_values[distances[node]],
            1 if node == varied_node else 0,
        )
        for node in graph
    }
    individual_derivatives[varied_node] = global_E_dual(
        graph, triangles, Dual(rho0), q_values
    ).d

shell_derivatives = []
for varied_shell in range(4):
    q_values = {
        node: Dual(
            base_shell_values[distances[node]],
            1 if distances[node] == varied_shell else 0,
        )
        for node in graph
    }
    shell_derivatives.append(global_E_dual(
        graph, triangles, Dual(rho0), q_values
    ).d)

orbit_constant = all(
    (individual_derivatives[node]-individual_derivatives[shell[0]]).contains(0)
    for shell in shells for node in shell
)
shell_sum_consistent = all(
    (
        shell_derivatives[index]
        - sum((individual_derivatives[node] for node in shell), arb(0))
    ).contains(0)
    for index, shell in enumerate(shells)
)
check(
    "all twelve boundary derivatives respect the four stabilizer orbits",
    orbit_constant,
)
check(
    "collective shell derivatives equal sums of individual edge responses",
    shell_sum_consistent,
    "E_q(shell)="+str([str(value) for value in shell_derivatives]),
)

individual_responses = {
    node: -individual_derivatives[node]/root_rho_derivative for node in graph
}
shell_responses = [
    -value/root_rho_derivative for value in shell_derivatives
]
check(
    "the implicit pole responses are finite on every boundary orbit",
    all(response.is_finite() for response in individual_responses.values())
    and all(response.is_finite() for response in shell_responses),
    "d rho/d q(shell)="+str([str(value) for value in shell_responses]),
)

# Causal margin at the root, reconstructed rather than inherited only by
# reference.  Every actual final tetrahedron is checked.
root_q_by_node = {
    node: base_shell_values[distances[node]] for node in graph
}
root_principal_minors = []
for triangle in triangles:
    root_principal_minors.extend(principal_minors(final_tetra_gram(
        [root_q_by_node[node] for node in triangle]
    )))
check(
    "the regular root remains strictly inside the spacelike boundary domain",
    all(value > 0 for value in root_principal_minors),
    f"smallest lower bound="
    f"{min(float(value.lower()) for value in root_principal_minors):.9g}",
)

# Independent high-precision finite differences of the complete E expression.
mp.mp.dps = 100
root_mp = (
    mp.mpf(root_lower.numerator)/root_lower.denominator
    + mp.mpf(root_upper.numerator)/root_upper.denominator
)/2
rho_mp = mp.mpf(1)/4
q_shell_mp = [root_mp, mp.mpf(3)/2, mp.mpf(4)/5, mp.mpf(3)/2]
q_node_mp = {node: q_shell_mp[distances[node]] for node in graph}


def theta_mp(rho_value, q1, q2, q3):
    c1 = (1-rho_value-q1)/2
    c2 = (1-rho_value-q2)/2
    c3 = (1-rho_value-q3)/2
    denominator = 4*(c1*c1+rho_value)
    p11 = (4*c1*c1-4*c1*c2+4*c2*c2+3*rho_value)/denominator
    p22 = (4*c1*c1-4*c1*c3+4*c3*c3+3*rho_value)/denominator
    p12 = (
        2*c1*c1-2*c1*c2-2*c1*c3+4*c2*c3+rho_value
    )/denominator
    return mp.acos(p12/mp.sqrt(p11*p22))


def global_E_mp(rho_value, q_values):
    deficits = {node: 2*mp.pi for node in graph}
    for triangle in triangles:
        for node in triangle:
            other = [vertex for vertex in triangle if vertex != node]
            deficits[node] -= theta_mp(
                rho_value,
                q_values[node],
                q_values[other[0]],
                q_values[other[1]],
            )
    return mp.fsum(
        deficits[node]
        *(1+q_values[node]+rho_value)
        /(4*mp.sqrt(4*rho_value+(1-q_values[node]-rho_value)**2))
        for node in graph
    )


h = mp.mpf("1e-25")
fd_rho = (
    global_E_mp(rho_mp+h, q_node_mp)-global_E_mp(rho_mp-h, q_node_mp)
)/(2*h)
fd_individual = {}
for varied_node in graph:
    plus = dict(q_node_mp)
    minus = dict(q_node_mp)
    plus[varied_node] += h
    minus[varied_node] -= h
    fd_individual[varied_node] = (
        global_E_mp(rho_mp, plus)-global_E_mp(rho_mp, minus)
    )/(2*h)

fd_shell = []
for shell in shells:
    plus = dict(q_node_mp)
    minus = dict(q_node_mp)
    for node in shell:
        plus[node] += h
        minus[node] -= h
    fd_shell.append((
        global_E_mp(rho_mp, plus)-global_E_mp(rho_mp, minus)
    )/(2*h))


def relative_error(reference, candidate):
    return abs(reference-candidate)/max(mp.mpf(1), abs(reference))


arb_rho_mid = mp.mpf(str(float(root_rho_derivative)))
arb_individual_mid = {
    node: mp.mpf(str(float(individual_derivatives[node]))) for node in graph
}
arb_shell_mid = [mp.mpf(str(float(value))) for value in shell_derivatives]
max_fd_error = max(
    [relative_error(fd_rho, arb_rho_mid)]
    + [relative_error(fd_individual[node], arb_individual_mid[node]) for node in graph]
    + [relative_error(fd_shell[index], arb_shell_mid[index]) for index in range(4)]
)
check(
    "independent centered differences reproduce the pole and boundary derivatives",
    max_fd_error < mp.mpf("2e-6"),
    f"maximum declared relative error={mp.nstr(max_fd_error, 8)}",
)

# If dS/dtau=2*tau*E(tau^2,q), then at E=0 and tau=1/2 the
# second tau derivative equals E_rho exactly.  This is a pole-only Hessian,
# not the full Lagrangian two-form.
tau0 = arb(1)/2
pole_hessian = 4*tau0*tau0*root_rho_derivative
check(
    "the pole-only action Hessian is nondegenerate at the stationary move",
    (pole_hessian-root_rho_derivative).contains(0)
    and strict_sign(pole_hessian) == certified_sign,
    f"d2S/dtau2 (up to the common action factor)={pole_hessian}",
)

root_decimal = mp.nstr(root_mp, 80)
result = {
    "protocol_commit": PROTOCOL_COMMIT,
    "witness_commit": WITNESS_COMMIT,
    "scope": {
        "equation": "single internal tent-pole Regge equation",
        "bulk_variable": "rho=tau^2/a^2",
        "boundary_variables": "twelve q_u, held as canonical boundary data",
        "not_computed": [
            "Lorentzian boundary/corner momenta",
            "full pre/post Legendre map",
            "mixed Lagrangian Hessian",
            "adjacent-move constraint matching",
        ],
    },
    "root": {
        "x": root_decimal,
        "dyadic_enclosure_width": str(float(root_upper-root_lower)),
        "rho": "1/4",
        "E": str(root_equation.v),
    },
    "regularity": {
        "certificate_method": method,
        "E_rho_whole_bracket": str(rho_derivative_whole),
        "certified_sign": "positive" if certified_sign > 0 else "negative",
        "E_rho_at_root": str(root_rho_derivative),
        "E_x_at_root": str(root_x_derivative),
        "pole_hessian_sign": "positive" if certified_sign > 0 else "negative",
        "finite_difference_max_relative_error": mp.nstr(max_fd_error, 12),
    },
    "boundary_response": {
        "individual_E_q_by_node": {
            str(node): str(individual_derivatives[node]) for node in sorted(graph)
        },
        "shell_E_q": [str(value) for value in shell_derivatives],
        "individual_d_rho_d_q_by_node": {
            str(node): str(individual_responses[node]) for node in sorted(graph)
        },
        "collective_shell_d_rho_d_q": [str(value) for value in shell_responses],
    },
    "verdict": {
        "internal_pole": "DERIVED REGULAR LOCAL POLE",
        "implicit_local_evolution_relation": "DERIVED",
        "unique_boundary_state": "NOT REQUIRED AND NOT SELECTED",
        "full_canonical_tent_move": "OPEN",
        "physical_clock": "OPEN",
    },
    "tests": {"passed": passed, "total": tests},
}
OUTPUT.write_text(json.dumps(result, indent=2)+"\n")

print("-"*78)
print(f"Result: {passed}/{tests} checks passed")
print(f"Wrote {OUTPUT}")
if passed != tests:
    raise SystemExit(1)
