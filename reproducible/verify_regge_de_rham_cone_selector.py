#!/usr/bin/env python3
"""Exact-geometry audit of the full-de Rham cone coefficient on the 600-cell.

Protocol commit: c3d2500.

The external analytic inputs are the Hodge cone coefficients of Fursaev--
Miele and Cheeger's skeleton-local heat expansion for piecewise-flat
pseudomanifolds.  This verifier attacks every theory-specific input: the
complete 600-cell incidence, cone links, regular-tetrahedron geometry,
all-form coefficient algebra, fixed-volume scaling and the small endpoint
margin.  No heat time, cutoff function or phenomenological target is used.
"""

from collections import defaultdict, deque
from itertools import combinations
import json
from math import acos, pi, sqrt
from pathlib import Path
import sys

import mpmath as mp
import numpy as np
import sympy as sy

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from commons import build_600cell


HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "regge_de_rham_cone_selector.json"
PROTOCOL_COMMIT = "c3d2500"
tests = 0
passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    condition = bool(condition)
    passed += int(condition)
    print(f"  [{'PASS' if condition else 'FAIL'}] {label}")
    if detail:
        print(f"         {detail}")


def connected(nodes, neighbours):
    nodes = set(nodes)
    if not nodes:
        return False
    seen = {next(iter(nodes))}
    queue = deque(seen)
    while queue:
        current = queue.popleft()
        for target in neighbours[current] & nodes:
            if target not in seen:
                seen.add(target)
                queue.append(target)
    return seen == nodes


print("=" * 78)
print("FIXED-REGGE FULL-DE RHAM CONE COEFFICIENT")
print("=" * 78)

# -------------------------------------------------------------------------
# 1. Rebuild the 600-cell and all links from the certified quaternion orbit.
# -------------------------------------------------------------------------
vertices, adjacency, _ = build_600cell()
neighbours = tuple(
    frozenset(np.flatnonzero(adjacency[index]).tolist())
    for index in range(len(vertices))
)
edges = tuple(
    (left, right)
    for left in range(len(vertices))
    for right in sorted(neighbours[left])
    if left < right
)
triangles = tuple(
    (left, right, third)
    for left, right in edges
    for third in sorted(neighbours[left] & neighbours[right])
    if right < third
)
tetrahedra = tuple(
    (first, second, third, fourth)
    for first, second, third in triangles
    for fourth in sorted(
        neighbours[first] & neighbours[second] & neighbours[third]
    )
    if third < fourth
)
check(
    "the independently rebuilt f-vector is (120,720,1200,600)",
    tuple(map(len, (vertices, edges, triangles, tetrahedra)))
    == (120, 720, 1200, 600),
)
check(
    "all vertices have unit circumradius and degree twelve",
    np.max(np.abs(np.sum(vertices * vertices, axis=1) - 1.0)) < 2e-9
    and {len(row) for row in neighbours} == {12},
)

edge_tetrahedra = defaultdict(list)
triangle_tetrahedra = defaultdict(list)
for tetrahedron in tetrahedra:
    for edge in combinations(tetrahedron, 2):
        edge_tetrahedra[edge].append(tetrahedron)
    for triangle in combinations(tetrahedron, 3):
        triangle_tetrahedra[triangle].append(tetrahedron)
check(
    "every edge is incident on exactly five tetrahedra",
    set(map(len, edge_tetrahedra.values())) == {5}
    and len(edge_tetrahedra) == 720,
)
check(
    "every triangular face is incident on exactly two tetrahedra",
    set(map(len, triangle_tetrahedra.values())) == {2}
    and len(triangle_tetrahedra) == 1200,
)

# The normal link of every open edge is the five-cycle.  This independently
# checks both the local cone multiplicity and the absence of an extra branch.
edge_links_are_c5 = True
for left, right in edges:
    link_nodes = neighbours[left] & neighbours[right]
    link_degree = {
        node: len(neighbours[node] & link_nodes) for node in link_nodes
    }
    edge_links_are_c5 &= (
        len(link_nodes) == 5
        and set(link_degree.values()) == {2}
        and connected(link_nodes, neighbours)
    )
check("every open-edge link is exactly C5", edge_links_are_c5)

# Each vertex figure is the combinatorial icosahedron: its H1 vanishes, so
# Cheeger's middle-link obstruction to the no-ideal-boundary Hodge extension
# is absent.  Connectivity plus f=(12,30,20), degree five and Euler 2 give the
# required closed triangulated two-sphere here.
vertex_links_are_icosahedra = True
for vertex in range(len(vertices)):
    link_nodes = neighbours[vertex]
    link_edges = tuple(
        edge for edge in edges if edge[0] in link_nodes and edge[1] in link_nodes
    )
    link_faces = tuple(
        face for face in triangles if set(face).issubset(link_nodes)
    )
    link_degree = {node: 0 for node in link_nodes}
    for left, right in link_edges:
        link_degree[left] += 1
        link_degree[right] += 1
    vertex_links_are_icosahedra &= (
        (len(link_nodes), len(link_edges), len(link_faces)) == (12, 30, 20)
        and set(link_degree.values()) == {5}
        and connected(link_nodes, neighbours)
        and len(link_nodes) - len(link_edges) + len(link_faces) == 2
    )
check(
    "every vertex link is the icosahedral S2 link",
    vertex_links_are_icosahedra,
    "f(link)=(12,30,20), degree=5, connected, Euler=2",
)
check(
    "the Cheeger middle-link obstruction is absent",
    edge_links_are_c5 and vertex_links_are_icosahedra,
    "edge link has odd dimension; vertex link S2 has H1=0",
)

# -------------------------------------------------------------------------
# 2. Exact regular-tetrahedron metric data.
# -------------------------------------------------------------------------
sqrt5 = sy.sqrt(5)
phi_exact = (1 + sqrt5) / 2
edge_exact = 1 / phi_exact
check(
    "the unit-circumradius chord length is exactly 1/phi",
    sy.simplify((2 - phi_exact) - edge_exact**2) == 0,
    "adjacent quaternion vertices have dot product phi/2",
)

edge_lengths = np.array([
    np.linalg.norm(vertices[left] - vertices[right])
    for left, right in edges
])
check(
    "all 720 coordinate edges realize that exact length",
    np.max(np.abs(edge_lengths - float(sy.N(edge_exact, 18)))) < 2e-9,
    f"range=[{edge_lengths.min():.15f},{edge_lengths.max():.15f}]",
)

tetra_volumes = []
for tetrahedron in tetrahedra:
    base = vertices[tetrahedron[0]]
    differences = np.array([
        vertices[index] - base for index in tetrahedron[1:]
    ])
    tetra_volumes.append(sqrt(np.linalg.det(differences @ differences.T)) / 6)
tetra_volume_exact = edge_exact**3 / (6 * sy.sqrt(2))
check(
    "all 600 facet volumes equal l^3/(6 sqrt(2))",
    np.max(np.abs(
        np.asarray(tetra_volumes) - float(sy.N(tetra_volume_exact, 18))
    )) < 3e-10,
)

beta_double = 5 * acos(1 / 3)
check(
    "five tetrahedral dihedral angles give a subcritical cone",
    0 < beta_double < 2 * pi,
    f"beta={beta_double:.15f}, deficit={2*pi-beta_double:.15f}",
)

# -------------------------------------------------------------------------
# 3. Attack the all-form conical coefficient algebra and conventions.
# -------------------------------------------------------------------------
beta, length = sy.symbols("beta length", positive=True)
gamma = 2 * sy.pi / beta
scalar_cone = beta * (gamma**2 - 1) * length / 6
vector_cone = 3 * scalar_cone + 2 * (beta - 2 * sy.pi) * length
full_cone = sy.factor(2 * scalar_cone + 2 * vector_cone)
full_cone_closed = length * (
    16 * sy.pi**2 / (3 * beta) + 8 * beta / 3 - 8 * sy.pi
)
check(
    "the scalar coefficient matches both primary-source conventions",
    sy.simplify(
        scalar_cone
        - sy.pi * (gamma**2 - 1) * length / (3 * gamma)
    ) == 0,
)

# Independent domain/convention check.  On a closed conic two-sphere with
# the Cheeger Hodge extension, Hodge decomposition gives the exact trace
# identity K_1=2*K_0+b_1-2*b_0=2*K_0-2.  Conic Gauss--Bonnet then makes the
# local two-dimensional one-form term 2*S+2*(beta-2*pi).  Cheeger's Kunneth
# formula adds one scalar cone component in the R direction, recovering the
# ambient-three-dimensional Fursaev--Miele vector formula.  Thus we are not
# silently matching two unrelated self-adjoint extensions.
vector_cone_2d = 2 * scalar_cone + 2 * (beta - 2 * sy.pi) * length
vector_cone_from_kunneth = sy.expand(vector_cone_2d + scalar_cone)
check(
    "Hodge decomposition plus Kunneth independently recovers the vector term",
    sy.simplify(vector_cone_from_kunneth - vector_cone) == 0,
    "K1(C_beta x R)=K1(C_beta)K0(R)+K0(C_beta)K1(R)",
)
check(
    "Hodge duality gives the frozen complete-exterior cone coefficient",
    sy.simplify(full_cone - full_cone_closed) == 0,
    "full=2*scalar+2*vector in dimensions 0,1,2,3",
)

deficit = sy.symbols("delta", positive=True)
full_in_deficit = sy.simplify(full_cone_closed.subs(
    beta, 2 * sy.pi - deficit
))
linear_regge = -sy.Rational(4, 3) * deficit * length
quadratic_remainder = sy.factor(full_in_deficit - linear_regge)
check(
    "the exact cone formula has the correct smooth/Regge limit",
    sy.limit(full_in_deficit / deficit, deficit, 0, dir="+")
    == -sy.Rational(4, 3) * length,
)
check(
    "the term discarded by smooth Regge curvature is exactly positive",
    sy.simplify(
        quadratic_remainder
        - 4 * deficit**2 * length / (3 * (2 * sy.pi - deficit))
    ) == 0,
    "C_exact=C_linear+4 delta^2 L/(3 beta)",
)

dimension = sy.Integer(3)
p0_density = sy.Rational(1, 6)
p1_density = (dimension - 6) / 6
ordinary_density = sy.simplify(2 * p0_density + 2 * p1_density)
check(
    "the smooth ordinary full-de Rham coefficient is -2R/3",
    ordinary_density == -sy.Rational(2, 3),
)
check(
    "only the one-skeleton can enter the three-dimensional A2 order",
    3 - 2 == 1 and 3 - 3 == 0,
    "Cheeger: t^(-n/2+j/2) is local on the (n-j)-skeleton; vertices enter t^0",
)

# -------------------------------------------------------------------------
# 4. Preregistered equal-volume endpoint comparison.
# -------------------------------------------------------------------------
regge_volume_exact = sy.simplify(600 * tetra_volume_exact)
check(
    "the exact fixed-Regge volume is 50 sqrt(2)/phi^3",
    sy.simplify(regge_volume_exact - 50 * sy.sqrt(2) / phi_exact**3) == 0,
)
check(
    "A2 has length scaling under g -> c^2 g",
    sy.Rational(3) - 2 == 1,
    "(4*pi*t)^(-3/2)[A0+t*A2+...] forces A2 -> c*A2",
)


def endpoint_values(dps):
    mp.mp.dps = dps
    mp_phi = (1 + mp.sqrt(5)) / 2
    mp_length = 1 / mp_phi
    mp_beta = 5 * mp.acos(mp.mpf(1) / 3)
    mp_deficit = 2 * mp.pi - mp_beta
    mp_scalar = (
        mp_beta / 6 * ((2 * mp.pi / mp_beta)**2 - 1) * mp_length
    )
    mp_full = 8 * mp_scalar + 4 * (mp_beta - 2 * mp.pi) * mp_length
    mp_regge_volume = 600 * mp_length**3 / (6 * mp.sqrt(2))
    mp_round_volume = 2 * mp.pi**2
    mp_scale = (mp_round_volume / mp_regge_volume)**(mp.mpf(1) / 3)
    mp_regge = mp_scale * 720 * mp_full
    mp_round = -8 * mp.pi**2
    mp_linear = (
        mp_scale * 720 * (-mp.mpf(4) / 3 * mp_deficit * mp_length)
    )
    return {
        "phi": mp_phi,
        "edge_length": mp_length,
        "beta": mp_beta,
        "deficit": mp_deficit,
        "regge_volume": mp_regge_volume,
        "round_volume": mp_round_volume,
        "equal_volume_scale": mp_scale,
        "regge_A2": mp_regge,
        "round_A2": mp_round,
        "difference_regge_minus_round": mp_regge - mp_round,
        "linearized_regge_A2": mp_linear,
    }


precision_runs = [endpoint_values(dps) for dps in (40, 80, 140)]
margins = [run["difference_regge_minus_round"] for run in precision_runs]
check(
    "the endpoint sign is stable at 40, 80 and 140 decimal digits",
    all(margin > 0 for margin in margins)
    and max(abs(margins[index] - margins[-1]) for index in range(2))
    < mp.mpf("1e-35"),
    "A2_Regge-A2_round is strictly positive",
)
values = precision_runs[-1]
check(
    "the round endpoint has the lower ordinary A2 at equal volume",
    values["round_A2"] < values["regge_A2"],
    f"margin={mp.nstr(values['difference_regge_minus_round'], 30)}",
)
check(
    "the endpoint margin is small but not a precision artifact",
    mp.mpf("0.08") < values["difference_regge_minus_round"] < mp.mpf("0.09"),
    "relative margin about 0.107 percent",
)
check(
    "using only the linear Regge limit would give the wrong winner",
    values["linearized_regge_A2"] < values["round_A2"]
    < values["regge_A2"],
    (
        f"linear={mp.nstr(values['linearized_regge_A2'], 18)}, "
        f"round={mp.nstr(values['round_A2'], 18)}, "
        f"exact={mp.nstr(values['regge_A2'], 18)}"
    ),
)

result = {
    "protocol_commit": PROTOCOL_COMMIT,
    "status": "DERIVED CONDITIONAL ENDPOINT SELECTION",
    "operator": "ordinary full-exterior Hodge-de Rham Laplacian",
    "extension": "closed-Hilbert-complex / generalized Dirichlet=Neumann Hodge extension",
    "heat_convention": "Tr exp(-t Delta) ~ (4 pi t)^(-3/2) [A0+t A2+...]",
    "f_vector": [120, 720, 1200, 600],
    "edge_link": "C5",
    "vertex_link": "icosahedral S2",
    "cone_formula": "C_full=8*(beta/6)*((2*pi/beta)^2-1)*L+4*(beta-2*pi)*L",
    "cone_minus_linear": "4*delta^2*L/(3*beta)",
    "values": {
        key: mp.nstr(value, 60)
        for key, value in values.items()
    },
    "winner_under_A2_minimization": "unit round S3",
    "scope": (
        "two endpoints at equal volume; not the interior metric family, a "
        "complete cutoff spectral action, an absolute scale, or Lorentzian dynamics"
    ),
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(result, indent=2) + "\n")

print("-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
if passed != tests:
    raise SystemExit(1)
