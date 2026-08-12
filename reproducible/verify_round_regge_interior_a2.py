#!/usr/bin/env python3
"""Preregistered full-path audit of ordinary de Rham A2.

Protocol commit: 033c935.

The metric is the frozen affine interpolation between the fixed regular
600-cell Regge metric and the radially pulled-back unit-round metric.  The
coefficient includes the bulk curvature, the de Rham transmittal face term,
and the exact complete-exterior conical edge term.  Numerical quadrature is
performed at every preregistered order and every u=j/200.  A finite grid can
refute round selection but cannot prove the continuum inequality.
"""

import json
from math import pi, sqrt
from pathlib import Path

import numpy as np
import sympy as sy


HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "round_regge_interior_a2.json"
PROTOCOL_COMMIT = "033c935"
ORDERS = (16, 24, 32, 40, 48)
GRID = np.linspace(0.0, 1.0, 201)
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


def gauss_unit(order):
    nodes, weights = np.polynomial.legendre.leggauss(order)
    return (nodes + 1.0) / 2.0, weights / 2.0


def tetra_quadrature(order, vertices):
    nodes, weights = gauss_unit(order)
    x, y, z = np.meshgrid(nodes, nodes, nodes, indexing="ij")
    wx, wy, wz = np.meshgrid(weights, weights, weights, indexing="ij")
    x = x.ravel()
    y = y.ravel()
    z = z.ravel()
    weight = (wx * wy * wz).ravel() * (1.0 - x) ** 2 * (1.0 - y)
    lambdas = np.column_stack((
        (1.0 - x) * (1.0 - y) * (1.0 - z),
        x,
        (1.0 - x) * y,
        (1.0 - x) * (1.0 - y) * z,
    ))
    points = lambdas @ vertices
    edge_matrix = np.column_stack((
        vertices[1] - vertices[0],
        vertices[2] - vertices[0],
        vertices[3] - vertices[0],
    ))
    weight *= abs(np.linalg.det(edge_matrix))
    return points, weight


def face_quadrature(order, face_vertices):
    nodes, weights = gauss_unit(order)
    x, y = np.meshgrid(nodes, nodes, indexing="ij")
    wx, wy = np.meshgrid(weights, weights, indexing="ij")
    x = x.ravel()
    y = y.ravel()
    weight = (wx * wy).ravel() * (1.0 - x)
    lambdas = np.column_stack((
        (1.0 - x) * (1.0 - y),
        x,
        (1.0 - x) * y,
    ))
    points = lambdas @ face_vertices
    first = face_vertices[1] - face_vertices[0]
    second = face_vertices[2] - face_vertices[0]
    return points, weight, first, second


def edge_quadrature(order, scale):
    t, weight = np.polynomial.legendre.leggauss(2 * order)
    points = scale * np.column_stack((np.ones_like(t), t, t))
    tangent = scale * np.array((0.0, 1.0, 1.0))
    return points, weight, tangent


def metric_scalars(points, u, a2):
    s2 = np.einsum("ni,ni->n", points, points)
    r2 = a2 + s2
    q = 1.0 - u + u / r2
    b = -u / r2**2
    p = q + b * s2
    return s2, r2, q, p, b


def inverse_inner(left, right, points, q, p, b):
    euclidean = np.einsum("...i,...i->...", left, right)
    y_left = np.einsum("ni,...i->n", points, left)
    y_right = np.einsum("ni,...i->n", points, right)
    return euclidean / q - b * y_left * y_right / (q * p)


def evaluate(order, u, vertices, a2, rho):
    tetra_points, tetra_weights = tetra_quadrature(order, vertices)
    face_points, face_weights, face_e1, face_e2 = face_quadrature(
        order, vertices[1:]
    )
    edge_points, edge_weights, edge_tangent = edge_quadrature(
        order, rho / sqrt(3.0)
    )

    # Bulk volume and scalar curvature.  In radial Cartesian coordinates
    # det(g)^(1/2)=q*sqrt(p).
    _, tetra_r2, tetra_q, tetra_p, _ = metric_scalars(
        tetra_points, u, a2
    )
    tetra_density = tetra_q * np.sqrt(tetra_p)
    one_tetra_volume = np.dot(tetra_weights, tetra_density)
    scalar_curvature = (
        8.0 * u * a2 / (tetra_q * tetra_p * tetra_r2**3)
        - 2.0 * u / (tetra_q**2 * tetra_r2**2)
    )
    volume = 600.0 * one_tetra_volume
    bulk = -(2.0 / 3.0) * 600.0 * np.dot(
        tetra_weights, scalar_curvature * tetra_density
    )

    # Face opposite v0.  The inward Euclidean covector is v0/rho and the
    # neighbouring facet contributes the reflected identical value.  The
    # GKV convention is L_ab=<nabla_ea eb,N_in>.
    _, face_r2, face_q, face_p, face_b = metric_scalars(
        face_points, u, a2
    )
    normal = vertices[0] / rho
    yn = face_points @ normal
    nn = (
        1.0 / face_q
        - face_b * yn**2 / (face_q * face_p)
    )
    ygi_y = np.einsum("ni,ni->n", face_points, face_points) / face_p
    ygi_n = yn / face_p
    projected_y2 = ygi_y - ygi_n**2 / nn
    gamma_connection = (
        2.0 * u * (1.0 - u)
        / (face_r2**3 * face_q * face_p)
    )
    mean_curvature = (
        gamma_connection * yn * projected_y2 / np.sqrt(nn)
    )

    # The intrinsic face area is sqrt(det(E^T g E)).
    e1 = np.broadcast_to(face_e1, face_points.shape)
    e2 = np.broadcast_to(face_e2, face_points.shape)
    e11 = (
        face_q * np.dot(face_e1, face_e1)
        + face_b * (face_points @ face_e1) ** 2
    )
    e22 = (
        face_q * np.dot(face_e2, face_e2)
        + face_b * (face_points @ face_e2) ** 2
    )
    e12 = (
        face_q * np.dot(face_e1, face_e2)
        + face_b * (face_points @ face_e1) * (face_points @ face_e2)
    )
    del e1, e2
    face_density = np.sqrt(e11 * e22 - e12**2)
    face = -(4.0 / 3.0) * 1200.0 * np.dot(
        face_weights, 2.0 * mean_curvature * face_density
    )

    # Edge v0-v1.  The two inward face covectors point toward v2 and v3.
    _, _, edge_q, edge_p, edge_b = metric_scalars(edge_points, u, a2)
    n1 = vertices[2] / rho
    n2 = vertices[3] / rho
    n11 = inverse_inner(n1, n1, edge_points, edge_q, edge_p, edge_b)
    n22 = inverse_inner(n2, n2, edge_points, edge_q, edge_p, edge_b)
    n12 = inverse_inner(n1, n2, edge_points, edge_q, edge_p, edge_b)
    cosine = -n12 / np.sqrt(n11 * n22)
    theta = np.arccos(np.clip(cosine, -1.0, 1.0))
    beta = 5.0 * theta
    tangent = np.broadcast_to(edge_tangent, edge_points.shape)
    tangent2 = (
        edge_q * np.dot(edge_tangent, edge_tangent)
        + edge_b * (edge_points @ edge_tangent) ** 2
    )
    del tangent
    edge_density = np.sqrt(tangent2)
    cone_density = (
        16.0 * pi**2 / (3.0 * beta)
        + 8.0 * beta / 3.0
        - 8.0 * pi
    )
    edge = 720.0 * np.dot(edge_weights, cone_density * edge_density)

    raw = bulk + face + edge
    normalized = (2.0 * pi**2 / volume) ** (1.0 / 3.0) * raw
    return {
        "volume": float(volume),
        "bulk": float(bulk),
        "face": float(face),
        "edge": float(edge),
        "raw": float(raw),
        "normalized": float(normalized),
        "beta_min": float(np.min(beta)),
        "beta_max": float(np.max(beta)),
    }


print("=" * 78)
print("FULL ROUND--REGGE ORDINARY DE RHAM A2 PATH")
print("=" * 78)

# -------------------------------------------------------------------------
# 1. Exact local identities and sign conventions.
# -------------------------------------------------------------------------
u, z, a2_symbol = sy.symbols("u z a2", positive=True)
r2_symbol = a2_symbol + z
q_symbol = 1 - u + u / r2_symbol
b_symbol = -u / r2_symbol**2
p_symbol = sy.factor(q_symbol + b_symbol * z)
gamma_unsimplified = (
    2 * u / r2_symbol**3 - 2 * b_symbol**2 / q_symbol
) / p_symbol
gamma_closed = 2 * u * (1 - u) / (
    r2_symbol**3 * q_symbol * p_symbol
)
check(
    "the connection coefficient has the exact u(1-u) factor",
    sy.simplify(gamma_unsimplified - gamma_closed) == 0,
)

s = sy.symbols("s", positive=True)
r2_s = a2_symbol + s**2
q_s = 1 - u + u / r2_s
p_s = 1 - u + u * a2_symbol / r2_s**2
warping = s * sy.sqrt(q_s)
warping_xi = sy.diff(warping, s) / sy.sqrt(p_s)
warping_xixi = sy.diff(warping_xi, s) / sy.sqrt(p_s)
curvature_warped = sy.simplify(
    -4 * warping_xixi / warping + 2 * (1 - warping_xi**2) / warping**2
)
curvature_closed = (
    8 * u * a2_symbol / (q_s * p_s * r2_s**3)
    - 2 * u / (q_s**2 * r2_s**2)
)
check(
    "the radial warped-product curvature reduces to the implemented formula",
    sy.simplify(curvature_warped - curvature_closed) == 0,
)
check(
    "the bulk curvature is zero at Regge and six at round",
    sy.simplify(curvature_closed.subs(u, 0)) == 0
    and sy.simplify(curvature_closed.subs(u, 1) - 6) == 0,
)

trace_identity = sy.simplify(
    sy.Rational(1, 6) * (2 * sy.Symbol("L") * 8 - 6 * 4 * sy.Symbol("L"))
)
check(
    "the full de Rham transmittal trace is exactly -4L/3",
    trace_identity == -sy.Rational(4, 3) * sy.Symbol("L"),
    "Tr(I)=8 and Tr(U)=4L in dimension three",
)

# -------------------------------------------------------------------------
# 2. Frozen regular tetrahedron and quadrature controls.
# -------------------------------------------------------------------------
a2 = float((7.0 + 3.0 * sqrt(5.0)) / 16.0)
rho = sqrt(1.0 - a2)
scale = rho / sqrt(3.0)
vertices = scale * np.array((
    (1.0, 1.0, 1.0),
    (1.0, -1.0, -1.0),
    (-1.0, 1.0, -1.0),
    (-1.0, -1.0, 1.0),
))
length = (sqrt(5.0) - 1.0) / 2.0
pair_lengths = [
    np.linalg.norm(vertices[left] - vertices[right])
    for left in range(4) for right in range(left + 1, 4)
]
check(
    "the frozen centered tetrahedron has all six edges 1/phi",
    max(abs(value - length) for value in pair_lengths) < 2e-15,
)

euclidean_volume = length**3 / (6.0 * sqrt(2.0))
euclidean_face_area = sqrt(3.0) * length**2 / 4.0
quadrature_controls = []
for order in ORDERS:
    tq_points, tq_weights = tetra_quadrature(order, vertices)
    fq_points, fq_weights, fq_e1, fq_e2 = face_quadrature(order, vertices[1:])
    eq_points, eq_weights, eq_tangent = edge_quadrature(order, scale)
    tetra_constant = float(np.sum(tq_weights))
    face_constant = float(
        np.sum(fq_weights)
        * np.linalg.norm(np.cross(fq_e1, fq_e2))
    )
    edge_constant = float(np.sum(eq_weights) * np.linalg.norm(eq_tangent))
    quadrature_controls.append({
        "order": order,
        "tetra_error": tetra_constant - euclidean_volume,
        "face_error": face_constant - euclidean_face_area,
        "edge_error": edge_constant - length,
    })
check(
    "all frozen Duffy rules integrate the Euclidean cell measures",
    max(
        abs(item[key])
        for item in quadrature_controls
        for key in ("tetra_error", "face_error", "edge_error")
    ) < 3e-15,
)

# -------------------------------------------------------------------------
# 3. Execute every frozen order on the complete preregistered grid.
# -------------------------------------------------------------------------
all_runs = {}
for order in ORDERS:
    print(f"  evaluating order {order} on 201 frozen u values")
    all_runs[order] = [
        evaluate(order, float(value), vertices, a2, rho) for value in GRID
    ]

round_exact = -8.0 * pi**2
regge_volume_exact = 50.0 * sqrt(2.0) / ((1.0 + sqrt(5.0)) / 2.0) ** 3
regge_beta_exact = 5.0 * np.arccos(1.0 / 3.0)
regge_edge_raw_exact = 720.0 * length * (
    16.0 * pi**2 / (3.0 * regge_beta_exact)
    + 8.0 * regge_beta_exact / 3.0
    - 8.0 * pi
)
regge_normalized_exact = (
    (2.0 * pi**2 / regge_volume_exact) ** (1.0 / 3.0)
    * regge_edge_raw_exact
)

endpoint_errors = []
for order in ORDERS:
    left = all_runs[order][0]
    right = all_runs[order][-1]
    endpoint_errors.append({
        "order": order,
        "regge_volume": left["volume"] - regge_volume_exact,
        "regge_bulk": left["bulk"],
        "regge_face": left["face"],
        "regge_A2": left["normalized"] - regge_normalized_exact,
        "round_volume": right["volume"] - 2.0 * pi**2,
        "round_face": right["face"],
        "round_edge": right["edge"],
        "round_A2": right["normalized"] - round_exact,
    })
check(
    "all five quadrature orders recover both exact endpoint decompositions",
    max(
        abs(item[key]) for item in endpoint_errors for key in item if key != "order"
    ) < 2e-11,
)

reference_run = all_runs[ORDERS[-1]]
check(
    "the mandatory transmittal face term is positive at every interior u",
    abs(reference_run[0]["face"]) < 1e-14
    and abs(reference_run[-1]["face"]) < 1e-14
    and min(item["face"] for item in reference_run[1:-1]) > 0.0,
    "it vanishes only at the two smooth/piecewise-flat endpoints",
)
check(
    "the edge audit uses the varying interior cone angle, not an endpoint interpolation",
    abs(reference_run[0]["beta_min"] - regge_beta_exact) < 2e-14
    and abs(reference_run[-1]["beta_max"] - 2.0 * pi) < 2e-14
    and all(
        regge_beta_exact < item["beta_min"] <= item["beta_max"] < 2.0 * pi
        for item in reference_run[1:-1]
    )
    and max(
        item["beta_max"] - item["beta_min"] for item in reference_run[1:-1]
    ) > 1e-3,
)

last_orders = ORDERS[-3:]
fields = ("volume", "bulk", "face", "edge", "raw", "normalized")
convergence_widths = {}
for field in fields:
    widths = []
    for index in range(len(GRID)):
        values = [all_runs[order][index][field] for order in last_orders]
        widths.append(max(values) - min(values))
    convergence_widths[field] = max(widths)
check(
    "the last three frozen orders enclose every reported component tightly",
    max(convergence_widths.values()) < 5e-10,
    ", ".join(f"{key}={value:.2e}" for key, value in convergence_widths.items()),
)

grid_differences = {}
grid_margins = {}
for order in ORDERS:
    normalized = np.array([item["normalized"] for item in all_runs[order]])
    grid_differences[order] = np.diff(normalized)
    grid_margins[order] = normalized[:-1] - normalized[-1]
check(
    "every frozen-order grid is strictly decreasing toward the round endpoint",
    all(np.max(grid_differences[order]) < 0.0 for order in ORDERS),
    "largest adjacent difference="
    f"{max(np.max(grid_differences[order]) for order in ORDERS):.3e}",
)
check(
    "no preregistered interior grid point lies below round",
    all(np.min(grid_margins[order]) > 0.0 for order in ORDERS),
    "smallest interior margin="
    f"{min(np.min(grid_margins[order]) for order in ORDERS):.3e}",
)

# Factoring the numerically observed quadratic approach is useful hostile
# reconnaissance, but neither a vanishing first variation nor positivity on
# the continuous interval is inferred from a finite grid.
final_run = all_runs[ORDERS[-1]]
round_margin = np.array([
    item["normalized"] - final_run[-1]["normalized"] for item in final_run
])
quadratic_quotient = round_margin[:-1] / (1.0 - GRID[:-1]) ** 2
check(
    "the sampled round-endpoint quadratic quotient remains positive",
    np.min(quadratic_quotient) > 0.0,
    f"range=[{np.min(quadratic_quotient):.6g},{np.max(quadratic_quotient):.6g}]",
)

sample_indices = (0, 20, 50, 100, 150, 180, 198, 199, 200)
samples = [
    {"u": float(GRID[index]), **final_run[index]} for index in sample_indices
]
result = {
    "protocol_commit": PROTOCOL_COMMIT,
    "status": "PATTERN ONLY -- frozen grid prefers round; no continuum certificate",
    "operator": "ordinary complete-exterior Hodge-de Rham Laplacian",
    "metric_path": "g_u=(1-u)I+u(I/r^2-yy^T/r^4)",
    "heat_convention": "Tr exp(-t Delta)~(4*pi*t)^(-3/2)[A0+t A2+...]",
    "components": {
        "bulk": "-(2/3)*600*integral R dV",
        "face": "-(4/3)*1200*integral tr(L+ + L-) dA",
        "edge": "720*integral[16*pi^2/(3*beta)+8*beta/3-8*pi] dl",
    },
    "orders": list(ORDERS),
    "grid": "u=j/200, j=0,...,200",
    "quadrature_controls": quadrature_controls,
    "endpoint_errors": endpoint_errors,
    "maximum_last_three_order_widths": convergence_widths,
    "smallest_interior_margin_above_round": float(np.min(round_margin[:-1])),
    "largest_adjacent_grid_difference": float(np.max(grid_differences[ORDERS[-1]])),
    "quadratic_quotient_range": [
        float(np.min(quadratic_quotient)),
        float(np.max(quadratic_quotient)),
    ],
    "samples_order_48": samples,
    "order_48_grid": [
        {"u": float(value), **entry} for value, entry in zip(GRID, final_run)
    ],
    "verdict": (
        "PATTERN: all preregistered converged grid values lie strictly above "
        "the round endpoint and decrease toward it.  The finite grid cannot "
        "exclude an unsampled dip, so global path selection remains OPEN."
    ),
    "scope": (
        "one frozen affine metric path and one ordinary A2 coefficient; no "
        "cutoff, Lorentzian dynamics, absolute scale, Newton constant, or "
        "whole-space H4-invariant metric theorem"
    ),
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(result, indent=2) + "\n")

print("-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
print(result["verdict"])
if passed != tests:
    raise SystemExit(1)
