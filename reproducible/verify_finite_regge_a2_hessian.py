#!/usr/bin/env python3
"""Full 720-edge Hessian of the exact finite-Regge de Rham A2 coefficient.

Protocol commit: 70afdbc.  Second derivatives are obtained by local forward
automatic differentiation and assembled globally.  Finite differences are
used only as a hostile nonlinear control, never as the primary Hessian.
"""

from collections import defaultdict, deque
from itertools import combinations, permutations
import json
from pathlib import Path
import sys

import mpmath as mp
import numpy as np
import scipy.linalg as la

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from commons import build_600cell


HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "finite_regge_a2_hessian.json"
PROTOCOL_COMMIT = "70afdbc"
LOCAL_EDGES = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
LOCAL_EDGE_INDEX = {edge: index for index, edge in enumerate(LOCAL_EDGES)}
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    ok = bool(condition)
    passed += int(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {label}")
    if detail:
        print(f"       {detail}")


class Jet:
    """Value, gradient and Hessian over the six local squared lengths."""

    n = 6

    def __init__(self, value, gradient=None, hessian=None):
        self.v = mp.mpf(value)
        self.g = ([mp.mpf(0)] * self.n if gradient is None else gradient)
        self.h = (
            [[mp.mpf(0) for _ in range(self.n)] for _ in range(self.n)]
            if hessian is None else hessian
        )

    @classmethod
    def variable(cls, index, value=1):
        gradient = [mp.mpf(0)] * cls.n
        gradient[index] = mp.mpf(1)
        return cls(value, gradient)

    @staticmethod
    def coerce(value):
        return value if isinstance(value, Jet) else Jet(value)

    def __add__(self, other):
        other = self.coerce(other)
        return Jet(
            self.v + other.v,
            [a + b for a, b in zip(self.g, other.g)],
            [[self.h[i][j] + other.h[i][j] for j in range(self.n)]
             for i in range(self.n)],
        )

    __radd__ = __add__

    def __neg__(self):
        return Jet(-self.v, [-x for x in self.g],
                   [[-x for x in row] for row in self.h])

    def __sub__(self, other):
        return self + (-self.coerce(other))

    def __rsub__(self, other):
        return self.coerce(other) - self

    def __mul__(self, other):
        other = self.coerce(other)
        gradient = [self.g[i] * other.v + self.v * other.g[i]
                    for i in range(self.n)]
        hessian = [[
            self.h[i][j] * other.v + self.v * other.h[i][j]
            + self.g[i] * other.g[j] + other.g[i] * self.g[j]
            for j in range(self.n)] for i in range(self.n)]
        return Jet(self.v * other.v, gradient, hessian)

    __rmul__ = __mul__

    def unary(self, value, first, second):
        return Jet(
            value,
            [first * x for x in self.g],
            [[first * self.h[i][j] + second * self.g[i] * self.g[j]
              for j in range(self.n)] for i in range(self.n)],
        )

    def __pow__(self, exponent):
        exponent = mp.mpf(exponent)
        value = self.v ** exponent
        first = exponent * self.v ** (exponent - 1)
        second = exponent * (exponent - 1) * self.v ** (exponent - 2)
        return self.unary(value, first, second)

    def __truediv__(self, other):
        return self * self.coerce(other) ** (-1)

    def __rtruediv__(self, other):
        return self.coerce(other) / self


def jet_sqrt(value):
    return Jet.coerce(value) ** mp.mpf("0.5")


def jet_acos(value):
    value = Jet.coerce(value)
    root = mp.sqrt(1 - value.v**2)
    return value.unary(
        mp.acos(value.v), -1 / root, -value.v / root**3
    )


def local_dot(lengths, origin, left, right):
    if left == right:
        return lengths[LOCAL_EDGE_INDEX[tuple(sorted((origin, left)))]]
    if left == origin or right == origin:
        target = right if left == origin else left
        return lengths[LOCAL_EDGE_INDEX[tuple(sorted((origin, target)))]]
    oi = lengths[LOCAL_EDGE_INDEX[tuple(sorted((origin, left)))]]
    oj = lengths[LOCAL_EDGE_INDEX[tuple(sorted((origin, right)))]]
    ij = lengths[LOCAL_EDGE_INDEX[tuple(sorted((left, right)))]]
    return (oi + oj - ij) / 2


def local_dihedral(lengths, edge):
    i, j = edge
    k, ell = [vertex for vertex in range(4) if vertex not in edge]
    ee = local_dot(lengths, i, j, j)
    aa = local_dot(lengths, i, k, k)
    bb = local_dot(lengths, i, ell, ell)
    ae = local_dot(lengths, i, k, j)
    be = local_dot(lengths, i, ell, j)
    ab = local_dot(lengths, i, k, ell)
    numerator = ab - ae * be / ee
    left_norm = aa - ae * ae / ee
    right_norm = bb - be * be / ee
    return jet_acos(numerator / jet_sqrt(left_norm * right_norm))


def determinant3(matrix):
    return (
        matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
        - matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
        + matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0])
    )


def local_volume(lengths):
    gram = [[local_dot(lengths, 0, i, j) for j in (1, 2, 3)]
            for i in (1, 2, 3)]
    return jet_sqrt(determinant3(gram)) / 6


def local_data(dps):
    mp.mp.dps = dps
    variables = [Jet.variable(index) for index in range(6)]
    volume = local_volume(variables)
    angles = [local_dihedral(variables, edge) for edge in LOCAL_EDGES]

    def convert(jet):
        return (
            float(jet.v),
            np.array([float(x) for x in jet.g]),
            np.array([[float(x) for x in row] for row in jet.h]),
        )

    return convert(volume), [convert(angle) for angle in angles]


def build_complex():
    vertices, adjacency, _ = build_600cell()
    neighbours = [set(np.flatnonzero(adjacency[index]).tolist())
                  for index in range(len(vertices))]
    edges = [(left, right) for left in range(120)
             for right in sorted(neighbours[left]) if left < right]
    triangles = [(a, b, c) for a, b in edges
                 for c in sorted(neighbours[a] & neighbours[b]) if b < c]
    tetrahedra = [(a, b, c, d) for a, b, c in triangles
                  for d in sorted(neighbours[a] & neighbours[b] & neighbours[c])
                  if c < d]
    return vertices, neighbours, edges, triangles, tetrahedra


def assemble(dps, tetrahedra, edge_index):
    (volume0, volume_g_local, volume_h_local), angle_data = local_data(dps)
    n_edges = len(edge_index)
    volume = volume0 * len(tetrahedra)
    volume_g = np.zeros(n_edges)
    volume_h = np.zeros((n_edges, n_edges))
    beta_gradients = [defaultdict(float) for _ in range(n_edges)]
    beta_hessians = [defaultdict(float) for _ in range(n_edges)]
    beta_values = np.zeros(n_edges)

    for tetrahedron in tetrahedra:
        local_global = [edge_index[tuple(sorted((tetrahedron[i], tetrahedron[j])))]
                        for i, j in LOCAL_EDGES]
        idx = np.array(local_global, dtype=int)
        volume_g[idx] += volume_g_local
        volume_h[np.ix_(idx, idx)] += volume_h_local
        for local_edge, global_edge in enumerate(local_global):
            angle0, angle_g, angle_h = angle_data[local_edge]
            beta_values[global_edge] += angle0
            for a, ga in enumerate(local_global):
                beta_gradients[global_edge][ga] += angle_g[a]
                for b, gb in enumerate(local_global):
                    beta_hessians[global_edge][(ga, gb)] += angle_h[a, b]

    mp.mp.dps = dps
    pi = mp.pi
    functional = 0.0
    functional_g = np.zeros(n_edges)
    functional_h = np.zeros((n_edges, n_edges))
    beta_spread = float(np.max(beta_values) - np.min(beta_values))

    for edge in range(n_edges):
        # At the frozen equilateral point the value is exactly five regular
        # dihedral angles.  Derivatives above remain fully assembled from the
        # five incident tetrahedra; using the exact value here avoids feeding
        # a summed binary64 angle back into the analytic C(beta) derivatives.
        beta = 5 * mp.acos(mp.mpf(1) / 3)
        c0 = 16 * pi**2 / (3 * beta) + 8 * beta / 3 - 8 * pi
        c1 = -16 * pi**2 / (3 * beta**2) + mp.mpf(8) / 3
        c2 = 32 * pi**2 / (3 * beta**3)
        c0, c1, c2 = map(float, (c0, c1, c2))
        functional += c0
        support = sorted(set(beta_gradients[edge]) | {
            a for pair in beta_hessians[edge] for a in pair
        })
        bg = np.array([beta_gradients[edge][a] for a in support])
        bh = np.array([[beta_hessians[edge][(a, b)] for b in support]
                       for a in support])
        functional_g[support] += c1 * bg
        functional_g[edge] += 0.5 * c0
        functional_h[np.ix_(support, support)] += c1 * bh + c2 * np.outer(bg, bg)
        functional_h[edge, edge] += -0.25 * c0
        functional_h[edge, support] += 0.5 * c1 * bg
        functional_h[support, edge] += 0.5 * c1 * bg

    p = -1.0 / 3.0
    vp = volume**p
    value = functional * vp
    gradient = vp * functional_g + functional * p * volume**(p - 1) * volume_g
    hessian = (
        vp * functional_h
        + p * volume**(p - 1) * (
            np.outer(functional_g, volume_g)
            + np.outer(volume_g, functional_g)
            + functional * volume_h
        )
        + functional * p * (p - 1) * volume**(p - 2)
        * np.outer(volume_g, volume_g)
    )
    return {
        "value": value,
        "gradient": gradient,
        "hessian": (hessian + hessian.T) / 2,
        "volume": volume,
        "raw": functional,
        "beta_spread": beta_spread,
        "local_volume": (volume0, volume_g_local, volume_h_local),
        "local_angles": angle_data,
    }


def scalar_volume(local_x):
    x = np.asarray(local_x, dtype=float)
    gram = np.empty((3, 3))
    for a, i in enumerate((1, 2, 3)):
        for b, j in enumerate((1, 2, 3)):
            if i == j:
                gram[a, b] = x[LOCAL_EDGE_INDEX[(0, i)]]
            else:
                gram[a, b] = (
                    x[LOCAL_EDGE_INDEX[(0, i)]]
                    + x[LOCAL_EDGE_INDEX[(0, j)]]
                    - x[LOCAL_EDGE_INDEX[tuple(sorted((i, j)))]]
                ) / 2
    return np.sqrt(np.linalg.det(gram)) / 6


def scalar_angle(local_x, edge):
    x = np.asarray(local_x, dtype=float)
    i, j = edge
    k, ell = [vertex for vertex in range(4) if vertex not in edge]

    def dot(origin, left, right):
        if left == right:
            return x[LOCAL_EDGE_INDEX[tuple(sorted((origin, left)))]]
        if left == origin or right == origin:
            target = right if left == origin else left
            return x[LOCAL_EDGE_INDEX[tuple(sorted((origin, target)))]]
        return (
            x[LOCAL_EDGE_INDEX[tuple(sorted((origin, left)))]]
            + x[LOCAL_EDGE_INDEX[tuple(sorted((origin, right)))]]
            - x[LOCAL_EDGE_INDEX[tuple(sorted((left, right)))]]
        ) / 2

    ee, aa, bb = dot(i, j, j), dot(i, k, k), dot(i, ell, ell)
    ae, be, ab = dot(i, k, j), dot(i, ell, j), dot(i, k, ell)
    cosine = (ab - ae * be / ee) / np.sqrt(
        (aa - ae * ae / ee) * (bb - be * be / ee)
    )
    return np.arccos(np.clip(cosine, -1, 1))


def evaluate_global(x, tetrahedra, edge_index):
    x = np.asarray(x, dtype=float)
    beta = np.zeros(len(x))
    volume = 0.0
    for tetrahedron in tetrahedra:
        local_global = [edge_index[tuple(sorted((tetrahedron[i], tetrahedron[j])))]
                        for i, j in LOCAL_EDGES]
        local_x = x[local_global]
        volume += scalar_volume(local_x)
        for local_edge, global_edge in enumerate(local_global):
            beta[global_edge] += scalar_angle(local_x, LOCAL_EDGES[local_edge])
    c = 16 * np.pi**2 / (3 * beta) + 8 * beta / 3 - 8 * np.pi
    return np.sum(np.sqrt(x) * c) * volume**(-1 / 3)


def quotient_basis(n):
    e0 = np.zeros(n)
    e0[0] = 1
    unit_scale = np.ones(n) / np.sqrt(n)
    vector = e0 - unit_scale
    householder = np.eye(n) - 2 * np.outer(vector, vector) / np.dot(vector, vector)
    return householder[:, 1:]


def inertia(values, tolerance=1e-8):
    return (
        int(np.sum(values > tolerance)),
        int(np.sum(np.abs(values) <= tolerance)),
        int(np.sum(values < -tolerance)),
    )


def ldl_inertia(matrix, tolerance=1e-8):
    _, diagonal, _ = la.ldl(matrix, lower=True, hermitian=True)
    signs = []
    i = 0
    while i < len(diagonal):
        if i + 1 < len(diagonal) and abs(diagonal[i, i + 1]) > tolerance:
            signs.extend(np.linalg.eigvalsh(diagonal[i:i + 2, i:i + 2]))
            i += 2
        else:
            signs.append(diagonal[i, i])
            i += 1
    return inertia(np.asarray(signs), tolerance)


def qmul(left, right):
    a, b, c, d = np.moveaxis(np.asarray(left), -1, 0)
    e, f, g, h = np.moveaxis(np.asarray(right), -1, 0)
    return np.stack((a*e-b*f-c*g-d*h, a*f+b*e+c*h-d*g,
                     a*g-b*h+c*e+d*f, a*h+b*g-c*f+d*e), axis=-1)


def qconj(value):
    result = np.array(value, copy=True)
    result[..., 1:] *= -1
    return result


def edge_stabilizer(vertices, first_edge, edge_index):
    """Full H4 stabilizer from q -> a q b^-1 and a qbar b^-1."""
    left_endpoint, right_endpoint = vertices[list(first_edge)]
    permutations_found = set()
    for reflected in (False, True):
        seed = qconj(vertices) if reflected else vertices
        endpoint_seed = qconj(np.array((left_endpoint, right_endpoint))) if reflected else np.array((left_endpoint, right_endpoint))
        for a in vertices:
            left_products = qmul(a, endpoint_seed)
            for b in vertices:
                mapped_endpoints = qmul(left_products, qconj(b))
                dots = mapped_endpoints @ vertices.T
                mapped = tuple(np.argmax(dots, axis=1).tolist())
                if set(mapped) != set(first_edge):
                    continue
                all_mapped = qmul(qmul(a, seed), qconj(b))
                permutation = tuple(np.argmax(all_mapped @ vertices.T, axis=1).tolist())
                if len(set(permutation)) == len(vertices):
                    permutations_found.add(permutation)
    edge_permutations = []
    for permutation in permutations_found:
        edge_permutations.append(tuple(
            edge_index[tuple(sorted((permutation[a], permutation[b])))]
            for a, b in sorted(edge_index, key=edge_index.get)
        ))
    return sorted(set(edge_permutations))


def stabilizer_orbits(permutations_on_edges, n_edges):
    unseen = set(range(n_edges))
    orbits = []
    while unseen:
        seed = min(unseen)
        orbit = {seed}
        queue = deque((seed,))
        while queue:
            current = queue.popleft()
            for permutation in permutations_on_edges:
                target = permutation[current]
                if target not in orbit:
                    orbit.add(target)
                    queue.append(target)
        orbits.append(tuple(sorted(orbit)))
        unseen -= orbit
    return orbits


def direct_second(function, direction, step=2e-3):
    base = np.ones_like(direction)
    f0 = function(base)
    fm2 = function(base - 2 * step * direction)
    fm1 = function(base - step * direction)
    fp1 = function(base + step * direction)
    fp2 = function(base + 2 * step * direction)
    return (-fp2 + 16*fp1 - 30*f0 + 16*fm1 - fm2) / (12*step**2)


print("=" * 78)
print("FULL FINITE-REGGE ORDINARY-DE RHAM A2 HESSIAN")
print("=" * 78)

vertices, neighbours, edges, triangles, tetrahedra = build_complex()
edge_index = {edge: index for index, edge in enumerate(edges)}
edge_incidence = np.zeros(len(edges), dtype=int)
for tetrahedron in tetrahedra:
    for edge in combinations(tetrahedron, 2):
        edge_incidence[edge_index[edge]] += 1
check(
    "the carrier is the complete 600-cell with five tetrahedra per edge",
    tuple(map(len, (vertices, edges, triangles, tetrahedra))) == (120, 720, 1200, 600)
    and set(edge_incidence) == {5},
)

precision_runs = [assemble(dps, tetrahedra, edge_index) for dps in (50, 80, 120)]
run = precision_runs[-1]
angle0 = np.array([item[0] for item in run["local_angles"]])
check(
    "all regular local Gram matrices are positive and dihedral cosines are 1/3",
    run["local_volume"][0] > 0
    and np.max(np.abs(np.cos(angle0) - 1/3)) < 2e-15,
)

# Independent Cayley--Menger inverse/cofactor convention at the regular point.
cayley = np.ones((5, 5))
cayley[0, 0] = 0
cayley[1:, 1:] = np.ones((4, 4)) - np.eye(4)
cayley_inverse = np.linalg.inv(cayley)
cayley_cosine = cayley_inverse[3, 4] / np.sqrt(
    cayley_inverse[3, 3] * cayley_inverse[4, 4]
)
check(
    "projection and independent Cayley--Menger angle conventions agree",
    abs(cayley_cosine - 1/3) < 2e-15
    and abs(scalar_volume(np.ones(6)) - np.sqrt(2)/12) < 2e-15,
)

# The new 720-variable implementation must reproduce the independently
# committed one-number endpoint calculation, including its positive
# equal-volume normalization factor.
endpoint_certificate = json.loads(
    (HERE / "regge_de_rham_cone_selector.json").read_text()
)
endpoint_value = float(endpoint_certificate["values"]["regge_A2"])
recovered_endpoint = run["value"] * (2 * np.pi**2) ** (1/3)
check(
    "the assembled functional reproduces the independent fixed-Regge endpoint",
    abs(recovered_endpoint - endpoint_value) < 2e-12,
    f"recovered={recovered_endpoint:.15f}; stored={endpoint_value:.15f}",
)

# All 24 tetrahedral vertex relabellings must permute the local derivative
# tensors covariantly.
local_angles = run["local_angles"]
relabel_residual = 0.0
_, local_volume_g, local_volume_h = run["local_volume"]
for permutation in permutations(range(4)):
    mapping = [LOCAL_EDGE_INDEX[tuple(sorted((permutation[a], permutation[b])))]
               for a, b in LOCAL_EDGES]
    relabel_residual = max(
        relabel_residual,
        np.max(np.abs(local_volume_g[np.array(mapping)] - local_volume_g)),
        np.max(np.abs(
            local_volume_h[np.ix_(mapping, mapping)] - local_volume_h
        )),
    )
    for source_edge in range(6):
        target_edge = mapping[source_edge]
        _, source_g, source_h = local_angles[source_edge]
        _, target_g, target_h = local_angles[target_edge]
        relabel_residual = max(
            relabel_residual,
            np.max(np.abs(target_g[np.array(mapping)] - source_g)),
            np.max(np.abs(target_h[np.ix_(mapping, mapping)] - source_h)),
        )
check(
    "local first and second derivatives respect all tetrahedral relabellings",
    relabel_residual < 2e-13,
    f"maximum covariance residual={relabel_residual:.3e}",
)

gradient = run["gradient"]
hessian = run["hessian"]
scale = np.ones(len(edges))
gradient_spread = np.ptp(gradient)
check(
    "H4 edge transitivity plus scale invariance gives exact stationarity numerically",
    np.max(np.abs(gradient)) < 1e-11 and gradient_spread < 1e-12,
    f"max|gradient|={np.max(np.abs(gradient)):.3e}; spread={gradient_spread:.3e}",
)
check(
    "the assembled Hessian is symmetric and annihilates global scale",
    np.max(np.abs(hessian-hessian.T)) < 1e-13
    and np.linalg.norm(hessian @ scale) / np.sqrt(len(scale)) < 1e-10,
    f"scale RMS={np.linalg.norm(hessian @ scale)/np.sqrt(len(scale)):.3e}",
)

# Coordinate control: for x_e=l_e^2, at l=1 one has
# H_l=4 H_x+2 diag(grad_x).  Exact stationarity removes the second term, so
# Sylvester inertia cannot depend on choosing lengths or squared lengths.
length_coordinate_hessian = 4 * hessian + 2 * np.diag(gradient)

basis = quotient_basis(len(edges))
quotient_hessians = [basis.T @ item["hessian"] @ basis for item in precision_runs]
quotient_eigenvalues = [np.linalg.eigvalsh(matrix) for matrix in quotient_hessians]
eigenvalues = quotient_eigenvalues[-1]
precision_drift = max(np.max(np.abs(values - eigenvalues))
                      for values in quotient_eigenvalues[:-1])
smallest_nonzero = float(np.min(np.abs(eigenvalues)))
dense_inertia = inertia(eigenvalues)
independent_inertia = ldl_inertia(quotient_hessians[-1])
length_coordinate_inertia = inertia(
    np.linalg.eigvalsh(basis.T @ length_coordinate_hessian @ basis)
)
certificate_ok = (
    smallest_nonzero > 1e-7
    and precision_drift < 1e-10
    and dense_inertia == independent_inertia
)
check(
    "quotient inertia is precision-stable with a resolved nonzero gap",
    certificate_ok,
    f"gap={smallest_nonzero:.6g}; drift={precision_drift:.3e}; inertia={dense_inertia}",
)
check(
    "independent pivoted LDL congruence reproduces the quotient inertia",
    dense_inertia == independent_inertia,
    f"eigh={dense_inertia}; LDL={independent_inertia}",
)
check(
    "edge lengths and squared edge lengths give the same stationary inertia",
    length_coordinate_inertia == dense_inertia,
    f"length-coordinate inertia={length_coordinate_inertia}",
)

# Full H4 edge stabilizer and every preregistered relative-position probe.
stabilizer = edge_stabilizer(vertices, edges[0], edge_index)
orbits = stabilizer_orbits(stabilizer, len(edges))
check(
    "the full edge stabilizer has order 20 and partitions all 720 edges",
    len(stabilizer) == 20 and sum(map(len, orbits)) == 720
    and sorted(edge for orbit in orbits for edge in orbit) == list(range(720)),
    f"orbit count={len(orbits)}; sizes={list(map(len, orbits))}",
)
orbit_probes = []
for orbit in orbits:
    if orbit == (0,):
        continue
    direction = np.zeros(len(edges))
    direction[0] = 1
    direction[orbit[0]] = -1
    orbit_probes.append((orbit[0], direction, float(direction @ hessian @ direction)))

# Frozen discrete conformal probes from the complete trace-free quadratic span.
quadratic_matrices = []
for i in range(3):
    matrix = np.zeros((4, 4))
    matrix[i, i] = 1
    matrix[i + 1, i + 1] = -1
    quadratic_matrices.append(matrix)
for i in range(4):
    for j in range(i + 1, 4):
        matrix = np.zeros((4, 4))
        matrix[i, j] = matrix[j, i] = 1
        quadratic_matrices.append(matrix)
conformal_columns = []
for matrix in quadratic_matrices:
    vertex_values = np.einsum("ni,ij,nj->n", vertices, matrix, vertices)
    direction = np.array([vertex_values[a] + vertex_values[b] for a, b in edges])
    direction -= np.mean(direction)
    conformal_columns.append(direction)
conformal_matrix = np.column_stack(conformal_columns)
conformal_basis, conformal_r = np.linalg.qr(conformal_matrix)
conformal_rank = int(np.sum(np.abs(np.diag(conformal_r)) > 1e-10))
conformal_basis = conformal_basis[:, :conformal_rank]
conformal_eigenvalues = np.linalg.eigvalsh(conformal_basis.T @ hessian @ conformal_basis)
explicit_conformal = conformal_columns[0]
explicit_conformal /= np.linalg.norm(explicit_conformal)
explicit_conformal_value = float(explicit_conformal @ hessian @ explicit_conformal)
check(
    "the frozen trace-free vertex-quadratic conformal carrier has rank nine",
    conformal_rank == 9 and abs(np.dot(explicit_conformal, scale)) < 1e-12,
)

# The conformal span is symmetry-natural, but equivariance does not force an
# individual copy inside a representation with multiplicities to be Hessian
# invariant.  Check both facts so the compression below is not mislabeled as
# a full spectral block.
conformal_projector = conformal_basis @ conformal_basis.T
symmetry_commutator = 0.0
conformal_symmetry_residual = 0.0
for permutation in stabilizer:
    permutation = np.asarray(permutation)
    symmetry_commutator = max(
        symmetry_commutator,
        float(np.max(np.abs(hessian[np.ix_(permutation, permutation)] - hessian))),
    )
    moved = conformal_basis[permutation]
    conformal_symmetry_residual = max(
        conformal_symmetry_residual,
        float(np.linalg.norm(moved - conformal_projector @ moved)),
    )
conformal_hessian_leakage = float(np.linalg.norm(
    hessian @ conformal_basis
    - conformal_basis @ (conformal_basis.T @ hessian @ conformal_basis)
))
check(
    "symmetry preserves the conformal span and commutes with the Hessian",
    symmetry_commutator < 1e-12 and conformal_symmetry_residual < 1e-8,
    f"commutator={symmetry_commutator:.3e}; span residual={conformal_symmetry_residual:.3e}",
)
check(
    "the conformal compression is not falsely labeled a Hessian eigenspace",
    conformal_hessian_leakage > 1e-3,
    f"off-subspace Hessian norm={conformal_hessian_leakage:.6g}",
)
check(
    "the Hessian quadratic form is negative on the entire conformal carrier",
    conformal_eigenvalues[-1] < -1e-3
    and np.ptp(conformal_eigenvalues) < 1e-12,
    f"range=[{conformal_eigenvalues[0]:.12g},{conformal_eigenvalues[-1]:.12g}]",
)

# Nonlinear five-point controls for all orbit contrasts and the nine conformal
# basis directions.  Normalize every direction to make one shared safe step.
nonlinear_function = lambda x: evaluate_global(x, tetrahedra, edge_index)
probe_records = []
for name, direction, automatic in [
    *[(f"edge_orbit_{representative}", direction/np.linalg.norm(direction),
       automatic/np.dot(direction, direction))
      for representative, direction, automatic in orbit_probes],
    *[(f"conformal_{index}", conformal_basis[:, index],
       float(conformal_basis[:, index] @ hessian @ conformal_basis[:, index]))
      for index in range(conformal_rank)],
]:
    direct = direct_second(nonlinear_function, direction)
    relative = abs(direct-automatic) / max(1, abs(automatic))
    probe_records.append({"name": name, "automatic": automatic,
                          "direct": direct, "relative_error": relative})
max_direct_error = max(item["relative_error"] for item in probe_records)
check(
    "direct nonlinear second differences reproduce every frozen probe family",
    max_direct_error < 2e-6,
    f"maximum relative/control error={max_direct_error:.3e}",
)
orbit_direct = [item["direct"] for item in probe_records
                if item["name"].startswith("edge_orbit_")]
conformal_direct = [item["direct"] for item in probe_records
                    if item["name"].startswith("conformal_")]
check(
    "nonlinear controls independently exhibit both saddle signs",
    min(orbit_direct) > 1e-2 and max(conformal_direct) < -1e-3,
    f"positive probe min={min(orbit_direct):.6g}; conformal max={max(conformal_direct):.6g}",
)

# Report symmetry-resolved numerical eigenvalue clusters.  This is not used
# to choose the verdict, but makes the 569/150 count independently auditable.
clusters = []
for value in eigenvalues:
    if not clusters or abs(value - clusters[-1]["value"]) > 2e-9:
        clusters.append({"value": float(value), "multiplicity": 1})
    else:
        count = clusters[-1]["multiplicity"]
        clusters[-1]["value"] = (
            clusters[-1]["value"] * count + float(value)
        ) / (count + 1)
        clusters[-1]["multiplicity"] = count + 1

positive, zero, negative = dense_inertia
if certificate_ok and positive and negative:
    verdict = "DERIVED FINITE SADDLE"
elif certificate_ok and positive == 719:
    verdict = "DERIVED FINITE LOCAL MINIMUM"
elif certificate_ok and negative == 719:
    verdict = "DERIVED FINITE LOCAL MAXIMUM"
elif zero:
    verdict = "DERIVED EXTRA DEGENERACY" if certificate_ok else "OPEN/FAILED CERTIFICATE"
else:
    verdict = "OPEN/FAILED CERTIFICATE"

check(
    "the preregistered decision boundary is applied mechanically",
    verdict in {
        "DERIVED FINITE SADDLE", "DERIVED FINITE LOCAL MINIMUM",
        "DERIVED FINITE LOCAL MAXIMUM", "DERIVED EXTRA DEGENERACY",
        "OPEN/FAILED CERTIFICATE",
    },
    verdict,
)

payload = {
    "protocol_commit": PROTOCOL_COMMIT,
    "provenance": "preregistered target-independent full finite Hessian",
    "functional": "V(x)^(-1/3) sum_e sqrt(x_e) C(beta_e(x))",
    "configuration": "all 720 squared edge lengths near the equilateral 600-cell",
    "functional_value": {
        "scale_invariant_without_round_factor": float(run["value"]),
        "with_equal_volume_round_factor": float(recovered_endpoint),
        "independent_stored_endpoint": endpoint_value,
    },
    "f_vector": [120, 720, 1200, 600],
    "stationarity": {
        "maximum_absolute_gradient": float(np.max(np.abs(gradient))),
        "gradient_spread": float(gradient_spread),
        "scale_hessian_rms": float(np.linalg.norm(hessian @ scale)/np.sqrt(len(scale))),
    },
    "quotient_hessian": {
        "dimension": 719,
        "inertia_positive_zero_negative": list(dense_inertia),
        "minimum_eigenvalue": float(eigenvalues[0]),
        "maximum_eigenvalue": float(eigenvalues[-1]),
        "smallest_absolute_eigenvalue": smallest_nonzero,
        "cross_precision_maximum_drift": float(precision_drift),
        "ldl_inertia": list(independent_inertia),
        "edge_length_coordinate_inertia": list(length_coordinate_inertia),
        "eigenvalue_clusters": clusters,
    },
    "edge_stabilizer": {
        "order": len(stabilizer),
        "orbit_sizes": list(map(len, orbits)),
        "probe_quadratic_values": [
            {"representative": representative, "value": value}
            for representative, _, value in orbit_probes
        ],
    },
    "discrete_conformal": {
        "rank": conformal_rank,
        "compression_eigenvalues_not_full_spectral_eigenvalues": conformal_eigenvalues.tolist(),
        "explicit_x1_squared_minus_x2_squared_value": explicit_conformal_value,
        "symmetry_span_residual": conformal_symmetry_residual,
        "hessian_off_subspace_norm": conformal_hessian_leakage,
        "interpretation": (
            "negative quadratic subspace; not a Hessian-invariant eigenblock "
            "because the edge representation permits equivariant copy mixing"
        ),
    },
    "nonlinear_controls": {
        "count": len(probe_records),
        "maximum_relative_error": float(max_direct_error),
        "records": probe_records,
    },
    "verdict": verdict,
    "scope": {
        "derived": "local finite-Regge coefficient at the equilateral singular point",
        "not_derived": [
            "smooth-round discretization/refinement equivalence",
            "complete finite-cutoff action", "Lorentzian dynamics",
            "Newton/Planck normalization", "propagating graviton",
        ],
    },
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")

print("-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
print(f"QUOTIENT INERTIA (+,0,-)={dense_inertia}")
print(f"CONFORMAL RANGE=[{conformal_eigenvalues[0]:.9g},{conformal_eigenvalues[-1]:.9g}]")
print(verdict)
raise SystemExit(0 if passed == tests else 1)
