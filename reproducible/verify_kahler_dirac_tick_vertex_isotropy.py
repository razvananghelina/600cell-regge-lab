#!/usr/bin/env python3
"""Vertex-centred metric isotropy of the signed Kähler--Dirac local tick.

Protocol commit 0015738 froze the state, ticks 0..8, geometric observable,
calibrations, thresholds and scope before this computation was performed.
Protocol commit 82ca6c9 separately froze the coordinate-free angular
multipoles and parity diagnostics before their values were computed.
"""

from collections import defaultdict, deque
from itertools import combinations
import json
from pathlib import Path
import sys

import numpy as np
import scipy.linalg as la
import scipy.sparse as sparse

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from commons import build_600cell


OUTPUT = Path(__file__).with_name("kahler_dirac_tick_vertex_isotropy.json")
PROTOCOL_COMMIT = "0015738"
MULTIPOLE_PROTOCOL_COMMIT = "82ca6c9"
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    condition = bool(condition)
    passed += int(condition)
    print(f"  [{'PASS' if condition else 'FAIL'}] {label}")
    if detail:
        print(f"         {detail}")


def sparse_max_abs(matrix):
    return 0.0 if matrix.nnz == 0 else float(np.max(np.abs(matrix.data)))


def normalized(vector):
    return vector / np.linalg.norm(vector)


def build_complex(vertices, adjacency):
    neighbours = tuple(
        frozenset(np.flatnonzero(adjacency[index]).tolist())
        for index in range(120)
    )
    edges = tuple(
        (left, right)
        for left in range(120)
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
    cells = (
        tuple((index,) for index in range(120)),
        edges,
        triangles,
        tetrahedra,
    )
    dimensions = tuple(map(len, cells))
    offsets = np.cumsum((0, *dimensions))
    indices = [
        {simplex: index for index, simplex in enumerate(layer)}
        for layer in cells
    ]
    coboundaries = []
    incidences = []
    for degree in range(3):
        rows = []
        columns = []
        values = []
        for high_index, simplex in enumerate(cells[degree + 1]):
            for omitted in range(degree + 2):
                face = simplex[:omitted] + simplex[omitted + 1:]
                low_index = indices[degree][face]
                sign = (-1) ** omitted
                rows.append(high_index)
                columns.append(low_index)
                values.append(sign)
                incidences.append((
                    int(offsets[degree] + low_index),
                    int(offsets[degree + 1] + high_index),
                    sign,
                ))
        coboundaries.append(sparse.csr_matrix(
            (values, (rows, columns)),
            shape=(dimensions[degree + 1], dimensions[degree]),
            dtype=np.int8,
        ))
    blocks = [[None] * 4 for _ in range(4)]
    for degree, coboundary in enumerate(coboundaries):
        blocks[degree + 1][degree] = coboundary
        blocks[degree][degree + 1] = coboundary.T
    kahler_dirac = sparse.bmat(blocks, format="csr", dtype=np.int8)
    return (
        cells,
        dimensions,
        offsets,
        tuple(coboundaries),
        tuple(incidences),
        kahler_dirac,
    )


def build_walk(simplex_count, incidences, degrees):
    arc_tail = []
    arc_head = []
    arc_phase = []
    lookup = {}
    for lower, higher, sign in incidences:
        for tail, head, phase in (
            (lower, higher, sign),
            (higher, lower, 1),
        ):
            lookup[(tail, head)] = len(arc_tail)
            arc_tail.append(tail)
            arc_head.append(head)
            arc_phase.append(phase)
    arc_tail = np.asarray(arc_tail, dtype=np.int32)
    arc_head = np.asarray(arc_head, dtype=np.int32)
    arc_phase = np.asarray(arc_phase, dtype=np.int8)
    arc_count = len(arc_tail)
    outgoing = [[] for _ in range(simplex_count)]
    for arc, tail in enumerate(arc_tail):
        outgoing[int(tail)].append(arc)

    reverse = np.asarray([
        lookup[(int(head), int(tail))]
        for tail, head in zip(arc_tail, arc_head)
    ], dtype=np.int32)
    shift = sparse.csr_matrix(
        (np.ones(arc_count), (reverse, np.arange(arc_count))),
        shape=(arc_count, arc_count),
    )

    rows = []
    columns = []
    data = []
    exact_blocks = True
    for arcs in outgoing:
        degree = len(arcs)
        phases = arc_phase[arcs].astype(np.int64)
        numerator = 2 * np.outer(phases, phases) - degree * np.eye(
            degree, dtype=np.int64
        )
        exact_blocks &= np.array_equal(
            numerator @ numerator,
            degree * degree * np.eye(degree, dtype=np.int64),
        )
        for local_row, row in enumerate(arcs):
            for local_column, column in enumerate(arcs):
                value = numerator[local_row, local_column] / degree
                if value:
                    rows.append(row)
                    columns.append(column)
                    data.append(value)
    coin = sparse.csr_matrix(
        (data, (rows, columns)), shape=(arc_count, arc_count)
    )
    embedding = sparse.csr_matrix(
        (
            arc_phase / np.sqrt(degrees[arc_tail]),
            (np.arange(arc_count), arc_tail),
        ),
        shape=(arc_count, simplex_count),
    )
    return {
        "arc_tail": arc_tail,
        "arc_head": arc_head,
        "outgoing": outgoing,
        "shift": shift,
        "coin": coin,
        "walk": shift @ coin,
        "embedding": embedding,
        "exact_coin_blocks": exact_blocks,
    }


def simplex_centres(cells, vertices):
    return np.asarray([
        normalized(vertices[list(simplex)].sum(axis=0))
        for layer in cells
        for simplex in layer
    ])


def tangent_coordinates(base, centres):
    tangent_basis = la.null_space(base.reshape(1, 4))
    cosines = np.clip(centres @ base, -1.0, 1.0)
    angles = np.arccos(cosines)
    tangents = centres - cosines[:, None] * base
    norms = np.linalg.norm(tangents, axis=1)
    logs = np.zeros_like(tangents)
    regular = (angles > 1e-14) & (norms > 1e-12)
    antipodal = (angles > np.pi - 1e-10) & (norms <= 1e-12)
    logs[regular] = (
        angles[regular, None] * tangents[regular]
        / norms[regular, None]
    )
    return logs @ tangent_basis, antipodal


def moment_audit(probabilities, coordinates, degree_masks):
    mean = probabilities @ coordinates
    raw = (coordinates.T * probabilities) @ coordinates
    covariance = raw - np.outer(mean, mean)
    eigenvalues = np.linalg.eigvalsh(covariance)
    trace = float(np.trace(covariance))
    if trace > 1e-20 and eigenvalues[0] > 1e-20:
        ratio = float(eigenvalues[-1] / eigenvalues[0])
        residual = float(
            np.linalg.norm(covariance - trace / 3 * np.eye(3)) / trace
        )
        gate = bool(
            np.linalg.norm(mean) < 1e-11
            and abs(ratio - 1.0) < 1e-10
            and residual < 1e-10
        )
    else:
        ratio = None
        residual = None
        gate = None
    radial = np.linalg.norm(coordinates, axis=1)
    support = np.flatnonzero(probabilities > 1e-24)
    return {
        "probability_sum": float(probabilities.sum()),
        "support_simplex_count": int(len(support)),
        "form_degree_probabilities": [
            float(probabilities[mask].sum()) for mask in degree_masks
        ],
        "mean_tangent_vector": mean.tolist(),
        "mean_norm": float(np.linalg.norm(mean)),
        "covariance": covariance.tolist(),
        "covariance_eigenvalues": eigenvalues.tolist(),
        "eigenvalue_ratio": ratio,
        "normalized_traceless_residual": residual,
        "isotropy_gate": gate,
        "mean_radius": float(probabilities @ radial),
        "rms_radius": float(np.sqrt(probabilities @ (radial ** 2))),
        "maximum_support_radius": float(radial[support].max()) if len(support) else 0.0,
    }


def lifted_polar_geometry(vertices, cells, base_vertex):
    """Lift rounded Q(sqrt(5)) vertices and form tangent polar data.

    commons.cell600 deliberately stores coordinates rounded to ten decimals.
    That is ample for adjacency and second moments but loses roughly ten
    digits under a theoretically exact harmonic cancellation followed by a
    square root.  Every coordinate belongs to a five-element Q(sqrt(5))
    alphabet, so restore that alphabet in long-double arithmetic.
    """
    scalar_type = np.longdouble
    phi = (scalar_type(1) + np.sqrt(scalar_type(5))) / scalar_type(2)
    alphabet = np.asarray(
        [scalar_type(0), scalar_type(0.5), scalar_type(1) / (2 * phi),
         phi / 2, scalar_type(1)],
        dtype=scalar_type,
    )
    lifted = np.empty(vertices.shape, dtype=scalar_type)
    maximum_lift = 0.0
    for index, value in np.ndenumerate(vertices):
        closest = alphabet[np.argmin(np.abs(alphabet - abs(value)))]
        restored = np.copysign(closest, scalar_type(value)) if value else closest
        lifted[index] = restored
        maximum_lift = max(maximum_lift, abs(float(restored) - value))
    centres = []
    for layer in cells:
        for simplex in layer:
            centre = lifted[list(simplex)].sum(axis=0)
            centre /= np.sqrt(np.sum(centre * centre))
            centres.append(centre)
    centres = np.asarray(centres, dtype=scalar_type)
    base = lifted[base_vertex]
    cosines = np.clip(centres @ base, scalar_type(-1), scalar_type(1))
    radii = np.arccos(cosines)
    tangent = centres - cosines[:, None] * base
    tangent_norms = np.sqrt(np.sum(tangent * tangent, axis=1))
    directions = np.zeros_like(tangent)
    regular = tangent_norms > scalar_type("1e-18")
    directions[regular] = tangent[regular] / tangent_norms[regular, None]
    return {
        "radii": radii,
        "directions": directions,
        "maximum_coordinate_lift": maximum_lift,
    }


def harmonic_multipoles(
    probabilities, coordinates, maximum_degree=12, polar_geometry=None
):
    """Return three rotation-invariant angular multipole normalizations.

    The base point has no angular direction and is removed.  Conditional and
    unconditional angular powers use the same numerator; the former divides
    by the moving probability while the latter deliberately does not.  The
    solid-harmonic variant weights each direction by p r**ell.
    """
    scalar_type = np.longdouble
    if polar_geometry is None:
        radii = np.linalg.norm(coordinates, axis=1).astype(scalar_type)
        all_directions = np.zeros_like(coordinates, dtype=scalar_type)
        regular = radii > scalar_type("1e-13")
        all_directions[regular] = (
            coordinates[regular] / radii[regular, None]
        )
    else:
        radii = polar_geometry["radii"]
        all_directions = polar_geometry["directions"]
    moving = (radii > 1e-13) & (probabilities > 1e-24)
    directions = all_directions[moving]
    angular_weights = probabilities[moving].astype(scalar_type)
    moving_mass = float(angular_weights.sum())
    dots = np.clip(
        directions @ directions.T, scalar_type(-1), scalar_type(1)
    )
    variants = {
        "conditional_angular": [],
        "unconditional_angular": [],
        "solid_harmonic_radial": [],
    }
    minimum_raw_squared = 0.0
    previous_kernel = np.ones_like(dots)
    kernel = dots.copy()
    for ell in range(1, maximum_degree + 1):
        if ell > 1:
            next_kernel = (
                (2 * ell - 1) * dots * kernel
                - (ell - 1) * previous_kernel
            ) / ell
            previous_kernel, kernel = kernel, next_kernel
        angular_squared = float(angular_weights @ kernel @ angular_weights)
        radial_weights = angular_weights * radii[moving] ** ell
        radial_mass = float(radial_weights.sum())
        radial_squared = float(radial_weights @ kernel @ radial_weights)
        minimum_raw_squared = min(
            minimum_raw_squared, angular_squared, radial_squared
        )
        angular_power = np.sqrt(max(0.0, angular_squared))
        radial_power = np.sqrt(max(0.0, radial_squared))
        variants["conditional_angular"].append(
            float(angular_power / moving_mass)
        )
        variants["unconditional_angular"].append(float(angular_power))
        variants["solid_harmonic_radial"].append(
            float(radial_power / radial_mass)
        )
    return {
        "moving_probability": moving_mass,
        "direction_count": int(moving.sum()),
        "minimum_raw_squared_power": minimum_raw_squared,
        "by_weighting": variants,
    }


def icosahedral_invariant_multiplicities(maximum_degree=12):
    """Multiplicity of the A5 trivial irrep in each SO(3) harmonic."""
    angles_and_counts = (
        (0.0, 1),
        (np.pi, 15),
        (2 * np.pi / 3, 20),
        (2 * np.pi / 5, 12),
        (4 * np.pi / 5, 12),
    )
    raw = []
    rounded = []
    for ell in range(maximum_degree + 1):
        total = 0.0
        for angle, count in angles_and_counts:
            character = (
                2 * ell + 1
                if angle == 0
                else np.sin((ell + 0.5) * angle) / np.sin(angle / 2)
            )
            total += count * character
        value = total / 60
        raw.append(float(value))
        rounded.append(int(np.rint(value)))
    return raw, rounded


print("=" * 78)
print("VERTEX-CENTRED ISOTROPY OF THE KAEHLER--DIRAC LOCAL TICK")
print("=" * 78)

vertices, adjacency, _ = build_600cell()
(
    cells,
    dimensions,
    offsets,
    coboundaries,
    incidences,
    kahler_dirac,
) = build_complex(vertices, adjacency)
check("the oriented 600-cell cochain f-vector is exact",
      dimensions == (120, 720, 1200, 600))
check("both consecutive coboundary products vanish exactly over Z",
      (coboundaries[1] @ coboundaries[0]).nnz == 0
      and (coboundaries[2] @ coboundaries[1]).nnz == 0)

degrees = np.asarray(np.abs(kahler_dirac).sum(axis=1)).ravel().astype(int)
walk_data = build_walk(kahler_dirac.shape[0], incidences, degrees)
shift = walk_data["shift"]
coin = walk_data["coin"]
walk = walk_data["walk"]
embedding = walk_data["embedding"]
check("all local reflection blocks square exactly before division",
      walk_data["exact_coin_blocks"])
check("the reversal and assembled tick are unitary",
      sparse_max_abs(shift.T @ shift - sparse.eye(shift.shape[0])) == 0
      and sparse_max_abs(coin.T @ coin - sparse.eye(coin.shape[0])) < 3e-16
      and sparse_max_abs(walk.T @ walk - sparse.eye(walk.shape[0])) < 3e-16)

inverse_sqrt_degree = sparse.diags(1 / np.sqrt(degrees))
discriminant = embedding.T @ shift @ embedding
normalized_dirac = inverse_sqrt_degree @ kahler_dirac @ inverse_sqrt_degree
check("the walk discriminant is the normalized signed Kähler--Dirac operator",
      sparse_max_abs(discriminant - normalized_dirac) < 3e-16)

first_intertwiner = walk @ embedding - shift @ embedding
second_intertwiner = (
    walk @ shift @ embedding
    - (2 * shift @ embedding @ discriminant - embedding)
)
check("both exact invariant-sector intertwining identities hold numerically",
      sparse_max_abs(first_intertwiner) < 3e-16
      and sparse_max_abs(second_intertwiner) < 8e-16)

centres = simplex_centres(cells, vertices)
base_vertex = 0
base = vertices[base_vertex]
coordinates, antipodal = tangent_coordinates(base, centres)
multipole_polar_geometry = lifted_polar_geometry(
    vertices, cells, base_vertex
)
check("the rounded coordinates lift uniquely back to Q(sqrt(5))",
      multipole_polar_geometry["maximum_coordinate_lift"] < 5e-11,
      "maximum lift={:.3e}".format(
          multipole_polar_geometry["maximum_coordinate_lift"]
      ))
degree_masks = [
    np.arange(int(offsets[degree]), int(offsets[degree + 1]))
    for degree in range(4)
]

incidence_lengths = defaultdict(list)
for lower, higher, _sign in incidences:
    lower_degree = int(np.searchsorted(offsets, lower, side="right") - 1)
    cosine = float(np.clip(centres[lower] @ centres[higher], -1.0, 1.0))
    incidence_lengths[f"{lower_degree}->{lower_degree + 1}"].append(
        float(np.arccos(cosine))
    )
incidence_length_audit = {
    rank_pair: {
        "count": len(values),
        "mean": float(np.mean(values)),
        "minimum": float(np.min(values)),
        "maximum": float(np.max(values)),
        "spread": float(np.ptp(values)),
    }
    for rank_pair, values in sorted(incidence_lengths.items())
}
check("each rank-pair incidence has one symmetry-fixed geodesic length",
      all(audit["spread"] < 2e-9
          for audit in incidence_length_audit.values()),
      "lengths=" + str({
          pair: audit["mean"] for pair, audit in incidence_length_audit.items()
      }))
rank_pair_means = [
    audit["mean"] for audit in incidence_length_audit.values()
]
check("the three Hasse microstep lengths are not one common physical length",
      np.ptp(rank_pair_means) > 1e-3,
      f"spread among rank means={np.ptp(rank_pair_means):.6g}")

# Positive and deliberately distorted known-answer controls on the 12-neighbour
# icosahedral tangent shell.
neighbour_vertices = np.flatnonzero(adjacency[base_vertex]).astype(int)
shell = coordinates[neighbour_vertices]
shell_lengths = np.linalg.norm(shell, axis=1)
shell_directions = shell / shell_lengths[:, None]
uniform = np.full(12, 1 / 12)
control_mean = uniform @ shell_directions
control_covariance = (shell_directions.T * uniform) @ shell_directions
control_eigenvalues = np.linalg.eigvalsh(control_covariance)
control_residual = np.linalg.norm(
    control_covariance - np.eye(3) / 3
) / np.trace(control_covariance)
check("the 12-neighbour vertex figure is an equal-length isotropic shell",
      np.ptp(shell_lengths) < 2e-10
      and np.linalg.norm(control_mean) < 1e-14
      and abs(control_eigenvalues[-1] / control_eigenvalues[0] - 1) < 1e-13
      and control_residual < 1e-13,
      f"length spread={np.ptp(shell_lengths):.3e}, residual={control_residual:.3e}")

distorted = np.ones(12)
distorted[0] = 2
distorted /= distorted.sum()
distorted_mean = distorted @ shell_directions
distorted_raw = (shell_directions.T * distorted) @ shell_directions
distorted_covariance = distorted_raw - np.outer(distorted_mean, distorted_mean)
distorted_eigenvalues = np.linalg.eigvalsh(distorted_covariance)
distorted_residual = np.linalg.norm(
    distorted_covariance - np.trace(distorted_covariance) / 3 * np.eye(3)
) / np.trace(distorted_covariance)
check("the frozen one-direction perturbation is correctly rejected",
      np.linalg.norm(distorted_mean) > 1e-3
      and abs(distorted_eigenvalues[-1] / distorted_eigenvalues[0] - 1) > 1e-3
      and distorted_residual > 1e-3)

# Group-theoretic and known-shell controls for the separately preregistered
# angular multipole audit.
invariant_raw, invariant_multiplicities = (
    icosahedral_invariant_multiplicities(12)
)
check("the A5 character average is integral through ell=12",
      max(abs(raw - rounded) for raw, rounded in zip(
          invariant_raw, invariant_multiplicities
      )) < 2e-14,
      f"multiplicities={invariant_multiplicities}")
check("ell=6 is the first nonconstant icosahedrally invariant harmonic",
      invariant_multiplicities[1:6] == [0] * 5
      and invariant_multiplicities[6] == 1)

shell_probabilities = np.zeros(kahler_dirac.shape[0])
shell_probabilities[neighbour_vertices] = uniform
shell_multipoles = harmonic_multipoles(
    shell_probabilities, coordinates,
    polar_geometry=multipole_polar_geometry,
)
distorted_shell_probabilities = np.zeros(kahler_dirac.shape[0])
distorted_shell_probabilities[neighbour_vertices] = distorted
distorted_shell_multipoles = harmonic_multipoles(
    distorted_shell_probabilities, coordinates,
    polar_geometry=multipole_polar_geometry,
)
shell_conditional = shell_multipoles["by_weighting"]["conditional_angular"]
distorted_conditional = distorted_shell_multipoles["by_weighting"][
    "conditional_angular"
]
check("the uniform icosahedral shell first appears at angular degree six",
      max(shell_conditional[:5]) < 1e-10
      and shell_conditional[5] > 1e-6,
      "A1..A6=" + str(shell_conditional[:6]))
check("the frozen distorted shell produces a detectable lower multipole",
      max(distorted_conditional[:5]) > 1e-3,
      "max(A1..A5)={:.6g}".format(max(distorted_conditional[:5])))

# Exact Hasse distances from the initial vertex for the strict cone audit.
hasse_neighbours = [[] for _ in range(kahler_dirac.shape[0])]
rows, columns = kahler_dirac.nonzero()
for row, column in zip(rows, columns):
    hasse_neighbours[int(row)].append(int(column))
distances = np.full(kahler_dirac.shape[0], -1, dtype=int)
distances[base_vertex] = 0
queue = deque((base_vertex,))
while queue:
    current = queue.popleft()
    for neighbour in hasse_neighbours[current]:
        if distances[neighbour] < 0:
            distances[neighbour] = distances[current] + 1
            queue.append(neighbour)

state = np.asarray(embedding[:, base_vertex].todense()).ravel()
audits = []
norms = []
cone_violations = []
antipodal_probability = []
for tick in range(9):
    norms.append(float(np.vdot(state, state).real))
    probabilities = np.bincount(
        walk_data["arc_tail"],
        weights=np.abs(state) ** 2,
        minlength=kahler_dirac.shape[0],
    )
    audit = moment_audit(probabilities, coordinates, degree_masks)
    if tick > 0:
        audit["angular_multipoles"] = harmonic_multipoles(
            probabilities, coordinates,
            polar_geometry=multipole_polar_geometry,
        )
    audit["tick"] = tick
    audits.append(audit)
    occupied_arcs = np.flatnonzero(np.abs(state) > 1e-13)
    occupied_tails = np.unique(walk_data["arc_tail"][occupied_arcs])
    cone_violations.append(int(np.sum(distances[occupied_tails] > tick)))
    antipodal_probability.append(float(probabilities[antipodal].sum()))
    state = walk @ state

check("all nine evolved norms and simplex probability sums equal one",
      max(abs(value - 1) for value in norms) < 2e-14
      and max(abs(audit["probability_sum"] - 1) for audit in audits) < 2e-14)
check("the strict Hasse cone holds at every preregistered tick",
      cone_violations == [0] * 9,
      f"violations by tick={cone_violations}")
check("no occupied probability reaches the antipodal log ambiguity",
      max(antipodal_probability) < 1e-24,
      f"antipodal probabilities={antipodal_probability}")
check("the form-degree parity flips at every micro-tick",
      all(
          sum(probability for degree, probability in enumerate(
              audit["form_degree_probabilities"]
          ) if degree % 2 != tick % 2) < 2e-14
          for tick, audit in enumerate(audits)
      ))

nonzero_time = audits[1:]
isotropy_hits = sum(audit["isotropy_gate"] is True for audit in nonzero_time)
check("all eight preregistered nonzero ticks pass tangent isotropy",
      isotropy_hits == 8,
      "hits={}/8; max mean={:.3e}; max residual={:.3e}".format(
          isotropy_hits,
          max(audit["mean_norm"] for audit in nonzero_time),
          max(audit["normalized_traceless_residual"] for audit in nonzero_time),
      ))

multipole_raw_minimum = min(
    audit["angular_multipoles"]["minimum_raw_squared_power"]
    for audit in nonzero_time
)
check("all raw squared multipole powers are nonnegative within tolerance",
      multipole_raw_minimum > -1e-12,
      f"minimum raw squared power={multipole_raw_minimum:.3e}")
lower_multipole_maxima = {}
for weighting in (
    "conditional_angular",
    "unconditional_angular",
    "solid_harmonic_radial",
):
    lower_multipole_maxima[weighting] = max(
        max(audit["angular_multipoles"]["by_weighting"][weighting][:5])
        for audit in nonzero_time
    )
check("A1 through A5 vanish for every tick and frozen weighting",
      max(lower_multipole_maxima.values()) < 1e-9,
      f"maxima={lower_multipole_maxima}")

parity_diagnostics = {}
for weighting in (
    "conditional_angular",
    "unconditional_angular",
    "solid_harmonic_radial",
):
    values = [
        audit["angular_multipoles"]["by_weighting"][weighting][5]
        for audit in nonzero_time
    ]
    odd = values[0::2]
    even = values[1::2]
    separation_ratio = float(max(even) / min(odd))
    separated = bool(max(even) < min(odd))
    if separation_ratio < 0.1:
        status = "DERIVED parity separation"
    elif separated:
        status = "PATTERN: weak parity separation"
    else:
        status = "DERIVED NEGATIVE: even and odd ranges overlap"
    parity_diagnostics[weighting] = {
        "A6_by_tick_1_to_8": values,
        "A6_squared_power_by_tick_1_to_8": [
            float(value ** 2) for value in values
        ],
        "odd_minimum": float(min(odd)),
        "odd_maximum": float(max(odd)),
        "even_minimum": float(min(even)),
        "even_maximum": float(max(even)),
        "separation_ratio_max_even_over_min_odd": separation_ratio,
        "even_strictly_below_odd": separated,
        "status": status,
    }

maximum_step_length = max(
    audit["maximum"] for audit in incidence_length_audit.values()
)
front_residuals = [
    abs(audit["maximum_support_radius"] - audit["tick"] * maximum_step_length)
    for audit in audits
]
check("the causal front numerically saturates the exact pi/10 angular bound",
      abs(maximum_step_length - np.pi / 10) < 2e-9
      and max(front_residuals) < 2e-9,
      f"c_angle={maximum_step_length:.15f}, max residual={max(front_residuals):.3e}")

payload = {
    "protocol_commit": PROTOCOL_COMMIT,
    "multipole_protocol_commit": MULTIPOLE_PROTOCOL_COMMIT,
    "phenomenological_target_used": False,
    "carrier": {
        "cochain_dimension": int(kahler_dirac.shape[0]),
        "directed_hasse_arc_dimension": int(walk.shape[0]),
        "initial_state": "A e_0",
        "initial_vertex": base_vertex,
    },
    "control": {
        "vertex_neighbour_count": len(neighbour_vertices),
        "shell_length_mean": float(shell_lengths.mean()),
        "shell_length_spread": float(np.ptp(shell_lengths)),
        "uniform_mean_norm": float(np.linalg.norm(control_mean)),
        "uniform_covariance_eigenvalues": control_eigenvalues.tolist(),
        "uniform_normalized_traceless_residual": float(control_residual),
        "distorted_mean_norm": float(np.linalg.norm(distorted_mean)),
        "distorted_covariance_eigenvalues": distorted_eigenvalues.tolist(),
        "distorted_normalized_traceless_residual": float(distorted_residual),
        "icosahedral_invariant_multiplicities_ell_0_to_12": (
            invariant_multiplicities
        ),
        "uniform_shell_conditional_multipoles_ell_1_to_12": (
            shell_conditional
        ),
        "distorted_shell_conditional_multipoles_ell_1_to_12": (
            distorted_conditional
        ),
    },
    "rank_pair_incidence_geodesic_lengths": incidence_length_audit,
    "maximum_single_tick_geodesic_length": maximum_step_length,
    "dimensionless_angular_front_speed": maximum_step_length,
    "angular_speed_bound_exact_identification": "pi/10 per micro-tick",
    "front_saturation_status": (
        "DERIVED NUMERICAL through tick 8 at stored-coordinate precision"
    ),
    "front_saturation_residuals_by_tick": front_residuals,
    "tick_count_nonzero_N": 8,
    "isotropy_hit_fraction": [isotropy_hits, 8],
    "norms": norms,
    "cone_violations_by_tick": cone_violations,
    "antipodal_simplex_count": int(antipodal.sum()),
    "antipodal_probability_by_tick": antipodal_probability,
    "tick_audits": audits,
    "multipole_definition": (
        "sqrt(sum_ij q_i q_j P_ell(u_i dot u_j)) / sum_i q_i; "
        "unconditional variant omits the denominator"
    ),
    "multipole_parity_diagnostics": parity_diagnostics,
    "verdict": (
        "DERIVED VERTEX-CENTRED ISOTROPIC TICK"
        if isotropy_hits == 8
        else "DERIVED NEGATIVE: VERTEX-CENTRED ISOTROPY FAILURE"
    ),
    "scope": (
        "Fixed 600-cell, canonical embedded vertex state, ticks 1..8; "
        "not a continuum, generic-state or physical-speed theorem."
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
check("the structured vertex-isotropy certificate was written", OUTPUT.exists())

print("\n" + "-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
print(f"VERTEX_ISOTROPY_HITS={isotropy_hits}/8")
for audit in audits[1:]:
    print(
        "n={tick}: support={support_simplex_count}, mean={mean_norm:.3e}, "
        "R={eigenvalue_ratio:.12f}, A={normalized_traceless_residual:.3e}, "
        "r_rms={rms_radius:.6f}".format(**audit)
    )
print("A6 PARITY DIAGNOSTICS (ticks 1..8)")
for weighting, diagnostic in parity_diagnostics.items():
    print(
        f"  {weighting}: A6={diagnostic['A6_by_tick_1_to_8']}, "
        f"A6^2={diagnostic['A6_squared_power_by_tick_1_to_8']}, "
        f"R_sep={diagnostic['separation_ratio_max_even_over_min_odd']:.6g}, "
        f"{diagnostic['status']}"
    )
print("SCOPE: vertex-centred fixed-complex kinematics only.")
raise SystemExit(0 if passed == tests else 1)
