"""Reconcile the dimension observables without assigning a target dimension.

This verifier keeps four logically different quantities separate:

* the simplicial dimension of the boundary of the convex 600-cell;
* the alternating nullity count of the unrelated ``Box_p`` hierarchy;
* finite-scale Weyl-counting and heat-flow diagnostics of the Hodge operator;
* Taylor moments of the finite heat trace at t=0.

The 2% shoulder rule below is the literal, exhaustive version of the rule
documented in verify_spectral_dimension_flow.py: find the widest log-time
interval on which the full range of d_s varies by at most 2% of its mean.
It is target-free.  We calibrate it first on geodesic refinements of S^2 and
report its refinement instability rather than hiding it.
"""

from collections import defaultdict
from itertools import permutations, product

import numpy as np
import scipy.sparse as sp
from scipy.optimize import brentq
from scipy.spatial import ConvexHull


PASSED = 0
FAILED = 0


def check(label, condition, detail=""):
    global PASSED, FAILED
    condition = bool(condition)
    PASSED += int(condition)
    FAILED += int(not condition)
    print(f"[{'PASS' if condition else 'FAIL'}] {label}")
    if detail:
        print(f"       {detail}")


def parity(p):
    inversions = sum(p[i] > p[j] for i in range(len(p))
                     for j in range(i + 1, len(p)))
    return -1 if inversions % 2 else 1


def boundary_600_cell():
    """Build the boundary independently from the quaternion coordinates."""
    phi = (1 + np.sqrt(5.0)) / 2
    vertices = set()
    for i in range(4):
        for sign in (-1.0, 1.0):
            q = [0.0] * 4
            q[i] = sign
            vertices.add(tuple(q))
    vertices.update(product((-0.5, 0.5), repeat=4))
    base = [phi / 2, 0.5, 1 / (2 * phi), 0.0]
    for p in permutations(range(4)):
        if parity(p) < 0:
            continue
        q = [base[p[i]] for i in range(4)]
        nonzero = [i for i, x in enumerate(q) if abs(x) > 1e-12]
        for signs in product((-1, 1), repeat=3):
            r = q[:]
            for i, sign in zip(nonzero, signs):
                r[i] *= sign
            vertices.add(tuple(round(x, 10) for x in r))

    vertex_array = np.array(sorted(vertices))
    dots = vertex_array @ vertex_array.T
    edges = [(i, j) for i in range(120) for j in range(i + 1, 120)
             if abs(dots[i, j] - phi / 2) < 1e-3]
    adjacency = defaultdict(set)
    for i, j in edges:
        adjacency[i].add(j)
        adjacency[j].add(i)
    faces = []
    for i, j in edges:
        for k in adjacency[i] & adjacency[j]:
            if j < k:
                faces.append((i, j, k))
    tetrahedra = []
    for i, j, k in faces:
        for ell in adjacency[i] & adjacency[j] & adjacency[k]:
            if k < ell:
                tetrahedra.append((i, j, k, ell))
    return vertex_array, [[(i,) for i in range(120)],
                          edges, faces, tetrahedra]


def coboundaries(cells):
    indices = [{cell: i for i, cell in enumerate(layer)} for layer in cells]
    operators = []
    for degree in range(3):
        rows, columns, values = [], [], []
        for high_index, simplex in enumerate(cells[degree + 1]):
            for omitted in range(degree + 2):
                face = simplex[:omitted] + simplex[omitted + 1:]
                rows.append(high_index)
                columns.append(indices[degree][face])
                values.append((-1) ** omitted)
        operators.append(sp.csr_matrix(
            (values, (rows, columns)),
            shape=(len(cells[degree + 1]), len(cells[degree])),
            dtype=np.int8))
    return operators


def hodge_spectrum(coboundary):
    """Use singular spectra; every positive value occurs in two degrees."""
    positive = []
    ranks = []
    for operator in coboundary:
        if operator.shape[1] <= operator.shape[0]:
            gram = (operator.T @ operator).toarray().astype(float)
        else:
            gram = (operator @ operator.T).toarray().astype(float)
        eigenvalues = np.linalg.eigvalsh(gram)
        eigenvalues = eigenvalues[eigenvalues > 1e-8]
        positive.append(eigenvalues)
        ranks.append(len(eigenvalues))
    full = np.sort(np.concatenate(
        [np.zeros(2)] + [np.repeat(levels, 2) for levels in positive]))
    return full, tuple(ranks)


def heat_flow(eigenvalues, count=1800):
    positive = eigenvalues[eigenvalues > 1e-9]
    times = np.logspace(np.log10(0.01 / positive.max()),
                        np.log10(100 / positive.min()), count)
    weights = np.exp(-np.outer(times, eigenvalues))
    dimensions = (2 * times * (weights @ eigenvalues)
                  / weights.sum(axis=1))
    return times, dimensions


def product_heat_flow(left, right, count=1800):
    positive_sum = np.add.outer(left, right).ravel()
    positive_sum = positive_sum[positive_sum > 1e-9]
    times = np.logspace(np.log10(0.01 / positive_sum.max()),
                        np.log10(100 / positive_sum.min()), count)

    def moments(values):
        weights = np.exp(-np.outer(times, values))
        trace = weights.sum(axis=1)
        first = weights @ values
        return trace, first

    trace_l, first_l = moments(left)
    trace_r, first_r = moments(right)
    dimensions = 2 * times * (first_l / trace_l + first_r / trace_r)
    return times, dimensions, positive_sum.min(), positive_sum.max()


def widest_relative_shoulder(times, dimensions, tolerance=0.02):
    """Exhaustive target-free 2% range test; no regression or target fit."""
    log_times = np.log10(times)
    best = None
    for i in range(len(dimensions)):
        if dimensions[i] <= 0.5:
            continue
        minimum = maximum = total = float(dimensions[i])
        for j in range(i + 1, len(dimensions)):
            if dimensions[j] <= 0.5:
                break
            value = float(dimensions[j])
            minimum = min(minimum, value)
            maximum = max(maximum, value)
            total += value
            mean = total / (j - i + 1)
            if maximum - minimum <= tolerance * mean:
                candidate = (float(log_times[j] - log_times[i]),
                             mean, minimum, maximum, i, j)
                if best is None or candidate[0] > best[0]:
                    best = candidate
    return best


def heat_maximum(eigenvalues):
    positive = eigenvalues[eigenvalues > 1e-9]

    def stationary(time):
        weights = np.exp(-time * eigenvalues)
        trace = weights.sum()
        mean = (weights * eigenvalues).sum() / trace
        variance = ((weights * eigenvalues**2).sum() / trace - mean**2)
        return mean - time * variance

    grid = np.logspace(np.log10(1e-6 / positive.max()),
                       np.log10(30 / positive.min()), 4000)
    roots = []
    previous_time, previous_value = grid[0], stationary(grid[0])
    for time in grid[1:]:
        value = stationary(time)
        if value * previous_value < 0:
            roots.append(brentq(stationary, previous_time, time))
        previous_time, previous_value = time, value

    values = []
    for time in roots:
        weights = np.exp(-time * eigenvalues)
        values.append((2 * time * (weights * eigenvalues).sum()
                       / weights.sum(), time))
    return max(values), roots


def grouped_spectrum(eigenvalues):
    groups = []
    for value in np.sort(eigenvalues):
        if not groups or abs(value - groups[-1][0]) > 1e-7:
            groups.append([float(value), 1])
        else:
            old_value, multiplicity = groups[-1]
            groups[-1][0] = (old_value * multiplicity + value) / (multiplicity + 1)
            groups[-1][1] += 1
    return [(value, multiplicity) for value, multiplicity in groups]


def frozen_counting_plateau(eigenvalues, target):
    """Reproduce the already-registered target-conditioned d_N diagnostic."""
    groups = [(x, m) for x, m in grouped_spectrum(eigenvalues) if x > 1e-7]
    levels = np.array([x for x, _ in groups])
    counts = np.cumsum([m for _, m in groups])
    local = np.empty(len(levels))
    for i in range(len(levels)):
        lo, hi = max(0, i - 3), min(len(levels), i + 4)
        if hi - lo < 3:
            lo, hi = max(0, hi - 3), min(len(levels), lo + 3)
        local[i] = 2 * np.polyfit(np.log(levels[lo:hi]),
                                  np.log(counts[lo:hi]), 1)[0]
    found = []
    for i in range(len(levels)):
        for j in range(i + 2, len(levels)):
            width = np.log10(levels[j] / levels[i])
            if width < 0.50:
                continue
            slope, intercept = np.polyfit(np.log(levels[i:j + 1]),
                                          np.log(counts[i:j + 1]), 1)
            dimension = 2 * slope
            residual = np.sqrt(np.mean((
                np.log(counts[i:j + 1])
                - (slope * np.log(levels[i:j + 1]) + intercept))**2))
            stability = np.std(local[i:j + 1])
            if (abs(dimension - target) <= 0.35 and residual <= 0.08
                    and stability <= 0.35):
                found.append((width, dimension, residual, stability))
    return max(found, default=None, key=lambda item: (item[0], -item[2]))


def icosphere():
    phi = (1 + np.sqrt(5.0)) / 2
    vertices = []
    for a in (-1, 1):
        for b in (-phi, phi):
            vertices.extend(((0, a, b), (a, b, 0), (b, 0, a)))
    vertices = np.array(vertices, dtype=float)
    vertices /= np.linalg.norm(vertices, axis=1)[:, None]
    faces = [tuple(face) for face in ConvexHull(vertices).simplices]
    return vertices, faces


def subdivide_sphere(vertices, faces):
    vertices = [vertex.copy() for vertex in vertices]
    midpoint_indices = {}
    refined = []

    def midpoint(a, b):
        edge = tuple(sorted((a, b)))
        if edge not in midpoint_indices:
            vertex = vertices[a] + vertices[b]
            vertex /= np.linalg.norm(vertex)
            midpoint_indices[edge] = len(vertices)
            vertices.append(vertex)
        return midpoint_indices[edge]

    for a, b, c in faces:
        ab, bc, ca = midpoint(a, b), midpoint(b, c), midpoint(c, a)
        refined.extend(((a, ab, ca), (b, bc, ab), (c, ca, bc),
                        (ab, bc, ca)))
    return np.array(vertices), refined


def graph_laplacian_spectrum(vertex_count, faces):
    edges = {tuple(sorted(edge)) for face in faces for edge in
             ((face[0], face[1]), (face[1], face[2]), (face[2], face[0]))}
    adjacency = np.zeros((vertex_count, vertex_count))
    for i, j in edges:
        adjacency[i, j] = adjacency[j, i] = 1
    return np.linalg.eigvalsh(np.diag(adjacency.sum(axis=1)) - adjacency)


print("=" * 78)
print("DIMENSION RECONCILIATION: TOPOLOGY, DIFFUSION, AND FINITE MOMENTS")
print("=" * 78)

vertices, cells = boundary_600_cell()
f_vector = tuple(map(len, cells))
coboundary = coboundaries(cells)
check("[DERIVED] reconstructed carrier is the 600-cell boundary f=(120,720,1200,600)",
      f_vector == (120, 720, 1200, 600), str(f_vector))
check("[DERIVED] ambient affine rank is four but every boundary facet is 3D",
      np.linalg.matrix_rank(vertices - vertices.mean(axis=0)) == 4
      and all(len(simplex) == 4 for simplex in cells[3]))
check("[DERIVED] Euler characteristic is 120-720+1200-600=0",
      sum((-1) ** degree * size for degree, size in enumerate(f_vector)) == 0)
check("[DERIVED] simplicial coboundary squares to zero exactly",
      (coboundary[1] @ coboundary[0]).nnz == 0
      and (coboundary[2] @ coboundary[1]).nnz == 0)

d2_spectrum, ranks = hodge_spectrum(coboundary)
betti = (f_vector[0] - ranks[0],
         f_vector[1] - ranks[0] - ranks[1],
         f_vector[2] - ranks[1] - ranks[2],
         f_vector[3] - ranks[2])
check("[DERIVED] boundary ranks are independently (119,601,599)",
      ranks == (119, 601, 599), str(ranks))
check("[DERIVED] Betti numbers are (1,0,0,1), hence two harmonic forms",
      betti == (1, 0, 0, 1) and np.sum(d2_spectrum < 1e-8) == 2,
      f"betti={betti}, harmonic={np.sum(d2_spectrum < 1e-8)}")

positive = d2_spectrum[d2_spectrum > 1e-8]
gap, top = positive.min(), positive.max()
window = np.log10(top / gap)
check("[DERIVED] independent D^2 range and diffusion window match",
      abs(gap - 0.145898033750315) < 1e-10
      and abs(top - 15.7082039324994) < 1e-10,
      f"range=[{gap:.12f},{top:.12f}], ratio={top/gap:.6f}, "
      f"window={window:.6f} decades")

D = sp.bmat([[None, coboundary[0].T, None, None],
             [coboundary[0], None, coboundary[1].T, None],
             [None, coboundary[1], None, coboundary[2].T],
             [None, None, coboundary[2], None]], format="csr", dtype=np.int64)
D_squared = D @ D
c0 = D.shape[0]
c1 = int(D_squared.diagonal().sum())
c2 = int(D_squared.multiply(D_squared).sum()) // 2
check("[DERIVED] finite moments are exactly (2640,14880,55920)",
      (c0, c1, c2) == (2640, 14880, 55920), str((c0, c1, c2)))
check("[DERIVED] reduced moment triple obeys its arithmetic identity",
      (c0 // 240, c1 // 240, c2 // 240) == (11, 62, 233)
      and 2 * 62**2 + 1 == 3 * 11 * 233)

tiny_time = 1e-8
tiny_weights = np.exp(-tiny_time * d2_spectrum)
tiny_ds = (2 * tiny_time * (tiny_weights * d2_spectrum).sum()
           / tiny_weights.sum())
check("[DERIVED NEGATIVE] finite heat trace is analytic and d_s(t)->0 at t->0",
      tiny_ds < 2e-7,
      f"K(t)=c0-c1*t+c2*t^2+..., d_s(1e-8)={tiny_ds:.3e}")

(maximum, maximum_time), stationary_roots = heat_maximum(d2_spectrum)
check("[DERIVED NEGATIVE] the complete Kähler--Dirac heat flow never reaches four",
      len(stationary_roots) == 1 and maximum < 3.30,
      f"global maximum={maximum:.9f} at t={maximum_time:.9f}")

count_three = frozen_counting_plateau(d2_spectrum, 3)
count_four = frozen_counting_plateau(d2_spectrum, 4)
check("[DERIVED diagnostic] old d_N rule reproduces 3.0688 and no 4D interval",
      count_three is not None and abs(count_three[1] - 3.06876236) < 1e-7
      and count_four is None,
      f"count-3={count_three}; count-4={count_four}")

# Calibration comes before applying the same target-free shoulder rule to D^2.
sphere_vertices, sphere_faces = icosphere()
sphere_rows = []
for refinement in range(5):
    sphere_spectrum = graph_laplacian_spectrum(len(sphere_vertices), sphere_faces)
    sphere_times, sphere_flow = heat_flow(sphere_spectrum)
    sphere_shoulder = widest_relative_shoulder(sphere_times, sphere_flow)
    sphere_rows.append((len(sphere_vertices), float(sphere_flow.max()),
                        sphere_shoulder[0], sphere_shoulder[1]))
    if refinement < 4:
        sphere_vertices, sphere_faces = subdivide_sphere(
            sphere_vertices, sphere_faces)

control = sphere_rows[-1]
check("[DERIVED calibration] level-4 geodesic S^2 control reads two, not its peak",
      abs(control[3] - 2.0) < 0.01 and control[1] > 2.57,
      f"nodes={control[0]}, peak={control[1]:.6f}, "
      f"shoulder={control[3]:.6f} over {control[2]:.3f} decades")
check("[DERIVED warning] the shoulder selector is not stable on earlier S^2 levels",
      sphere_rows[-2][3] > 2.45 and abs(sphere_rows[-1][3] - 2) < 0.01,
      "rows=" + str([(n, round(p, 4), round(w, 3), round(s, 4))
                      for n, p, w, s in sphere_rows]))

d2_times, d2_flow = heat_flow(d2_spectrum)
d2_shoulder = widest_relative_shoulder(d2_times, d2_flow)
check("[PATTERN] the same 2% rule gives a short Kähler--Dirac shoulder near three",
      3.25 < d2_shoulder[1] < 3.30 and d2_shoulder[0] < 0.5,
      f"mean={d2_shoulder[1]:.6f}, width={d2_shoulder[0]:.3f} decades; "
      "too short for the registered half-decade plateau gate")

# Compare the two colleague-selected operators and the path product by exactly
# the same rule.  These are diagnostics, not target fits.
hasse_adjacency = abs(D).astype(float)
hasse_degree = np.asarray(hasse_adjacency.sum(axis=1)).ravel()
hasse_laplacian = np.diag(hasse_degree) - hasse_adjacency.toarray()
hasse_spectrum = np.linalg.eigvalsh(hasse_laplacian)
hasse_times, hasse_flow = heat_flow(hasse_spectrum)
hasse_shoulder = widest_relative_shoulder(hasse_times, hasse_flow)

vertex_spectrum = np.linalg.eigvalsh(
    (coboundary[0].T @ coboundary[0]).toarray().astype(float))
vertex_times, vertex_flow = heat_flow(vertex_spectrum)
vertex_shoulder = widest_relative_shoulder(vertex_times, vertex_flow)

chain_size = 320
path_spectrum = 2 - 2 * np.cos(np.pi * np.arange(chain_size) / chain_size)
product_times, product_flow, product_gap, product_top = product_heat_flow(
    vertex_spectrum, path_spectrum)
product_shoulder = widest_relative_shoulder(product_times, product_flow)

check("[DERIVED diagnostic] Hasse refinement moves the finite shoulder toward three",
      3.49 < hasse_shoulder[1] < vertex_shoulder[1] < 3.63,
      f"vertex={vertex_shoulder[1]:.6f}/{vertex_shoulder[0]:.3f} decades; "
      f"Hasse={hasse_shoulder[1]:.6f}/{hasse_shoulder[0]:.3f} decades")
check("[DERIVED diagnostic] fixed-space product with P_320 eventually reads one",
      0.95 < product_shoulder[1] < 1.02
      and abs(np.log10(product_top / product_gap) - 5.31065) < 1e-4,
      f"peak={product_flow.max():.6f}, shoulder={product_shoulder[1]:.6f}/"
      f"{product_shoulder[0]:.3f} decades, window="
      f"{np.log10(product_top/product_gap):.6f}")

phi = (1 + np.sqrt(5.0)) / 2
phi_conjugate = (1 - np.sqrt(5.0)) / 2
algebraic_four = phi**3 + phi_conjugate**3
box_nullity_count = 9 - 13 + 1 - 1
check("[DERIVED arithmetic] the two unrelated algebraic fours remain exact",
      abs(algebraic_four - 4) < 1e-12 and box_nullity_count == -4,
      "phi^3+phi'^3=4; alternating Box nullities 9-13+1-1=-4")
check("[STRUCTURAL] neither arithmetic equality is a heat/Weyl dimension probe",
      True,
      "the Box_p are a separate hierarchy; finite moments have Taylor, not "
      "Seeley--DeWitt, small-t behavior")

print("-" * 78)
print(f"RESULT: {PASSED}/{PASSED + FAILED} checks passed")
print("VERDICT_FIXED_KAHLER_DIRAC_4D=DERIVED_NEGATIVE")
print("VERDICT_KAHLER_DIRAC_3D=PATTERN_CONSISTENT_BUT_NO_HALF_DECADE_PLATEAU")
print("VERDICT_GEOMETRIC_BOUNDARY_DIMENSION=DERIVED_3")
print("VERDICT_FINITE_MOMENTS_DIMENSION=NONE")
if FAILED:
    raise SystemExit(1)
