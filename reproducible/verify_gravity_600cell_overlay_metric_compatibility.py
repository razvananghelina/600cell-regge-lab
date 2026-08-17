#!/usr/bin/env python3
"""Exact metric-compatibility audit for the universal staircase overlay."""

from collections import Counter
from functools import reduce
from itertools import combinations, permutations
import hashlib
import json
from pathlib import Path

import sympy as sp


HERE = Path(__file__).resolve().parent
OVERLAY_SOURCE = HERE / "gravity_600cell_universal_staircase_overlay.json"
FACE_SOURCE = HERE / "gravity_600cell_overlay_face_poset.json"
OUTPUT = HERE / "gravity_600cell_overlay_metric_compatibility.json"
OVERLAY_SHA256 = "0dd03eed878f599463a44160484c74ddeaa0511fc70c8b2e77bc05a2f36dd3dc"
FACE_SHA256 = "439a3d067d50415f0a47c79091ec746c12dd7975b2246b6143f3f7a70847ce13"
PRIOR_ART_COMMIT = "5443238"
PROTOCOL_COMMIT = "5735a15"
PROTOCOL_CLARIFICATION_COMMIT = "0d2b8b4"

VERTICES = tuple(range(4))
ORDERS = tuple(permutations(VERTICES))
FULL_MASK = 15
MASKS = tuple(range(1, FULL_MASK))
MASK_INDEX = {mask: index for index, mask in enumerate(MASKS)}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def mask_of(items):
    result = 0
    for item in items:
        result |= 1 << item
    return result


def parse_word(word):
    return tuple(1 if character == "+" else -1 for character in word)


def word_of(pattern):
    return "".join("+" if sign > 0 else "-" for sign in pattern)


def expression_key(expression):
    return sp.srepr(sp.expand(expression))


def matrix_key(matrix):
    return tuple(expression_key(entry) for entry in matrix)


def matrix_equal(left, right):
    return left.shape == right.shape and all(
        sp.expand(a-b) == 0 for a, b in zip(left, right)
    )


def affine_coefficients(expression, coordinates):
    zero = {coordinate: 0 for coordinate in coordinates}
    coefficients = [sp.expand(expression.subs(zero))]
    coefficients.extend(sp.expand(sp.diff(expression, coordinate))
                        for coordinate in coordinates)
    reconstructed = coefficients[0]+sum(
        coefficient*coordinate
        for coefficient, coordinate in zip(coefficients[1:], coordinates)
    )
    if sp.expand(expression-reconstructed) != 0:
        raise ValueError("expression is not affine")
    return coefficients


tests = []


def check(label, condition):
    tests.append((label, bool(condition)))


overlay = json.loads(OVERLAY_SOURCE.read_text())
faces = json.loads(FACE_SOURCE.read_text())
words = tuple(overlay.get("feasible_sign_words", ()))
patterns = tuple(parse_word(word) for word in words)
face_words = tuple(record["sign_word"] for record in faces.get("chambers", ()))

overlay_source_ok = bool(
    digest(OVERLAY_SOURCE) == OVERLAY_SHA256
    and overlay.get("outcome") == "UNIVERSAL_STAIRCASE_OVERLAY_CERTIFIED"
    and overlay.get("passed") == overlay.get("tests") == 12
    and overlay.get("prior_art_commit") == "4fb243b"
    and overlay.get("protocol_commit") == "4dc2feb"
    and overlay.get("subset_masks") == list(MASKS)
    and overlay.get("full_dimensional_chamber_count") == 148
    and len(words) == len(set(words)) == 148
)
face_source_ok = bool(
    digest(FACE_SOURCE) == FACE_SHA256
    and faces.get("outcome") == "UNIVERSAL_OVERLAY_FACE_POSET_CERTIFIED"
    and faces.get("passed") == faces.get("tests") == 13
    and faces.get("arrangement_vertex_count") == 33
    and faces.get("polyhedral_f_vector_local") == [33, 206, 468, 442, 148]
    and set(face_words) == set(words)
)
check("the two frozen overlay artifacts and hashes are reproduced",
      overlay_source_ok and face_source_ok)


# Exact symbols and the regular 600-cell tetrahedron Gram matrix.
R_MINUS, R_PLUS, T, ratio = sp.symbols(
    "R_minus R_plus T r", positive=True
)
SQRT5 = sp.sqrt(5)
PHI = (1+SQRT5)/2
ADJACENT_DOT = PHI/2
U = sp.Matrix(4, 4, lambda i, j: 1 if i == j else ADJACENT_DOT)
ETA = sp.diag(1, 1, 1, 1, -1)
ETA[:4, :4] = U

lambda0, lambda1, lambda2, time = sp.symbols(
    "lambda_0 lambda_1 lambda_2 t"
)
COORDINATES = (lambda0, lambda1, lambda2, time)
LAMBDAS = (
    lambda0,
    lambda1,
    lambda2,
    1-lambda0-lambda1-lambda2,
)


def parameter_vertex(vertex, top):
    lambdas = [sp.Integer(0)]*4
    lambdas[vertex] = sp.Integer(1)
    return (sp.Integer(1), *lambdas[:3], sp.Integer(1 if top else 0))


def target_vertex(vertex, top):
    scale = R_PLUS if top else R_MINUS
    spatial = [sp.Integer(0)]*4
    spatial[vertex] = scale
    return (*spatial, T if top else sp.Integer(0))


def simplex_vertices(order, split):
    return tuple(
        [(vertex, False) for vertex in order[:split+1]]
        + [(vertex, True) for vertex in order[split:]]
    )


def solved_affine_map(order, split):
    labelled_vertices = simplex_vertices(order, split)
    domain = sp.Matrix([
        parameter_vertex(vertex, top) for vertex, top in labelled_vertices
    ])
    target = sp.Matrix([
        target_vertex(vertex, top) for vertex, top in labelled_vertices
    ])
    determinant = sp.det(domain)
    if determinant == 0:
        raise ValueError("singular staircase simplex")
    # Rows are target coordinates; columns are (constant, lambda0, ..., t).
    forms = (domain.inv()*target).T.applyfunc(sp.expand)
    return determinant, forms


def closed_affine_map(order, split):
    pivot = order[split]
    before = order[:split]
    after = order[split+1:]
    tail_after = sum((LAMBDAS[index] for index in after), sp.Integer(0))
    beta_pivot = time-tail_after
    alpha_pivot = LAMBDAS[pivot]-beta_pivot
    spatial = [sp.Integer(0)]*4
    for vertex in before:
        spatial[vertex] = R_MINUS*LAMBDAS[vertex]
    for vertex in after:
        spatial[vertex] = R_PLUS*LAMBDAS[vertex]
    spatial[pivot] = R_MINUS*alpha_pivot+R_PLUS*beta_pivot
    expressions = (*spatial, T*time)
    return sp.Matrix([
        affine_coefficients(expression, COORDINATES)
        for expression in expressions
    ])


affine_maps = {}
metrics = {}
determinants = {}
closed_form_matches = []
metric_symmetry = []
for order in ORDERS:
    for split in range(4):
        determinant, forms = solved_affine_map(order, split)
        closed = closed_affine_map(order, split)
        key = (order, split)
        affine_maps[key] = forms
        determinants[key] = determinant
        closed_form_matches.append(matrix_equal(forms, closed))
        jacobian = forms[:, 1:]
        metric = (jacobian.T*ETA*jacobian).applyfunc(sp.expand)
        metrics[key] = metric
        metric_symmetry.append(matrix_equal(metric, metric.T))

affine_reconstruction_ok = bool(
    len(affine_maps) == len(metrics) == 96
    and all(value != 0 for value in determinants.values())
    and all(closed_form_matches)
)
check("all 96 affine maps are independently reconstructed exactly",
      affine_reconstruction_ok)


# Adjacent staircase maps must agree on their common internal facet.
continuity_failures = []
for order in ORDERS:
    for split in range(3):
        wall = sum(
            (LAMBDAS[index] for index in order[split+1:]),
            sp.Integer(0),
        )
        left = affine_maps[(order, split)]
        right = affine_maps[(order, split+1)]
        for target_index in range(5):
            difference = left[target_index, :]-right[target_index, :]
            expression = difference[0]+sum(
                coefficient*coordinate
                for coefficient, coordinate in zip(
                    difference[1:], COORDINATES
                )
            )
            restricted = sp.expand(expression.subs(time, wall))
            if restricted != 0:
                continuity_failures.append({
                    "order": list(order),
                    "split": split,
                    "target_coordinate": target_index,
                    "residual": str(restricted),
                })
continuity_ok = not continuity_failures
check("all 72 staircase internal facets pass exact PL continuity",
      continuity_ok)


# Every schedule has the same time map T*t, so lapse-dependent terms must
# cancel from all inter-schedule metric differences.  Check this globally.
metric_types = {}
for key, metric in metrics.items():
    metric_types.setdefault(matrix_key(metric), (key, metric))
metric_type_records = tuple(metric_types.values())
t_cancellation_failures = []
for (_, left), (_, right) in combinations(metric_type_records, 2):
    for row in range(4):
        for column in range(4):
            difference = sp.expand(left[row, column]-right[row, column])
            if difference.has(T):
                t_cancellation_failures.append((row, column, str(difference)))
metric_controls_ok = bool(all(metric_symmetry) and not t_cancellation_failures)
check("all pullback metrics are symmetric and lapse cancels from differences",
      metric_controls_ok)


# Reconstruct each chamber/order assignment without trusting recorded counts.
def assignments(pattern, order):
    result = []
    for split in range(4):
        lower_mask = mask_of(order[split+1:])
        upper_mask = mask_of(order[split:])
        lower_ok = lower_mask == 0 or pattern[MASK_INDEX[lower_mask]] > 0
        upper_ok = upper_mask == FULL_MASK or pattern[MASK_INDEX[upper_mask]] < 0
        if lower_ok and upper_ok:
            result.append(split)
    return tuple(result)


assignment_table = []
assignment_ok = True
computed_order_records = []
for order_index, order in enumerate(ORDERS):
    split_counts = Counter()
    multiplicities = Counter()
    for pattern in patterns:
        found = assignments(pattern, order)
        multiplicities[len(found)] += 1
        if len(found) == 1:
            split_counts[found[0]] += 1
        else:
            assignment_ok = False
    recorded = overlay["staircase_orders"][order_index]
    assignment_ok &= bool(
        tuple(recorded["order"]) == order
        and recorded["simplex_assignment_counts"]
        == {str(split): split_counts[split] for split in range(4)}
        and recorded["assignment_multiplicity_distribution"]
        == {str(key): value for key, value in sorted(multiplicities.items())}
    )
    computed_order_records.append({
        "order": list(order),
        "simplex_assignment_counts": [split_counts[index] for index in range(4)],
    })

for pattern in patterns:
    row = []
    for order in ORDERS:
        found = assignments(pattern, order)
        row.append(found[0] if len(found) == 1 else None)
    assignment_table.append(tuple(row))

assignment_ok &= bool(
    len(assignment_table) == 148
    and all(split is not None for row in assignment_table for split in row)
)
check("all 148 x 24 chamber assignments are independently unique",
      assignment_ok)


# Static global-affine control.
global_expressions = tuple(R_MINUS*value for value in LAMBDAS)+(T*time,)
global_forms = sp.Matrix([
    affine_coefficients(expression, COORDINATES)
    for expression in global_expressions
])
global_jacobian = global_forms[:, 1:]
global_metric = (global_jacobian.T*ETA*global_jacobian).applyfunc(sp.expand)
static_maps_ok = all(
    matrix_equal(forms.subs(R_PLUS, R_MINUS), global_forms)
    for forms in affine_maps.values()
)
static_metrics_ok = all(
    matrix_equal(metric.subs(R_PLUS, R_MINUS), global_metric)
    for metric in metrics.values()
)
static_control_ok = bool(static_maps_ok and static_metrics_ok)
check("the static limit is one global affine map and one metric",
      static_control_ok)


# Chamber-by-chamber exact equality classes.
chamber_records = []
map_count_distribution = Counter()
metric_count_distribution = Counter()
pair_agreement = [[0 for _ in ORDERS] for _ in ORDERS]
mismatch_witness = None
for word, row in zip(words, assignment_table):
    selected_keys = [(order, split) for order, split in zip(ORDERS, row)]
    selected_maps = [affine_maps[key] for key in selected_keys]
    selected_metrics = [metrics[key] for key in selected_keys]
    map_keys = [matrix_key(item) for item in selected_maps]
    metric_keys = [matrix_key(item) for item in selected_metrics]
    map_count = len(set(map_keys))
    metric_count = len(set(metric_keys))
    map_count_distribution[map_count] += 1
    metric_count_distribution[metric_count] += 1
    chamber_records.append({
        "sign_word": word,
        "distinct_affine_maps": map_count,
        "distinct_pullback_metrics": metric_count,
        "splits_by_lexicographic_order": list(row),
    })
    for left_index in range(len(ORDERS)):
        for right_index in range(len(ORDERS)):
            if metric_keys[left_index] == metric_keys[right_index]:
                pair_agreement[left_index][right_index] += 1
    if mismatch_witness is None and metric_count > 1:
        for left_index in range(len(ORDERS)):
            for right_index in range(left_index+1, len(ORDERS)):
                if metric_keys[left_index] == metric_keys[right_index]:
                    continue
                difference = selected_metrics[left_index]-selected_metrics[right_index]
                first = next(
                    (index for index, entry in enumerate(difference)
                     if sp.expand(entry) != 0),
                    None,
                )
                mismatch_witness = {
                    "sign_word": word,
                    "left_order": list(ORDERS[left_index]),
                    "left_split": row[left_index],
                    "right_order": list(ORDERS[right_index]),
                    "right_split": row[right_index],
                    "first_differing_entry": [first//4, first % 4],
                    "first_difference_factored": str(
                        sp.factor(difference[first], extension=SQRT5)
                    ),
                    "difference_matrix_factored": [
                        [str(sp.factor(difference[i, j], extension=SQRT5))
                         for j in range(4)]
                        for i in range(4)
                    ],
                }
                break
            if mismatch_witness is not None:
                break

identically_compatible_chambers = metric_count_distribution[1]
chamber_census_ok = bool(
    len(chamber_records) == 148
    and sum(map_count_distribution.values()) == 148
    and sum(metric_count_distribution.values()) == 148
    and all(pair_agreement[index][index] == 148 for index in range(24))
    and all(pair_agreement[i][j] == pair_agreement[j][i]
            for i in range(24) for j in range(24))
)
check("the complete affine-map and intrinsic-metric census is closed",
      chamber_census_ok)


# Exact S4 x C2 chamber action and orbit-invariance of equality-class counts.
def permute_mask(mask, permutation):
    return mask_of(
        permutation[vertex]
        for vertex in VERTICES
        if mask & (1 << vertex)
    )


def transform_pattern(pattern, permutation, reflect_time):
    result = [None]*len(MASKS)
    for mask, sign in zip(MASKS, pattern):
        target = permute_mask(mask, permutation)
        if reflect_time:
            target = FULL_MASK ^ target
            sign = -sign
        result[MASK_INDEX[target]] = sign
    return tuple(result)


pattern_to_index = {pattern: index for index, pattern in enumerate(patterns)}
transforms = tuple(
    (permutation, reflect_time)
    for permutation in ORDERS
    for reflect_time in (False, True)
)
symmetry_ok = True
unseen = set(patterns)
computed_orbits = []
while unseen:
    representative = min(unseen)
    orbit = {
        transform_pattern(representative, permutation, reflect_time)
        for permutation, reflect_time in transforms
    }
    symmetry_ok &= orbit <= set(patterns)
    indices = [pattern_to_index[item] for item in orbit]
    map_counts = {chamber_records[index]["distinct_affine_maps"] for index in indices}
    metric_counts = {
        chamber_records[index]["distinct_pullback_metrics"] for index in indices
    }
    symmetry_ok &= len(map_counts) == len(metric_counts) == 1
    computed_orbits.append({
        "representative": word_of(representative),
        "size": len(orbit),
    })
    unseen -= orbit

recorded_orbits = sorted(
    overlay["symmetry_orbits"], key=lambda item: item["representative"]
)
computed_orbits = sorted(
    computed_orbits, key=lambda item: item["representative"]
)
symmetry_ok &= computed_orbits == recorded_orbits
check("the metric-count census is invariant under the full S4 x C2 action",
      symmetry_ok)


# Complete global compatibility polynomial after R_minus=1, R_plus=r.
polynomial_expressions = {}
polynomial_static_ok = True
for row in assignment_table:
    reference = metrics[(ORDERS[0], row[0])]
    for order, split in zip(ORDERS[1:], row[1:]):
        difference = metrics[(order, split)]-reference
        for entry in difference:
            expression = sp.expand(entry.subs({R_MINUS: 1, R_PLUS: ratio}))
            polynomial_static_ok &= bool(
                not expression.has(T)
                and sp.expand(expression.subs(ratio, 1)) == 0
            )
            if expression != 0:
                polynomial_expressions.setdefault(expression_key(expression), expression)

polynomials = [
    sp.Poly(expression, ratio, extension=SQRT5)
    for expression in polynomial_expressions.values()
]

gcd_factor_control_ok = True
gcd_poly = None
radical_gcd = None
common_factor_polys = {}
compatible_roots = []
positive_roots = []
unknown_sign_roots = []
direct_root_checks = {}

if polynomials:
    gcd_poly = reduce(sp.gcd, polynomials).monic()
    radical_gcd = gcd_poly.sqf_part().monic()
    common_factor_keys = None
    factor_catalog = {}
    for polynomial in polynomials:
        factor_keys = set()
        _, factors_list = sp.factor_list(
            polynomial.as_expr(), ratio, extension=SQRT5
        )
        for factor, _multiplicity in factors_list:
            factor_poly = sp.Poly(factor, ratio, extension=SQRT5).monic()
            key = expression_key(factor_poly.as_expr())
            factor_keys.add(key)
            factor_catalog[key] = factor_poly
        common_factor_keys = (
            factor_keys if common_factor_keys is None
            else common_factor_keys & factor_keys
        )
    common_factor_polys = {
        key: factor_catalog[key] for key in sorted(common_factor_keys or ())
    }
    independent_radical = sp.Poly(1, ratio, extension=SQRT5)
    for factor_poly in common_factor_polys.values():
        independent_radical *= factor_poly
    independent_radical = independent_radical.monic()
    gcd_factor_control_ok = independent_radical == radical_gcd

    degree = radical_gcd.degree()
    if degree == 1:
        compatible_roots = [sp.simplify(-radical_gcd.nth(0)/radical_gcd.nth(1))]
    elif degree == 2:
        a = radical_gcd.nth(2)
        b = radical_gcd.nth(1)
        c = radical_gcd.nth(0)
        discriminant = sp.expand(b*b-4*a*c)
        compatible_roots = [
            sp.simplify((-b-sp.sqrt(discriminant))/(2*a)),
            sp.simplify((-b+sp.sqrt(discriminant))/(2*a)),
        ]
    elif degree > 2:
        compatible_roots = [
            sp.simplify(root)
            for root in sp.solve(radical_gcd.as_expr(), ratio)
        ]
        gcd_factor_control_ok &= len(compatible_roots) == degree

    unique_roots = {}
    for root in compatible_roots:
        unique_roots.setdefault(expression_key(root), root)
    compatible_roots = list(unique_roots.values())
    for root in compatible_roots:
        real_status = sp.simplify(sp.im(root)) == 0 or root.is_real is True
        positive_status = root.is_positive
        if real_status and positive_status is True:
            positive_roots.append(root)
        elif real_status and positive_status is None:
            unknown_sign_roots.append(root)

    for root in positive_roots:
        all_equal = True
        for row in assignment_table:
            reference = metrics[(ORDERS[0], row[0])].subs(
                {R_MINUS: 1, R_PLUS: root}
            )
            for order, split in zip(ORDERS[1:], row[1:]):
                candidate = metrics[(order, split)].subs(
                    {R_MINUS: 1, R_PLUS: root}
                )
                all_equal &= matrix_equal(reference, candidate)
        direct_root_checks[str(root)] = bool(all_equal)
else:
    # Frozen convention: the zero polynomial means compatibility at every
    # positive ratio, not a failed gcd calculation.
    direct_root_checks["ALL_POSITIVE_RATIOS"] = bool(
        identically_compatible_chambers == len(words)
    )

root_control_ok = bool(
    polynomial_static_ok
    and gcd_factor_control_ok
    and not unknown_sign_roots
    and all(direct_root_checks.values())
    and (
        (not polynomials and direct_root_checks["ALL_POSITIVE_RATIOS"])
        or (polynomials and any(sp.simplify(root-1) == 0 for root in positive_roots))
    )
)
check("gcd, independent factor intersection and direct root substitutions agree",
      root_control_ok)


controls_ok = all(condition for _, condition in tests)
if not controls_ok:
    outcome = "OVERLAY_METRIC_COMPATIBILITY_CONTROL_FAILED"
elif not polynomials:
    outcome = "OVERLAY_INHERITS_DYNAMIC_REGGE_METRIC"
else:
    dynamic_roots = [
        root for root in positive_roots if sp.simplify(root-1) != 0
    ]
    if dynamic_roots:
        outcome = "OVERLAY_INHERITS_DYNAMIC_REGGE_METRIC"
    elif identically_compatible_chambers:
        outcome = "OVERLAY_METRIC_COMPATIBILITY_PARTIAL_ONLY"
    else:
        outcome = "OVERLAY_INHERITS_STATIC_METRIC_ONLY"


def polynomial_text(polynomial):
    return None if polynomial is None else str(sp.factor(
        polynomial.as_expr(), extension=SQRT5
    ))


payload = {
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "protocol_clarification_commit": PROTOCOL_CLARIFICATION_COMMIT,
    "input_sha256": {
        "universal_overlay": digest(OVERLAY_SOURCE),
        "overlay_face_poset": digest(FACE_SOURCE),
    },
    "symbols": {
        "spatial_scales": [str(R_MINUS), str(R_PLUS)],
        "outer_time": str(T),
        "normalized_ratio": str(ratio),
        "regular_tetrahedron_off_diagonal_gram": str(ADJACENT_DOT),
    },
    "labelled_staircase_simplex_count": len(affine_maps),
    "distinct_affine_map_types_global": len({
        matrix_key(value) for value in affine_maps.values()
    }),
    "distinct_pullback_metric_types_global": len(metric_types),
    "simplex_vertex_matrix_determinant_distribution": {
        str(key): value for key, value in sorted(
            Counter(str(value) for value in determinants.values()).items()
        )
    },
    "continuity_failure_count": len(continuity_failures),
    "continuity_failures": continuity_failures,
    "time_cancellation_failure_count": len(t_cancellation_failures),
    "chamber_count": len(chamber_records),
    "chambers": chamber_records,
    "distinct_affine_map_count_distribution": {
        str(key): value for key, value in sorted(map_count_distribution.items())
    },
    "distinct_pullback_metric_count_distribution": {
        str(key): value for key, value in sorted(metric_count_distribution.items())
    },
    "identically_metric_compatible_chambers": identically_compatible_chambers,
    "schedule_pair_metric_agreement_counts": pair_agreement,
    "first_metric_mismatch": mismatch_witness,
    "symmetry_orbits": computed_orbits,
    "unique_nonzero_compatibility_polynomial_count": len(polynomials),
    "global_compatibility_gcd": polynomial_text(gcd_poly),
    "global_compatibility_radical_gcd": polynomial_text(radical_gcd),
    "independent_common_irreducible_factors": [
        str(polynomial.as_expr()) for polynomial in common_factor_polys.values()
    ],
    "all_exact_common_roots": [str(root) for root in compatible_roots],
    "positive_compatible_ratios": (
        "ALL_POSITIVE_RATIOS" if not polynomials
        else [str(root) for root in positive_roots]
    ),
    "unknown_sign_roots": [str(root) for root in unknown_sign_roots],
    "direct_root_substitution_checks": direct_root_checks,
    "static_global_affine_control": static_maps_ok,
    "static_global_metric_control": static_metrics_ok,
    "regge_action_evaluations": 0,
    "dust_action_evaluations": 0,
    "continuum_target_parsed": False,
    "tests": len(tests),
    "passed": sum(condition for _, condition in tests),
    "outcome": outcome,
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")

for label, condition in tests:
    print(f"{'PASS' if condition else 'FAIL'}: {label}")
print(f"global affine-map types={payload['distinct_affine_map_types_global']}")
print(f"global pullback-metric types={payload['distinct_pullback_metric_types_global']}")
print(f"per-chamber affine-map distribution={dict(sorted(map_count_distribution.items()))}")
print(f"per-chamber metric distribution={dict(sorted(metric_count_distribution.items()))}")
print(f"identically compatible chambers={identically_compatible_chambers}/148")
print(f"compatibility gcd={payload['global_compatibility_gcd']}")
print(f"positive compatible ratios={payload['positive_compatible_ratios']}")
print(f"OUTCOME: {outcome}")
print(f"{payload['passed']}/{payload['tests']} tests passed")

raise SystemExit(0 if controls_ok else 1)
