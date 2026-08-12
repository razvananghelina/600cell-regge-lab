#!/usr/bin/env python3
"""Symbolic certificate for all-q locality of the rank-edgewise tower.

Protocol commit 8d0c557 froze the theorem, hypotheses, sharp bounds and
falsification gates before q=3 or the finite partition certificate was run.

The infinite step uses the fixed-dimension link classification of Jojic and
Papaz (arXiv:2408.12756v3).  This verifier checks all finite consequences of
that theorem, the relative gluing to sd(boundary Delta^4), and the q<=4
saturation/control data.  It does not infer an infinite claim from a fit.
"""

from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations, permutations
from math import factorial
import json
from pathlib import Path

import sympy as sy

from whitney_trace_refinement_tools import (
    barycentric_refine,
    edgewise_local_facets,
    make_base_level,
)


OUTPUT = Path(__file__).with_name("whitney_dual_resolution_all_k.json")
PROTOCOL_COMMIT = "8d0c557"
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    condition = bool(condition)
    passed += int(condition)
    print(f"  [{'PASS' if condition else 'FAIL'}] {label}")
    if detail:
        print(f"         {detail}")


def integer_partitions(total, maximum=None):
    """Unordered positive partitions in nonincreasing order."""
    if total == 0:
        yield ()
        return
    if maximum is None or maximum > total:
        maximum = total
    for first in range(maximum, 0, -1):
        for rest in integer_partitions(total - first, first):
            yield (first,) + rest


def multinomial(parts):
    total = sum(parts)
    value = factorial(total)
    for part in parts:
        value //= factorial(part)
    return value


def cyclic_gap_partition(ranks):
    ranks = tuple(sorted(ranks))
    gaps = [right - left for left, right in zip(ranks, ranks[1:])]
    gaps.append(4 - ranks[-1] + ranks[0])
    return tuple(sorted(gaps, reverse=True))


def base_chamber_multiplicity(ranks):
    """Number of S_5 flags containing a fixed chain of these ranks."""
    ranks = tuple(sorted(ranks))
    gaps = [ranks[0]]
    gaps.extend(right - left for left, right in zip(ranks, ranks[1:]))
    gaps.append(5 - ranks[-1])
    value = 1
    for gap in gaps:
        value *= factorial(gap)
    return value


def enumerated_base_chamber_multiplicity(ranks):
    """Independent S_5 flag count for one canonical nested face chain."""
    targets = {
        rank: frozenset(range(rank)) for rank in ranks
    }
    return sum(
        all(frozenset(ordering[:rank]) == targets[rank] for rank in ranks)
        for ordering in permutations(range(5))
    )


def expected_f_vector(q):
    f0 = 20 * q**3 + 10 * q
    return (f0, 140 * q**3 + 10 * q, 240 * q**3, 120 * q**3)


def expected_vertex_histogram(q):
    values = {
        12: 20,
        16: 30 * (q - 1),
        24: 20 * q**3 - 20 * q + 10,
    }
    return {key: value for key, value in values.items() if value}


def expected_edge_histogram(q):
    return {
        4: 60 * q**3 + 30 * q,
        6: 80 * q**3 - 20 * q,
    }


def expected_degree_histogram(q):
    vertex = expected_vertex_histogram(q)
    degree_by_tetrahedra = {12: 8, 16: 10, 24: 14}
    return {
        degree_by_tetrahedra[key]: value for key, value in vertex.items()
    }


def keyed_rank_edgewise(q):
    """Build K_q while retaining exact base carriers and parent chambers."""
    regular_vertices = tuple(map(sy.Matrix, (
        (1, 1, 1),
        (1, -1, -1),
        (-1, 1, -1),
        (-1, -1, 1),
    )))
    coarse = make_base_level(regular_vertices)
    ranked = barycentric_refine(coarse)
    coarse_vertex_cells = tuple(
        cell for layer in coarse["cells"] for cell in layer
    )
    coarse_vertex_ranks = tuple(map(len, coarse_vertex_cells))

    local_facets = edgewise_local_facets(q)
    key_to_vertex = {}
    vertex_keys = []
    top_records = []
    for parent, top in enumerate(ranked["top"]):
        for facet in local_facets:
            child = []
            for numerator in facet:
                key = tuple(
                    (top[rank], Fraction(weight, q))
                    for rank, weight in enumerate(numerator)
                    if weight
                )
                if key not in key_to_vertex:
                    key_to_vertex[key] = len(vertex_keys)
                    vertex_keys.append(key)
                child.append(key_to_vertex[key])
            top_records.append((tuple(sorted(child)), parent))

    top = tuple(record[0] for record in top_records)
    edges = tuple(sorted({
        tuple(sorted(edge))
        for tetrahedron in top
        for edge in combinations(tetrahedron, 2)
    }))
    triangles = tuple(sorted({
        tuple(sorted(triangle))
        for tetrahedron in top
        for triangle in combinations(tetrahedron, 3)
    }))
    edge_index = {edge: index for index, edge in enumerate(edges)}

    vertex_parent_counts = [Counter() for _ in vertex_keys]
    edge_parent_counts = [Counter() for _ in edges]
    neighbours = [set() for _ in vertex_keys]
    for tetrahedron, parent in top_records:
        for vertex in tetrahedron:
            vertex_parent_counts[vertex][parent] += 1
        for raw_edge in combinations(tetrahedron, 2):
            edge = tuple(sorted(raw_edge))
            index = edge_index[edge]
            edge_parent_counts[index][parent] += 1
            left, right = edge
            neighbours[left].add(right)
            neighbours[right].add(left)

    parent_sets = []
    for support in (
        {base_vertex for base_vertex, _ in key} for key in vertex_keys
    ):
        parent_sets.append({
            parent for parent, tetrahedron in enumerate(ranked["top"])
            if support.issubset(tetrahedron)
        })

    vertex_records = []
    for vertex, key in enumerate(vertex_keys):
        support = {base_vertex for base_vertex, _ in key}
        ranks = tuple(sorted(coarse_vertex_ranks[item] for item in support))
        counts = vertex_parent_counts[vertex]
        vertex_records.append({
            "ranks": ranks,
            "parent_set_exact": set(counts) == parent_sets[vertex],
            "parent_count": len(counts),
            "local_counts": tuple(sorted(counts.values())),
            "global_tetrahedra": sum(counts.values()),
            "degree": len(neighbours[vertex]),
        })

    edge_records = []
    for edge, counts in zip(edges, edge_parent_counts):
        support = {
            base_vertex
            for vertex in edge
            for base_vertex, _ in vertex_keys[vertex]
        }
        ranks = tuple(sorted(coarse_vertex_ranks[item] for item in support))
        expected_parents = {
            parent for parent, tetrahedron in enumerate(ranked["top"])
            if support.issubset(tetrahedron)
        }
        edge_records.append({
            "ranks": ranks,
            "parent_set_exact": set(counts) == expected_parents,
            "parent_count": len(counts),
            "local_counts": tuple(sorted(counts.values())),
            "global_tetrahedra": sum(counts.values()),
        })

    return {
        "f_vector": (
            len(vertex_keys), len(edges), len(triangles), len(top_records)
        ),
        "vertex_records": vertex_records,
        "edge_records": edge_records,
    }


def local_saturation(q=4):
    facets = edgewise_local_facets(q)
    vertex_counts = Counter(vertex for top in facets for vertex in top)
    edge_tops = defaultdict(list)
    for top in facets:
        for edge in combinations(top, 2):
            edge_tops[tuple(sorted(edge))].append(top)
    vertex_types = set()
    for vertex, count in vertex_counts.items():
        support = tuple(index + 1 for index, value in enumerate(vertex) if value)
        vertex_types.add((cyclic_gap_partition(support), count))
    edge_types = set()
    edge_links_have_expected_topology = True
    for edge, incident_tops in edge_tops.items():
        count = len(incident_tops)
        carrier_size = sum(
            bool(edge[0][index] or edge[1][index]) for index in range(4)
        )
        edge_types.add((carrier_size, count))
        link_edges = [
            tuple(vertex for vertex in top if vertex not in edge)
            for top in incident_tops
        ]
        link_degrees = Counter(vertex for item in link_edges for vertex in item)
        adjacency = defaultdict(set)
        for left, right in link_edges:
            adjacency[left].add(right)
            adjacency[right].add(left)
        reached = set()
        if adjacency:
            frontier = [next(iter(adjacency))]
            while frontier:
                vertex = frontier.pop()
                if vertex in reached:
                    continue
                reached.add(vertex)
                frontier.extend(adjacency[vertex] - reached)
        connected = reached == set(adjacency)
        if carrier_size < 4:
            expected = (
                connected
                and sorted(link_degrees.values()).count(1) == 2
                and all(value in (1, 2) for value in link_degrees.values())
            )
        else:
            expected = connected and all(
                value == 2 for value in link_degrees.values()
            )
        edge_links_have_expected_topology &= expected
    return vertex_types, edge_types, edge_links_have_expected_topology


print("=" * 78)
print("ALL-RESOLUTION EDGEWISE DUAL-LOCALITY CERTIFICATE")
print("=" * 78)


# Theorem 3.3/3.5: vertex links K_lambda for every partition of four.
theoretical_vertex_types = {
    (partition, multinomial(partition))
    for partition in integer_partitions(4)
}
expected_vertex_types = {
    ((4,), 1),
    ((3, 1), 4),
    ((2, 2), 6),
    ((2, 1, 1), 12),
    ((1, 1, 1, 1), 24),
}
check(
    "the partition classification gives exactly five local vertex-link types",
    theoretical_vertex_types == expected_vertex_types,
    str(sorted(theoretical_vertex_types)),
)


# Theorem 4.2/Corollary 4.3: an edge link is a join of two K_sigma.
theoretical_edge_types = set()
edge_partition_records = []
for lam in (partition for partition in integer_partitions(4) if len(partition) == 2):
    for sigma_left in integer_partitions(lam[0]):
        for sigma_right in integer_partitions(lam[1]):
            carrier_size = len(sigma_left) + len(sigma_right)
            local_facets = multinomial(sigma_left) * multinomial(sigma_right)
            theoretical_edge_types.add((carrier_size, local_facets))
            edge_partition_records.append({
                "lambda": lam,
                "sigma_left": sigma_left,
                "sigma_right": sigma_right,
                "carrier_size": carrier_size,
                "local_facets": local_facets,
            })
expected_edge_types = {(2, 1), (3, 2), (3, 3), (4, 4), (4, 6)}
check(
    "the face-link partition theorem leaves exactly five relative edge types",
    theoretical_edge_types == expected_edge_types,
    str(sorted(theoretical_edge_types)),
)


# Every nonempty rank face of one chamber, with exact S_5 coface count.
rank_records = []
analytic_global_vertex_values = set()
for size in range(1, 5):
    for ranks in combinations(range(1, 5), size):
        partition = cyclic_gap_partition(ranks)
        local = multinomial(partition)
        base = base_chamber_multiplicity(ranks)
        global_value = base * local
        analytic_global_vertex_values.add(global_value)
        rank_records.append({
            "ranks": ranks,
            "base_chambers": base,
            "cyclic_gap_partition": partition,
            "local_tetrahedra": local,
            "global_tetrahedra": global_value,
        })

check(
    "all fifteen base rank chains obey the flag/cyclic-gap formula",
    len(rank_records) == 15
    and all(
        record["base_chambers"]
        == enumerated_base_chamber_multiplicity(record["ranks"])
        for record in rank_records
    ),
)
check(
    "the analytic all-q vertex incidences are exactly 12, 16 and 24",
    analytic_global_vertex_values == {12, 16, 24},
    str(sorted(analytic_global_vertex_values)),
)


local_by_carrier_size = defaultdict(set)
for carrier_size, local_count in theoretical_edge_types:
    local_by_carrier_size[carrier_size].add(local_count)
analytic_global_edge_values = set()
for size in (2, 3, 4):
    for ranks in combinations(range(1, 5), size):
        base = base_chamber_multiplicity(ranks)
        for local in local_by_carrier_size[size]:
            analytic_global_edge_values.add(base * local)
check(
    "relative gluing turns every legal local edge type into C4 or C6",
    analytic_global_edge_values == {4, 6},
    str(sorted(analytic_global_edge_values)),
)


# Symbolic all-q census.  A base face of rank-set S occurs 120/m_B(S)
# times, and its relative interior contains binomial(q-1, |S|-1) vertices.
q_symbol = sy.symbols("q", integer=True, positive=True)
symbolic_f0 = sy.Integer(0)
symbolic_vertex_histogram = defaultdict(lambda: sy.Integer(0))
for record in rank_records:
    base_face_count = 120 // record["base_chambers"]
    refined_vertices = base_face_count * sy.binomial(
        q_symbol - 1, len(record["ranks"]) - 1
    )
    symbolic_f0 += refined_vertices
    symbolic_vertex_histogram[record["global_tetrahedra"]] += refined_vertices
symbolic_f = (
    sy.expand_func(symbolic_f0).expand(),
    None,
    240 * q_symbol**3,
    120 * q_symbol**3,
)
symbolic_f = (
    symbolic_f[0],
    sy.expand(symbolic_f[0] + symbolic_f[2] - symbolic_f[3]),
    symbolic_f[2],
    symbolic_f[3],
)
symbolic_vertex_histogram = {
    key: sy.expand_func(value).expand()
    for key, value in symbolic_vertex_histogram.items()
}
check(
    "the all-q f-vector and vertex histogram follow symbolically",
    symbolic_f == (
        20 * q_symbol**3 + 10 * q_symbol,
        140 * q_symbol**3 + 10 * q_symbol,
        240 * q_symbol**3,
        120 * q_symbol**3,
    )
    and symbolic_vertex_histogram == {
        12: sy.Integer(20),
        16: 30 * q_symbol - 30,
        24: 20 * q_symbol**3 - 20 * q_symbol + 10,
    },
)

# Since every edge has degree four or six, E4+E6=f1 and
# 4E4+6E6=6f3 determine the complete histogram without a fit.
edge_four = sy.expand(3 * symbolic_f[1] - 3 * symbolic_f[3])
edge_six = sy.expand(3 * symbolic_f[3] - 2 * symbolic_f[1])
check(
    "double incidence derives the all-q C4/C6 edge histogram",
    edge_four == 60 * q_symbol**3 + 30 * q_symbol
    and edge_six == 80 * q_symbol**3 - 20 * q_symbol,
)


(
    observed_vertex_types,
    observed_edge_types,
    local_edge_topology_exact,
) = local_saturation(4)
check(
    "q=4 realizes every partition-classified local vertex-link type",
    observed_vertex_types == theoretical_vertex_types,
    f"observed={sorted(observed_vertex_types)}",
)
check(
    "q=4 realizes every relative edge-link isomorphism type",
    observed_edge_types == theoretical_edge_types,
    f"observed={sorted(observed_edge_types)}",
)
check(
    "each saturated local edge link is the predicted path or cycle",
    local_edge_topology_exact,
)


control_records = []
all_f_vectors = True
all_histograms = True
all_relative_gluing = True
for q in (1, 2, 3, 4):
    audit = keyed_rank_edgewise(q)
    vertex_hist = Counter(
        record["global_tetrahedra"] for record in audit["vertex_records"]
    )
    edge_hist = Counter(
        record["global_tetrahedra"] for record in audit["edge_records"]
    )
    degree_hist = Counter(record["degree"] for record in audit["vertex_records"])
    all_f_vectors &= audit["f_vector"] == expected_f_vector(q)
    all_histograms &= (
        dict(sorted(vertex_hist.items())) == expected_vertex_histogram(q)
        and dict(sorted(edge_hist.items())) == expected_edge_histogram(q)
        and dict(sorted(degree_hist.items())) == expected_degree_histogram(q)
    )

    vertex_gluing = all(
        record["parent_set_exact"]
        and record["parent_count"] == base_chamber_multiplicity(record["ranks"])
        and len(set(record["local_counts"])) == 1
        and record["local_counts"][0]
            == multinomial(cyclic_gap_partition(record["ranks"]))
        for record in audit["vertex_records"]
    )
    edge_gluing = all(
        record["parent_set_exact"]
        and record["parent_count"] == base_chamber_multiplicity(record["ranks"])
        and len(set(record["local_counts"])) == 1
        and record["local_counts"][0]
            in local_by_carrier_size[len(record["ranks"])]
        for record in audit["edge_records"]
    )
    all_relative_gluing &= vertex_gluing and edge_gluing
    control_records.append({
        "q": q,
        "f_vector": audit["f_vector"],
        "vertex_tetrahedron_histogram": dict(sorted(vertex_hist.items())),
        "edge_tetrahedron_histogram": dict(sorted(edge_hist.items())),
        "vertex_degree_histogram": dict(sorted(degree_hist.items())),
        "vertex_relative_gluing_exact": vertex_gluing,
        "edge_relative_gluing_exact": edge_gluing,
    })
    print(
        f"q={q}: f={audit['f_vector']} "
        f"vertex={dict(sorted(vertex_hist.items()))} "
        f"edge={dict(sorted(edge_hist.items()))}"
    )

check("the exact all-q f-vector formula holds on q=1,2,3,4", all_f_vectors)
check(
    "the closed-form vertex, edge and degree histograms hold on all controls",
    all_histograms,
)
check(
    "parentwise relative gluing agrees with the analytic certificate",
    all_relative_gluing,
)
check(
    "the preregistered previously untested q=3 control is exact",
    control_records[2]["f_vector"] == (570, 3810, 6480, 3240)
    and control_records[2]["vertex_tetrahedron_histogram"]
        == {12: 20, 16: 60, 24: 490}
    and control_records[2]["edge_tetrahedron_histogram"]
        == {4: 1710, 6: 2100}
    and control_records[2]["vertex_degree_histogram"]
        == {8: 20, 10: 60, 14: 490},
)


# The degree formula is a symbolic consequence of a closed triangulated
# 2-sphere link: 3F=2E and V-E+F=2 imply V=F/2+2.
sphere_degree_map = {facets: facets // 2 + 2 for facets in (12, 16, 24)}
check(
    "sphere Euler identities give vertex degrees 8, 10 and 14",
    sphere_degree_map == {12: 8, 16: 10, 24: 14},
)


certificate = {
    "protocol_commit": PROTOCOL_COMMIT,
    "external_hypothesis": {
        "source": "Jojic--Papaz, arXiv:2408.12756v3",
        "fixed_simplex_parameter": 4,
        "claim_used": "vertex and face links are classified by partitions",
    },
    "theoretical_local_vertex_types": [
        {"partition": list(partition), "facets": facets}
        for partition, facets in sorted(theoretical_vertex_types)
    ],
    "theoretical_relative_edge_types": [
        {"carrier_size": size, "local_facets": facets}
        for size, facets in sorted(theoretical_edge_types)
    ],
    "edge_partition_records": [
        {
            key: list(value) if isinstance(value, tuple) else value
            for key, value in record.items()
        }
        for record in edge_partition_records
    ],
    "rank_chain_records": [
        {
            key: list(value) if isinstance(value, tuple) else value
            for key, value in record.items()
        }
        for record in rank_records
    ],
    "all_q_formulas": {
        "f_vector": ["20*q^3+10*q", "140*q^3+10*q", "240*q^3", "120*q^3"],
        "vertex_tetrahedron_histogram": {
            "12": "20",
            "16": "30*(q-1)",
            "24": "20*q^3-20*q+10",
        },
        "edge_tetrahedron_histogram": {
            "4": "60*q^3+30*q",
            "6": "80*q^3-20*q",
        },
        "sharp_maxima": {"a0": 24, "a1": 6, "r3": 14},
    },
    "controls": control_records,
    "proof_scope": (
        "all positive q, conditional on the cited fixed-dimension partition "
        "classification; no finite extrapolation"
    ),
    "physical_scope": "kinematic locality only",
}
OUTPUT.write_text(json.dumps(certificate, indent=2) + "\n")


print("-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
if passed == tests:
    print("DERIVED UNIFORM LOCALITY under the stated partition-classification hypotheses")
    print("a0(q)=24, a1(q)=6, r3(q)=14 for every positive integer q")
else:
    print("ALL-q THEOREM NOT ESTABLISHED")

raise SystemExit(0 if passed == tests else 1)
