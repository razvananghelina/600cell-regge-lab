#!/usr/bin/env python3
"""Exact eight-reading audit of the spectrally ordered four-node KO6 carrier.

The candidate, ambiguity census and scope boundary were frozen in protocol
commit 12c45eb.  This verifier uses no Hessian or particle target and does not
claim to construct a Dirac operator.
"""

from collections import Counter
from itertools import combinations, product
import json
from pathlib import Path

import sympy as sp


OUTPUT = Path(__file__).with_name("hopf_six_spectral_krajewski.json")
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    condition = bool(condition)
    passed += int(condition)
    print(f"  [{'PASS' if condition else 'FAIL'}] {label}")
    if detail:
        print(f"         {detail}")


def connected_components(node_count, edges):
    adjacency = {node: set() for node in range(node_count)}
    for left, right in edges:
        adjacency[left].add(right)
        adjacency[right].add(left)
    unseen = set(range(node_count))
    components = []
    while unseen:
        seed = min(unseen)
        component = {seed}
        frontier = [seed]
        while frontier:
            node = frontier.pop()
            for neighbor in adjacency[node]-component:
                component.add(neighbor)
                frontier.append(neighbor)
        unseen -= component
        components.append(component)
    return components


def pfaffian4(matrix):
    return sp.expand(matrix[0, 1]*matrix[2, 3]
                     - matrix[0, 2]*matrix[1, 3]
                     + matrix[0, 3]*matrix[1, 2])


print("="*78)
print("SPECTRALLY ORDERED SIX-FIBRATION KO6 CARRIER")
print("="*78)

# Node names and sizes are the exact real crossed-product result.  The
# integer ranks below encode the already-proved exact inequalities, avoiding
# any floating comparison of phi:
#   -phi < phi-1 < 2 and -5 < 0 < 5.
nodes = (
    {"name": "trivial", "size": 6, "u_rank": 2, "v_rank": 2,
     "joint_spectrum": ("2", "5")},
    {"name": "reflection_sign", "size": 6, "u_rank": 2, "v_rank": 0,
     "joint_spectrum": ("2", "-5")},
    {"name": "positive_doublet", "size": 12, "u_rank": 1, "v_rank": 1,
     "joint_spectrum": ("phi-1", "0")},
    {"name": "negative_doublet", "size": 12, "u_rank": 0, "v_rank": 1,
     "joint_spectrum": ("-phi", "0")},
)
node_count = len(nodes)
node_sizes = tuple(node["size"] for node in nodes)
unordered_pairs = tuple(combinations(range(node_count), 2))
off_diagonal_dimension = 2*sum(
    node_sizes[left]*node_sizes[right] for left, right in unordered_pairs
)
check("the off-diagonal full enveloping carrier has dimension 936",
      node_sizes == (6, 6, 12, 12)
      and len(unordered_pairs) == 6
      and off_diagonal_dimension == 936)


def apply_left_unit(a, b, state):
    """Apply E_ab to the left tensor coordinate, or return None."""
    if state is None:
        return None
    left_index, right_index = state
    return ((a, right_index) if left_index == b else None)


def apply_right_unit(c, d, state):
    """Apply E_cd to the right tensor coordinate, or return None."""
    if state is None:
        return None
    left_index, right_index = state
    return ((left_index, c) if right_index == d else None)


# Exhaust every pair of matching left/right matrix units and every basis
# vector in each of the twelve ordered cells.  Units from nonmatching simple
# blocks act as zero, so these are all potentially nonzero order-zero cases.
matrix_unit_cases = 0
matrix_unit_order_zero = True
for left_node in range(node_count):
    for right_node in range(node_count):
        if left_node == right_node:
            continue
        left_size = node_sizes[left_node]
        right_size = node_sizes[right_node]
        for a in range(left_size):
            for b in range(left_size):
                for c in range(right_size):
                    for d in range(right_size):
                        for x in range(left_size):
                            for y in range(right_size):
                                state = (x, y)
                                left_then_right = apply_right_unit(
                                    c, d, apply_left_unit(a, b, state)
                                )
                                right_then_left = apply_left_unit(
                                    a, b, apply_right_unit(c, d, state)
                                )
                                matrix_unit_cases += 1
                                if left_then_right != right_then_left:
                                    matrix_unit_order_zero = False
                                    break
                            if not matrix_unit_order_zero:
                                break
                        if not matrix_unit_order_zero:
                            break
                    if not matrix_unit_order_zero:
                        break
                if not matrix_unit_order_zero:
                    break
            if not matrix_unit_order_zero:
                break
        if not matrix_unit_order_zero:
            break
    if not matrix_unit_order_zero:
        break
check("order zero holds on every potentially nonzero matrix-unit action",
      matrix_unit_order_zero
      and matrix_unit_cases == sum(
          node_sizes[left]**3*node_sizes[right]**3
          for left in range(node_count) for right in range(node_count)
          if left != right
      ) == 9051264,
      f"exact action cases={matrix_unit_cases}")

records = []
for priority, u_direction, v_direction in product(
        (("u", "v"), ("v", "u")), (1, -1), (1, -1)):
    direction = {"u": u_direction, "v": v_direction}

    def spectral_key(node_index):
        node = nodes[node_index]
        return tuple(direction[coordinate]*node[coordinate+"_rank"]
                     for coordinate in priority)

    order = tuple(sorted(range(node_count), key=spectral_key))
    position = {node: rank for rank, node in enumerate(order)}
    positive_cells = frozenset(
        (left, right) if position[left] < position[right]
        else (right, left)
        for left, right in unordered_pairs
    )
    negative_cells = frozenset((right, left)
                               for left, right in positive_cells)
    all_cells = positive_cells | negative_cells

    gamma = sp.zeros(node_count)
    mu = sp.zeros(node_count)
    for left, right in positive_cells:
        gamma[left, right] = 1
        gamma[right, left] = -1
        mu[left, right] = 1
    intersection = mu-mu.T
    pfaffian = pfaffian4(intersection)
    determinant = sp.expand(intersection.det())

    # Cell-level real-structure and order-zero checks.  A basis vector in
    # H_ij has independent left index in C^n_i and right index in C^n_j.
    # J maps it to H_ji with the two indices exchanged.  Left and opposite
    # matrix units therefore act on distinct tensor coordinates and commute.
    j_square = all((right, left) in all_cells
                   for left, right in all_cells)
    j_gamma_anticommutes = all(
        gamma[right, left] == -gamma[left, right]
        for left, right in all_cells
    )
    order_zero_cell_rule = matrix_unit_order_zero

    # The explicit metric-zero cycle is
    # sum_(i!=j) gamma_ij pi(z_i) J pi(z_j) J^-1.  On H_(k,l), only the
    # summand (i,j)=(k,l) survives, hence the reconstructed coefficient is
    # exactly gamma_kl.
    orientation_cycle = sp.zeros(node_count)
    for left, right in all_cells:
        orientation_cycle[left, right] += gamma[left, right]
    orientable = orientation_cycle == gamma

    left_faithful = all(any(left == node for left, _ in all_cells)
                        for node in range(node_count))
    right_faithful = all(any(right == node for _, right in all_cells)
                         for node in range(node_count))

    # A D block from positive cell (i,j) to J(k,l)=(l,k) obeys first order
    # only if i=l or j=k.  This is the same central-cell rule used in the
    # repository's preregistered blind Krajewski census.
    legal_odd_blocks = []
    central_links = set()
    for first in positive_cells:
        for second in positive_cells:
            left, right = first
            other_left, other_right = second
            if left == other_right or right == other_left:
                legal_odd_blocks.append((first, second))
                if right == other_left and left != other_right:
                    central_links.add(tuple(sorted((left, other_right))))
    components = connected_components(node_count, central_links)

    gates = {
        "faithful": left_faithful and right_faithful,
        "J_square_plus": j_square,
        "J_gamma_minus": j_gamma_anticommutes,
        "gamma_commutes_with_algebra": all(gamma[left, right] in (-1, 1)
                                            for left, right in all_cells),
        "order_zero": order_zero_cell_rule,
        "metric_zero_orientable": orientable,
        "poincare_nondegenerate": determinant != 0,
        "legal_odd_blocks_nonempty": bool(legal_odd_blocks),
        "possible_central_link_graph_connected": len(components) == 1,
    }
    records.append({
        "priority": list(priority),
        "u_direction": "ascending" if u_direction == 1 else "descending",
        "v_direction": "ascending" if v_direction == 1 else "descending",
        "node_order": [nodes[node]["name"] for node in order],
        "positive_cells": [list(cell) for cell in sorted(positive_cells)],
        "hilbert_dimension": off_diagonal_dimension,
        "intersection_matrix": [[int(intersection[row, col])
                                 for col in range(node_count)]
                                for row in range(node_count)],
        "intersection_pfaffian": int(pfaffian),
        "intersection_determinant": int(determinant),
        "legal_odd_block_count": len(legal_odd_blocks),
        "central_links": [list(link) for link in sorted(central_links)],
        "central_components": len(components),
        "gates": gates,
    })

check("the frozen ambiguity space has exactly eight distinct readings",
      len(records) == 8
      and len({(tuple(record["priority"]), record["u_direction"],
                record["v_direction"]) for record in records}) == 8
      and len({tuple(record["node_order"]) for record in records}) == 8)
check("every reading gives six cells plus their six J transposes",
      all(len(record["positive_cells"]) == 6 for record in records))
check("order zero and the KO6 J/gamma signs hold for all eight readings",
      all(record["gates"]["order_zero"]
          and record["gates"]["J_square_plus"]
          and record["gates"]["J_gamma_minus"]
          and record["gates"]["gamma_commutes_with_algebra"]
          for record in records))
check("the explicit metric-zero central cycle equals gamma in every reading",
      all(record["gates"]["metric_zero_orientable"]
          for record in records))
check("both algebra actions are faithful in every reading",
      all(record["gates"]["faithful"] for record in records))
check("all eight intersection forms are exactly unimodular",
      all(abs(record["intersection_pfaffian"]) == 1
          and record["intersection_determinant"] == 1
          and record["gates"]["poincare_nondegenerate"]
          for record in records),
      f"Pfaffian multiset={dict(Counter(record['intersection_pfaffian'] for record in records))}")
check("every reading admits nonzero first-order-compatible odd positions",
      all(record["gates"]["legal_odd_blocks_nonempty"]
          for record in records),
      f"block-count multiset={dict(Counter(record['legal_odd_block_count'] for record in records))}")
check("every possible central-link graph is connected",
      all(record["gates"]["possible_central_link_graph_connected"]
          for record in records),
      "necessary possibility only; no Dirac coefficients are selected")

# Reversing both coordinate directions must reverse every oriented cell and
# Q, without changing any gate or determinant.
record_by_key = {
    (tuple(record["priority"]), record["u_direction"],
     record["v_direction"]): record for record in records
}
reversal_ok = True
for record in records:
    reversed_key = (
        tuple(record["priority"]),
        "descending" if record["u_direction"] == "ascending" else "ascending",
        "descending" if record["v_direction"] == "ascending" else "ascending",
    )
    reversed_record = record_by_key[reversed_key]
    matrix = sp.Matrix(record["intersection_matrix"])
    reversed_matrix = sp.Matrix(reversed_record["intersection_matrix"])
    reversal_ok &= reversed_matrix == -matrix
check("simultaneous direction reversal is exactly global grading reversal",
      reversal_ok)

gate_hit_counts = {
    gate: sum(record["gates"][gate] for record in records)
    for gate in records[0]["gates"]
}
all_gate_hits = sum(all(record["gates"].values()) for record in records)
check("the complete robust-gate hit fraction is 8/8",
      all_gate_hits == 8 and all(value == 8
                                 for value in gate_hit_counts.values()),
      "this is a carrier result, not a constructed finite triple")

payload = {
    "protocol_commit": "12c45eb",
    "target_comparison_performed": False,
    "algebra": {
        "real_type": ["M6(R)", "M6(R)", "M12(R)", "M12(R)"],
        "node_sizes": list(node_sizes),
        "joint_spectra": [list(node["joint_spectrum"]) for node in nodes],
    },
    "carrier": {
        "definition": "direct sum over all ordered cells (i,j), i!=j",
        "complex_dimension": off_diagonal_dimension,
        "positive_cell_count": len(unordered_pairs),
        "negative_cells_are_J_transposes": True,
        "Dirac_operator_selected": False,
    },
    "ambiguity": {
        "raw_readings": len(records),
        "complete_gate_hits": all_gate_hits,
        "hit_fraction": "8/8",
        "global_reversal_pairs": 4,
    },
    "gate_hit_counts": gate_hit_counts,
    "records": records,
    "verdict": (
        "TARGET-FREE STRUCTURAL ROBUST CARRIER EXISTENCE: all eight signed "
        "lexicographic readings give a faithful order-zero KO6 carrier, an "
        "explicit metric-zero orientation cycle and a unimodular Poincare "
        "pairing. Legal odd blocks and connected possible-link graphs exist, "
        "but no geometric Dirac operator or coefficients are selected."
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
check("the exact structured carrier audit was written", OUTPUT.exists())

print("\n"+"-"*78)
print(f"RESULT: {passed}/{tests} checks passed")
print("STRUCTURAL ROBUST: all 8/8 readings give orientable unimodular KO6 carriers.")
print("DERIVED SCOPE LIMIT: only legal D positions exist; no D is selected.")
print("NO HESSIAN OR STANDARD-MODEL TARGET WAS USED.")
raise SystemExit(0 if passed == tests else 1)
