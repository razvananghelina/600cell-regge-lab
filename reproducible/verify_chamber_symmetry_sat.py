#!/usr/bin/env python3
"""SMT decision: can a first-order orientable chamber partition retain symmetry?

The query ranges over every set partition of the 120 oriented barycentric
chambers, encoded by a restricted-growth string.  It is stronger than a
search in the 90-bit contraction family.  For each prime-order conjugacy
class of A5 it asks whether the partition equivalence relation is invariant.
"""

import contextlib
import io
import os
import runpy
import time

import z3


HERE = os.path.dirname(os.path.abspath(__file__))
with contextlib.redirect_stdout(io.StringIO()):
    data = runpy.run_path(os.path.join(HERE, "verify_oriented_chamber_double.py"))

edges = data["chamber_edges"]
reflection = data["reflection"]
gamma = data["gamma"]
rotations = data["chamber_rotations"]
N = 120


def compose(p, q):
    return tuple(p[q[i]] for i in range(N))


identity = tuple(range(N))


def permutation_order(p):
    current = identity
    for order in range(1, 61):
        current = compose(p, current)
        if current == identity:
            return order
    raise AssertionError("rotation order exceeds 60")


representative = {}
for p in rotations:
    order = permutation_order(p)
    if order in (2, 3, 5) and order not in representative:
        representative[order] = p

# In A5 the two order-five classes are represented by g and g^2.
g5 = representative[5]
queries = (
    ("C2", representative[2]),
    ("C3", representative[3]),
    ("C5A", g5),
    ("C5B", compose(g5, g5)),
)

positive = [i for i in range(N) if int(gamma[i]) == 1]
negative = [i for i in range(N) if int(gamma[i]) == -1]
timeout_ms = int(os.environ.get("CHAMBER_SYMMETRY_TIMEOUT_MS", "120000"))


def solve_query(name, permutation):
    labels = [z3.Int(f"{name}_p_{i}") for i in range(N)]
    maxima = [z3.Int(f"{name}_m_{i}") for i in range(N)]
    solver = z3.Solver()
    solver.set(timeout=timeout_ms)

    # Every set partition exactly once.
    solver.add(labels[0] == 0, maxima[0] == 0)
    for i in range(1, N):
        solver.add(labels[i] >= 0)
        solver.add(labels[i] <= maxima[i-1] + 1)
        solver.add(maxima[i] == z3.If(labels[i] > maxima[i-1],
                                      labels[i], maxima[i-1]))

    # First order and metric-dimension-zero orientability.
    for x, y in edges:
        solver.add(z3.Or(labels[x] == labels[y],
                         labels[reflection[x]] == labels[reflection[y]]))
    for x in positive:
        for y in negative:
            solver.add(z3.Or(labels[x] != labels[y],
                             labels[reflection[x]] != labels[reflection[y]]))

    # The symmetry may permute block labels.  What must be invariant is the
    # equivalence relation, not the restricted-growth integer names.
    for x in range(N):
        for y in range(x+1, N):
            solver.add((labels[x] == labels[y])
                       == (labels[permutation[x]] == labels[permutation[y]]))

    start = time.monotonic()
    result = solver.check()
    elapsed = time.monotonic() - start
    if result == z3.sat:
        model = solver.model()
        values = tuple(model.eval(label).as_long() for label in labels)
        return "sat", elapsed, values, ""
    if result == z3.unsat:
        return "unsat", elapsed, None, ""
    return "unknown", elapsed, None, solver.reason_unknown()


print("=" * 78, flush=True)
print("CHAMBER PARTITION SYMMETRY SAT AUDIT", flush=True)
print("scope: arbitrary partitions; first order + orientability", flush=True)
results = {}
for name, permutation in queries:
    status, elapsed, labels, reason = solve_query(name, permutation)
    results[name] = (status, labels)
    print(f"{name}: {status} in {elapsed:.2f}s", flush=True)
    if reason:
        print(f"  reason={reason}", flush=True)
    if labels is not None:
        print(f"  blocks={len(set(labels))}", flush=True)
        print("  labels=" + ",".join(map(str, labels)), flush=True)

unknown = [name for name, (status, _) in results.items() if status == "unknown"]
sat = [name for name, (status, _) in results.items() if status == "sat"]
if unknown:
    print("STATUS=INCOMPLETE", flush=True)
    print("unknown=" + ",".join(unknown), flush=True)
    raise SystemExit(2)
if sat:
    print("STATUS=SYMMETRIC_LEGAL_PARTITION_EXISTS", flush=True)
    print("sat=" + ",".join(sat), flush=True)
else:
    print("STATUS=DERIVED_NO_GO", flush=True)
    print("No first-order orientable partition has a nontrivial A5 stabilizer.",
          flush=True)

print("RESULT: 4/4 prime-order A5 conjugacy-class queries decided", flush=True)
