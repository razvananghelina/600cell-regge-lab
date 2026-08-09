#!/usr/bin/env python3
"""
SMT certificate attempt for the global lower bound dim Omega^1 >= 126.

Variables encode an arbitrary partition of all 120 oriented chambers using
a restricted-growth string, so every set partition occurs exactly once.
Constraints impose:

  * first order on every chamber edge;
  * metric-dimension-zero orientability on every opposite-orientation pair.

The objective is encoded as the number of distinct unordered quotient-graph
edges.  Since D is self-adjoint, dim_C Omega^1 is twice this number.

UNSAT for at most 62 unordered quotient edges proves the global lower bound
dim Omega^1 >= 126 for every orientable first-order partition algebra.
"""

import os
import runpy
import time

import z3


HERE = os.path.dirname(os.path.abspath(__file__))
data = runpy.run_path(os.path.join(HERE, "verify_oriented_chamber_double.py"))

chamber_edges = data["chamber_edges"]
reflection = data["reflection"]
gamma = data["gamma"]
N = 120

labels = [z3.Int(f"p_{i}") for i in range(N)]
maxima = [z3.Int(f"m_{i}") for i in range(N)]
s = z3.Solver()
s.set(timeout=int(os.environ.get("CHAMBER_SMT_TIMEOUT_MS", "600000")))

# Restricted-growth canonical labeling of an arbitrary set partition.
s.add(labels[0] == 0, maxima[0] == 0)
for i in range(1, N):
    s.add(labels[i] >= 0)
    s.add(labels[i] <= maxima[i - 1] + 1)
    s.add(maxima[i] == z3.If(labels[i] > maxima[i - 1],
                              labels[i], maxima[i - 1]))

# First order for a partition algebra: each D-edge is collapsed either by
# the represented partition or by its J-reflected opposite partition.
for x, y in chamber_edges:
    s.add(z3.Or(labels[x] == labels[y],
                labels[reflection[x]] == labels[reflection[y]]))

# Orientability: (P(x),P(Jx)) must determine the sign gamma_x.  Thus no
# positive and negative chamber may have the same ordered block pair.
positive = [i for i in range(N) if int(gamma[i]) == 1]
negative = [i for i in range(N) if int(gamma[i]) == -1]
for x in positive:
    for y in negative:
        s.add(z3.Or(labels[x] != labels[y],
                    labels[reflection[x]] != labels[reflection[y]]))

# Distinct unordered quotient edges.  Internal chamber edges do not
# contribute to represented one-forms.  The first occurrence of each
# unordered label pair contributes one.
active = []
low = []
high = []
new_pair = []
for e, (x, y) in enumerate(chamber_edges):
    active_e = labels[x] != labels[y]
    low_e = z3.If(labels[x] < labels[y], labels[x], labels[y])
    high_e = z3.If(labels[x] < labels[y], labels[y], labels[x])
    previous_different = [
        z3.Not(z3.And(active[k],
                      low_e == low[k],
                      high_e == high[k]))
        for k in range(e)
    ]
    new_e = z3.And(active_e, *previous_different)
    active.append(active_e)
    low.append(low_e)
    high.append(high_e)
    new_pair.append(new_e)

unordered_quotient_edges = z3.Sum(
    [z3.If(flag, 1, 0) for flag in new_pair])

# A counterexample to the proposed 126 lower bound would have Omega <=124,
# equivalently at most 62 unordered quotient edges.
s.add(unordered_quotient_edges <= 62)

print("=" * 78)
print("GLOBAL SMT MINIMUM SCREEN: seeking dim Omega <= 124")
print(f"partition variables={N}, chamber edges={len(chamber_edges)}")
print(f"orientability cross-sign constraints={len(positive)*len(negative)}")
start = time.monotonic()
result = s.check()
elapsed = time.monotonic() - start
print(f"SMT_RESULT={result}")
print(f"elapsed={elapsed:.2f}s")
if result == z3.unsat:
    print("GLOBAL_LOWER_BOUND=126")
    print("STATUS=CERTIFIED_GLOBAL_MINIMUM (combined with the Omega=126 witness)")
elif result == z3.sat:
    model = s.model()
    values = [model.eval(label).as_long() for label in labels]
    quotient_edges = model.eval(unordered_quotient_edges).as_long()
    print(f"COUNTEREXAMPLE_OMEGA={2*quotient_edges}")
    print("COUNTEREXAMPLE_LABELS=" + ",".join(map(str, values)))
    print("STATUS=126_REFUTED")
    raise SystemExit(1)
else:
    print("STATUS=INCOMPLETE -- solver returned unknown")
    print(f"REASON={s.reason_unknown()}")
    raise SystemExit(2)
