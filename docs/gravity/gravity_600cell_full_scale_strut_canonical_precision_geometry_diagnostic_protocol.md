# Preregistered geometry-control diagnostic

Date: 2026-08-20  
Status: **preregistered diagnostic; no scientific criterion changes**

## Question

Which exact conjunct makes the P160 multiprecision geometry control fail in the
first completed run frozen at commit
`bc343e117df1eccc921e3884ee4c2f488314fc4b`?

## Frozen hypotheses and exclusions

The mathematical object, source matrices, precision levels, finite-difference
steps, interval arithmetic, singular-value calculations, candidate vector,
thresholds, classification hierarchy, and outcome logic remain unchanged.
This diagnostic may not repair a failing condition and may not introduce a new
criterion.

## Instrumentation to add

For each pair in `{P100,P160} x {even,odd}`, serialize the raw value and Boolean
result for every existing geometry-control conjunct:

1. logical carrier columns cover exactly the 120 vertices;
2. the sparse carrier has exactly 4440 entries;
3. sector dimensions equal `[1,1,1,2,2,2,3]`;
4. the branch-entry control passes;
5. base negative-count histogram equals `{1:2400}`;
6. displaced negative-count histogram equals `{1:1600}`;
7. the maximum imaginary kernel residue is below the already frozen arithmetic
   floor for that precision level.

The console check will include a compact diagnostic summary, and the JSON artifact
will contain the complete records. No threshold or pass/fail expression changes.

## Preregistered interpretation

- If exactly one or more conjuncts fail, report them mechanically; do not repair
  them in the same execution.
- If the serialized conjuncts all pass while the aggregate control fails, the
  instrumentation is inconsistent and the resolver remains **OPEN**.
- A correction, if justified, requires a separate written protocol and commit
  after this diagnostic artifact is frozen.
