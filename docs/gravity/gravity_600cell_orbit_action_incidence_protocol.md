# Protocol: exact flag-incidence audit of the reduced slab action

Date: 2026-08-21

Prior-art gate: commit `c14b5ac`.

Frozen discrepancy artifact: commit `af862ab`.

Status: frozen before enumerating any triangle--simplex flag orbit.

Only the new targeted verifier and the static registry audit may be run.  No
existing evaluator or artifact may be changed during this audit.

## 1. Carrier and coordinate reconciliation

Load definitions only from the 100-decimal orbit-action source and from the
binary64 direct-action source.  Do not run either main program.  For each
parity require exact equality of:

- the 2400 labelled simplices;
- the order-24 stabilizer permutations;
- old, internal and final edge sets;
- triangle sets and boundary-triangle sets;
- old/internal/final coordinate lookup maps after removing the direct source's
  fixed final-index offset;
- every edge Jacobian/sign convention.

Reconstruct the frozen regular and off-shell states independently.  Expand
each orbit coordinate to every labelled edge and require exact invariance under
all 24 stabilizer elements.  A failure here forbids an incidence verdict.

## 2. Exact flag enumeration

Let

```text
F = {(triangle,simplex): triangle is a three-vertex face of simplex}.
```

Require `|F|=24000`.  Enumerate all orbits of `F` under the stabilizer by exact
integer permutations.  Also enumerate triangle and simplex orbits independently
from the labelled sets rather than trusting stored orbit lists.

For every flag orbit `F_a`, identify its triangle orbit `T_i` and simplex orbit
`S_j` and compute exactly

```text
c_a = |F_a|/|T_i|.
```

Require `c_a` to be a positive integer by direct divisibility.  Verify by
choosing every triangle in `T_i` and counting incident flags from `F_a`; all
counts must equal `c_a`.

Record before any action comparison:

- numbers and size distributions of triangle, simplex and flag orbits;
- the complete coefficient multiset;
- the complete list of `(T_i,S_j,c_a)` records;
- a SHA-256 digest of the canonical JSON representation.

## 3. Audit the existing shortcut

Reproduce mechanically the shortcut's incidence rule: for each stored simplex
orbit choose its lexicographic representative and include each of its ten
local triangle flags once.  Map each such flag to its exact flag orbit and
count the resulting shortcut coefficient `s_a`.

Compare the full vectors `(s_a)` and `(c_a)` exactly.  No action value enters
this comparison.  Report all mismatching flag orbits, not only their number or
maximum discrepancy.

As independent counting controls require

```text
sum_a |F_a| = 24000,
sum_a c_a |T_i(a)| = 24000,
sum_a s_a = 10 * number_of_simplex_orbits.
```

## 4. Action discriminator

At 80 decimal digits implement two reduced actions using one locally coded
Lorentzian angle routine:

1. `S_flag`, with the exact `c_a` coefficients and triangle-orbit sizes;
2. `S_shortcut`, with the existing representative-simplex rule.

Evaluate both at the regular state and the already frozen off-shell state for
both parities.  No root, perturbation, scale or tolerance may be selected from
output.

Only after all eight reduced action values are complete, load the three frozen
artifacts.  Require:

- at the off-shell state, `S_flag` agrees with the stored 80-decimal direct
  action within `1e-45`;
- `S_shortcut` agrees with the stored primary orbit action within `1e-45`;
- the two reproduce the stored direct/orbit disagreement within `1e-45` in
  each parity;
- at the regular state, both reduced actions agree with each other and the
  stored arbitrary-precision published control within `1e-45`.

Every evaluated representative simplex must retain one Lorentzian negative
direction, minimum nonzero leading-minor modulus above `1e-20`, minimum angle
argument modulus above `1e-6`, and imaginary action below `1e-60`.

## 5. Outcome hierarchy

Assign exactly one:

1. `ORBIT_ACTION_INCIDENCE_CONTROL_FAILED` for a carrier, invariance,
   enumeration, branch, provenance or known-control failure;
2. `ORBIT_ACTION_DIRECT_CONSTRUCTION_SUSPECT` if `s_a=c_a` exactly but the
   existing disagreement remains;
3. `ORBIT_ACTION_FLAG_MULTIPLICITY_BUG_DERIVED` if `s_a!=c_a`, `S_flag`
   matches the direct action, `S_shortcut` matches the primary orbit action,
   and all regular controls pass;
4. `ORBIT_ACTION_DISAGREEMENT_UNRESOLVED` otherwise.

Outcome 3 establishes an exact implementation bug and requires a separate
scope audit before any affected nonhomogeneous dynamics result is reused.  It
does not invalidate direct Regge calculus or the independent scale-homogeneity
theorem.  No existing verifier may be silently patched as part of this result.

