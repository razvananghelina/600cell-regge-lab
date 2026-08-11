# Preregistration: complete orbit census of the B1 chamber embedding

Date: 2026-08-11

## Why this test is being run

The exact counterexample in `chamber_b1_refutation.md` proves that the fixed
120-chamber geometry admits a noncommutative all-gate algebra

```text
A = M2(C) + C + C + C.
```

That witness is not yet physical evidence: its positive-sheet Krajewski
cells and their dimensions were obtained during a falsification search, and
the displayed chamber colouring has trivial `A5` stabilizer.  Calling one
such colouring geometrically selected would therefore be fitting.

There is one sharply falsifiable rescue.  The 60 rotations carry the known
colouring through an orbit of 60 valid colourings.  If the fixed graph and
the already frozen cell data admit exactly that orbit and nothing else, then
the **embedding**, conditional on those cell data, is selected up to the
geometry's symmetry.  Choosing one member could then be read as spontaneous
symmetry breaking.  If further orbits exist, the embedding remains a
look-elsewhere construction.

The known witness and its 60-element orbit were inspected before this
protocol.  This is an adversarial uniqueness test, not a blind discovery.

## Frozen geometry

Reconstruct independently:

1. the exact 120 complete flags of the icosahedron;
2. chamber adjacency `D`, orientation grading `gamma`, and central
   reflection `J`;
3. the positive-sheet graph

   ```text
   S = (D J)|H_plus,
   ```

   with 60 vertices and 90 edges;
4. the full graph automorphism group `Aut(S)`, not merely the 60 rotations
   supplied in advance, and its relation to the geometric `A5` action.

No floating-point graph matching is allowed after the exact combinatorial
carrier has been constructed.

## Frozen algebra and colouring constraints

The node sizes and oriented cells are taken verbatim from the committed B1
certificate:

```text
node sizes = (2,1,1,1),
cells      = ((0,1)x2, (1,2)x25, (3,1)x12, (2,3)x19),
capacities = (4,25,12,19).
```

Colour the 60 vertices by the four cells.  Every graph edge must have one of
the exact unordered colour pairs

```text
{(0,1), (1,2), (1,3), (2,3)}.
```

The capacities are pairwise distinct, so no nontrivial colour relabelling is
quotiented out.  Only graph automorphisms count as equivalences.

These capacities and the algebra type are **STRUCTURAL INPUTS**, not outputs
of the chamber geometry.  Even a unique embedding orbit cannot retroactively
derive them.

## Complete enumeration and exact gates

Enumerate every colouring satisfying the capacity and edge constraints,
with no solution-count limit and no solver timeout in the accepted run.
Record:

- `N_support`, the total number of support-compatible colourings;
- the complete orbit decomposition under `Aut(S)`;
- every orbit size and stabilizer order;
- whether `Aut(S)` is exactly the geometric `A5` action;
- whether the committed witness orbit is the only orbit.

For one representative of every orbit, reconstruct the full represented
algebra on `C^120` and check exactly:

1. faithfulness, unitality, star closure and noncommutativity;
2. order zero and full first order on the seven complex algebra basis
   elements;
3. `[gamma,A]=0` and nonzero inner one-forms;
4. the explicit metric-dimension-zero orientation cycle;
5. the exact antisymmetric intersection matrix, rank and determinant;
6. connectedness by the exact rational rank of
   `a -> [D,pi(a)]`.

Record separately `N_all_gate` and its orbit decomposition.  A support
colouring is not silently promoted to a valid triple if connectedness or a
full tensor equation fails.

## Decision boundary

- **DERIVED CONDITIONAL EMBEDDING ORBIT:** `N_all_gate=60`, the full graph
  automorphism group is geometric `A5`, and all 60 solutions form exactly the
  committed free orbit.  This selects the embedding only conditional on the
  structural algebra/cell/multiplicity input.
- **PATTERN / LOOK-ELSEWHERE:** more than one all-gate orbit exists.  Report
  the exact hit fraction of the committed orbit and do not call it selected.
- **REFUTATION OF THE OLD CERTIFICATE:** the committed orbit fails any exact
  gate when rebuilt independently.
- **INCOMPLETE:** a resource limit is reached before a complete census.  A
  lower bound on solutions is not a uniqueness result.

No Standard-Model character, mass, coupling, generation count or
phenomenological target will be examined.  Only the targeted verifier will
be run; the full suite remains excluded by user instruction.

## Hostile scope statement

Even the strongest positive outcome would not derive the algebra type
`M2(C)+C^3`, its Krajewski support, or the capacity vector.  It would show
only that, once those previously searched structural inputs are fixed, their
embedding into the chamber geometry has no further discrete ambiguity beyond
the exact symmetry orbit.  This distinction must remain in the verdict.
