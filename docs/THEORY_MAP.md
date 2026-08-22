# Authoritative theory map

Updated: 2026-08-22

This is the navigation and duplication-control map for the active public
600-cell Regge programme.  It is not a promotional summary.  Every route has
a complete scope in [`theory_map.json`](theory_map.json), including its
hypotheses, evidence, dependencies, search aliases and, for a bounded no-go,
the condition required to reopen it.

Read this file, [`gravity/CURRENT_STATUS.md`](gravity/CURRENT_STATUS.md) and
the binding [`../CLAUDE.md`](../CLAUDE.md) before starting a calculation.

## How to use the map

Before a new calculation:

1. formulate the proposed object, carrier, operator and hypotheses;
2. search `theory_map.json` using all technical aliases of that object;
3. search the referenced result notes and verifier names before treating the
   route as new;
4. classify the proposal as an existing result, a reusable control, a
   genuinely open gate or a route already closed under the same hypotheses;
5. only then perform the literature gate and preregister a new test;
6. keep exactly one route marked `ACTIVE_GATE`;
7. when the result is consolidated, update this map and the JSON registry in
   the same commit.

Absence from this map is not evidence of novelty.  It means the repository
search must be widened before a new route is authorized.

The scientific label and route state answer different questions:

- `DERIVED`, `STRUCTURAL`, `PATTERN`, `OPEN` say how strong the evidence is;
- `ACCEPTED`, `BOUNDED_NO_GO`, `OPEN_GATE`, `ACTIVE_GATE`, `METHOD_CONTROL`
  and `PATTERN_CONTROL` say what to do with the route next.

A `BOUNDED_NO_GO` is never an unrestricted impossibility claim.  Its exact
kill scope and reopening condition are in the JSON registry.

## Main map

```mermaid
flowchart TB
    F0["F-CARRIER<br/>600-cell spatial S3 carrier"]
    F1["F-ACTION<br/>Lorentzian Regge + conserved dust"]

    subgraph H["Homogeneous finite-height dynamics"]
        H1["H-ONE-SLAB<br/>classified finite-height slabs"]
        H2["H-NONUNIQUE<br/>two physical second branches"]
        H3["H-LOCAL-SELECTORS<br/>standard local selectors fail"]
        H4["H-EXTENDIBILITY<br/>branch B survives"]
        H5["H-INVARIANT-REGION<br/>unique successor at every finite step"]
        H6["H-INFINITE-LIMIT<br/>convergence / total duration open"]
        H7["H-UNIVERSAL-LIMIT<br/>continuous fixed family"]
        H8["H-GENERIC-BASIN<br/>generic incoming basin open"]
    end

    subgraph NH["Nonhomogeneous canonical response"]
        N1["NH-CARRIER<br/>rank-240 geometry-selected carrier"]
        N2["NH-INTERNAL-RANK<br/>rank 239; one homogeneous line"]
        N3["NH-FREE-TICK<br/>line is not a free tick"]
        N4["NH-FIRST-TANGENT<br/>first 1440D symplectic response"]
        N5["NH-SECOND-TANGENT<br/>accepted two-step response"]
        N6["P-CONSTRAINT-QUOTIENT<br/>no exact coarse quotient"]
        N6R["P-REFINEMENT-CONSTRAINT<br/>ACTIVE: gauge restoration"]
        N7["P-WAVES<br/>local waves / spatial operator open"]
        N8["P-C<br/>limiting speed open"]
    end

    subgraph CTRL["Repository controls: do not rediscover"]
        C1["CTRL-OLD-TWO-STEP<br/>old near-static two-step cocycle"]
        C2["CTRL-OLD-JACOBI<br/>old acceleration-drift-stiffness form"]
        C3["R-FRIEDMANN<br/>three-level refinement pattern"]
        C4["CTRL-REFINED-H4-CONSTRAINT<br/>single refined constrained control"]
        C5["CTRL-NESTED-TANGENTIAL-CARRIER<br/>canonical K0/K1 displacement map"]
    end

    subgraph SCALE["Scale and downstream physics"]
        S1["S-SCALE-NOGO<br/>no absolute classical tick"]
        S2["P-G-PLANCK<br/>G / Planck units open"]
        S3["P-MATTER-SM<br/>particles / SM open"]
    end

    F0 --> F1 --> H1 --> H2 --> H3 --> H4 --> H5
    H5 --> H6
    H4 --> H7
    H5 --> H8

    H1 --> N1 --> N2 --> N3 --> N4 --> N5 --> N6 --> N6R --> N7 --> N8
    H4 --> N5

    F1 --> C1 --> C2
    C1 -. method control .-> N5
    F1 --> C3
    F1 --> C4
    C3 -. refinement pattern .-> N6R
    C4 -. singular control .-> N6R
    C5 -. matched configuration carrier .-> N6R

    F1 --> S1
    S1 -. physical scale lift .-> N5
    S1 --> S2
    C3 --> S2
    N7 --> S2
    N6R --> S3

    classDef foundation fill:#d9edf7,stroke:#31708f,color:#111;
    classDef accepted fill:#dff0d8,stroke:#3c763d,color:#111;
    classDef noGo fill:#f2dede,stroke:#a94442,color:#111;
    classDef open fill:#fcf8e3,stroke:#8a6d3b,color:#111;
    classDef active fill:#eadcff,stroke:#6f42c1,stroke-width:4px,color:#111;
    classDef control fill:#eeeeee,stroke:#666,color:#111;
    classDef pattern fill:#fff1b8,stroke:#9a7b00,color:#111;

    class F0,F1 foundation;
    class H1,H4,H5,N1,N2,N4,N5 accepted;
    class H2,H3,H7,N3,N6,S1 noGo;
    class H6,H8,N7,N8,S2,S3 open;
    class N6R active;
    class C1,C2,C4,C5 control;
    class C3 pattern;
```

The red nodes are bounded negatives, not failed projects.  They tell us
which exact argument must not be rerun under unchanged hypotheses.  The
purple node is the only active calculation.

## Route ledger

| Route ID | Evidence label | Route state | Current statement | Primary evidence |
|---|---|---|---|---|
| `F-CARRIER` | STRUCTURAL | FOUNDATION | Fixed regular 600-cell spatial `S^3` carrier, f-vector `(120,720,1200,600)` | [face-poset result](gravity/gravity_600cell_overlay_face_poset_result.md) |
| `F-ACTION` | DERIVED | FOUNDATION | Reproducible zero-`Lambda` Lorentzian Regge-plus-dust slab action | [published-action control](gravity/gravity_600cell_published_dust_control_result.md) |
| `H-ONE-SLAB` | DERIVED | ACCEPTED | Nonzero positive-height one-slab roots are exactly classified in the homogeneous sector | [classification](gravity/gravity_600cell_finite_height_classification_result.md) |
| `H-NONUNIQUE` | DERIVED | BOUNDED_NO_GO | At `v=3/2`, the action alone gives two physical second slabs | [composition result](gravity/gravity_600cell_finite_height_composition_result.md) |
| `H-LOCAL-SELECTORS` | DERIVED | BOUNDED_NO_GO | Causality, orientation, real branch and local regularity do not select between them | [selector audit](gravity/gravity_600cell_finite_height_selector_result.md) |
| `H-EXTENDIBILITY` | STRUCTURAL | ACCEPTED | Branch A dies; branch B has the accepted third, fourth and fifth continuations | [third-slab result](gravity/gravity_600cell_finite_height_third_slab_result.md) |
| `H-INVARIANT-REGION` | DERIVED | ACCEPTED | The accepted branch has one successor at every prescribed finite step | [invariant-region theorem](gravity/gravity_600cell_finite_height_invariant_region_result.md) |
| `H-INFINITE-LIMIT` | OPEN | OPEN_GATE | Convergence, infinite total proper duration and completeness are not proved | [invariant-region limits](gravity/gravity_600cell_finite_height_invariant_region_result.md) |
| `H-UNIVERSAL-LIMIT` | DERIVED | BOUNDED_NO_GO | The compactified boundary contains a continuous fixed family, not one universal endpoint | [asymptotic map](gravity/gravity_600cell_finite_height_asymptotic_map_result.md) |
| `H-GENERIC-BASIN` | OPEN | OPEN_GATE | The finite census has 36 signatures; only a local neighbourhood of the representative tree is certified | [basin census](gravity/gravity_600cell_finite_height_incoming_basin_discovery_result.md) |
| `NH-CARRIER` | STRUCTURAL | ACCEPTED | A geometry-selected rank-240 scale-plus-strut tangent carrier survives the canonicity audit | [quadratic carrier](gravity/gravity_600cell_finite_height_carrier_quadratic_result.md) |
| `NH-INTERNAL-RANK` | DERIVED | ACCEPTED | Internal stationarity has rank 239 and kills every nonhomogeneous carrier direction | [rank result](gravity/gravity_600cell_finite_height_internal_carrier_rank_result.md) |
| `NH-FREE-TICK` | DERIVED | BOUNDED_NO_GO | The survivor is the lapse-constraint tangent and fixed incoming momentum removes it | [kernel reconciliation](gravity/gravity_600cell_finite_height_internal_kernel_canonical_reconciliation_result.md) |
| `NH-FIRST-TANGENT` | DERIVED | ACCEPTED | The first accepted slab has an adversarially replicated regular symplectic 1440D forced response | [first full tangent](gravity/gravity_600cell_finite_height_full_boundary_tangent_result.md) |
| `NH-SECOND-TANGENT` | DERIVED | ACCEPTED | Both complete routes establish a regular canonical two-slab response, exact physical scale lift and schedule-robust four-way composition on the unreduced carrier | [consolidated result](gravity/gravity_600cell_second_full_boundary_tangent_result.md); [preserved adversarial failure](gravity/gravity_600cell_second_full_boundary_tangent_adversarial_first_run_failure.md) |
| `CTRL-OLD-TWO-STEP` | DERIVED | METHOD_CONTROL | Two-step tangent machinery already exists on the old near-static `tau=0.0102` background | [old two-step result](gravity/gravity_600cell_dust_two_step_full_tangent_result.md) |
| `CTRL-OLD-JACOBI` | DERIVED | METHOD_CONTROL | Three-slice acceleration-drift-stiffness machinery already exists on that old background | [old Jacobi result](gravity/gravity_600cell_dust_three_slice_jacobi_result.md) |
| `CTRL-REFINED-H4-CONSTRAINT` | DERIVED | METHOD_CONTROL | A different stationary barycentric `H4` sector has one exact internal null, a compatibility hyperplane and a schedule-independent constrained response | [null coupling](gravity/gravity_600cell_refined_h4_null_coupling_result.md); [constrained response](gravity/gravity_600cell_refined_h4_constrained_response_result.md) |
| `CTRL-NESTED-TANGENTIAL-CARRIER` | DERIVED | METHOD_CONTROL | The normalized `K0 -> K1` refinement canonically prolongs all 7920 spatial tangential vertex-displacement directions with an exact old-vertex left inverse | [nested carrier result](gravity/gravity_600cell_nested_vertex_displacement_result.md) |
| `R-FRIEDMANN` | PATTERN | PATTERN_CONTROL | Three refinement levels approach the homogeneous Friedmann acceleration; no convergence theorem | [refinement comparison](gravity/gravity_600cell_projected_refinement_acceleration_comparison_result.md) |
| `S-SCALE-NOGO` | DERIVED | BOUNDED_NO_GO | The scale-covariant classical action cannot select an absolute nonzero tick | [scale-covariance theorem](gravity/gravity_600cell_tick_scale_covariance_result.md) |
| `P-CONSTRAINT-QUOTIENT` | DERIVED | BOUNDED_NO_GO | The regular two-slab map has open pre/post images and no exact local quotient selected by coarse Legendre degeneracy | [bounded no-go](gravity/gravity_600cell_finite_height_constraint_quotient_result.md) |
| `P-REFINEMENT-CONSTRAINT` | OPEN | ACTIVE_GATE | A canonical K0/K1 tangential carrier exists, but no matched on-shell two-resolution family yet shows pseudo-gauge directions converging to an exact constraint kernel | [nested carrier result](gravity/gravity_600cell_nested_vertex_displacement_result.md) |
| `P-WAVES` | OPEN | OPEN_GATE | No physical graviton, wave equation or spatial tensor operator is derived | [old Jacobi firewall](gravity/gravity_600cell_dust_three_slice_jacobi_result.md) |
| `P-C` | OPEN | OPEN_GATE | No effective limiting speed is derived | [scale and Jacobi limits](gravity/gravity_600cell_tick_scale_covariance_result.md) |
| `P-G-PLANCK` | OPEN | OPEN_GATE | `G` and Planck units require both physical local dynamics and independently derived scale breaking | [scale no-go](gravity/gravity_600cell_tick_scale_covariance_result.md) |
| `P-MATTER-SM` | OPEN | OPEN_GATE | No particles, masses, gauge algebra or Standard-Model sector are derived in the active public theory | [repository scope](../README.md) |

## No-go ledger and reopening rules

| Route | Closed claim under the stated hypotheses | What would legitimately reopen it |
|---|---|---|
| `H-NONUNIQUE` | The frozen homogeneous action plus its existing physical inequalities defines a unique second tick | An independently motivated selector, or a proof that nonhomogeneous/refined equations remove the ambiguity |
| `H-LOCAL-SELECTORS` | The four tested local conditions select branch A or B | A different preregistered physical selector; do not rerun the same four tests |
| `H-UNIVERSAL-LIMIT` | The current compactified map selects one universal scale ratio | A separately derived condition that breaks the continuous fixed family |
| `NH-FREE-TICK` | The rank-239 internal kernel line is a free mode or a tick | A different, derived physical carrier or boundary-data problem; not a relabelling of this line |
| `S-SCALE-NOGO` | The current scale-free classical action derives seconds, a Planck length or an absolute tick | An independently derived dimensionful scale-breaking input; inserting a measured target after inspection is fitting |
| `P-CONSTRAINT-QUOTIENT` | The fixed coarse branch-B first/second slabs select a nontrivial exact local pre/post-constraint quotient through Legendre degeneracy | A later singular move with a propagated constraint, an exact continuous symmetry of an enlarged action, a coherent refinement limit to an exact kernel, or a different/perfect action |

## Current decision path

The only active route is `P-REFINEMENT-CONSTRAINT`.  The former active gate,
`P-CONSTRAINT-QUOTIENT`, is now a bounded no-go: all accepted first/second
pre-Legendre systems are regular by strict margins, and two independent
finite-dimensional proofs give zero local pre/post-constraint codimension.
Resolved nonzero pseudo-constraint directions cannot be removed by a chosen
threshold without changing the finite theory.

The `CTRL-NESTED-TANGENTIAL-CARRIER` feasibility gate has passed.  Two
constructions give the same `K0 -> K1` normalized-midpoint map, its old-vertex
restriction is an exact left inverse and the spatial tangential rank is 7920.
This removes post-hoc eigenvector matching for that declared carrier, but it
does not transport momenta or the normal/lapse component.

The existing `CTRL-REFINED-H4-CONSTRAINT` result remains a different
stationary, homogeneous, single-level control.  The active route must now
construct a matched on-shell `K1` finite-height seed with the same action,
physical radius and conserved-matter normalization as `K0`.  Only after that
may it preregister singular-value/coupling scaling toward an exact kernel.  If
the matched seed does not exist or temporal schedule ambiguity survives at
the effective boundary level, the refinement route closes; comparing
unmatched Hessians or post-hoc mode overlaps is forbidden.

Until this gate passes, `P-WAVES`, `P-C`, `P-G-PLANCK` and `P-MATTER-SM`
remain blocked.  The accepted tangent is not a graviton, a value of `c`, a
physical tick or local general relativity.

## Update invariant

The machine verifier requires:

- unique route IDs and exactly one `ACTIVE_GATE`;
- valid evidence labels and route states;
- existing evidence files and valid dependency IDs;
- no dependency cycle;
- hypotheses, scope and duplicate-search terms for every route;
- kill scopes and reopening conditions for every bounded no-go;
- a next test for every open or active gate;
- every route ID to occur in this Markdown map;
- links to this map from the curated repository indexes and binding rules;
- exactly one registry entry for the map verifier.

Run only:

```bash
/home/razvan/science/.venv/bin/python reproducible/verify_theory_map.py
```

This governance verifier does not run any scientific calculation.
