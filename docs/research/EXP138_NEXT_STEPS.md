# Post-EXP138: Suggested Next Experiments

**Context**: Exp138 showed that fermion loops on 600-cell vertices CANNOT generate Yang-Mills gauge structure. Gauge bosons need different geometric realization.

---

## EXP139: Edge-Based Gauge Theory (HIGHEST PRIORITY)

**Hypothesis**: Gauge bosons live on EDGES, fermions on VERTICES (standard lattice gauge theory).

### Setup
- 600-cell: 120 vertices, 720 edges
- **Fermions**: Type-C vertices (96) - KEEP vertex assignment (works well)
- **Gauge bosons**: Edges (720 total)
- Edge types: A-C (96), B-C (192), C-C (432)

### Tasks
1. **Edge classification**
   - Each edge connects two vertices of types {A,B,C}
   - A-C edges (96): candidate SU(3) gluons?
   - B-C edges (192): candidate SU(2) W/Z bosons?
   - C-C edges (432): fermion propagators? Yukawa?

2. **Plaquette action**
   - Triangular faces (1200 from exp136)
   - Each triangle = 3 edges forming plaquette
   - Wilson loop around triangle: U = g_1 g_2 g_3
   - Field strength: F ~ Tr(1 - U)
   - Action: S = sum_triangles Tr(1 - U_triangle)

3. **Gauge transformations**
   - Gauge transformation at vertex v: g_v in SU(3) or SU(2)
   - Edge variable: U_ij → g_i U_ij g_j^†
   - Check gauge invariance of triangle action

4. **Counting**
   - ACC triangles (240): involve 1 A-vertex, 2 C-vertices, 3 edges
   - What gauge group per edge type?
   - Do 96 A-C edges form SU(3) adjoint (8-dimensional)?

5. **Fermion-gauge coupling**
   - Fermion on vertex i
   - Gauge field on edge ij
   - Fermion hopping: ψ_i U_ij ψ_j
   - Standard lattice QCD form!

**Expected outcome**:
- A-C edges: 96 = 12 * 8, possibly 8-component SU(3) gluons at each of 12 "sites"?
- B-C edges: 192 = 12 * 16, need to understand decomposition
- May naturally get gauge self-interaction (edges form closed loops)

**Category**: DERIVED (if we compute exact edge adjacency structure)

---

## EXP140: Spectral Gauge Theory (HIGH PRIORITY)

**Hypothesis**: Gauge bosons are Laplacian eigenmodes (spectral approach).

### Setup
From exp110-111:
- 600-cell Laplacian has 9 eigenvalues (exact phi-algebraic)
- Multiplicities: 1, 4, 9, 16, 25, 36, 9, 16, 4 (all perfect squares = H_4 irrep dimensions)
- Eigenspaces could encode gauge symmetry

### Tasks
1. **Gauge field as eigenfunction**
   - Laplacian eigenvector φ_k (vector on 120 vertices)
   - Interpret as gauge field configuration
   - Different eigenspaces → different gauge groups?

2. **Dirac operator** (antisymmetric structure)
   - Define Dirac operator D on 600-cell (graph Dirac operator)
   - D^† = -D (automatically antisymmetric!)
   - Commutator [D, f] for function f encodes gauge structure
   - Check if [D, φ] eigenspaces form Lie algebra

3. **Non-commutative geometry**
   - Algebra A = functions on vertices
   - Inner derivation: [D, a] for a in A
   - Check if this is isomorphic to SU(3) or SU(2) Lie algebra

4. **Eigenspace decomposition**
   - 8-dimensional subspaces? (SU(3) adjoint)
   - 3-dimensional subspaces? (SU(2) adjoint)
   - From exp111: Phi-sector dimension = 26 = 4+9+9+4
   - Does 26 decompose as 8 + 3 + 15? (SU(3) + SU(2) + ?)

5. **Connection to E8**
   - From exp120d: 600-cell embeds in E8 (S union T = 240 roots)
   - E8 → H4 branching rule: 248 = 2*120 + 8 (exp120)
   - Can we extract SU(3) x SU(2) from E8 structure?

**Expected outcome**:
- Gauge bosons = collective modes (like phonons in crystal)
- Antisymmetric Dirac operator provides Lie algebra structure
- Eigenspace multiplicities may match gauge group dimensions

**Category**: DERIVED (eigenvalues) + PATTERN (if we find group matches)

---

## EXP141: Dual Lattice (120-cell) Gauge Theory (MEDIUM PRIORITY)

**Hypothesis**: Fermions on 600-cell, gauge bosons on 120-cell (dual polytope).

### Setup
From exp117:
- 120-cell: 600 vertices, 1200 edges, degree 4
- Dual to 600-cell (vertices ↔ cells)
- Each 120-cell vertex = center of 600-cell tetrahedron

### Tasks
1. **Primal-dual correspondence**
   - 600-cell vertex ↔ 120-cell cell (tetrahedron)
   - 600-cell edge ↔ 120-cell face (triangle)
   - 600-cell face ↔ 120-cell edge
   - 600-cell cell ↔ 120-cell vertex

2. **Field assignment**
   - Fermions: 600-cell vertices (120)
   - Gauge bosons: 120-cell vertices (600)
   - Check dimensions: 600 = 8 * 75? Factor of 75 mysterious.
   - Or: 600 = 24 * 25 (24 gauge bosons, 25 polarizations/sites?)

3. **Interaction vertices**
   - 120-cell vertex (gauge boson) in center of 600-cell tetrahedron
   - Tetrahedron has 4 vertices (4 fermions)
   - Gauge-fermion coupling: 4-point vertex (unusual!)
   - Or restrict to face (triangle): 3-point vertex (QED/QCD-like)

4. **Spectral structure**
   - From exp117: 120-cell has 9 algebraic fields (not pure phi)
   - Discriminants are Fibonacci numbers (2, 5, 13, 21)
   - Less "pure" than 600-cell (which has only Q(sqrt(5)))
   - Does spectral impurity encode gauge group complexity?

**Expected outcome**:
- May naturally separate fermions (120 on 600-cell) from gauge (600 on 120-cell)
- Factor of 5 discrepancy (600 vs 120) mysterious but may be geometric
- 120-cell degree 4 (vs 600-cell degree 12) may simplify dynamics

**Category**: DERIVED (polytope structure) + SPECULATIVE (field assignment)

---

## EXP142: Holonomy and Wilson Loops (MEDIUM PRIORITY)

**Hypothesis**: Gauge bosons as holonomy (parallel transport) on 600-cell graph.

### Setup
- Connection 1-form ω on edges (SU(3) or SU(2) valued)
- Parallel transport: U(path) = P exp(∫ ω)
- Wilson loop: W(C) = Tr[U(closed path C)]

### Tasks
1. **Define connection**
   - Each edge ij: assign matrix U_ij in SU(3) or SU(2)
   - Constraint: U_ji = U_ij^† (consistency)
   - Question: What determines U_ij? (Free choice or derived from geometry?)

2. **Minimal loops**
   - Triangles (1200): smallest non-trivial loops
   - Icosahedral cycles (vertex figure): 5-gons, 6-gons?
   - Diameter-5 paths: maximal loops

3. **Curvature**
   - Field strength: F_ij = [D_i, D_j] = [d + ω, d + ω]
   - For triangle ijk: F ~ (U_ij U_jk U_ki - 1)
   - Sum over triangles: total curvature = ?
   - Compare to Euler characteristic χ = 0 (S^3)

4. **Holonomy groups**
   - For each vertex v, consider all loops through v
   - Holonomy group Hol(v) = group generated by all loop holonomies
   - Check if Hol(v) = SU(3) or SU(2) or other

5. **Gauge-invariant observables**
   - Wilson loops (already gauge invariant)
   - 't Hooft loops (magnetic dual)
   - Polyakov loops (if we have temporal direction on S^3?)

**Expected outcome**:
- Connection approach is standard in gauge theory
- Question is whether 600-cell geometry FIXES the connection or leaves freedom
- If geometry fixes it: connection is DERIVED (good!)
- If freedom remains: connection is FITTING (bad)

**Category**: Depends on whether connection is derived or chosen

---

## EXP143: Fermion Self-Energy and Emergent Gauge Bosons (LOW PRIORITY)

**Hypothesis**: Gauge bosons emerge from fermion-fermion interactions (C-C edges).

### Setup
- Start with ONLY fermions (96 Type-C vertices)
- C-C adjacency: 432 edges (fermion-fermion interactions)
- Gauge bosons = collective excitations

### Tasks
1. **Fermion propagator**
   - Free fermion on C vertices
   - Propagator: G_0(i,j) ~ (E - L)^{-1} where L = Laplacian
   - Fermion self-energy: Σ(E) from C-C interactions

2. **Collective modes**
   - Look for poles in fermion 2-point function
   - Pole at E = E_gauge → gauge boson
   - Residue → gauge coupling strength

3. **Gauge symmetry emergence**
   - Check if low-energy effective theory has gauge invariance
   - Like in spin liquids (emergent U(1) photon from spin interactions)
   - Or in strongly coupled systems (ρ meson emerges from pions)

4. **Connection to Type-A/B vertices**
   - Do collective modes "live" at Type-A/B locations?
   - Are Type-A/B vertices the SOURCES of fermions?
   - Is gauge boson = bound state of fermions from A/B sources?

**Expected outcome**:
- Highly speculative
- Would be emergent gauge theory (like in condensed matter)
- Requires solving strongly coupled fermion system (hard!)

**Category**: SPECULATIVE (emergent gauge symmetry not guaranteed)

---

## EXP144: 24-Cell as Pure Gauge Sector (WILD IDEA)

**Hypothesis**: Type-A + Type-B (24 vertices) form 24-cell, which IS the gauge manifold.

### Setup
From exp136-137:
- 24 gauge vertices form D_4 root system (24-cell)
- 24-cell is self-dual regular polytope
- Has 24 vertices, 96 edges, 96 faces, 24 cells
- Symmetry group: F_4 (Weyl group of D_4)

### Tasks
1. **Isolate 24-cell subgraph**
   - Remove Type-C vertices
   - Consider only effective 2-step gauge-gauge adjacency (via fermions)
   - This IS the 24-cell graph (from exp137)

2. **Gauge group from 24-cell symmetry**
   - Symmetry group of 24-cell: F_4 has order 1152
   - F_4 ⊃ SO(9)? (24-cell in some sense related to SO(9))
   - Check if F_4 contains SU(3) x SU(2) as subgroup

3. **Root system interpretation**
   - 24-cell vertices are D_4 roots (exactly!)
   - D_4 Lie algebra = so(8)
   - SU(3) x SU(2) ⊂ SO(8)? Check embedding.

4. **Fermions as rep of F_4**
   - 96 Type-C vertices
   - Is 96 an irrep dimension of F_4?
   - F_4 irreps: 1, 26, 52, 273, 324, 1053, ...
   - 96 is NOT an irrep of F_4 (problem!)

**Expected outcome**:
- 24-cell structure is beautiful but may not match SM gauge group
- D_4 = so(8) is close to SU(3) x SU(2) x U(1) but not exact
- Might need to break F_4 → SU(3) x SU(2) x U(1) (spontaneous?)

**Category**: PATTERN (24-cell = D_4) + SPECULATIVE (gauge group match)

---

## Which Experiment to Do FIRST?

### RECOMMENDATION: **EXP139 (Edge-Based Gauge Theory)**

**Reasons**:
1. **Standard framework**: Lattice gauge theory is well-established
2. **Natural fit**: Fermions on vertices, gauge on edges (conventional)
3. **Exact counts**: 96 A-C edges, 192 B-C edges (can be derived exactly)
4. **Plaquettes exist**: 1200 triangles (from exp136) provide field strength
5. **Minimal speculation**: Uses known structures (edges, faces)

**Concrete first steps**:
1. Classify all 720 edges by vertex types
2. Build edge-edge adjacency (edges sharing a vertex)
3. Enumerate triangles with each edge type (A-C, B-C, C-C)
4. Compute "plaquette matrix" (which edges form triangles together)
5. Check if 96 A-C edges have SU(3)-like structure

If EXP139 fails → try EXP140 (spectral).
If both fail → try EXP141 (dual lattice).
All three fail → maybe gauge structure is emergent (EXP143) or doesn't directly embed at all.

---

## Summary Table

| Exp | Priority | Framework | Key Idea | Main Risk |
|-----|----------|-----------|----------|-----------|
| 139 | HIGHEST | Lattice gauge | Gauge on edges | Edge count mismatch |
| 140 | HIGH | Spectral/NCG | Gauge as eigenmodes | No clean group match |
| 141 | MEDIUM | Dual lattice | Gauge on 120-cell | Factor 5 mismatch |
| 142 | MEDIUM | Holonomy | Gauge as connection | Freedom in choice |
| 143 | LOW | Emergent | Gauge from fermions | Strongly coupled |
| 144 | WILD | F_4 symmetry | 24-cell as gauge | Not SM group |

**Next action**: Implement EXP139 and see if edge-based gauge theory works.

If successful, this would be a MAJOR positive result to balance the EXP138 negative result.
