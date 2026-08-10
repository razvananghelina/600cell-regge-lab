# B1 is false: an exact noncommutative all-gate chamber algebra

Date: 2026-08-10

## Headline verdict

**DERIVED REFUTATION.** On the fixed chamber carrier `H=C^120`, with the
derived three-regular chamber adjacency `D`, orientation grading `gamma`, and
geometric reflection `J`, there is a faithful unital representation of

`A=M2(C) direct-sum C direct-sum C direct-sum C`

which is noncommutative and satisfies every B1 hypothesis:

1. order zero;
2. first order;
3. `[gamma,A]=0`;
4. nonzero represented inner one-forms;
5. a nondegenerate KO6 intersection form;
6. metric-dimension-zero orientability;
7. connectedness, `[D,pi(a)]=0` only for `a` in `C 1_A`.

Thus B1 is false. This is an exact certificate on the unchanged chamber
matrices. There is no fitted `D`, isospectral replacement, numerical rank
tolerance, or assumption that a Krajewski block mask is the full first-order
condition.

The intersection determinant is `1444=38^2`, so it is nondegenerate but not
unimodular. The stronger conjecture obtained by replacing “nondegenerate” by
“integrally unimodular” remains **OPEN**.

## Complete algebraic convention

Here `A` is a finite-dimensional complex star algebra, the repository's
convention for a represented complex star subalgebra of `M120(C)`. The
antiunitary is `J=P K`, where `P` is the exact reflection permutation and `K`
is coordinate conjugation. Hence

`pi(b)^0 = J pi(b)^* J^-1 = P pi(b)^T P`.

The verifier checks complex linearity, multiplication, and the star
operation, not merely real matrix units.

## Attack on the finite-type framing

### The general `k<=7` bound is false

Total faithfulness says only

`ker(pi_+) intersect ker(pi_-)=0`.

It does not make either sheet representation faithful. The one-cell
bimodule `(M60(C),C)` gives a faithful unital two-sheet representation of
`M60(C)+C`: `M60` acts on `H+`, the scalar summand acts on `H-`, and order
zero is automatic. Thus order zero, grading, faithfulness, and unitality
permit `k=60`. **DERIVED NEGATIVE.**

The argument `k*m=60`, `M_k subset M_m`, hence `k^2<=60`, is valid only when
the same summand fills both sheets in the stated way. Under that complete
extra hypothesis it remains **DERIVED**.

### The corrected bound is `k<=30` with orientability and first order

Order zero decomposes `H+` into joint `A-A^op` cells

`H_ij=C^(k_i) tensor C^(k_j) tensor C^(mu_ij)`.

Metric-dimension-zero orientability forbids diagonal cells and reverse-paired
cells. First order then makes the compression of

`S=(D J)|H+`

to each occupied cell zero: a block from `(i,j)` to its own `J`-paired copy
would change both Krajewski indices. Since fixed `S` is invertible, its
restriction injects a cell of dimension `w_ij` into its complement of
dimension `60-w_ij`. Therefore `w_ij<=30`.

Every simple summand is incident to an occupied cell by total faithfulness,
so `k_i*k_j<=w_ij<=30` for some `j`, and hence `k_i<=30`. **DERIVED.** This
is finite, but it does not make a Wedderburn-type census decisive.

### The even-summand condition stands

For `A=direct_sum_i M_(k_i)(C)`, full `K0(A)` has one generator per simple
summand. The KO6 intersection form is antisymmetric. An odd-order
antisymmetric matrix over characteristic zero is singular, so nondegeneracy
requires an even number of summands. **DERIVED necessary condition**, not a
sufficiency result.

### The decisive defect: the gates are not type-level

Orientability is determined by the directed bimodule support `mu_ij`, not by
the tuple `(k_i)`. Connectedness depends further on the fixed-`D` blocks and
their embedding. The repository now contains two representations of exactly
the same type `M2(C)+C^3` on the same carrier:

- the older witness passes order zero, first order, and nondegenerate PD but
  fails orientability and connectedness;
- the witness below passes both additional gates.

Thus neither axiom is a Wedderburn-type invariant. A type-only enumeration
cannot decide B1. One would need signed Krajewski multiplicities and then a
fixed-`D` embedding calculation with continuous unitary data in general.
The proposed enumeration strategy is **DERIVED INSUFFICIENT**, even though
the corrected bound makes the abstract type set finite.

## Exact counterexample

### Fixed graph and cells

Order `H-` by the reflection of sorted `H+`. Then `S=(D J)|H+` is an exact
symmetric invertible, zero-diagonal, three-regular integer matrix on 60
vertices, with 90 edges.

Index the summands by sizes `(k0,k1,k2,k3)=(2,1,1,1)`. The positive sheet has
four cells:

| cell | multiplicity | dimension |
|---|---:|---:|
| `(0,1)` | `2` | `4` |
| `(1,2)` | `25` | `25` |
| `(3,1)` | `12` | `12` |
| `(2,3)` | `19` | `19` |

Their dimensions sum to 60. The `M2` node is a pure source.

An exact colouring of the 60 vertices has capacities `(4,25,12,19)` and
sends every edge of `S` into one of

`(cell0,cell1), (cell1,cell2), (cell1,cell3), (cell2,cell3)`.

These are composable Krajewski pairs. More importantly, every shared index
is a scalar node. Arbitrary entries in these blocks therefore satisfy the
full tensor first-order condition. The verifier nevertheless checks the full
condition on all seven complex basis elements.

The 60-entry certificate is frozen in
`reproducible/verify_chamber_b1_counterexample.py` and checked against all 90
edges. It was found by an exact finite constraint search; the proof is the
independently checked certificate, not a statistical claim about a search.

### Representation and order conditions

On a positive cell `(i,j)`, `A` acts on the left `C^(k_i)` factor. On its
reflected negative cell it acts on the `C^(k_j)` factor. The four `M2` matrix
units and three scalar units form a complex basis. Exactly:

- the action is complex-linear, multiplicative, star preserving, faithful,
  unital, and noncommutative;
- every represented element commutes with `gamma`;
- all `7*7` order-zero commutators vanish entrywise;
- all `7*7` first-order double commutators vanish entrywise;
- at least one `[D,pi(a)]` is nonzero.

No center-only test is used for either order condition.

### Orientability

Let `p_i` be the central summand projectors and `p_j^0` their opposite
actions. The explicit Hochschild zero-cycle

`sum_(i,j in cells) (p_i p_j^0-p_j p_i^0)`

is represented by exactly `gamma`: `+I` on all 60 positive chambers and
`-I` on all 60 negative chambers. Metric-dimension-zero orientability is
**DERIVED**.

### Intersection form

Using a minimal rank-one projection for `M2` and the scalar projectors gives

```text
Cap = [[  0,   2,   0,   0],
       [ -2,   0,  25, -12],
       [  0, -25,   0,  19],
       [  0,  12, -19,   0]].
```

It is antisymmetric, has rank four, Pfaffian `38`, and determinant `1444`.
Poincare nondegeneracy is **DERIVED**.

### Connectedness

Flatten the seven matrices `[D,pi(e_q)]` into the commutator map. Its exact
rational rank is six. Its complex kernel is therefore one-dimensional and
is spanned by the represented identity:

`[D,pi(a)]=0 iff a in C 1_A`.

Connectedness is **DERIVED**, not inferred only from a quotient graph.

## Correction to the previous numerical detector

`verify_math_to_physics_bridge.py` previously called its isospectral
block-mask projection a “first-order detector.” That was too strong. A
shared-index mask is only the support part of first order. If the shared
factor is `M2` or `M3`, the block must also be an identity intertwiner on
that factor. The script did not test those equations.

Its exact enumerations remain valid as cell-support screens, and its
machine-precision output is a **PATTERN for the coarser mask only**. It is not
evidence for a full first-order physical algebra. The new B1 witness avoids
the issue structurally and checks the full equations anyway.

## Status ledger

### DERIVED

- an exact noncommutative all-gate witness on unchanged `D,J,gamma`;
- B1 is false;
- the general order-zero `k<=7` inference is false;
- the corrected `k<=30` bound under orientability, first order, and
  invertible fixed `S`;
- odd simple-summand count is incompatible with nondegenerate KO6 pairing;
- orientability and connectedness are not Wedderburn-type invariants;
- the old numerical mask omitted full tensor first order.

### STRUCTURAL

- the choice of `M2(C)+C^3`, its cells, and the displayed colouring;
- treating its non-abelian unitary factor as a candidate gauge ingredient.

### OPEN

- B1 strengthened to require an integrally unimodular form;
- a geometry-selected rather than existence-selected noncommutative algebra;
- a color `M3(C)` sector, Standard-Model matter, hypercharge, anomalies,
  Higgs/Yukawa data, and continuum physics.

### REJECTED

- B1 as stated;
- exhaustive type-only enumeration as a decision procedure;
- the claim that a Krajewski block mask alone verifies full first order.

## Physics boundary

The chamber axioms do not force every admissible gauge algebra to be a torus.
But this does not open the Standard-Model gate. The witness was selected to
falsify a universal theorem, not selected uniquely by the geometry; it has an
`M2` factor but no derived `M3` color factor or matter representation. It
proves possibility, not physical selection.

Registered verifier:

`/home/razvan/science/.venv/bin/python reproducible/verify_chamber_b1_counterexample.py`

Focused result: `17/17 PASS`. No PDF build was attempted.

Implementation and registration commit: `4c05edb`.

Final unique registered suite, run on that commit from the repository root:

`/home/razvan/science/.venv/bin/python reproducible/run_all.py`

Result: `80/80 scripts completed successfully` in `803.2 s`, with process
exit `0`. The bidirectional coverage and duplicate-registration guards
passed. No PDF build was attempted.
