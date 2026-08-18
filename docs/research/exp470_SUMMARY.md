# exp470: E8 Rotations and 600-Cell Decompositions -- SUMMARY

## Status: IN PAPER (Remark in Section 7.7, full S7.5 in Supplementary)
## Classification: STRUCTURAL THEOREM about E8 (DERIVED)
## Quantum Gravity: NEGATIVE (mechanism doesn't work)

---

## The Question

How many ways can the 240 E8 roots be projected to give a 600-cell?
Can different projections parameterize quantum gravity configurations?

## The Answer

### Task 1: Number of 600-cell decompositions

**ANSWER: Exactly 40 distinct decompositions.** 40 = rank(E8) * a1 = 8 * 5.

The Coxeter projection R^8 -> R^4 (onto the H4 eigenspace with exponents {1,11,19,29}) maps 240 E8 roots to two concentric 600-cells:

- **Inner 600-cell**: 120 roots at projected norm^2 = 1 - 1/sqrt(a1)
- **Outer 600-cell**: 120 roots at projected norm^2 = 1 + 1/sqrt(a1)
- **Norm ratio** = phi (golden ratio, EXACT)
- **Product of norm^2** = 4/a1 = 4/5

All 8! = 40320 orderings of simple reflections yield exactly 40 distinct partitions of the 240 roots into inner/outer sets.

### Task 2: Graph isomorphism

**ALL 40 decompositions are isomorphic** as root subsets:
- Identical Gram matrix eigenvalues
- Same E8-adjacency spectra (degree 32 for all 120 roots in each half)
- Same inner product distribution: -2(60), -1(1920), 0(3240), +1(1920)

### Task 3: Grassmannian structure

**EXACT THEOREM** (verified for all 780 pairs):

> **overlap(D_i, D_j) = h(E8) * sum_k cos^2(theta_k)**
>
> where theta_1,...,theta_4 are the principal angles between the two 4D subspaces, and h(E8) = 30 is the Coxeter number.

The sum of squared cosines takes exactly 4 discrete values:
- sum cos^2 = 2m/a1 for m = a1, a1+1, ..., rank(E8)  (i.e., m = 5, 6, 7, 8)

This gives overlaps: 60, 72, 84, 96 (always multiples of 12 = vertex degree).

| m | sum cos^2 | overlap | # exchanged | # pairs |
|---|-----------|---------|-------------|---------|
| 5 | 10/5 = 2  | 60      | 60 = 5*12   | 176     |
| 6 | 12/5 = 2.4| 72      | 48 = 4*12   | 296     |
| 7 | 14/5 = 2.8| 84      | 36 = 3*12   | 184     |
| 8 | 16/5 = 3.2| 96      | 24 = 2*12   | 124     |

Number of overlap levels = rank(E8) - a1 + 1 = 4 = dim(projected space).

### Grassmannian angles

The "clean" principal angle cosines are:

**cos(theta) = sqrt(k/a1) for k = 0, 1, 2, 3, 4, 5**

Notably:
- sqrt(1/a1) = 1/sqrt(5) = 0.4472
- **sqrt(2/a1) = c_bare** (the bare self-energy coefficient from the mass formula!)
- sqrt(3/a1) = sqrt(3/5) = 0.7746
- sqrt(4/a1) = 2/sqrt(5) = 0.8944

Some pairs also have "mixed" cosines that are NOT of the sqrt(k/a1) form, but the sum of squared cosines is ALWAYS 2m/a1 with integer m.

Product of all 4 SVs for special patterns:
- a1^(-1/2), a1^(-1), a1^(-3/2), a1^(-2), 0

### Task 4: Path integral

Since all 40 decompositions are isomorphic, Z = 40 * exp(-S_0) is trivial.
The 40-element discrete structure is too small for continuum gravity.
**STATUS: NEGATIVE.**

### Task 5: Regge calculus connection

The 40 decompositions define 40 distinct 4D subspaces, but all give the same
600-cell geometry. The Grassmannian angles between subspaces are determined
by a1 = 5, but this doesn't connect to Regge edge lengths.
**STATUS: NEGATIVE.**

---

## Key Structural Results (potential paper material)

1. **40 = rank(E8) * a1**: the number of Coxeter 600-cell decompositions of E8.
2. **overlap = h * Tr(cos^2 theta)**: exact relation between combinatorial overlap and Grassmannian metric. This appears to be a new result in root system theory.
3. **Norm ratio = phi**: the two concentric 600-cells have radius ratio phi.
4. **n1^2 * n2^2 = 4/a1**: product of squared radii.
5. **cos(theta) = sqrt(k/a1)**: principal angles in the Grassmannian.
6. **c_bare = sqrt(2/a1) appears as a Grassmannian angle**: potential structural origin for the mass formula coefficient.

---

## Framework Connections

| Quantity | Value | Expression |
|----------|-------|------------|
| # decompositions | 40 | rank(E8) * a1 |
| # overlap levels | 4 | rank(E8) - a1 + 1 |
| overlap formula | 30 * sum cos^2 | h(E8) * Tr(cos^2) |
| inner norm^2 | 0.5528 | 1 - 1/sqrt(a1) |
| outer norm^2 | 1.4472 | 1 + 1/sqrt(a1) |
| norm ratio | phi | (1+sqrt(a1))/2 |
| product | 0.8 | 4/a1 |
| Grassmannian angle | 50.77 deg | arccos(sqrt(2/a1)) = arccos(c_bare) |

---

## c_bare = sqrt(2/a1) Grassmannian Analysis (exp470h-i)

Three convergent routes to c_bare = sqrt(2/a1):
1. **McKay CG vertex**: ||v_k|| = sqrt(C_2/d_k) = sqrt(2/5) [amplitude norm]
2. **S-matrix prefactor**: sqrt(2/(k+2)) = sqrt(2/5) [modular unitarity]
3. **Grassmannian angle**: cos(theta) = sqrt(2/5) [principal angle in Gr(4,8)]

All three trace to the **Plancherel measure** of the Verlinde ring Z[phi].
Classification: STRUCTURAL CONVERGENCE (not independent, but new geometric content).
The Grassmannian provides an inherent justification for sqrt (amplitude, not probability).

Proof: For the 132 pairs with shared Coxeter plane, conjugate symmetry of
exponents {11,19} forces cos^2_3 = cos^2_4 = (14/a1 - 2)/2 = 2/a1. QED.

D6 -> H3 comparison: FAILED (5 shells). E8 is special (all exponents distinct).

## Scripts
- `exp470_quantum_gravity_e8.py` - Initial construction (superseded)
- `exp470b_e8_analysis.py` - Wrong simple roots (superseded)
- `exp470c_correct_e8.py` - First correct results
- `exp470d_counting.py` - 40 decompositions found
- `exp470e_final_analysis.py` - Full structural analysis
- `exp470f_verification.py` - Verification of exact formulas
- `exp470g_deep_analysis.py` - Weyl group counting, D6 test, h=a1*b1
- `exp470h_cbare_grassmannian.py` - c_bare as Grassmannian angle (full analysis)
- `exp470i_cbare_connection.py` - Three routes comparison, theorem proof
