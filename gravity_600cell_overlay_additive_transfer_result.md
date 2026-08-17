# Additive common-refinement transfer is highly nonunique

Date: 2026-08-17

Only the targeted exact transfer verifier was run.  No gravity action and no
full suite were run.

## 1. Provenance

- prior-art gate: `c77a87c`
- preregistered protocol: `b0b30dc`
- verifier registered before calculation: `7061647`
- verifier:
  `reproducible/verify_gravity_600cell_overlay_additive_transfer.py`
- source overlay artifact SHA-256:
  `0dd03eed878f599463a44160484c74ddeaa0511fc70c8b2e77bc05a2f36dd3dc`
- result artifact:
  `reproducible/gravity_600cell_overlay_additive_transfer.json`
- result artifact SHA-256:
  `63fb1c57e99a0e10747a026fea7e6d9514fa519e56d132ff0f2fb91f195a5808`
- aggregation-matrix SHA-256:
  `cdd3dc6e63e6e8f0b9d2a9421f0e9c3607d0850dd76371a91f899a78aa871832`

The registered verifier passed `14/14` controls.  Rational ranks were checked
independently modulo `1000003`, `1000033` and `1000037`.

## 2. Complete hypothesis

The result concerns one additive scalar weight per full-dimensional fine
chamber.  For each of the 96 staircase four-simplices, the coarse value is the
sum of the weights of the contained fine chambers.  Positivity and invariance
under the full local `S4 x C2` are imposed.

This is not the nonlinear Regge action.  In Regge calculus curvature lives on
hinges and deficit angles couple neighbouring top cells.  All conclusions
below are restricted to piecewise-constant additive top-cell transfer.

## 3. Exact rank result

For the exact `96 x 148` aggregation matrix `R`,

```text
rank_Q(R)    = 15,
nullity_Q(R) = 133.
```

All three modular ranks are also 15.  Thus all 24 coarse triangulations
together observe only 15 independent combinations of 148 fine weights; 133
directions are invisible.

The small rank has an exact structural explanation.  For every subset `A`,
let `b_A(C)` indicate whether the chamber lies above the wall `h_A=0`, with
`b_empty=1` and `b_full=0`.  A staircase simplex is the interval between two
consecutive subsets `L subset U` on a maximal chain of `B4`, so

```text
R[(o,k),C] = b_L(C)-b_U(C).
```

The 96 labelled rows are repetitions of the 32 Hasse edges of `B4`: 24 edge
rows occur twice and eight boundary-edge rows occur six times.  They lie in
the span of the constant and 14 internal subset indicators, giving rank at
most 15.  Exact elimination shows that this upper bound is saturated.  This
is **DERIVED**, not a numerical pattern.

## 4. Symmetry does not select a lift

The 148 fine chambers form 14 `S4 x C2` orbits.  The 96 coarse labels form two
48-element orbits, corresponding to outer and inner staircase positions.  On
the invariant subspace,

```text
rank_Q(R_inv)    = 2,
nullity_Q(R_inv) = 12.
```

All three modular ranks are also 2.  Therefore symmetry retains only two
coarse totals and leaves 12 of the 14 invariant fine weights undetermined.

## 5. Constructive positive witness

The verifier gives a particularly small invariant kernel vector in the 14
fine-orbit coordinates:

```text
v = (-2,0,1,0,0,0,0,0,0,0,0,0,0,0).
```

The two affected fine orbits have sizes two and four, so the total variation
is `2*(-2)+4*(1)=0`.  With `epsilon=1/4`, define

```text
x_plus  = 1 + epsilon v,
x_minus = 1 - epsilon v.
```

After lifting orbit coordinates to all 148 chambers:

- both vectors are strictly positive;
- their minimum entries are respectively `1/2` and `3/4`;
- both are exactly `S4 x C2` invariant;
- both have total weight 148;
- both give exactly `(19,55,55,19)` on every staircase order.

They are distinct but have identical values in all 96 coarse cells.  Varying
the rational coefficient in a sufficiently small interval gives infinitely
many such positive invariant lifts.  The preregistered outcome is therefore

```text
POSITIVE_INVARIANT_ADDITIVE_TRANSFER_NONUNIQUE.
```

## 6. Hostile interpretation

- **DERIVED NEGATIVE:** incidence, conservation, positivity and full finite
  symmetry do not select additive fine weights.
- **DERIVED:** the fine-to-coarse aggregation map itself is canonical; its
  inverse/lift is not.
- **STRUCTURAL:** the common overlay solved schedule choice at the carrier
  level but relocated the freedom to the dynamical transfer law.
- **NOT DERIVED:** that no canonical Regge--dust action exists.  Hinge-local,
  metric-volume, Galerkin, nonlocal and dynamically perfect constructions are
  outside this theorem.

Adding Euclidean/Lorentzian volumes, an `L2` norm or a path-integral measure
could select a lift, but each is extra structure and must be derived or openly
declared.  Conservation alone cannot justify it.

## 7. Status ledger

| Claim | Status |
|---|---|
| Common-refinement aggregation `R` is canonical | **DERIVED COMBINATORIAL** |
| `rank(R)=15`, `nullity(R)=133` | **DERIVED EXACT** |
| Boolean-lattice edge-difference factorization | **DERIVED EXACT** |
| Invariant rank 2 and nullity 12 | **DERIVED EXACT** |
| Positive invariant nonuniqueness | **DERIVED EXACT / explicit witness** |
| Carrier alone selects an additive action density | **REFUTED** |
| Carrier alone selects nonlinear Regge dynamics | **OPEN, not tested** |
| Metric-volume or Galerkin transfer is physically selected | **OPEN** |
| Perfect 4D Regge--dust action exists on this overlay | **OPEN** |

## 8. Consequence and next step

The correct continuation is not to fit one of the 12 invariant weights.  It is
to derive a new action directly from geometry.  The minimal deterministic
route is:

1. enumerate the complete overlay face poset;
2. take its functorial barycentric subdivision;
3. test exact Lorentzian nondegeneracy at the certified slab geometry;
4. only then evaluate the ordinary Regge--dust action on that new carrier.

This would define a third, geometry-derived discretization rather than infer
fine dynamics from the old schedules.  Its agreement with either old schedule
or with a continuum limit would remain a subsequent falsifier.

