# Second finite-height full-boundary tangent: consolidated result

Date: 2026-08-22

## Verdict

**DERIVED COMPUTATIONAL / ADVERSARIALLY REPLICATED, WITH THE SCOPE BELOW.**

On the frozen branch-B history at incoming representative `v=3/2`, both
staircase parities have regular complete pre-Legendre systems on the first and
second finite-height slabs. Their directly solved 1440-dimensional boundary
tangents are canonical within the propagated numerical bounds. The normalized
and direct-physical second-slab actions realize the exact degree-two scale
lift, and all four parity-ordered two-step products agree under the frozen
10/100 classifier.

This is a two-step linearized canonical response on the unreduced boundary
phase carrier. It is not yet a physical perturbation spectrum.

Primary artifact commit: `2a80690` (`31/31`).  
Adversarial protocol commit: `eceda30`.  
Adversarial implementation commit: `977f3ba`.  
First failed run preserved at commit: `94d6b3b` (`27/28`).  
Correction preregistration commit: `88e3ef3`.  
Corrected implementation commit: `8668c0e`.  
Accepted adversarial artifact SHA-256:
`1355f8cf339d18c1cf2855ecb1228e97e868d73f7a1ef739e4c11ce9521fcd4b`.  
Accepted adversarial verifier SHA-256:
`a0793b5cbe865ed1b96c7ab1411474057cc275d4b3bcb158f01d2f92680b2d9d`.

Only the targeted adversarial verifier was run. It returned `28/28 PASS` and

```text
TWO_STEP_FULL_BOUNDARY_TANGENT_SCHEDULE_ROBUST_ADVERSARIALLY_REPLICATED
```

## Frozen hypotheses

The claim requires all of the following:

1. the regular 600-cell one-slab carrier with 720 old, 840 internal and 720
   new edge variables;
2. the frozen zero-`Lambda` Lorentzian Regge-plus-conserved-dust action;
3. the accepted branch-B history reconstructed from the fixed brackets
   `(9,10)` and `(31,32)`;
4. the two frozen even/odd staircase schedules;
5. centered logarithmic differentiation at `1e-18`, `5e-19`, `2.5e-19` and
   `1.25e-19`, followed by the three adjacent Richardson levels;
6. the implicit boundary Legendre map with the actual `K_NO` term and the
   lexicographic 720-edge old-to-new identification;
7. the preregistered 10/100 uncertainty classifier.

No other incoming state, triangulation family, matter model or continuum
limit is covered.

## Mechanically different replication

The primary route used 95-by-95 orbit-convolution kernels and seven minimal
binary-tetrahedral sectors. The adversarial route did not use these objects
for its scientific decision. At 120 decimal digits it independently:

- assembled six complete real 2280-by-2280 Hessian families from all 2400
  four-simplices and 6240 triangles;
- checked raw reciprocity before the one licensed symmetrization;
- formed eighteen complete 1560-by-1560 pre-Legendre systems;
- compared both `gesvd` and `gesdd` at every Richardson level;
- solved eighteen systems with 1440 right-hand sides;
- checked the three real symplectic block identities on each map and product;
- froze the dense scale, rank, canonicality and schedule labels before opening
  a primary tangent entry.

The independently recomputed first slab reproduced all six previously
accepted dense tangent hashes exactly. This was the positive known-answer
control. Omitting `K_NO`, using the identity or `r1` momentum lift, and adding
a common `1e-3` product corruption were the negative controls; all were
refuted.

## Quantitative result

- Both parities returned `SCALE_LIFT_CONFIRMED` on all four raw and three
  Richardson Hessian comparisons. The largest normalized distance was
  `3.19e-15`, at less than `3.51e-5` of its uncertainty.
- The normalized and directly assembled physical tangents agreed with
  `D_c T D_c^-1`, `D_c=diag(I,cI)` and `c=r1^2`; the largest normalized
  distance was `2.99e-12`.
- The second-slab even/odd distance was at most `6.64e-13` against uncertainty
  `9.37e-11`.
- All six pairs among the four two-step products agreed. Their largest
  normalized distance was `5.63e-12`.
- The `1e-3` hostile corruption had normalized distance `5.55e-9` against
  uncertainty `5.14e-12` and was therefore refuted above the frozen
  100-times boundary.
- After the dense labels were frozen, all 14 first-slab, 14 second-slab and 28
  product minimal blocks agreed entrywise with the primary ball archives:
  `56/56 AGREES`. The largest distance-to-uncertainty ratio was `4.76`, below
  the frozen agreement boundary of 10.

The normalized smallest pre-Legendre singular values were approximately
`3.34e-4` to `4.09e-4`; their strict numerical gates were `4.56e-9` to
`9.13e-8`. Thus regularity is not a marginal threshold event in this
calculation.

## Preserved failed run and correction

The first adversarial execution returned `27/28` and
`SECOND_FULL_TANGENT_DENSE_CONTROL_FAILED`. It is preserved rather than
silently overwritten. The delayed closure had constructed its coordinate
bases at 120 digits, whereas the primary archive had constructed them at 180.
A target-free 24-element carrier diagnostic found order-one unitary frame
changes in every 2D/3D sector and only `1e-17` changes in the 1D sectors. This
exactly explained why all scalar blocks agreed while every non-scalar block
was initially refuted. The corrected closure reconstructs the disclosed
primary coordinate convention at 180 digits only after the dense firewall.

The same failed run exposed an overly loose product error term. For already
solved maps, the first-order perturbation identity

```text
delta(T2 T1) = delta(T2) T1 + T2 delta(T1) + delta(T2) delta(T1)
```

gives additive, not multiplicative, condition amplification. The corrected
term `kappa1+kappa2+u*kappa1*kappa2` was preregistered before rerunning the
unchanged `1e-3` control. These were implementation corrections, not changes
to a scientific threshold or target.

## What this establishes and what it does not

**DERIVED COMPUTATIONAL:** the frozen nonlinear history admits a regular,
canonical, schedule-robust two-slab linear response on the complete boundary
carrier, and the action's exact degree-two scale covariance is realized by
the directly assembled physical systems.

**STRUCTURAL LIMIT:** canonical here means symplectic for the unreduced
boundary Legendre data. Constraint, pseudo-constraint, gauge and dust
directions have not been separated. The 1440 coordinates are not 1440
physical degrees of freedom.

**STRUCTURAL LIMIT:** robustness was tested for the two derived staircase
parities, not for every triangulation or refinement. It was tested at the
representative branch-B history, not on the complete incoming basin.

**OPEN:** the external novelty of the branch-specific coefficients remains
open. Two-step/Jacobi machinery itself already existed elsewhere in this
repository on the old near-static background, so this result is not a claim
to a new canonical formalism.

**NOT DERIVED:** no tangent eigenvalue or singular-value spectrum was opened;
there is no graviton, wave equation, limiting speed, absolute tick, `G`,
Planck scale, particle mass or Standard-Model result here.

## Next gate

The single next gate is `P-CONSTRAINT-QUOTIENT`: derive from the action the
constraint and pseudo-constraint covectors, freeze a quotient criterion that
does not choose a numerical threshold after seeing a spectrum, and prove that
the two-step maps descend to the same nondegenerate physical phase space.
Only after that gate passes is a finite-height Jacobi/wave analysis physically
interpretable.
