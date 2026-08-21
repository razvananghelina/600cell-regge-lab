# Protocol: complete nonhomogeneous internal Hessian of the refined slab

Date: 2026-08-21

This protocol is frozen before any complete local Hessian is assembled or any
nonhomogeneous eigenvalue is observed.

Prior-art gate commit: `d4dc6c7`.

## 1. Frozen inputs

```text
commons/cell600.py
  ea5bce4b6c52e0834539ca4b1df9c6a67a3a5ed4da32f4e0298a493fc5315c7f

reproducible/verify_gravity_600cell_refined_h4_stationary_fill.py
  89aab727792e20a81e7577e0425f8fa4b1e84e2a7ae66caa9e79a4aebf3581e7

reproducible/gravity_600cell_refined_local_curvature_mass.json
  180010a79177ba16620ebea9847443c57a7a6d2d8a3df71ad6ecb83f454ef091

reproducible/verify_gravity_600cell_refined_local_curvature_mass.py
  c54f17708a2678b925cfce96fcfc7d6baaeeaf0577bedbf22b5d0435c069fae6

reproducible/gravity_600cell_refined_h4_null_coupling_adversarial.json
  5c1f596958f9d878c8d9d3ccb6ecc8359f72164e8f36dd9930fb71ddc1351ce9

docs/gravity/gravity_600cell_refined_h4_null_coupling_result.md
  660a3707f24f44d0393e6a1804e407fa45aa4782a98960438a296da50c35825a

docs/gravity/gravity_600cell_refined_nonhomogeneous_internal_hessian_prior_art.md
  6e766435caa20404fcc9d30403cc27a969b827ab4d3fc232ce3ee58b1f90cd38
```

The implementation must stop before scientific assembly if any digest differs.

## 2. Geometry and coordinates

Rebuild the 600-cell tetrahedra from adjacency, barycentrically subdivide
them, and generate all 24 rank-colour staircase slabs directly.  For every
schedule require exactly

```text
spatial f-vector         (2640,17040,28800,14400)
slab pentachora          57600
boundary edges           34080
internal cross edges     17040
internal vertical edges   2640
all internal edges       19680
all slab edges           53760.
```

Order boundary and internal edges lexicographically by their labelled vertex
pairs.  This ordering is bookkeeping only; no schedule identification other
than the explicit layer reversal is allowed.

At the fixed product seed, assign signed squared lengths from the exact
rank-pair geometry and `tau0=0.0102`.  The local variables are logarithms of
the absolute signed squared lengths.  A vertical edge is negative and every
cross diagonal is positive on the accepted branch.

Use the per-vertex masses from the frozen curvature-mass artifact, indexed by
the barycentric rank.  No density variable is differentiated.

## 3. Sparse local Hessian

For each distinct local pentachoron signed-length pattern, compute all ten
dihedral angles and their `10 x 10` first derivative matrix.  Use centred
differences at

```text
h, h/2, h/4 with h=1e-10
```

and fourth-order Richardson extrapolation.  Repeat at 100 and 140 decimal
digits.  The selected stencil is the finest Richardson value at 140 digits;
the local derivative envelope is

```text
100 * (successive-Richardson difference
       + 100-vs-140 digit difference
       + 1e-70).
```

No observed spectrum enters this envelope.  Every displaced simplex must
retain one Lorentzian Gram negative direction, a nonzero logarithm argument,
and smooth angle continuation.

Assemble the internal block with the Regge--Schlaefli formula

```text
C = -i sum_h [epsilon_h Hess(A_h)
              + grad(A_h) tensor grad(epsilon_h)] + Hess(S_dust).
```

Area gradients and Hessians are analytic in the three signed squared hinge
lengths.  The dust contribution on vertical edge `v` is the diagonal value

```text
-2*pi*m_v*tau0.
```

Accumulate a parallel sparse absolute-error matrix from the local derivative
envelopes.  Add a standard forward binary64 summation envelope
`gamma_(k+32)` times the absolute term sum for every entry.  The decisive
operator-error bound is the maximum absolute row sum of the symmetric error
matrix.  The raw matrix may be replaced by `(C+C^T)/2` only after its
antisymmetric row norm is inside that bound.

## 4. Required controls

For every schedule:

1. every individual internal gradient must be zero within its assembled
   forward-error envelope;
2. the complete action and the assembled gradient/Hessian must be real within
   their branch and arithmetic envelopes;
3. with `U_s` the ten-column indicator matrix of the six cross-rank and four
   vertical-rank edge orbits, `U_s^T C_s U_s` must agree entrywise with a
   freshly differentiated `10 x 10` aggregate action block built from the
   frozen action source;
4. the analytic product tangent

   ```text
   (n_s)_vertical = 1,
   (n_s)_cross(ra,rb) = -rho0/x_cross(ra,rb)
   ```

   normalized in the Euclidean local coordinate norm, must satisfy `C_s n_s`
   inside `100` times the combined operator and multiplication envelope;
5. the explicitly rebuilt reverse schedule must be congruent under layer
   reversal to the forward matrix inside the combined sparse envelope;
6. changing one selected local angle-derivative entry by `1e-4` before
   assembly must be detected by either the aggregate pullback or null-line
   control.  The corrupted matrix is never used scientifically.

The aggregate comparison is a different incidence path from the local
assembly: the aggregate evaluator groups simplices and hinges by rank-layer
state before differentiation; the local route differentiates stencils and
then sums actual labelled incidences.

## 5. Complete-kernel decision

For each time-reversal pair form, for one representative, the bordered sparse
matrix

```text
K_s = [[C_s,n_s],[n_s^T,0]].
```

Run two shift-invert symmetric eigensolves with different deterministic start
vectors, tolerances `1e-10` and `1e-12`, and request the eight eigenvalues
nearest zero.  Also perform sparse LU solves of `K_s` on the all-ones vector
and seven fixed signed integer probe vectors.  Record Ritz residuals, solve
residuals, factorization permutations and the two smallest-spectrum lists.

Define

```text
uncertainty = operator row-error bound
              + maximum Ritz residual
              + maximum normalized sparse-solve residual * ||K_s||_inf.
```

The bordered matrix is `NUMERICALLY_RESOLVED_NONSINGULAR` only if:

- both eigensolves converge and pair their eight nearest-zero values within
  `100*uncertainty`;
- every solve residual is at most `100*uncertainty/max(1,||K_s||_inf)`;
- the smallest absolute returned eigenvalue exceeds `100*uncertainty`;
- the reverse schedule has the same list within the combined uncertainty.

This is a reproducible computational certificate, not interval arithmetic.
If shift-invert convergence, spectrum coverage or separation is doubtful, the
outcome is **OPEN**, not nonsingular.

## 6. Frozen verdict hierarchy

```text
LOCAL_EXTENSION_INVALID
```

if a provenance, geometry, stationarity, branch, reciprocity, aggregate,
null-line, time-reversal or corruption control fails.

```text
COMPLETE_INTERNAL_KERNEL_NUMERICALLY_OPEN
```

if the local construction passes but at least one bordered system is not
resolved under the frozen numerical gates.

```text
COMPLETE_INTERNAL_KERNEL_IS_PRODUCT_DURATION_LINE
```

only if all 12 representatives and their 12 explicit reverses pass and every
bordered matrix is resolved nonsingular.  The conclusion is then **DERIVED
COMPUTATIONAL** pending the mechanically different adversarial replication
required by Rule 4.

No eigenvector may be compared with a graph Laplacian, continuum harmonic,
graviton polarization or desired dispersion relation in this mission.  No
schedule is averaged or selected.  No full verifier suite and no deferred
nonlinear root census will run.

## 7. Artifact and registration

Write

```text
reproducible/gravity_600cell_refined_nonhomogeneous_internal_hessian.json
```

containing the input digests, complete counts, all 24 schedule diagnostics,
all 12 bordered spectra, error envelopes, control outcomes, CSR matrix hashes
and the frozen verdict.  Do not store dense matrices.  Register exactly one
new verifier in `reproducible/run_all.py`, strengthen no unrelated test, and
run only that verifier plus the static registry audit.
