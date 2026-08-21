# Protocol: constrained refined H4 linearized boundary response

Date: 2026-08-21

Prior-art commit: `8ecbd2a`.

This protocol is frozen before constructing or inspecting any constrained
effective response matrix.

## 1. Frozen inputs

Require these exact files and hashes:

```text
reproducible/verify_gravity_600cell_refined_h4_stationary_fill.py
  89aab727792e20a81e7577e0425f8fa4b1e84e2a7ae66caa9e79a4aebf3581e7
reproducible/gravity_600cell_refined_local_curvature_mass.json
  180010a79177ba16620ebea9847443c57a7a6d2d8a3df71ad6ecb83f454ef091
reproducible/gravity_600cell_refined_local_curvature_mass_adversarial.json
  c59890d12bf929c4677dffed1b932ad8c05ab0ac00980be15ba780e62744c28e
reproducible/gravity_600cell_refined_boundary_cotangent.json
  4e7bf0beb0327a3ee1bddbec13126fbef99380970e62cecf74eb24ce8d6dafaa
reproducible/gravity_600cell_refined_effective_h4_hessian.json
  56e08db9a840b95e686fadb2763e89400b09220e88b80e9d35c17c1e73eef0a3
reproducible/gravity_600cell_refined_h4_null_coupling.json
  6b6fbd95b07f365b3fcac332fa3546021e8d756a510af0184bc974e52d5efa79
reproducible/gravity_600cell_refined_h4_null_coupling_adversarial.json
  5c1f596958f9d878c8d9d3ccb6ecc8359f72164e8f36dd9930fb71ddc1351ce9
docs/gravity/gravity_600cell_refined_h4_null_coupling_result.md
  660a3707f24f44d0393e6a1804e407fa45aa4782a98960438a296da50c35825a
docs/gravity/gravity_600cell_refined_h4_constrained_hessian_prior_art.md
  222f31862e911e03a1a7740696618948e370e43164812120120d85e834f0f639
commons/cell600.py
  ea5bce4b6c52e0834539ca4b1df9c6a67a3a5ed4da32f4e0298a493fc5315c7f
```

Require every accepted upstream outcome, all 24 internal inertias `(9,1,0)`,
the analytic product null line and the accepted nonzero rank-one compatibility
row.  Load only function definitions from the frozen action source.  Do not
execute its top-level P1 calculation or import either Hessian verifier.

## 2. Coordinates, null line and compatibility space

Use the fixed total-orbit log squared-edge ordering

```text
B=(old_01,old_02,old_03,old_12,old_13,old_23,
   new_01,new_02,new_03,new_12,new_13,new_23),
I=(cross_01,cross_02,cross_03,cross_12,cross_13,cross_23,
   rho_0,rho_1,rho_2,rho_3).
```

Rebuild the analytic internal tangent

```text
n_cross,rs=-tau0^2/q_cross,rs,
n_rho,r=1.                                         (1)
```

Use the twelve-component adversarial row already frozen in the accepted
artifact as `c`; do not refit it from the new Hessians.  Recheck
`H_ii n=0` and `H_bi n=c` for all schedules using

```text
e_n=100*10*e_H*max(1,||n||max)+1e-65.             (2a)
```

Require both maximum component errors below `e_n` and require
`||c||max>10^6 e_n`; the latter prevents a numerically unresolved coupling
from licensing the reduction.

For a nonzero vector `a` and pivot `p`, define the algebraic kernel basis

```text
E(a,p)_j = e_j-(a_j/a_p)e_p,  j != p.             (2)
```

Freeze

```text
P     =E(c,3),   P_alt=E(c,9),
Q     =E(n,9),   Q_alt=E(n,6).
```

Thus `P,P_alt` span `ker(c^T)` and `Q,Q_alt` span `ker(n^T)`.  Require their
dimensions and annihilation identities explicitly.  Let `T` be the exact
change of boundary basis fixed by `P_alt=P T`; construct it from the eleven
non-pivot rows of `P_alt`, where the corresponding rows of `P` form the
identity.

Let `R` exchange the six old and six new boundary coordinates.  Require
`c^T R P=0`; no rank-colour relabelling is permitted.

## 3. High-precision Hessian ladders

At the exact static product `tau0=0.0102`, set the legacy P1 mass to zero and
add the selected matter Hessian only through

```text
H_dust[rho_r,rho_r]=-2*pi*m_r*tau0.               (3)
```

Before differentiation, recheck all ten total internal gradients below
`1e-60` for all 24 schedules.

At `100` and `140` decimal digits, differentiate the analytic 22-component
gravitational log-gradient by centred log-coordinate differences.  Freeze

```text
h0=1e-10, h1=5e-11, h2=2.5e-11,
H100a=(4*D(h1)-D(h0))/3,
H100b=(4*D(h2)-D(h1))/3,
H140b=(4*D(h2)-D(h1))/3 at 140 digits.
```

Add (3) independently to each and use the real symmetric part for reduction.
For every schedule define

```text
e_H=100*max(||H100a-H100b||max,
            ||H100b-H140b||max)
    +1e-50*max(1,||H140b||max).                   (4)
```

Require raw imaginary part and antisymmetry no larger than `e_H`.

## 4. Constrained reduction without a pseudoinverse

For each symmetric ladder matrix split

```text
H=[[A,B],[B^T,C]].
```

For each boundary basis `U` in `{P,R P,P_alt}` and internal complement
`V` in `{Q,Q_alt}` when requested, solve

```text
C_V=V^T C V,
C_V Y=-V^T B^T U,
X=V Y,
K(U,V)=U^T(A U+B X).                              (5)
```

Use linear solves, not an inverse or pseudoinverse.  For the primary
`K(P,Q)`, compute `C_Q` at all three ladder levels and define

```text
e_C=100*max(||C_Q,100a-C_Q,100b||max,
            ||C_Q,100b-C_Q,140b||max)
    +1e-50*max(1,||C_Q,140b||max).                (6)
```

Require every eigenvalue of the real symmetric `C_Q,140b` to exceed `100e_C`.
Require reduced solve residual below `1e-60`.  The full internal residual
`B^T U+C X`, which is zero only after compatibility is used, must be no larger
than

```text
e_r=100*max(||r100a-r100b||max,
            ||r100b-r140b||max)
    +1e-45*max(1,||r140b||max).                   (7)
```

For every requested `(U,V)` define

```text
e_K=100*max(||K100a-K100b||max,
            ||K100b-K140b||max)
    +1e-45*max(1,||K140b||max).                   (8)
```

Require antisymmetry no larger than `e_K`.  Store every primary `11 x 11`
matrix, internal reduced spectrum, lifts, residuals and envelopes.

The evidential object is the bilinear form on `ker(c^T)`, not any `12 x 12`
extension.  It is the restricted linearized boundary-momentum response modulo
the conormal `c`; no nonlinear admissible boundary surface is assumed.

## 5. Basis invariance, reversal and complete class census

For all 24 schedules require

```text
K(P_alt,Q) = T^T K(P,Q) T,
K(P,Q_alt) = K(P,Q),                               (9)
```

with maximum error no larger than ten times the larger propagated `e_K`
(including the factor `max(1,||T||max^2)` for the first comparison).  Also
require both internal representatives to solve the full compatible equation
within (7).

For schedule order `o` and its reversed order `rev(o)`, require

```text
K_o(P,Q)=K_rev(o)(R P,Q)                          (10)
```

inside the larger of their envelopes.

Canonicalize only temporal orientation:

```text
K_can,o = K_o(P,Q)       if o <= rev(o),
          K_o(R P,Q)     otherwise.               (11)
```

In lexicographic schedule order, assign a matrix to the first prior class
representative within the larger envelope (8), or start a new class.  Store
and print the complete class membership before assigning an outcome.  No
schedule selection, averaging or additional permutation is allowed.

## 6. Controls

1. For

   ```text
   A=[[5,1],[1,4]], C=diag(2,0),
   B=[[1,3],[0,6]], n=(0,1), c=(3,6), p=(2,-1),
   ```

   the constrained scalar response must equal exactly `18`; the incompatible
   direction `(1,0)` must fail the internal solvability condition.
2. The selected-matter minus gravity-only Hessian must have precisely the four
   entries (3), with error `<1e-80`.
3. For schedules `0,1,22,23`, lift three frozen coefficient directions in the
   `P` basis: the first unit vector, all ones, and alternating signs.  Centred
   complete-action second differences at `1e-10` and `5e-11`, Richardson
   combined, must reproduce `y^T K(P,Q)y` with relative error `<1e-28`.
4. Adding `1e-6*max(1,||K||max)` to one component of a copied canonical matrix
   must be detected as inequivalent.
5. Setting the largest component of `c` to zero must make the frozen basis or
   compatibility controls fail by a resolved amount.
6. No Moore--Penrose inverse, root search, nested census, nonhomogeneous
   spectrum, continuum comparison or physical-constant extraction is run.

## 7. Frozen outcomes

Use the first applicable outcome:

1. `REFINED_H4_CONSTRAINED_RESPONSE_CONTROL_FAILED` for provenance, topology,
   on-shell, precision, null/coupling, basis, solve, residual, symmetry,
   directional, corruption or scope failure;
2. `REFINED_H4_CONSTRAINED_RESPONSE_INTERNAL_COMPLEMENT_SINGULAR` if any
   reduced `9 x 9` internal block fails (6);
3. `REFINED_H4_CONSTRAINED_RESPONSE_TIME_REVERSAL_FAILED` if (10) fails;
4. `REFINED_H4_CONSTRAINED_RESPONSE_MULTIPLE_SCHEDULE_CLASSES` if more than
   one class remains under (11);
5. `REFINED_H4_CONSTRAINED_RESPONSE_SINGLE_SCHEDULE_CLASS` otherwise.

Outcome 4 is a **DERIVED COMPUTATIONAL NEGATIVE** for a canonical bare
staircase response in the invariant sector, subject to adversarial
replication.  Outcome 5 advances only to constructing a nonhomogeneous
quadratic response; it does not establish an exact nonlinear constraint
surface, propagation, a tick, `c` or `G`.

## 8. Execution and acceptance boundary

Register the verifier before its first execution.  Run it twice and require a
byte-identical artifact.  Perform only the static registry audit.  Do not run
the full suite or the deferred nonlinear census.  No outcome with physical or
mathematical weight is accepted before a mechanically different adversarial
reconstruction.
