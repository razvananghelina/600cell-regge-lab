# Protocol adversarial: direct-action constrained H4 response

Date: 2026-08-21

Status: frozen before evaluating any new direct-action matrix in this audit.

Prior-art gate: `8ecbd2a`.  Primary result: `d9f3abc`.

## 1. Question and complete hypotheses

On the frozen internally stationary Lorentzian product over

```text
K0=P(sd K_600), f=(2640,17040,28800,14400), tau0=0.0102,
```

use the corrected Regge action with its boundary terms, the curvature-selected
rank masses, all 24 colour-ordered staircase schedules, the accepted analytic
internal null line `n`, and the accepted nonzero compatibility covector
`c=B n`.  In the fixed total-orbit log-squared-edge coordinates, the object is
the constrained linearized boundary-momentum form on

```text
S=ker(c^T),
```

after stationarity in the nine-dimensional internal slice `ker(n^T)`.

The primary calculation differentiated the analytic gradient to construct a
`22 x 22` Hessian and then performed a restricted block elimination.  This
audit asks whether a mechanically different construction, using only centred
second differences of the **complete scalar action**, independently gives:

1. one response class across all 24 schedules;
2. time-reversal covariance; and
3. the same 24 response matrices as the frozen primary result.

No nonlinear boundary constraint surface, nonhomogeneous mode, continuum
limit, speed, coupling constant, Planck scale or particle observable is part
of the claim.

## 2. Independence and frozen inputs

The decisive direct construction must not import, execute, or call either the
primary Hessian builder or its gradient-differentiation routine.  Load only
function definitions from

```text
reproducible/verify_gravity_600cell_refined_h4_stationary_fill.py
  89aab727792e20a81e7577e0425f8fa4b1e84e2a7ae66caa9e79a4aebf3581e7
```

to obtain the scalar action and exact geometry.  Freeze also:

```text
reproducible/gravity_600cell_refined_local_curvature_mass.json
  180010a79177ba16620ebea9847443c57a7a6d2d8a3df71ad6ecb83f454ef091
reproducible/gravity_600cell_refined_h4_null_coupling_adversarial.json
  5c1f596958f9d878c8d9d3ccb6ecc8359f72164e8f36dd9930fb71ddc1351ce9
reproducible/gravity_600cell_refined_h4_constrained_response_corrected.json
  85adea23f6a19153f61f3ed066137a5e40ab77b8901d4cc81cfc4f864e0bc093
docs/gravity/gravity_600cell_refined_h4_constrained_response_primary_result.md
  fd07977bdb2e45bf3170d1ba98919690e57bc8ee476ce48506859bbffc0253ad
docs/gravity/gravity_600cell_refined_h4_constrained_hessian_prior_art.md
  222f31862e911e03a1a7740696618948e370e43164812120120d85e834f0f639
commons/cell600.py
  ea5bce4b6c52e0834539ca4b1df9c6a67a3a5ed4da32f4e0298a493fc5315c7f
```

The primary artifact may be read only after every direct matrix and direct
class has been constructed.  It supplies the comparison target, not a lift,
Hessian block, matrix entry, class representative or error estimate to the
direct route.

## 3. Algebraic slices

Reconstruct the analytic product tangent from the base coordinates,

```text
n_cross,rs=-tau0^2/q_cross,rs,  n_rho,r=1,
```

and read `c` only from the accepted adversarial null-coupling artifact.  For a
nonzero vector `a` and pivot `p`, use

```text
E(a,p)_j=e_j-(a_j/a_p)e_p, j != p.
```

Freeze

```text
P=E(c,3) in R^(12 x 11),  Q=E(n,9) in R^(10 x 9).
```

The `22 x 20` direct slice is

```text
W=diag(P,Q).
```

Recheck `c^T P=0`, `n^T Q=0`, full column rank, the static internal equations,
and the same real Lorentzian branch before accepting any result.  Alternative
pivots `P_alt=E(c,9)` and `Q_alt=E(n,6)` must be related to the primary bases
by exact algebraic change-of-basis identities; this is a coordinate control,
not a second numerical fit.

Let `R` exchange the six old and six new boundary coordinates.  Solve the
algebraic identity `R P=P T_R` from the eleven identity rows of `P` and require
it componentwise before using `T_R`.

## 4. Direct scalar-action reconstruction

At 180 decimal digits use the complete gravitational plus selected-dust
action.  For a 22-component log-coordinate direction `v`, define

```text
D_h(v)=[S(exp(hv))+S(exp(-hv))-2S(1)]/h^2
```

at the five frozen steps

```text
h=(1e-10, 5e-11, 2.5e-11, 1.25e-11, 6.25e-12).
```

Apply the fixed extrapolation ladder

```text
R_j=(4 D_(j+1)-D_j)/3,
U_j=(16 R_(j+1)-R_j)/15,
V_j=(64 U_(j+1)-U_j)/63.
```

For each schedule evaluate all 210 directions

```text
e_a                         (20 directions),
e_a+e_b, a<b                (190 directions).
```

At each of the two eighth-order levels `V_0,V_1`, reconstruct a symmetric
`20 x 20` matrix `G_l` by polarization:

```text
(G_l)_aa=q_l(e_a),
(G_l)_ab=[q_l(e_a+e_b)-q_l(e_a)-q_l(e_b)]/2.
```

This exhaustively determines the restricted second variation; no random
directions, fitted coefficients or primary lift is allowed.  Track every
action's imaginary part, angle-identity residual and minimum branch argument.

Split each `G_l` as

```text
G_l=[[A_l,B_l],[B_l^T,C_l]], dimensions 11+9,
```

solve with a linear solve

```text
C_l Y_l=-B_l^T,
K_l=A_l+B_l Y_l,
```

and take `K=K_1`.  No inverse, pseudoinverse or analytic Hessian is allowed.
Freeze the direct envelopes

```text
e_G=100 ||G_1-G_0||max + 1e-80 max(1,||G_1||max),
e_K=100 ||K_1-K_0||max + 1e-60 max(1,||K_1||max).
```

Require the smallest eigenvalue of real-symmetric `C_1` to exceed `100 e_G`,
the solve residual divided by `max(1,||B_1||max)` to be below `1e-80`, and
the largest coordinate displacement at the coarsest step to be below `2e-5`.

## 5. Class, reversal and primary comparisons

For an order `o`, let `rev(o)` reverse the four entries.  The direct
time-reversal test is

```text
K_o = T_R^T K_rev(o) T_R.
```

Canonicalize only temporal orientation:

```text
K_can,o=K_o                         if o <= rev(o),
        T_R^T K_o T_R               otherwise.
```

For transformed matrices propagate a max-entry error exactly as

```text
e_can = [max_j sum_i |(T_R)_ij|]^2 e_K.
```

Assign schedules in lexicographic order to the first prior representative
whose max-entry difference is at most the sum of their direct envelopes;
otherwise open a new class.  Store the full membership before assigning an
outcome.  No averaging or schedule selection is permitted.

Only after that census, load the frozen primary artifact.  For every schedule
require

```text
||K_direct-K_primary||max
 <= 10 (e_K,direct + e_K,primary).
```

The factor ten is frozen for cross-method propagation; the corruption control
below must demonstrate that it cannot hide a resolved false match.

## 6. Precision and falsification controls

1. Repeat the complete 210-direction reconstruction for schedules `0` and
   `23` at 220 decimal digits, using the same steps.  Their direct matrices
   must agree with the 180-digit matrices inside ten times the sum of the two
   direct envelopes.
2. On a known symmetric polynomial action
   `S(x)=x^T M x/2+11 sum x_i^4+13 sum x_i^6+17 sum x_i^8+19 sum x_i^10`,
   with

   ```text
   M=[[7,2,-1],[2,5,3],[-1,3,11]],
   ```

   the same five-step and polarization implementation must recover `M` within
   its frozen envelope.
3. Omitting the factor `1/2` in the off-diagonal polarization of that control
   must fail by a resolved amount.
4. Add `1e-6 max(1,||K||max)` to one entry of a copied schedule matrix.  It
   must form a new class and must fail the primary-comparison gate.
5. Reverse the time convention without `T_R`; if that convention happens to
   pass, report it but do not use it as evidence.  The preregistered covariant
   convention remains the evidential test.
6. Run the completed verifier twice and require byte-identical JSON.  This
   checks reproducibility only; mechanical independence comes from the direct
   scalar-action reconstruction.

## 7. Frozen outcomes

Use the first applicable label:

1. `ADVERSARIAL_REFINED_H4_CONSTRAINED_RESPONSE_CONTROL_FAILED` if provenance,
   topology, branch, basis, known-answer, precision or corruption controls
   fail.  No physics conclusion follows.
2. `ADVERSARIAL_REFINED_H4_CONSTRAINED_RESPONSE_DISAGREEMENT` if the controls
   pass but the direct census has more than one class, reversal fails, or any
   direct matrix disagrees with the frozen primary matrix.  The primary
   single-class claim remains **OPEN**, and the disagreement is the result.
3. `ADVERSARIAL_REFINED_H4_CONSTRAINED_RESPONSE_CORROBORATED` only if every
   control passes, all 24 direct matrices form one class, reversal is
   covariant, and every direct matrix matches the primary calculation.

Even outcome 3 establishes only the finite `H4` constrained linear response.
It neither supplies a physical tick nor proves nonhomogeneous gravitational
propagation.
