# Preregistration: finite-height carrier quadratic canonicity

Date: 2026-08-22.

Prior-art gate commit: `5fab3f5`.

Status: **FROZEN BEFORE THE FIRST FINITE-HEIGHT FULL-HESSIAN ASSEMBLY.**

No quadratic parity difference, spectrum, eigenvector or continuum target has
been inspected for the background specified here.

## 1. Frozen inputs

The verifier must reject any input hash mismatch.

```text
reproducible/gravity_600cell_finite_height_fourth_slab.json
  cf322cf0d60668d8f3f58e251425c9ad6bf43b112f22f9f3aebbc28f86212468

reproducible/gravity_600cell_homothetic_frustum_action.json
  c0226a47607113930a31259d0cbee8ea33df2f7b0ba9416f9dbe5d647cede52d

reproducible/gravity_600cell_full_scale_strut_canonical_precision.json
  75351ae4dfde26dd75ed8faa927b0a49cd725d83c7629d4545268030b54e2706

reproducible/gravity_600cell_full_scale_strut_symbolic_gap_resolution.json
  ea2c52f0cd227516734defc509330e528b140f71bfd0f50e87036f3fa9832179

reproducible/verify_gravity_global_regge_orbits.py
  ad93cdd08fabeeee56b009f23936696837c4362f88ae23f92a36d0395e61ffaf

reproducible/verify_gravity_600cell_dust_full_anisotropic_legendre_rank.py
  834b97c85e386def853b6308e65e831c52d62d7cbcc4b23118602120d6c676e5

reproducible/verify_gravity_600cell_full_scale_strut_canonical_intersection.py
  a2d5390d39c725a5fb586fefce9da34cede3a1fb84bbe36791f8b0599b3eae42
```

The last two source files are implementation libraries only.  Their old
background constants, old target matrices and bottom-level verifier bodies
must not be executed or imported as scientific inputs.

## 2. Independently reconstructed background

At 120 decimal digits define

```text
z(q)       = (q^2+2)/(2(q^2+3)),
epsilon(q) = 2*pi-5*acos(z(q)),
mu(q)      = 180*epsilon(q)/(pi*sqrt(q^2+4)),
p(q)       = 180*q*epsilon(q)/sqrt(q^2+4)
             -600*sqrt(3)*asinh(q/sqrt(8(q^2+3))).
```

For `v=3/2`, solve the exact scalar equation

```text
E(q)=4*pi*(mu(q)-mu(v))+q*(p(q)-p(v))=0
```

by deterministic bisection on the frozen bracket `(9,10)`.  Then reconstruct

```text
h      = (p(q)-p(v))/(2*pi*mu(q)),
L0     = 1,
lambda = 1+h*q,
rho    = h^2,
M      = mu(v).
```

Require `E(q)` below `1e-90`, agreement with the committed representative
through 70 decimal digits and the positive physical gates

```text
h>0, lambda>0, rho>0.
```

For the local carrier conversion set

```text
tau^2 = 3*(lambda-1)^2+8*rho/L0^2
```

and require exactly the generic-domain exclusions

```text
lambda != 1,
tau != 0,
(lambda-1)^2-3*tau^2 != 0,
q_diag=lambda*L0^2-rho != 0.
```

No root outside `(9,10)` may replace this one.

## 3. Complete action derivatives

For each parity representative `even` and `odd`, independently reconstruct
the full 2,400-simplex slab and its 2,280 logarithmic signed-squared-edge
variables.  Use the exact same Lorentzian branch, boundary term, zero
cosmological constant and conserved dust term as the certified homogeneous
cellular action.

Set the base squared lengths to

```text
old boundary       L0^2,
artificial diagonal lambda*L0^2-rho,
same-vertex strut  -rho,
new boundary       lambda^2*L0^2.
```

The complete gradient is evaluated directly at arbitrary precision.  The
Hessian is assembled without numerical differentiation of the full action:
triangle-area derivatives are analytic and dihedral-angle derivatives use
the four frozen centered logarithmic steps

```text
operational_primary  1e-20,
operational_shadow   1e-15,
validation_primary   3e-20,
validation_shadow    3e-15.
```

All 2,400 simplices and displaced local patterns must retain Lorentzian
inertia `(3,1)`, positive leading-minor and angle-argument margins, and
imaginary contamination below `1e-60`.

Required first-derivative gates are:

```text
maximum artificial/internal gradient              < 1e-25,
maximum even/odd old-boundary edgewise difference  < 1e-25,
maximum even/odd new-boundary edgewise difference  < 1e-25,
maximum binary64/direct-gradient assembly error    < 2e-11.
```

The edgewise comparisons use physical edge labels, not orbit positions.
Failure of any gate is `CONTROL_FAILED` and forbids interpreting a Hessian
difference.

For every Hessian approximation report raw antisymmetry.  It must lie inside
ten times the sum of the cross-step derivative discrepancy and the existing
forward-summation roundoff envelope.  Do not symmetrize a matrix to make this
gate pass.

## 4. Exact geometry-selected tangent carrier

Build a fresh `1560 x 240` matrix for each parity.  Rows are the 840 internal
edges followed by 720 new-boundary edges, identified from the complete model
by physical edge labels.  Columns are fixed globally as

```text
sigma_0,...,sigma_119,c_0,...,c_119.
```

With

```text
A = -16*rho/(L0^2*(lambda-1)^2),
B = 8+16*rho/(L0^2*(lambda-1)^2),
a = L0^2*A/(8*q_diag),
b = L0^2*B/(8*q_diag),
k = rho/((lambda-1)*q_diag),
```

the nonzero log-length responses are frozen as follows.

```text
pole at u:
  G[pole,c_u]=1

oriented diagonal u(lower)->w(upper):
  G[diagonal,sigma_u]=a
  G[diagonal,sigma_w]=b
  G[diagonal,c_u]=k
  G[diagonal,c_w]=-lambda*k

new edge {u,w}:
  G[new edge,sigma_u]=1/lambda
  G[new edge,sigma_w]=1/lambda.
```

Require the exact support census

```text
120 rows of support 1,
720 rows of support 2,
720 rows of support 4,
scale-column support 24 at every vertex,
strut-column support 13 at every vertex.
```

Require the pole identity and connected non-bipartite upper graph proof of
exact rank 240.  Floating SVD is diagnostic only.  Also require the two
collective identities

```text
sum sigma: pole=0, diagonal=L0^2/q_diag, new=2/lambda,
sum c:     pole=1, diagonal=-rho/q_diag, new=0.
```

## 5. Frozen quadratic comparison

Restrict each complete Hessian approximation to its internal-plus-new rows
and columns.  For every derivative scheme `s`, compute

```text
Q_p,s     = G_p^T H_p,s G_p,
B_p,s     = (Q_p,s+Q_p,s^T)/2,
Delta_s   = B_even,s-B_odd,s.
```

No spectrum or eigenvector of `B_p,s` may be computed before the parity
classification is frozen.  Matrix two-norms are allowed only for numerical
error and outcome classification.

Use

```text
N = max(1, ||B_even,op||_2, ||B_odd,op||_2),
d_s = ||Delta_s||_2/N.
```

The target-independent uncertainty is

```text
e_step = max over p,s ||B_p,s-B_p,op||_2/N,

e_round = [1560*entry_bound_even*||G_even||_2^2
           +1560*entry_bound_odd *||G_odd ||_2^2]/N,

e_total = e_step+e_round+100*machine_epsilon.
```

Here `entry_bound_p` is the previously defined forward-summation maximum
absolute Hessian-entry error envelope, recomputed at this background without
using an observed parity difference.

The representative comparison covers the two H4 schedule orbits only after
the verifier checks the committed census of 60 schedules per parity and the
within-orbit edge-permutation covariance.  Otherwise the conclusion remains
representative-scoped.

## 6. Controls that must pass and fail

1. Reordering each parity Hessian and carrier independently by lexicographic
   physical edge labels must leave its `B_p,op` unchanged inside
   `10*e_total*N`.
2. On the two-dimensional span of the uniform `sum sigma` and `sum c`
   columns, the even/odd difference must lie inside `10*e_total*N`.  This is
   the known homogeneous subdivision control, not the main result.
3. Add the synthetic symmetric perturbation

   ```text
   H_odd -> H_odd + u*u^T,
   u=G_odd[:,0]/||G_odd[:,0]||_2.
   ```

   Its induced rank-one carrier perturbation is known algebraically and must
   be detected above `100*e_total` relatively.
4. Change the lexicographically first oriented-diagonal source-scale
   coefficient `a` by `+1/10`.  The corresponding pulled-back form must
   change by more than `100*e_total` relatively.  This is a carrier
   corruption control, not an alternate admissible geometry.

If a hostile control is accidentally below the main error floor, the outcome
is `CONTROL_FAILED`; no coefficient or threshold may be changed after seeing
the result.

## 7. Frozen outcome hierarchy

### `FINITE_HEIGHT_QUADRATIC_CONTROL_FAILED`

Use if provenance, root, branch, gradient, reciprocity, carrier, covariance
or hostile-control gates fail.

### `FINITE_HEIGHT_QUADRATIC_PARITY_INDEPENDENT_PRIMARY`

Use only if

```text
max_s d_s <= 10*e_total
```

and every control passes.  Label the result **DERIVED COMPUTATIONAL /
STRUCTURAL, PRIMARY ONLY** pending the mandatory mechanically different
adversarial replication.  It is a necessary quadratic canonicity result, not
a propagation law.

### `FINITE_HEIGHT_QUADRATIC_PARITY_DEPENDENT_PRIMARY`

Use only if

```text
min_s d_s > 100*e_total
```

and every control passes.  Also require the four normalized difference
matrices to agree pairwise within `10*e_total` in two-norm.  Label the result
**DERIVED COMPUTATIONAL NEGATIVE, PRIMARY ONLY** pending adversarial
replication.

### `FINITE_HEIGHT_QUADRATIC_PARITY_OPEN`

Use for the gap between the two thresholds or disagreement among derivative
schemes.

No primary outcome is consolidated before a separately preregistered,
mechanically different high-precision replication.  Two executions of this
same verifier establish reproducibility only.

## 8. Interpretation firewall

This is a one-sided, old-boundary-fixed quadratic action test on an exact
infinitesimal carrier.  A positive result does not establish nonlinear
integrability, a complete canonical boundary map, gauge reduction, two
graviton polarizations, stability, locality, a continuum limit, a limiting
speed, `c`, `G`, a Planck scale or particle masses.

A negative result says only that the present bare length-Regge action plus
the complete geometry-selected scale-plus-strut carrier does not remove the
two-schedule ambiguity at quadratic order on this finite-height background.
Adding a connection, area-angle variables or a perfect action would be a new
theory input and must pass a separate prior-art and preregistration gate.

Only the new targeted verifier and static registry checks may be run.  The
full suite is forbidden for this mission unless the user explicitly asks for
it.
