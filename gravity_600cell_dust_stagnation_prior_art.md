# Prior-art gate: nonlinear defect stagnation on the 600-cell dust slab

Date: 2026-08-13

Status: **completed before computing a displaced local Jacobian or trying a
local-Jacobian step at any stored nonlinear state**.

This search is not a novelty proof.  External novelty of the carrier-specific
matrices remains **OPEN**.

## Exact object and hypotheses

On the already certified De Felice--Fabri time-symmetric 600-cell dust slab,
retain the frozen order-24 invariant carrier, both schedule parities, the
four committed zero-sum Helmert phase contrasts, both signs, amplitude
`eta=1e-4` and the five collective values in `|t|<=0.1`.

For each of the 80 final states stored by the first nonlinear run, the second
frozen verifier used

```text
F(z;t,b) = Q^T E(log u_path(t)+Qz,b),
p_H      = H_Q^-1 F,
```

where `H_Q` is the quotient Hessian at the regular background.  The method
stopped at 72 states because no damping in `1,1/2,...,1/1024` reduced
`norm(F)`.  The proposed diagnostic compares this frozen linear model with
the actual local Jacobian

```text
J(z;t,b) = partial F / partial z
```

at the same stored states.  It does not change the boundary data, amplitude,
carrier or equations.

## KNOWN: fixed-Jacobian failure is not root nonexistence

The residual condition underlying an inexact Newton step is standard.  The
linear step must make

```text
norm(F + J s) <= eta * norm(F),  eta < 1,
```

and the forcing term controls convergence and robustness.  See the primary
paper:

[Choosing the forcing terms in an inexact Newton method](https://www.osti.gov/biblio/218521)

The present fixed-Hessian defect step uses the background `H_Q` in place of
the displaced `J`.  Therefore its stagnation can result from an inaccurate
local linear model even when a nearby root exists.  This numerical mechanism
is **KNOWN**; observing it here would not be a new method.

## KNOWN: nonlinear Regge pseudo-constraints can change the soft sector

On curved Regge solutions exact discrete gauge symmetry is generically
broken and constraints are replaced by pseudo-constraints:

[(Broken) Gauge Symmetries and Constraints in Regge Calculus](https://arxiv.org/abs/0905.1670)

Higher-order Regge equations can impose consistency conditions on background
gauge parameters, and the quadratic constraints become background-dependent:

[From covariant to canonical formulations of discrete gravity](https://arxiv.org/abs/0912.1817)

Canonical simplicial evolution likewise permits initially free data to be
fixed by later constraints:

[Canonical simplicial gravity](https://arxiv.org/abs/1108.1974)

Thus displacement or rotation of the four soft relative-lapse modes is a
**KNOWN possible mechanism**, not evidence for new dynamics.

## CONTROL

- the complete analytic logarithmic equation and the independently certified
  complete-action implementation;
- the exact base complement `Q` and committed precision `H_Q` matrices;
- all 16 boundary cases and all 80 stored terminal states;
- the 100-decimal three-step action derivative estimator already used for
  candidate rejection;
- the certified Lorentzian branch gates.

## OPEN

- whether the displaced local Jacobian is numerically resolved at the soft
  scale;
- whether it differs materially from `H_Q`, in particular in the four soft
  directions;
- whether a local-Jacobian step supplies descent where the fixed-Hessian step
  did not;
- whether the binary analytic equation is above its independently estimated
  complete-action error at representative displaced states;
- existence or nonexistence of a stationary root at any fixed nonzero
  boundary deformation.

## Framing correction

The proposed diagnostic can distinguish an unresolved numerical floor from
a failed fixed linear model.  It cannot establish the absence of a root.
Even a well-resolved nonsingular local Jacobian and a failed line search are
only local numerical facts.  A genuine no-root result would require a
separately specified domain plus an interval-Newton/Krawczyk exclusion,
topological degree argument, or another exhaustive certificate.

Accordingly, the third verbal possibility--"real absence of a root"--must
remain **OPEN** after this diagnostic unless such a certificate is built.  No
solver failure will be cited as evidence for a claim its search space could
not falsify.

## Decision

Preregister a non-solving diagnostic which:

1. computes converged centered-difference local Jacobians at all 80 stored
   states;
2. measures the inexact-Newton forcing factor of the frozen `H_Q` step;
3. tests exactly one local-Jacobian descent step without continuing to a
   root;
4. uses deterministic complete-action anchors to quantify the binary
   equation floor;
5. reports `PRECISION_LIMITED`, `FIXED_MODEL_MISMATCH` or
   `LOCAL_MODEL_INCONCLUSIVE`, never `NO_ROOT`.
