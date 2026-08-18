# Boundary-coboundary test: OPEN because the Dirichlet product was invalid

Date: 2026-08-17

Only the targeted verifier was run.  The full suite was not run.

## 1. Provenance

- prior-art gate: `5ccb29b`
- preregistered protocol: `84c8fd8`
- verifier registered before its first launch: `fea7462`
- control-only implementation correction: `3584ab7`
- result artifact:
  `reproducible/gravity_600cell_dust_boundary_coboundary.json`
- result SHA-256:
  `accebd761c59d7bfda5aa41857505ce2e39e6e4f4def67625017289874e7785f`

The first launch evaluated zero actions because of a missing standard-library
import.  That failure and the direction-serialization tolerance correction are
recorded in the protocol.  No case, amplitude or decision rule changed before
the scientific launch.

## 2. Complete hypotheses of the attempted claim

For each of the two order-24 600-cell dust schedules, the attempted
construction fixed arbitrary nearby old and new logarithmic boundary data,
solved all 35 internal equations, and intended to define an on-shell action

```text
W_p(o,n)=S_p(o,x_p(o,n),n).
```

It then intended to test whether

```text
W_odd(Po,Pn)-W_even(o,n)=B_old(o)+B_new(n)+constant
```

by 128 finite mixed rectangles.  The directions, signs, amplitudes, physical
permutation, action and branch were frozen.  Both operational and validation
fixed-J solvers were required at 161 unique boundary points, for 644 internal
solves in total.

The missing hypothesis was decisive: the 35 by 35 internal Dirichlet Hessian
must be invertible, or an independently established constrained local
principal function must replace the product-neighbourhood function above.
Neither holds in the attempted formulation.

## 3. Mechanical output

The targeted verifier passed its `10/10` bookkeeping and control assertions
and assigned the preregistered non-scientific outcome

```text
BOUNDARY_COBOUNDARY_OPEN.
```

The solve census was

```text
CONVERGED       4 / 644  (the base point, two parities, two calibrations)
NO_ARMIJO_STEP 640 / 644
```

All 128 rectangles were consequently `OPEN_SOLVE`; no mixed rectangle was
evaluated and no separability classification exists.  Among failed solves:

- 506 stopped after the first attempted Newton step;
- the rest stopped after two to six accepted steps;
- residual infinity norms at failure ranged from `3.5095e-3` to `8.3537e-2`;
- the fixed-base-J corrections ranged from `6.7845` to `1609.12` in internal
  log coordinates.

The four base solves converged to residuals below `7e-71`.  Thus the action and
base implementation work; the failure is tied to the arbitrary displaced
two-boundary data.

## 4. Why the framing fails

This repository had already derived, at the same background:

```text
rank(K_xx)=34,
one collective internal lapse null,
one nonzero boundary compatibility equation,
compatible final-boundary tangent space = 29-dimensional zero-sum space.
```

Those statements are in
`gravity_600cell_dust_gauge_quotient_precision_result.md` and
`gravity_600cell_dust_exact_lapse_path_result.md`.  They imply that the usual
implicit-function theorem does not define `x_p(o,n)` on an arbitrary
60-dimensional product neighbourhood.  The canonical 65 by 65 pre-Legendre
problem becomes regular only after initial momentum is included; that fact
does not repair the fixed-two-boundary Dirichlet problem.

The preregistered rectangles varied old and new shape directions independently
and therefore did not enforce the missing compatibility equation.  The huge
fixed-J corrections are consistent with attempting to invert the collective
null direction.  They are not evidence that an endpoint term exists or does
not exist.

This was an avoidable framing error: the prior-art gate mentioned the regular
canonical map but failed to propagate the repository's explicit rank-34
Dirichlet result into the hypotheses.  The error is recorded rather than
hidden by changing the solver after seeing the failures.

## 5. Post-result primary-source reconciliation

- Dittrich and Hoehn, [*Constraint analysis for variational discrete
  systems*](https://arxiv.org/abs/1303.4294), formulate discrete variational
  evolution as constrained pre/post canonical relations.  A regular
  unconstrained two-point generating function cannot be assumed when the
  relevant Legendre or Hessian block is singular.
- Dittrich and Hoehn, [*From covariant to canonical formulations of discrete
  gravity*](https://arxiv.org/abs/0912.1817), show how nonlinear Regge dynamics
  produces consistency conditions and pseudo-constraints.
- Dittrich and Hoehn, [*Canonical simplicial
  gravity*](https://arxiv.org/abs/1108.1974), use Hamilton's principal function
  within the full pre/post constraint formalism, not as permission to prescribe
  arbitrary incompatible boundary configurations.

These are known structural issues.  No external novelty is claimed.

## 6. Status ledger

| Claim | Status |
|---|---|
| The targeted artifact is reproducible and registered | **DERIVED / 10 of 10 controls** |
| The base internal stationary solution is recovered | **DERIVED COMPUTATIONAL** |
| Arbitrary nearby `(o,n)` admit internal stationary roots | **REFUTED AS AN ASSUMPTION / not established** |
| A smooth scalar `W_p(o,n)` exists on the tested product neighbourhood | **OPEN; the attempted hypothesis is invalid** |
| The two schedule actions differ only by endpoint terms | **OPEN** |
| Endpoint separability is refuted by this run | **NO; no rectangle was available** |
| The canonical initial-value map remains locally regular | **DERIVED UPSTREAM; unaffected** |
| The nonlinear schedule defect is physical | **OPEN** |

## 7. Correct continuation

The endpoint-term question must be reformulated for constrained generating
families or canonical relations.  Two legitimate routes remain:

1. compare the configuration projections of the two constrained canonical
   relations, enforcing the derived compatibility condition and quotienting
   the collective lapse;
2. compare the regular initial-value maps modulo vertical momentum shifts
   generated by endpoint functions, with every integrability condition stated
   before inspecting component defects.

Simply switching to an adaptive 35-dimensional Newton method would not repair
the missing hypothesis and could select a distant pseudo-constraint branch.
The present 128-case artifact must therefore remain **OPEN**, not be rerun with
a more permissive solver.

## 8. Consequence for the theory

There is no new physical result from this attempted calculation.  What is
gained is a sharper boundary: the dust slab is naturally a constrained
canonical relation, not an ordinary unconstrained two-point propagator on all
boundary geometries.  Any future claim about a physical tick, schedule
equivalence or perfect action must respect that constraint structure.
