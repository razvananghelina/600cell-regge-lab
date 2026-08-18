# Prior-art gate: repairing the complete-action transverse solver

Date: 2026-08-14

Frozen numerical-boundary result: `64a13f6`

Status: written before evaluating any new complete-action row after result
`64a13f6`.

## 1. Exact object and hypotheses

At each of the same 80 frozen `(parity,direction,sign,t,boundary)` states,
the mathematical problem is to solve

```text
F(z)=Q^T E_action(t,z,boundary)=0,   z in R^34,
```

where `E_action` is reconstructed from 100-decimal central differences of
the complete Lorentzian Regge-plus-dust action.  The collective scalar

```text
g=w(t)^T E_action
```

may be classified only after an independently validated transverse root.

The available local matrix `J_b` is the Richardson-extrapolated Jacobian of
the binary analytic equation, not a derived Jacobian of the complete-action
row.  It may therefore be used as an invertible proposal generator or norm
preconditioner, but its entries are not licensed as the physical residual.

Result `64a13f6` supplies the new controls:

- the complete-action solver and validation rows agree to
  `3.47e-16--5.72e-16`;
- their old error proxies differ by `80.978--80.982`, the expected
  fourth-order step ratio near `3^4=81`;
- all 63 solver-zero labels fail independent validation;
- 11 of 17 no-descent states lie in an algebraic dead zone of the old
  residual-norm acceptance inequality;
- raw preconditioned trials nevertheless reduce the residual strongly in
  most of those states.

These facts diagnose the numerical method.  They do not establish a root or
its absence.

## 2. Primary prior art

### KNOWN: inexact Newton forcing terms

Dembo--Eisenstat--Steihaug introduced inexact Newton iterations in which the
linearized residual is reduced only to a controlled fraction of the current
nonlinear residual.  Eisenstat--Walker later made the forcing term adaptive
and scale-independent, explicitly warning that solving the local Newton
equation more accurately than the nonlinear model warrants can give little
or no decrease in the true residual.

Primary sources:

- Dembo, Eisenstat and Steihaug, *Inexact Newton Methods*, SIAM Journal on
  Numerical Analysis 19 (1982), 400--408,
  <https://doi.org/10.1137/0719025>.
- Eisenstat and Walker, *Choosing the Forcing Terms in an Inexact Newton
  Method*, SIAM Journal on Scientific Computing 17 (1996), 16--32,
  <https://doi.org/10.1137/0917003>.

The latter paper also gives a globalized inexact Newton algorithm with
backtracking on the actual nonlinear residual.  Therefore neither damping,
residual backtracking nor an accuracy forcing term is novel here.

### KNOWN: line search and trust-region globalization

Line-search and trust-region modifications are standard devices for making a
Newton proposal conditional on decrease of the actual nonlinear problem.
General inexact-Newton dogleg methods and their convergence analysis already
exist.

Primary source:

- Pawlowski, Simonis, Walker and Shadid, *Inexact Newton Dogleg Methods*,
  SIAM Journal on Numerical Analysis 46 (2008), 2112--2132,
  <https://doi.org/10.1137/050632166>.

### KNOWN: affine-covariant Newton corrections

The Newton correction `J(x)^-1 F(x)` and affine-covariant/natural level
functions are standard ways to account for badly scaled nonlinear systems.
This is the relevant mathematical precedent for measuring progress in a
preconditioned residual norm rather than treating every component of `F` as
equally scaled.

Reference:

- Deuflhard, *Newton Methods for Nonlinear Problems: Affine Invariance and
  Adaptive Algorithms*, Springer Series in Computational Mathematics 35,
  <https://doi.org/10.1007/978-3-642-23899-4>.

Using a fixed nonsingular `J_b` at one iterate makes
`norm(J_b^-1 F)` an ordinary norm of the complete-action residual.  Descent
in that norm is mathematically meaningful even if `J_b` is only a
preconditioner.  Updating the norm at the next accepted iterate does not turn
`J_b` into the physical equation.

## 3. What the literature does not supply

The cited methods do not answer this repository-specific problem:

1. whether the De Felice--Fabri 600-cell dust sandwich has a transverse root
   on the frozen order-24 carrier;
2. whether the reduced scalar vanishes at such a root;
3. whether the binary analytic Jacobian is a sufficiently good proposal
   generator for the complete Regge action near `1e-12` residuals;
4. how the deterministic Richardson truncation proxy should be separated
   from the root-acceptance and final-validation windows on this action;
5. whether either schedule parity defines a physical tick.

All five remain **OPEN**.

## 4. Framing attack

The previous proposed repair, "use a smaller window", is insufficiently
specified.  A smaller solver window followed by a still smaller validation
window can reproduce the same false-stop hierarchy at a new scale.  Merely
relaxing the ten-error threshold would be worse: it would convert an observed
failure into a fitted tolerance.

Likewise, `norm(J_b^-1 F)` cannot be called a physical action or a new
equation.  It is only a scale-sensitive merit norm for the already fixed
complete-action equation `F=0`.

The safe repair must therefore separate three roles:

1. the complete-action row is the only nonlinear residual;
2. `J_b` proposes a step and defines a temporary invertible norm, while
   actual action rows decide acceptance;
3. final zero and scalar labels use a disjoint, smaller action window and are
   never fed back into the final reported validation.

## 5. Proposed difference and evidence label

The repository-specific continuation will combine known ingredients:

- complete-action residuals at high precision;
- an inexact binary-Jacobian proposal;
- backtracking on an active complete-action merit;
- a separate validation window;
- fixed look-elsewhere accounting over the same 80 states.

The only project-specific contribution is the exact placement of these known
tools around the already derived 600-cell action and the observed
Richardson-scale failure.  That is **STRUCTURAL / numerical method design**,
not a new Newton theorem.

Because the positive raw scalar was inspected before this repair, any future
common-sign outcome retains **PATTERN-informed provenance** even after a new
preregistration.  A validated scalar zero would be a
**DERIVED COMPUTATIONAL LOCAL** hit on the frozen carrier, subject to its
predeclared 16/8/4 look-elsewhere count.

External novelty is **OPEN**; this search was a focused prior-art gate, not a
systematic review.
