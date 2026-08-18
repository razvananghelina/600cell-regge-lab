# Preregistration: stationary roots of the 35-variable global Regge slab

Date: 2026-08-13

Certified evaluator commit: `d9fe159`

Status: **preregistered before any stationary-root iteration has been run**.

No value from a particle target, `a1=5`, a coupling, a speed or a Planck unit
is used.  This protocol tests existence of vacuum stationary points in the two
already derived phase-parity fixed subspaces.

## 1. Complete scope

For each of the two `H4`-inequivalent ordered phase parities:

1. use its certified staircase slab of 2400 Lorentzian four-simplices;
2. fix all 1440 initial/final boundary squared edge lengths to one;
3. use all 30 diagonal-length orbits and all five positive pole-magnitude
   orbits under the exact order-24 schedule stabilizer;
4. assign diagonal squared lengths `q_j>0` and pole squared lengths
   `-rho_k<0`;
5. use the corrected plus complex-angle branch and zero-volume Lorentzian
   Regge action already certified in commit `d9fe159`;
6. require every four-simplex to retain inertia `(-,+,+,+)` and every angle
   argument to stay off its logarithmic cut;
7. solve all 35 restricted internal-edge equations.  Since every orbit has
   size 24, a restricted zero implies all 840 individual derivatives vanish.

This is a fixed-symmetry search.  A root is a root of the complete slab action;
failure is not a no-go outside the 35-dimensional invariant subspace.

## 2. Frozen starts

Exactly six starts are used for each parity:

```text
S0: q_j=1,                     rho_k=1/4
S1: q_j=1+(j+1)/1000,          rho_k=1/4+(k+1)/1000
S2: q_j=1-(j+1)/2000,          rho_k=1/4+(5-k)/1500
S3: q_j=1+0.02 cos(2 pi(j+1)/31),
    rho_k=1/4+0.01 cos(2 pi(k+1)/5)
S4: q_j=1+0.02 sin(2 pi(j+1)/31),
    rho_k=1/4+0.01 sin(2 pi(k+1)/5)
S5: q_j=1+0.01(-1)^j,          rho_k=1/4+0.005(-1)^k.
```

Here `j=0,...,29`, `k=0,...,4`.  The phase/orbit ordering is the
lexicographic ordering frozen by the evaluator.  If a start is noncausal it is
recorded and rejected, not modified.

These starts test the regular neighborhood.  They do not constitute a global
exhaustion of the positive cone.

## 3. Variables and residual

Optimize logarithmic variables

```text
y_i = log x_i,
x = (q_0,...,q_29,rho_0,...,rho_4).
```

The residual is the per-edge action gradient

```text
r_i(y) = (1/24) dS/dx_i.
```

The factor 24 removes the common orbit multiplicity and changes no zero.  No
component-dependent rescaling or fitted weighting is permitted.  The
Jacobian with respect to `y` is

```text
J_ij = (1/24) H_ij x_j,
```

where `H` is obtained by centered differences of the certified analytic
Schlaefli gradient.

## 4. Frozen causal Levenberg--Marquardt iteration

At every iteration:

1. compute `r` and a centered Jacobian with relative `y` step `2e-5`;
2. use the symmetrized `x`-Hessian only as a diagnostic; solve the residual
   least-squares step

   ```text
   (J^T J + mu I) delta = -J^T r;
   ```

3. try `mu/sigma_max(J)^2` in the fixed list

   ```text
   0, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 1, 100;
   ```

4. for each candidate use backtracking factors `2^-b`, `b=0,...,20`;
5. accept the first candidate satisfying all causal gates and

   ```text
   ||r_new||_2 <= (1-1e-4*alpha) ||r_old||_2;
   ```

6. stop successfully when both

   ```text
   ||r||_inf < 1e-10,
   ||r||_2   < 3e-10;
   ```

7. stop unsuccessfully after 80 accepted iterations, no acceptable step, a
   nonfinite value, or relative Jacobian condition number above `1e14` for
   three consecutive iterations.

The causal gates are:

- all 100 simplex-orbit representatives have exactly one negative Gram
  eigenvalue;
- minimum absolute Gram eigenvalue exceeds `1e-8`;
- minimum angle-argument modulus exceeds `1e-6`;
- action and gradient imaginary residuals are below `2e-7` relative;
- every `q_j` and `rho_k` lies strictly between `exp(-6)` and `exp(6)`.

A candidate that reaches the artificial box boundary is rejected and called
inconclusive, not a root or no-go.

## 5. Frozen post-root validation

Every converged candidate must pass all of the following.

1. Re-evaluate the unreduced 2400-simplex action and all 840 individual
   internal-edge gradients.
2. Require

   ```text
   max individual |dS/ds_e| < 2e-10.
   ```

3. Recompute the 35 Jacobian/Hessian at relative steps `1e-4`, `3e-5`, and
   `1e-5`; require stable rank at `1e-7,1e-9,1e-11` and record all spectra.
4. Compare centered differences of the complete action in all 35 orbit
   directions with the analytic gradient at steps `1e-4` and `3e-5`, requiring
   relative agreement below `3e-5`.
5. Round the candidate to 14 decimal digits, restart the frozen solver, and
   require convergence to the same root within relative `2e-10`.
6. Cluster candidates using relative Euclidean distance `1e-8`; report the
   exact number of distinct roots found and which starts reach each root.
7. Verify roots of the two phase parities independently.  No parity is called
   selected merely because one numerical start failed in the other.

A float64 candidate passing these controls is labeled **DERIVED NUMERICAL
WITNESS**, not an exact theorem.  Before any physical claim, it must receive a
separate arbitrary-precision or interval certificate whose protocol is frozen
after publishing the candidate digits.

## 6. Acceptance and negative boundaries

### Acceptance

The invariant vacuum route advances for a parity if at least one candidate
passes every post-root validation.  If exactly one parity has validated roots
from at least two distinct starts while all six searches in the other parity
terminate regularly without boundary contact, this is a **PATTERN candidate
for dynamical parity selection**, not yet a theorem of uniqueness.

### Honest negative

If no root is found, report only:

> no stationary point was found from the six frozen starts in the stated
> causal box and invariant subspace.

That is not evidence of global nonexistence and does not close the full route.

### Kill boundaries

- failure of the full 840-gradient reconstruction kills a candidate;
- loss of Lorentzian inertia or an angle-branch crossing kills a candidate;
- convergence only on the artificial box boundary is inconclusive;
- a parity conclusion based on unequal search effort or different starts is
  forbidden;
- an exact/physical claim without subsequent high-precision certification is
  forbidden.

Only the targeted search verifier will be run.  The full suite remains out of
scope.
