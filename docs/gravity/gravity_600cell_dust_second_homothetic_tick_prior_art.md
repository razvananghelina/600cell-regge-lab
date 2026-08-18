# Prior-art gate: second canonical homothetic 600-cell dust tick

Date: 2026-08-16

Upstream accepted first tick: `46a7361`.

Status: **completed before evaluating the second-tick seed, Jacobian or
equations**.

This is a targeted primary-source map, not proof of external novelty.

## 1. Exact object and hypotheses

Retain the complete Lorentzian Regge+dust action, the fixed total dust mass,
the regular 600-cell spatial carrier, both derived order-24 staircase
schedules and the exact geometric orbit maps used by the accepted first tick.

For each parity, let the committed first tick provide

```text
a1 = log(L1/L0),
r1 = log(rho1/rho0),
p_post,1 on the thirty final boundary orbits.
```

The second slab has fixed lower geometry `L1`, fixed total mass and fixed
incoming canonical datum `p_post,1`.  Its only two unknowns are

```text
u = log(L2/L1),
v = log(rho2/rho1).
```

The exact homothetic squared lengths are

```text
q_old       = exp(2*a1)*L0^2,
q_new       = exp(2*(a1+u))*L0^2,
pole        = exp(r1+v)*rho0,
diagonal    = exp(2*a1+u)*L0^2 - exp(r1+v)*rho0.
```

The diagonal formula is the exact identity
`L_lower*L_upper-rho2`, not an interpolation.

The old-orbit target is obtained mechanically from the committed
old-to-final orbit map:

```text
target[old] = p_post,1[old_to_final[old]].
```

The equations are

```text
F0(u,v) = mean of all five complete pole equations,
F1(u,v) = mean[p_pre,2-target].
```

Every candidate must subsequently pass all 35 internal equations and all 30
canonical matching components.  The mass, action coefficient, incoming
momentum and lower geometry are not adjustable.

## 2. KNOWN

Iteration of pre/post canonical data is the standard structure of discrete
canonical evolution.  Hamilton's principal function generates the canonical
move, and data introduced at one move may be fixed by later constraints:

- Dittrich and Hoehn, *Canonical simplicial gravity*,
  <https://arxiv.org/abs/1108.1974>.

Consistent discretizations can turn continuum multiplier-like variables into
variables fixed by discrete equations and implement evolution as an
unconstrained canonical transformation:

- Gambini and Pullin, *Consistent discretization and canonical classical and
  quantum Regge calculus*, <https://arxiv.org/abs/gr-qc/0511096>.

Curved Regge solutions need not retain exact lapse gauge; constraints can
become background-dependent pseudo-constraints:

- Bahr and Dittrich, *(Broken) Gauge Symmetries and Constraints in Regge
  Calculus*, <https://arxiv.org/abs/0905.1670>.

The 600-cell has previously been evolved with Sorkin-type schemes, including
dust cosmology and enlarged variable sets:

- De Felice and Fabri, <https://arxiv.org/abs/gr-qc/0009093>;
- De Felice and Fabri, <https://arxiv.org/abs/gr-qc/0106077>.

Therefore neither iterating a Regge step nor finding a lapse-like variable
from discrete consistency is externally new.

## 3. CONTROL

Before any second-tick solve, the new verifier must:

1. load the committed first-tick artifact and require its exact passing
   outcome and `7/7` status;
2. reconstruct the first tick with a generalized lower/upper homothetic
   evaluator and reproduce all 35 residuals and both boundary momenta;
3. load, validate and apply the vertex-derived old-to-final orbit map rather
   than rely on its observed identity form;
4. require a complete, finite and parity-consistent 30-component post-momentum
   target;
5. retain the Lorentzian branch and 100-decimal arithmetic.

Failure of a control forbids evaluation of a scientific outcome.

## 4. OPEN difference

No located primary source supplies this exact second canonical step on this
carrier with this action, fixed mass and inherited numerical boundary datum.
The following remain open before calculation:

1. whether the extrapolated seed lies on the Lorentzian branch;
2. whether the two-equation Jacobian retains calibrated rank two or its weak
   singular direction collapses toward gauge;
3. whether a connected deterministic Newton correction reaches a root;
4. whether the root passes all 65 equations;
5. whether both staircase parities give the same endpoint;
6. whether the second scale increment continues contraction, stalls or turns
   toward expansion;
7. whether the lapse is stable, changes, or becomes numerically undetermined;
8. whether further iteration and refinement converge.

External novelty of the exact numerical result remains **OPEN** pending a
dedicated review.

## 5. Framing attack

One first tick does not yet define an evolution law.  A second accepted tick
would demonstrate local iterability on the homogeneous subspace, but would
still not establish continuum dynamics, general relativity, inflation or
emergent time.

In particular, the initial absolute `tau0` remains supplied by hand.  The
calculation can at most select subsequent *relative* lapses from canonical
consistency.  If the weak singular value decreases under iteration, the
apparent selection may be a finite-carrier pseudo-constraint tending back to
gauge; that is a negative outcome, not a clock.
