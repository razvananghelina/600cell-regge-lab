# Prior-art gate: out-of-sample fourth canonical dust tick

Date: 2026-08-16

## Exact object and hypotheses

The object is the fourth iterate of the same two-variable homothetic canonical
map used in the committed weak-lapse calculation.  The carrier is the fixed
600-cell staircase, the Lorentzian Regge action and dust term are unchanged,
the dust mass is conserved, and the order-24 even/odd schedules remain
independent controls.

For each

```text
lambda in {1/2, 1/4, 1/8},
```

the lower fourth-slab state and incoming momentum are the committed third
output at that same `lambda` and parity.  The new unknowns are only the fourth
absolute log scale `a4` and absolute relative-lapse log `r4`.  No fixed-lapse
target and no `lambda=1` fourth state are admitted.

The out-of-sample prediction implied by the committed three-step weak-lapse
jet is

```text
u4/u1           -> 4,
a4/u1           -> 10,
v4/v1           -> 7,
r4/v1           -> 16,
p_post,4/k      -> 9,
```

where `u_n=a_n-a_(n-1)`, `v_n=r_n-r_(n-1)` and
`k=lambda*epsilon3*L0*tau0/4`.  Tick four was not used to identify any of
these integers.

## Primary prior art

- Gambini and Pullin derive consistent Regge dynamics as a canonical
  transformation: <https://arxiv.org/abs/gr-qc/0511096>.
- Dittrich and Hoehn formulate action-generated simplicial evolution in terms
  of pre/post data and consistency conditions:
  <https://arxiv.org/abs/1108.1974>.
- De Felice and Fabri give multi-step dust evolution of the 600-cell:
  <https://arxiv.org/abs/gr-qc/0009093> and
  <https://arxiv.org/abs/gr-qc/0106077>.
- Variational error analysis relates a discrete Lagrangian's approximation
  order to the order of its variational integrator:
  <https://arxiv.org/abs/1102.2685> and
  <https://arxiv.org/abs/1609.02309>.

The literature establishes the broad framework and makes an odd-momentum /
triangular-position law unsurprising for a constant leading force.  The search
did not locate the five specific fourth-step limits above for this fixed
600-cell staircase action.  External novelty is **OPEN**; search is not proof.

## KNOWN / CONTROL / OPEN

- **KNOWN:** action-generated discrete canonical maps and successive 600-cell
  dust evolution exist in the literature.
- **DERIVED COMPUTATIONAL:** the first three iterates on this branch have a
  quadratic weak-lapse limit with momentum ratios `(3,5,7)`, scale increments
  `(1,2,3)`, lapse increments `(1,3,5)`, triangular cumulative scale logs and
  square cumulative lapse logs.
- **CONTROL:** the committed weak-lapse artifact has SHA-256
  `500be1c4e2d7ec4104b9773bc1cfc71065c9d930607119eb616367d18fa5d8f9`
  and records `tick4_target_parsed=false`.
- **CONTROL:** the fourth incoming momentum is freshly taken from the
  committed third post-momentum and the independently derived orbit map.
- **CONTROL:** both schedule parities and all complete action/branch/Jacobian
  gates remain mandatory.
- **OPEN:** whether the five fourth-step observables converge quadratically to
  `(4,10,7,16,9)`.
- **OPEN:** exact all-order recurrence, spatial refinement, anisotropic
  stability, continuum Einstein convergence and emergent time.

## Framing attack

One successful fourth step would validate extrapolation in the iteration
index `n`, which the three-step weak-lapse calculation did not test.  It would
not prove the recurrence for all `n`: all four iterates are still local in
time and live on one homogeneous fixed carrier.

The three `lambda` values are not three independent physical universes; they
are a convergence control.  A hit establishes the leading temporal jet of
this discrete map through `n=4`, not a spatial continuum limit.  It also does
not select an absolute lapse or make time emergent.

Only the new targeted verifier will be run.  The full suite will not be run.
