# Result: no same-state analytic half-step at the dust turning point

Date: 2026-08-21

## Provenance

```text
prior-art gate                                  4f4b5c2
frozen primary protocol                        a025951
registered primary implementation              55cb3a2
primary result and adversarial protocol         177232b
registered adversarial implementation           18f69cd
preserved first adversarial control failure     1b0a024
first exact-reduction correction                19cc044
preserved second adversarial control failure    e548944
final complex-branch correction                 7a1dff3
```

Targeted verifiers:

```text
reproducible/verify_gravity_600cell_relational_step_composition.py
  8/8 PASS

reproducible/verify_gravity_600cell_relational_step_composition_adversarial.py
  9/9 PASS
```

Accepted adversarial artifact:

```text
reproducible/gravity_600cell_relational_step_composition_adversarial.json
SHA-256 aec6012991d0a4c86a237a7122244939352bf4aa7e4f3c66ebca0cc49b67d9a5
```

The two failed adversarial artifacts remain in the repository.  No nonlinear
root search and no full suite were run.  The exact leading branch count made
the nonlinear stage inapplicable under the frozen protocol.

## Headline

```text
SAME_STATE_HALF_STEP_ABSENCE_ADVERSARIALLY_CORROBORATED
```

> **DERIVED NEGATIVE, scoped and adversarially corroborated:** on the exact
> homogeneous cellular 600-cell Regge-plus-conserved-dust action, the
> time-symmetric weak-lapse state family has no analytic first half-step with
> the same incoming canonical state, actual proper duration asymptotic to
> half the coarse duration, `log L=O(e^2)` and all lapse equations retained.

The old scaled-lapse family does possess such a branch only because it also
halves the incoming momentum.  It therefore changes the initial physical
state and cannot serve as a temporal-refinement test.

## Exact obstruction

Let

```text
epsilon = 2*pi-5*acos(1/3),
D       = 5*sqrt(2)/3-epsilon > 0,
p0(e)   = 180*epsilon*e.
```

For the coarse ansatz

```text
L_plus=exp(A e^2+O(e^4)),
rho=e^2 exp(O(e^2)),
```

the leading lapse and same-state momentum equations factor as

```text
45 A(D A+4 epsilon)=0,
90   (D A+4 epsilon)=0.
```

They share the accepted contracting root

```text
A_coarse=-4 epsilon/D
        =-12 epsilon/(5 sqrt(2)-3 epsilon),
```

which independently reproduces the blind cellular weak-lapse theorem.

For a nominal half-step from the literal same state,

```text
rho=(e^2/4) exp(O(e^2)),
P_pre=p0(e),
```

the first-slab leading equations instead become

```text
90 A(D A+epsilon)=0,
180  (D A+3 epsilon/2)=0.
```

The lapse roots are

```text
A=0,  -epsilon/D,
```

whereas the momentum equation has the sole root

```text
A=-3 epsilon/(2D).
```

There is no common root.  The exact resultant is nonzero:

```text
-729000 epsilon^2(3 epsilon-5 sqrt(2)) != 0.
```

Because the first fine slab is absent already at leading order, no second
fine slab or nonlinear series seed exists in the preregistered branch class.
A numerical root search after this result would be an unregistered search
for a different branch, not a continuation of the failed one.

## Why the old half-lapse trajectory survives

If the fine incoming momentum is changed from `p0(e)` to its own static value

```text
p0(e)/2=90*epsilon*e,
```

the fine momentum equation becomes

```text
180(D A+epsilon)=0.
```

It now shares the nonzero lapse root `A=-epsilon/D`.  Thus the obstruction is
not a sign error or a verifier unable to find a fine branch.  The branch is
restored exactly by the one change that invalidates a same-state comparison.

This explains the earlier successful `lambda` hierarchy: each value of
`lambda` starts from a different canonical momentum.  Its quadratic integer
law remains a valid family of trajectories, but it is not evidence that one
coarse evolution step equals two fine steps.

## Adversarial route

The primary verifier expanded `S/e` and then differentiated its exact finite
jet.  The adversarial verifier used the opposite order:

1. differentiate the complete unexpanded cellular action;
2. substitute the half-step path;
3. take the exact rescaled limits;
4. eliminate the scale coefficient only afterward.

It recovered the same two polynomials and nonzero resultant.  Direct
100-decimal evaluations on the literal exponential path gave quadratic
convergence for all three defining zero residuals.  At the two incompatible
candidate roots, the nonzero momentum and lapse obstructions dominated their
finest halving drifts by respectively

```text
563740.6 and 155072.9.
```

The first two adversarial runs returned `7/9` because generic SymPy
simplification did not reduce two exact complex-branch identities.  Both
failures and hashes were preserved.  Explicit polynomial reduction and
`expand_complex` repaired only those control booleans; no equation, root,
precision point, tolerance or outcome changed.

## Physical status ledger

| Claim | Status |
|---|---|
| Existing coarse weak-lapse root | **DERIVED** |
| Same-state analytic half-step in the registered turning-point class | **DERIVED NEGATIVE / ADVERSARIALLY CORROBORATED** |
| Old `lambda` family is a same-state refinement | **REFUTED** |
| Finite homogeneous lapse acts as a state-dependent pseudo-constraint | **STRUCTURAL, strongly supported** |
| One carrier slab is a fundamental indivisible tick | **NOT DERIVED** |
| Absolute tick from the current classical theory | **DERIVED NEGATIVE under the scale-covariance hypotheses** |
| Generic nonzero-velocity half-step | **OPEN** |
| Nonanalytic or different-carrier half-step | **OPEN** |
| Perfect/improved action restores refinement | **OPEN** |
| Explicit Brown--Kuchar relational clock in this discrete model | **OPEN** |

## Framing verdict

This result does **not** promote the surviving coarse lapse to a fundamental
time quantum.  A discretization can fix a multiplier because continuum gauge
symmetry is broken; that is the standard pseudo-constraint mechanism.  A
fundamental tick would require a separately justified principle selecting the
carrier and forbidding temporal refinement.  The repository has no such
principle.

Combined with the global scale theorem, the honest current statement is:

> The finite homogeneous action selects a relative lapse from canonical data,
> but its turning-point map is not locally divisible into two analytic
> same-state half-steps, and no absolute duration is selected.

## Post-result prior-art audit

The learned terms `same-state half-step`, `lapse pseudo-constraint`,
`temporal refinement` and `perfect action` recover the established mechanism:

- [Gambini--Pullin](https://arxiv.org/abs/gr-qc/0511096) explain how
  consistent discretization can determine continuum multiplier-like
  variables and yield a constraint-free canonical transformation.
- [Bahr--Dittrich](https://arxiv.org/abs/0907.4323) construct improved/perfect
  actions by solving refined equations subject to coarse boundary data and
  connect refinement invariance to restored gauge symmetry.
- [Dittrich--Hoehn](https://arxiv.org/abs/1108.1974) provide the pre/post
  canonical simplicial framework.
- [Brown--Kuchar](https://arxiv.org/abs/gr-qc/9409001) show what an explicit
  dust clock requires: dust proper time and its conjugate momentum as phase
  variables.

No located source prints the present two polynomial factors for this
600-cell action.  That negative search is not evidence of novelty.
**External novelty remains OPEN.**

## Next load-bearing test

The result is only at the time-symmetric scaling `p0=O(e)` and
`log L=O(e^2)`.  Before modifying the theory, test the generic-velocity
scaling

```text
p0=O(1),
log L_plus=V e+A e^2+...,
rho=(e^2/4) exp(...),
```

with the same incoming `(L0,p0,M)` in coarse and fine histories.  Enumerate
all leading branches before comparing endpoints.

- If a unique same-state fine branch exists and composes, the present no-go
  is a turning-point staggering effect.
- If it also fails, the current cellular action has no local relational
  temporal refinement in either natural velocity class; the next honest
  route is an explicit dust-clock/deparameterized or improved/perfect action,
  not another fitted lapse search.

