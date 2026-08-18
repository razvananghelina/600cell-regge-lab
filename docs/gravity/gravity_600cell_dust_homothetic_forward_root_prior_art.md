# Prior-art gate: non-static homothetic forward root

Date: 2026-08-16

Status: **completed after the frozen local-response result and before any
evaluation inside its newly predicted root bracket**.

Upstream local-response result: `086009a`.

This is a targeted primary-source map, not proof of external novelty.

## 1. Exact object and complete hypotheses

Retain both derived order-24 staircase schedules, the complete Lorentzian
Regge plus fixed-mass dust action, the published lower regular boundary
`L0`, the chosen proper strut length `tau0=0.0102`, and

```text
M = (90/pi)*(2*pi-5*acos(1/3))*L0.
```

On the exact homothetic embedding parameterized by
`s=log(L_+/L0)`, use

```text
q_old       = L0^2,
q_new       = exp(2s)*L0^2,
pole        = tau0^2,
diagonal    = exp(s)*L0^2-tau0^2.
```

The new object is the nonzero negative root predicted after the frozen local
test, together with two independent substitutions:

```text
g_internal(s)[35] = 0,
p_pre(s)[30] = P p_post(static)[30].
```

The orbit map `P` and target momentum must be loaded from the already
committed two-slab gluing artifact.  The root is not accepted from the
globally summed lapse equation alone.

This mission fixes `tau0`; it does not ask whether the theory selects a lapse
or an absolute time unit.

## 2. KNOWN

The broad physical construction is old.  De Felice and Fabri build the
initial time-symmetric 600-cell sandwich, then evolve it with five Sorkin
stages, fixed total dust mass and a chosen timelike lapse.  Their
five-dimensional embedding supplies homothetic trial lengths, and their
subsequent spatial sections contract from maximum size:

- <https://arxiv.org/abs/gr-qc/0009093>;
- <https://arxiv.org/abs/gr-qc/0106077>.

The paper explicitly chooses `tau_k/l_k` approximately constant, rather than
deriving `tau_k`.  It also solves five restricted `4 x 4` systems after
fixing lapse and shift and omitting the corresponding Bianchi-related
equations.  Therefore a contracting first step is **KNOWN**, but it is not
the same object as the present complete `35+30` internal-and-canonical
substitution.

The continuum closed dust solution near its maximum is also a control, not a
new prediction.  In the outer-time embedding,

```text
R(v)=Rmax-v^2/(4 Rmax).
```

If the equal-edge initial sandwich has endpoints at `v=-tau0/2` and
`v=+tau0/2`, the next equal-duration endpoint is `v=3 tau0/2`.  To leading
order this predicts

```text
Delta log L approximately -tau0^2/(2 Rmax^2)
                         approximately -zeta^2*tau0^2/(2 L0^2).
```

This comparison is only a continuum control: the finite Regge root need not
equal it.

Canonical evolution by an action-generated pre/post momentum match is known
structure:

- Dittrich and Hoehn, <https://arxiv.org/abs/1108.1974>.

Global versus local Regge variation is known to be non-equivalent in general
closed FLRW discretizations:

- Liu and Williams, <https://arxiv.org/abs/1501.07614>.

## 3. CONTROL

The root verifier must first reproduce, without reevaluating a search space,
the committed facts:

- the local-response artifact passes `10/10`;
- `dE_lapse/ds` is resolved positive;
- every tested non-static diagonal residual is below `1e-90`;
- the two-slab gluing artifact passes and supplies 30-component pre/post
  momenta and the forward orbit identification;
- the exact homothetic diagonal formula remains unchanged.

The static root `s=0` is excluded explicitly.  A solver that returns it has
failed the mission.

## 4. OPEN difference

The repository has not yet established any of the following:

1. a certified nonzero bracket and root for the pole/lapse equation;
2. vanishing of all 35 local equations at that root;
3. equality of all 30 pre-momenta with the committed forward target;
4. parity independence of the root and canonical data;
5. uniqueness of the nonzero root in a stated bracket;
6. continuation to a second dynamically produced slab;
7. lapse or clock selection.

No located primary source prints this exact root and complete substitution
for the repository carrier.  External novelty remains **OPEN**.

## 5. Framing attack before calculation

Even if the root and junction both pass, this will not yet show that the
theory has a unique evolution map on its full 65-variable carrier.  The
homothetic embedding is an additional symmetry/geometry selection, and the
earlier artificial momentum homotopy followed a different, zero-lapse
connected branch.  The result must therefore be called a **geometrically
selected homogeneous forward root**, not a proof of global uniqueness.

Nor would it supply emergent time.  `tau0` is input.  The possible advance is
more modest and concrete: the first action-stationary, canonically glued
non-static 600-cell slab for one chosen proper interval.
