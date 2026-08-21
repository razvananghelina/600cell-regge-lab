# Prior-art gate: on-shell effective H4 boundary Hessian across schedules

Date: 2026-08-21

Status: completed before evaluating any new second derivative.

## Exact object and complete hypotheses

Use the static Lorentzian product over

```text
K0=P(sd K_600),  f=(2640,17040,28800,14400),
tau0=0.0102,
```

with the exact projected rank geometry, all 24 colour-ordered staircase
triangulations, the repository's corrected Regge action and boundary terms,
and the independently selected rank masses `m_r` obtained from
`m_v=K_v/(8*pi)`.  This point is stationary in all ten internal H4 variables
for every schedule.

Let

```text
b = 12 old/new boundary log squared-edge orbit coordinates,
i =  6 cross-diagonal + 4 lapse log coordinates.
```

For each schedule, form the total-orbit action Hessian at the same labelled
boundary data,

```text
H = [[H_bb,H_bi],
     [H_ib,H_ii]].
```

If and only if `H_ii` is invertible within a certified numerical envelope,
eliminate the linearized internal equations and define the Hessian of
Hamilton's principal function

```text
K_eff = H_bb-H_bi H_ii^(-1) H_ib.                (1)
```

No pseudoinverse, schedule average, selected colour order, continuum target,
mode identification, speed, `G` or Planck scale is allowed.

## The matter Hessian is fixed, not fitted

Evaluate the gravitational Hessian with the dust mass disabled, then add the
already selected dust term analytically.  In logarithmic lapse-square
coordinate `z_r=log(rho_r)`,

```text
S_dust=-8*pi*sum_r m_r sqrt(rho_r),
dS_dust/dz_r=-4*pi*m_r*sqrt(rho_r),
d2S_dust/dz_r2=-2*pi*m_r*sqrt(rho_r).             (2)
```

At the static product `sqrt(rho_r)=tau0`; all mixed and boundary dust
derivatives vanish.  Equation (2), including its factor and sign, is frozen
before the calculation.  The old conditional `P1` distribution is forbidden.

## KNOWN from primary literature

- Hamilton's principal function generates the canonical simplicial evolution
  from boundary data after the bulk equations are imposed: Dittrich and
  Hoehn, *Canonical simplicial gravity*,
  [arXiv:1108.1974](https://arxiv.org/abs/1108.1974), DOI
  `10.1088/0264-9381/29/11/115009`.
- Linearized Regge dynamics and its Hessian can depend on bulk triangulation
  in four dimensions; triangulation independence is a substantive condition,
  not automatic: Dittrich and Steinhaus,
  [arXiv:1110.6866](https://arxiv.org/abs/1110.6866), DOI
  `10.1103/PhysRevD.85.044032`.
- Curved Regge backgrounds generically break exact discrete gauge symmetry
  and replace continuum constraints by background-dependent
  pseudo-constraints: Bahr and Dittrich,
  [arXiv:0905.1670](https://arxiv.org/abs/0905.1670), DOI
  `10.1088/0264-9381/26/22/225011`.
- Improved/perfect actions are the established route for recovering
  discretization-independent continuum dynamics; a bare finite Regge action
  need not have that property: Bahr and Dittrich,
  [arXiv:0907.4323](https://arxiv.org/abs/0907.4323), DOI
  `10.1103/PhysRevD.80.124030`.

No located source evaluates (1) for the projected 600-cell carrier, its 24
staircase schedules or the curvature-selected matter branch.  Search absence
does not prove novelty; external novelty remains **OPEN**.

## What the repository already knows, and why it is insufficient

- **DERIVED:** the 24 schedules use 24 distinct sets of temporal diagonals.
- **STRUCTURAL:** at the inherited *off-shell P1* fill, the ten-dimensional
  internal Hessians form 12 time-reversal classes.
- **DERIVED:** at that same point, eliminating only the six cross variables
  gives coincident four-dimensional lapse Schur complements.
- **DERIVED, later:** replacing P1 by curvature-selected masses makes the
  static point internally on shell and selects one schedule-independent
  first boundary derivative.

None of these statements computes (1).  Raw internal Hessian disagreement
does not imply effective boundary disagreement, and first-derivative
agreement does not imply equality of quadratic response.

## Framing attack

This calculation has a deliberately narrow meaning.

1. The twelve boundary variables are only the `H4`-invariant rank-orbit
   sector.  Even exact schedule independence here would not establish
   nonhomogeneous graviton propagation or determine `c`.
2. The background has equal spatial boundaries.  Its Hessian is a linear
   response about a static on-shell product, not a demonstrated nonzero tick.
3. A disagreement among schedules would close the **bare unselected
   staircase construction** as a canonical quadratic evolution.  It would
   not refute Regge calculus, improved/perfect actions, or a theory that
   supplies an independent schedule sum/selection rule.
4. A singular `H_ii` cannot be repaired by choosing an arbitrary
   pseudoinverse.  Gauge-independent elimination would require a separately
   proved null-space condition; absent that proof, the outcome is unresolved.
5. Time reversal exchanges old and new boundary variables.  A schedule and
   its reversed order must be compared after this fixed layer-swap
   congruence.  No other rank-colour permutation is licensed because the four
   ranks have inequivalent geometry and matter weights.

## OPEN and exact next gate

Preregister a target-free high-precision construction of all 24 full `22x22`
Hessians, add (2), certify every internal block, compute (1), and enumerate
the complete multiset of effective matrices before interpreting it.

- One class after the fixed time-reversal identification advances only the
  invariant quadratic sector.
- More than one class is a **DERIVED NEGATIVE** for canonical bare staircase
  evolution on this carrier.
- Any unresolved internal null direction leaves the result **OPEN**.

Only after a single invariant class survives may a larger nonhomogeneous
operator be built and compared with the spatial Laplacian.  No dispersion
relation or speed is to be inferred from the present `12x12` test.

