# Two-slab 600-cell dust gluing: corrected result and status ledger

Date: 2026-08-16

## Provenance

- prior-art gate: `8c45290`;
- corrected consecutive-schedule framing: `620461d`, `6c4a377`;
- frozen two-slab protocol: `29dcfa5`;
- first implementation/result: `1148dbf`;
- post-result branch-precision correction preregistration: `ab75d91`;
- corrected implementation committed before re-evaluation: `037f49b`;
- corrected result artifact: `a766740`;
- corrected JSON SHA-256:
  `a5a22d219b71e49c154c1ef80ed9da93b1aef0b93cd2d6ed22f041b71f62db77`.

The result is in
`reproducible/gravity_600cell_dust_two_slab_gluing.json`; the registered
verifier is
`reproducible/verify_gravity_600cell_dust_two_slab_gluing.py`.

Only this verifier was run.  The full suite was deliberately not rerun at the
user's request.

## Status ledger

| Claim | Status | Evidence / boundary |
|---|---|---|
| Three-layer combinatorial carrier | **DERIVED COMPUTATIONAL** | 360 vertices, 4800 distinct four-simplices and one shared 600-cell |
| Consecutive-slab orbit map `P` | **DERIVED** | logical endpoint incidence; bijective; no permutation search |
| Time-reversal quotient map `R` | **DERIVED** | 24 reversing automorphisms induce one and the same 30-orbit map |
| Lorentzian branch | **DERIVED COMPUTATIONAL** | every arbitrary-precision action evaluation has one timelike Gram direction |
| Direct action gluing | **DERIVED CONTROL** | direct union equals the sum of factors at 61 frozen points |
| Pre/post sign convention | **DERIVED CONTROL** | `p_pre + R p_post = 0` orbitwise |
| Shared variation identity | **DERIVED CONTROL** | `dS_total/dq_shared = p_post - P p_pre` orbitwise |
| Repeated identical sandwich is stationary | **DERIVED NEGATIVE** | shared residual is nonzero by `2.72e19` uncertainty norms |
| Direction of the repeated-sandwich cusp | **DERIVED DIAGNOSTIC** | scale projection contains the full vector within numerical error |
| A next spatial frame | **OPEN / NOT SOLVED** | requires canonical inversion, not action gluing |
| Unique physical tick duration | **OPEN** | temporal edge squares are variables, but the collective lapse/gauge status is unsettled |
| Expansion, contraction or inflation | **OPEN** | no next-frame scale has been solved |
| Full gravitational perturbations | **OPEN** | present carrier retains only the order-24-invariant 30-orbit boundary sector |

## Exact object and result

For each of the two schedule parities, the verifier constructs

```text
S_total(q0,x1,q1,x2,q2)
  = S_sigma(q0,x1,q1) + S_sigma(q1,x2,q2)
```

both as two separately evaluated one-slab factors and as a direct Regge action
on their 4800-simplex union.  The published regular dust sandwich is used on
both factors.  There are exactly two attempts, `even` and `odd`; no phase,
orbit, branch or tolerance search is performed.

Both parities return

```text
TWO_SLAB_GLUING_CONTROL_PASSED
```

and the corrected run reports `25/25` checks.

The common geometry counts are

```text
vertices                         360
four-simplices                  4800
shared tetrahedra                600, each with multiplicity 2
triangles                      11280
triangle orbits                  470
four-simplex orbits              200
```

All direct simplex and triangle orbits have size 24.  The consecutive map
`P` is the identity in the independently sorted logical-endpoint orbit
coordinates.  This is an output of the incidence construction, not an assumed
momentum pairing.  For each parity, exactly 24 vertex automorphisms reverse
the five phase cells and all 24 induce one quotient map `R`.

## Corrected Lorentzian branch certificate

The first result exposed a real audit limitation: the auxiliary branch check
converted the `1e-20` derivative steps to binary64, where `exp(1e-20)` rounds
to one.  The actions and derivatives were still at 100 decimal digits, but
the advertised branch coverage was too strong.  The result was therefore
made provisional and the correction was preregistered before re-evaluation.

The corrected verifier retains a binary64 audit over the 61-point
`+/-1e-6` envelope and additionally applies Jacobi's signature criterion to
every 100-decimal simplex Gram matrix in every action evaluation.  Per parity
it certifies

```text
arbitrary-precision action evaluations        903
evaluated four-simplex orbit representatives 120400
negative Gram directions in every simplex       1
minimum absolute leading principal minor   1.0404e-4
minimum complex angle-argument modulus       0.99534
```

The maximum action imaginary contamination is `5.73e-97` in the even case
and `2.30e-96` in the odd case, below the frozen `1e-70` bound.

## Action and momentum identities

The maximum direct-versus-factor action error over the base and all sixty
shared-orbit audit points is `1.143e-97` in both parities.  All 90 boundary
derivatives per parity pass the two frozen operational/validation pairs.

Orbitwise,

```text
max |p_pre + R p_post|                 8.32e-80  even
                                      6.41e-76  odd

max |dS_total/dq_shared
      - (p_post - P p_pre)|            1.20e-81  both
```

At the published regular sandwich every component is, to the stored
precision,

```text
p_pre,j   = -0.00090810444890653157621005256593...
p_post,j  = +0.00090810444890653157621005256593...
cusp_j    = +0.00181620889781306315242010513186...
```

Thus

```text
||cusp||                 = 0.00994778582473809847980212619916...
||cusp uncertainty||     = 3.6513653962e-22
signal / uncertainty     = 2.7244e19
```

The repeated copy is therefore not a stationary two-slab history.  It exits
the symmetric sandwich and then restarts its time reverse, producing a
momentum cusp.

## Scale versus anisotropic direction

Each of the 30 boundary variables represents an orbit of 24 physical
600-cell edges.  Since all orbit weights are equal, the normalized constant
vector is the natural global scale direction in this reduced carrier.
Projecting the computed cusp gives, in both parities,

```text
mean component                  0.00181620889781306315242010513186...
max component spread            1.713e-39
||cusp_perp|| / ||cusp||        3.205e-37
relative cusp uncertainty       3.671e-20
```

Hence the anisotropic projection is seventeen orders of magnitude below the
already tiny relative uncertainty: the cusp is entirely in the homogeneous
scale direction for this control.

This is **DERIVED DIAGNOSTIC**, not yet a derivation of expansion.  Full
600-cell symmetry already strongly constrains the direction at the regular
background.  The sign and magnitude of the actual change in spatial scale
must be obtained by solving the next slab with

```text
p_pre(next) = p_post(current).
```

A proposed scaling such as `cusp ~ rho*tau^2` cannot be inferred from the
single frozen value `tau=0.0102`; it requires a separate preregistered family.

## Framing attack

The action-gluing and shared-variation identities follow from the variational
composition law once the geometry, hinge constants and signs are correct.
Their successful reproduction is therefore an implementation and canonical-
structure control, not new dynamics.  The nonzero cusp says that repeating
the symmetric sandwich is invalid; it does not by itself say whether the next
frame expands or contracts.

The present `30+35` model is not a one-variable minisuperspace, but it is also
not the full 600-cell edge phase space.  The 720 boundary edges are quotient
to 30 orbits of 24 by the ordered-schedule stabilizer.  Anisotropic modes
outside this invariant sector have not been tested.

The temporal pole squares are varied in the complete one-slab action and have
stationarity equations, so their metric lengths are not merely fixed action
parameters.  Nevertheless the combinatorial existence/order of a slab is
supplied structurally, and this control freezes the poles at the published
`tau`.  Whether the collective lapse is exact gauge, a weak pseudo-constraint
or dynamically selected remains **OPEN**.

## Post-result prior-art check

The technical terms learned in the calculation were searched again after the
result.  Dittrich and Hoehn's [canonical simplicial gravity
formalism](https://arxiv.org/abs/1108.1974) already establishes the use of a
discrete action as Hamilton's principal function, pre/post momenta and
evolution by gluing simplices.  Their earlier [covariant-to-canonical Regge
analysis](https://arxiv.org/abs/0912.1817) explains how nonlinear curvature
can turn exact constraints into point-dependent pseudo-constraints.  General
[variational-integrator constructions](https://arxiv.org/abs/1102.2685) also
make action-generated implicit evolution and momentum preservation known
structure.

No located primary source prints the present 30-orbit maps or this explicit
two-slab 600-cell dust audit.  A search cannot prove absence; external novelty
remains **OPEN**.  The mathematical composition law itself is **KNOWN**, not
a discovery of this project.

## Next preregistered question

Before solving a next frame, compute the complete singular spectrum of the
`65 x 65` canonical-inversion Jacobian at the published point:

```text
unknowns:   35 internal x + 30 final-boundary q_new
equations:  35 internal stationarity equations
          + 30 equations p_pre = p_target
```

The protocol must freeze the following outcome interpretation in advance:

- rank 65: locally unique map; an exact lapse gauge is refuted;
- rank 64 with a geometrically identified collective null vector: one
  lapse/gauge direction;
- additional null vectors only at the symmetric solution: symmetry-frozen
  modes, not automatically gauge;
- small but nonzero singular values: pseudo-constraints, reported as a
  spectrum rather than rounded to zero;
- stable additional nulls with explicit geometric generators: candidate
  constraints/gauge directions.

A generic off-shell perturbation cannot decide physical gauge.  The spectrum
should be recomputed at a genuinely stationary neighboring solution if and
when canonical inversion produces one.  Only after this rank/gauge census
passes should the published pre-momentum be used as a reproduction control,
followed by `p_target=p_post(published)` as the first candidate forward tick.
