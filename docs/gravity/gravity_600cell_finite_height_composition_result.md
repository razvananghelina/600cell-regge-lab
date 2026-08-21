# Result: finite-height composition is not a unique evolution map

Date: 2026-08-21.

## Provenance

```text
composition prior-art gate                         c9eb996
composition protocol frozen before calculation    386466e
primary verifier registered before first run       c7bffc3
first symbolic timeout preserved                   1fe35f4
local junction-certificate correction              ac6527a
primary nonuniqueness artifact                     0f832c6
adversarial protocol frozen                        e8628bc
adversarial verifier registered before first run   34fc405
first adversarial comparison failure preserved     c1371a6
serialization-comparison correction frozen         25c4fed
corrected comparison implementation                0764442
accepted adversarial artifact                      90cb5b7
```

Targeted verifiers:

```text
reproducible/verify_gravity_600cell_finite_height_composition.py
  10/10 PASS

reproducible/verify_gravity_600cell_finite_height_composition_adversarial.py
  9/9 PASS
```

Accepted artifacts:

```text
reproducible/gravity_600cell_finite_height_composition.json
SHA-256 d4e36141863bd2ae515b96eeeff4f50eb087016cca8cfb6f4b1e3355d6fba447

reproducible/gravity_600cell_finite_height_composition_adversarial.json
SHA-256 d50e87f736e51585596aa1d7778238febaf7422840d668499878d8bd917f99e9
```

Both targeted verifiers were rerun together and reproduced these artifacts
byte for byte.  No full-suite run was performed.  The static registry audit
found `407/407` distinct registrations, two deliberate exclusions and zero
duplicates, unregistered verifiers or missing files.

## Headline

```text
FINITE_HEIGHT_TWO_SLAB_NONUNIQUE_ADVERSARIALLY_CORROBORATED
```

> **DERIVED COMPUTATIONAL NEGATIVE FOR UNIQUENESS, REPRESENTATIVE-SCOPED /
> STRUCTURAL OPEN SELECTION / ADVERSARIALLY CORROBORATED:** under the frozen
> homogeneous 600-cell Regge-plus-conserved-dust action and canonical
> momentum-matching convention, an admitted first slab at `v=3/2` has at
> least two distinct isolated physical second slabs.  Therefore the action
> and the stated physical inequalities do not define a single-valued
> two-slab evolution map.

This does not negate the exact one-slab existence theorem.  It shows that its
finite state-dependent update is not, by itself, a deterministic tick.  A
further branch-selection principle may exist, but none is derived here.

## Complete hypotheses

Use the previously certified positive-Lorentzian homogeneous cellular slab
at zero cosmological constant, with conserved global dust and the exact
finite-height first-slab branch.  Define canonical boundary momenta by

```text
p_pre  = -(L_minus/2) * partial S/partial L_minus,
p_post =  (L_plus /2) * partial S/partial L_plus.
```

At the shared slice retain the physical mass and normalize the outgoing
canonical data as

```text
m1  = M/L1,
pi1 = p_post,1/L1^2.
```

The second slab is required to have positive proper height, positive next
scale and a nonzero two-equation Jacobian.  The calculation does not cover
nonhomogeneous edge data, another Regge action, refinement, a cosmological
constant, a different matter clock or an additional branch-selection axiom.

## Canonical composition equations

The complete action obeys

```text
S(lambda*L_minus,lambda*L_plus,lambda^2*rho;lambda*M)
  = lambda^2*S(L_minus,L_plus,rho;M).
```

Consequently the outgoing momentum has degree two.  Direct differentiation
of two joined slabs gives the shared-slice equation

```text
(L1/2) * partial(S1+S2)/partial L1 = p_post,1-p_pre,2 = 0.
```

For normalized second-slab variables `(h2,q2)`, the two residuals reduce to

```text
C2 = 8*pi[mu(q2)-m1] + 4*pi*h2*q2*mu(q2),
P2 = p(q2)-pi1       - 2*pi*h2*mu(q2).
```

Their scalar elimination is

```text
E2 = 4*pi[mu(q2)-m1] + q2[p(q2)-pi1].
```

The primary route used this elimination and monotone intervals to count all
real roots.  The adversarial route instead redifferentiated the full action
and solved the original two equations directly, without using `E2` or reading
the primary roots until construction was complete.

## Frozen representative census

The representative inputs were fixed before evaluation:

| Incoming `v` | All real second roots | Physical second roots | Result |
|---:|---:|---:|---|
| `3/2` | 3 | 2 | nonunique continuation |
| `3` | 3 | 2 | nonunique continuation |
| `20` | 1 | 0 | no forward continuation |

For the load-bearing witness `v=3/2`, the first slab produces

```text
L1  = 2.96260125838050816101450426467548...,
m1  = 3.96692448160445187605754289398116...,
pi1 = -172.883601601732563658264525133347....
```

The two physical second slabs are

```text
branch A:
  q2       = 0.0212594337109692113846320849415857...,
  h2       = 7.28930126260954098473864773192187...,
  L2/L1    = 1.15496641699173171194831360136404...;

branch B:
  q2       = 31.2792236208252636198289982042168...,
  h2       = 0.0689338063134504452159117927059578...,
  L2/L1    = 3.15619594271307285222389331068446....
```

Both have positive height, positive scale ratio and nonzero Jacobian.  The
third solution is the time-reversed algebraic branch and has negative height.

## Adversarial confirmation

The independent verifier:

- reconstructed and differentiated the complete action;
- solved in `(h2,q2)` at 80, 120 and 180 decimal digits;
- recovered both roots from all eight preregistered `+-5%` seed
  perturbations;
- checked the direct constraint, momentum and shared-slice residuals;
- rejected the reverse root by its negative height;
- showed that the wrong momentum scaling, reversed sign and mass reset each
  fail distinctly;
- compared with the primary artifact only after both roots existed.

The first adversarial execution passed every physical test but failed an
impossible `1e-70` comparison against values stored to only 60 significant
digits.  That failure was preserved before a correction was frozen.  The
corrected check compares the independent roots with the actual serialized
precision; every printed 60-digit value agrees exactly.  No root, equation,
seed, physical inequality or outcome rule was changed.

## Scope of the negative result

| Claim | Status |
|---|---|
| The certified one-slab finite update exists | **DERIVED EXACT** |
| The fixed action defines a unique second slab at every admitted state | **DERIVED NEGATIVE** |
| At `v=3/2` there are two isolated physical continuations | **DERIVED COMPUTATIONAL / ADVERSARIALLY CORROBORATED** |
| At `v=3` there are two physical continuations | **DERIVED COMPUTATIONAL** |
| At `v=20` there is no physical continuation | **DERIVED COMPUTATIONAL** |
| The complete first-state composition domain is classified | **OPEN** |
| A canonical physical principle selects one branch | **OPEN** |
| The finite update is a deterministic or fundamental tick | **NOT DERIVED** |
| Nonhomogeneous data remove the ambiguity | **OPEN** |
| Refinement removes the ambiguity | **OPEN** |
| `c`, `G` or Planck time follow | **NOT DERIVED** |
| External novelty of the coefficient-level result | **OPEN** |

The negative statement needs only one valid state with two physical outputs;
therefore the incomplete global first-state classification does not weaken
the refutation of uniqueness.  Conversely, the representative census cannot
support the stronger claim that no deterministic refinement or additional
physical selector exists.

## Interpretation

The present homogeneous action supplies a **relation**, not a single-valued
flow: some incoming states have multiple allowed successors and at least one
tested state has none.  Calling the isolated one-slab solution a `tick` would
therefore be too strong.  The scientifically honest next question is not to
pick the preferred numerical branch after seeing it, but to identify and
preregister an independently motivated selector, or to move to the
nonhomogeneous canonical system and test whether the ambiguity is lifted.
