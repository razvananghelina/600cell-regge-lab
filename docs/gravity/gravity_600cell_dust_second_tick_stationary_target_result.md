# Stationary roots versus the second canonical target

Date: 2026-08-16

## Provenance

- target-independent roots committed first: `caaf1f1`;
- prior-art gate and mass-framing audit: `fcc4d7c`;
- frozen comparison protocol: `9127b04`;
- pre-evaluation artifact-control correction: `28c4dd1`;
- registered verifier before evaluation: `3eb77bb`;
- targeted verifier:
  `reproducible/verify_gravity_600cell_dust_second_tick_stationary_target.py`;
- artifact:
  `reproducible/gravity_600cell_dust_second_tick_stationary_target.json`.

Only this targeted verifier was run.  It returned **5/5**.  The full suite was
not run.

## Mechanical verdict

```text
STATIONARY_SECOND_TICK_NO_HIT
hits = 0/2
```

Both preregistered roots were compared on all 30 components and both derived
schedule parities with the canonically mapped post-momentum of the accepted
first tick.

| committed root | structural label | residual per component | norm | certified bound | hit |
|---|---:|---:|---:|---:|---:|
| 0 | contracting | `+1.6161034538904e-9` | `8.8517631695780e-9` | `3.6513653962011e-21` | no |
| 1 | time reversal | `-5.4485975915103e-3` | `2.9843198076385e-2` | `3.6513653962011e-21` | no |

The contracting residual is uniform but about `2.42e12` times the inherited
norm bound.  It is close in scale and completely nonzero by the frozen gate;
it is not an approximate pass.

## What the mass observation changes

**DERIVED:** the action keeps

```text
M=(90/pi)*epsilon3*L0
```

fixed at the original scale.  It does not impose `M proportional to L` anew
on each evolving boundary.  The exact all-lapse cancellation was proved only
for the static product family `L_lower=L_upper=L0`.

**DERIVED:** when the two boundary scales differ while `M` is conserved, the
pole/lapse equation ceases to be an identity.  That breaking already supplied
one of the two equations that selected the accepted first tick.  This is the
interesting branch of the proposed dichotomy, not the static-only branch.

**DERIVED NEGATIVE AT FIXED INHERITED LAPSE:** mass conservation does not make
the contracting stationary slab glue when its lapse is held exactly at the
first-tick value.  A nonzero lapse correction is required if the second tick
exists locally.

**OPEN:** mass conservation has not supplied an absolute clock.  The initial
`tau0=0.0102` is still external.  The established mechanism selects relative
scale/lapse changes through the joint pole equation and canonical seam.

## Next falsification test

Use the contracting root chosen by geometry (`b<a1`), not by target proximity,
as the sole seed of a separately preregistered two-variable correction in
absolute `(b,r)`.  Solve simultaneously

```text
G(b,r)=0,
P_pre(b,r)=mapped P_post,first.
```

The stationary enumeration already certifies a resolved nonzero determinant
there.  Acceptance still requires all 35 internal equations, all 30 momentum
components, both schedule parities and the original Lorentzian branch gates.
Failure remains a genuine local negative, not permission to select another
seed.
