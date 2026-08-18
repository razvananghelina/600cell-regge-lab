# Result: the complete third-slab tangent is regular and canonical

Date: 2026-08-18

## Headline

**DERIVED COMPUTATIONAL, BLIND.**  The complete action-derived tangent of the
third accepted non-static dust slab exists on the full `1440`-dimensional
position/momentum carrier:

```text
T_3 : Phase(slice 2) -> Phase(slice 3).
```

Every pre-Legendre determinant ball excludes zero, every minimal block of
`T_3` is canonical, and the rigorous shifted product

```text
C_32 = T_3 T_2
```

is canonical.  All fourteen primary comparisons between the independently
derived schedule parities are robust.  The preregistered outcome is

```text
LATER_SLAB_TANGENT_CERTIFIED.
```

This upgrades the project from a tangent known only over its first two slabs
to a verified later-background derivative.  It is not yet a graviton,
instability or Lyapunov result.

## Provenance and reproduction

| stage | commit |
|---|---|
| primary-literature and framing gate | `920ce5d` |
| blind protocol | `7e4f47f` |
| verifier registered before first execution | `c45cd2e` |

Verifier:

```text
reproducible/verify_gravity_600cell_dust_later_slab_tangent.py
```

Artifacts:

```text
reproducible/gravity_600cell_dust_later_slab_tangent.json
SHA-256 58a95d90d569b25a3aa396346f5198472c8aed706846b4057182b04c9f7480c4

reproducible/gravity_600cell_dust_later_slab_tangent.npz
SHA-256 77b4dd54a5dcba9d1aa12870b361c9d7d7dde11ccaaa558361b9c1dc24768196
```

The first and second targeted runs both reported `16/16`; both artifacts were
byte-identical.  The full suite was deliberately not run.

## Background and seam

The tangent is evaluated on the committed, unequal later geometry

```text
(a_2,r_2) -> (a_3,r_3),
```

with the original conserved dust mass.  The mass is not recomputed from a
later scale.  The complete thirty-component canonical seam between slabs 2
and 3 has, in both schedules,

```text
maximum residual  3.11959e-41,
inherited bound   3.65137e-21.
```

The old-to-new orbit map is the literal identity on the `30` boundary orbit
labels.  The full carrier retains seven minimal sectors of dimensions
`3,2,2,2,1,1,1`, restoring

```text
sum 60 d^2 = 1440
```

real phase directions.

## Structural controls

For both schedule parities:

- all four derivative variants reproduce exactly `3305` nonzero local
  Hessian entries;
- all `28` third-slab pre-Legendre determinant balls exclude zero;
- all seven `T_3` blocks satisfy the Flint-ball symplectic identities;
- all seven directly multiplied `T_3T_2` blocks satisfy them as well;
- the deterministic archive contains the preregistered `448` arrays.

The largest symplectic midpoint/error ratio across all recorded map blocks is
below `2.68e-5`.  The maps are highly conditioned—the largest tangent
condition number is approximately `1.0631e12`—but the complete propagated
balls, rather than raw binary residuals, decide canonicality.

All schedule comparisons give

```text
ordered singular spectra  14/14 SCHEDULE_ROBUST,
matched eigen spectra      14/14 SCHEDULE_ROBUST.
```

## Blind branch diagnostic

After restoring representation multiplicity, each schedule gives for `T_3`

```text
120 resolved contracting,
120 resolved expanding,
240 total resolved off-unit,
the remaining 1200 unit-consistent or open.
```

The resolved moduli span approximately

```text
0.0210854 ... 47.4262.
```

For `C_32=T_3T_2`, each schedule again gives `120+120`, with resolved moduli

```text
0.00510785 ... 195.777.
```

The different open-versus-unit-consistent counts between schedules are
conservative uncertainty labels for near-unit, nonnormal eigenvalues.  They
do not signal a schedule-dependent spectrum; all preregistered invariant
comparisons are robust.

**STRUCTURAL DIAGNOSTIC:** these counts describe finite canonical maps in the
declared logarithmic coordinates.  They are not a physical instability count
and are not Lyapunov exponents.

## Hostile interpretation audit

1. A regular symplectic tangent is mandatory for a variational theory and is
   not by itself evidence for general relativity.
2. `T_3T_2` contains only two updates on a changing background.  Its
   eigenvalues cannot be promoted to asymptotic growth rates.
3. The large condition numbers make raw eigenvalue labels less trustworthy
   than the ball canonicality and schedule-invariant singular data.
4. No earlier conformal/shape basis, negative eigenspace or spatial spectrum
   was loaded.  Therefore this artifact does not yet test persistence of the
   previously found `30` negative-stiffness directions.
5. Dust perturbations are not independent canonical variables, and no exact
   constraint quotient or physical polarization count has been supplied.

## Status ledger

- **DERIVED COMPUTATIONAL:** the third-slab full tangent exists, is regular
  and is canonical on all `1440` declared phase directions.
- **DERIVED COMPUTATIONAL:** `T_3T_2` is canonical and schedule robust.
- **DERIVED COMPUTATIONAL:** the complete slab-2/slab-3 seam passes by about
  twenty orders of magnitude inside its inherited bound.
- **STRUCTURAL:** the blind finite-map census contains `120` resolved
  reciprocal contracting/expanding pairs.
- **OPEN:** persistence of the conformal/shape split and the `30` old negative
  modes on the shifted Jacobi operator.
- **OPEN:** physical quotient, polarization, long-time growth, dispersion,
  limiting speed, refinement and external novelty.

## Next load-bearing gate

Use only the now committed `T_2,T_3,T_3T_2` balls to reconstruct the shifted
three-slice Jacobi equation

```text
K_-^(2) delta q_1 + K_0^(2) delta q_2 + K_+^(2) delta q_3 = 0.
```

Commit that operator before loading the previous conformal/shape or negative
carriers.  Only the subsequent target-disclosed comparison may decide whether
the same action-selected subsystem persists along the evolving background.

