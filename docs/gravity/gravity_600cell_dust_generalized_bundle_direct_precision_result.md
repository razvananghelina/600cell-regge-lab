# Result: direct generalized bundle is common only within certified error

Date: 2026-08-18

## Headline

The direct high-precision reconstruction preserves the rank-`15`
Hermitian-definite generalized-mode bundle in each of the two disclosed
symmetry sectors and at both centered times.  Under the frozen `10/100`
classification, every old/shifted comparison is common and every local and
cross-time `Gamma/Omega` leakage is zero-consistent:

```text
16/16 projector comparisons     GENERALIZED_COMMON_FIBER_RESOLVED
64/64 local leakages             LEAKAGE_ZERO_CONSISTENT
64/64 cross-time leakages        LEAKAGE_ZERO_CONSISTENT
mechanical outcome               DIRECT_GENERALIZED_COMMON_BUNDLE_RESOLVED
```

The scientifically careful statement is weaker than the mechanical label.
The observed old/shifted projector distance is about `3.960493e-7`, while the
complete comparison error is about `8.886344e-4`.  The calculation therefore
certifies consistency with one common bundle and excludes no sufficiently
small rotation inside that error.  It does **not** prove analytic equality of
the two fibers.

## Provenance and reproducibility

| stage | commit |
|---|---|
| primary-literature and framing gate | `6fdddc6` |
| preregistered direct-precision protocol | `d177762` |
| registered verifier before first scientific execution | `c125d84` |
| immutable first-result artifact | `6901ebd` |
| disclosed diagnostic-only protocol amendment | `e7db568` |
| diagnostic artifact | `67d9fa6` |

The first-result artifact has SHA-256

```text
a7f6f915b9284905ad1931131edaa5cd2402dd3b13d1161be12e4201252641a7.
```

The diagnostic artifact has SHA-256

```text
8ded406366dbf291da02dfbf995c4e37036cc6ce745d9240d14905664ba6042a.
```

The amended verifier source has SHA-256

```text
01479fcaa7e5354ea3bb72306ac8cd433a87b11a539f912075d69273a014b510.
```

The first result was executed twice byte-identically.  After the diagnostic
amendment, the verifier was again executed twice byte-identically and passed
`13/13` both times.  Removing only the newly disclosed diagnostic fields from
the amended artifact reproduces the complete first scientific payload
exactly.  No operator, projector, threshold, label or outcome changed.  Only
this targeted verifier was run; the full suite was not run.

## What was computed

For each of

```text
2 times * 2 parities * 2 sectors * 4 derivative variants = 32 cells,
```

the three slab Hessians were reconstructed directly at high precision rather
than loaded through the old binary tangent serialization.  On the fixed
25-dimensional shape carrier, the generalized Hermitian-definite problem was

```text
A_t = -V_t,       B_t = -M_t > 0,
A_t x = lambda B_t x.
```

Every cell resolved `15` negative and `10` positive generalized modes.  The
rank-`15` Euclidean projectors were compared under the literal edge identity;
no Procrustes alignment, polar fit or post-result transport was introduced.
Each projector was then tested against both `Gamma=M^-1 N` and `Omega=M^-1 V`
at its own time and at the other time.

All six slab reconstructions pass, all `48` raw Hessian symmetry defects lie
inside the four-variant derivative-family variation, and all `160` direct
matrix balls overlap the broader serialized controls.  The direct stiffness
radius improves by factors from approximately `1.86e7` to `1.97e9`.

## Error attribution

The post-first-result amendment decomposed, without changing the calculation,
the projector bound into a shape-carrier term, a generalized-eigenspace term,
and a kinetic-metric term.  Across all `32` cells:

| contribution | range |
|---|---:|
| shape-carrier term | `2.24e-9 ... 9.08e-9` |
| generalized-eigenspace term | `2.43e-4 ... 6.46e-4` |
| kinetic-metric term | `1.96e-8 ... 7.88e-8` |
| complete projector error | `2.43e-4 ... 6.46e-4` |

The generalized-eigenspace term is the largest contribution in all `32/32`
cells.  The geometry/shape basis and the kinetic metric are not the present
precision bottlenecks.  The dominant bound comes from propagating the pencil
error across the approximately `3.998e-5` generalized spectral gap.

This is useful negative localization: rebuilding the incidence geometry or
arbitrarily raising arithmetic precision would attack the wrong term.  The
next calculation must use a preregistered residual-based invariant-subspace
bound for the already reconstructed Hermitian-definite pencil, or else leave
the tiny-rotation question open.

## Status ledger

- **DERIVED COMPUTATIONAL, conditional on the frozen derivative family:** all
  `32` generalized pencils have inertia `15 negative + 10 positive`.
- **DERIVED COMPUTATIONAL:** all `16` old/shifted projector comparisons are
  common under the repository's frozen certified-error classification.
- **DERIVED COMPUTATIONAL:** all `128` local and cross-time closure residuals
  are zero-consistent.
- **DERIVED COMPUTATIONAL:** the conservative projector uncertainty is
  dominated in every cell by the generalized-eigenspace perturbation term.
- **STRUCTURAL:** the literal edge identity is the only transport used in this
  audit.  The result says it is compatible with the data; it does not derive a
  physical bundle connection.
- **OPEN:** analytic equality versus a projector rotation smaller than about
  `8.9e-4` in the two-time comparison.
- **OPEN:** a reduced non-autonomous propagator, its characteristic roots,
  long-time stability and continuum/refinement behavior.
- **OPEN:** graviton, wave, particle-inertia, mass or limiting-speed
  interpretations.

Here “negative” and “inertia” refer to the signature of the quadratic Regge
action.  They do not yet mean particle mass or ordinary inertial resistance.

## Scientific meaning and next gate

This is a real step toward dynamics: the same action-selected rank-`30`
subsystem (two rank-`15` sectors) is compatible with two consecutive centered
ticks and is closed, within certified error, under the matrices that enter the
linear recurrence.  It is not yet a propagation law.

The next gate is to tighten only the dominant generalized-eigenspace bound.
If a residual-based bound falls below the observed `3.96e-7` displacement,
the two fibers can be classified as resolved-rotated or genuinely common
under the same frozen thresholds.  If it cannot, identity versus tiny rotation
remains **OPEN**, and any reduced product must carry that uncertainty rather
than silently choosing a common basis.
