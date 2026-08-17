# Result: the weak lapse stiffness is a gravity--dust near cancellation

Date: 2026-08-17  
Targeted verifier: `7/7`; no full-suite run.  
Artifact SHA-256:
`caf44bad3b5677a3ab92c49fb4e736597c55861e9627e0cfae11d94ebd8d27ce`.

## Provenance and honesty boundary

- prior-art gate: `475f534`;
- disclosed post-result protocol: `a0a5491`;
- registered verifier before its execution: `ba61c87`;
- machine result: `937f67a`.

The scalar appearance of the pole Schur matrix and a rough cancellation ratio
were already seen before this verifier.  Consequently this is not a blind
prediction.  The affine dust decomposition is **DERIVED ALGEBRAIC**; scalarity
and near cancellation remain **PATTERN**.

## Exact decomposition

The uniform dust action contributes only to the 120 pole-pole entries.  In
logarithmic pole coordinates,

```text
h_dust = -(2 pi M/120) sqrt(rho).
```

It does not alter the strong block or the strong-pole couplings.  Therefore
the Schur family under a formal dust multiplier `mu` is exactly

```text
S(mu) = S_gravity + mu h_dust I_120.
```

This statement is independent of the observed small eigenvalue.

## Numerical result

**DERIVED COMPUTATIONAL / PATTERN:** all fourteen minimal `2T` blocks, covering
120 full directions in each parity, are consistent with one common scalar
under the stored uncertainty.  The largest non-scalar deviation is only
`1.312 epsilon`, versus the frozen acceptance boundary `10 epsilon`.

At the physical dust multiplier `mu=1`,

```text
gravity Schur scalar = +5.448612752365180e-3
dust Schur scalar    = -5.448616996926997e-3
------------------------------------------------
total residual       = -4.244561817270933e-9
```

Thus

```text
|residual| / max(|gravity|,|dust|) = 7.790163668e-7.
```

The affine zero occurs at

```text
mu_star = 0.9999992209836331,
1-mu_star = 7.790163669e-7.
```

The fixed physical mass is therefore about `0.779 ppm` above the mass that
would make this local pole Schur stiffness vanish at the accepted geometry.

## What this means

**DERIVED:** the previously mysterious `4.24456e-9` is not an arbitrary
regularizer and not an unresolved numerical zero.  It is the residual of two
separately order-`5.45e-3` terms with opposite signs.

**STRUCTURAL / PATTERN:** this is exactly the form expected when a continuum
lapse constraint is almost, but not exactly, restored on a curved discrete
background.  It is consistent with a Regge pseudo-constraint.

**OPEN:** the calculation does not show why the committed mass normalization
lands within `0.779 ppm` of the zero.  It may follow from the same static
Regge balance used to define the mass, from the small non-static deformation,
or from a carrier-specific cancellation.  Refinement or an analytic incidence
derivation is required to distinguish these.

It does not derive a dust clock, a Hamiltonian constraint, gauge symmetry,
gravitons, `c`, an absolute tick or Planck units.

## Consequence for the next calculation

The pole block is regular but extremely soft for a now-understood reason.
Sectorwise Schur elimination is therefore preferable to a raw binary64 inverse
when constructing the complete `1440 x 1440` boundary tangent map.  The next
physical gates are symplecticity, stability and geometric mode content; only
after those pass is a dispersion comparison meaningful.
