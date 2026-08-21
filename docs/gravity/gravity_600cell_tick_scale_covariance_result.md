# Classical tick scale covariance: consolidated result

Date: 2026-08-21

## Accepted provenance

- prior-art gate: `3275575`;
- original frozen protocol: `fd1f8a8`;
- primary implementation before evaluation: `127a4af`;
- primary result: `0ac0aba` (`12/12`);
- first adversarial implementation: `b94d7ec`;
- retained first disagreement: `9951c7a` (`10/12`);
- retained precision disagreement: `af862ab` (`7/10`);
- exact flag-incidence diagnosis: `0fa8947` (`9/12`, zero incidence mismatches);
- corrected-state protocol: `8e0f8a6`;
- corrected implementation before evaluation: `4561579`;
- accepted artifact:
  `reproducible/gravity_600cell_tick_scale_covariance_state_correction.json`;
- accepted artifact SHA-256:
  `5a607577a6db295fd1af03a1e8693053d8abd86a864cc1375b9cf4bba8855ac7`.

The accepted targeted verifier returned **21/21**.  No full-suite or deferred
nonlinear census was run.

## Headline

```text
ABSOLUTE_CLASSICAL_TICK_NO_GO_ADVERSARIALLY_CORROBORATED
```

> **DERIVED EXACT / ADVERSARIALLY CORROBORATED:** on the fixed 600-cell
> staircase carrier, the zero-cosmological-constant classical Lorentzian
> Regge-plus-dust action is globally scale covariant when all geometrized dust
> masses scale with the geometry.  Its stationary equations and canonical data
> occur in continuous global scale families.  Under these complete hypotheses,
> the theory cannot select an absolute nonzero tick.

## Exact theorem and hypotheses

For every positive `alpha`, scale every signed squared edge length and dust
proper-length square by `alpha^2`, and every geometrized dust mass by `alpha`:

```text
q_e   -> alpha^2 q_e,
rho_v -> alpha^2 rho_v,
m_v   -> alpha m_v.
```

Assume flat-simplex four-dimensional Regge gravity with boundary terms,
pressureless point dust, zero cosmological constant, and no fixed external
length, mass, or quantum scale.

A simplex Gram matrix has squared-length degree one.  Dihedral angles and
deficits are invariant; every triangular hinge area has length degree two.
The Regge action therefore has length degree two.  Each dust term
`m_v sqrt(rho_v)` has the same degree under the simultaneous mass scaling.
Thus

```text
S(alpha^2 q, alpha m) = alpha^2 S(q,m).
```

Derivatives with respect to logarithmic squared lengths scale by `alpha^2`;
derivatives with respect to raw squared lengths are invariant.  Consequently
zeros of the internal equations map to zeros, while pre/post canonical momenta
transform covariantly.

The theorem covers both current selected mass rules:

```text
M_coarse = (90/pi) epsilon3 L,
mu_r     = K_r/(8*pi),
```

because both have length degree one.

## Computational corroboration

The frozen state is nonhomogeneous, off shell, branch-valid, and has all 35
internal plus all 60 boundary derivative components resolved as nonzero.

The 100-decimal orbit action passed all 95 logarithmic-derivative identities
for both parities and `alpha in {3/5,7/4}`, with maximum errors between about
`1e-93` and `1e-98`.

After correcting the independently diagnosed state-construction error, the
direct 2400-simplex binary64 action and all 95 raw derivatives pass within
conditioning-aware bounds of `1.31e-6`--`3.64e-6`; observed errors are at most
`9.16e-8` for current/final derivatives and `1.96e-10` for old-boundary
derivatives.

The literal direct 2400-simplex action at 80 decimals passes the two scale
identities in both parities with errors below `9e-94`.  Its corrected base
actions agree with both a fresh orbit evaluation and the previously stored
primary values within `5e-56`.

## Falsification controls

Holding the mass fixed while rescaling geometry breaks covariance in every
case.  The normalized action defects are `0.6277` or `1.0`, and the raw pole
derivative defects are resolved.  Thus the positive result is conditional on
mass scaling and the verifier can falsify the theorem when one hypothesis is
removed.

The exact finite incidence audit independently found, per parity,

```text
260 triangle orbits,
100 simplex orbits,
1000 incident-flag orbits,
24000 labelled flags,
0 shortcut/exact coefficient mismatches.
```

Thus the representative orbit reduction is not hiding a flag-multiplicity
fit.

## Retained failure history

The first direct adversarial and precision runs compared different off-shell
states.  The primary correctly perturbed the new boundary by
`exp(1e-6*((i mod 11)-5))`; the adversarial concatenated internal and new
variables and applied the internal modulus-five rule to both.  Both states were
separately group-invariant and scale covariant, so ordinary invariance gates did
not catch the mismatch.

The error was exposed by the exact flag audit and then proved structurally from
the frozen AST.  No failed artifact or threshold was overwritten.  The episode
adds a reusable methodological requirement: adversarial implementations must
compare fully expanded labelled input states, not only carrier maps and summary
parameters.

## Physical status ledger

| Claim | Status |
|---|---|
| Complete classical action has global scale covariance | **DERIVED EXACT / ADVERSARIALLY CORROBORATED** |
| Current equations select an absolute nonzero tick | **DERIVED NEGATIVE under stated hypotheses** |
| `tau=0.0102` is derived | **REFUTED; it is inherited input** |
| Canonical consistency can select a next-lapse ratio | **DERIVED COMPUTATIONAL LOCAL** |
| `tau_next/tau0=0.999998220375...` is an absolute clock | **REFUTED framing; it is relative** |
| A fixed external mass can break the symmetry | **DERIVED by hostile control** |
| Such a fixed mass is itself derived by current theory | **OPEN / unsupported** |
| Dust can provide relational proper time | **KNOWN STRUCTURE; model realization still partial** |
| `tau/L` can be selected dynamically | **OPEN** |
| Seconds, `c`, `G`, `hbar`, Planck time or Planck mass are derived | **OPEN / unsupported** |

## What is required for a physical tick

The classical geometry may still derive a dimensionless clock rule, such as
`tau/L` or a relative lapse recurrence.  To convert that rule into an absolute
time, at least one independently justified dimensionful structure must break
the scale family.  Examples include a fixed physical matter mass, a
cosmological scale, a boundary scale, or a genuinely quantum scale.  Supplying
one after inspecting the desired answer is fitting, not derivation.

The next valid question is therefore not another search for a bare value of
`tau`.  It is whether the canonical dust dynamics selects a nonzero,
refinement-stable **dimensionless** `tau/L` without inheriting `tau0`.  If it
does not, time remains relational/gauge in this classical sector.  If it does,
the separate scale-origin problem remains.

