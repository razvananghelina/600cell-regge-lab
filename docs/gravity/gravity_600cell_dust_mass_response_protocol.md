# Preregistered protocol: conserved inhomogeneous dust-mass response

Date: 2026-08-18

Prior-art gate commit: `6b587ea`

Status: frozen before constructing or evaluating any mass-response matrix.
The known `119+1` tangent count and the previously observed proximity of the
expanding branch to the canonical weak Schur lift are disclosed in advance.
This is a confirmatory, target-disclosed test, not a blind discovery search.

## 1. Frozen inputs and exclusions

Require exact hashes and passing outcomes for:

- the accepted first dust tick;
- the complete full-rank pre-Legendre artifact and its source;
- the complete 1,440-dimensional boundary tangent JSON and NPZ;
- the full 120-dimensional pole-Schur artifact;
- the complete internal-curvature response artifact and source;
- the direct 600-cell slab geometry source.

Use both frozen five-stage schedules and all four derivative variants:

```text
operational primary, operational shadow,
validation primary, validation shadow.
```

No continuum harmonic, desired mode count, wave speed, Planck scale, particle
mass or Standard-Model target is loaded.  No singular vector may be used to
choose the source carrier.

## 2. Matter parameter and analytic source

At each vertex world-line write

```text
m_v = (M/120) (1 + eta_v).
```

The background is `eta_v=0`.  The dust action remains exactly

```text
S_dust = -8 pi sum_v m_v sqrt(rho_v).
```

For the logarithmic pole magnitude `z_v=log rho_v`, its mixed derivative is

```text
b = partial^2 S_dust / partial z_v partial eta_v
  = -4 pi (M/120) sqrt(rho_v).
```

Construct `B_m` literally from pole-edge incidence: one entry `b` in the
internal pole-equation row for the same world-line and zero everywhere else.
The old-momentum rows have no direct mass derivative.  Require exactly 120
nonzero entries and rank 120 before symmetry reduction.

The new-boundary momentum has no direct mass term because the declared dust
action depends only on the pole magnitudes.  This is a property of the fixed
point-particle discretization, not a claim about general dust.

## 3. Forced canonical response

For each minimal binary-tetrahedral sector of irrep dimension `d`, reconstruct
the full action Hessian block and form the already certified canonical matrix

```text
J_d : (delta x_internal[35d], delta q_new[30d])
      -> (delta E_internal[35d], delta p_old[30d]).
```

Project the five free pole orbits into the same sector to obtain a source
matrix `B_d` with `5d` columns.  Solve with complex Flint balls at 80 decimal
digits:

```text
Y_d = -J_d^-1 B_d.
```

Form the outgoing response

```text
delta q_new = Y_d,new,
delta p_new = K_nx Y_d,internal + K_nn Y_d,new,
```

and apply the already frozen final-to-next boundary-orbit permutation.  Call
the resulting `60d x 5d` ball matrix `R_d`.

The rank of `Y_d` is a mandatory control but carries no physical evidence:
it follows algebraically from regular `J_d` and injective `B_d`.

## 4. Conserved zero-sum source

Every nontrivial irreducible sector is automatically orthogonal to the global
constant mass vector, so retain all `5d` source columns there.  In the unique
trivial `d=1` sector, restrict the five orbit weights to

```text
eta_1 + eta_2 + eta_3 + eta_4 + eta_5 = 0
```

using the deterministic QR orthonormalization of the four columns
`e_i-e_5`, `i=1,...,4`.  The total restored source dimension is fixed before
calculation:

```text
4 + sum_nontrivial 5 d^2 = 119.
```

This equality with the already known 119 strong tangent pairs is disclosed;
it is not evidence.  The calculation must decide the subspace relation.

## 5. Rank and curvature gates

For every schedule, sector and derivative variant:

1. classify the singular values of the zero-sum outgoing phase response;
2. embed `Y_d` as the full slab-edge variation
   `(delta q_old=0, delta x_internal, delta q_new)`;
3. apply the independently reconstructed Jacobian of all 3,840 internal
   causal Regge deficits;
4. classify the singular values of that curvature response.

Use the same global hierarchy as the upstream precision calculations.  For an
operational singular triplet, set `epsilon` to the sum of:

- operational/shadow and validation/shadow step differences;
- operational/validation difference;
- propagated Flint-ball radius;
- direct SVD residual plus the `gesdd`/`gesvd` discrepancy;
- `1e-70` arithmetic floor.

Classify a singular value as:

```text
NONZERO_RESOLVED  if s > 100 epsilon,
ZERO_CONSISTENT   if every variant is < 10 epsilon,
OPEN              otherwise.
```

Report all minimal-sector counts and restore them with irrep multiplicity.
No rank threshold may be changed after execution.

## 6. Target-disclosed tangent comparison

Load each frozen tangent midpoint and radius.  Select invariant spaces by
ordered modulus alone:

- in nontrivial sectors, the `5d` largest-modulus and `5d`
  smallest-modulus Schur spaces;
- in the trivial sector, the four largest and four smallest spaces, leaving
  the previously known near-unit reciprocal pair unselected.

Require a mechanical spectral gap above `2` between each selected and
unselected set.  Compare the zero-sum mass-response image separately with the
expanding and contracting spaces.  For every comparison compute the maximum
principal-angle sine/projector distance.

Its uncertainty is the sum of:

- all four derivative-variant distance changes;
- the propagated mass-response ball radius divided by its smallest singular
  value;
- the tangent-ball radius times the eigenvector condition divided by spectral
  separation;
- the ordered-Schur versus direct-eigenvector subspace discrepancy;
- `10 eps_machine` times the largest participating condition number;
- `1e-70`.

Assign exactly:

```text
IDENTIFIED       if distance <= 10 epsilon,
SEPARATED        if distance > 100 epsilon,
NUMERICALLY_OPEN otherwise.
```

The complete look-elsewhere ledger is fixed:

```text
N = 2 schedules x 7 sectors x 2 branches = 28 comparisons.
```

Report the full hit fractions, not only a favorable sector.  Combining the
two branches or choosing a different count after seeing the spectrum is
forbidden.

## 7. Controls and falsification attempts

Mandatory controls:

1. all input hashes and upstream outcomes match;
2. the direct geometry import retains all 43 certificates;
3. the seven minimal sectors have dimensions `1,1,1,2,2,2,3`, exhaust the
   full regular representation and have high-precision basis residual below
   `1e-70`;
4. all Lorentzian branch and Hessian-reality gates pass;
5. every `J_d` determinant ball excludes zero;
6. current `J_d` singular spectra reproduce the frozen rank artifact;
7. the source has exactly one incidence-selected entry per pole world-line;
8. direct differentiation of the analytic dust term at relative-mass steps
   `1e-20` and `3e-20` reproduces `b` within the high-precision arithmetic
   floor;
9. a zero source gives exactly zero response;
10. the uniform trivial-sector source is retained as a reported control but
    excluded from the physical zero-sum comparison;
11. QR and null-space constructions of the trivial zero-sum domain give the
    same image;
12. Schur-ordered and direct eigenvector tangent spaces agree within their
    reported numerical uncertainty;
13. even/odd response singular spectra receive a calibrated schedule label.

Falsification explicitly varies derivative steps, schedule parity, SVD
driver, tangent-space construction and the basis of the zero-sum source.  It
does not vary the physical source definition.

## 8. Outcome hierarchy

After controls:

- `DUST_MASS_RESPONSE_EXPANDING_IDENTIFIED` only if all 14 expanding
  comparisons are `IDENTIFIED`, all 14 contracting comparisons are
  `SEPARATED`, and both phase and curvature maps have restored rank 119;
- `DUST_MASS_RESPONSE_CONTRACTING_IDENTIFIED` under the reversed branch
  condition;
- `DUST_MASS_RESPONSE_BOTH_BRANCHES_SEPARATED` if all 28 comparisons are
  `SEPARATED` while both response ranks are 119;
- `DUST_MASS_RESPONSE_MIXED_OR_OPEN` for every controlled partial/open
  comparison or rank result;
- `DUST_MASS_RESPONSE_CONTROL_FAILED` if a mandatory control fails.

Even an identified branch remains **STRUCTURAL** as a source sensitivity on
the frozen finite carrier.  It is not a propagating dust degree of freedom or
a graviton.

## 9. Required independent audit

Before consolidation, a separate registered verifier must avoid the decisive
sectorwise Flint solve.  It will use the frozen binary representative matrices
as exact dyadic inputs or reconstruct a literal full-carrier response, use a
different complement/projector algorithm, and include synthetic identified
and separated controls.  Disagreement leaves the result **OPEN** under Rule 4.
