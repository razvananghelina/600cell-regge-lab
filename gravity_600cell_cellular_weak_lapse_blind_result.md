# Blind result: analytic weak-lapse jet of the cellular 600-cell map

Date: 2026-08-17

Prior-art commit: `b77856a`  
Preregistered protocol commit: `71d10b4`  
Rank-one protocol correction: `35f37d4`  
Registered verifier commit: `80f3164`

Stage-A artifact:

```text
reproducible/gravity_600cell_cellular_weak_lapse_blind.json
SHA-256 6d39e9a4594d9c9ead102f94cf9115d8474132ecce511fe7359826dcc73b9de0
```

Provenance status: **no committed tick artifact was parsed or compared**.
The artifact records `tick_artifacts_parsed=false`.  This note and artifact
are to be committed before Stage B opens the stored weak-lapse results.

## Verdict

**DERIVED:** `CELLULAR_WEAK_LAPSE_JET_DERIVED` (`8/8` targeted checks).

The exact cellular Regge action and canonical seam equations select one real
contracting weak-lapse coefficient branch through four slabs.  At each step
the leading system has rank one and fixes the scale coefficient `A_n`; the
next system in `(B_n,R_n)` has rank two with the same exact determinant

```text
det = 16200 epsilon_3^2 != 0.
```

No integer sequence was assumed in obtaining the solution.

## Exact result

Put

```text
epsilon_3 = 2*pi-5*acos(1/3),
q         = 5*sqrt(2)-3*epsilon_3 > 0,
x         = e^2.
```

For

```text
log L_n       = A_n x+B_n x^2+O(x^3),
log(rho_n/x)  = R_n x+O(x^2),
```

the leading coefficients are

```text
A_1 = -12 epsilon_3/q,
A_n = n(n+1) A_1/2,                  n=1,2,3,4,

R_1 = -10 epsilon_3(7 sqrt(2) epsilon_3+60)/q^3,
R_n = n^2 R_1,                       n=1,2,3,4.
```

Numerically,

```text
A = (-0.230433886199,
     -0.691301658598,
     -1.38260331720,
     -2.30433886199),

R = (-0.263208404079,
     -1.05283361632,
     -2.36887563671,
     -4.21133446526).
```

The exact target-blind normalized sequences printed by the verifier are

```text
(A_n-A_(n-1))/(A_1-A_0) = 1,2,3,4,
A_n/(A_1-A_0)           = 1,3,6,10,

(R_n-R_(n-1))/(R_1-R_0) = 1,3,5,7,
R_n/(R_1-R_0)           = 1,4,9,16,

p_out,n/k                = 3,5,7,9.
```

`B_n` is also derived exactly in the artifact.  It was retained as the
registered nuisance coefficient needed to resolve `R_n`; it was not used as
a sequence target.

## Independent controls

- The closed action independently reproduces the static zero action and
  collective pre/post momentum.
- Symbolic derivatives agree with 100-digit centered finite differences at
  three generic nonstatic points to maximum relative error
  `4.5445182e-80`.
- The calculation discovers rather than assumes that `F/e` first survives
  at `x^1`, whereas the momenta start at `x^0`.
- The exact leading equations vanish on the selected solution; the affine
  next-order equations are solved by the stated nonzero-determinant Cramer
  identity.
- Substitution of the truncated jet into the full, unexpanded equations at
  `e=1/100,1/200,1/400` gives halving orders approaching `7` for the lapse
  residual and `5` for the seam residual for every `n=1..4`.

Under the preregistered volume map and half-step continuum convention,

```text
A_1(discrete)/A_1(closed-FLRW) = 1.0789794680413509...
```

Thus the fixed 600-cell coefficient is about `7.90%` larger in magnitude
than that continuum control.  This is a finite-carrier comparison, not a
spatial-refinement result.

## Preserved implementation history

The first registered blind run found that the nominal leading `2x2` system
was actually rank one and stopped before producing an artifact.  The
protocol was corrected in commit `35f37d4` before continuing.  Later runs
were interrupted only after stack traces showed generic SymPy
`simplify/series/expand` expression swell.  Commits `1c97935`, `4a52e1e`,
`664f65b`, `8da272d`, `a7883da` and `c708413` replace those heuristic
operations by the same exact finite jet, Cramer identity and arithmetic in
`Q(sqrt(2))(epsilon_3)`.  No tick target entered these corrections.

## Meaning and limits

**DERIVED:** the homogeneous cellular action has a unique local contracting
canonical jet through four slabs, and its leading motion is discrete
constant-acceleration kinematics.

**STRUCTURAL:** this is the local minisuperspace gravitational evolution of
the regular 600-cell carrier.  It is the discrete analogue of the
turning-point expansion of closed dust Friedmann cosmology.

**OPEN:** comparison with the previously committed numerical ticks.  That is
Stage B and has deliberately not been performed here.

**OPEN:** convergence under spatial refinement, anisotropic stability,
propagating tensor modes, and whether the selected finite-carrier lapse is a
physical clock rather than a pseudo-constraint.

**OPEN:** an absolute time unit, limiting speed, Planck scales and particle
masses.  The present result supplies a dimensionless local evolution law,
not those dimensional constants.

**KNOWN / not a novelty claim:** regular 600-cell Regge dust cosmology and
quadratic turning-point motion are established in the Collins--Williams
line of work.  The narrow exact reconciliation with this repository's four
canonical ticks is not assessed until Stage B and post-result literature
search.

