# STEP 3: Whitney/Hopf comparison with `a_1=5`

Date: 2026-08-10

## Provenance

STEP 1 was committed before this comparison:

- preregistration commit: `9884c95`;
- message: `STEP 1 preregister Whitney Hopf spectra -- NO target comparison`;
- blind artifact SHA-256:
  `7052c9c12c817a40ad463aafe80bcdca32d2b487fac87abdb43151e28c045f64`.

The committed JSON contained the full raw and generalized spectra for all six
Hopf fibrations before this file compared any ratio with the bootstrap seed.

## Result

For every one of the six fibrations, the preregistered data give

```text
raw gap ratio cross/fiber          = 5.0
generalized Whitney ratio          = 5.0
fiber generalized gap multiplicity = 4
cross generalized gap multiplicity = 4
```

The two four-dimensional gap eigenspaces coincide exactly to numerical
precision in the Whitney metric: all four principal cosines are one.  On that
common subspace,

`K_cross=5 K_fiber`.

Moreover,

`dim ker(K_cross-5 K_fiber)=9`,

matching the already-derived vertex `Box` kernel dimension.

## Exact mechanism

Every regular tetrahedral facet has the same volume and stiffness.  Assembly
therefore gives

`K_W=w L_full`,

and the Hopf support split gives

`K_fiber=w L_fiber`, `K_cross=w L_cross`.

The consistent zero-form mass has the symmetry-forced form

`M_0=2V I+(V/4)A=(V/4)(20I-L_full)`.

On the common gap subspace the graph eigenvalues are

`lambda_f=2-phi`,

`lambda_c=5(2-phi)`,

`lambda_full=6(2-phi)`.

Consequently `M_0` is one scalar on that subspace and cancels from the ratio.
The generalized Whitney ratio is therefore exactly five, not a numerical
coincidence.

## Hostile interpretation

This is a **DERIVED PRESERVATION RESULT**, not an independent second selection
of `a_1=5`.  The regular-facet Whitney stiffness is proportional to the graph
Laplacian, so the old Hopf gap ratio is inherited.  The nontrivial content is
that the consistent geometric mass does not destroy the ratio and that the
same four-dimensional mode realizes both gaps across all fibrations.

It does not yet derive a physical speed:

1. both component forms are positive spatial kinetic forms;
2. the Hopf fiber has not been identified with Lorentzian time;
3. the ratio is dimensionless and fixes no conversion to seconds;
4. no refinement-mode extension has yet been certified.

Thus the present status of `c^2=a_1=5` advances from an isolated graph reading
to a **STRUCTURAL DYNAMICAL COMPATIBILITY** with the Whitney kinetic metric,
not to a complete physical derivation.

## New local datum

The blind enumeration found that every one of the 600 coarse tetrahedra has
exactly one fiber edge, for all six fibrations.  This is potentially the
missing refinement rule: it selects one local kinetic direction per
tetrahedron without choosing a new fine edge by hand.

The next gate is to convert that single edge into a positive local tensor,
restrict the tensor to all 24 barycentric children, and compute the new-mode
fiber/cross spectra.  If the ratio is not stable, the physical `c^2=5` route
is closed even though the base-level compatibility is exact.
