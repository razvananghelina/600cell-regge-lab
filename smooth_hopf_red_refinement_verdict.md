# Shape-controlled smooth Hopf tower: comparison and verdict

Date: 2026-08-10

## Provenance and incidental correction

- preregistered definition: `6b8a80b`;
- source-coordinate tolerance correction, before mode inspection: `cde7ad2`;
- blind result: `390e38b`;
- blind JSON SHA-256:
  `efc09b7b54b987b53d8b0a3086f3846020888362d2f12b7f33c12cda7d0485e0`.

The first blind execution exited `9/10` solely because the original
`build_600cell()` coordinates have norm residual `2.80e-11`, while the guard
used `3e-15`.  The correction changed the audit tolerance to `5e-10`, printed
the residual, and did not touch the mesh, matrices, spectra, comparison gates,
or output schema.  It was committed before rerunning and before opening any
mode values.

## Topology and shape

All three meshes are closed Euler-zero tetrahedral complexes and every face
has incidence two:

```text
level     V       E       F       T      chi
0       120     720    1200     600       0
1       840    5640    9600    4800       0
2      6480   44880   76800   38400       0
```

The preregistered shape data are

```text
level     h_max       q_min       q_median     min |P_t X|
0       0.618034     1.000000      1.000000       1.000000
1       0.459506     0.857143      0.927839       0.996718
2       0.232846     0.851101      0.860163       0.999129
```

Both refinement levels exceed `q_min=0.5`, and

`q_min(2)/q_min(1)=0.99295 > 0.8`.

Maximum chord length decreases strictly.  The Hopf field remains far above
the preregistered `0.98` projection threshold.  The stiffness split residuals
are at most `2.42e-16`.  Every topology, shape, field and algebra gate passes.

This is a two-level finite result.  It does not prove a level-independent
shape bound for the infinite tower.

## Charged coordinate modes: target `(1,2,3)`

The maximum Ritz error in each four-dimensional coordinate space is

```text
operator       level 0       level 1       level 2
vertical       0.129420      0.035187      0.008907
horizontal     0.258840      0.070356      0.017812
full           0.388260      0.105036      0.026623
```

All three errors decrease strictly at both steps.

## Fiber-invariant base modes: target `(0,8,8)`

```text
operator       level 0       level 1       level 2
vertical       0.211146      0.073372      0.018673
horizontal     1.788854      0.510559      0.129929
full           2.000000      0.583693      0.148555
```

All three errors also decrease strictly at both steps.  In particular, the
three explicit base modes are moving into the analytically known infinite
vertical kernel rather than remaining at a finite lattice artifact.

Between levels 1 and 2, `h_max` decreases by a factor `1.973`, while all six
mode errors decrease by factors between approximately `3.929` and `3.950`.
This is numerically consistent with `O(h^2)` eigenvalue convergence.  Because
only one fine-level ratio is available and no order estimate was an acceptance
gate, this is **PATTERN**, not a convergence theorem.

## Low spectra

The full operator has exactly one resolved zero at every level.  Its first
four positive modes move

`3.38826 -> approximately 3.104 -> approximately 3.02644`,

toward the exact multiplicity-four value `3`.

The separated vertical low spectrum begins with a constant followed by bands
of three and five modes.  Their upper edges decrease as

```text
three-mode band: 0.21115 -> 0.06007 -> 0.01566
five-mode band:  0.82918 -> 0.48819 -> 0.13736
```

This is the expected beginning `1,3,5,...` of the base-harmonic tower, but the
band identification was not an acceptance criterion.  It remains **PATTERN**.
The load-bearing evidence is the preregistered three-dimensional Hopf-map
coordinate space.

## Verdict

- **DERIVED POSITIVE:** the projected `1 -> 8` refinement passes every
  preregistered gate through `6480` vertices and `38400` tetrahedra.
- **DERIVED POSITIVE:** the true rank-one/rank-two Hopf split reproduces the
  canonical charged and base modes increasingly accurately at both steps.
- **DERIVED:** the combined operator retains only the constant zero mode while
  explicit base modes move toward the vertical zero sector.
- **PATTERN:** the last-step rate is consistent with second order and the low
  vertical bands display `1,3,5,...`.
- **OPEN:** uniform infinite-level shape regularity and a full convergence
  proof.
- **DERIVED NEGATIVE:** none of this restores a continuum role for the old
  factor five or selects a non-round kinetic coefficient.

This is now a credible numerical geometry result rather than a one-mesh
coincidence.  Its physical content remains limited: it validates the round
Hopf carrier, not matter, Lorentzian time, inertia, or a value of `r`.

