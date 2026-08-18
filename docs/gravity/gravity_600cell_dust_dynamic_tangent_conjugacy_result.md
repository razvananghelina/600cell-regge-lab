# Geometric conjugacy audit of the dynamic tangent maps

Date: 2026-08-17

## Verdict

> **DERIVED COMPUTATIONAL / STRUCTURAL:** the two non-isomorphic ordered slabs
> induce the same full 60-dimensional linear canonical boundary map, within the
> frozen calibration, after the preregistered identification of identical
> physical boundary-edge orbits.  All 60 carrier-derived direct boundary lifts
> pass and all 60 time-reversed lifts fail.

Frozen combined outcome:

```text
DYNAMIC_TANGENT_BOUNDARY_COVARIANCE_ONLY
```

This explains the previously certified characteristic isospectrality at the
level of the finite boundary tangent maps.  It is not a complete-slab
isomorphism and not general triangulation independence.

## Provenance and commit ordering

- prior-art gate: `0129053`
- protocol, before candidate enumeration: `e006b8c`
- enumeration verifier registered before first run: `0efa93a`
- target-blind candidate artifact committed: `858512c`
- comparison verifier registered only after that commit: `72629dd`

Enumeration artifact SHA-256:

```text
51b52457eba84ca1e41926b6e4fb1c51032f788b70bde916a3fb755d0323cb3e
```

Comparison artifact:
`reproducible/gravity_600cell_dust_dynamic_tangent_conjugacy.json`

Comparison artifact SHA-256:

```text
be8919a33029a6aa65c2393e117bc36be6afb0b193a76ac317994c18752ddeae
```

The enumeration passed `7/7`; the comparison passed `8/8`.  Only these targeted
verifiers were run.  The full suite was not run.

## Stage A: the geometry before the matrices

All 14,400 `H4` actions were enumerated without loading a tangent matrix or
spectrum.  The committed counts were:

```text
complete direct even-to-odd slab candidates       0
complete time-reversed slab candidates             0
distinct boundary permutations                    60
supporting boundary H4 actions                  1,440.
```

Thus the strongest explanation is **DERIVED NEGATIVE**: the even and odd
four-dimensional triangulated slabs are not related by an `H4` carrier
isomorphism, either directly or by layer reversal.

The unique identification of identical physical edge sets was present before
comparison.  It is one of the 60 boundary maps and is supported by 24 `H4`
actions.  The remaining maps enlarge the look-elsewhere denominator; none was
selected after seeing a residual.

## Stage B: full canonical-map identities

For every boundary permutation `Q`, the verifier tested exactly

```text
direct:   T_odd C(Q) - C(Q) T_even,
reversed: T_odd K(Q) T_even - K(Q),

C(Q)=diag(Q,Q),
K(Q)=diag(Q,-Q).
```

The complete result is

```text
direct boundary PASS     60/60
direct boundary OPEN      0/60
direct boundary FAIL      0/60

reversed boundary PASS    0/60
reversed boundary OPEN    0/60
reversed boundary FAIL   60/60.
```

For the preregistered identical-physical-edge map, the direct identity has

```text
Frobenius residual       5.0850204375e-29
calibrated uncertainty   2.6174930701e-16
residual / uncertainty   1.9427063611e-13.
```

Across all 60 direct maps, residuals range from `5.085e-29` to `2.122e-28`,
and every ratio is below `8.108e-13`.  The reversed residual is
`8.5225628816e7`, against uncertainty `1.1766615368e-10`, a resolved failure by
ratio `7.243e17` for every candidate.

## What the positive result means

- **DERIVED COMPUTATIONAL:** after eliminating all 35 internal variations and
  passing to pre/post canonical boundary variables, the schedule-dependent
  difference seen in the unreduced static Hessians disappears at linear order
  about the accepted dynamic homothetic tick.
- **DERIVED COMPUTATIONAL:** the physical-edge identification gives a genuine
  similarity of the two finite tangent maps within calibration.  Their shared
  spectrum is therefore explained, not merely observed.
- **STRUCTURAL:** this is boundary covariance between two different bulk
  triangulations.  There is no complete carrier isomorphism connecting the
  bulks.
- **PATTERN:** all 60 direct carrier variants pass and all 60 reversed variants
  fail.  The other 59 direct passes are symmetry-related attempts, not 59
  independent discoveries; the load-bearing hit is the separately fixed
  physical-edge identification.
- **OPEN:** whether the equality is exact, whether it persists away from the
  homogeneous background, and whether it survives refinement or the full
  720-edge phase space.

The result is stronger than equality of eigenvalues but narrower than a
perfect discretization.  It concerns one derivative of one finite canonical
map at one highly symmetric, dust-filled solution.  It does not establish
nonlinear schedule independence, lattice gravitons, a causal cone, a limiting
speed or a physical duration for the tick.

## Post-result primary-source audit

Dittrich and Steinhaus study linearized Regge actions under changes of
triangulation and find triangulation independence to be special rather than
automatic, especially in four dimensions:
<https://arxiv.org/abs/1110.6866>.  Dittrich, Kaminski and Steinhaus further
show strong restrictions and non-locality for 4D discretization-independent
measures: <https://arxiv.org/abs/1404.5288>.  Therefore the present cancellation
after internal elimination is not licensed by a general theorem and remains a
specific finite result.

No located primary source reports the same two dust-filled 600-cell schedule
maps or this boundary covariance.  External novelty is **OPEN**, not proved by
search.

## Next decisive test

The positive equality is measured at the maximally homogeneous dynamic point.
The cheapest falsifier is nonlinear and anisotropic: use the already resolved
canonical inversion to evolve the same preregistered small zero-sum boundary
perturbations through both schedules, identify the physical edge orbits, and
compare the complete outputs.  If the difference begins at quadratic order,
the present result is a symmetry-point cancellation.  If the nonlinear maps
continue to agree across a frozen amplitude family, there is evidence for a
genuine effective boundary discretization independence.
