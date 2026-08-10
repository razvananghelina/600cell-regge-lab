# Dimension reconciliation: what the three probes actually measure

Date: 2026-08-09

## Verdict

The apparent `4` versus `3` disagreement was mostly a category error.

- **DERIVED:** the carrier used by `D=d+d*` is the three-dimensional
  simplicial boundary of a convex four-dimensional 600-cell.  “Four” here is
  the ambient/polytope dimension; the cochain degrees are `0,1,2,3`.
- **DERIVED:** the four `Box_p` nullities are `(9,13,1,1)`, so their
  alternating count is `-4`.  The `Box_p` do not intertwine the simplicial
  coboundaries: the exact maximum residuals of
  `d_p Box_p - Box_(p+1) d_p` are `(30,44,22)`.  Thus this is neither the
  Euler characteristic nor a Fredholm/cohomological index of that complex.
- **DERIVED:** the Kähler--Dirac operator has two harmonic modes, exactly
  matching `(b0,b1,b2,b3)=(1,0,0,1)`.  Its complete heat flow has a unique
  maximum `3.295663771`; it never reaches four.
- **PATTERN, not a resolved dimension measurement:** a target-free exhaustive
  2% shoulder rule gives `3.274268` over only `0.376` decade for `D^2`.  This
  is consistent with a three-dimensional carrier but fails the repository's
  already-registered half-decade heat-plateau gate.
- **DERIVED NEGATIVE:** the exact numbers `(2640,14880,55920)` are finite
  matrix moments, not Seeley--DeWitt coefficients.  They have no dimensional
  exponent.

There is no audited spectral, heat, Weyl, topological, or index probe in the
repository that supports a four-dimensional **static carrier**.  The exact
integer `4` survives in several arithmetic identities and as the ambient
dimension of the polytope, but neither is a spacetime-dimension measurement.
If the theory retains a `3+1` spacetime, the time/fourth direction is an
additional dynamical structure that has not been derived here.

## Registry repair completed first

Commit `12952c9` adds the missing duplicate-registration guard and removes the
duplicate entry.

At parent `c9205c0` the registry contained `79` entries and `78` distinct
names.  The sole duplicate was
`verify_incidence_operator_enumeration.py`.  The four names cited in the
mission (`verify_alpha_spectral.py`, `verify_gravity.py`,
`verify_oriented_chamber_double.py`, and `verify_rg_bootstrap.py`) each occur
once at that commit.  The old printed denominator `78/79` itself implies one
duplicate, not four.

The repaired coverage guard now exits `2` for any of:

- an on-disk verifier missing from the registry;
- a registered name missing on disk;
- a name registered more than once.

The first run with the truthful unique denominator completed at `77/78` in
`658.3 s`: every unique verifier passed except
`verify_chamber_symmetry_sat.py`, which exceeded the suite's `300 s` timeout.
That is not reported as “all passing.”

Commit `7246e9d` contains the dimension reconstruction, corrected labels and
registered verifier.  Commit `f7adc8c` gives the declared exhaustive chamber
census its measured `600 s` allowance and makes any non-PASS suite result a
nonzero process exit.  On that committed tree the final unique registered
suite completed `79/79 PASS` in `901.0 s`; the chamber census itself completed
in `216.9 s`.  Thus the earlier timeout was an inadequate suite allowance,
not a verifier failure.

## Complete hypotheses and observable definitions

### Static cochain carrier

The vertices are the 120 unit-quaternion vertices of the 600-cell.  Edges,
triangles and tetrahedral facets are reconstructed from their adjacency.  The
verified `f`-vector is

`(120,720,1200,600)`.

The ambient affine span has rank four, while every maximal boundary simplex
has four vertices and hence dimension three.  The boundary of a convex
4-polytope is a 3-sphere; this explicit construction, not the Betti numbers
alone, identifies `S^3`.

### Kähler--Dirac diffusion observable

For the total oriented cochain space,

`D=d+d*`,

`D^2=Delta_0 direct_sum Delta_1 direct_sum Delta_2 direct_sum Delta_3`,

and

`d_s(t)=2t Tr(D^2 exp(-tD^2))/Tr(exp(-tD^2))`.

This is the return-probability exponent of the theory's own operator.  On a
finite matrix it tends to zero at both ends.  Only an intermediate plateau
or a controlled refinement scaling can approximate a continuum dimension.

### Static counting observable

The registered `d_N=3.0688` is the fitted finite-scale Weyl exponent

`d_N(Lambda)=2 d log N_D2(Lambda)/d log Lambda`

on a selected contiguous set of distinct positive `D^2` shells, with shell
multiplicities retained.  It measures a finite spectral-counting slope, not a
topological invariant.  The registered target-conditioned rule finds a
`0.794911`-decade 3D interval with value `3.068762`, and no 4D interval.
Because its positive controls were too small, its geometric interpretation
remains finite-size limited.

### The `Box` alternating nullity

The value called a “spectral index” was defined as

`dim ker(Box_0)-dim ker(Box_1)+dim ker(Box_2)-dim ker(Box_3)`

`=9-13+1-1=-4`.

The four `Box_p` are separately constructed signed Laplacians on the four
cell layers.  They are not the components of `D`, and they do not form a
cochain map for the natural simplicial differential.  Consequently their
alternating nullity has no Euler, Witten, APS, or Fredholm index theorem behind
it.  It is an exact algebraic invariant of that four-operator hierarchy.

The actual simplicial Euler characteristic is independently

`120-720+1200-600=0`,

as required for an odd-dimensional sphere.

## Independent Kähler--Dirac reconstruction

The new verifier does not import the spectrum from
`verify_kahler_dirac.py`.  It reconstructs the cells and integer
coboundaries, then obtains the full `D^2` spectrum from the three singular
spectra.  Every positive singular value occurs in the adjacent Hodge degrees.

The results are:

| quantity | independent result |
|---|---:|
| cochain dimension | `2640` |
| ranks `(d0,d1,d2)` | `(119,601,599)` |
| Betti numbers | `(1,0,0,1)` |
| harmonic modes of `D^2` | `2` |
| smallest positive eigenvalue | `0.145898033750` |
| largest eigenvalue | `15.708203932499` |
| ratio | `107.665631` |
| log window | `2.032077` decades |
| unique heat-flow maximum | `3.295663771` at `t=1.774404055` |

The two harmonic modes are precisely the harmonic 0-form and harmonic
3-form.  This is an independent Hodge consistency check of the boundary
construction.

### Attack on the reported `3.287`

Commit `c9205c0` contains only 31 comment lines with the calibration and
dimension table.  It contains no executable geodesic-sphere construction,
Hasse spectrum, product spectrum, plateau-value function, or checks.  The
exact claim `3.287 over 0.285 decades` therefore has no reproducible algorithm
in that commit.

The function that does exist in `verify_spectral_dimension_flow.py` claims to
find the widest 2% span, but greedily skips overlapping starting points and
returns only a width.  It cannot produce the claimed plateau value.

The registered reconciliation verifier implements the literal rule
exhaustively and without a target:

> Choose the widest contiguous log-time interval with `d_s>0.5` for which
> `(max d_s-min d_s)/mean d_s <= 0.02`; report its mean.

Applied to the same exact `D^2` spectrum, it gives

`3.274268 over 0.376 decades`,

not `3.287 over 0.285 decades`.  Thus the exact colleague number is
**REFUTED AS REPRODUCED PROVENANCE**, although its qualitative conclusion
“near three and not four” survives.  Most importantly, the exact global
maximum `3.295664` proves that no alternative plateau rule can find a 4D
region in this heat flow.

### Calibration and its limitation

The same rule is applied first to successive midpoint-and-radial geodesic
refinements of an icosahedral `S^2`:

| vertices | peak | widest 2% span | reported shoulder |
|---:|---:|---:|---:|
| 12 | 2.0447 | 0.145 | 2.0308 |
| 42 | 2.4155 | 0.173 | 2.3993 |
| 162 | 2.5355 | 0.184 | 2.5180 |
| 642 | 2.5704 | 0.184 | 2.5533 |
| 2562 | 2.5794 | 0.250 | 2.0045 |

The finest control reproduces the known value two and demonstrates why a
peak must not be read as dimension.  But the selector jumps discontinuously
from the peak shoulder to the geometric shoulder only at the last tested
level.  A single successful control level therefore does not make every
sub-half-decade shoulder a trustworthy dimension estimator.  This is why the
Kähler--Dirac value is labelled **PATTERN consistent with three**, while the
absence of four is **DERIVED** from the full maximum.

For comparison, the same fully specified rule gives:

| finite operator | global peak | 2% shoulder | width |
|---|---:|---:|---:|
| 600-cell vertex Laplacian | `3.6323` | `3.6084` | `0.140` |
| barycentric Hasse carrier | `3.5624` | `3.5378` | `0.320` |
| Kähler--Dirac `D^2` | `3.2957` | `3.2743` | `0.376` |
| vertex graph times path `P_320` | `4.4842` | `0.9804` | `0.455` |

The first three move monotonically toward three as more of the actual
cochain geometry enters.  This is a meaningful finite-size pattern, not a
continuum-limit proof.  The product result says that at the late scales where
the fixed spatial factor is saturated, diffusion sees the path alone.  It
closes that finite unrescaled product route; it does not refute every possible
continuum time construction.

## Why the finite “Seeley--DeWitt coefficients” contain no dimension

The exact arithmetic is confirmed:

`c0=Tr(I)=2640`,

`c1=Tr(D^2)=14880`,

`c2=(1/2)Tr(D^4)=55920`,

with reduced triple `(11,62,233)` and

`2*62^2+1=3*11*233`.

For a finite matrix, however,

`Tr exp(-tD^2)=c0-c1 t+c2 t^2+O(t^3)`.

It is analytic at `t=0`; there are no powers `t^((k-d)/2)` with a negative
leading exponent.  Numerically and analytically,

`d_s(10^-8)=1.127e-7 -> 0`.

Therefore:

- **DERIVED:** the moment triple and its Diophantine identity are exact;
- **DERIVED NEGATIVE:** they contain no finite-complex dimension;
- **REJECTED:** calling them Seeley--DeWitt coefficients or assigning them
  cosmological, Einstein--Hilbert and Yang--Mills terms without a controlled
  continuum asymptotic;
- **OPEN:** a refinement family with a proved heat-kernel asymptotic from which
  continuum Seeley--DeWitt coefficients could be extracted.

## Does any remaining repository probe support four?

No, after separating unlike meanings of the integer four:

| occurrence of `4` | what it measures | dimension evidence? |
|---|---|---|
| ambient coordinates / regular 4-polytope | affine dimension of the solid 600-cell | no; its used boundary is 3D |
| `dim_R(C^2/2I)=4` | real dimension of the singular orbifold whose link is the 600-cell boundary | not by itself; identifying that bulk with spacetime is an extra bridge |
| `9-13+1-1=-4` | alternating nullities of non-cochain `Box_p` | no |
| `phi^3+phi'^3=4` | algebraic trace in `Q(sqrt(5))` | no |
| `p=4` or `Tr(D^4)` in downstream formulae | chosen fourth moment/Sobolev power | no; identifying it with `d_ST` was circular |
| old `d_N=3.9951` window | selected finite-shell regression | rejected by the registered full-curve test |
| exact finite moments | Taylor coefficients of a finite heat trace | no |
| Kähler--Dirac heat flow | diffusion of the theory's own cochain operator | never reaches four |
| tower and warped products | tested finite dynamical extensions | no stable 4D plateau |

Several registered phenomenology scripts still substitute `d_ST=4` into
downstream formulae, sometimes renaming the arithmetic trace `4` as a
spacetime derivation.  Those substitutions are assumptions or patterns, not
independent probes.  Their numerical identities may remain true after the
rename, but every physical conclusion that needs a derived four-dimensional
Sobolev/Seeley exponent is now **OPEN** and requires a separate audit.

Accordingly `logic_chain_map.md` section 5 has been corrected: it no longer
lists `spectral index -> spacetime dimension` as a structural bridge.  The
bridge is rejected.  Section 16 now calls `(c0,c1,c2)` finite moments and
explicitly rejects the Seeley--DeWitt label on the fixed complex.

## Status ledger

### DERIVED

- registry uniqueness guard and the actual one-duplicate provenance;
- 600-cell boundary `f`-vector, ambient rank four and simplicial dimension
  three;
- `d^2=0`, ranks `(119,601,599)`, Betti numbers `(1,0,0,1)`, and two harmonic
  modes;
- complete `D^2` range, moment triple, and unique heat-flow maximum;
- reproduction of the old target-conditioned `d_N=3.068762` interval and
  absence of a 4D counting interval;
- exact non-intertwining of the `Box_p` hierarchy with the simplicial
  coboundary;
- finite heat-trace Taylor behavior and absence of a Seeley exponent;
- scoped negative for the tested path-product route.

### STRUCTURAL

- reading the boundary carrier as the theory's spatial geometry;
- using a finite heat-flow shoulder as an approximation to continuum spectral
  dimension.

### PATTERN

- convergence of vertex, Hasse and Kähler--Dirac shoulders toward three;
- the short `D^2` shoulder near `3.27` as a positive dimension measurement;
- downstream appearances of the algebraic integer four as spacetime
  quantities.

### OPEN

- a preregistered refinement sequence with a stable half-decade-or-longer
  Kähler--Dirac plateau and positive higher-dimensional controls;
- a derived Lorentzian/time direction and a scale-matched `3+1` operator;
- continuum Seeley--DeWitt coefficients from a proved refinement asymptotic;
- the consequences for every mass/coupling formula that currently uses
  `d_ST=4` as an input.

### REJECTED

- `Box` alternating nullity as spacetime dimension;
- finite moments as Seeley--DeWitt coefficients;
- the exact unregistered `3.287/0.285` claim as a reproducible result;
- any assertion that the fixed cochain geometry itself supplies four
  dimensions.

## Reproducibility

Primary verifier:

`/home/razvan/science/.venv/bin/python reproducible/verify_dimension_reconciliation.py`

The `Box_p`/coboundary non-intertwining checks are registered in
`reproducible/verify_edge_gauge_spectrum.py`.  The new verifier is registered
in `reproducible/run_all.py`.  No PDF build was attempted.

Final command, run from the repository root on commit `f7adc8c`:

`/home/razvan/science/.venv/bin/python reproducible/run_all.py`

Result: `79/79 scripts completed successfully` in `901.0 s`.  The registry
contained 79 distinct verifier names, the bidirectional coverage and
duplicate guards passed, and the process exited `0`.  No PDF build was
attempted.
