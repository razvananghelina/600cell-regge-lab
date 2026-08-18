# Adversarial audit of the 2026-07-22/23/24 theory fixes

Date: 2026-07-24

## Executive result

Baseline: `PYTHONPATH=/tmp/science-python-deps python3 reproducible/run_all.py`
completed in 243.8 seconds with `46/46` scripts successful.

This audit found one broken provenance claim and four materially overstated
claim/certificate pairs.  It did not find a counterexample to the
central-parity theorem, the subgroup classification, the four-parameter
bracket classification, either dimension-group calculation, the scoped
matter no-gos, the JUNO arithmetic, or the Kähler--Dirac module
decomposition.

Summary by finding:

- **BROKEN: 1.** The holographic/tower/warped criteria cannot be called
  verifiably preregistered from this repository.
- **WEAKENED: 4.** C3 “Dirac existence”; canonicity of the sampled A5 color
  metric; exactness of the invariant spectral gap; exactness of the
  mass-ratio negative.
- **SURVIVED: 8 principal cores.** Each target retained at least one
  substantive theorem after attack; details follow.

Severity ranking:

1. **High — BROKEN provenance:** no immutable record establishes that the
   numerical plateau thresholds, controls, tower weights, and warp choices
   preceded inspection of results.
2. **High — WEAKENED C3 claim:** the segregation verifier had no real algebra
   representation, opposite action, `J`, matrix `D`, or first-order double
   commutator.  It certified only a shared-index Krajewski block type.
3. **Medium — WEAKENED exact spectrum:** `gap=phi^-4` and the ratio no-go were
   inferred from double-precision clusters.  Exact moments survive through a
   new independent integer-matrix check.
4. **Medium — WEAKENED canonical metric:** the A5 verifier constructs an
   abstract icosahedral sampling embedding but does not compare it with the
   actual 720-edge kernel intertwiner.
5. **Low — verifier fidelity defect fixed:** the advertised
   `188,908,396` and `1,362,811,872,984` Krajewski counts were printed but
   previously asserted only as `>0` and increasing.  Exact equality is now
   asserted; both totals pass.

## Target-by-target attack log

### 1. Central-parity segregation

**ATTACK ->** Separate the theorem from the claimed C3 first-order Dirac
existence.

**METHOD ->** Read the note and verifier line by line; reran the verifier;
checked `SL(2,5)` directly; inspected the C3 cell indices and compared the
verifier's `row equal or column equal` Boolean with the matrix first-order
condition.  Added an independent residual-character intersection check for
the two proposed cells.

**OUTCOME -> WEAKENED.** The central lemma, unique involution, odd subgroup
types `C1,C3,C5`, subgroup counts, and odd-endomorphism dimensions survived.
The C3 witness has the advertised parity assignment, a shared left
Krajewski index, and common C3 character support, so a nonzero residual C3
intertwiner exists.  But no finite real spectral triple or matrix-level
first-order equation was constructed.  The note and verifier now call this a
residual-equivariant Krajewski-legal block candidate, not a Dirac theorem.

Concrete evasion boundary: a different residual-C3 algebra/opposite-action
allocation can make the double commutator differ from the shared-index
combinatorics.  The old wording silently treated that allocation as already
defined.

### 2. A5-equivariant brackets

**ATTACK ->** Challenge completeness of the four-parameter family and the
metric exclusion of noncompact branches.

**METHOD ->** Recomputed the character dimension encoded in the verifier;
checked that the four tensor maps are linearly independent; inspected the
full quadratic Jacobi coefficient row space and the real branch reduction;
tested the metric-invariance constraint against the constructed sampling
Gram matrix.

**OUTCOME -> WEAKENED.** The four-dimensional equivariant-map count, Jacobi
variety, real branch classification, compact/split Killing signatures, and
rigidity under block rescaling survived.  Positive-definite metric
invariance does exclude the split, semidirect, and nilpotent representatives
for the displayed metric.  What did not survive is the unqualified word
“canonical”: the verifier builds the 12-point sampling map but never
reconstructs the actual 720-edge gauge kernel and identifies its isotypic
intertwiner with that map.  The note now makes the metric conclusion
conditional on that embedding.

Concrete evasion: an independently normalized A5 intertwiner from the
abstract `3'+5` to the edge kernel changes the two Schur scalars.  It does not
rescue a noncompact algebra from *some* positive invariant metric, but it can
change the claimed geometric ratio and selected representative.

### 3. Bratteli inflation

**ATTACK ->** Treat finite floors as insufficient evidence for an AF limit
and challenge the golden mass-lattice embedding.

**METHOD ->** Inspected the unimodular Fibonacci bonding map, the telescoped
rooted McKay stable image, Smith forms, trace row, primitive order, and parity
argument.  Reran the exact verifier.

**OUTCOME -> SURVIVED.** The Fibonacci limit does not rest on extrapolating
12 floors: `det(F)=-1` makes every bonding map a lattice automorphism, and the
PF functional identifies the ordered group with `Z[phi]`.  The mass lattice
identification is canonical only as the same abstract rooted fusion
polynomial/state; the note already labels physical mass/scale interpretation
structural.  The McKay result is honestly the stationary presentation
`lim(Z^4,M)` with dyadic trace and infinitesimals.  The all-level matter
no-go follows from bipartite support, not tested depth.

Concrete evasion already admitted by the note: a different AF
representation/opposite action or cumulative-floor matter module can evade
the fixed-floor/canonical-endpoint no-go.  It is not covered by the theorem.

### 4. Invariant spectrum

**ATTACK ->** Demand exact arithmetic for the exact gap/moments and test
whether the mass hunt is genuinely exhaustive.

**METHOD ->** Read the clustering and joint diagonalization path; tightened
the moment check through a different code path using integer Laplacian
entries and `Tr(Delta^2)=sum_ij Delta_ij^2`; inspected all 1326 numerical
ratios and the `1e-10` log tolerance.

**OUTCOME -> WEAKENED.** The exact integer checks now independently give
`Tr(D^2)=14880` and `Tr(D^4)/2=55920`; `c0=2640` is a dimension identity.
The gap and unresolved roots are still double-precision eigenvalues.
Therefore `gap=phi^-4` is numerically identified, not exactly certified.
The mass hunt is exhaustive over the 52 clustered numerical levels at its
stated tolerance, but is not an exact-algebraic negative for unresolved
roots.  The note has been corrected accordingly.

### 5. Holographic dimension, tower spacetime, warped spacetime

**ATTACK ->** Audit whether criteria were genuinely frozen and whether
controls validate the dimension verdict.

**METHOD ->** Compared the March protocol with the July verifier, checked
filesystem timestamps, searched for repository history, reran the baseline
outputs, and inspected the decision gates and controls.

**OUTCOME -> BROKEN provenance; SURVIVED scoped numerics.** There is no Git
history or external timestamped registry.  The older protocol specifies no
numerical thresholds.  A verifier constant appearing textually before its
spectral function prevents mutation during one run; it does not prove
pre-result registration.  All three notes now say “fixed for rerun” and mark
pre-result provenance unverified.

The negative results themselves survive conditionally: under the displayed
criteria the 600-cell has no calibrated 4D plateau, and the derived golden,
McKay, and warped finite constructions have no N-stable 4D plateau.  The
T4 control also fails to demonstrate positive sensitivity, so the correct
interpretation remains a scoped negative/inconclusive calibration—not a
measurement of dimension.

Concrete evasion: other tower coefficients, boundary conditions, larger
truncations, or plateau criteria are outside these no-gos.  The notes largely
say this already.

### 6. Matter no-go chain

**ATTACK ->** Challenge the `1.9e8` enumeration, representation-support
constraint, and each `D=0` scope.

**METHOD ->** Inspected the generating-function DP, all 25/50 cell fusion
vectors, missing irrep support, endpoint Krajewski rule, preprojective
diagonal bimodule, and functor-route availability statements.  Replaced weak
count assertions with exact equality assertions.

**OUTCOME -> SURVIVED after verifier repair.** Both advertised enumeration
totals are reproduced exactly.  More importantly, the no-go does not require
enumerating all dimension-only matrices: the derived diagonal restriction
omits required irreps.  The edge and preprojective `D=0` results are valid for
the full canonical endpoint/diagonal algebras.  The functor note correctly
calls missing representations undefined rather than inconsistent.

Concrete evasions are real but already scoped out: a separately derived
smaller algebra, non-diagonal action, different opposite representation, or
different matter Hilbert space could allow nonzero odd blocks.  None is
ruled out by this chain.

### 7. JUNO comparison

**ATTACK ->** Recheck external constants and look for undisclosed variant
selection.

**METHOD ->** Recomputed the exact SymPy observables and asymmetric pulls;
opened the cited primary arXiv records online.  JUNO's
`0.3092+-0.0087` and `(7.50+-0.12)e-5`, DESI's `0.0642 eV`, the
dynamical-DE `0.098(+0.016,-0.037) eV`, and KATRIN's `0.45 eV` are confirmed.
The NuFIT arXiv record was available, but its official parameter-table PDF
returned an upstream error during this audit.

**OUTCOME -> SURVIVED with source limitation.** The arithmetic, directional
error choices, and two-input chi-square ranking survive.  Variant I is not a
derived selection: the note explicitly calls both the `1/45` correction and
its scope a pattern.  That disclosure prevents the ranking from becoming a
hidden model-selection claim.  The note now records that the NuFIT table
numbers were not independently re-fetched on 2026-07-24.

Garden-of-forking-paths warning: the three variants do not exhaust possible
state-specific corrections, correlations are ignored, and the preferred
variant was compared after the correction rule existed.  Therefore
`chi^2=0.609` is descriptive, not discovery significance; the note already
states this.

### 8. Kähler--Dirac matter

**ATTACK ->** Test `22 Reg`, grading independence, and whether `2640=c0`
changes Hilbert spaces.

**METHOD ->** Inspected all 120 signed cell actions, exact coboundary
equivariance, layer characters, spin/form parity matrices, and the definition
of `c0` used by the spectral-action verifier.  Reran the standalone verifier.

**OUTCOME -> SURVIVED.** The layer characters give
`Reg,6Reg,10Reg,5Reg`, hence `22Reg`, with dimension closure.  Spin parity is
an antipodal signed permutation; form parity is degree sign.  They commute,
are not equal up to sign, and `D` is spin-even/form-odd.  `2640=c0` is the
definition `Tr(I)` on the same oriented-cell cochain space, not a transported
identity from another Hilbert space.  The note is already careful that
matter-carrier availability does not select a finite algebra.

## Corrections made

- `segregation_theorem.md`: downgraded C3 Dirac existence to a
  residual-equivariant Krajewski candidate and stated the missing
  matrix-level data.
- `verify_segregation_theorem.py`: renamed the test and added a separate C3
  character-intersection certificate.
- `a5_equivariant_bracket_theorem.md`: made the edge metric conclusion
  conditional on the displayed sampling embedding.
- `invariant_spectrum.md`: downgraded exact gap and exact ratio-negative
  language to numerical certification.
- `verify_invariant_spectrum.py`: added independent exact integer moment
  checks.
- `holographic_dimension.md`, `tower_spacetime.md`,
  `warped_spacetime.md`: corrected unverifiable preregistration language.
- `verify_bimodule_krajewski.py`: now asserts both advertised enumeration
  totals exactly.
- `juno_2026_comparison.md`: recorded the 2026-07-24 source recheck and the
  NuFIT-table retrieval limitation.

No PDF build was attempted.

## Final suite evidence

After all corrections:

`PYTHONPATH=/tmp/science-python-deps python3 reproducible/run_all.py`

completed in 236.9 seconds with
`Result: 46/46 scripts completed successfully.`
