# Hostile audit of the Q8 factor-swap transplant to the `m=120` arena

Date: 2026-07-28

## Decision

The proposed matter carrier is **REFUTED**.  The narrower
configuration-dependent zero-fluctuation claim **SURVIVES**.

The right-algebra/left-symmetry configuration on

`H0 = C[2I]_x tensor C[2I]_y`

does evade the factor-preserving no-go.  With flip-star `J0`, the subgroup
average `X`, and `D0=X+J0XJ0^-1`, it has exact order zero, first order,
reality, equivariance, and nonzero commutators.  Thus the failure reported
for the left-algebra/adjoint-symmetry configuration was not an intrinsic
impossibility theorem for every configuration on the same vector-space
arena.  **DERIVED.**

However, the chirality-doubled construction is not a manifold-like finite
spectral triple.  Orientability fails, the intersection form is identically
zero, Poincare duality fails maximally, and connectedness fails.  The result
is therefore only algebraic KO-dimension-six data with nonzero one-forms.
It does **not** reopen the physical matter gate.  **DERIVED negative.**

## Verdicts

### 1. Degeneracy — WEAKENED

For a subgroup `C10`,

`X^2=10X`.

Thus `X=10P`, where `P` is an orthogonal rank-12 projection.  The two
factor projections commute, so `D0` has exact spectrum

| eigenvalue | multiplicity |
|---:|---:|
| `0` | `11664` |
| `10` | `2592` |
| `20` | `144` |

After chirality doubling, the kernel has dimension `23328`.

A kernel, even a large one, does not violate compact resolvent in finite
dimension.  The fatal manifold-like defects are instead:

- **DERIVED:** connectedness fails because every represented element of
  `C[C10]`, not only scalars, commutes with `D`;
- **DERIVED:** metric-dimension-zero orientability fails.  Every operator
  `sum pi(a) J pi(b) J^-1` is identical on the two added sheets, while
  `Gamma=diag(1,-1)` has opposite sheet signs;
- **DERIVED:** every entry
  `Tr(Gamma pi(p_i) J pi(p_j) J^-1)` cancels between sheets.  The `9 x 9`
  intersection matrix is zero and has rank zero, so Poincare duality fails.

Concrete evasion boundary: the listed bilinear/reality axioms survive, but
the construction is not an orientable, Poincare-dual, connected finite
geometry.

### 2. Chirality doubling — WEAKENED

The block formulas do give the KO-6 signs exactly.  **DERIVED.**

The claim that such doubling is familiar from NCG particle models is only
an analogy.  It does not derive these particular sheets, their identical
algebra representation, or the off-diagonal `D`.  **PATTERN.**

No identification with either repository double is available:

- the Galois double is a 60-dimensional weighted McKay-node construction
  with a Galois permutation and a scoped node algebra;
- the primal-dual double is a 5280-dimensional cellular construction with
  form parity and a cellular Hodge map;
- the present doubled space has dimension 28800 and no supplied
  intertwiner transporting either algebra, grading, `J`, and `D`.

Dimension mismatch alone does not prove that no embedding can ever exist,
but no canonical carrier map is present.  Moreover, the currently available
`C^9` representation on the Galois double also has an identically zero
intersection form: its Galois permutation preserves both node dimensions
and McKay parity, so paired graded traces cancel.  The SM algebra action is
not defined there, hence no other intersection form can honestly be tested.
The present doubling therefore remains **STRUCTURAL**, and a viable
nontrivial identification is **OPEN**.

Concrete evasion boundary: it is legitimate as an explicit algebraic
construction, but illegitimate as a claimed derived grading or as a repair
of the canonical-arena matter program.

### 3. Reading of the canonical arena — SURVIVED

The verifier `verify_canonical_bimodule_arena.py` makes separate choices:

- lines 88–98 construct flip-star and the canonical left/right order-zero
  placement;
- lines 100–146 introduce diagonal adjoint symmetry and decompose it;
- lines 148–176 test the central adjacency/class-sum candidates;
- lines 178–194 exclude only the listed derived grading candidates;
- lines 196–202 check the nonunital rank-14 corners.

No line proves that diagonal adjoint symmetry is forced.  No prior theorem
cited by the note forces it either.  The note itself says “used here” and
limits its conclusion to its derived candidate list.  **DERIVED audit.**

The new configuration changes both entries relevant to the old failure:
the symmetry is left translation on `x`, and the represented algebra is
right multiplication on `x`.  Left/right regular actions commute, while a
nonnormal subgroup average commutes with the former but not the latter.
That is the exact missed horn.  **DERIVED.**

Concrete evasion boundary: this refutes an intrinsic no-go across all
configurations on the `14400`-dimensional undoubled vector space.  It does
not refute the negative result for the canonical note's fixed
left-algebra/adjoint-symmetry configuration.

### 4. Physics — SURVIVED as a negative

For

`C[2I] = direct_sum M_d(C)`,
`d=(1,2,2,3,3,4,4,5,6)`,

the unitary group of the full complex algebra is

`U(1) x U(2)^2 x U(3)^2 x U(4)^2 x U(5) x U(6)`.

Calling this directly “the gauge group” is slightly too strong: the
spectral-triple gauge group is an image/quotient of unitaries and may also
be subject to unimodularity.  The raw unitary group is nevertheless nothing
like the Standard Model group.  **STRUCTURAL.**

The proposed `C+H+M3(C)` support has complex regular rank
`1^2+2^2+3^2=14`.  Its corner unit is nonunital on the full arena.
**DERIVED.**

No Standard-Model representation, hypercharge, anomaly condition,
generations, or Yukawa sector follows.  The algebraic gate opens while the
physics gate stays closed.  **DERIVED negative.**

### 5. Exactness and tensor reduction — SURVIVED

The group is constructed over `Q(phi)` and the operators used in the axiom
checks are integer sparse matrices.  No floating-point comparison or
eigenvalue computation is used.  **DERIVED.**

The reduction is airtight:

`pi(a)=R_a tensor I`,

`J pi(b) J^-1=I tensor L_(b^-1)`,

and

`[D0,pi(a)]=[X,R_a] tensor I`.

Therefore order zero and first order are identities of the form
`[A tensor I,I tensor B]=0`.  Linearity extends them from group elements to
all of `C[2I]`.  The step-grid loops are only redundant controls and are not
the proof.  **DERIVED.**

Concrete evasion boundary: the reduction proves the bilinear axioms but
does not prove orientability, Poincare duality, connectedness, or physical
adequacy; those were separately audited and fail.

### 6. `C10`, centrality, and canonicity — WEAKENED

Exact enumeration gives six order-ten subgroups in `2I`, all conjugate.
Each is inverse-closed and nonnormal.  For a subgroup `H`,

`sum_(h in H) h`

is central exactly when `H` is invariant under conjugation.  Hence every
one of these six subgroup sums is genuinely noncentral.  The verifier also
finds exactly 100 of 120 group elements with nonzero commutator.  **DERIVED.**

The individual choice is not canonical: `next(g for ... order(g)==10)`
depends on the coordinate ordering.  The geometry supplies a conjugacy
orbit of Hopf-fiber partitions, not a distinguished member in the current
data.  Moreover, the same mechanism works for any inverse-closed nonnormal
subgroup sum.  The selected `X` is therefore a **STRUCTURAL carrier choice**
from a derived family, and its claimed special physical role is only a
**PATTERN**.  Uniqueness or an independent selector is **OPEN**.

## Status ledger

### DERIVED

- the right-algebra/left-symmetry factor-swap data satisfy the stated
  algebraic KO-6, order-zero, first-order, and nonzero-form identities;
- the tensor-factor reduction is exhaustive;
- `X=10P`, `rank(P)=12`, and the exact spectrum/multiplicities above;
- exactly six conjugate `C10` subgroups exist and their sums are noncentral;
- orientability, connectedness, and Poincare duality fail;
- the actual `128 x 128` Q8 positive control has the same trivial-doubling
  failure, with all 64 group-element graded traces equal to zero;
- the scoped `C^9` Galois-double intersection form is also rank zero;
- the literal SM corner is rank 14 and nonunital.

### STRUCTURAL

- the added chirality double;
- selection of one `C10` carrier from its conjugacy orbit;
- the full complex group algebra as a physical internal algebra.

### PATTERN

- analogy with particle/antiparticle-chirality doubling in Standard-Model
  NCG;
- interpreting the selected subgroup average as a distinguished matter
  carrier.

### OPEN

- a canonical identification with the Galois or primal-dual double;
- an independent geometric selector for one `C10` carrier;
- any orientable, Poincare-dual repair that does not enlarge the arena or
  fit `D` or `gamma`.

## Final boundary

The honest conclusion is:

> The earlier `m=120` failure was configuration-dependent at the level of
> the tested algebraic axioms.  A missed configuration has nonzero inner
> fluctuations.  But its manufactured grading makes the intersection form
> identically zero, so it is not a manifold-like finite spectral triple and
> it supplies no Standard-Model physics.
