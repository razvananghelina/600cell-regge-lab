# The missing-link audit

Date: 2026-08-10

## Decision

There is no evidence for one missing numerical identity between the exact
600-cell geometry and a physical theory.  The current construction stops at
several logically independent selections.  The shortest honest name for the
missing object is a **refinement-natural dynamical spectral functor**: a rule
which assigns the algebra, chiral module, Dirac data and state to every
refinement and intertwines them without inspecting a Standard-Model target.

This is a specification of the missing object, not its construction.
Merely naming the functor does not advance the physics gate.  **OPEN.**

The theory's existing finite spectral action is the strongest internal
candidate for the selection rule, but it has no projective scale in its
domain and therefore cannot vary or select one.  The exact calculation below
tests its most direct local extension to the projective cometric.  That
extension does not select even the first new Dirac scale unless its cutoff
function or polynomial coefficients already contain the choice.  Together
with the independent finite-Dirac critical-circle result, this rules out the
spectral-action prescription **as presently defined and directly extended**
as the missing link.  It does not rule out a future coupled action with new
dynamical terms.  **DERIVED NEGATIVE under the complete hypotheses below.**

Verifier: `reproducible/verify_missing_link_selection.py`.

## 1. What the established geometry supplies

The following inventory is already certified and is not re-derived here.

1. The boundary of the 600-cell is a simplicial `S^3` with exact cochain
   complex, Betti numbers `(1,0,0,1)` and Kähler--Dirac `d+d*`.
2. The unit round embedding selects a Levi--Civita `SU(2)` spin transport
   once the round metric and short geodesics are accepted.  Its face
   holonomy is nontrivial and its spin-fibre algebra is `M2(C)`.
3. Normalised Haar measure is projectively consistent under edge
   subdivision.
4. The projective configuration cometric is not selected: already on one
   regular tetrahedron and its barycentric subdivision there is an exact
   family

   `K_f(t) = H K_c H^T + t Q`, `t>0`, `rank(Q)=44`,

   with the same coarse pullback and all tetrahedral symmetries.
5. The natural diagonal holonomy action on `S+S` remains `M2(C)` and does not
   separate the copy label.  A larger finite factor requires a representation
   choice.
6. The canonical `(2,3,5)` line-incidence census does not produce the derived
   `M16` matter target.
7. The static geometry contains no fourth direction.  Time, Lorentzian
   signature and a controlled refinement limit remain absent.

Items 4--7 can be changed while items 1--3 are kept fixed.  Therefore the
present geometric data do not logically determine those outputs.  This is
stronger than a failure to find the right formula: explicit alternative
extensions already exist.  **DERIVED/STRUCTURAL independence statement,
with each component scoped by its cited audit.**

## 2. Test of the direct spectral-action extension

### Complete hypotheses

The no-selection statement assumes all and only the following.

1. The local refinement is the complete barycentric subdivision of one
   regular tetrahedron, with 6 coarse and 50 fine edges.
2. `A` is the linearised composition of the two oriented fine half-edges on
   every coarse edge.
3. `H=A^T(AA^T)^-1` and `Q=I-HA`.
4. The fine configuration cometric is
   `K_f(t)=K_h+tQ`, where `K_h=H K_c H^T`, `K_c` is positive and `t>0`.
5. A proposed spectral functional sees this sector through
   `S_f(t)=Tr f(D_t^2)` with `D_t^2=K_f(t)` locally.

Hypothesis 5 is the most favourable possible reading for spectral selection:
it assumes that a global Clifford module and Dirac operator realising this
square have already been supplied.  Those objects are themselves still
open.

### Exact reduction

The identities

`Q^2=Q`, `rank(Q)=44`, `K_h Q=Q K_h=0`

give, for every polynomial `f`,

`Tr f(K_f(t)) = Tr_horizontal f(K_h) + 44 f(t)`.

The same block decomposition gives the corresponding functional-calculus
identity for continuous spectral cutoffs on this finite block.  Hence a
stationary positive scale exists exactly when the externally specified
scalar function `f` has a stationary point there.

Three relevant cases are exact:

- for `f(x)=a x+b x^2`, with `a>=0`, `b>0`, the action is strictly increasing
  for `t>0`; its infimum is the degenerate boundary `t=0`;
- for the heat cutoff `f(x)=exp(-s x)`, `s>0`, it is strictly decreasing and
  has no finite critical point;
- for `f(x)=b x^2-a x`, `a,b>0`, criticality gives
  `t=a/(2b)`.  The scale is exactly the unselected coefficient ratio.

Thus the spectral action either does not select a finite `t`, or merely moves
the free parameter from the cometric into the cutoff.  **DERIVED NEGATIVE.**

The conclusion is unchanged for the natural constant-coefficient torus
completion.  A Fourier mode `k` has

`lambda_k(t)=k^T K_h k+t k^T Q k`,

and `k^T Q k=||Qk||^2>=0`.  Every non-increasing positive cutoff therefore
makes every mode contribution non-increasing in `t`; the heat action has no
finite interior minimum.  This extends the monotonic part of the result, not
the special finite formula `44 f(t)`, to the local Fourier tower.

This local result is decisive only against a direct sum of the present local
spectral terms.  Couplings between tetrahedra could lift the local mode, but
specifying those couplings would be new dynamical data and must itself pass
the same target-independent selection test.  No no-go against every possible
nonlocal dynamics is claimed.

### Independent finite-sector control

The failure is not peculiar to the projective connection space.  On the
existing `C3` finite-Dirac arena, the spectral-triple axioms leave a generic
122-real-dimensional gauge-and-scale quotient.  Positive quadratic/quartic
functionals select only zero, while the first symmetry-breaking quartic has
an exact gauge-inequivalent critical `S1`.  Therefore polynomial spectral
criticality has already failed to select the internal Dirac independently of
the new parameter `t`.  **DERIVED NEGATIVE**, certified by
`reproducible/verify_dirac_selection.py`.

## 3. Why a scalar action cannot yet select all missing objects

`Tr f(D^2)` is a functional only after its Hilbert space and operator have
been declared.  In the current repository it does not range over:

- real associative algebra types and their faithful bimodules;
- chiral matter functors or higher-rank equivariant bundles;
- the connection group `SU(2)` versus `U(2)`;
- Lorentzian time evolutions, refinement embeddings or continuum states.

Consequently it cannot select these objects by variation.  Adding them to a
larger configuration space is possible, but defining that space, its measure
and its comparison of inequivalent Hilbert spaces is precisely the missing
construction.  **STRUCTURAL logical boundary.**

The current evidence therefore refutes the framing “one coefficient is
missing”.  It does not prove that no single deeper principle could determine
all sectors.  It proves only that no such principle is presently defined and
that the existing spectral action is not it.

## 4. Minimal specification of the missing link

Let `Ref(X)` denote the refinement category generated from the 600-cell
spatial carrier.  A candidate link must assign, at every level `n`, data

`F(X_n)=(G_n,A_n,H_n,D_n,J_n,gamma_n,Omega_n)`

and refinement morphisms which satisfy all of the following before any
particle target is inspected:

1. **Functoriality:** composition of refinements gives composition of the
   Hilbert/state embeddings and intertwines inherited observables.
2. **Kinematic selection:** `G_n`, the vertical cometric scales and the
   Clifford module are unique up to declared gauge/unitary equivalence.
3. **Finite selection:** `A_n` and its faithful chiral bimodule are selected,
   not merely embeddable, and satisfy the real spectral-triple gates.
4. **Matter content:** kernels or cohomology of selected operators supply the
   matter module; Euler characteristics or virtual identities do not count.
5. **Dynamics:** the state/action selects the remaining Dirac data rather
   than importing a cutoff function with equivalent free coefficients.
6. **Spacetime:** a Lorentzian evolution or reconstruction supplies the
   missing time direction and has a controlled continuum/EFT limit.
7. **Rigidity and provenance:** the admissible family and look-elsewhere
   count are fixed before comparison with Standard-Model representations or
   measured constants.

No current construction satisfies this list.  It is nevertheless a useful
acceptance boundary: any proposed “missing link” can now be killed at the
first item it fails, rather than being credited for compatibility with one
desired output.

## 5. Current logical cut

The repository's current status map already displays three decisive gates:

- selected algebra/representation;
- canonical matter functor;
- time/refining spacetime dynamics.

The HD audit adds the projective Dirac scale beneath the third gate.  The
finite spectral action does not cross either the projective-Dirac gate or the
finite-Dirac gate.  Hence the shortest current verdict is:

> The missing link is not yet a theorem or an equation.  It is an absent
> target-independent functorial dynamics that must select the spectral data
> themselves across refinement.

This diagnosis narrows the next search, but it is not a physical theory.
**OPEN.**

## 6. Can the missing object be only a state?

A distinguished state is an attractive smaller candidate.  It could provide
a GNS representation, a covariance scale and, on a noncommutative von
Neumann algebra, a modular automorphism flow.  The thermal-time proposal is
the relevant precedent: it reads a state-dependent time flow from
Tomita--Takesaki modular theory ([Connes--Rovelli,
1994](https://arxiv.org/abs/gr-qc/9406019)).  This precedent motivates the
test; it does not identify the repository's state.

The finite audit gives three exact boundaries.

1. The normalized trace on `M2(C)` has density `rho=I/2`; its modular
   generator `[log(rho),a]` vanishes for every `a`.
2. The configuration multiplication algebra is commutative.  Every state on
   a commutative algebra is tracial, so its modular automorphism group is
   trivial.  The projectively consistent Haar state therefore supplies no
   time flow on this algebra.
3. A faithful nontracial state
   `rho_p=diag(p,1-p)` on `M2(C)` gives

   `sigma_s(E12)=exp(i s log(p/(1-p))) E12`.

   The modular frequency is the arbitrary density-eigenvalue ratio.  For
   example `p=2/3` gives `log(2)` and `p=3/4` gives `log(3)`; any positive
   ratio can be inserted by choosing `p`.

There is a second obstruction.  The faithful GNS representation of `M2(C)`
on its four-dimensional Hilbert--Schmidt space is left multiplication
`A -> A tensor I2`.  Its image and commutant both have complex dimension
four.  It is unitarily the same multiplicity structure which already failed
pure-state separation; a state on `M2` cannot manufacture the 16-dimensional
algebra `M4` or select a different algebra type.

Adding finite geometric symmetries does not change the canonical-state
conclusion.  For a measure-preserving finite action, the coefficient-of-the-
identity state on the crossed product is tracial.  The verifier exhausts the
control `C(Z3) crossed_product Z3 = M3(C)`: all 81 ordered basis products pass
traciality, so even this full matrix algebra has trivial modular flow in its
canonical invariant state.  Noncommutativity alone is insufficient; a
nontracial state must also be selected.

**DERIVED NEGATIVE:** neither the canonical Haar/trace state nor an arbitrary
faithful state is the missing link.  The former has no modular time; the
latter imports its clock rate through an unselected density matrix and cannot
change the algebraic type.  A modular-time route would first require a
geometrically selected noncommutative algebra and a selected nontracial state
on it.  Those are components of the functor in Section 4, not a replacement
for it.

Exact verifier: `reproducible/verify_missing_link_modular_state.py`.

## 7. What the recently revised convergence theorem contributes

The March 2026 revision of Farsi--Latremoliere--Packer proves convergence of
inductive sequences of metric spectral triples in spectral propinquity under
a full-quantum-isometry bridge condition.  Such convergence controls state
spaces, spectra, bounded functional calculus and the quantum dynamics
generated by the already supplied Dirac operators
([arXiv:2301.00274](https://arxiv.org/abs/2301.00274)).

This is directly useful as an **acceptance topology** for Section 4: it can
make “the refinement converges” a mathematical statement rather than a
visual plateau.  It is not a selection principle.  Its hypotheses begin with
an inductive sequence of algebras, Hilbert spaces and Dirac operators; its
examples add proper length functions with bounded doubling.  It tests data
of exactly the type which are still unselected here.

The 600-cell refinement quantifies the remaining freedom.  For the full
barycentric subdivision, not the 7440-edge Hasse control graph, the exact
f-vector map is

```text
[V']   [1 1 1  1] [V]
[E'] = [0 2 6 14] [E]
[F']   [0 0 6 36] [F]
[T']   [0 0 0 24] [T].
```

Starting from `(120,720,1200,600)`, the first level is
`(2640,17040,28800,14400)`.  The eigenvalues of the subdivision map are
`1,2,6,24`, and the edge count has the exact form

`E_n=(120/11) 2^n+(7800/11) 24^n`.

Consider the minimal level-weight model in which a new connection mode born
at level `n` receives Dirac eigenvalue `a_n=b^n`, `b>1`.  This is a
**STRUCTURAL model**, not the still-missing global Dirac.  Its zeta series is
an exact pair of geometric series and has convergence abscissa

`d(b)=log(24)/log(b)`.

Thus compact resolvent and exponential refinement do not determine a
dimension:

- `b=2` gives `d=4.5849625...`;
- `b=3` gives `d=2.8927893...`;
- `b=24^(1/3)` gives `d=3`;
- `b=24^(1/4)` gives `d=4`.

The 3-dimensional volume split supplies a structural reason for the third
choice: every top simplex has volume `1/24` of its parent, so inverse isotropic
linear scaling is `24^(1/3)`.  This recovers three only after the spatial
three-dimensional scaling law is assumed; it is a consistency check, not an
independent dimension derivation.  Choosing `24^(1/4)` to obtain four would
be target fitting and conflicts with that spatial scaling interpretation.

**DERIVED CONDITIONAL NEGATIVE:** convergence/compact-resolvent criteria do
not select the Dirac scaling law.  Spectral propinquity gives the right tool
for validating a proposed functor, but a separate dynamical principle must
still construct it and fix `b`.  Exact verifier:
`reproducible/verify_missing_link_refinement_scaling.py`.

## 8. Attack on the singular “missing link” framing

The functor in Section 4 is an acceptance envelope, but calling it *the* link
risks hiding the problem by putting every missing object into one tuple.  No
current theorem says the same principle must select external spacetime,
internal matter and their physical matching.  The dependency graph instead
has three contracts which can fail independently:

1. **External/refinement contract:** construct a Lorentzian or modular
   inductive spectral geometry, including a derived scale law `b`, a
   nontracial state if modular time is used, and convergence in a topology
   such as spectral propinquity.
2. **Finite/matter contract:** select a real algebra, faithful chiral bimodule,
   finite Dirac and matter kernel/cohomology.  The line-incidence route is
   closed; higher-rank bundles remain open only with a preregistered
   canonicity principle.
3. **Product/matching contract:** prove how the external and finite triples
   combine, how charges and anomalies are fixed, and how the resulting
   effective dynamics is matched and renormalized.

The fixed spatial geometry can be held constant while the refinement weight
`b` changes.  The same fixed chamber carrier admits thousands of finite
algebra designs and at least one exact noncommutative all-gate witness.  These
are explicit independence witnesses, not merely missing calculations.

Therefore the singular framing is **REFUTED as a claim about the current
logical graph**.  A future deeper theory could unify the three contracts, but
the repository contains no evidence that it does.  The most advanced route
is presently Contract 1, because its carrier, connection maps, Haar measure
and a rigorous convergence topology are available.  Its sharp missing datum
is now:

> a geometrically selected nontracial dynamics which fixes the refinement
> Dirac scale law and supplies Lorentzian time.

Even success there would leave Contract 2 open.  This is a smaller and more
falsifiable target than “derive the Standard Model from the 600-cell.”

## Status ledger

- **DERIVED:** the projective cometric has a rank-44 free vertical sector.
- **DERIVED:** its spectral functional reduces to a fixed horizontal term
  plus `44 f(t)`.
- **DERIVED NEGATIVE:** positive moments and heat cutoffs have no finite
  positive stationary scale.
- **STRUCTURAL:** a symmetry-breaking polynomial can choose a scale only
  through an unselected coefficient ratio.
- **DERIVED NEGATIVE:** polynomial spectral criticality independently fails
  to select the finite Dirac.
- **DERIVED NEGATIVE:** the Haar/trace state has trivial modular flow; a
  nontracial state moves the time scale into an arbitrary density ratio.
- **DERIVED NEGATIVE:** adjoining finite measure-preserving symmetries gives a
  noncommutative crossed product but leaves the canonical state tracial.
- **DERIVED NEGATIVE:** faithful GNS of `M2` remains an `M2` image with a
  four-dimensional commutant and cannot select `M4`.
- **DERIVED:** full barycentric edge growth is controlled by the exact factor
  `24`.
- **DERIVED CONDITIONAL:** level weights `b^n` have spectral abscissa
  `log(24)/log(b)`; convergence alone leaves `b` free.
- **STRUCTURAL:** inverse spatial linear scaling gives three, not four.
- **STRUCTURAL ADVANCE:** spectral propinquity supplies a rigorous convergence
  test after a sequence of spectral triples has been constructed.
- **STRUCTURAL:** the current diagram has multiple independent selection
  gates; this is a dependency audit, not a no-go theorem against every future
  unifying principle.
- **OPEN:** construct the refinement-natural dynamical spectral functor.
- **OPEN, SHARP EXTERNAL TARGET:** select a nontracial refinement dynamics and
  its Dirac scale law before using spectral propinquity to test convergence.
- **NOT CLAIMED:** Standard-Model matter, four-dimensional spacetime, a
  unique continuum Dirac operator or a parameter-free spectral action.
