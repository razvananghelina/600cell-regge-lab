# Prior-art gate: two-slab 600-cell dust composition

Date: 2026-08-16

Status: **written before evaluating a new two-slab action, shared-slice
derivative or orbitwise momentum match**.

Upstream one-slab result: `94675be`

Upstream one-slab audit: `fd81cb5`

This is a targeted primary-source map, not a novelty proof.  External novelty
remains **OPEN** pending a dedicated review.

## 0. Framing correction before protocol

The first draft proposed gluing two copies of the published time-symmetric
sandwich and expecting their shared momenta to match.  That expectation is
**RETRACTED before preregistration and before any new two-slab evaluation**.

The already certified one-slab boundary data show

```text
sorted(p_pre) = sorted(-p_post)
```

to relative error `2.043e-15` in both parities.  This is exactly the expected
time-reversal behavior: the sandwich enters and leaves equal spatial
geometries with opposite momenta.  Repeating it would restart the motion and
create a momentum cusp.  A time-reversed copy also reverses the momentum and
does not provide the forward continuation.

Therefore the symmetric double is a **negative sign control**, not a known
stationary two-slab solution.  The corrected first protocol must test action
gluing and the formula for the shared derivative, while expecting the
repeated-sandwich shared residual to be nonzero.  Actual forward evolution
requires solving a new slab whose pre-momentum equals the first slab's
post-momentum.

## 1. Exact proposed object, operator, carrier and hypotheses

For an ordered five-phase schedule `sigma`, write the already certified
one-slab complete action as

```text
S_sigma(q_old, x, q_new),
```

with

```text
q_old in R_+^30   old 600-cell boundary edge-orbit squares,
x     in R_+^35   30 staircase diagonals + 5 pole squares,
q_new in R_+^30   new 600-cell boundary edge-orbit squares.
```

The action is the Lorentzian Regge curvature action with its boundary hinge
terms plus the published De Felice--Fabri dust world-line action on the five
pole orbits.  The mass and the symmetric control lengths are inputs from the
published solution, not fitted targets.

The first new object is the composed action

```text
S_2 = S_sigma(q0, x1, q1) + S_sigma(q1, x2, q2).
```

The shared 600-cell slice `q1` is internal to the composite.  With the
already certified conventions

```text
p_pre  = - partial S / partial q_old,
p_post = + partial S / partial q_new,
```

its equation is

```text
partial S_2 / partial q1 = 0
    iff
p_post(first slab) = p_pre(second slab).
```

This is a 30-component orbit equation.  Equality of sorted momentum
multisets is insufficient.  The calculation must derive the explicit
geometric permutation identifying every final-boundary edge orbit of the
first slab with the corresponding old-boundary orbit of the reversed second
slab, then test all 30 signed components.

The structural control carrier will use only the two canonical schedule
parity representatives, repeated on consecutive slabs, together with their
derived old/new layer-identification and time-reversal audit maps:

```text
even,
odd.
```

Arbitrary pairings or permutation searches are not candidates and must not
enlarge the look-elsewhere count.

At the repeated-sandwich sign control:

```text
q0 = q1 = q2 = published regular 600-cell boundary,
x1 = x2 = published time-symmetric internal dust sandwich,
```

with the second slab identified by the derived layer and orbit map.  The
shared derivative is **not expected to vanish**.  It must instead agree
componentwise with `p_post(first)-p_pre(second)`, with the nonzero sign and
factor fixed before evaluation.  No deformation, optimizer, fitted tolerance
or branch selection is allowed in this first control.

## 2. KNOWN structure

### KNOWN: implicit multi-step Regge evolution

Barrett, Galassi, Miller, Sorkin, Tuckey and Williams construct implicit tent
evolution by adjoining simplices to a spatial hypersurface and explicitly
apply the scheme to a dust-filled 600-cell.  Repeated tent moves and the
staircase subdivision are therefore established structure, not inventions of
this repository:

- [*A Parallelizable Implicit Evolution Scheme for Regge
  Calculus*](https://arxiv.org/abs/gr-qc/9411008).

De Felice and Fabri correct the 600-cell schedule to five classes of 24,
evolve those classes successively and publish the time-symmetric dust
sandwich used here:

- [*The Friedmann universe of dust by Regge Calculus: study of its ending
  point*](https://arxiv.org/abs/gr-qc/0009093).
- [*Singularities of the closed RW metric in Regge Calculus: a generalized
  evolution of the 600-cell*](https://arxiv.org/abs/gr-qc/0106077).

### KNOWN: action composition and momentum matching

Canonical simplicial gravity treats the discrete action as Hamilton's
principal function.  Gluing/removing simplices gives canonical evolution;
pre- and post-constraints and later consistency conditions control data
introduced by earlier moves:

- B. Dittrich and P. A. Hoehn, [*Canonical simplicial
  gravity*](https://arxiv.org/abs/1108.1974).

Thus adding adjacent discrete actions and extremizing their shared boundary
is the standard variational composition law.  The equation
`p_post=p_pre` is **KNOWN**, not a new dynamical principle.

### KNOWN: curvature can turn constraints into pseudo-constraints

On curved Regge solutions exact vertex-displacement gauge symmetry can be
broken, leaving pseudo-constraints that depend on next-step data:

- B. Bahr and B. Dittrich, [*(Broken) Gauge Symmetries and Constraints in
  Regge Calculus*](https://arxiv.org/abs/0905.1670).

This is a known reason why arbitrary boundary data need not admit a
stationary filling.  It motivates the two-slab consistency test but does not
predict its 600-cell result.

### KNOWN: causal cell type and boundary admissibility matter

In a different Lorentzian Regge cosmology, Jercher and Steinhaus find that a
height solution depends on causal regularity and an inequality relating
matter and geometric boundary data:

- A. F. Jercher and S. Steinhaus, [*Cosmology in Lorentzian Regge calculus:
  causality violations, massless scalar field and discrete
  dynamics*](https://arxiv.org/abs/2312.11639).

This blocks the inference that any prescribed pair of spatial boundaries
must have a Regge filling.

## 3. Repository controls already established

The following are **CONTROL**, not new hypotheses:

1. `reproducible/verify_gravity_global_boundary_legendre.py` certifies a real
   65-variable action with 35 internal equations, 30 final post-momenta and
   30 old pre-momenta (`33/33`).
2. Its regular pure-gravity pre/post momentum **multisets** obey time reversal:
   `p_pre` matches `-p_post`, not `p_post`.  This refutes the initial repeated-
   sandwich stationary-control proposal but does not yet give the required
   orbitwise dust-sandwich map.
3. `reproducible/verify_gravity_600cell_published_dust_control.py` reproduces
   all 35 complete internal dust equations in both parities (`14/14`).
4. The published symmetric control supplies a stationary one-slab factor, a
   known pre-momentum and a known opposite post-momentum.  It is a calibration
   point for canonical inversion, not a stationary repeated double.
5. The calibrated continuation result `94675be` finds 63 nonstationary
   deformed boundary-value states and leaves 17 numerically open.  It does not
   test composition.

## 4. OPEN object

The following have not been established in the repository or located in the
primary sources above:

- the explicit orbitwise time-reversal permutation for all 30 shared
  600-cell boundary edge orbits;
- signed orbitwise relation between all 30 dust-sandwich pre/post momenta
  under that map;
- direct cancellation of the two one-slab boundary hinge terms into the
  internal `2*pi` shared-hinge term;
- equality between the summed one-slab action and an explicitly reconstructed
  two-slab Regge action on the same Lorentzian branch;
- the rank and gauge structure of the full shared-slice equation;
- the adjacent non-time-symmetric slab whose pre-momentum equals the published
  slab's post-momentum;
- any multi-step scale sequence, continuum limit, causal speed or physical
  clock.

No located source prints the present order-24 orbit permutation or the
30-component momentum match.  That absence does not establish novelty.

## 5. Framing attack

A successful action-gluing and shared-derivative identity would be only a
**DERIVED CONTROL**.  It would prove that the repository can compose its
one-slab action with the correct signs and orbit identification.  Correctly
recovering a nonzero cusp for the repeated sandwich is part of that control,
not a failure.

A failure is load-bearing.  Until the orbitwise momentum signs and action
gluing pass, no second-slab evolution result may be physically interpreted.

The first test must also distinguish:

1. a wrong orbit permutation;
2. a boundary-sign or extrinsic-curvature mismatch;
3. a Lorentzian branch mismatch;
4. a dust normalization mismatch;
5. genuine failure of the published factors to compose.

It must not use an optimizer or search over permutations to repair a failure.
The permutation must be derived from vertex incidence and schedule reversal
before momenta are compared.

## 6. Decision and kill boundary

The next protocol will preregister only the structural composition and sign
control.  It must freeze:
It must freeze:

- the two canonical schedule parity representatives;
- the vertex-level reversal map;
- the induced 30-orbit permutation;
- the complete one-slab action and dust term;
- arbitrary-precision derivative steps and branch gates;
- the signed nonzero repeated-sandwich residual predicted from the separately
  computed one-slab momenta;
- orbitwise shared-derivative and action-gluing tolerances;
- exactly two attempts and no permutation search.

**Acceptance boundary:** both parity controls pass the direct action-gluing
identity, all branch gates and all 30 identities
`dS_2/dq_shared = p_post(first)-p_pre(second)`.  The repeated-sandwich
residual must be nonzero with the preregistered sign relation.

**Kill boundary:** if either parity fails after the geometric map, sign
conventions and complete action are independently audited, the present
one-slab objects cannot yet be composed into dynamics.  Stop before any
canonical continuation solve.

Only after this control passes may a separate preregistration perform
canonical continuation:

```text
given q_old and p_target = p_post(published slab),
solve the next slab for x and q_new such that
internal equations vanish and p_pre(next)=p_target.
```

That solve must include the known collective lapse/gauge issue and first
reproduce the published slab when `p_target=p_pre(published slab)`.  The
canonical continuation, not the structural gluing control, is the first
candidate for an actual next frame.
