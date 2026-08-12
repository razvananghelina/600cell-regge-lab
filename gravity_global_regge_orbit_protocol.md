# Preregistration: symmetry-reduced Regge system on the five-phase slab

Date: 2026-08-12

Upstream certified schedule commit: `d439b07`

Status: **preregistered before evaluating any global Regge gradient, Hessian or
root on this slab**.

The orbit counts in section 2 were seen in a scratch enumeration before this
commit and are disclosed.  No action value, gradient, Hessian rank or
stationary point has been computed.

No comparison with `a1=5`, a particle target, a coupling, `c`, `G`, Planck
time or Planck mass is permitted in this audit.

## 1. Complete hypotheses

1. Use the certified regular 600-cell as both initial and final boundary.
2. Choose one representative of the even phase-order orbit from the certified
   ten unordered covers.  Representatives inside that orbit are equivalent
   by `H4`; no individual labeling is physical.
3. Use its certified staircase triangulation with 2400 four-simplices.
4. Fix all 720 old and all 720 final boundary squared edge lengths to one.
5. Vary all internal squared lengths only in the fixed subspace of the exact
   order-24 automorphism group of this ordered slab.
6. Parameterize diagonal squared lengths as positive and pole squared lengths
   as negative.  Every four-simplex must retain inertia `(-,+,+,+)`.
7. Use the corrected plus complex-angle branch of Borissova--Dittrich already
   certified in `verify_gravity_lorentzian_tent_legendre.py`.
8. Use the ordinary zero-volume Lorentzian Regge action.  No cosmological,
   matter, higher-curvature or fitted coefficient is added.
9. Every non-boundary triangle has curvature constant `2*pi`; boundary
   triangles have `pi`.  The gradient with respect to internal edges is
   evaluated through the Schlaefli-reduced area derivative only after the
   full angle branch is independently checked.

This is a symmetry-restricted existence route.  A root in the fixed subspace
is a full root by finite-group invariance: gradients are constant on each edge
orbit, so a zero derivative with respect to every orbit variable makes every
one of the 840 internal edge derivatives zero.  Failure to find a root in
this subspace cannot exclude nonsymmetric roots in the full 840-variable
system.

## 2. Exploratory orbit counts already seen

The pointwise phase stabilizer appeared to act freely on every relevant
simplex layer.  The disclosed counts are:

```text
internal edges       840 = 35 orbits x 24
  poles              120 =  5 orbits x 24
  diagonals          720 = 30 orbits x 24
internal triangles  3840 = 160 orbits x 24
internal tetrahedra 5400 = 225 orbits x 24
four-simplices      2400 = 100 orbits x 24.
```

There is one pole orbit for each of the five phases.  For each of the ten
unordered phase pairs, its 72 diagonals split into three orbits of 24.  Thus
the honest invariant system has 35 variables, not the tempting 15 obtained by
silently equating those three diagonal orbits.

The registered verifier must reconstruct these counts.  No additional
equalities between the 35 variables are allowed unless an automorphism proves
them.

## 3. Frozen construction of the reduced action

For each orbit choose the lexicographically first representative only after
the complete orbit is constructed.  At invariant edge data:

- evaluate one representative from each of the 100 four-simplex orbits;
- assign its ten complex dihedral angles to the corresponding triangle
  orbits, with exact incidence multiplicity;
- assemble each internal triangle curvature from all incident simplex-orbit
  contributions;
- multiply each representative triangle contribution by its orbit size 24;
- differentiate its signed Heron area with respect to the 35 orbit variables.

No angle, area or edge coefficient may be fitted.  The reduced gradient is
the derivative of the action restricted to the invariant subspace; because
all edge orbits have size 24, division by 24 must reproduce the derivative at
each individual edge.

## 4. Frozen hostile controls

### A. Combinatorial reconstruction

- reconstruct both ordered parity representatives independently of the JSON;
- verify their automorphism stabilizers have order 24;
- reproduce every count in section 2 for both parity classes;
- verify all orbits have size 24.

### B. Reduced versus unreduced action

At each deterministic control below, evaluate all 2400 four-simplices and the
100-orbit reduction independently:

```text
R0: every diagonal = 1, every rho_phase = 1/4
R1: diagonal orbit j = 1 + (j+1)/1000,
    rho_phase k = 1/4 + (k+1)/1000
R2: diagonal orbit j = 1 - (j+1)/2000,
    rho_phase k = 1/4 + (5-k)/1500.
```

Every control must keep all simplices Lorentzian.  Compare:

- the full complex action;
- all 35 restricted derivatives;
- every internal triangle curvature grouped by orbit;
- simplex inertia and complex-angle argument margins.

Agreement tolerance is `2e-8` relative for action/gradient and `2e-9` for
orbit curvature values.  Failure is fatal to the reduced evaluator.

### C. Schlaefli and direct differentiation

At `R0`, `R1`, and `R2`, use centered differences of the complete complex
action in every one of the 35 orbit directions.  The result must agree with
the Schlaefli-reduced gradient to relative tolerance `2e-5`.  In addition,
test the per-simplex complex Schlaefli identity in deterministic orbit
directions, using the already corrected mixed-causal cosine sign.

### D. Known regular negative

At `R0`, recover the exact regular pole deficit and verify that none of the
five pole-orbit derivatives vanishes.  This is a negative control, not a root
search.

### E. Linearized system only

Compute the centered-difference `35 x 35` Hessian of the restricted action at
`R0`, separately for the two phase-parity representatives.  Report:

- real/complex residual;
- symmetry residual;
- singular values and numerical rank at frozen relative thresholds
  `1e-7`, `1e-9`, `1e-11`;
- the dimension and orbit support of any nullspace.

The Hessian is evaluated at a nonstationary point and therefore is not called
a physical propagator.  Its role is only to determine whether a subsequent
root search is locally well posed and whether phase parity changes the
linearized algebra.

## 5. Acceptance and kill boundaries

The **orbit reduction** advances only if A--C pass.  The regular negative in D
must reproduce the existing local theorem.  The linearized audit advances if
E gives a stable, interpretable rank under all three thresholds.

This protocol does not search for a root.  A root-search protocol may be
written only after the evaluator and its conditioning are committed.

Kill or downgrade conditions:

- mismatch between full and orbit-reduced gradients kills the reduction;
- a complex action/gradient at the real causal controls kills the claimed
  real Lorentzian branch there;
- loss of Lorentzian inertia kills that control, not the entire slab;
- singular or ill-conditioned linearization does not prove no dynamics, but
  kills a naive Newton continuation;
- parity-dependent spectra are a structural difference, not automatically
  physical chirality;
- no negative result in 35 variables may be promoted to a no-go for all 840
  variables.

Only the targeted verifier will be run.  The full suite remains outside this
mission.
