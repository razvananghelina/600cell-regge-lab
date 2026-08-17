# Prior-art gate: universal staircase overlay of a tetrahedral prism

Date: 2026-08-17

Status: written before enumerating any chamber of the proposed arrangement.

## 1. Exact object and hypotheses

Let `Delta^3` have barycentric coordinates `lambda_0,...,lambda_3`, with

```text
lambda_i >= 0,  sum_i lambda_i = 1,
```

and let `t in [0,1]` be the interval coordinate.  The local carrier is the
tetrahedral prism `Delta^3 x I`.  For every nonempty proper subset
`A subset {0,1,2,3}`, introduce the rational hyperplane

```text
h_A = t - sum_(i in A) lambda_i = 0.
```

There are exactly 14 such hyperplanes.  Their induced polyhedral arrangement
inside the prism is the proposed **universal staircase overlay**.  The global
carrier would be obtained by applying this construction functorially to every
tetrahedron of the 600-cell boundary and gluing along common faces.  A
barycentric subdivision of the resulting polyhedral complex would be a
simplicial carrier.

This mission contains no Regge action, metric, dust distribution, continuum
target, physical constant or fitted parameter.

## 2. Why these hyperplanes occur

For a vertex order `(v_0,v_1,v_2,v_3)`, the four standard staircase
four-simplices are the regions indexed by `k=0,...,3` satisfying

```text
sum_(i>k) lambda_(v_i) <= t <= sum_(i>=k) lambda_(v_i).
```

Thus every internal facet of every staircase triangulation lies on an
`h_A=0` with `A` nonempty and proper.  Conversely, every such subset occurs as
a tail of some order.  Therefore the arrangement of all 14 hyperplanes is,
structurally, the polyhedral overlay of all 24 vertex-order staircases.  The
planned exact computation audits rather than assumes that its open chambers
are assigned uniquely to one simplex of every staircase.

## 3. Primary prior art

- Staircase triangulations belong to the general theory of triangulations of
  products of simplices and the Cayley trick; see Francisco Santos,
  [*The Cayley trick and triangulations of products of
  simplices*](https://arxiv.org/abs/math/0312069).  This supplies the known
  product/staircase setting, not the requested chamber census.
- Arrangements whose normals are all nonzero `0/1` vectors are known as the
  resonance or all-subsets arrangement; see Lukas Kühne,
  [*The Universality of the Resonance Arrangement and its Betti
  Numbers*](https://arxiv.org/abs/2008.10553), and Zachary Chroman and Mihir
  Singhal, [*Computations associated with the resonance
  arrangement*](https://arxiv.org/abs/2106.09940).  The present arrangement is
  an affine restriction inside a simplex prism, not the central resonance
  arrangement itself.
- Karim Adiprasito and Igor Pak prove that PL-homeomorphic triangulations have
  a common stellar subdivision in [*All triangulations have a common stellar
  subdivision*](https://arxiv.org/abs/2404.05930).  This is an existence
  theorem and does not select this overlay or a minimal refinement.
- An improved or perfect Regge action requires a dynamical coarse-graining
  construction, not merely a common carrier; see Benjamin Bahr and Bianca
  Dittrich, [*Improved and Perfect Actions in Discrete
  Gravity*](https://arxiv.org/abs/0907.4323).

No primary source located in this gate gives the exact chamber and symmetry
orbit census of this 14-hyperplane restriction in `Delta^3 x I`.  That absence
does not prove novelty; external novelty remains **OPEN** pending a dedicated
review.

## 4. KNOWN / CONTROL / OPEN

### KNOWN

- Each of the 24 vertex orders gives a standard four-simplex staircase
  triangulation of `Delta^3 x I`.
- All its internal facets are among the 14 subset-sum hyperplanes above.
- Hyperplane arrangements cut a convex polytope into convex polyhedral cells.
- Any two PL-equivalent triangulations admit some common stellar subdivision.

### CONTROL

- Exhaust all `2^14` strict sign patterns exactly; no sampling or floating
  tolerance is allowed.
- Check the full local `S4 x C2` action, where time reflection sends
  `sign(h_A)` to `-sign(h_(A^c))`.
- For every feasible chamber and every one of the 24 orders, require a unique
  containing staircase simplex.
- On each spatial face `lambda_j=0`, require the hyperplanes to restrict to
  the corresponding six nontrivial subset hyperplanes for a triangular prism,
  plus only boundary duplicates `t=0,1`.  This is the gluing control.

### OPEN

- The exact number of full-dimensional chambers and their `S4 x C2` orbits.
- Whether all exact assignment and face-restriction controls pass.
- The full face vector of the global overlay or its barycentric subdivision.
- A Lorentzian realization, a dust/Regge action, coarse-graining weights,
  canonical momenta and any nontrivial evolution.
- Whether the resulting carrier is physically preferred over another
  functorial common refinement.

## 5. Framing attack

Even a successful census would establish only a canonical **common carrier**.
It would not make the old even and odd actions equal and would not be a perfect
action.  The construction uses all 24 local schedules symmetrically, so it
removes the arbitrary choice of one schedule, but it may be combinatorially
large.  A large chamber count would be a practical negative for the proposed
next dynamical calculation and must be reported without relabelling it as
physical progress.

## 6. Post-result terminology correction

The exact census exposed a sharper classical identification that the initial
gate missed.  A chamber sign word defines a Boolean function on subsets by

```text
f(A)=1  iff  sum_(i in A) lambda_i > t.
```

It is therefore a nonconstant positive (monotone) threshold Boolean function
on four variables.  Conversely, strict positive weights can be normalized to
`sum lambda_i=1`, and zero weights of inessential variables can be perturbed
positively inside a strict separation margin.  Hence the local chamber census
is exactly the classical census of positive four-variable threshold functions,
with the two constant functions removed.

The classical count is 150 including the constants, so the computed 148 is an
external combinatorial cross-check, not a new count.  Early enumerations are
summarized by Saburo Muroga, *Threshold Logic and Its Applications* (Wiley,
1971), Table 2.3.2, and the associated references S. Muroga, I. Toda and
M. Kondo, *Majority decision functions of up to six variables*, Mathematics
of Computation 16 (1962), 459--472.  The current tabulation and references are
also indexed by [OEIS A002078](https://oeis.org/A002078).

This missed keyword does not change the preregistered common-refinement test,
but it changes the novelty ledger: the number 148 is **KNOWN**, while its use
as the local universal staircase overlay in this repository is
**STRUCTURAL / external novelty OPEN**.
