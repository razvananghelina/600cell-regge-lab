# Compatible canonical-data projection census

Date: 2026-08-19

## Complete hypotheses

Use the frozen complete variable-face flat-frustum compatibility matrix
`[F E S]` on the fixed lower regular 600-cell, with 3600 cell-flex columns,
720 upper spatial squared-edge columns, and 120 strut squared-length columns.
Evaluate the two nonstatic rational representatives, two primes, two exact
right-inverse graphs, reversed face orientation, odd canonical relabelling,
and reversed metric sign.  This is a kinematic linear compatibility test, not
an action, Hessian, constraint algebra, or evolution equation.

## Frozen result

The registered verifier passed 11/11 checks.  Its first JSON was committed as
`9b97775` with the explicit message that no carrier comparison had been
performed.  The artifact SHA-256 is
`f011ef9848a6139408a9f8495a12e0d8e0050e04f39aa5c00ca88c02dde26beb`.

Every prime and all seven complete constructions give

```text
rank(F)       = 3600
rank([F E])   = 4200
rank([F S])   = 3600
rank([F E S]) = 4200

dim compatible data              = 240
dim edge-only compatible data    = 120
dim strut-only compatible data   = 120
dim edge projection              = 120
dim strut projection             = 120
```

## Exact modular consequences

Let `K_p` be the compatible data space over either frozen finite field.  Let
`K_E=K_p intersect (E direct-sum 0)` and
`K_S=K_p intersect (0 direct-sum S)`.  Their intersection is zero and

```text
dim(K_E) + dim(K_S) = 120 + 120 = dim(K_p).
```

Therefore

\[
K_p=K_E\oplus K_S.
\]

**DERIVED (modular).** The compatible space splits into independent edge-only
and strut-only sectors.  Since the full strut ambient space has dimension 120,
`K_S` is that entire ambient space: every infinitesimal strut datum is
kinematically compatible modulo both primes.  The edge sector is a
120-dimensional subspace of 720 edge data.

**OPEN.** No equality over the rationals follows from two primes alone.  No
specific basis of `K_E`, rational lift through `F`, symplectic form, or
dynamical interpretation has yet been proved.

## Post-result prior-art check and framing attack

Free lapse/shift-type vertex data on a flat linearized Regge background are
not new.  Hoehn shows that the linearized 1-4 Pachner move generates four
lapse-and-shift variables and their vertex-displacement generators; see
[arXiv:1411.5672](https://arxiv.org/abs/1411.5672), especially the 1-4 move
analysis.  Thus the full 120-dimensional strut sector is structurally
consistent with known flat-background gauge freedom; it is not evidence for a
new physical clock.

The warning is equally important in the other direction.  Brewin's Regge
Friedmann analysis explains that in constrained nonlinear cosmological models
struts and diagonals need not remain arbitrary and may have to be evolved with
the spatial geometry; see DOI
[10.1088/0264-9381/4/4/023](https://doi.org/10.1088/0264-9381/4/4/023),
especially the discussion of constraints versus variation.  Our result is
only flat-frustum compatibility and cannot settle the action-level lapse
question.

Vertex-based conformal edge scalings are established mathematical prior art,
including piecewise-flat three-manifolds and the Regge scalar-curvature
functional; see Glickenstein
[arXiv:0906.1560](https://arxiv.org/abs/0906.1560).  This makes a
120-vertex conformal candidate natural, but dimension 120 does not prove that
the computed `K_E` is that subspace.

External novelty of the exact 600-cell projection census remains **OPEN**.

## Protocol deviation

Before the JSON was written, but after the full source and criteria were
committed, the assistant disclosed the speculative target `119+121`.  Commit
`a3fe2b9` records this while the artifact was still absent.  The result instead
returned `120+120`, so the disclosed numerical guess was falsified.  Even so,
the strict blindness claim is not restored: the census is **STRUCTURAL**, and
any carrier identification requires mechanically independent exact evidence.

## Next falsifiable distinction

The old failed vertex scale/lapse construction specified both boundary data
and a particular local lift into the 3600 cell-flex variables.  Its failure
refutes that augmented lift.  It does not by itself refute the possibility
that its 240-dimensional *data image* equals `K_p` with a different, globally
solved lift.

The next test must therefore distinguish these claims:

1. construct the 720-by-120 unsigned vertex-edge incidence image dictated by
   the derived squared-length variation on an edge;
2. test whether this exact data subspace equals `K_E`, using image inclusion
   plus the already frozen equal dimension, not dimension alone;
3. treat the 120 strut coordinates as the full `K_S` candidate;
4. only after data-space equality, construct and verify the unique rational
   cell-flex lifts through the full-column-rank block `F`.

Until steps 2 and 4 pass independently, “scale plus lapse” remains **OPEN** at
the data level and **DERIVED NEGATIVE** for the original local augmented lift.
