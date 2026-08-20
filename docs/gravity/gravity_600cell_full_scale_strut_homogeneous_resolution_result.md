# Primary result: unique homogeneous weak-pole canonical line

Date: 2026-08-20  
Status: **DERIVED STRUCTURAL/COMPUTATIONAL, primary only; adversarial replication required**

## Outcome

The preregistered targeted verifier reports

```text
HOMOGENEOUS_WEAK_POLE_LINE_UNIQUE
10/10 tests passed
```

Artifact SHA-256:

```text
70d7583756acdbee77893f98d57054ab074d9353a86247840cc1eb2c7b6be931
```

Verifier source SHA-256:

```text
e1e1ae5b0244837bd59f550b1c7a9be6bbd42e5882951617b4cc44a7fc13ea48
```

Only this targeted verifier was run.  The full suite was not run.

## Exact statement

Under the hypotheses frozen in the protocol—the real nondegenerate complete
scale+strut carrier, the closed homogeneous Regge action, the frozen canonical
endpoint and the weak-pole canonical graph—the homogeneous generator is

```text
sigma = -lambda p_z,
c     =  p_s,
```

where `sigma=delta lambda`, `delta s=sigma/lambda`, and `p_s,p_z` are the two
logarithmic derivatives of the old collective momentum.  It obeys identically

```text
p_s (sigma/lambda) + p_z c = 0.
```

The exact carrier responses are the corresponding homothetic differentials.
The frozen exhaustive minors give rank 9 for the ten-column D block and rank 14
for the fifteen-column K block in both parities.  Consequently both kernels are
exactly one-dimensional.  The even/odd analytic generator projectors agree
exactly in the evaluated representation.

The physical carrier ratio is

```text
sigma/c = -2.3370425140082096076325345290206321744897033430589695e-6.
```

It agrees with the independently stored endpoint-Jacobian construction to
`5.34e-42`, inside its frozen error bound `2.29e-25`, and with the disclosed P100
near-null vector to `2.33e-42`.  Neither comparison supplies the generator or a
fitted coefficient; both identify the already derived line.

## Disclosed method failure

The first execution of source commit `5889c2e` produced a bundled symbolic FAIL
and was interrupted after about 11 minutes in a global SymPy factorization.  No
artifact was written.  The failure and an implementation-only correction were
frozen in commit `3ca01d0` before source modification.  The correction proves
the rational identities in algebraically independent derivative symbols and
uses the full action derivatives only for the automatic P160/P220 bridge.  No
formula, input, tolerance or outcome criterion changed.

## What this does and does not establish

- **DERIVED:** every already replicated nonhomogeneous canonical intersection is
  zero; the primary homogeneous intersection has exactly one line.
- **STRUCTURAL:** the line is the infinitesimal weak-lapse response that preserves
  the old collective momentum in the reduced canonical graph.
- **OPEN:** a mechanically different adversarial replication of the homogeneous
  result.
- **OPEN:** whether the omitted pole/lapse equation preserves or eliminates the
  line.
- **NOT ESTABLISHED:** gauge freedom, a solution tangent, propagation, a physical
  clock, a tick, `c`, `G` or a Planck scale.

Thus this is a kinematic/canonical uniqueness result, not yet dynamics.

## Post-result primary-literature check

The search was repeated after the outcome.  Dittrich and Höhn derive the general
generating-function framework and the possibility that discrete constraints fix
formerly free data, while their nonlinear analysis explains pseudo-constraints
from broken gauge symmetry:

- B. Dittrich and P. A. Höhn, *From covariant to canonical formulations of
  discrete gravity*, <https://arxiv.org/abs/0912.1817>.
- B. Dittrich and P. A. Höhn, *Canonical simplicial gravity*,
  <https://arxiv.org/abs/1108.1974>.
- P. A. Höhn, *Canonical linearized Regge Calculus: counting lattice gravitons
  with Pachner moves*, <https://arxiv.org/abs/1411.5672>.

The 600-cell evolution literature studies homogeneous dust evolution and its
causality-breaking endpoint, including generalized variables, but the searched
primary sources did not expose this exact complete-carrier intersection or its
one-dimensional kernel:

- A. De Felice and E. Fabri, *The Friedmann universe of dust by Regge Calculus:
  study of its ending point*, <https://arxiv.org/abs/gr-qc/0009093>.
- A. De Felice and E. Fabri, *Singularities of the closed RW metric in Regge
  Calculus: a generalized evolution of the 600-cell*,
  <https://arxiv.org/abs/gr-qc/0106077>.
- R. G. Liu and R. M. Williams, *Regge calculus models of closed lattice
  universes*, <https://arxiv.org/abs/1502.03000>.

The inference “no exact prior result was found” is search-limited.  **External
novelty remains OPEN.**

