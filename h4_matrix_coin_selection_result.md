# H4 symmetry does not uniquely select a local matrix coin

Date: 2026-08-11  
Preregistration commit: `d19c453`

## Result

On the natural four-colour chamber carrier, H4 symmetry and the local metric
reduce the matrix freedom but do not select one coin.

> **DERIVED NONSELECTION:** chamber equivariance alone leaves the complete
> group `U(4)`, with 16 real parameters.

> **DERIVED METRIC NONSELECTION:** requiring the coin to commute with either
> preregistered local Gram matrix leaves `U(1)^4`, with four independent
> phases rather than one global phase.

There is nevertheless one useful coefficient-free candidate: reflection in
the one-dimensional zero mode of the Gram matrix,

\[
C_0=2P_{\ker G}-I_4.
\]

It is real, symmetric, exactly unitary to numerical precision and commutes
with `G`.  It is a **canonical candidate**, not a theorem that the dynamics
must choose it.

The targeted verifier passes `17/17` in about 1.4 seconds.  No isotropy or
continuum target was evaluated, and no full suite was run.

## Why H4 equivariance leaves `U(4)`

The 14,400 complete flags form one chamber orbit.  The finite Coxeter group
H4 also has order 14,400, so the chamber action is free and transitive.  A
pointwise coin field equivariant under this action is constant over the
orbit, but its value on the four-dimensional multiplicity fibre remains an
arbitrary unitary matrix.

Thus

\[
\mathcal H\simeq \mathbb C[H_4]\otimes\mathbb C^4
\]

contains a full multiplicity algebra `M_4(C)`.  Constancy is not uniqueness.

This conclusion is not evaded by merely declaring the fibre a spinor.  On a
free orbit there is no nontrivial little group at a chamber; an internal
group cocycle can be trivialized by transporting a fibre frame from one
base chamber.  A genuinely selecting spin constraint must therefore come
from additional structure such as a connection, Clifford action or
compatibility between neighbouring frames—not from chamber equivariance
alone.

## What the metric adds

For both frozen metric conventions, the four-by-four Gram spectrum consists
of one zero mode and three distinct positive modes:

| metric | eigenvalues |
|---|---|
| literal geodesic steps | `0, 0.00469401, 0.00778918, 0.02683574` |
| unit directions | `0, 0.59254899, 1.41020150, 1.99724952` |

All eigenspaces are one-dimensional.  Therefore

\[
\{X:[X,G]=0\}\cong\mathbb C^4,
\]

and its unitary subgroup is

\[
U(1)^4.
\]

The commutant dimension was obtained independently from the eigenvalue
multiplicities and from the nullity of the explicit 16-column commutator
map.  Both give complex dimension four.

The labelled Gram matrices were checked over all 14,400 chambers.  Their
maximum spreads are `3.75e-12` for literal steps and `1.47e-10` for unit
directions, so the census is not tied to a selected chamber.

## Complete discrete spectral census

With four distinct eigenspaces, every Hermitian unitary spectral function is

\[
C_\epsilon=\sum_{j=0}^3\epsilon_jP_j,
\qquad \epsilon_j=\pm1.
\]

After identifying `C` and `-C` as global-phase equivalents, there are

\[
N_{\rm spectral}=2^{4-1}=8
\]

classes for each metric convention.  Hence even the coefficient-free
involution class is a census, not a unique answer.

The preregistered member `C_0` is distinguished from the other seven because
zero is an algebraically distinguished eigenvalue.  That makes it a
defensible next candidate without looking at its dynamics.  It does not make
the other allowed phase choices disappear.

## Framing correction: the paper's block coin

Under the literal identity identification between the four Coxeter colours
and the four active component indices, imposing both

\[
[C,G]=0
\]

and the paper's restricted form

\[
C=I_2\otimes Y
\]

leaves only scalar matrices.  The same complex dimension one occurs for all
24 permutations of the four labels.  The zero-mode reflection does not have
that block form; its least-squares block residual is about `1.49` for literal
steps and `1.46` for unit directions.

This is **STRUCTURAL**, not a derived bridge to the published robust walk.
The Coxeter colours label geometric neighbour directions, whereas the
paper's four active components are algorithmic spin/ancilla channels whose
meaning changes across the three bond macros.  No canonical intertwiner
between these two four-dimensional spaces has been derived.  In a continuous
change to the Gram eigenbasis the same block-commutant dimension becomes two,
which demonstrates the basis dependence directly.

Therefore the block-collapse calculation cannot be used to claim that the
published walk has a uniquely selected coin.

## Physical verdict

The extra matrix freedom is real:

- H4 symmetry alone permits 16 continuous real parameters;
- adding the local scalar metric still permits four phases;
- the existing robust-component architecture has no derived identification
  with the metric colour space.

Consequently, choosing a coin because it produces a desirable cone would be
fitting.  The geometry does, however, provide the zero-mode reflection as a
specific, coefficient-free proposal on a **new four-colour flip-flop
carrier**.  It can be tested without tuning, provided that carrier change is
declared rather than presented as the old robust walk.

The next honest dynamical test is therefore narrow: preregister the
four-colour shift

\[
S\lvert k,i\rangle=\lvert s_i k,i\rangle
\]

and the designated literal-metric coin `C_0`, then measure connectivity and
metric propagation without comparing or optimizing the other seven spectral
sign classes first.

## Status ledger

- **DERIVED:** the chamber orbit is transitive and the fibre multiplicity is
  four.
- **DERIVED:** H4-equivariant constant coins form `U(4)`.
- **DERIVED:** both local Gram matrices have four simple eigenspaces.
- **DERIVED:** Gram-compatible unitary coins form `U(1)^4`.
- **DERIVED:** there are eight spectral involution classes modulo global
  sign for each metric.
- **DERIVED:** the zero-mode reflection is coefficient free, real, symmetric
  and unitary.
- **STRUCTURAL:** promoting that distinguished mathematical reflection to the
  physical coin.
- **STRUCTURAL:** identifying Coxeter colour space with the robust paper
  component space.
- **OPEN:** dynamics of the zero-mode reflection on the natural colour
  carrier.
- **NOT CLAIMED:** isotropy, a Dirac cone, a speed of light, mass or Planck
  scale.

## Reproduction

```bash
/home/razvan/science/.venv/bin/python \
  reproducible/verify_h4_matrix_coin_selection.py
```

Expected result: `17/17`.
