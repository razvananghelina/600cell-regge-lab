# Minimal first-class conversion relocates rather than removes nonlocality

Date: 2026-08-11

Preregistration commit: `7256f8f`  
Pre-computation sign correction: `9bae8d8`

Targeted verifier:
`reproducible/verify_whitney_first_class_conversion.py`

Targeted result: **11/11 PASS**.  The verifier is registered.  The full suite
was not run by explicit user request.

## Result

The minimal linear auxiliary conversion successfully changes the exact copy
constraints from second class to first class.  However, every exact
gauge-invariant physical coordinate and every quadratic gauge-invariant
Hamiltonian retaining the Whitney physical block requires the inverse of the
multiplier Gram matrix.

> **DERIVED RELOCATION NO-GO:** first-class conversion removes the nonzero
> constraint bracket but relocates the global solve into the unique physical
> dressing and Hamiltonian.  It does not produce a local Whitney tick.

This closes the minimal linear conversion only.  It is not a theorem against
all possible gauge theories.

## Exact conversion

For

\[
G=CM^{-1}C^*>0,

\]

add one complex auxiliary coordinate $\eta\in\mathbb C^r$ with bracket

\[
\{\eta,\eta^*\}=+iG.

\]

Then

\[
\Phi=Cu+\eta

\]

is exactly first class because

\[
\{\Phi,\Phi^*\}=-iCM^{-1}C^*+iG=0.

\]

The verifier checks this as an exact rational matrix identity.

The dimension count is also correct:

| level | local (n) | auxiliaries/constraints (r) | extended (n+r) | physical (n-r) |
|---|---:|---:|---:|---:|
| base | 9,000 | 6,360 | 15,360 | 2,640 |
| first barycentric | 216,000 | 153,120 | 369,120 | 62,880 |

Thus the conversion does not create physical modes after quotienting both the
first-class constraint surface and its gauge orbits.

## Sign correction

The original preregistration wrote the wrong sign in the dressed coordinate.
This was detected and corrected in commit `9bae8d8` before any computation.
With the frozen brackets, the gauge transformations are

\[
\delta u=-iM^{-1}C^*\epsilon,
\qquad
\delta\eta=+iG\epsilon.

\]

Therefore the invariant coordinate is

\[
\widetilde u=u+X\eta,

\]

not $u-X\eta$.

## Uniqueness of the nonlocal dressing

Gauge invariance requires

\[
XG=M^{-1}C^*.

\]

Because $G>0$, the solution is unique:

\[
X=M^{-1}C^*G^{-1}.

\]

On the first-class surface $\eta=-Cu$,

\[
\widetilde u
=(I-M^{-1}C^*G^{-1}C)u
=P_Du,

\]

exactly the previously derived Dirac projector.  The verifier checks both the
gauge-generator annihilation and this projector reduction exactly.

Consequently the inverse is not an artefact of a poor Hamiltonian choice.  It
is forced already by the requirement that a physical coordinate be
gauge-invariant and agree with $u$ on the auxiliary-zero slice.

## Unique quadratic Hamiltonian blocks

Let the extended Hermitian quadratic Hessian retain the fixed Whitney block:

\[
H_{\rm ext}=
\begin{pmatrix}A&B\\B^*&D\end{pmatrix}.

\]

Gauge invariance fixes

\[
BG=AM^{-1}C^*,

\]

so

\[
B=AM^{-1}C^*G^{-1}=AX.

\]

The remaining block is

\[
D=X^*AX.

\]

and hence

\[
H_{\rm ext}
=
\begin{pmatrix}I\\X^*\end{pmatrix}
A
\begin{pmatrix}I&X\end{pmatrix}.

\]

The complete gauge generator lies in its kernel exactly.  These equations are
linear and unique under the stated fixed-top-block hypothesis; no coefficient
search was performed.

## Exact small-control support

On the boundary-of-4-simplex control, $(n,r)=(75,45)$:

- exact dressing entries beyond the two endpoint tetrahedra of their
  constraint: **460**;
- exact Hamiltonian cross-block entries beyond those endpoints: **665**;
- exact physical-projector entries between different tetrahedra: **1,000**.

This is a support calibration, not a full-complex support census.  The general
load-bearing result is the unique appearance of (G^{-1}); the earlier
full-complex/refinement certificates independently show that the corresponding
reduced solve grows in depth.

## Canonicity caveat

The minimal conversion uses an independent constraint basis.  No canonical
independent basis exists without choosing spanning trees in the occurrence
graphs.  Retaining every canonical neighbour row instead produces a
reducible first-class system with multiplier relations and would require a
ghost-for-ghost hierarchy.

That hierarchy may be mathematically legitimate, but its BRST complex and
Hamiltonian are additional structures.  They are not selected merely by the
fact that the constraint kernel is conforming assembly.

## Physical verdict

The sequence of results is now sharp:

1. exact spectral geometry has a uniformly local constrained pencil;
2. the original constraints are second class;
3. a positive multiplier metric changes the spectrum and inserts a scale;
4. minimal first-class conversion restores gauge algebra but makes the exact
   physical dressing and Hamiltonian nonlocal.

So neither “make the multiplier dynamical” nor “call it gauge” produces the
missing causal tick.  A successful continuation needs genuinely new local
degrees of freedom and a Hamiltonian selected independently of the target
Whitney reduction.

## Status ledger

- **DERIVED:** exact minimal first-class constraint algebra.
- **DERIVED:** correct physical dimensions at both levels.
- **DERIVED:** unique gauge-invariant dressing containing (G^{-1}).
- **DERIVED:** unique quadratic cross block containing (G^{-1}).
- **DERIVED CONTROL:** exact remote support counts `(460,665,1000)`.
- **DERIVED NEGATIVE:** no local bypass in the minimal linear conversion.
- **OPEN:** reducible local BRST/ghost-for-ghost construction with a selected
  Hamiltonian.
- **OPEN:** nonlinear or differently embedded first-class theories.
- **OPEN:** finite-stiffness approximate conformity.
- **NOT CLAIMED:** physical time, mass, $c$, $\hbar$, Newton's $G$, or a Planck
  scale.

## Reproduction

```bash
/home/razvan/science/.venv/bin/python \
  reproducible/verify_whitney_first_class_conversion.py
```

Expected result: `11/11`.
