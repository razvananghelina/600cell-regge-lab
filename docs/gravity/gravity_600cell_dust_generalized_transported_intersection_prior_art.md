# Prior-art and framing gate: transported generalized phase intersection

Date: 2026-08-18

## Exact object, carrier and complete hypotheses

For each of the two disclosed one-dimensional symmetry sectors, both schedule
parities and all four derivative schedules, let

```text
F0 = E_old direct-sum E_old*       subset C^60,
F1 = E_shifted direct-sum E_shifted* subset C^60
```

be the already residual-certified rank-30 phase fibers, and let `T_2` be the
regular, symplectic second-slab tangent derived from the fixed Regge--dust
action and identity seam map.

The full inclusion `T_2 F0 subset F1` is now a replicated computational
refutation.  The next canonical object is

```text
K0 = F0 intersection T_2^-1(F1)
   = kernel((I-Q1) T_2 |_F0).
```

This mission asks only for the exact dimension of `K0`.  No desired dimension,
configuration graph, Lagrangian condition or continuum target may be inspected
before the complete rank census is committed.

## Primary literature

**KNOWN.** Principal angles and their sine/cosine SVD characterizations detect
intersection and separation of finite-dimensional subspaces; zero principal
angles correspond to common directions.  Knyazev and Argentati also give
perturbation estimates:
<https://doi.org/10.1137/S1064827500377332>.

**KNOWN.** Knyazev, Jujunashvili and Argentati characterize subspace angles,
projectors and corresponding Ritz values:
<https://doi.org/10.1016/j.jfa.2010.05.018> and
<https://arxiv.org/abs/0705.1023>.

**KNOWN.** In discrete variational systems with constraints, Dittrich and
Hoehn derive evolution between pre- and post-constraint surfaces as a
pre-symplectic map and show that propagating degrees of freedom depend on both
endpoints:
<https://doi.org/10.1063/1.4818895> and
<https://arxiv.org/abs/1303.4294>.

**KNOWN.** Discrete action principles generate the endpoint canonical map;
this continues to justify using `T_2`, not a separately selected transport:
<https://doi.org/10.1017/S096249290100006X>.

These results justify the object and rank diagnostic.  They do not predict the
rank in this 600-cell calculation.  External novelty is **OPEN**.

## Rank logic

Let

```text
R = (I-Q1) T_2 Q0.
```

The exact right factor `Q0` has rank 30, hence `rank(R)<=30`.  Moreover

```text
dim(K0) = 30-rank(R).
```

Therefore resolving 30 nonzero singular values of `R` certifies `K0={0}`
without needing to certify any numerical zero among the remaining 30 ambient
singular values.  Fewer than 30 resolved singular values do not by themselves
certify a positive intersection: small singular values consistent with zero
leave the exact nullity **OPEN** unless an additional exact identity supplies
an upper rank bound.

The error bound for `R` is already supplied by the accepted Flint/projector
calculation.  Weyl's singular-value perturbation bound permits every midpoint
singular value larger than that operator error to certify a nonzero exact
singular value.  The protocol will retain the established conservative factor
100 rather than choose a threshold after seeing the spectrum.

## Attack on the framing

The intersection is basis-independent and canonical once `F0,F1,T_2` are
fixed.  It is nevertheless conditional on the generalized-pencil definition
of `E_t` and on the four finite derivative schedules.  It is not the complete
pre/post-constraint surface of Regge gravity unless later derived from the
action itself.

- If `K0=0` in all cells, the present generalized-fiber phase route is closed;
  this does not kill other Regge perturbations.
- If `K0` is nonzero, that alone is not a propagator or a physical mode.  Its
  graph, symplectic/isotropic/Lagrangian character, schedule stability and
  next-slab transport must be tested separately.
- Choosing a lower-dimensional graph after seeing singular vectors would be
  fitted and is forbidden.

## KNOWN, CONTROL, OPEN

- **DERIVED UPSTREAM:** 16 canonical phase tangents and 32 residual-certified
  configuration projectors.
- **DERIVED UPSTREAM:** full phase transport fails in all 16 cells and is
  adversarially corroborated.
- **CONTROL:** replay the accepted high-precision phase verifier byte-for-byte
  and use its complete error bound.
- **OPEN:** exact rank of `R` and dimension of `K0`.
- **OPEN:** any action-selected constraint interpretation of a nonzero `K0`.
- **FORBIDDEN HERE:** inspecting desired ranks, fitting graphs, propagators,
  dispersion, mass, inertia, limiting speed or particle claims.

