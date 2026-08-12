# The metric phase space is canonical; its Hamiltonian is not

Date: 2026-08-12

Preregistered protocol commit: `3531d9a`

Registered verifier:
`reproducible/verify_gravity_metric_phase_space_canonicity.py`

Machine-readable result:
`reproducible/gravity_metric_phase_space_canonicity.json`

## Headline

> **DERIVED H4 KINETIC CANONICITY OBSTRUCTION.** The 720 Regge edge lengths
> have a canonical 1440-dimensional cotangent phase space, but the full `H4`
> symmetry leaves 47 independent real symmetric quadratic kinetic forms.
> Nearest-neighbour locality leaves four parameters, and locality within one
> tetrahedron leaves three. Positivity leaves open families in both classes.

Thus the current geometry supplies the arena for metric dynamics but does not
select the dynamics. The targeted verifier passes `16/16`. No full-suite run
was performed.

The value 47 was seen in an exploratory calculation before preregistration and
was disclosed in the protocol. It is a certified result, not a blind
prediction and not a phenomenological target hit.

## 1. What is canonical

Let `M_Regge` be the open set of admissible positive squared lengths `x_e` on
the 720 edges, locally around the equilateral 600-cell. Its cotangent bundle
has coordinates `(x_e,p_e)` and the coordinate-independent canonical form

```text
omega = sum_e d p_e wedge d x_e.
```

Therefore the following are **DERIVED STRUCTURALLY**:

- 720 metric configuration variables;
- 720 conjugate covectors;
- a 1440-dimensional phase-space arena;
- the cotangent lift of the `H4` edge-permutation action, preserving `omega`.

Changing from squared lengths to lengths induces the corresponding canonical
change of momenta; it is not a different symplectic theory.

This construction supplies no Hamiltonian, lapse, constraint, time unit,
Lorentzian signature, `c`, `G` or Planck normalization.

## 2. Exact invariant-operator count

At the equilateral point write the momentum-quadratic term as

```text
T_K(p) = (1/2) p^T K p,       K=K^T.
```

For a transitive finite group action on a set `X`, equivariant matrices are
constant on the group orbits of ordered pairs `X x X`. Equivalently, their
dimension is the number of stabilizer orbits on `X`.

Here `X` is the set of 720 unoriented edges. The verifier independently
rebuilds the exact finite permutations from the quaternion action and obtains

```text
edge orbit size          = 720
edge stabilizer order    = 20
implied |H4|             = 14,400
ordered-pair orbitals    = 62
```

The stabilizer-orbit size multiset is

```text
size 1  :  2 orbits
size 2  :  4 orbits
size 5  :  6 orbits
size 10 : 32 orbits
size 20 : 18 orbits.
```

Transpose cannot be inferred from the enumeration order. For each orbital the
verifier constructs an explicit group element sending the base edge to an
orbital representative and uses its inverse to identify the transposed
orbital. The involution has

```text
32 self-transpose orbitals
15 pairs of distinct transpose orbitals,
```

so

```text
dim Sym(End_H4(R^720)) = 32+15 = 47.
```

This is the exact linear freedom in an `H4`-invariant quadratic kinetic term
at the equilateral point.

## 3. Locality does not select a ray

Using distance in the line graph of the 600-cell edges, the 47 symmetric
parameters split as

| exact edge-line distance | parameters |
|---:|---:|
| 0 | 1 |
| 1 | 3 |
| 2 | 9 |
| 3 | 14 |
| 4 | 13 |
| 5 | 7 |

Consequently a nearest-neighbour kinetic form, supported only when two edges
coincide or share a vertex, still has dimension

```text
1+3 = 4.
```

An alternative cell-local requirement, support only when the two edges lie in
a common tetrahedron, has dimension 3. These two locality notions are not
nested: opposite edges of a tetrahedron do not share a vertex, while some
edges sharing a vertex do not lie in a common tetrahedron. The word “further”
in item 9 of the protocol should therefore be read as “alternative stronger
cell support”, not as a subspace inclusion. The registered predicates were
the explicit support relations and are unchanged by this wording correction.

Positivity does not restore uniqueness. Let `A_line` be the line-graph
adjacency and `A_tet` the adjacency of distinct edges lying in a common
tetrahedron. Their exact degrees are 22 and 15. Hence

```text
I + (1/44) A_line
I + (1/30) A_tet
```

are symmetric, invariant, nonproportional to `I`, and positive definite by
strict diagonal dominance. Varying either coefficient in an interval around
zero gives a continuum of positive rays.

Only the ultralocal condition `K_ef=0` for every `e != f` reduces the class to
`K proportional to I`. The geometry has not derived that condition. Imposing
it now because it creates uniqueness would be an ansatz.

## 4. Attack on the conclusion

This obstruction is deliberately limited.

1. It proves non-selection by the current spatial symmetry, incidence
   locality and positivity. It does not prove non-selection after a temporal
   slab or four-dimensional Regge action is supplied.
2. A first-class constraint algebra could restrict the 47 parameters. That
   algebra is precisely what the current theory lacks, so it cannot yet be
   used as a selector.
3. A DeWitt/Regge supermetric could be imported from general relativity, but
   then the dynamics comes from that external input rather than from the
   present 600-cell data.
4. The finite `A2` Hessian is a candidate potential curvature, not by itself a
   kinetic supermetric. It is also indefinite on the scale quotient, so
   silently using it as `K` would create a new physical assumption.

The obstruction is therefore real but not terminal: the missing selector
must involve inter-slice geometry, constraint closure or another independent
principle that is absent from a single spatial 600-cell.

## 5. Status ledger

| Claim | Status |
|---|---|
| `T*M_Regge` and its canonical symplectic form exist | **DERIVED STRUCTURALLY** |
| `H4` selects a unique quadratic kinetic form | **REFUTED** |
| Full symmetric invariant dimension is 47 | **DERIVED COMPUTATIONAL** |
| Nearest-neighbour invariant dimension is 4 | **DERIVED COMPUTATIONAL** |
| Common-tetrahedron invariant dimension is 3 | **DERIVED COMPUTATIONAL** |
| Positivity selects one of those local forms | **REFUTED** |
| `K=I` is the unique ultralocal invariant ray | **DERIVED CONDITIONAL** |
| Ultralocality is selected by the current geometry | **OPEN / NOT DERIVED** |
| A temporal/slab construction can select `K` | **OPEN** |
| Constraint closure can select `K` | **OPEN** |
| Lorentzian gravity, `c`, `G` or Planck units follow | **OPEN** |

## 6. What comes next

The next admissible construction must add data that a spatial symmetry cannot
see. Candidate classes must be enumerated before testing their physical
consequences:

1. a one-step four-dimensional slab between two 600-cell slices;
2. the allowed lapse/shift carriers on that slab;
3. a discrete action whose Legendre transform produces `K`;
4. primary/secondary constraints and their exact Poisson algebra;
5. only then, whether the 150 negative Euclidean `A2` directions are gauge,
   constrained or propagating.

This is a sharper target than “find dynamics”: derive the extra principle
that cuts a 47-dimensional kinetic family to a physical ray without looking
at the desired continuum answer.

## 7. Reproduction history

The first targeted run passed every combinatorial count but failed a deliberately
strict `1e-12` coordinate-residual guard at `1.12e-10`. The shared 600-cell
utility rounds its vertices to ten decimals and already has unit-norm residual
`5.60e-11`. The verifier was made more independent by rebuilding the standard
`Q(sqrt(5))` coordinate set without decimal rounding. The residual became
`2.22e-16`, while the nearest-wrong-vertex gap remained `0.190983`. No orbit,
transpose, locality or dimension count changed. The final targeted run passes
`16/16`.

No full suite and no PDF build were run.
