# Adversarial protocol: global prism shift gluing

Date: 2026-08-19

Primary result commit: `4d4377b`.

Primary artifact SHA-256:

```text
1ab6654ae57c83a49dd4f427154b891c0b8ae613631773ab6733a1227b9999fa
```

The primary verifier returned a 119-dimensional face-matched kernel and a
finite nonuniqueness family.  This audit is designed to refute that result,
not reproduce its decisive matrix.

## 1. Independence requirement

The audit must not import the primary verifier, its local `2400 x 1800`
matching matrix, its row builder, its `F_2` elimination or its disclosed
one-vertex potential.

Instead it must use the graph's fundamental-cycle basis:

1. construct a deterministic spanning tree of the 600-cell edge graph;
2. use the `E-V+1=601` non-tree edges as coordinates on graph cycles;
3. project every oriented triangular boundary into those chord coordinates;
4. compute the triangle-boundary span over both `F_3` and
   `F_1000003` with a sparse elimination independent of the primary bitset
   code.

If the triangle boundaries fail to span all 601 cycle coordinates over
either field, the claimed vanishing of first cohomology is not accepted.

## 2. Controls

- The boundary of one tetrahedron must have cycle dimension three and
  triangle-boundary rank three over both primes.
- The periodic `3 x 3` triangulated torus must have cycle dimension 19 and
  triangle-boundary rank 17 over both primes, leaving `b1=2`.
- Removing all triangular relations from the 600-cell must leave all 601
  graph cycles.  This guards against a rank routine that silently assumes
  exactness.

## 3. Independent reconstruction of the 119 modes

Do not form local tetrahedral covector unknowns.  Assign an edge number

```text
x_(u,v)=phi(v)-phi(u)
```

from a vertex potential.  Use path integration on the spanning tree to
reconstruct `phi` from `x`, then verify every non-tree edge.  Repeat for at
least five deterministic integer potentials, including:

- graph distance from a seed;
- squared graph distance;
- a signed coordinate-order potential;
- two deterministic modular polynomials of the vertex label.

The reconstruction must be exact over integers.  Adding a constant to a
potential must leave every edge value fixed.

## 4. Different finite Lorentzian family

Use the squared-distance potential, not the primary one-vertex potential,
and set `rho=7/5`, not `rho=1`.

For every tetrahedron independently form

```text
H=[G a; a^T -rho].
```

Verify its `(3,1)` signature through the positive spatial Gram matrix and
negative Schur complement.  For every shared triangle, independently build
the complete labelled `6 x 6` squared-interval matrix of its lateral prism:

```text
d^2(bottom_i,bottom_j)=d^2(top_i,top_j)=ell_ij^2,
d^2(bottom_i,top_j)=ell_ij^2+2(phi_j-phi_i)-rho.
```

The two incident four-cells must produce the same matrix to `1e-10`.
At least one four-volume must differ from the zero-potential value while the
16 natural lengths of each cell remain unchanged.

## 5. Attack on the matching convention

Also evaluate the stronger ambient-vector convention rejected before the
primary protocol.  Embed the regular 600-cell in `R4`, compute the
intersection of all 600 tetrahedral direction spaces, and report its
dimension.

If that intersection is zero, record explicitly:

- full ambient-vector matching would kill all modes;
- it is not equivalent to labelled lateral-face metric matching because it
  supplies a common ambient frame absent from an intrinsic Regge complex.

If the face metric actually determines that missing normal component, the
primary definition is too weak and the 119-mode conclusion is refuted.

## 6. Verdict

Return

```text
GLOBAL_SHIFT_GLUING_CORROBORATED
```

only if both prime-field cycle calculations, all controls, exact potential
reconstructions and the independent finite face metrics pass.  Return
`PRIMARY_RESULT_REFUTED` on an exact disagreement and `PRECISION_OPEN` if
only floating geometry is unresolved.

Even a passing audit establishes only a **DERIVED KINEMATIC** ambiguity.  It
does not classify the modes as gauge and does not authorize an ADM or
graviton claim.

Only this targeted audit and static guards may run; the full suite is
excluded.
