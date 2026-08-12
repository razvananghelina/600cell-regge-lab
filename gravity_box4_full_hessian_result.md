# Corrective result: the full Hessian of `Tr(Box^4)`

Date: 2026-08-12  
Protocol commit: `d909513`  
Registered verifier: `reproducible/verify_gravity_box4_full_hessian.py`  
Machine-readable result: `reproducible/gravity_box4_full_hessian.json`

## Headline

**DERIVED CORRECTION.**  The matrix previously reported as the Hessian of
`Tr(Box^4)` is not the full Hessian.  It is a diagonal-eigenvalue-sensitivity
Gram matrix.  It depends on the arbitrary choice of eigenvectors inside
degenerate eigenspaces.

For the complete hypothesis list frozen in protocol commit `d909513`, the
actual `720 x 720` Hessian has exact rational rank `720`.  The spectral
Hessian formula makes it positive semidefinite, so exact full rank makes it
positive definite on the 720-dimensional edge-weight perturbation space:

```text
positive = 720, zero = 0, negative = 0.
```

The former `101 positive + 619 zero` result is reproduced, diagnosed, and
withdrawn as a Hessian or graviton statement.

## Complete hypotheses

The calculation uses the fixed real symmetric chamber-independent vertex
operator

```text
B = 5 A_fiber - A_cross
```

on the 120 vertices of the 600-cell.  The twelve derived Hopf decagons contain
120 undirected fiber edges; the remaining 600 graph edges are cross edges.
For every undirected edge `e={i,j}`,

```text
E_e = c_e (|i><j| + |j><i|),
c_e = 5  for a fiber edge,
c_e = -1 for a cross edge,
B(w) = B + sum_e w_e E_e.
```

All 720 real edge weights are varied independently.  No gauge quotient,
metric dictionary, physical target, graviton ansatz, continuum limit, or
phenomenological number enters the computation.

## Exact derivative

Let `F(w)=Tr(B(w)^4)`.  Noncommutative differentiation gives

```text
H_ef
 = 4 Tr(E_f B^2 E_e + B E_f B E_e + B^2 E_f E_e)
 = 8 Tr(B^2 E_f E_e) + 4 Tr(B E_f B E_e).
```

The second equality uses cyclicity of trace and transposition of the real
symmetric factors.  Every entry is an integer.  Direct dense trace checks and
three exact degree-four polynomial directional checks have zero residual.

## Exact rank and inertia

The frozen matrix `B` is singular:

```text
rank_Q(B) = 111,  nullity_Q(B) = 9.
```

This is not inferred from a floating tolerance.  Nine independent exact
rational kernel vectors certify the upper bound, while row reduction modulo
1009 gives rank 111 and certifies the matching lower bound.

For the full edge-restricted Hessian:

```text
rank_F101(H)  = 708,
rank_F1009(H) = 720,
hence rank_Q(H) = 720.
```

The unsuccessful first prime is reported to avoid hiding prime selection.
Full rank modulo one prime is an exact nonzero-determinant certificate over
the integers and hence over the rationals.

For any real symmetric perturbation `E`, in an eigenbasis of `B`,

```text
d^2 Tr(B^4)[E,E]
 = sum_i 12 lambda_i^2 E_ii^2
   + sum_{i<j} 8(lambda_i^2 + lambda_i lambda_j + lambda_j^2) E_ij^2.
```

Every coefficient is nonnegative because

```text
x^2 + xy + y^2 = (x+y/2)^2 + 3y^2/4 >= 0.
```

Thus the Hessian is positive semidefinite on the full space of symmetric
matrices.  Its exact rank 720 on the linearly independent edge directions
makes the restricted Hessian positive definite.  As a numerical control only,

```text
lambda_min(H) = 434.491379676...
lambda_max(H) = 51366.5394399...
```

No numerical tolerance is needed for the rank or inertia conclusion.

## Why the old matrix failed

The legacy script constructed

```text
G_ef = 12 sum_k lambda_k^2 S_ke S_kf,
S_ke = <psi_k, E_e psi_k>.
```

For a simple eigenvalue this is only the contribution
`12 lambda_k^2 (lambda'_k)^2` to the second derivative of
`sum_k lambda_k^4`.  It omits the `4 lambda_k^3 lambda''_k` contribution,
which contains eigenvector mixing.  At degeneracies, its individual diagonal
sensitivities additionally depend on the chosen orthonormal eigenbasis.

Three quantitative controls expose the problem:

```text
legacy numerical inertia                 = 101 + 619 + 0
||H-G||_F / ||H||_F                      = 0.949678740588...
||G_rotated-G||_F / ||G||_F              = 0.416082162046...
```

The rotation is block-orthogonal inside the thirteen degenerate eigenspaces
of the same `B`; it changes no operator or geometry.  Therefore `G` is not a
canonical function of `B` and the edge perturbations.

There is also a decisive search-space objection.  `G` factors through only
120 eigenvalue-sensitivity rows.  Because nine eigenvalues of `B` vanish,

```text
rank(G) <= 111, hence nullity(G) >= 609
```

by construction.  The old calculation could never have falsified a huge
zero-mode sector.  Its 619 zero modes therefore were not evidence for gauge
freedom or spectral degeneracy; the search space had already forced at least
609 of them.

## Physical interpretation

**DERIVED.**  `Tr(B(w)^4)` has a non-diagonal, positive-definite quadratic
stiffness at the frozen background in the chosen 720 edge-weight variables.

**DERIVED NEGATIVE.**  The statements “101 positive modes”, “619 zero modes”,
and “the 619 modes contain the 119 gauge directions” do not describe the full
Hessian and are withdrawn.

**DERIVED NEGATIVE.**  This calculation does not establish a propagating
graviton.  A Hessian of a static finite action is a stiffness matrix, not a
propagator by itself.  Propagation requires a derived time evolution or
Lorentzian kinetic operator; a graviton additionally requires a field
dictionary, gauge constraints and quotient, spin-2 transformation law,
continuum limit, and coupling to a conserved stress tensor.  Those ingredients
are not present in this calculation.

**STRUCTURAL.**  Separate cochain calculations may still exhibit coexact
mode multiplicities analogous to two continuum transverse form families.
Without the missing field dictionary, calling those modes the two physical
spin-2 graviton polarizations is not derived.

**OPEN.**  Whether another theory-derived action and constraint structure
turn any part of the edge/cochain carrier into gravity remains open.  This
correction removes evidence; it neither proves nor disproves such a future
construction.

**OPEN.**  No Newton constant, Planck length, Planck mass, or absolute physical
scale follows from this Hessian.  Positive stiffness fixes a dimensionless
shape only after the normalization of the action and the physical meaning of
the edge variables are supplied.

## Consequent documentation corrections

The Hessian error exposed a broader category boundary in the same gravity
section.  The active summaries and paper variants now state the following.

- **DERIVED:** `Im(d_0)` and `Im(d_1^T)` have dimensions 119 and 601 and form
  exact/coexact Hodge sectors.  Calling them gauge and physical sectors needs
  a derived metric-field and gauge-action dictionary.
- **DERIVED:** the chosen uniform-neighbor Ollivier curvature is zero on all
  720 graph edges.  This graph/measure statement is not continuum Ricci
  flatness and does not by itself define an Einstein vacuum.
- **STRUCTURAL:** `C h_T = 8 pi G T_T` is a candidate response ansatz.  Only
  `C=d_1^T d_1` is constructed; `h_T`, `T_T`, `G`, and the equation are not
  selected by the finite data.
- **STRUCTURAL:** `c_lat^2=5` is a dimensionless lattice coefficient under a
  chosen update reading.  Physical length/tick units and identification with
  light are additional dictionaries, so the SI speed of light is not derived.
- **STRUCTURAL/OPEN:** the product
  `lambda_1^(exact) 4 phi^2=24=|2T|` is exact arithmetic, but its former
  “gap--Planck” and graviton-mass labels have no derived scale dictionary.

These corrections do not prove that gravity cannot emerge.  They make clear
that any such emergence must be demonstrated in collective/refinement
dynamics with a spin-2 field map and universal source coupling, rather than
declared from finite Hodge dimensions or a static Hessian.

## Status ledger

| Claim | Status | Evidence |
|---|---|---|
| 600-cell edge split is `120+600` | **DERIVED** | exact integer reconstruction |
| `rank_Q(B)=111`, `nullity_Q(B)=9` | **DERIVED** | exact kernel plus modular lower bound |
| frozen full-Hessian formula | **DERIVED** | noncommutative polynomial derivative |
| `rank_Q(H)=720` | **DERIVED** | rank 720 modulo 1009 |
| restricted inertia is `(720,0,0)` | **DERIVED** | analytic PSD plus exact full rank |
| legacy `101+619` is reproducible | **STRUCTURAL** | floating eigendecomposition control |
| legacy matrix equals the Hessian | **REFUTED** | 0.95 relative discrepancy |
| legacy matrix is canonical | **REFUTED** | 0.42 change under legal basis rotation |
| Hessian supplies a graviton propagator | **REFUTED AS CLAIMED** | static stiffness lacks required dynamics/dictionary |
| theory derives gravity by another route | **OPEN** | not decided by this audit |
| theory selects a Planck scale here | **OPEN / unsupported** | no dimensionful normalization or `G` |

## Reproduction

Run only the targeted verifier:

```bash
/home/razvan/science/.venv/bin/python \
  reproducible/verify_gravity_box4_full_hessian.py
```

Expected result at this commit: `14/14 checks passed`.  The full repository
suite is deliberately not run for this mission, following the user's current
instruction.
