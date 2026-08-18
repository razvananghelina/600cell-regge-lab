# Preregistered protocol: canonicity of the first metric phase space

Date: 2026-08-12

Status at registration: **PROTOCOL ONLY -- NO HAMILTONIAN SELECTED**

## 1. Question and scope

The preceding inventory found that the current repository has no Hamiltonian
constraint acting on the 600-cell metric.  This protocol asks the next,
strictly prior question:

> Does the fixed 600-cell geometry and its full `H4` symmetry select a
> quadratic kinetic Hamiltonian on the natural Regge metric phase space, or
> does it leave a multi-parameter family?

The configuration arena is the open set `M_Regge` of admissible positive
squared lengths on the 720 edges, near the equilateral point.  The phase-space
arena is its cotangent bundle

```text
T* M_Regge = {(x_e,p_e): e an edge},
omega = sum_e d p_e wedge d x_e.
```

This cotangent symplectic form is canonical.  It does not by itself select a
Hamiltonian, lapse carrier, constraint or time normalization.

## 2. Frozen operator class

At the equilateral point consider real quadratic momentum terms

```text
T_K(p) = (1/2) p^T K p,
K=K^T,
```

invariant under the full `H4` permutation action on the 720 unoriented edges.
Thus `K` must lie in the symmetric part of the permutation commutant.

This is a necessary local linearization of any smooth invariant kinetic
Hamiltonian.  A failure of uniqueness here refutes selection by `H4` alone;
it is not a no-go against a future extra geometric principle.

The earlier finite-Hessian audit already reported an order-20 stabilizer and
62 stabilizer orbits on edges.  Before this registration an exploratory
uncommitted calculation suggested that transpose pairs those 62 orbitals into
47 symmetric parameters.  That number is therefore **not blind evidence**;
the registered verifier must rebuild it and expose the complete mechanical
certificate.  No physical target or desired kinetic matrix has been compared.

## 3. Frozen exact/combinatorial checks

The verifier must independently rebuild the 600-cell and check:

1. f-vector `(120,720,1200,600)`;
2. the full quaternionic `H4` action is transitive on unoriented edges;
3. the stabilizer of one edge has order 20;
4. its edge orbits form the orbital basis of the full equivariant
   endomorphism algebra;
5. transpose of an orbital is computed using an explicit group transporter,
   not guessed from enumeration order;
6. the number of transpose-fixed orbitals and paired orbitals, hence the exact
   dimension of the symmetric commutant;
7. the number of symmetric parameters supported at each line-graph distance;
8. the parameter count after the preregistered nearest-neighbour condition
   `K_ef=0` whenever the two edges neither coincide nor share a vertex;
9. the parameter count when support is restricted further to edge pairs
   contained in a common tetrahedron;
10. positivity does not reduce a linear family to a ray: because the identity
    is positive definite, every symmetric invariant support space containing
    it has a relatively open positive cone around the identity.

All orbit partitions and transpose pairings must be written to a
machine-readable result.  Floating coordinates may identify the exact finite
permutations only if bijectivity, incidence preservation, group-order/stabilizer
counts and orbit partitions all pass exactly after identification.

## 4. Decision boundary

- **DERIVED CANONICAL KINETIC RAY:** the symmetric invariant space, after a
  locality condition already fixed above and itself supplied by the incidence
  geometry, is one-dimensional up to overall time/energy scale.
- **DERIVED CANONICITY OBSTRUCTION:** the full invariant space has dimension
  greater than one and even the nearest-neighbour or common-tetrahedron class
  contains more than one positive-definite ray.  Then `H4`, incidence locality
  and positivity do not select a kinetic Hamiltonian.
- **OPEN/INCOMPLETE:** the finite action, transporter pairing or locality
  counts cannot be certified.

The diagonal ultralocal choice `K proportional to I` is a canonical-looking
member of the family, but it counts as selected only if ultralocality is
derived independently.  Declaring all off-diagonal couplings zero because
that yields uniqueness would be an ansatz, not a result.

## 5. Claims explicitly excluded

This audit does not claim that:

- no canonical Hamiltonian can arise after adding a time/slab geometry;
- the Regge/DeWitt supermetric is invalid;
- any invariant quadratic form is the physical one;
- first-class constraints or Lorentzian signature have been derived;
- a Planck scale or the value of `c` follows.

Only the targeted verifier and a static registry check may run.  No full suite
and no PDF build.
