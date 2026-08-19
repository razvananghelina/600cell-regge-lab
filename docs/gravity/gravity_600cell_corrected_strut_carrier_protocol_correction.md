# Pre-execution correction: the uniform strut line need not be singular

Date: 2026-08-19

The frozen corrected-carrier protocol asked for “the singular value carried
by the normalized uniform vertex vector” and for 119 complement singular
values.  That wording assumed without proof that the uniform line is an
invariant right-singular subspace.

The schedule stabilizer has order 24 and acts freely on 120 vertices.  Its
vertex representation therefore contains five copies of the trivial
representation.  Equivariance permits mixing inside that five-dimensional
isotypic component; it does not force the particular globally uniform vector
to diagonalize `G^T G`.

This correction is committed before the verifier exists and before any
accepted-slab carrier spectrum is evaluated.  It changes no carrier formula,
incidence condition, target firewall or mechanical outcome.  The verifier
will instead record:

```text
||G u||,
spectrum(Q^T G^T G Q),
||Q^T G^T G u||,
```

where `u` is normalized uniform and `Q` is a deterministic Householder basis
of its Euclidean orthogonal complement.  The last quantity explicitly tests,
rather than assumes, invariance.

No Hessian, Schur or tangent target has been loaded.  This is a correction of
an invalid representation-theoretic inference, not a response to data.

