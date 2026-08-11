# B1 chamber embedding: exact non-selection certificate

Date: 2026-08-11  
Preregistration commit: `2c16049`

## Question and scope

The fixed 120-chamber carrier admits the exact noncommutative B1 witness

\[
A=M_2(\mathbb C)\oplus\mathbb C^3.
\]

The witness already satisfies the complete stated hypothesis list:
faithful unital star representation, order zero, first order,
`[gamma,A]=0`, nonzero inner one-forms, metric-dimension-zero
orientability, nondegenerate intersection form and connectedness.  The
remaining question here is narrower: does the chamber geometry **select**
the displayed embedding, or was it one result of a large search?

The algebra type, oriented cells and capacities remain STRUCTURAL INPUTS:

```text
cells      = ((0,1)x2, (1,2)x25, (3,1)x12, (2,3)x19)
capacities = (4,25,12,19).
```

No Standard-Model or phenomenological target is used in this test.

## Change from the preregistered protocol

The preregistration requested a complete orbit census.  The exhaustive
CP-SAT run was deliberately stopped before completion because a finite
countercertificate had already made the exact total irrelevant to the
selection claim.  Therefore:

- **INCOMPLETE:** the complete number of admissible colourings and orbits is
  not known and is not reported;
- no partial count or hit fraction is promoted to a theorem;
- the accepted result below uses only two explicit certificates and exact
  finite checks.

This does not satisfy the preregistered census acceptance boundary.  It does
settle the weaker but physically decisive yes/no question of uniqueness:
two inequivalent valid orbits are enough to refute it.

## Exact geometry

The verifier reconstructs the icosahedron over
\(\mathbb Q(\sqrt5)\), its 120 complete flags, chamber adjacency \(D\),
orientation \(\gamma\), central reflection \(J\), and

\[
S=(DJ)|_{H_+}.
\]

It derives that `S` is a connected, loopless, invertible, 3-regular graph
with 60 vertices and 90 edges.  Its complete graph automorphism group has
order 60 and is exactly the geometric rotation group \(A_5\).

For a colouring \(c\) and automorphism \(g\), the orbit action is

\[
(g\cdot c)(g x)=c(x).
\]

The committed witness and a second explicit colouring both have the frozen
capacities and only the allowed edge-colour pairs

\[
\{01,12,13,23\}.
\]

Each has a free orbit of size 60.  The two complete `Aut(S)=A5` orbits are
disjoint.

## Full gate reconstruction

For **each** of the two colourings, the verifier rebuilds the full
120-dimensional representation rather than inferring validity from the
support alone.  On the seven complex basis elements of
\(M_2(\mathbb C)\oplus\mathbb C^3\), it checks:

1. faithful, unital, star-closed and noncommutative representation;
2. commutation with \(\gamma\);
3. order zero and full first order;
4. nonzero inner one-forms;
5. the exact orientation cycle;
6. the antisymmetric intersection matrix

   \[
   \cap=
   \begin{pmatrix}
   0&2&0&0\\
   -2&0&25&-12\\
   0&-25&0&19\\
   0&12&-19&0
   \end{pmatrix},
   \qquad \operatorname{rank}\cap=4,
   \qquad \det\cap=1444;
   \]

7. connectedness: the complex-linear commutator map has rank 6 on the
   seven-dimensional algebra, leaving exactly the scalar kernel.

Both certificates pass every gate exactly.

## Independent continuous ambiguity

A central support colouring determines only the four-dimensional carrier of
the noncommutative cell.  It does not choose a particular copy

\[
M_2(\mathbb C)\otimes I_2\subset M_4(\mathbb C).
\]

Its unitary conjugacy orbit has real dimension

\[
\dim U(4)-\dim\frac{U(2)\times U(2)}{U(1)}
=16-(4+4-1)=9.
\]

Thus even a unique central support would not have selected the internal
matrix factor.

## Status ledger

- **DERIVED:** `Aut(S)` is precisely geometric \(A_5\), of order 60.
- **DERIVED:** at least two disjoint free \(A_5\)-orbits of supports satisfy
  every B1 gate.
- **DERIVED NEGATIVE:** the fixed geometry does not select the B1 embedding.
- **DERIVED:** every such central support retains a 9-real-dimensional
  internal \(M_2\) embedding ambiguity.
- **INCOMPLETE:** the complete support-orbit census and its exact total.
- **STRUCTURAL:** the algebra type, cells, multiplicities and capacities.
- **DERIVED:** B1 itself remains refuted as a mathematical commutativity
  claim by the exact noncommutative witness.
- **DERIVED NEGATIVE (physics):** that existence theorem does not open a
  Standard-Model gate; the noncommutative algebra was not selected by the
  geometry.

## Reproduction

Run only the targeted verifier:

```bash
/home/razvan/science/.venv/bin/python \
  reproducible/verify_chamber_b1_embedding_orbits.py
```

Expected result: `17/17` with `COMPLETE_CENSUS=False` and
`INEQUIVALENT_ALL_GATE_ORBITS>=2`.

