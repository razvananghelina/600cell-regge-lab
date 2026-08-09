# Galois doubling of the McKay node space: exact scope and obstruction

## Result

The proposed doubling removes the **dimension/sign obstruction** to a KO-dimension-six real structure, but it does **not** yet construct the finite Standard-Model triple.  More sharply, its proposed linchpin separates into two different operations:

1. **DERIVED:** the nontrivial outer automorphism of (A_5) and the field automorphism
   \(\sigma:\sqrt5\mapsto-\sqrt5\) have the same action on the character table;
2. **DERIVED (negative):** complex conjugation of the derived color embedding does not realize that action.  The embedding (3':A_5\to SO(3)\subset SU(3)) is real, so both the defining (3) and its complex conjugate \(\bar3\) restrict to the same (A_5)-module (3').  The outer twist is the inequivalent module (3).

Thus Galois doubling is a viable structural source of a second sheet and KO signs, but it does not by itself select color conjugation, the finite algebra action, or hypercharge.

All finite statements below are checked by `reproducible/verify_galois_doubling.py` using exact arithmetic in \(\mathbb Q(\sqrt5)\).

## 1. The linchpin character check

Use the (A_5) class order

\[
(1A,2A,3A,5A,5B).
\]

The outer automorphism induced by an odd permutation in (S_5) exchanges (5A\leftrightarrow5B).  With the repository convention

\[
\chi_3=(3,-1,0,\phi,\phi'),\qquad
\chi_{3'}=(3,-1,0,\phi',\phi),
\]

this gives

\[
1\mapsto1,quad3\leftrightarrow3',\quad4\mapsto4,quad5\mapsto5,
\]

exactly as does \(\sigma(\phi)=\phi'\).  **DERIVED.**

The outer automorphism lifts to (2I\cong SL(2,5)): conjugation by an element of (GL(2,5)) of nonsquare determinant preserves (SL(2,5)) and induces the non-inner automorphism of (PSL(2,5)\cong A_5).  On the repository's irreps its character action is

\[
\rho_2\leftrightarrow\rho_3,qquad
\rho_4\leftrightarrow\rho_5,qquad
\rho_1,\rho_6,\rho_7,\rho_8,\rho_9\ \hbox{fixed}.
\]

Equivalently,

\[
2\leftrightarrow2',\quad3\leftrightarrow3',\quad
1,4,4_s,5,6\ \hbox{fixed}.
\]

On the (2I) classes it exchanges (10A\leftrightarrow10B) and (5A\leftrightarrow5B), fixing the other five classes.  Direct comparison of every exact character proves that pullback by this outer automorphism equals coefficientwise Galois conjugation.  **DERIVED.**

### Why this is not (3\leftrightarrow\bar3)

For the already-derived embedding \(\iota_{3'}:A_5\to SO(3)\to SU(3)\), all matrices are real.  Therefore

\[
\overline{\iota_{3'}(g)}=\iota_{3'}(g),
\qquad
\bar3\!\downarrow_{A_5}\cong3\!\downarrow_{A_5}\cong3'.
\]

By contrast, \(\iota_{3'}\circ\alpha\) has character (3), not (3').  The two embeddings cannot be conjugate in (SU(3)), since conjugate representations have equal characters and their values already differ on (5A).  **DERIVED (negative).**

Consequently the finite subgroup forgets the (3)-versus-\(\bar3\) orientation.  Interpreting the outer/Galois sheet as charge conjugation is an additional choice, not a consequence of the derived real embedding.  **OPEN.**

## 2. The doubled sign table

Let

\[
W=\bigoplus_{i=1}^9\rho_i,qquad \dim W=30,
\]

and let (W^\sigma) carry the outer-twisted action.  The permutation above preserves integer versus half-integer spin, hence commutes with the McKay bipartite grading (g).  Its weighted split remains (16+14).  **DERIVED.**

Set

\[
H=W\oplus W^\sigma,qquad
\Gamma=\begin{pmatrix}g&0\\0&-g\end{pmatrix}.
\]

Choose anti-linear isometric intertwiners between each irrep and its twisted partner and denote their direct sum by (S).  Then

\[
J_+(v,w)=(S^{-1}w,Sv)
\]

satisfies

\[
J_+^2=+1,qquad J_+\Gamma=-\Gamma J_+.
\]

For a real McKay adjacency (D) and its correctly relabelled twist (D^\sigma),

\[
D_H=\operatorname{diag}(D,D^\sigma)
\]

gives (J_+D_H=D_HJ_+), and (D_H\Gamma=-\Gamma D_H).  These are precisely the finite KO-dimension-six signs

\[
(J^2,JD/DJ,J\Gamma/\Gamma J)=(+,+,-).
\]

The (+/-) eigenspaces of \(\Gamma\) both have dimension (30), so the old (16\ne14) same-space obstruction is absent.  **DERIVED, conditional on the doubled space and chosen anti-linear isometries.**

Other exact variants are also possible: putting (+g) on both sheets changes the last sign to (+), while changing the sign on the second Dirac block changes (JD=DJ) to (JD=-DJ); a signed sheet swap gives (J^2=-1).  **DERIVED.**

Two qualifications prevent promotion to a derived spectral triple:

- Arithmetic \(\sigma\) is not itself an antiunitary operation for the chosen real embedding of \(\mathbb Q(\sqrt5)\); it changes positive real magnitudes.  Antiunitarity requires the separately chosen isometric intertwiners (S).  **STRUCTURAL choice / OPEN derivation.**
- The second sheet (W^\sigma), although natural and representation-theoretically exact, is not selected by an existing finite-action theorem in the repository.  **STRUCTURAL candidate.**

## 3. Algebra action and first-order gate

No derived representation

\[
\rho:\mathbb C\oplus\mathbb H\oplus M_3(\mathbb C)\longrightarrow\operatorname{End}(W)
\]

currently exists.  The statements “color on integer-spin nodes” and “weak action on spinorial nodes” do not define a unital associative-algebra representation: a (2I) irrep is not thereby an (M_3(\mathbb C))- or \(\mathbb H\)-module, nor is a block assignment forced.  **OPEN.**

Therefore the target order-zero equation

\[
[\rho(a),J\rho(b)^*J^{-1}]=0
\]

cannot honestly be evaluated, and neither can its first-order equation or the surviving commutant (u(1)).  Defining unrelated left and right actions on separate sheets would make order zero tautological but would insert the missing bimodule by hand.  It is not counted as a solution.  **DERIVED scope limit.**

There is one useful route-specific screen.  The canonical vertex-scalar algebra \(\mathbb C^9\) acts diagonally and satisfies order zero, but independent endpoint projectors on any McKay edge give a nonzero

\[
[[D,a],b^{\mathrm o}].
\]

The verifier exhibits this exactly on the edge \(\rho_1-\rho_2\); global constants survive.  Thus McKay adjacency plus the full independent node algebra is not a finite first-order triple.  **DERIVED (scoped negative).**  This does not rule out a smaller, non-diagonal, or genuinely derived Standard-Model bimodule action.

Because the construction stops at this gate, there is no derived (Y) to feed into the Route C anomaly equations.  **OPEN.**

## 4. Secondary doubled-quiver Hom test

For the eight McKay edges, both orientations give

\[
2\sum_{i-j}\dim\rho_i\dim\rho_j=240,
\]

equal to the number of (E_8) roots.  **DERIVED equality; PATTERN as a physical identification.**

Under the diagonal (2I) action,

\[
\bigoplus_{i\leftrightarrow j}\operatorname{Hom}(V_i,V_j)
\cong
16\rho_2\oplus6\rho_3\oplus16\rho_7\oplus22\rho_9.
\]

The dimension is (32+12+64+132=240).  In particular the space has no integer-spin (1,3,3',4,5) summands.  **DERIVED.**  Hence its direct diagonal decomposition does not expose the color and singlet blocks required by (M_{15}) or (M_{16}).  This is a route-specific negative, not a proof against other quiver functors.  The weak tensor-by-(2) incidence is canonical, but extracting chiral matter would require an additional projection or module functor.  **OPEN.**

## Classification ledger

- **DERIVED:** outer twisting equals \(\mathbb Q(\sqrt5)\) Galois conjugation on all (A_5) and (2I) characters.
- **DERIVED (negative):** it is not ordinary complex conjugation of the real (3') color embedding and does not distinguish (3) from \(\bar3\).
- **DERIVED, conditional:** doubling realizes the KO6 sign table and removes the (16/14) dimension obstruction.
- **STRUCTURAL:** (W\oplus W^\sigma) and chosen anti-linear isometric intertwiners are coherent candidate data, not yet forced data.
- **DERIVED (scoped negative):** full vertex scalars with McKay adjacency fail first order.
- **DERIVED:** the bidirected Hom space has dimension (240) and the exact decomposition displayed above.
- **PATTERN:** identifying that (240) with the (E_8) root count or the Galois sheet with physical charge conjugation.
- **OPEN:** a derived unital \(\mathbb C\oplus\mathbb H\oplus M_3(\mathbb C)\) bimodule action, an antiunitary selected by the discrete geometry, color orientation, a first-order Dirac operator, and the resulting anomaly-free generation-blind (Y).

## Sharp remaining target

A successful continuation must construct, rather than assign, a bimodule representation on (H) for which order zero holds, solve first order for a derived (D), and produce a commutant generator constant on weak doublets and blind to generation.  Only then does the already-derived anomaly factorization decide whether its charges are the Standard-Model tuple.
