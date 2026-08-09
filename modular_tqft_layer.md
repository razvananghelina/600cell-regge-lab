# Modular/TQFT foundational-layer audit

Date: 2026-07-27 (sixth session).

This audit obeys RULE ZERO and the binding post-audit vocabulary:
**DERIVED** means an exact consequence of stated objects; **STRUCTURAL** means
a precise compatible identification requiring an additional dictionary;
**PATTERN** means an observed concordance that does not constrain the physical
formula; **OPEN** means that a required construction is absent.  Modular data
below are computed by `reproducible/verify_modular_tqft_layer.py`.

## 1. What the bootstrap actually says

The repository has used two different definitions.

1. The original paper/requested bootstrap is
   \[
   d_{1/2}(n)=2\cos(\pi/n)=(1+\sqrt n)/2,\qquad n=k+2.
   \]
   Here \(d_{1/2}\), not \(d_1\), is the fundamental \(SU(2)_{n-2}\)
   quantum dimension.
2. `s02_bootstrap_closure.md` later defines
   \[
   d_1(n)={\sin(3\pi/n)\over\sin(\pi/n)}
          =1+2\cos(2\pi/n)
   \]
   and matches it to the already-fixed constant \(\phi\).  This is a
   different conditional problem.  At \(n=5\), both \(d_{1/2}\) and \(d_1\)
   happen to equal \(\phi\), but the definitions must not be conflated.

For the original equation, if \(n\ge9\), its right side is at least \(2\)
while \(2\cos(\pi/n)<2\), so there are no solutions.  Exact evaluation for
\(3\le n\le8\) gives only \(n=5\).  At that point the equation is precisely
\[
2\cos(\pi/5)=(1+\sqrt5)/2.
\]
Thus uniqueness is **DERIVED**, conditional on choosing this matching
equation.  The choice of matching equation is **STRUCTURAL**, not a physical
vacuum principle.  Its content is the classical pentagon identity, with the
integer five appearing both as root-of-unity order and in the surd.  Calling
it dynamical selection is **PATTERN/overstatement**.  No nature-selection or
vacuum functional follows.

## 2. Exact modular data

We use \(T=e^{-2\pi i c/24}\operatorname{diag}(\theta_a)\); listing twists
separately removes the conventional framing phase.

### Fibonacci category \(\mathrm{Fib}=(G_2)_1\)

Objects are \(1,\tau\), with
\[
\tau^2=1+\tau,\quad d=(1,\phi),\quad
\mathcal D=\sqrt{\phi+2},
\]
\[
S={1\over\sqrt{\phi+2}}
\begin{pmatrix}1&\phi\\ \phi&-1\end{pmatrix},\qquad
\theta=(1,e^{4\pi i/5}),\qquad c=14/5\pmod8.
\]
The displayed \(S\) squares to the identity and Verlinde diagonalization
gives \(N_\tau=\left(\begin{smallmatrix}0&1\\1&1\end{smallmatrix}\right)\).

### \((F_4)_1\)

Objects are \(1,x\), with the same fusion rule and unitary dimensions as Fib:
\[
x^2=1+x,\quad d=(1,\phi),\quad\mathcal D=\sqrt{\phi+2},\quad S=S_{\rm Fib}.
\]
The distinction is ribbon data:
\[
\theta=(1,e^{6\pi i/5}),\qquad c=26/5\pmod8.
\]

### \((E_8)_1\)

There is only the vacuum:
\[
1^2=1,\quad d=\mathcal D=1,\quad S=(1),\quad\theta=(1),\quad c=8=0\pmod8.
\]
Its character \(T\)-phase is \(e^{-2\pi i/3}\), despite its trivial twist.

### \(SU(2)_3\)

Label objects by \(a=0,1,2,3\) (spin \(a/2\)).  Fusion is truncated
Clebsch--Gordan:
\[
a\otimes b=\!\!\sum_{\substack{c=|a-b|\\c\equiv a+b\ (2)}}^{
\min(a+b,6-a-b)}c.
\]
In particular \(1^2=0+2,\ 1\,2=1+3,\ 1\,3=2,\ 2^2=0+2,\ 2\,3=1,\
3^2=0\).  Moreover
\[
d=(1,\phi,\phi,1),\quad \mathcal D=\sqrt{5+\sqrt5},
\]
\[
S_{ab}=\sqrt{2/5}\sin((a+1)(b+1)\pi/5),
\]
\[
\theta=(1,e^{3\pi i/10},e^{4\pi i/5},-i),\qquad c=9/5\pmod8.
\]
The full Verlinde ring has rank four and is not \(\mathbb Z[\phi]\).  Its
even subcategory \(\{0,2\}\) is Fib and has based ring
\(\mathbb Z[t]/(t^2-t-1)\cong\mathbb Z[\phi]\).

All statements in this section are **DERIVED**.

## 3. The exceptional conformal-embedding bridge

The central charges close exactly:
\[
c(G_2{}_1)+c(F_4{}_1)=14/5+26/5=8=c(E_8{}_1).
\]
The \(E_8{}_1\) vacuum branches as
\[
\chi^{E_8}_0=\chi^{G_2}_1\chi^{F_4}_1+
              \chi^{G_2}_\tau\chi^{F_4}_x.
\]
This is certified directly by modular data.  In product basis
\((11,1x,\tau1,\tau x)\), the branching vector \(b=(1,0,0,1)^T\) obeys
\((S_{\rm Fib}\otimes S_{\rm Fib})b=b\).  Also
\(h_\tau+h_x=2/5+3/5=1\), so both summands have the same \(T\)-phase
\(e^{-2\pi i/3}\) as the \(E_8{}_1\) vacuum.  Equivalently the extension
algebra is \(1\boxtimes1\oplus\tau\boxtimes x\).

The embedding and compatibility are **DERIVED** mathematics.  Using this
extension as the bridge between the paper's Fibonacci and McKay/\(E_8\)
pillars is **STRUCTURAL**: the repository supplies no functor from this
chiral conformal extension to its 600-cell operators or mass assignments.

## 4. Galois conjugation, \(F_4\), and entropy

The non-unitary Galois conjugate of Fib is the Yang--Lee category:
\[
\tau^2=1+\tau,\quad d_\tau=-1/\phi,\quad
\mathcal D_{\rm YL}=\sqrt{3-\phi},
\]
\[
S_{\rm YL}={1\over\sqrt{3-\phi}}
\begin{pmatrix}1&-1/\phi\\-1/\phi&-1\end{pmatrix},\quad
\theta=(1,e^{-2\pi i/5}),\quad c=-22/5=18/5\pmod8.
\]
It has the same based fusion ring as Fib and \((F_4)_1\), but it is neither
unitary Fib nor \((F_4)_1\).  In particular the nontrivial \(F_4{}_1\) twist
is \(e^{6\pi i/5}=e^{-4\pi i/5}\), not the Yang--Lee
\(e^{-2\pi i/5}\), and its positive dimension is \(\phi\), not
\(-1/\phi\).  Therefore the dark sector is not the complementary \(F_4\)
factor.  That proposed identification is **DERIVED negative**.

Formally,
\[
\log{\mathcal D_{\rm Fib}\over\mathcal D_{\rm YL}}
={1\over2}\log{\phi+2\over3-\phi}=\log\phi.
\]
For the four-object \(SU(2)_3\) Galois pair the squared dimensions are twice
these values, \(5\pm\sqrt5\), so the same ratio results.  The paper's
\(\ln\phi\) arithmetic is therefore **DERIVED, right for the wrong reason**:
it follows from Galois-conjugate global dimensions.  Calling
\(\log\mathcal D_{\rm YL}\) a physical topological entanglement entropy is
**OPEN/PATTERN**, because Yang--Lee is non-unitary and the paper constructs
no dark-sector ground state, spatial cut, or non-unitary entropy
prescription.

## 5. Load-bearing tests

### Radiative claims

For \(SU(2)_3\),
\[
N_{1/2}=\begin{pmatrix}0&1&0&0\\1&0&1&0\\0&1&0&1\\0&0&1&0\end{pmatrix},
\qquad\|N_{1/2}\|_F^2=6.
\]
For Fib, \(\|N_\tau\|_F^2=3\).  Thus “\(=6\)” is **DERIVED** only in the
four-object \(SU(2)_3\) convention.  Its equality to the separately named
repository integer \(b_1=6\) is a **PATTERN** unless a physical vertex map
is constructed.

The prefactor \(\sqrt{2/5}\) is visibly present in \(S^{SU(2)_3}\), but an
overall normalization is basis/convention data constrained by unitarity.
In the purported self-energy
\(S_{jl}S_{1l}/S_{0l}\), nonconstant sine/eigenvalue factors remain.
No propagator, label assignment, contraction, or map to \(\delta_1\) is
defined.  Hence neither \(c_{\rm bare}=\sqrt{2/5}\) nor
\[
c_{\rm eff}=\sqrt{2/5}(1-6\alpha^2)
\]
is derived by modular data.  Both are **PATTERN** mass corrections.  Ward
identities from QED are not consequences of this modular category.

### Three generations versus level three

\(k=3\) counts the truncation level of \(SU(2)_3\); \(N_{\rm gen}=3\) counts
putative fermion generations.  Modular data provide no representation,
index, anomaly, branching multiplicity, or functor identifying them.
Their equality is **PATTERN**, not a derivation.  A real derivation would
need a canonical functor from category/simple-object or conformal-embedding
data to three chiral copies of a constructed SM matter module, with the
generation multiplicity proved to equal an appropriate categorical index.
That remains **OPEN**.

### The two appearances of \(\mathbb Z[\phi]\)

This is the genuine positive result.  The Fib based fusion ring has basis
\((1,\tau)\), multiplication matrix \(N_\tau\), and positive character
\(\tau\mapsto\phi\).  The repository's stationary Fibonacci tower uses
\(F=\left(\begin{smallmatrix}1&1\\1&0\end{smallmatrix}\right)\), which is
\(N_\tau\) after swapping the basis.  Because \(\det F=-1\), its ordered
\(K_0\) is canonically
\[
(\mathbb Z[\phi],\ \mathbb Z[\phi]\cap\mathbb R_{\ge0},\ 1),
\]
with the same Perron--Frobenius positive character.  Thus these are the same
ordered based ring up to the explicit vacuum/nontrivial basis swap:
**DERIVED**, not two unrelated sightings.  Embedding the paper's selected
mass labels into that ring is algebraic **DERIVED** as already audited;
interpreting order or tower level as physical mass/energy is
**STRUCTURAL**.

## 6. Verdict

- **DERIVED:** all displayed modular data; bootstrap uniqueness conditional
  on the chosen equation; the even \(SU(2)_3\) Fib subring; the
  \(G_2{}_1\times F_4{}_1\subset E_8{}_1\) modular branching certificate;
  the Yang--Lee conjugate; the formal \(\ln\phi\) dimension ratio; and the
  canonical ordered-ring identification with the Fibonacci tower.
- **STRUCTURAL:** choosing the bootstrap matching; using the exceptional
  conformal embedding as a bridge to the 600-cell/McKay operators; and
  reading ordered \(K_0\) levels as physical scales.
- **PATTERN/decoration:** TQFT derivation of the radiative mass coefficient,
  the equality \(6=b_1\) as physical vertex phase space, \(N_{\rm gen}=k\),
  “dark sector = Galois category,” and physical entropy language.
- **OPEN:** a functor connecting the conformal embedding to the repository's
  discrete operators/matter module; a physical non-unitary entropy
  prescription; and any categorical derivation of radiative corrections or
  generations.

The modular/TQFT layer contains rigid and worthwhile mathematics, especially
the exceptional conformal embedding and the ordered-ring identity.  It is
not presently load-bearing for the paper's masses, generations, or dark
sector.
