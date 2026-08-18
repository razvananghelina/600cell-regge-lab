# S06 Hopf Closure

Scop: sa inchidem `S06` exact pe cele patru propoziții stabilite dinainte.

## Formal Setting

Nu vorbim despre fibrări Hopf continue ale lui `S^3` în general.

Obiectul exact este:

- o `discrete Hopf fibration` a 600-cell-ului,
- adică o partiție a celor 120 vârfuri în 12 decagoane,
- obținută din clase laterale ale unui subgrup `C_{10} \le 2I`.

Pentru o astfel de fibrară `F`, definim:

\[
A_{\mathrm{fiber}}(F),
\qquad
A_{\mathrm{cross}}(F)=A-A_{\mathrm{fiber}}(F),
\]

\[
L_{\mathrm{fiber}}(F)=2I-A_{\mathrm{fiber}}(F),
\qquad
L_{\mathrm{cross}}(F)=10I-A_{\mathrm{cross}}(F),
\]

\[
\Box_F(c)=cA_{\mathrm{fiber}}(F)-A.
\]

## Verified Propositions

Scriptul focalizat

- [reproducible/verify_hopf_fibration_invariants.py](D:\infinity\ToE\science\reproducible\verify_hopf_fibration_invariants.py)

verifică exact următoarele patru propoziții și nimic mai mult.

### P1. Count

Există exact

\[
6
\]

fibrări discrete distincte de acest tip.

### P2. Unique nontrivial coefficient

Pentru fiecare fibrară `F`, singurul coeficient netrivial pentru care

\[
\ker(\Box_F(c)) \neq 0
\]

este

\[
c=6.
\]

### P3. Stable kernel sector

Pentru fiecare fibrară `F`,

\[
\dim \ker(\Box_F(6)) = 9,
\]

iar kernelul coincide cu același sector spectral al lui `A`:

\[
\ker(\Box_F(6)) = E_A(12)\oplus E_A(6\phi)\oplus E_A(6\phi').
\]

Echivalent, sub acțiunea lui `2I`,

\[
\ker(\Box_F(6)) = \rho_0 \oplus \rho_1 \oplus \rho_8.
\]

### P4. Stable gap ratio

Pentru fiecare fibrară `F`,

\[
\frac{\lambda_1(L_{\mathrm{cross}}(F))}{\lambda_1(L_{\mathrm{fiber}}(F))}
= 5.
\]

## Closure Decision

`S06` este închis ca:

- `Computational fact` for the exact count of the 6 discrete fibrations;
- `Theorem-level stable statement inside the verified six-fibration class`
  for `P2`, `P3`, and `P4`.

În limbaj practic al exact-core:

- alegerea fibrării nu schimbă nucleul spectral relevant al construcției;
- rezultatele-cheie folosite mai departe sunt uniforme pe toată clasa discretă
  de fibrări compatibile cu `2I`.

## What Is Not Claimed

Nu se pretinde aici:

- că orice mărime derivată din fibrară este invariantă;
- că semnătura completă a lui `\Box_F` este aceeași pe toate fibrările;
- că orice fibrară Hopf continuă a lui `S^3` intră în această afirmație.

Se pretinde doar uniformitatea exactă a celor patru propoziții de mai sus pe
clasa discretă relevantă.
