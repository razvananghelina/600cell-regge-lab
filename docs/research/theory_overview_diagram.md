# Teoria 600-Cell: Diagrama de Ansamblu
## "One Integer to Rule Them All" - Flow Diagram (March 2026, v5)

```
                        +=======================+
                        |   a1! = 4a1(a1+1)     |
                        |   SOLUTIE UNICA: a1=5  |
                        +===========+===========+
                                    |
            +-----------------------+-----------------------+
            |                       |                       |
            v                       v                       v
   +================+    +==================+    +==================+
   | phi=(1+sqrt5)/2|    |  N=a1!/(a1-2)!  |    |  b1=a1+1=6      |
   | Golden Ratio    |    |  =120 vertices   |    |  First Betti     |
   +=======+========+    |  600-cell / 2I   |    +========+=========+
           |              +========+=========+             |
           |                       |                       |
           |    +------------------+------------------+    |
           |    |                  |                   |    |
           |    v                  v                   v    |
           |  +===========+  +=================+  +===================+
           |  | McKay      |  | Spectral        |  | TQFT             |
           |  | Graph      |  | Action          |  | k+2=a1           |
           |  | (aff E8)   |  | D=d+d*          |  | Verlinde ring    |
           |  +=====+=====+  | 2640 eigenvalues |  | =Z[phi]          |
           |        |         | Lam_max=b1*phi^2 |  +========+==========+
           |        |         +=======+==========+           |
           |        |                 |                      |
     ======+========+=========+===+==+===+==================+===========
     GAUGE |        |         |   |      |                   |
           |        v         |   |      |                   |
           |  +-------------+ |   |      |                   |
           |  | 12 = 1+3'+  | |   |      |                   |
           |  |     3+5     | |   |      |                   |
           |  | A5 reps     | |   |      |                   |
           |  +------+------+ |   |      |                   |
           |         |        |   |      |                   |
           |         v        |   |      |                   |
           |  +===============+   |      |                   |
           |  | SU(3)xSU(2)  |   |      |                   |
           |  | xU(1)        |   |      |                   |
           |  | dim(G)=12    |   |      |                   |
           |  +======+========+   |      |                   |
           |         |            |      |                   |
     ======+=========+============+======+==================+===========
     COUPLINGS       |            |      |                   |
           |         |            |      |                   |
     +-----+--+  +---+----+   +--+------+--+    +-----------+-----+
     |alpha_s |  |sin^2tW |   |  alpha     |    | sqrt(2/a1)      |
     |=1/     |  |=b1/    |   | 2*pi*a^2-  |    | one-loop QED    |
     |(2phi^3)|  |(a1^2+1)|   | (N/b1)*    |    | coeff from      |
     |        |  |=6/26   |   | phi^4*a    |    | Verlinde        |
     |0.11%   |  |0.19%   |   | +1=0       |    | S-matrix        |
     +----+---+  +---+----+   |0.0001%     |    +-------+---------+
          |          |         +----+-------+            |
     =====+============+============+===================+================
     MASSES          |              |                    |
          |          |              v                    |
          |          |   +===========================+   |
          |          |   |  m_f = m_e * phi^(5a+6b   |   |
          |          |   |  + d1*(1+sqrt(2/a1)*alpha))|  |
          |          |   |  chi^2/dof = 1.26          |   |
          |          |   +============+===============+   |
          |          |                |                   |
          |          |    +-----------+----------+        |
          |          |    |           |          |        |
          |          |    v           v          v        |
          |          |  +------+  +------+  +------+     |
          |          |  |(a,b) |  |ln|N| |  |-2/a1 |     |
          |          |  |Z[phi]|  |Arake- |  |char   |     |
          |          |  |norm  |  |lov    |  |cross-  |     |
          |          |  |selec-|  |height |  |term    |     |
          |          |  |tion  |  |       |  |(d qk)  |     |
          |          |  +------+  +------+  +------+     |
          |          |    ^           ^          ^        |
          |          |    +---------- +----------+        |
          |          |         GALOIS CASCADE              |
          |          |    (3 levels, 7 branches,           |
          |          |     ONE mechanism: Galois            |
          |          |     obstruction on W=(z,rho))        |
          |          |                                      |
     =====+==========+=====================================+=========
     EW BOSONS       |                                      |
          |          |                                      |
          |    +-----+------+                               |
          |    |  m_Z = m_e  |                               |
          |    |  *phi^25    |                               |
          |    |  *16/15     |                               |
          |    |  (0.28%)    |                               |
          |    +-----+------+                               |
          |          |                                      |
          |          v                                      |
          |    +--------------+    +---------------------+  |
          |    | m_W = m_Z *  |    | m_H^2/m_W^2 =      |  |
          |    | cos(theta_W) +--->| phi^2-16*alpha*phi  |  |
          |    +--------------+    | m_H = 125.0 GeV     |  |
          |                       | (0.20%)              |  |
          |                       +---------------------+  |
          |                                                |
     =====+================================================+=========
     MIXING                                                |
          |                                                |
     +----+--------------------+  +---------------------+  |
     | CKM: n12=3, n13=12     |  | PMNS: TBM from      |  |
     | n23=7                   |  | Galois eigenvalues   |  |
     | delta_CKM=arctan(sqrt5) |  | {0,3,5}             |  |
     | Vacuum offset DERIVED   |  | sin^2(th13)=1/45    |  |
     | (14th characterization) |  | delta_PMNS=197.7 deg |  |
     +-------------------------+  +----------------------+  |
                                                           |
     ======================================================+=========
     NEUTRINOS & CP                                        |
          |                                                |
     +----+--------------+  +---------------------+        |
     | m3=2*m_e/phi^35   |  | theta_QCD = 0 EXACT |        |
     | = 49.53 meV        |  | (Galois parity)     |        |
     | n_seesaw=35        |  +---------------------+        |
     | =N/2-a1^2=60-25   |                                  |
     | M_R ~ 5.3 TeV      |                                  |
     | m2/m3=sqrt(a*phi^3)|                                  |
     +--------------------+                                  |
                                                             |
     ========================================================+=========
     GRAVITY & COSMOLOGY                                     |
          |                                                  |
     +----+---------------+  +--------------+  +-------------+--------+
     | Seeley-DeWitt      |  | Starobinsky  |  | Spectral Data       |
     | A0=11, A1=62       |  | R^2 inflation|  | lambda_max=b1*phi^2 |
     | A2=233=F13         |  | n_s=0.9654   |  | lam1*lam_max=b1^2   |
     | Graviton emergent  |  | r=0.003      |  | Mults=d_i^2         |
     | R(l)->2 (TT exact) |  |              |  | (Peter-Weyl, 2I)    |
     +--------------------+  +--------------+  | Alt zeta=0 (ALL s)  |
          |                                    | z_Pl=2*D^2-d_ST     |
     +----+---------------------------------------------------+------+
     |                    GALOIS DUALITY                               |
     |  sigma: phi -> phi' = (1-sqrt5)/2                               |
     |                                                                 |
     |  VIZIBIL (phi)          |  DARK (phi')                          |
     |  alpha real -> EM       |  alpha' complex -> no EM              |
     |  alpha_s > 0 -> confine |  alpha_s'< 0 -> no confinement       |
     |  m_f                    |  m_f' = m_e^2/m_f                     |
     |  Omega_b                |  Omega_DM ~ (7-phi)*Omega_b ?         |
     |                                                                 |
     |  Norm Sum Rules: SUM N(z_f)=b1=6, SUM N^3=126 (EXACT)         |
     |  CC: alpha^(57-alpha_s) <- Galois cancellation (PATTERN)        |
     +-----------------------------------------------------------------+

     =================================================================
     WHY 600-CELL? (exp415/420)
     +---------------------------------------------------------------+
     |  Tr(A^2)/|G| MINIMIZED by 2I (binary icosahedral)             |
     |  rho(q) = (q+3)(6-q)/(12q), f'<0 for all q                   |
     |  All 6 regular 4D polytopes tested: 600-cell 7/7,             |
     |  others 0/7. Texas Sharpshooter REFUTED.                      |
     +---------------------------------------------------------------+

     =================================================================
     OPEN PROBLEMS (7 remaining)
     +---------------------------------------------------------------+
     |  #1 CC functional form (SERIOS)                                |
     |  #2 Newton's constant G / m_e absolute (STRUCTURAL)           |
     |  #3 Scalaron mass / f4 (SERIOS)                               |
     |  #4 Omega_DM/Omega_b mechanism (PARTIAL)                      |
     |  #5 Baryogenesis factor 1.8x (MINOR)                          |
     |  #6 Seeley-DeWitt c1 discrepancy (MINOR)                      |
     +---------------------------------------------------------------+
```

## Logica "Highway" (fluxul principal al paper-ului)

```
a1=5  -->  phi, N=120, b1=6
              |
              +-->  Gauge: A5 reps -> SU(3)xSU(2)xU(1)
              |
              +-->  Couplings: alpha, alpha_s, sin^2(tW) (ALL from a1)
              |
              +-->  Masses: m_f = m_e*phi^(5a+6b+delta) (Galois Cascade)
              |       |
              |       +-- (a,b) from Z[phi] norm selection
              |       +-- corrections: Arakelov + character + Morrey
              |       +-- one-loop: sqrt(2/a1)*alpha from Verlinde
              |
              +-->  Mixing: CKM (n=3,12,7) + PMNS (TBM + 1/45)
              |       |
              |       +-- Cayley spectral: lam_t/lam_u=phi^4, lam_b/lam_d=phi^2
              |       +-- Vacuum offset "-1" DERIVED (14th char.)
              |
              +-->  EW: m_Z -> m_W -> m_H chain (0.20%)
              |
              +-->  Neutrinos: seesaw n=35, m3=49.53 meV
              |
              +-->  CP: theta_QCD=0, delta_CKM, delta_PMNS
              |
              +-->  Gravity: emergent graviton, Starobinsky inflation
              |       |
              |       +-- Spectral: lambda_max=b1*phi^2, lam1*lam_max=b1^2
              |       +-- z_Planck=2*D^2-d_ST DERIVED
              |
              +-->  Galois: dark sector, norm sum rules, CC (pattern)
```

## Unificarea Mass-Gauge-TQFT (exp445-446)

```
         Wilson lines z_f in Z[phi]
                    |
         +----------+----------+
         |          |          |
         v          v          v
    Galois norm  TQFT vacuum   Gauge
    N(z_f)       E_vac=phi^3   sin^2(tW)
         |          |          |
         v          v          v
    SUM_N=b1=6  =1/(2*alpha_s) =b1/(a1^2+1)
         |          |          |
         +----------+----------+
                    |
              SAME b1 = 6
         (First Betti number)
```

## Spectral Data Summary (exp463-464)

```
    D^2 on 600-cell (2640 eigenvalues)
    |
    +-- Delta_0 (120 dim): 9 levels, mults = d_i^2 (Peter-Weyl)
    |     lambda_1 * lambda_max = b1^2 = 36 EXACT
    |     Eigenvalues = degree - Adj_eigenvalues (Cayley graph)
    |
    +-- Delta_1 (720 dim): 28 levels
    |     min eigenvalue: N(z) = a1 = 5
    |
    +-- Delta_2 (1200 dim): 45 levels
    |     min eigenvalue: N(z) = 1 (UNIT in Z[phi])
    |
    +-- Delta_3 (600 dim): 27 levels
    |     min eigenvalue: N(z) = 1 (UNIT in Z[phi])
    |
    +-- Galois: ALL eigenvalues in Z[phi], sigma permutes them
    +-- Alternating zeta: SUM (-1)^k zeta_k(s) = 0 (ALL s)
    +-- Poincare: log det'(D0-D3) = log det'(D1-D2) = -433.66
    +-- Even=Odd: log det'(D0+D2) = log det'(D1+D3) = 2013.68
    +-- lambda_max = b1*phi^2 = 15.7082 EXACT
    +-- zeta_0(1) = degree = 12, zeta_0(2) = 2
```

## Cele 14 Caracterizari ale lui a1=5

1.  a1! = 4*a1*(a1+1) -- Diophantine
2.  Regular 4-polytope cu decoupling -> 600-cell
3.  McKay correspondence -> affine E8
4.  N_gen = (a1-1)/2 + 1 = 3 integer
5.  sin^2(tW) = b1/(a1^2+1) in (0,1/4)
6.  a1*(a1-5)=0 from moment-Weinberg
7.  Fibonacci TQFT at k+2=a1
8.  Z[phi] is a PID (class number 1)
9.  Tr(A^2)/|G| minimized (variational)
10. E8 root system unique at rank 8
11. Bootstrap: d1=phi iff a1=5 (TQFT vacuum)
12. SUM N(z_f)^3 = 126 iff a1+1=(a1-2)! (norm sum)
13. CKM T-matrix: I2*(a1-2)(a1-5)=0, I3*(a1-1)(a1-5)=0 (concordance)
14. Cayley-CKM vacuum offset: uniform delta=1 <=> a1=5
