# Prior-art gate: relational temporal composition of the homogeneous dust map

Date: 2026-08-21

Status: written before deriving a same-state half-step branch, before running
any new root solve, and before comparing one coarse slab with two fine slabs.

## 1. Exact object and complete hypotheses

Use the already certified homogeneous cellular 600-cell action

```text
S(L_minus,L_plus,rho;M)
 = 360(L_minus+L_plus)h[2*pi-5*acos(c)]
 + 600*sqrt(3)(L_minus^2-L_plus^2)asinh(b)
 - 8*pi*M*sqrt(rho),
```

where

```text
Delta = L_plus-L_minus,
h     = sqrt(rho+Delta^2/4),
c     = (Delta^2+2*rho)/(2(Delta^2+3*rho)),
b     = Delta/sqrt(8(Delta^2+3*rho)).
```

Fix `L0=1`, zero cosmological constant, the conserved mass

```text
M=(90/pi)[2*pi-5*acos(1/3)],
```

the positive Lorentzian branch, homothetic boundaries, and the action-derived
pre/post sign convention.  The current boundary phase state contains the
spatial scale and its canonical momentum.  The squared strut proper duration
`rho` is an internal variable and obeys `partial S/partial rho=0`.

The new question is deliberately narrower than “does the theory derive a
tick?”  Ask whether the current action admits a local weak-lapse comparison
from one and the same initial canonical state:

```text
one nominal slab of duration e
versus
two nominal slabs of duration e/2,
```

with every slab satisfying its lapse equation and every shared boundary
satisfying its canonical seam equation.  The compared observables are final
scale, final momentum and accumulated dust-worldline proper length.

No absolute unit, continuum target, experimental datum or desired convergence
order is part of the definition.

## 2. Why the existing lambda family is not this test

The committed weak-lapse family starts at

```text
p_initial(lambda)=lambda*k0.
```

Therefore its `lambda=1/2` history and its `lambda=1` history do not have the
same initial canonical momentum.  Comparing them as a half-step and a full
step would change the physical initial state while changing the nominal
duration.  This is **DERIVED FROM THE REGISTERED DEFINITIONS** and makes that
old family inadmissible as a step-composition test.

For the new test, if the coarse incoming momentum is `p0`, the first fine slab
must receive exactly the same `p0`, not the static momentum associated with
its own half lapse.  A fine branch must be established by its equations; it
may not be chosen by proximity after inspecting the coarse answer.

## 3. Primary prior art

- Brown and Kuchar, *Dust as a Standard of Space and Time in Canonical Quantum
  Gravity*, arXiv:`gr-qc/9409001`, make dust proper time and its conjugate
  momentum canonical variables and resolve the Hamiltonian constraint with
  respect to the dust-time momentum:
  <https://arxiv.org/abs/gr-qc/9409001>.
- Dittrich and Hoehn, *Canonical simplicial gravity*, arXiv:`1108.1974`, use
  Hamilton's principal function and pre/post Legendre data to generate
  composable simplicial evolution:
  <https://arxiv.org/abs/1108.1974>.
- Marrero, Martin de Diego and Martinez, *On the exact discrete Lagrangian
  function for variational integrators*, arXiv:`1608.01586`, show that an
  exact discrete Lagrangian generates the fixed-time Hamiltonian flow and
  provides the reference object for variational error analysis:
  <https://arxiv.org/abs/1608.01586>.
- Bahr and Dittrich, *Improved and Perfect Actions in Discrete Gravity*,
  arXiv:`0907.4323`, define improved actions by solving refined equations
  subject to coarse boundary data and identify refinement independence with
  restored reparametrization symmetry:
  <https://arxiv.org/abs/0907.4323>.
- De Felice and Fabri evolve a dust-filled 600-cell with a Sorkin/Regge
  scheme; this establishes multi-step 600-cell dust evolution, not the
  same-state composition certificate below:
  <https://arxiv.org/abs/gr-qc/0009093>.

## 4. KNOWN / CONTROL / OPEN

### KNOWN

- The cellular action above is exactly the common restriction of all
  admissible homogeneous staircase subdivisions.
- It generates the existing one-step and multi-step canonical data.
- The static mass-balanced family leaves `rho` arbitrary.
- The all-index weak-lapse trajectory has constant-acceleration leading
  coefficients when each nominal lapse is accompanied by its own static
  incoming momentum.
- Exact/perfect discrete Lagrangians and continuum flows have a composition
  law; generic finite discretizations need only converge to it.

### CONTROL

- Hold `(L0,p0,M)` literally identical between the coarse and fine histories.
- Require all lapse and seam equations; fixing `rho` while dropping its
  equation is a different, gauge-fixed system and does not count.
- Enumerate all real weak-lapse coefficient branches permitted by the frozen
  ansatz.  Do not select a half-step branch by a numerical seed alone.
- Calibrate the state alignment and composition code on an exact
  constant-force discrete Lagrangian.
- Include the hostile old construction with the fine static momentum and
  require the same-state gate to reject it.

### OPEN

- Whether a real same-state half-step branch exists on the current action.
- If it exists, whether two fine steps converge to one coarse step in final
  scale, final momentum and accumulated proper time.
- The convergence order and whether a refinement limit restores a freely
  chosen relational interval.
- Spatial-refinement stability, anisotropic propagation and external novelty.

## 5. Framing attack

The notation `Phi_h` is not currently justified merely by writing
`rho=h^2 exp(r)`: because `rho` is varied, changing the nominal `h` can be
only a coordinate change in `r`, or it can lead to a different solution
branch.  The branch must be proved to have actual proper duration asymptotic
to `h`.

If no same-state half-step branch exists, the correct result is not that a
fundamental tick has been discovered.  It is that the finite action supplies
a lapse-dependent pseudo-constraint and fails the present temporal
refinement gate.  Calling that selected lapse fundamental would require an
additional principle excluding refined carriers.

Conversely, if the one-versus-two comparison converges, that is evidence for
a relational continuum map and evidence *against* interpreting a single
carrier slab as a fundamental indivisible tick.  It still does not select an
absolute time unit because the already proved global scale covariance
remains in force.

## 6. Proposed evidential difference

The only possible new result is an internally preregistered, same-state
temporal-composition certificate or a bounded no-go for that certificate on
the homogeneous cellular action.  The general mechanism is **KNOWN**.
External novelty of any narrow coefficient identity is **OPEN**; a search
cannot establish novelty.

