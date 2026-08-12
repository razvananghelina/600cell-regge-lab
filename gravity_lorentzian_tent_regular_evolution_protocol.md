# Preregistered protocol: regularity of the asymmetric Lorentzian tent equation

Date: 2026-08-12

Status at registration:
**PROTOCOL ONLY -- THE POLE-DERIVATIVE SIGN HAS NOT BEEN COMPUTED**

## 1. Framing correction

The previous result correctly states that bare icosahedral symmetry does not
select the target-found four-shell boundary profile. It is too strong,
however, to require a gravitational evolution law to select one unique
boundary state without initial data.

In canonical Regge calculus, the action is Hamilton's principal function.
New edge lengths can be free data at one move and can be restricted by
pre-constraints from later moves. An invariant action may also admit a full
symmetry orbit of asymmetric states. Therefore:

- failure of bare geometry to choose `q_u` is still a non-prediction result;
- it is not, by itself, failure of a dynamical law;
- the immediate dynamical test is whether the internal pole equation is
  locally solvable for the pole when boundary data are supplied.

This distinction agrees with the established canonical tent-move framework
of Dittrich and Hoehn,
[arXiv:1108.1974](https://arxiv.org/abs/1108.1974). That framework is external
theory input; it does not certify this repository's Lorentzian conventions.

## 2. Complete frozen hypotheses

Use exactly the carrier, causal branch, action convention and boundary
witness certified in protocol/result commits `4a63f66` and `cc71574`:

1. `T_v=[v,v']*L_v`, with `L_v` the combinatorial icosahedron;
2. old cone and link squared lengths equal `a^2`;
3. the pole squared interval is `-rho a^2`, `rho>0`;
4. final cone squared lengths are `q_u a^2`;
5. all final tetrahedra are strictly spacelike;
6. all four-simplices have signature `(-,+,+,+)`;
7. internal timelike hinges use the ordinary real Lorentzian Regge angle and
   real area magnitude already checked against explicit Minkowski normals;
8. the volume/cosmological coefficient, matter and higher-curvature terms
   vanish;
9. the only bulk variable is the pole; `q_u` are boundary data.

The frozen point is

```text
rho_0=1/4,
(q_0,q_1,q_2,q_3)=(x*,3/2,4/5,3/2),
x* in [11/25,9/20],
E(rho_0,q)=0.
```

The root `x*` is defined by the certified rational bracket and unique-zero
theorem, not by treating its decimal expansion as exact.

## 3. Equation and regularity criterion

For every internal hinge `u`, let

```text
epsilon_u=2*pi-sum_(five incident simplices) theta,
w_u=d(A_u/a^2)/d rho>0,
E(rho,q)=sum_u epsilon_u w_u.
```

After removing the common positive conversion factor between `rho` and pole
proper length, `E=0` is the internal Regge equation.

The load-bearing test is

```text
partial E/partial rho != 0
```

at the root. If it holds, the implicit-function theorem gives a unique smooth
local pole `rho=rho(q)` for all sufficiently nearby boundary data. This is a
regular local evolution relation, not a unique vacuum and not a global time
step.

If `tau=sqrt(rho)` is used, the physical stationarity equation differs by a
nonzero factor `2 tau`. At a root its second derivative has the same
nondegeneracy/sign information as `partial E/partial rho` because
`tau_0=1/2`.

## 4. Frozen verification procedure

The verifier must rebuild the combinatorial icosahedron and the exact angle,
area and pole equation rather than differentiating stored decimals.

Use first-order dual numbers over Arb balls to compute derivatives. It must:

1. reproduce the certified endpoint signs and `partial E/partial x>0`;
2. evaluate `partial E/partial rho` on the entire closed `x` bracket at
   `rho=1/4`;
3. if dependency inflation prevents one-box certification, subdivide the
   rational bracket into exactly `16` equal closed rational boxes and require
   the same strict sign on all of them;
4. compute high-precision values at the isolated root for
   `E_rho`, all twelve individual `E_(q_u)`, and the four shell-coordinate
   derivatives;
5. verify that individual derivatives are constant on the stabilizer orbits
   and that a shell derivative equals the sum of its individual derivatives;
6. report the implicit responses
   `d rho/d q_u=-E_(q_u)/E_rho` and the four collective-shell responses;
7. verify these derivatives independently by centered finite differences of
   the complete pole equation, with a declared relative tolerance `2e-6`;
8. verify that the root remains strictly inside the already certified causal
   domain; no new causal claim may be inferred merely from the derivative.

No desired sign, magnitude, simple constant or Planck-scale comparison is a
target. A derivative of either sign passes if it is rigorously separated from
zero.

## 5. Canonicity and scope attack

The test can prove only regularity of the single bulk equation. It cannot
claim the full pre/post Legendre map because that additionally requires:

- the complete Lorentzian boundary/corner action with fixed branch and sign;
- pre- and post-momenta for all boundary edges;
- the mixed Lagrangian Hessian and its constraint/null directions;
- compatibility with adjacent tent moves.

The existing Euclidean full-boundary-action control is not silently
analytically continued. Lorentzian spacelike boundary hinges have boost/corner
branches that must be fixed independently before computing momenta.

Likewise, the theorem does not select the four-shell boundary state. It says
that, once nearby boundary data are supplied, the internal pole is or is not
locally determined.

## 6. Decision boundary

- **DERIVED REGULAR LOCAL POLE:** Arb separates `E_rho` from zero on the root
  bracket and independent derivatives agree. Then the implicit-function
  theorem licenses a unique local `rho(q)` near the root.
- **DERIVED DEGENERATE POLE:** the exact derivative vanishes or a certified
  enclosure forces zero.
- **OPEN/INCONCLUSIVE:** interval dependency contains zero and the frozen
  16-box subdivision cannot decide the sign.
- **REFUTED WITNESS:** rebuilding the prior equation, causal branch or root
  signs fails.

Even the strongest passing result is **STRUCTURAL LOCAL DYNAMICS**, not a
selected physical clock, constraint algebra, global foliation, light cone or
Planck unit.

Only the targeted verifier and a static registry audit may run. No full suite
and no PDF build.
