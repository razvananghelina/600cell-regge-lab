"""Exact trace-index test for the gauge-prefactor ratio 8:5:2.

Only spaces with an already-defined action of all three gauge factors are
eligible.  The repository inventory is recorded explicitly below; spaces
carrying only 2I, A5, Hodge, or particle-label data are not promoted to gauge
modules by this script.
"""

from fractions import Fraction as F
import sys

failed = False
run = passed = 0


def check(name, condition, detail=""):
    global failed, run, passed
    run += 1
    if condition:
        passed += 1
        print(f"  [PASS] {name}")
    else:
        failed = True
        print(f"  [FAIL] {name}")
    if detail:
        print(f"         {detail}")


def target_ratio(indices):
    t1, t2, t3 = indices
    return 5*t1 == 8*t2 and 2*t2 == 5*t3 and t3 != 0


print("="*72)
print("MATTER TRACE-INDEX TEST")
print("="*72)

target = (F(8), F(5), F(2))
check("prefactors reduce exactly to trace-index ratio 8:5:2",
      target_ratio(target))
check("ratio condition is 5*T1=8*T2 and 2*T2=5*T3",
      5*target[0] == 8*target[1] and 2*target[1] == 5*target[2])

# One SM generation, expressed as left-handed Weyl fields:
# Q=(3,2)_{1/6}, u^c=(bar3,1)_{-2/3}, d^c=(bar3,1)_{1/3},
# L=(1,2)_{-1/2}, e^c=(1,1)_1.  A sterile nu^c contributes zero.
# Convention: T(fundamental SU(n))=1/2 and T1=Tr(Y^2).
sm_terms = {
    "Q":   (F(6)*F(1, 6)**2, F(3)*F(1, 2), F(2)*F(1, 2)),
    "uc":  (F(3)*F(2, 3)**2, F(0),          F(1, 2)),
    "dc":  (F(3)*F(1, 3)**2, F(0),          F(1, 2)),
    "L":   (F(2)*F(1, 2)**2, F(1, 2),       F(0)),
    "ec":  (F(1),             F(0),          F(0)),
}
sm = tuple(sum((term[i] for term in sm_terms.values()), F(0)) for i in range(3))
check("one SM Weyl generation has indices (10/3,2,2)",
      sm == (F(10, 3), F(2), F(2)), f"indices={sm}")
check("ordinary-hypercharge SM ratio is 5:3:3, not 8:5:2",
      not target_ratio(sm) and tuple(x*F(3, 2) for x in sm) == (F(5), F(3), F(3)))
sm_gut = (F(3, 5)*sm[0], sm[1], sm[2])
check("GUT-normalized SM hypercharge gives equal indices (2,2,2)",
      sm_gut == (F(2), F(2), F(2)))

# Finite inventory.  `common_action` means that explicit actions of u(1),
# su(2), and su(3) on the same vector space have been constructed, not merely
# that its dimension or labels resemble matter content.
inventory = {
    "120 vertices / 2I regular sectors": {
        "common_action": False, "derived_action": "left/right 2I"},
    "nine affine-E8 McKay node representations": {
        "common_action": False, "derived_action": "2I representation ring"},
    "nine (a,b) fermion slots": {
        "common_action": False, "derived_action": "labels and mass grading only"},
    "12 Hopf fiber amplitudes before brackets": {
        "common_action": False, "derived_action": "A5 permutation/Hodge data"},
    "other scalar/vector/Hodge spectral sectors": {
        "common_action": False, "derived_action": "spectral and A5/2I data"},
    "1+3+8 gauge algebra in its adjoint module": {
        "common_action": True, "indices": (F(0), F(2), F(3))},
}
check("inventory contains all five requested candidate families",
      len(inventory) == 6)
eligible = {name: item["indices"] for name, item in inventory.items()
            if item["common_action"]}
check("only the conditional gauge adjoint has all three explicit actions",
      eligible == {"1+3+8 gauge algebra in its adjoint module": (F(0), F(2), F(3))})
check("adjoint indices are (0,2,3)",
      next(iter(eligible.values())) == (F(0), F(2), F(3)))
check("no derived common module has index ratio 8:5:2",
      not any(target_ratio(indices) for indices in eligible.values()))

# If the ordinary one-generation SM content is retained and normalized so
# T3=2 (hence scale k=1), any additional color-singlet content would have to
# supply precisely these missing indices.  This is a condition, not a claim
# that such extra multiplets are present.
deficit = tuple(target[i]-sm[i] for i in range(3))
check("SM-like target deficit at T3=2 is exactly (14/3,3,0)",
      deficit == (F(14, 3), F(3), F(0)), f"deficit={deficit}")

print("\nClassification:")
print("  DERIVED: target equations and SM/adjoin trace indices")
print("  DERIVED (negative): no derived matter module exists; conditional adjoint fails 8:5:2")
print("  OPEN: construct a common matter action; particle labels alone do not define one")
print(f"\nTOTAL: {passed}/{run} tests PASSED")
sys.exit(1 if failed else 0)
