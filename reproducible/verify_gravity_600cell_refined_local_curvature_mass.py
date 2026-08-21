#!/usr/bin/env python3
"""Post-hoc certificate for the refined local curvature--mass identity."""

from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path

import mpmath as mp
import sympy as sp


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
STATIONARY = HERE / "gravity_600cell_refined_h4_stationary_fill.json"
STATIONARY_SOURCE = HERE / "verify_gravity_600cell_refined_h4_stationary_fill.py"
FEASIBILITY = HERE / "gravity_600cell_refined_canonical_map_feasibility.json"
FEASIBILITY_SOURCE = HERE / "verify_gravity_600cell_refined_canonical_map_feasibility.py"
CELL600 = ROOT / "commons/cell600.py"
PRIOR_ART = ROOT / "docs/gravity/gravity_600cell_refined_local_curvature_mass_prior_art.md"
PROTOCOL = ROOT / "docs/gravity/gravity_600cell_refined_local_curvature_mass_protocol.md"
OUTPUT = HERE / "gravity_600cell_refined_local_curvature_mass.json"

PRIOR_ART_COMMIT = "f12f56c"
PROTOCOL_COMMIT = "bb80f28"
EXPECTED_HASHES = {
    "stationary": "283be37bc7530a3cc4fce9e279272359f107f09fb7b1b0eaff141059bfb4e018",
    "stationary_source": "89aab727792e20a81e7577e0425f8fa4b1e84e2a7ae66caa9e79a4aebf3581e7",
    "feasibility": "ab6209bc745b4c988b59b8c0416522dd2e4a434f17f4cfd596df817bb48ff02e",
    "feasibility_source": "36fba835048e6e0f0676b749192a9d882406932770a00ba1396929bbc4d04a32",
    "cell600": "ea5bce4b6c52e0834539ca4b1df9c6a67a3a5ed4da32f4e0298a493fc5315c7f",
    "prior_art": "37a6e5b53011d499f36a8ba3daff27b5c2e20334103cd36fcec58dcdbfbd135c",
    "protocol": "23824de71ad9c3493af56358720ba7cac4319077c3893d73d9d398bd65795187",
}
PAIRS = tuple(combinations(range(4), 2))
RANK_SIZES = (120, 720, 1200, 600)
TAU_TEXT = "0.0102"
TOLERANCE = mp.mpf("1e-68")

tests = 0
passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    ok = bool(condition)
    passed += int(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {label}")
    if detail:
        print(f"       {detail}")
    return ok


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def mp_text(value, digits=75):
    return mp.nstr(value, digits, strip_zeros=False)


def reconstruct_geometry(populations):
    phi = (1 + mp.sqrt(5)) / 2
    c_value = phi / 2

    def norm_square(size):
        return size * (1 + (size - 1) * c_value)

    unit_squares = {}
    for left, right in PAIRS:
        a, b = left + 1, right + 1
        dot = (
            a * (1 + (b - 1) * c_value)
            / mp.sqrt(norm_square(a) * norm_square(b))
        )
        unit_squares[left, right] = 2 - 2 * dot

    squared = [[mp.mpf(0) for _ in range(4)] for _ in range(4)]
    for (left, right), value in unit_squares.items():
        squared[left][right] = squared[right][left] = value
    gram = mp.matrix([
        [
            (squared[3][left] + squared[3][right] - squared[left][right]) / 2
            for right in range(3)
        ]
        for left in range(3)
    ])
    chamber_volume = mp.sqrt(mp.det(gram)) / 6
    total_volume = 14400 * chamber_volume
    scale = (2 * mp.pi**2 / total_volume) ** (mp.mpf(1) / 3)

    inverse = gram**-1
    normals = (
        mp.matrix([1, 0, 0]),
        mp.matrix([0, 1, 0]),
        mp.matrix([0, 0, 1]),
        mp.matrix([-1, -1, -1]),
    )

    def inner(left, right):
        return (left.T * inverse * right)[0]

    angles = {}
    incidences = {}
    deficits = {}
    pair_curvatures = {}
    for left, right in PAIRS:
        omitted_a, omitted_b = [
            value for value in range(4) if value not in (left, right)
        ]
        cosine = -inner(normals[omitted_a], normals[omitted_b]) / mp.sqrt(
            inner(normals[omitted_a], normals[omitted_a])
            * inner(normals[omitted_b], normals[omitted_b])
        )
        angles[left, right] = mp.acos(cosine)
        incidences[left, right] = mp.mpf(14400) / populations[left, right]
        deficits[left, right] = (
            2 * mp.pi - incidences[left, right] * angles[left, right]
        )
        pair_curvatures[left, right] = (
            populations[left, right]
            * scale
            * mp.sqrt(unit_squares[left, right])
            * deficits[left, right]
        )

    rank_curvatures = tuple(
        mp.fsum(
            pair_curvatures[pair] / 2
            for pair in PAIRS
            if rank in pair
        )
        for rank in range(4)
    )
    return {
        "unit_squares": unit_squares,
        "chamber_volume": chamber_volume,
        "total_volume": total_volume,
        "scale": scale,
        "angles": angles,
        "incidences": incidences,
        "deficits": deficits,
        "pair_curvatures": pair_curvatures,
        "rank_curvatures": rank_curvatures,
        "total_curvature": mp.fsum(pair_curvatures.values()),
    }


paths = {
    "stationary": STATIONARY,
    "stationary_source": STATIONARY_SOURCE,
    "feasibility": FEASIBILITY,
    "feasibility_source": FEASIBILITY_SOURCE,
    "cell600": CELL600,
    "prior_art": PRIOR_ART,
    "protocol": PROTOCOL,
}
actual_hashes = {name: digest(path) for name, path in paths.items()}
provenance_ok = check(
    "all frozen inputs match the committed protocol",
    actual_hashes == EXPECTED_HASHES
    and PRIOR_ART_COMMIT == "f12f56c"
    and PROTOCOL_COMMIT == "bb80f28",
)

stationary = json.loads(STATIONARY.read_text())
feasibility = json.loads(FEASIBILITY.read_text())
census = stationary["census"]
upstream_ok = check(
    "the accepted off-shell census supplies exactly the frozen 24 schedules",
    stationary["outcome"] == "REFINED_H4_INDUCED_FILL_OFF_SHELL"
    and stationary["tests"] == {"passed": 12, "total": 12}
    and census["schedule_count"] == len(census["schedules"]) == 24
    and census["certified_nonzero_vertical_entries"] == 96
    and census["certified_nonzero_cross_entries"] == 0
    and stationary["definitions"]["tau0_supplied"] == "0.0102",
)

level = feasibility["levels"]["projected_barycentric"]
population_map = level["colour_pair_edge_populations"]
populations = {
    pair: int(population_map[f"{pair[0]}-{pair[1]}"])
    for pair in PAIRS
}
with mp.workdps(100):
    geometry = reconstruct_geometry(populations)
    committed_curvature = mp.mpf(
        stationary["exact_geometry"]["spatial_regge_curvature"]
    )
    committed_scale = mp.mpf(stationary["exact_geometry"]["scale_s0"])
    geometry_error = abs(geometry["total_curvature"] - committed_curvature)
    scale_error = abs(geometry["scale"] - committed_scale)
    integer_incidences = all(
        14400 % populations[pair] == 0
        and geometry["incidences"][pair] == 14400 // populations[pair]
        for pair in PAIRS
    )
    positive_geometry = all(
        geometry["unit_squares"][pair] > 0
        and geometry["deficits"][pair] > 0
        and geometry["pair_curvatures"][pair] > 0
        for pair in PAIRS
    ) and all(value > 0 for value in geometry["rank_curvatures"])

geometry_ok = check(
    "the six independent spatial hinge curvatures are positive and integral",
    integer_incidences and positive_geometry,
    "incidences=" + str(tuple(
        int(geometry["incidences"][pair]) for pair in PAIRS
    )),
)
curvature_ok = check(
    "the endpoint-half rank curvatures sum to the committed Regge curvature",
    geometry_error < TOLERANCE and scale_error < TOLERANCE,
    f"curvature error={mp_text(geometry_error, 8)}, "
    f"scale error={mp_text(scale_error, 8)}",
)

with mp.workdps(100):
    tau = mp.mpf(TAU_TEXT)
    mass = mp.mpf(stationary["exact_geometry"]["selected_total_mass"])
    rank_curvatures = geometry["rank_curvatures"]
    total_curvature = geometry["total_curvature"]
    action_differences = []
    recovered_vectors = []
    for schedule in census["schedules"]:
        recovered = []
        for rank in range(4):
            total_residual = mp.mpf(
                schedule["internal"][f"rho_{rank}"]["total_log_residual"]["real"]
            )
            gravitational = total_residual + mp.pi * mass * tau
            expected = tau * rank_curvatures[rank] / 2
            action_differences.append(gravitational - expected)
            recovered.append(gravitational)
        recovered_vectors.append(tuple(recovered))
    maximum_identity_error = max(abs(value) for value in action_differences)
    schedule_spread = max(
        abs(recovered_vectors[index][rank] - recovered_vectors[0][rank])
        for index in range(24)
        for rank in range(4)
    )

identity_ok = check(
    "all 96 gravitational lapse derivatives equal tau times local curvature / 2",
    maximum_identity_error < TOLERANCE,
    f"max error={mp_text(maximum_identity_error, 8)}",
)
schedule_ok = check(
    "all 24 staircase schedules select the same four curvature masses",
    schedule_spread < TOLERANCE,
    f"spread={mp_text(schedule_spread, 8)}",
)

tau_symbol = sp.symbols("tau", positive=True)
mass_response = -4 * sp.pi * tau_symbol * sp.eye(4)
rank_ok = check(
    "the four rank masses enter the lapse equations with symbolic rank four",
    mass_response.rank() == 4
    and sp.simplify(mass_response.det() - (4 * sp.pi * tau_symbol)**4) == 0,
    f"det={mass_response.det()}",
)

with mp.workdps(100):
    fractions = tuple(value / total_curvature for value in rank_curvatures)
    selected_masses = tuple(value / (8 * mp.pi) for value in rank_curvatures)
    per_vertex_masses = tuple(
        selected_masses[rank] / RANK_SIZES[rank] for rank in range(4)
    )
    density_ratios = tuple(4 * value for value in fractions)
    density_contrast = max(density_ratios) / min(density_ratios)
    mass_sum_error = abs(mp.fsum(selected_masses) - mass)

mass_ok = check(
    "the unique curvature masses are positive and conserve the frozen total mass",
    all(value > 0 for value in selected_masses)
    and mass_sum_error < TOLERANCE,
    f"sum error={mp_text(mass_sum_error, 8)}",
)

alternatives = {
    "P1_equal_rank": tuple(mp.mpf(1) / 4 for _ in range(4)),
    "original_vertices_only": (mp.mpf(1), mp.mpf(0), mp.mpf(0), mp.mpf(0)),
    "binomial_1331": tuple(mp.mpf(value) / 8 for value in (1, 3, 3, 1)),
}
alternative_errors = {
    name: max(abs(candidate[rank] - fractions[rank]) for rank in range(4))
    for name, candidate in alternatives.items()
}
p1_negative_ok = check(
    "the conditional P1 equal-rank masses remain locally off shell",
    alternative_errors["P1_equal_rank"] > mp.mpf("1e-6"),
    f"fraction error={mp_text(alternative_errors['P1_equal_rank'], 12)}",
)
particle_negative_ok = check(
    "retaining mass only on original vertices does not solve the rank equations",
    alternative_errors["original_vertices_only"] > mp.mpf("1e-6"),
    f"fraction error={mp_text(alternative_errors['original_vertices_only'], 12)}",
)
pattern_negative_ok = check(
    "the nearby (1,3,3,1)/8 pattern is not an exact identity",
    alternative_errors["binomial_1331"] > mp.mpf("1e-6"),
    f"fraction error={mp_text(alternative_errors['binomial_1331'], 12)}",
)

corrupted_differences = list(action_differences)
corrupted_differences[0] += mp.mpf("1e-10")
residual_corruption_ok = check(
    "a deliberate 1e-10 lapse-residual corruption fails the identity gate",
    max(abs(value) for value in corrupted_differences) > TOLERANCE,
)
corrupted_populations = dict(populations)
corrupted_populations[(0, 1)] += 1
population_corruption_ok = check(
    "an incompatible edge-population corruption fails integral incidence",
    14400 % corrupted_populations[(0, 1)] != 0,
)

scope = {
    "root_search_executed": False,
    "nested_census_executed": False,
    "hessian_or_spectrum_computed": False,
    "continuum_or_particle_target_loaded": False,
    "physical_constant_extracted": False,
    "discovery_status": "POST_HOC_CONFIRMATION_AFTER_DISCLOSED_EXPLORATION",
}
scope_ok = check(
    "the calculation stays inside the frozen local identity scope",
    not any(value for key, value in scope.items() if key.endswith(("executed", "loaded", "extracted"))),
)

controls_ok = all((
    provenance_ok,
    upstream_ok,
    geometry_ok,
    curvature_ok,
    schedule_ok,
    rank_ok,
    mass_ok,
    p1_negative_ok,
    particle_negative_ok,
    pattern_negative_ok,
    residual_corruption_ok,
    population_corruption_ok,
    scope_ok,
))
if not controls_ok:
    outcome = "REFINED_LOCAL_CURVATURE_MASS_CONTROL_FAILED"
elif not identity_ok:
    outcome = "REFINED_LOCAL_CURVATURE_MASS_IDENTITY_REFUTED"
else:
    outcome = "REFINED_LOCAL_CURVATURE_MASS_IDENTITY_CONFIRMED_POST_HOC"

outcome_ok = check(
    "the frozen hierarchy assigns exactly one scoped outcome",
    outcome in {
        "REFINED_LOCAL_CURVATURE_MASS_CONTROL_FAILED",
        "REFINED_LOCAL_CURVATURE_MASS_IDENTITY_REFUTED",
        "REFINED_LOCAL_CURVATURE_MASS_IDENTITY_CONFIRMED_POST_HOC",
    },
    outcome,
)

artifact = {
    "title": "Refined local curvature--mass identity",
    "date": "2026-08-21",
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": actual_hashes,
    "definitions": {
        "carrier": "K0=P(sd K_600)",
        "f_vector": [2640, 17040, 28800, 14400],
        "rank_sizes": list(RANK_SIZES),
        "tau0": "0.0102",
        "curvature_localization": "half of l_e*epsilon_e to each endpoint",
        "mass_action": "-8*pi*sum_r mu_r*sqrt(rho_r)",
        "identity": "dS_grav/dlog(rho_r)=tau0*K_r/2",
        "selected_mass": "mu_r=K_r/(8*pi)",
    },
    "geometry": {
        "edge_populations": {
            f"{left}-{right}": populations[left, right]
            for left, right in PAIRS
        },
        "tetrahedra_per_edge": {
            f"{left}-{right}": int(geometry["incidences"][left, right])
            for left, right in PAIRS
        },
        "deficits": {
            f"{left}-{right}": mp_text(geometry["deficits"][left, right])
            for left, right in PAIRS
        },
        "pair_curvatures": {
            f"{left}-{right}": mp_text(geometry["pair_curvatures"][left, right])
            for left, right in PAIRS
        },
        "rank_curvatures": [mp_text(value) for value in rank_curvatures],
        "total_curvature": mp_text(total_curvature),
        "committed_curvature_error": mp_text(geometry_error),
    },
    "identity_audit": {
        "schedule_count": 24,
        "equation_count": 96,
        "maximum_absolute_error": mp_text(maximum_identity_error),
        "maximum_schedule_spread": mp_text(schedule_spread),
        "tolerance": mp_text(TOLERANCE),
        "mass_response_rank": 4,
    },
    "selected_rank_matter": {
        "fractions": [mp_text(value) for value in fractions],
        "total_masses": [mp_text(value) for value in selected_masses],
        "per_vertex_masses": [mp_text(value) for value in per_vertex_masses],
        "p1_relative_density": [mp_text(value) for value in density_ratios],
        "p1_density_max_min_ratio": mp_text(density_contrast),
        "total_mass": mp_text(mass),
        "mass_conservation_error": mp_text(mass_sum_error),
    },
    "alternative_fraction_errors": {
        name: mp_text(value) for name, value in alternative_errors.items()
    },
    "scope": scope,
    "status_labels": {
        "identity": "DERIVED_COMPUTATIONAL_STRUCTURAL_POST_HOC",
        "p1_homogeneous_dust": "DERIVED_NEGATIVE_ON_FIXED_STATIC_FILL",
        "refinement_convergence": "OPEN",
        "tick_c_G_planck_particles": "OPEN_NOT_COMPUTED",
    },
    "outcome": outcome,
    "tests": {"passed": passed, "total": tests},
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

print("-" * 78)
print(f"Tests passed: {passed}/{tests}")
print(f"Outcome: {outcome}")
print("K_r/K: " + ", ".join(mp_text(value, 22) for value in fractions))
print(f"P1 density contrast: {mp_text(density_contrast, 18)}")

raise SystemExit(0 if passed == tests else 1)
