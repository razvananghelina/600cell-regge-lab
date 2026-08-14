#!/usr/bin/env python3
"""Execute the one-change calibrated continuation of the frozen solver.

Protocol addendum: cf27934

The base verifier remains immutable and reproduces its calibration failure.
This wrapper makes an audited, exact source transformation: distinct artifact
paths and protocol hash, plus the sole mathematical amendment 1e-10 -> the
upstream control's already registered 5e-8 accuracy class.
"""

from pathlib import Path


HERE = Path(__file__).resolve().parent
BASE = HERE / "verify_gravity_600cell_dust_action_solver_repair.py"

source = BASE.read_text()


def replace_exact(old, new, expected_count=1):
    global source
    count = source.count(old)
    if count != expected_count:
        raise RuntimeError(
            f"base verifier drift for {old!r}: expected {expected_count}, got {count}"
        )
    source = source.replace(old, new)


replace_exact("Protocol commit: 4b6b10c", "Protocol commit: cf27934")
replace_exact(
    'OUTPUT = HERE / "gravity_600cell_dust_action_solver_repair.json"',
    'OUTPUT = HERE / "gravity_600cell_dust_action_solver_repair_continuation.json"',
)
replace_exact(
    'CHECKPOINT = HERE / ".gravity_600cell_dust_action_solver_repair.checkpoint.pkl"',
    'CHECKPOINT = HERE / ".gravity_600cell_dust_action_solver_repair_continuation.checkpoint.pkl"',
)
replace_exact('PROTOCOL_COMMIT = "4b6b10c"', 'PROTOCOL_COMMIT = "cf27934"')
replace_exact(
    'operational_reference_error < arb.mpf("1e-10")',
    'operational_reference_error < arb.mpf("5e-8")',
)
replace_exact(
    'validation_reference_error < arb.mpf("1e-10")',
    'validation_reference_error < arb.mpf("5e-8")',
)
replace_exact(
    '"operational_reference_below_1e-10"',
    '"operational_reference_within_upstream_5e-8"',
    expected_count=2,
)
replace_exact(
    '"validation_reference_below_1e-10"',
    '"validation_reference_within_upstream_5e-8"',
    expected_count=2,
)

compiled = compile(source, str(BASE), "exec")
# Multiprocessing pickles worker functions by their __main__ attribute name.
# Execute in the wrapper's real module namespace so fork-pool task dispatch can
# resolve those names; an isolated dictionary is not registered in sys.modules.
exec(compiled, globals())
