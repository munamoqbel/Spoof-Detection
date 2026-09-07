# Octave regression shim (NOT for the simulation)

Lets `run_recovery.m` and `run_scheduler_test.m` run headless under GNU
Octave (tested 8.4 + statistics package) without MATLAB:

- `CST_spfMode.m` — Octave has no `enumeration` classes; this exposes the
  same three `uint8` constants as Constant properties.
- `pf.m` — Octave cannot call a script-local function defined after its use.
- plotting stubs (`figure`, `plot`, `subplot`, ...) — headless no-ops.
- `run_recovery_headless.m`, `run_sched_headless.m` — drivers.

Run from the repo root:

    SHIM=$PWD/tools/octave_shim octave --no-gui --quiet tools/octave_shim/run_sched_headless.m

Octave's `rng(42)`/`randn` stream differs from MATLAB's, so the random
geometry and noise differ: expect the same design constants, the same
qualitative events (latch on all 3 axes at attack onset, 4 probation
attempts / 3 vetoes / 1 commit, handback at epoch 208) but not bit-exact
epochs. Bit-exact regression against `reference/*.txt` needs MATLAB.

Do NOT copy this folder into the target simulation.
