# Spoof-Detection — Kujur et al. (2024) INS monitor with Section-5 recovery

MATLAB implementation of the CPI + solution-separation spoofing monitor from
Kujur, Khanafseh & Pervan, *Optimal INS Monitor for GNSS Spoofer Tracking
Error Detection*, NAVIGATION 71(1) 2024 (navi.629), with the post-detection
recovery FSM (NOMINAL -> COAST -> PROBATION -> handback), packaged for a
two-rate host: a 100 Hz navigation loop and a 2 Hz GNSS/monitor task.

## Runtime path (what the target simulation needs)

Host contract: `docs/HOST_2HZ_WIRING.m`. The gate works in update-increment
space (`x+ - x_bar`), so it is independent of the host's feedback/gain
bookkeeping of `KF.states`.

| file | role |
|---|---|
| `spoofMonitor2hz.m` | persistent-state 2 Hz gate; wraps `protectedNav` |
| `protectedNav.m` | FSM in increment form; host owns its filters, gate returns corrections |
| `monitorPool.m` | overlapping-window bank (N slots) |
| `cpiMonitor.m`, `ssMonitor.m` | per-window CPI (Eq. 33/35) and SS (Eq. 49-52) tests |
| `insCoast.m`, `revalidation.m` | INS-only propagation (E28/E31), chi-square re-validation |
| `STRUCT_SPF.m` | all struct constructors + window open/close |
| `CST_spfParam.m` | solved monitor constants (single source of truth at runtime) |
| `CST_spfMode.m` | mode enumeration |
| `CST_gnssHybrid.m` | **stand-in** (`NO_STATES = 60`); the simulation has its own |

## Harness / offline design (not deployed)

`run_recovery.m` (flat 2 Hz demo + offline design), `run_scheduler_test.m` +
`scheduler_sim.m` (100 Hz + 2 Hz split-rate architecture test),
`generate_test_data.m`, `build_Phi_Q.m`, `kalman_update_step.m`,
`kujur_params.m`, `solve_N_min.m`, `compute_PMD_eq38.m`,
`recovery_nav_sim.m` (a host following `docs/HOST_2HZ_WIRING.m`),
`ins_100hz_template.m` (what the host's 100 Hz function must add).

## Tests

- `run_scheduler_test.m` — 100 Hz + 2 Hz architecture (T1-T6)
- `tests/test_monitors.m` — engine unit + Monte-Carlo checks (replaces the
  dual_monitor-era `test_all.m`, which cannot drive `cpiMonitor` since the
  window length / thresholds moved into `CST_spfParam`)
- `tests/test_variable_numMeas.m` — varying satellite count, padded vs exact
- `run_recovery.m` — end-to-end demo with offline design and plots

## Regression

`reference/` holds the console output of both scripts from the author's
MATLAB run. `tools/octave_shim/` lets both scripts run under GNU Octave
(constants and architecture tests reproduce exactly; event epochs differ by
the RNG stream).
