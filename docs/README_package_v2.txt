SPOOFING DETECTION + RECOVERY PACKAGE  (complete set from this session)
=======================================================================

UNZIP EVERYTHING INTO ONE MATLAB FOLDER. Filenames are already correct
(underscores) - do not rename anything.

CONTENTS (22 files)
-------------------
Struct classes (definitions only):
  KujurConfig.m   NavState.m   PoolState.m   NavMode.m
  NavInfo.m       PoolReport.m SsResult.m

Logic functions:
  protected_nav_step.m        the FSM (NOMINAL / COAST / PROBATION)
  monitor_pool_step.m         the windows (SS every epoch, CPI at close)
  monitor_pool_open_window.m  monitor_pool_clear.m
  ss_monitor_step.m           SS test, Eq. 49-52
  kalman_update_step.m        one KF epoch (harness; replace with YOUR EKF)
  ins_coast_step.m            INS-only propagation (E28/E31)
  revalidation_test.m         GNSS-vs-coast chi-square (end-of-attack)

Integration templates (your two-function architecture):
  spoof_monitor_2hz.m         persistent wrapper for your 2 Hz function
  ins_100hz_template.m        the three additions to your 100 Hz function

Harness:
  recovery_nav_sim.m  run_recovery.m  generate_test_data.m

Docs:
  QUICKSTART_integration.m    Implementation_Guide_v2.pdf

KEEP FROM YOUR CURRENT FOLDER (8 files - already there, do NOT delete)
---------------------------------------------------------------------
  cpi_monitor_window.m   kujur_params.m    solve_N_min.m
  compute_PMD_eq38.m     build_Phi_Q.m     test_all.m
  ss_monitor_update.m    dual_monitor.m

These are YOUR original files; your validated results were produced
with them, so they are the reference copies. The new code calls
cpi_monitor_window and the design tools directly.

TOTAL FOLDER = 22 (this zip) + 8 (yours) = 30 files.

TEST ORDER
----------
  1. test_all        -> all 10 tests pass (validates the shared math)
  2. run_recovery    -> full end-to-end run; compare one plot against
                        your validated diagonal-attack results
  3. which protected_nav_step monitor_pool_step ss_monitor_step ...
     (any 'not found' = a file is missing or misnamed)
