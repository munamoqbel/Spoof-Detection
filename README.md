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
| `SPF_gate.m` | persistent-state 2 Hz gate; wraps `SPF_protectedNav` |
| `SPF_protectedNav.m` | FSM in increment form; host owns its filters, gate returns corrections |
| `SPF_monitorPool.m` | overlapping-window bank (N slots) |
| `SPF_cpiMonitor.m`, `SPF_ssMonitor.m` | per-window CPI (Eq. 33/35) and SS (Eq. 49-52) tests |
| `SPF_insCoast.m`, `SPF_revalidation.m` | INS-only propagation (E28/E31), chi-square re-validation |
| `matrixInv.m` | STAND-IN for the simulation's SVD inverse (`[inv, invalid] = matrixInv(A)`), used by the CPI and re-validation on the zero-padded covariances; replace with the simulation's own file |
| `STRUCT_SPF.m` | all struct constructors + window open/close |
| `CST_spfParam.m` | solved monitor constants (single source of truth at runtime); `NOMINAL_BIAS` is the nominal GNSS excursion the host tolerates per axis, added to the SS gate and the PL and a deadband of the CPI, 0 in the repository |
| `CST_spfMode.m` | mode enumeration |
| `CST_gnssHybrid.m` | **stand-in** (`NO_STATES = 60`); the simulation has its own |

## Minimalist copies for the simulation

`deploy/` holds the same ten runtime files with the commentary stripped to
a few lines each (identical logic, signatures and struct layouts; bitwise
identical telemetry on the harness scenarios). Copy `deploy/` into the
simulation and keep the root files as the annotated reference. See
`deploy/README.md`. After editing either set, `tools/deployEquivalence.m`
re-runs the harness scenarios with both and reports whether the telemetry
is still bitwise identical.

## Harness / offline design (not deployed)

`run_recovery.m` (flat 2 Hz demo + offline design), `run_scheduler_test.m` +
`scheduler_sim.m` (100 Hz + 2 Hz split-rate architecture test),
`generate_test_data.m`, `build_Phi_Q.m`, `kalman_update_step.m`,
`kujur_params.m`, `solve_N_min.m`, `compute_PMD_eq38.m`,
`recovery_nav_sim.m` (a host following `docs/HOST_2HZ_WIRING.m`),
`ins_100hz_template.m` (what the host's 100 Hz function must add).

## Tests

Run `runAllTests` from the repo root before integrating (MATLAB; for Octave see
`tools/octave_shim/README.md`). Do not add `tools/octave_shim` to a MATLAB path. It runs the four scripts below in
order and prints one PASS/FAIL line each.

- `run_scheduler_test.m` — 100 Hz + 2 Hz architecture (T1-T6)
- `tests/test_monitors.m` — engine unit + Monte-Carlo checks (replaces the
  dual_monitor-era `test_all.m`, which cannot drive `SPF_cpiMonitor` since the
  window length / thresholds moved into `CST_spfParam`)
- `tests/test_variable_numMeas.m` — varying satellite count, padded vs exact
- `tests/test_feedback_split.m` — host setKF feedback split (GAIN), drain in COAST
- `tests/test_error_handlers.m` — unusable covariances via matrixInv, numMeas clamp, NaN input, anchor stamp
- `tests/test_matrixInv.m` — contract of the host's SVD inverse: zero padding, accuracy on code+carrier S, non-finite flag, rank deficiency; run it with the simulation's own `matrixInv` on the path
- `tests/test_shadow_mode.m` — guide Step 5: shadow mode on clean data (zero alarms),
  then a scripted attack with authority (latch, vetoes, one commit)
- `run_recovery.m` — end-to-end demo with offline design and plots

## Host decisions

`docs/HOST_DECISIONS.md` tabulates, per mode and event, what `setKF` sets
(`stateFB`, `KF.states`, `KF.covariance`), how the trial is handled, what
each mode compares against, and the gate outputs. `docs/HOST_2HZ_WIRING.m`
is the code-level contract.

## Logging and diagnostics

`tools/SPF_logInit.m`, `tools/SPF_logAppend.m` and `tools/SPF_plotLog.m`
(MATLAB desktop, not code generation) record the gate telemetry once per
2 Hz epoch and draw the standard figures: mode timeline with events,
SS/CPI margins per axis (alarm at 1), protection level, host NIS,
re-validation statistic against its threshold, dwell, coast time, faults.
The header of `SPF_plotLog.m` says how to read them for a jamming or
false-latch investigation. The mode figure shades the warm-up: after a
gate initialisation the monitors stay disarmed (`info.armed = false`, no
latch possible) until the filter has had `ARM_EPOCHS` consecutive epochs
with enough rows and a protection level below `ARM_PL_MAX`
(`CST_spfParam`).

`tools/reacqChecks.m` is the re-acquisition diagnostic: given four `.mat`
captures of the host's `kfUpdate` (last epoch before an outage, first call
back, first accepted update, a healthy epoch) it reproduces the host update,
compares h(x) with H, tests the innovation and the prior against the truth,
splits the update by row type, and traces the prior back through the
attitude residual. The header lists the variables each capture must hold.

`tools/runCompare.m` decomposes an INS-only outage drift from the logs of one
or two full runs, with nothing re-run: it estimates the truth time-tag lag,
gives the aided accuracy against the truth, the velocity, tilt and bias error
at the last aided epoch with their sigmas, fits the outage error trace as
a + b t + c t^2 + d t^3 to separate velocity, bias/tilt and gyro drift, and
predicts the drift from the state errors. With two runs it prints their
difference at the outage start in sigmas, and the aided position and
velocity error rms per aided window with the ratio between the runs (a
ratio that is the same in every window is a filter-wide change, one
concentrated in a window or two is an episode inside them). Its time-alignment section first
checks the truth arrays against the differenced truth position (a velocity
logged at another point of the host step, in another frame or missing, and a
yaw with another convention show there) and the nav velocity against the
differenced nav position, then finds the lag of the nav solution behind the
truth separately on the position, the velocity and the yaw. The attitude lag
is the reference, because a gyro-integrated attitude cannot trail the truth;
what the position and velocity trail beyond it is error of the nav solution,
the signature of GNSS measurements applied late, and the regression of the
errors on the truth velocity and acceleration quantifies it, and the position
error in the body frame regressed on the speed separates a reference-point
offset (a constant along or across the heading) from a latency (a speed
term). The slope of the
position error over the first 10 s of the outage is printed as the model-free
value of the velocity error at the outage start, and the budget includes the
yaw error acting on the manoeuvres of the outage. A clean run without an
outage is analysed with tOut and tReacq empty: the alignment and aided-accuracy
sections then run over the whole run and the outage sections are skipped. The
header lists the arrays each run struct must hold (nav, truth, bias estimates,
optional P diagonal).
`tools/timeOffset.m` is the filter-free companion: the time shift between
two position records of the same vehicle, for instance the receiver's own
position output against the reference, which says whether the GNSS tags run
late against the reference before anything in the host is changed.
`tools/latencyFromInnov.m` measures the same latency from inside the 2 Hz
function: per epoch it fits the pseudorange innovations to a clock term plus
dlt times the line-of-sight velocity, so observables one epoch older than the
state they are compared with show up as dlt of 0.5 s, and a split by epoch
type tells which of the two 2 Hz epochs is paired wrongly.
`tools/pairingCheck.m` needs no filter quantity at all: from the measured
pseudoranges, the satellite positions the host uses in h(x) and the
reference position, it finds the shift of the reference instant at which the
ranges fit, so observables one epoch older than the pulse they are paired
with show a minimum at −0.5 s. It takes the host's per-epoch log as logged
(prLog / satPosLog / sysIdLog, one block of slots per constellation) through
`tools/pairingLogUnpack.m`, and prints the receiver clock it removed per
system with a warning when the residual says the ranges or the positions are
not the ones of h(x). `tests/test_pairingCheck.m` covers the log format.
`tools/dopplerCheck.m` reads the spfMeas records saved per 2 Hz epoch (the
struct the host already builds for the gate, which carries y, S, H, R and the
state before and after the update) and says what the range-rate rows do in
the update: per row type the measurement sigma, what the filter expects from
its own state, the innovation rms, the gain fraction applied, the NIS, the
velocity correction per epoch, and whether the velocity partial of the
range-rate rows shares the sign and frame of the position partial of the range
rows. Verdict lines name an R far above the state uncertainty (rows switched
off in effect), a range rate built from the nav solution itself, and an
opposite-sign or wrong-frame derivative. `tests/test_dopplerCheck.m` covers
the five cases.
`tools/feedbackCheck.m` asks whether the KF update reaches the mechanisation,
when, and in full. Per nav sample the change of the logged position minus the
integral of the logged velocity is the part of the motion the velocity does
not explain; summed over the half epoch either side of tKf it is the
correction the update put into the position, whenever and however it was
applied, and its running sum regressed on the increment postState minus
priorState across the epochs is the step profile, from which the instant of
the update is read (10, 50 and 90 % of the step). The change is regressed
jointly on the increments of the records k-1, k and k+1 (-1 at k an increment
applied in full at its own epoch, 0 one that never reaches the logged
solution, a coefficient at k-1 or k+1 records offset by one epoch) and on the
state x+ of record k (-GAIN when the host applies GAIN x+). The velocity is
read the same way against lines fitted either side of the step. From the
records alone, x- of record k+1 regressed on x+ of record k gives the share of
each state the bookkeeping keeps after the feedback (1 - GAIN for the host's
split), the velocity residual integrated into the position residual (0.5 s
times the kept share when the residual is propagated with Phi) and the
position residual leaking into the velocity residual (Phi transposed). With
biasEst the sum of the accelerometer bias increments over the window is set
against the change of biasEst. Verdict lines name an axis whose increment
does not reach the solution, one applied in part, over more than one epoch or
with the opposite sign, offset records, a slewed step, a logged solution that
receives less of the state than the bookkeeping removes (part of each
correction lost between setKF and the mechanisation), and a residual
propagated with Phi transposed or without its coupling.
`tests/test_feedbackCheck.m` covers the applied, half-applied, unapplied,
offset and slewed cases, the host's split with GAIN 0.4, and a split that
books 0.4 but applies 0.16 with Phi transposed.
`tools/covPropCheck.m` reads the same spfMeas records and asks whether the
host's P keeps its cross terms from one 2 Hz epoch to the next: the previous
P+ is propagated over the interval with the INS error dynamics (no Q), the
epoch's own rows are applied, and the result is compared with the logged P+
block by block (pos-vel, vel-tilt, vel-accel bias, tilt-gyro bias). A block
that the propagation rebuilds and the stored P+ has lost is removed somewhere
in the host's 2 Hz path, which no Q can repair; a chain from tOut shows what
a P with every cross term would carry through the outage next to the logged
one, and the pDiag given to runCompare is checked to be diag(postCov). Given
the propTel the gate received per epoch as well, it reads the interval the
accumulator covers from accumPhi(1,4) (one epoch when the reset follows the
gate call, zero when it precedes it, more when a reset is missed), compares
the coupling blocks of accumPhi with the error dynamics, and repeats the
one-step prediction with the host's own accumPhi and accumQ.

`tools/spoofProfile.m` reads the spfMeas records of a spoofed run and says
what the spoof did to the measurements when its profile is not known: per
epoch it splits the range rows where a gap in their innovations is wider
than the group spread, takes the group farther from zero as the offset
group, and prints onset, end, size and rows of the offset and its course in
time; with the system of each row it names the spoofed constellation. The
offset is read whole while the gate coasts and as a lower bound while the
filter follows a slow spoof. The split needs rows of both kinds in one
record, so the tool also prints the mean innovation per system epoch by
epoch over the stretch where some system is off, with the host's row set
next to it, and with the pairing log and the reference position it reads
the profile with no filter in the way: each range minus the geometric
range from the reference leaves the receiver clock, read where two systems
agree and interpolated elsewhere, the inter-system biases as run medians,
and the offset per system, with its onset, end, size and course. A spoof
that moves the position gives every satellite its own offset, so the tool
also fits, per epoch, the position and clock the off system's rows ask for
(in the innovations, relative to the host's prediction) and the position
and clock its ranges describe (against the reference), per system, with
the fix in the minute before and after the stretch (a sky that is clean
again reads near zero after), and prints a timeline of the record's row
set around the spoof with the log's used slots and the gate mode under it
when the gate log is given. Off epochs closer than 10 s form one stretch;
every stretch is listed, the tables and the timeline follow the strongest,
and the timeline says whether the log holds ranges the host did not use or
flagged with corrValid false. `docs/SPOOFPROFILE_PRINTOUT.md` walks through
the printout block by block. The Doppler rows get the same reading: the
velocity and clock drift the range-rate rows ask for, and against a
stationary reference the log's range rates minus the line-of-sight
satellite velocity, so a spoof that moves the position without the
matching Doppler shows up as a velocity in the ranges and none in the
rates.

`tools/spoofInject.m` puts a spoof into one constellation's slot arrays of
the host for one epoch, the ranges and range rates a receiver standing at
the true position plus a position offset, moving with a velocity offset
and carrying a clock offset and drift would measure, exact in the
geometry; unused slots are left alone. `tools/spoofScenario.m` gives the
offsets at time t for a position step, a position ramp with its matching
velocity, a clock step or a clock ramp between tOn and tOff, and zero
outside. Together they run the matrix on the clean twin: the host calls
spoofScenario once per epoch and spoofInject once per constellation to
spoof, after the corrected ranges are formed and before the rows are
built, and spoofProfile reads the injection back from the ranges.

`runCompare` ends the outage at the first logged P that falls by more than
30% inside the given interval (the first update back, or a reset): the
drift, the trace fit and the budget stop there, so a tReacq given too late
does not pull post-hand-back samples into the INS-only fit.

`tools/designConstants.m` is the offline design step (quickstart step 2) on
recorded data: from the same spfMeas records it takes Eq. 20 per monitored
axis on every usable epoch, keeps as the design sky the 5th percentile over
the epochs the gate would be armed on (enough rows and a protection level
inside `ARM_PL_MAX`, from the P diagonal; the whole-run minimum sits on the
first epoch after an outage and is printed for information only), solves the
window length per axis where the missed-detection requirement is reachable,
and prints the CPI threshold, SS gate and PL multiplier at the window used
next to the values in `CST_spfParam`. When the requirement is out of reach,
as it is for a code-based filter against a 0.1 m spoofer, the current window
is kept, the constants come out equal to the published ones, and the tool
prints what the data does set: the constant offset per axis the CPI sees
within one window on the median sky, at the threshold and at the
missed-detection probability, `NOMINAL_BIAS` included (given as an argument
or read from the `CST_spfParam` on the path). With the P diagonal it also
reports the epoch from which the protection level allows arming. S carries
the host's P and R, so it is rerun after every change of Q or R.

## Printable copies

`docs/pdf/` holds PDF renderings of this README, `docs/HOST_DECISIONS.md`,
`docs/REVIEW_FINDINGS.md` and `docs/CONCEPTS.md` (separation, coast, P_C /
P_KF, anchor and warm-up explained in question-and-answer form), plus
`ALL_DOCS.pdf` with all of them in one file.
Regenerate after editing the Markdown with `docsToPdf` (a MATLAB
function in `tools/`): it writes the HTML files itself and calls Edge or
Chrome headless for the PDF step; without a browser the HTML files can be
printed from any browser.

## Code generation

The runtime path (`SPF_gate` and everything it calls) is written for
MATLAB Coder: fixed-size padded buffers, one struct layout per type from
`STRUCT_SPF`, `uint8` enumeration modes, no cells, varargin, dynamic
fields, function handles or try/catch. `tools/codegenSpoofMonitor.m`
(MATLAB only) runs `coder.screener` on the call tree and builds a C library
with the entry-point types taken from the `STRUCT_SPF` zero constructors.
The gate must only run in navigation mode (see the alignment section of
`docs/HOST_2HZ_WIRING.m`).

## Regression

`reference/` holds the console output of both scripts from the author's
MATLAB run; the current engine reproduces it exactly on MATLAB. `tools/octave_shim/`
lets everything run under GNU Octave (event epochs differ by the RNG stream).
