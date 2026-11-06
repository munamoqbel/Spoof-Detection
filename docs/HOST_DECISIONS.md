# Host integration decisions

Companion to `HOST_2HZ_WIRING.m`. Records what the host's 2 Hz function does
with the operational `KF`, the `trialKF` and the gate outputs in every case,
and what each gate mode compares against. `GAIN` is the host's constant
feedback gain (0.4). "Mode" is the gate mode **after** this epoch's call to
`SPF_gate`; "in probation" is the mode **before** it.

## 1. Definitions

| Symbol | Meaning |
|---|---|
| `xPrior` | `activeKF.states` captured before calling `kfUpdate` (extrapolated at 100 Hz) |
| `xPost` | `kfUpdate` output: `activeKF.states + K*y` |
| `PPrior`, `PPost` | `activeKF.covariance` before, `kfUpdate` covariance after |
| `y`, `H`, `R`, `numMeas` | innovation `z - h(xPrior)`, Jacobian, the NOMINAL measurement noise, all the rows the host built (pressure row included; `numMeas = 1` on a GNSS outage with the baro row alone, 0 with none). A host that scales R on large innovations or rejects rows for its own update passes those beside the nominal pair (`setKfMeasBoth`, `kfMeasFromUpdateBoth`): the monitors read the nominal pair only (section 8) |
| increment | `xPost - xPrior = K*y`: the only thing the monitors consume |
| separation | sum of Phi-propagated increments since a reference epoch = host solution minus the coast of that reference |
| anchor | last alarm-free closed window: `separation` (kept live) and coast covariance |
| `probSep` | sum of the trial's increments since probation opened |
| `nav.state`, `nav.covar` | a complete `(x+, P+)` for the operational KF, same definition as the `kfUpdate` outputs |
| active filter | operational `KF` in NOMINAL and COAST, `trialKF` in PROBATION |

The gate never stores an absolute state. Every quantity it keeps is a
difference between two solutions that share the same mechanization history,
which is what makes it invariant to how `setKF` splits `xPost` between
`stateFB` and `KF.states` (verified: `tests/test_feedback_split.m`, identical
estimate for `GAIN` 0.1, 0.4 and 1.0).

## 2. Operational KF: what `setKF` sets in every case

Exactly one `setKF` call per epoch on the operational KF. `stateFB` rows are
position/velocity/attitude; bias and clock feedback follow the existing
cumulative rule in every case.

| Case | Gate flags | `stateFB` | `KF.states` | `KF.covariance` | Scratch `xPost`, `PPost` |
|---|---|---|---|---|---|
| NOMINAL, normal update | none | `GAIN * xPost` | `xPost - GAIN * xPost` | `PPost` | used |
| NOMINAL, no satellite | none | zero except biases | `xPost` except pressure-alt `(1-GAIN)*xPost` | `PPost` | used (existing rule) |
| NOMINAL, invalid update | none | zero | zero | as today | to the gate: `numMeas = 0`, `xPost = xPrior`, `PPost = PPrior` |
| LATCH (NOMINAL -> COAST) | `nav.applyCorrection`, `info.eventLatched` | `GAIN * nav.state` | `nav.state - GAIN * nav.state` | `nav.covar` (anchor coast covariance) | discarded |
| COAST (stays) | none | `GAIN * states_extrap` (drain) | `states_extrap - GAIN * states_extrap` | input `KF.covariance` unchanged (extrapolated) | discarded; only `y, H, R, numMeas` feed the gate |
| COAST -> PROBATION | `kfCommand.startTrial`, `info.eventProbationStarted` | drain, as COAST | drain, as COAST | input unchanged | discarded; next epoch `trialKF = KF` (states and covariance, already extrapolated) |
| PROBATION (stays) | none | drain, as COAST | drain, as COAST | input unchanged | belongs to the trial (see section 3) |
| VETO (PROBATION -> COAST) | `info.eventProbationVetoed` | drain, as COAST | drain, as COAST | input unchanged | trial discarded |
| COMMIT (PROBATION -> NOMINAL) | `nav.applyCorrection`, `info.eventHandback` | `GAIN * nav.state` | `nav.state - GAIN * nav.state` | `nav.covar` (trial `PPost`) | `nav.state` is the trial's `xPost` |
| No usable anchor at alarm | `info.anchorMissing`, no correction | as COAST | as COAST | input unchanged | discarded |
| Alignment (`navActive = false`) | gate called, returns `zeroTel` and re-arms itself | as today | as today | as today | as today |
| GNSS outage in navigation (pressure row only) | gate called with `numMeas = 1`; cannot pass re-validation (`REVAL_MIN_MEAS = 4`) | as today (no-satellite case) | as today | as today | used |

Implementation: `setKF` itself is unchanged. The host selects its input
before the one call per epoch: `nav.state` / `nav.covar` on latch and commit,
the extrapolated `KF.states` / `KF.covariance` on COAST and PROBATION (the
normal case on the prior is the drain), `kfPost` otherwise. On the two
overriding branches the fields that select the case (`failed`, satellite
count) are forced to the normal case, so a failed or no-satellite scratch
update can neither zero the clean state nor trigger the 100 Hz reset.

Rules behind the table:

- Whatever is fed to the mechanization is removed from `KF.states` in the same
  call. The drain applies the normal split to the extrapolated state, so the
  leftover of the latch correction reaches the output smoothly, by `(1-GAIN)`
  per epoch, instead of sitting in `KF.states` for the whole coast.
- The covariance in COAST and PROBATION is never the scratch `PPost` (it was
  reduced by a measurement being refused) and never the gate's output (only
  meaningful on latch and commit epochs; equal to the extrapolated covariance
  otherwise).
- The coast is pure INS for now: no baro update in COAST or PROBATION. A
  baro-only aiding loop during the coast is a later gate extension.
- The latch correction is bounded by the protection level, since the monitors
  fire before the spoofer gets further; it is not the size of the attack.

## 3. Trial KF (PROBATION only)

| Step | Action |
|---|---|
| Probation opens (`startTrial`, acted on at the start of the next epoch) | `trialKF = KF`, states and covariance as the 100 Hz side extrapolated them; no propagation on that epoch. The operational KF keeps draining meanwhile (negligible by then) |
| Start of each PROBATION epoch | `trialKF.states = accumPhi * trialKF.states`; `trialKF.covariance = accumPhi * P * accumPhi' + accumQ`, symmetrised (the trial has no 100 Hz extrapolation of its own) |
| Update | `kfUpdate(measurement, trialKF)`; `h(x)` must use the states of the struct passed in |
| After the gate call | `trialKF = kfUpdate result` (states and covariance). No `setKF`, no feedback, no drain: the trial has no mechanization of its own and its states must keep the full increments |
| Veto | trial dropped |
| Commit | trial `xPost`, `PPost` go through the normal `setKF` as `nav.state`, `nav.covar`; trial no longer used |

## 4. What each mode compares, and against what

| Mode | Active filter | Test | Compares | Against | Statistic and threshold | Result |
|---|---|---|---|---|---|---|
| NOMINAL, warm-up (`info.armed = false`) | operational KF, updated normally | none: the monitors arm after `ARM_EPOCHS` (10) consecutive epochs with `numMeas >= REVAL_MIN_MEAS` and `K_MISSED_DETECTION * sigma_pos(P+) < ARM_PL_MAX` (100 m); a PL above the limit restarts the count, a row dropout with the PL inside pauses it; sticky once armed (until a gate re-initialisation) | | | | no windows, no alarms, no latch; the anchor follows the current solution so it is fresh when the monitors start |
| NOMINAL | operational KF, updated normally | SS, every open window, every monitored axis, every epoch | window separation (increments since the window opened) | that window's INS-only coast | `abs(d_axis) > K_FALSE_ALERT * sqrt(P_C - P_KF) + NOMINAL_BIAS(axis)` (the allowance is 0 in the repository; the host sets the nominal GNSS excursion it tolerates) | any alarm: LATCH to the anchor |
| NOMINAL | same | CPI, when a window reaches N = 10 | normalised innovation projections on the axis beyond the allowance, `xi = max(abs(f'S^-1 y) / sqrt(f'S^-1 f) - NOMINAL_BIAS(axis) * sqrt(f'S^-1 f), 0)` | N(0,1) under no attack | `sum(xi^2) > CPI_THRESHOLD` (Gamma, 45.64) | any alarm: LATCH; a window alarm-free for its whole life refreshes the anchor |
| COAST | operational KF **not** updated (scratch update feeds `y, H, R` only) | re-validation, every epoch | raw innovation `y` (all rows) | the INS-only coast (the operational KF's own prior) with its grown covariance | `y' (H P_C H' + R)^-1 y < chi2inv(1 - 1e-3, numMeas)` and `numMeas >= REVAL_MIN_MEAS` (4); baro-only epochs cannot pass | 10 consecutive passes: PROBATION; a fail resets the dwell; an epoch with fewer than `REVAL_MIN_MEAS` rows (no GNSS, baro only, a repeated message the host dropped) holds it. Monitors do not run |
| PROBATION | trial KF, updated | SS and CPI on the trial, same as NOMINAL | trial increments (`probSep`) | the coast (`P_C - P_trial`) | as NOMINAL | any alarm: VETO to COAST (dwell reset); 18 quiet epochs: COMMIT |
| PROBATION | same | re-validation, diagnostic only | `y - H * (-Phi * probSep)` | the coast | same statistic, logged as `info.qReval` | no decision |

Why two stages: the COAST test is coarse and its gate widens as the coast
covariance grows, so a spoofer resuming with an offset inside that band passes
it. PROBATION lets a filter follow the GNSS and runs the full monitor bank on
it, with NOMINAL sensitivity, without touching the operational solution. None
of these tests reads `KF.states`; the residual shrinking during the coast is
the estimate having no information, not a sign that the attack ended.

## 5. Gate outputs per event

| Event | Flags | `nav.state` | `nav.covar` | `kfCommand` |
|---|---|---|---|---|
| Latch | `applyCorrection`, `eventLatched`, `eventAnchorEpoch` | operational `xPost - anchor.separation` | `anchor.covariance` | none |
| Probation opens | `eventProbationStarted` | none | none | `startTrial = true` |
| Veto | `eventProbationVetoed` | none | none | none |
| Commit | `applyCorrection`, `eventHandback` | trial `xPost` | trial `PPost` | none |
| Alarm without usable anchor | `eventLatched`, `anchorMissing` | none | none | none; mode still COAST |
| Every epoch | `mode`, `ssAlarm`, `cpiAlarm`, `alarmPerAxis`, `maxProtectionLevel`, `qReval`, `dwellCount`, `coastEpochs`, `ssRatio` (SS margin per axis, alarm > 1), `cpiRatio` (CPI margin per axis, alarm > 1), `armed` (monitors active; false during the warm-up after init), `resetInhibit` (gate in COAST / PROBATION or alarming: the host must not apply its protective reset on the GNSS-vs-INS discrepancy), `coastBudgetExceeded` (coast PL `K_MD * sigma_C` above `ARM_PL_MAX`, or `COAST_BUDGET_EPOCHS` passed: the host may reset, integrity not assured) | | `coastCov` (diagnostic) | |
| Fault telemetry | `inputFault` (non-finite input: epoch dropped, gate re-initialises on the next good epoch), `numMeasClamped` (host passed more than `MAX_MEAS` rows), `solveFault` (an S or residual covariance was not positive definite: that CPI epoch counts `xi = 0`, re-validation cannot pass) | | | |

After a commit the anchor is the committed solution (separation zero,
covariance the trial's `PPost`), the probation windows carry over, and a new
attack goes through the same cycle with no limit on the number of cycles.

## 6. Exception handlers (audit of the runtime path)

| Hazard | Where | Handling |
|---|---|---|
| `S^-1 y`, `S^-1 f` on the padded S | `SPF_cpiMonitor` | live-row variances must be positive; padding diagonal filled with the largest live variance so the matrix is non-singular; then the host's `matrixInv` (flags non-finite / singular input); a flagged epoch is dropped (`xi = 0`), `solveFault` |
| `r' S_r^-1 r` on an unusable residual covariance | `SPF_revalidation` | live-row variances must be positive, then the host's `matrixInv`; `q` must be finite and `>= 0` to pass; flagged: no pass, `solveFault` |
| `numMeas > MAX_MEAS` slicing fixed buffers | `setKfMeas`, `kfMeasFromUpdate`, `SPF_cpiMonitor`, `SPF_revalidation` | clamped at ingress and defensively at use; `numMeasClamped` |
| NaN / Inf in any input (failed host update, uninitialised `propTel`) | `SPF_gate` | `isfinite` check on every input; epoch dropped with `inputFault`, state and epoch counter reset, re-init on the next good epoch (mirrors the host's own failed-update reset) |
| `sqrt` of a negative variance (`P_C - P_KF`, `P_C`, `coastCov` diagonal) | `SPF_ssMonitor`, `SPF_protectedNav` | guarded: test reported undefined (no alarm), sigma 0 |
| division `gamma / sqrt(sigma2)` with `sigma2 = 0` (axis unobservable) | `SPF_cpiMonitor` | guarded: `xi = 0` |
| startup anchor stamped with epoch 0 (the "no anchor" sentinel) and propagated once too often | `SPF_protectedNav` | during the warm-up the anchor is re-seeded every epoch with that epoch's `(x+, P+)` and stamp, so the first armed epoch starts from a fresh anchor |
| false latch on an unconverged filter right after (re)initialisation (few rows, large P, interval matrices just reset) | `SPF_protectedNav` | warm-up arming: no window and no latch until `ARM_EPOCHS` consecutive epochs qualify (`numMeas >= REVAL_MIN_MEAS`, PL < `ARM_PL_MAX`); `info.armed` |
| `uint32` epoch differences, `uint8` counters and loop variables | all | checked: no wrap possible (`anchor.epoch <= epoch`), classes consistent |
| `REVAL_THRESHOLD_TABLE(numMeas)` index | `SPF_revalidation` | `numMeas` clamped to `MAX_MEAS` (= `CST_gnssHybrid.MAX_MEASURES`), table has 60 entries and the index is clamped to its length |

Covered by `tests/test_error_handlers.m`. On a host update flagged `failed`, pass the epoch as `numMeas = 0`, `xPost = xPrior`, `PPost = PPrior`, or set `resetRequest = true` on that call; a NaN slipping through is caught by the input check.

## 7. Deferred

- Baro-only aiding during COAST and PROBATION (vertical channel).
- Coast budget (Implementation Guide rule 4): reported only
  (`K_MD * sigma_C > ARM_PL_MAX`, with `COAST_BUDGET_EPOCHS` as a time belt,
  `info.coastBudgetExceeded`); the host decides what
  to do beyond it (its protective reset). The gate does not declare dead
  reckoning itself. Each cycle still costs a coast and widens the coarse
  re-validation gate.
- `MAX_ANCHOR_AGE` is a startup-only guard today.
- Warm-up arming: `ARM_EPOCHS` qualifying epochs (rows, a PL below
  `ARM_PL_MAX`) on which the monitors, run in shadow, would not have
  alarmed (section 10). `ARM_EPOCHS` is 10 in the repository (simulation
  references); the host runs with 120 (60 s). A row dropout pauses the
  count; a PL above `ARM_PL_MAX` or a shadow alarm restarts it.
- Outage re-arm (disarm after a long GNSS outage so re-acquisition goes
  through the warm-up): not needed on the jammed data set, where the host's
  own reset already re-initialises the gate; relevant only for outages that
  end below the host's reset limit.

## 8. The measurement noise the gate reads

The CPI normalises each epoch's innovation by S = H P H' + R and the
re-validation compares the GNSS rows with the coast through H P_C H' + R.
Both are tests against the null hypothesis, so their R is the noise the
design assumes under no spoof, not the weight the host chose for its own
update. A host that scales R up on large innovations normalises a spoof
innovation to about one sigma before the gate sees it, and a host that
drops a satellite whose range or range rate jumps between epochs (this
host: for 30 epochs) removes, at a step onset, the rows that see the step
most, for 15 s, and the same rows again for 15 s at the step back. Measured
on the vehicle data (`tools/cpiReplay`, `tools/matrixRow`): a 10 m step on
every constellation read 0.09 of the CPI threshold with the host's R and
0.30 with the nominal one (no latch either way, as designed); a 20 m step
read 0.85 with the host's R and crossed the threshold at its second epoch
with the nominal one; a 30 m step latched after 4.8 s with the host's R,
and in COAST the re-validation passed on 48 of 48 epochs with q at 0.42 of
its threshold on a nearly full row set, handing the nav back into the spoof
24 s after the latch; with the nominal R the same residuals give q near
twice the threshold. Decision: the gate's `innovationCov` and
`measNoiseCov` carry the nominal R on all the rows the host built, taken
before the host's spike handling; the host's own S and R travel beside them
(`innovationCovHost`, `measNoiseCovHost`) for the record. The host's update
keeps its scaling and its rejections. With the nominal R the gate also sees
the receiver's own events whole; on 55 minutes of the clean vehicle run the
CPI then needs an epoch floor of 8 range rows (reflections in a partial
outage read as tens of metres on a one-axis projection) and reads two
recording slips of one epoch as 16 m jumps; the real near misses stay
below 0.65 of the threshold (`tools/cpiReplay` with the mode `nominal`).

## 9. Re-admission bounded by COAST_SIGMA_MAX

Run 5 of the vehicle matrix (30 m step on every constellation, the gate
handed S with the nominal R): latch at 600.78 s, 0.8 s after onset, the
coast pure (the innovation projection stayed at 31 to 35 m), and yet the
re-validation passed on all rows at 644.87 s and the probation committed
at 658.01 s, 42 s before the spoof ended. The chi-square compares the rows
with the coast through H P_C H' + R, and the coast sigma had grown from
1.0 m at the latch to about 6 m; 30 m then reads as 1.5 sigma per row. The
probation did not veto because the trial, starting from that covariance,
absorbed 85 % of the offset in its first update, before any window was
open, and sat still afterwards. In the paper's terms the offset was inside
the protection level, about 60 m at that moment, and the recovery admits
an error within the PL. For a gate whose purpose is to keep the nav off a
spoof that is not acceptable, so the re-admission tests use the coast
covariance bounded by `COAST_SIGMA_MAX` (3 m, `STRUCT_SPF.boundCoastCov`,
a congruence scaling that keeps the matrix positive semi-definite) on the
`COAST_BOUND_STATES`, the NED position and the receiver clock bias: every
range row carries the clock, so an offset common to all rows is absorbed
by a clock variance that grows freely through the coast (run 4 of the
matrix, a 20 m D step, with the position bounded alone: the chi-square
passed on 60 % of the coast epochs and the coast CPI fell from 6 to 2
over the 100 s, the spoof held only by the CPI against the coast), and
the first probation window opens on the coast with that covariance so the
trial's first update counts. With 5 m range rows the re-validation then
refuses any offset above about 24 m for as long as it lasts, the probation
vetoes a trial that jumps by more than k_FA sigma_SS + B, about 17 m, and
clean GNSS is still re-admitted after a real coast whose drift stays below
those figures, two to three minutes on this INS. Beyond that the gate
stays in COAST until the coast budget (unbounded P_C, as before) and the
host's reset, which is what the host does after a jam today. Offsets
below about 17 m would still be admitted by the chi-square, which is one
epoch at a time at 1e-3 on all rows, while the CPI latches from about
20 m. So the dwell also requires the CPI statistic of the innovations
against the coast over the last N epochs (the pool's windows, S with the
bounded covariance, a zero increment so the SS stays silent, B and T_N as
in NOMINAL) to stay below T_N: GNSS the monitors would latch on is not
re-admitted. With 5 m range rows and B = 5 m this refuses offsets from
about 14 m, below what the gate detects, and admits clean GNSS after a
real coast whose drift stays below about 13 m, some 100 s on this INS;
beyond that the gate stays in COAST until the budget and the host's reset.
The CPI margin is logged through the coast (cpiRatio), so the plots show
how far the coast is from re-admission.

## 10. Arming on the monitors' silence

On every vehicle run the gate latched at 83.1 s, 7 s after it armed, on
the start transient of 82.6 to 84.6 s (both constellations 12 m off for
five epochs while the heading was still 4 degrees off; it converges at the
first manoeuvre near 200 s). The host's alignment flag is true from the
start and cannot carry the heading, so `navActive` cannot hold the gate
off, and a longer `ARM_EPOCHS` would be tied to this trajectory. The
monitors therefore run in shadow through the warm-up: windows open and
close as in NOMINAL, the margins are logged (the grey part of the plots),
nothing latches and the alarm flags stay down, and a shadow alarm restarts
the arming count. The gate arms only after `ARM_EPOCHS` epochs on which it
would not have alarmed, whatever the cause of a transient, and the pool
carries over into NOMINAL so the first windows after the arming are
already old. On this data the 83 s transient restarts the count and the
gate arms near 143 s; the manoeuvre at 200 s alarmed in no run.

## 11. The host's spike filter and the gate

The host drops a satellite for 30 epochs (15 s) when its range or range
rate jumps between consecutive epochs beyond a limit. It runs on the raw
measurements, before the innovation, in every mode. For the gate this
means: a ramp never trips it, so the gate sees every row during a ramp;
at a step onset the satellites whose line of sight projects the step most
are gone for 15 s and come back carrying the spoof, which the gate then
reads whole (run 5: 6 range rows at the first spoofed epoch, 14 again
from 618 s); at the step back the same happens to the clean rows, so a
hand-back comes about 15 to 20 s after the spoof ends (13.6 and 18.5 s
measured). The one risk is inside those 15 s in COAST: the chi-square
runs on the rows that saw the step least and can pass on a spoofed epoch
(one isolated pass at 604.30 s in three runs); the dwell of ten
consecutive passes and the CPI against the coast cover it. The gate takes
the rows as the host gives them; passing the dropped rows as well is not
needed.
