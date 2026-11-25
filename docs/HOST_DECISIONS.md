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
| Return at a latch | `applyCorrection`, `eventLatched`, `eventReturn`, `eventAnchorEpoch` (the long coast's start) | operational `xPost - longCoast.separation` | the long coast's covariance | none |
| Return in COAST | `applyCorrection`, `eventReturn`, `eventAnchorEpoch` (the long coast's start) | coasting `xPrior - longCoast.separation` | the long coast's covariance | none; the mode stays COAST, the dwell restarts |
| Probation opens | `eventProbationStarted` | none | none | `startTrial = true` |
| Veto | `eventProbationVetoed` | none | none | none |
| Commit | `applyCorrection`, `eventHandback`, `handbackUncertified` (the probation opened past the hold limit `COAST_BUDGET_EPOCHS`: re-admitted on the coast's own covariance, integrity not assured across it) | trial `xPost` | trial `PPost` | none |
| Alarm without usable anchor | `eventLatched`, `anchorMissing` | none | none | none; mode still COAST |
| Every epoch | `mode`, `ssAlarm`, `cpiAlarm`, `alarmPerAxis`, `maxProtectionLevel`, `qReval`, `dwellCount`, `coastEpochs`, `ssRatio` (SS margin per axis, alarm > 1), `ssSeparation` and `ssSigma` (the separation `d` and its `sigma_SS` per axis of the window with the largest SS margin; the margin is `|d| / (K_FALSE_ALERT sigma_SS + B)`, so the largest clean `|d|` of a run is the `B` the host sets), `cpiRatio` (CPI margin per axis, alarm > 1), `longSeparation`, `longSigma` and `longAge` (the oldest long INS coast: nav minus coast on N and E, its sigma, its age in epochs; while the GNSS is trusted the separation is the INS's own wander over the last 80 to 120 s as the GNSS corrected it), `armed` (monitors active; false during the warm-up after init), `resetInhibit` (gate in COAST / PROBATION or alarming: the host must not apply its protective reset on the GNSS-vs-INS discrepancy), `coastBudgetExceeded` (coast PL `K_MD * sigma_C` above `ARM_PL_MAX`, or `COAST_BUDGET_EPOCHS` passed: the host may reset, integrity not assured) | | `coastCov` (diagnostic) | |
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
- Coast budget (Implementation Guide rule 4): `K_MD * sigma_C > ARM_PL_MAX`,
  with `COAST_BUDGET_EPOCHS` as a time belt, `info.coastBudgetExceeded`,
  reported only; the host may reset on it. Past the time belt itself the
  gate re-admits as the paper does, on the coast's own covariance
  (section 9), and flags the hand-back `handbackUncertified`. A PL-driven
  re-admission was measured and withdrawn: on run-wide floors the PL
  budget comes at 50 s against 3 m of real drift and would hand a 100 s
  step spoof to the spoofer at 70 s. The gate does not declare dead
  reckoning itself.
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
over the 100 s, the spoof held only by the CPI against the coast). With
the clock bias bounded as well the run was repeated and the chi-square
still passed on 57 % of the coast epochs: a D offset is nearly common to
all rows, and a single-direction signal of about five sigma is diluted
over the 25 degrees of freedom of the chi-square, which then sits at its
threshold. The chi-square is the paper's test and stays as a condition of
the dwell; what holds a common-mode offset is the CPI against the coast,
the directional test, at 2 to 7 times T_N on that run. The clock bound
stays for a host whose clock variance grows through the coast (the
synthetic case D of test_coastBound), and
the first probation window opens on the coast with that covariance so the
trial's first update counts. With 5 m range rows the re-validation then
refuses any offset above about 24 m for as long as it lasts, the probation
vetoes a trial that jumps by more than k_FA sigma_SS + B, about 17 m, and
clean GNSS is still re-admitted after a real coast whose drift stays below
those figures, two to three minutes on this INS. Beyond that the gate
stays in COAST until the hold limit `COAST_BUDGET_EPOCHS`; past it the
re-validation, the CPI hold and the probation's windows run on the
unbounded P_C, the paper's recovery, so the GNSS is re-admitted through
the normal dwell and probation and the hand-back is flagged
`handbackUncertified` (measured on the long run: a followed 0.3 m/s
drag whose 30 m step back latched left the clean GNSS refused to the
end, the nav INS-only and 1.5 km off after 15 minutes, the host's reset
bound to a 5 km gap; the limit bounds that loss to the drift over its
length). Offsets
below about 17 m would still be admitted by the chi-square, which is one
epoch at a time at 1e-3 on all rows, while the CPI latches from about
20 m. So the dwell also requires the CPI statistic of the innovations
against the coast over the last N epochs (the pool's windows, S with the
bounded covariance, a zero increment so the SS stays silent, B and T_N as
in NOMINAL) to stay below T_N: GNSS the monitors would latch on is not
re-admitted. With 5 m range rows and B = 5 m this refuses offsets from
about 14 m, below what the gate detects, and admits clean GNSS after a
real coast whose drift stays below about 13 m, some 100 s on this INS;
beyond that the gate stays in COAST until the hold limit, then re-admits as above.
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

## 12. The return path

Run 5 of the vehicle matrix, a 0.3 m/s east ramp for 100 s: followed, as
every ramp below about 3.6 m/s is on this host (1.5 m per 5 s window
against a gate near 18 m, an innovation offset of the rate times the 3 s
filter lag under the 5 m deadband), 30 m by the end. At the step back the
gate latched, 2 s later on the rows the spike filter had left, and the
bound refused the clean GNSS for the rest of the run: the coast, which
had inherited the ramp's velocity, moved away at 0.3 m/s, 57 m by the end.
"The GNSS is 30 m from the nav" reads both as a spoof that has jumped and
as a spoof that has ended; the paper re-admits both after about 60 s, the
bound refuses both. The inertial sensors tell them apart: in a followed
spoof the nav has wandered from where the INS alone would have carried it
over the last minutes, and the returning GNSS lands on that INS coast
within its drift (7 to 11 m over two minutes on this host); an onset lands
away from both the nav and the INS coast, and the nav had not wandered.

The gate therefore keeps `RETURN_BANK` long INS-only coasts beside the
anchor (host solution minus a coast, propagated as the anchor is, started
in turn every `RETURN_RESTART` epochs so that the oldest spans 80 to 120 s;
they start with the arming and again after a return and a hand-back). At
a latch, and on every COAST epoch with rows (the clean GNSS reaches the
gate only once the spike filter has released its rows, and a followed drag
can end while the gate already coasts on a later onset), it reads the
position the rows ask for, the CPI projection with the full S, against the
oldest long coast: a return when the nav before the event (the anchor at
a latch, so that a latch a few epochs late still reads right; the coasting
nav in COAST) had wandered from it horizontally by more than the yardstick
and the GNSS lands within the yardstick of it, the yardstick being that
coast's own horizontal protection level at its age, `K_MISSED_DETECTION *
sigma_long + B` from the long coast's covariance, never below
`RETURN_WANDER` / `RETURN_RADIUS` (`returnYardstick` in the telemetry). On a return the nav is corrected back to the long
coast, velocity included (`nav.applyCorrection`, in COAST relative to the
coasting KF), the coast's uncertainty becomes that coast's own and is the
floor of the bound (`boundCoastCov(P_C, sigmaLatch)`: no sigma above the
larger of `COAST_SIGMA_MAX` and the coast's sigma at the latch), so the
clean GNSS re-validates and the hand-back follows in about 15 s. The first
probation window keeps the short bound whatever opened the probation, so a
trial that jumps more than about 20 m is still vetoed. Otherwise nothing
changes: an onset latches on the anchor as before, and the long coasts
follow the nav through the correction.

The two floors, 15 m, sit below the smallest step the gate latches on
(about 20 m), so no onset qualifies through them; the yardstick above
them keeps the INS's own drift out. Fixed radii alone were not enough:
sized above a two-minute INS drift of 7 to 11 m, they met a vehicle whose
long coast drifted 15 m in 100 s, where a 20 m onset landed 15 m from that
coast while the nav had wandered 15 m from it, read as a return, loosened
the re-validation to the long coast's own uncertainty and handed the nav
to the spoofer (the matrix on the long run, step N 20 m). With the
yardstick the return path reaches only followed spoofs larger than the
INS's own drift uncertainty at that age, so its reach is the IMU's grade
through `sigma_long`; a smaller followed spoof that ends with a latch is
held on the clean GNSS until the hold limit, where the paper's recovery
re-admits it, flagged. What a spoofer can get
through the door is a hand-back within the yardstick of the long coast,
after a drag of at least the yardstick the spoofer could make anyway, and
the first window caps the probation's jump. Vertical returns are not
read (the GLONASS clock-height pattern reads 7 to 11 m on single epochs
in D); a ramp longer than the horizon has its last two minutes undone and
the rest must fit the loosened re-admission; an INS drift beyond 15 m in
two minutes makes the gate miss a return, which is the stranding of
before, not a new failure. `test_coastBound` cases E and G hold a
followed ramp's end and a drag-then-step's end on the synthetic host.

## 13. The gate's own coast covariance

Every gate threshold is k times a coast sigma: the separation gate, the
protection level and the coast budget, the re-validation's S, the long
coasts' covariance that floors the bound after a return. Until now the
gate propagated its coasts with the host's accumQ, and the host's Q is a
filter tuning: inflated to keep the filter responsive under the GNSS, it
overstated the drift of an INS-only coast about three times in sigma on
the vehicle host (section 9, `tools/coastDriftCheck`), which is why the
re-admission had to be bounded by hand and why the budget flag came at
68 s. The IMU model was then brought to the datasheet and its Allan plots
(`tools/imuModelCheck`, `tools/imuModelTerms`), and the host's Q to the
sensor on the attitude and bias blocks with its remaining terms named
(the position block for the receiver's 0.11 s offset, the horizontal
velocity for the heading's unobserved drift on a wheeled vehicle). What
the host cannot give up, the gate does not need: it builds its own Q.

With `COAST_Q_DATASHEET` the gate replaces accumQ, once per epoch and for
every coast it propagates, by `SPF_coastQ`: on the inertial states the
datasheet, the velocity from the VRW, the attitude from the ARW, the two
bias blocks from their in-run stability read as a walk over a coast plus
their bias walk, the position block without noise of its own, integrated
over the interval the host's accumPhi covers as dt/2 (Phi Qc Phi' + Qc);
on every other state, clock, pressure, GNSS states, the host's accumQ as
it comes. The host's propTel is unchanged on the wire; the host sets the
flag and the `IMU_*` and `COAST_*_STATES` constants in its CST_spfParam,
the states in its own layout. At two minutes the velocity block alone then
contributes 0.2 m of position sigma where the host's tuning gave 23 m,
and the coast sigma is the sensor's: about 5 m from the ARW plus what the
filter's errors at the latch grow into.

Those errors at the latch are the other half, and the host's P does not
state them honestly on every state: on a wheeled vehicle without heading
aiding the heading is unobserved on straight legs and drifts a degree in
a quarter of an hour with a P of 0.03 degree (`tools/navConsistency`,
the heading row at 0.1 with the mean out; the tilt and the horizontal
velocity at 0.4 through the turns). The gate therefore reads the host's P
through `COAST_SIGMA_FLOOR_STATES / _VALUES`: the listed diagonals are
raised to the sigma the consistency read measured, once per epoch on
kfMeas.postCov, before the P starts a window, an anchor or a long coast
and before the separation test subtracts it, so the floor cancels in the
separation and only its growth through the dynamics remains. At a commit
the host takes its own P back, not the floored one. The lists are empty in
the repository; a host whose P reads consistent leaves them empty.

The order of honesty is therefore: the sensor from its datasheet and its
Allan plots, the host's Q from the sensor plus named terms, the gate's
coasts from the sensor, and the host's P floored where its own
consistency read says so. `test_coastQ` checks the Q against the formula
and the floor against its definition; the test references, run with the
flag off and empty floors, stand unchanged.

## 14. The monitors' S carries the gate's covariance

Two latches survived on the clean 2000 s run with B at [8 8 10] and the
floors in, at 446 s on D and at 1262 s on E, and the host's R was scaled
on neither (`tools/measNoiseCheck`; `tools/cpiReplay` with the mode
`host` reads the same sums): what the gate reads is the host's S with
the nominal R, normalised by the host's P, and the host's P is honest
about the clock only right after a gap. At 446 s two satellites dropped
and the clock shift of both systems stepped by 16 m together for two
epochs, back by 447 s, a common mode the CPI projection, holding the
clocks at the host's P after seven minutes of aiding, read as 26 m on D;
with the clocks free the rows put the shift at 1.4 m. A 58 m step at
322 s, right after a gap, latched nothing, because the host's clock P had
grown there. At 1262 s the rows agreed on a 16 m east shift of their own
solution while the nav sat 9 m east of the truth: a real disagreement,
held in the innovations for nine epochs by a host that takes about a
tenth of it per epoch.

Decision: the S the monitors read is the host's S plus H dP H'
(`SPF_floorInnovCov`), dP the increase of the floored P over the host's
on every state (`COAST_SIGMA_FLOOR_*`) and, on the receiver clock bias
states (`CPI_CLOCK_STATES`), the variance of `CPI_CLOCK_SIGMA` over the
host's. It is applied once per epoch at the top of the gate, so the CPI
in NOMINAL and the monitors on the trial in PROBATION read it, and to
the S the CPI against the coast reads in COAST; the re-validation keeps
the bounded clock of section 9, since a vertical spoof is mostly a
common mode and must not be re-admitted. A shift common to every range
is a clock by definition: with the clock free the projection reads the
position shift the rows agree on, the same quantity at every epoch of
every run, where before it read a common mode as position by a share
that followed the aiding history. That share is sigma^2 over sigma^2 +
R / n, the rows' own clock resolution: 56 % at the host's 2 m, 89 % at
5 m, 97 % at 10 m, 99.7 % at 30 m; above ten times
sqrt(R_range / REVAL_MIN_MEAS) the value is immaterial, 25 m with 5 m
rows and four of them, 30 m for the host, and no run is involved in it.
The host lists one clock bias state per constellation it ranges on
([16 18 21] in its layout): each system's rows carry their own column,
so a floor on one system alone would still read the other's common mode
as position; a state whose rows are absent adds nothing, and the clock
drift states belong to the rate rows and stay out. The same holds for
`COAST_BOUND_STATES` of section 9, which bounds one clock per
constellation on such a host.

Measured on the clean run (`tools/cpiReplay` with the floors as its
ninth argument): the two windows read 0.81 and 0.00 against 1.42 and
1.84, no other window rises, the sigma of the reading goes from 2.84,
2.31, 2.26 m to 3.92, 3.37, 5.47 m on N, E, D, and the offset one window
latches on from 14.1, 12.9, 14.8 m to 16.4, 15.2, 21.7 m: 2 m on the
horizontal axes, 7 m on the vertical, where the clock and the altitude
share a signature and the SS and the baro row keep the vertical. The
repository leaves the clock states empty and the sigma at 0, and the
floor lists empty, so the references stand; `test_floorInnovCov` checks
the matrix against its formula and the two CPI readings, a clock jump
and a position step.

## 15. The monitors return to shadow above the protection level limit

On the long vehicle run the receiver is gone for 450 s (2400 to 2850 s).
The INS drifts 383 m north and 328 m east in the gap, the heading being
2.5 degrees off, while the host's position sigma reaches about 25 m and
the protection level passes 100 m for 342 s (`coastBudgetExceeded` on
684 epochs). At the reacquisition the first window reads the 350 m jump
at 14 sigma on N and latches (`tools/matrixRow`: CPI 2.77 at 2870.48 s).
The run was in shadow, so the host went on updating and the gate's coast
was the host's own prior, back on the GNSS: the re-validation then read
the host's innovations, 0.49 of its threshold on 31 rows once the spike
filter released them, and handed back at 2905.52 s. An acting gate would
have held the nav on the drifted INS, 500 m from the truth, and the
re-validation against that coast would have failed until the hold
limit, where the paper's recovery re-admits the GNSS, flagged.

The arming rule already says what the gate can certify: the monitors arm
only when the protection level is below `ARM_PL_MAX`, because above it
the INS-only solution cannot judge the GNSS. Arming was sticky, so a gap
in which the PL passes the limit left the monitors armed on a solution
they could not judge. Decision: on an armed NOMINAL epoch whose
protection level exceeds `ARM_PL_MAX` the monitors return to shadow, the
count from zero, exactly as after a re-initialisation; the windows and
margins go on being logged, nothing latches, a shadow alarm restarts the
count, and the gate re-arms after `ARM_EPOCHS` qualifying epochs on the
re-converged host. `info.armed` says so through the gap and for the
60 s after it, which is the honest statement: after a gap that long the
INS cannot certify the first fix, and a spoofer who jams and then
transmits is met by the host's own filter, not by a coast that is 500 m
wrong. `test_shadowArm` case D: a gap long enough for the PL to pass the
limit disarms the monitors, an 80 m jump at the reacquisition raises no
latch, the gate re-arms `ARM_EPOCHS` epochs after the last shadow alarm
and a 30 m step afterwards latches within two epochs. No constant was
added.

## 16. The host's R rule and the row screen

What the records showed (long run, the braking before the tunnel at 2355
to 2396 s, and the quiet legs at 1000 to 1500 s). The host's measurement
setup sets each row's R from the row's misclosure against the filter:
nominal within one nominal sigma, the misclosure squared above it, a cap
at three sigma on the ranges and at fifteen on the range rates. Read from
the records with both R on them (`updateStepCheck`, `dopplerCheck`): the
state takes exactly the Kalman step for the R the gain is built with, and
the covariance is updated with the same R, so the host is a consistent
filter, but two to seven times slower than the nominal R would make it,
and slowest where the INS is most wrong. A row that disagrees with the
filter by more than one sigma is weighted down until it counts as one
sigma, so the filter never experiences a disagreement larger than one
sigma on any row and corrects an error at the rate of a one-sigma
innovation whatever its size. In the braking, 35 to 9 m/s under 21 mrad
of heading error, the honest Doppler rows sat at four to six sigma, were
weighted down 25 to 42 times, and the north velocity error grew to
0.55 m/s, the heading error times the speed lost, against fifteen rows
that pointed at it every epoch. The manoeuvres that should correct a
heading are where the rule bites hardest; with the heading and gyro bias
process noise small, P shrinks over the hour to 0.34 mrad against a 21 mrad
error. The lock-in is these two together, both in the host.

The rule confuses two things that need opposite treatment. One row
disagreeing with the other rows is an outlier, multipath, a cycle slip, a
bad satellite, and deserves to be weighted down or dropped. The rows
agreeing among themselves and disagreeing with the filter together are the
filter's own error, or under attack a spoof, and must pass at full weight:
following consistent rows is the host's job, deciding whether they are
honest is the gate's, and the gate reads the nominal R by design (sections
11 and 14, `STRUCT_SPF.setKfMeasBoth`).

Two steps, in order. First, constants only: the host's thresholds moved
from one to three nominal sigma and the caps brought down to meet them, so
the rule reads nominal weight within three sigma of the filter and a ninth
beyond, on both row types. That restores the quiet legs and bounds the loss
in a manoeuvre; it still cannot tell an outlier from a filter error.
Second, when the manoeuvres still come up short, `SPF_rowScreen`: the rows
of each type fitted alone, ranges to a position shift and one clock per
constellation, range rates to a velocity and one drift per constellation,
the worst row beyond `SCREEN_K` (3) times its residual's own sigma dropped
(the standardised residual of Baarda's w-test, so a low satellite that
pulls the fit toward itself is judged on what it hides), the fit repeated
while a row of redundancy remains; kept rows at the nominal R whatever
their misclosure, dropped rows at the nominal R times `DROP_FACTOR`. It is
one call in the host's measurement setup after the rows are assembled and
before the pressure row, with the loops left as they are; the host's
override keeps its own branch. The same screen is the row screen of
`cpiReplay`. `test_rowScreen`: a 60 m single row dropped; a 40 m row on a
10 degree satellite dropped behind an 80 m one; 20 m and 0.5 m/s on every
row, five sigma against the filter, nothing dropped and every row nominal;
a 30 m jump on one constellation read as its clock; no redundancy, no
decision; the override untouched.

## 17. The Doppler episode screen

What the records showed (long run, against the reference). The receiver's
Doppler goes wrong on every row at once under strong transient dynamics:
at the bend before the tunnel the rows sat 0.4 to 0.5 m/s from the
reference for 20 s, the error turning from east to north with the heading,
and the hour holds 36 such episodes of 1.5 to 7 s, most in the urban
stretches, one 40 s stretch with the rows reversing by 0.4 m/s every epoch.
No consensus test can see them, since the rows agree among themselves, and
the row screen of section 16 passes them by construction. The filter
believed the portal and put the disagreement on the tilt, and the tilt is
most of what the 475 s coast integrated into the coast error of runs 9 and
10. With the heading brought down by the process-noise changes, this one
receiver fault is what is left of that error.

The screen (`SPF_dopplerScreen`). The host's innovation is total-state,
predicted minus measured minus H times the prior error state, so the
velocity the rows are measured against is the mechanisation velocity less
the prior state's velocity, and that difference does not move at a
feedback: a hand-over moves the mechanisation and the state by the same
amount. The range-rate rows are fitted alone to a velocity and one clock
drift per constellation, and the velocity they say, that difference less
the fit, is held against a coast that follows the change of the same
difference every epoch with the previous update's own step added back,
which is the IMU's increment less the error growth the filter's model
expects; no hand-over enters and nothing is assumed about how the
mechanisation applies one. Every accepted epoch re-anchors the coast on the rows
through a slow filter, a tenth of the difference, so the anchor is the
recent honest rows and not one epoch, and a slow drift of the rows passes
by design. The tolerance is the floor, 0.2 m/s, grown with the time the
coast has run on the IMU alone, 0.005 m/s per second for the tilt and bias
drift, and with a 5 mrad heading error on the speed change since the
anchor. Two consecutive epochs beyond it start an episode; while it lasts
the rows are held out and the coast runs on the IMU; four consecutive
epochs back within the tolerance end it, plus one more for every ten
epochs the rows have been out, so rows that were wrong for long must agree
for longer; the coast is kept at a release, so rows that came only half
way back are flagged again at once. The IMU alone carries the coast for
forty epochs, flagged or without rows alike. At forty epochs of a hold the
episode is closed and the coast is re-anchored on the filter's own
velocity, nav less prior, which has run without the range-rate rows for
the whole hold and carries the ranges' aiding, a Doppler-free second
opinion the pure-IMU coast is not; the tolerance is kept as grown, so rows
within it are accepted again and rows still beyond it start a new episode,
and a long fault is held out for as long as the filter's velocity
disagrees with it. After forty epochs without rows the coast is lost and
the next rows anchor it afresh, so a tunnel does not leave a coast metres
per second off to flag the honest rows at its exit. The end of an episode
is a return to the coast, and nothing fires there. With the two rules in,
run 10's records read 18.7 percent of the Doppler flagged, 77 false flags
in 12 runs and 318 release tails against 31.6 percent, 531 and 501 before,
every episode caught; the anchor's share was then set from the same
records, a twentieth, a tenth and three twentieths giving 1041, 903 and
896 flagged epochs with 40, 40 and 39 episodes caught, so a tenth, the
largest share that loses none.

Run 11, the first with the screen in the host. The host caught the portal,
its three episodes flagged at 83, 100 and 100 percent in its own log, and
the coast came out at 950 m north and 500 m east against 1033 and 487 in
runs 9 and 10. So the portal Doppler was not what set the tilt at the
tunnel entry; what is left is the filter's own tilt accuracy, 0.03 to 0.05
degrees by its pitch and roll sigma, and 0.04 degrees over 475 s is 800 to
1000 m on its own. That coast moves only with a better tilt before the
entry or with motion information inside the tunnel, for a land vehicle the
non-holonomic constraints or an odometer, a host design question outside
the screen. The screen itself held 60 percent of the Doppler in the host
against 17 for the emulation on the same records, and the heading wandered
to 0.7 degrees mid-run with the filter starved of rows. Where host and
emulation judged the same epoch they agreed to 0.01 m/s, so the inputs
were right; the difference was the epoch they refused to judge. The host
fits only when the rate rows outnumber the unknowns, three velocity axes
plus one drift per constellation present, seven rows with three
constellations, while the emulation fitted with four; in the urban
stretches the host left many epochs unjudged, and an episode in progress
then stood with its hold frozen, since the hold advanced only on judged
epochs, so holds that should end in 20 s ran for as long as the rows stayed
few. Hence: the hold is time, advancing on every call while an episode is
on, judged or not, with the hold limit and its re-anchor applying on an
unjudged call too, since the re-anchor needs no rows; an unjudged call
neither confirms nor releases; and the emulation fits under the host's
redundancy rule, so it reproduces the host instead of flattering it. Under
that rule run 11's records read 15.9 percent flagged, 21 false flags in 6
runs and 260 release tails, 42 of 42 episodes caught, the number to expect
once the hold is time.

What run 11 settled about the coast. The horizontal acceleration error at
the tunnel entry, the slope of the velocity error over the first 40 s of
the coast, is the same in runs 10 and 11 to a few percent, about 0.009 m/s²
in the coast's direction: holding the portal out did not move it, so the
residual is systematic and the portal was the wrong suspect. It is tilt,
inside the filter's own sigma, and its source is the yaw-rate slope the
rows have carried from the start, about 0.4 m per rad/s: in every turn the
rows are 0.2 to 0.3 m/s off the same way in the body frame, which the
filter believes and puts into tilt, the bias steps at hard turns being its
signature. The host's own log of its lever-arm term shows the term formed
from a real rate, so the slope is not an uncompensated arm; the other thing
with that signature is a pairing latency, whose error in a turn is the
centripetal acceleration times the latency, and the time-pairing read has
given 30 ms on every run, which at 15 m/s is 0.45 m per rad/s. The two
part with speed, an arm's slope constant, a latency's in proportion, and
the tool reads them jointly against the reference. On run 11 the across
slope read the same at 6 and at 11 m/s, a point offset and not a latency,
and a ten-minute replay with the arm input at zero gave the arm as the data
see it: the configured length to within a few centimetres, the direction
7 degrees off in the IMU's horizontal axes, the antenna 0.3 m less to one
side than the drawing says, by the host's own view and against the
reference alike. The residual run 11 carried against the reference, 0.36 m
per rad/s, is that geometry difference; the host's own view of it, 0.55 m,
is inflated by the filter's response in turns, as the synthetic had shown;
and the 30 ms the time-pairing read gave on every run was this residual
arm, since an arm under a centripetal acceleration reads as a latency of
its length over the speed, which the zero-arm replay confirmed by reading
150 ms for 2.47 m. A transient episode the screen holds out; a consistent
error at every turn no screen can, and that fix is the arm's components in
the IMU's axes, from the data for the next run and from a survey of the
mounting for the product. The two coast forms,
compared on run 11 with the hand-over logged, read 0.0053 against
0.0052 m/s per epoch and 21 against 5 false flags in the hour, so the host
function stays on the update's step. The
flag is one input to the measurement setup: the range-rate rows' R goes to
the host's maximum while it is on, the layout stays, the ranges and the
pressure row update as before, and the gate reads every row at the nominal
R, as for the host's own reject (section 11). A wrong Doppler row in a
turn does not help the heading, it damages it; in an honest turn the rows
and the IMU agree and nothing is flagged, so the heading keeps the turns.

Emulated offline on run 10 (`dopplerEpisodeCheck`, the same logic on the
records and the run; this first form took the nav velocity alone and the
host's velocity correction from a log): all 36 episodes caught, the portal
half a second after the rows leave the reference and held to the tunnel;
false flags 1.2 percent of the hour's Doppler, release tails 4 percent,
and 5 percent in the band where the rows are 0.15 to 0.25 m/s off, which
the screen holds out and the scoring counts as neither clean nor wrong.
That emulation needed four reads on its inputs before it could be
trusted, each a tool option now: a lever-arm term in the rows' velocity,
the reference point's own offset from the IMU, how much of the host's
velocity correction the logged velocity carries and when, and the time
pairing of the rows, which is sound. The total-state form above reads the
prior and the update's step from the records and needs no correction log;
the third read becomes a check, the share of the state left in the
increment, which is zero when the mechanisation takes the whole hand-over
and the prior is the state the rows were predicted with. In the host the
screen reads the mechanisation velocity, the prior and the update's step
directly and needs none of the reads.

Run 10 under the total-state form, no correction log: every one of 40
episodes caught, the portal half a second after the rows leave the
reference; the hand-over cancels on the quiet epochs, the share of the
state left in the increment reading 0.001, so the mechanisation takes the
whole hand-over and the delivery question is closed. But 31.6 percent of
the Doppler flagged, false flags 531 epochs in 35 runs of 10 to 20 s and
release tails 501, against 1.2 and 4 percent under the first form. Runs of
that length are sparse steps in nav minus prior that the trimmed share
read does not see, and the suspect is a feedback branch in which the state
loses more than the mechanisation receives, so the difference steps by
what was lost and not taken. The tool now classes every epoch by the
state's loss between an update and the next prior, as a share of the
posterior, and reads per class how much of the loss the nav took; a class
that lost more than the nav took is named LOSS NOT TAKEN, and the per-epoch
table carries the loss and the increment's error. The read on run 10:
no class lost more than the nav took, the limited branch included, so the
host function's increment from the update's step is exact at every epoch,
and the increment's error against the reference is 0.001 to 0.006 m/s per
epoch through the false-flag runs. The cost is plain IMU coast drift: the
rows sat in the grey band from 468 s, the coast lost its anchor and ran on
the IMU alone for 27 s, arrived 0.27 m/s off, and from 498 s, the rows
clean again, drifted on at about the tolerance's own growth rate, so the
flags lasted to 545 s. The tool now carries the coast's own error against
the reference as a table column and reads its drift rate over every
flagged stretch of 10 epochs or more against DRIFT_RATE. On run 10 that
drift reads a median near zero over the 50 long holds, which mostly sit on
true episodes where the coast is right, a 90th percentile of 0.006 and a
worst of 0.013 m/s per second, and the table through 465 to 500 s shows
the mechanism: the anchor filter at a twentieth followed a drifting IMU with a
standing lag of the drift per epoch over ALPHA, about 0.1 m/s there, a
21-epoch stretch without rows then adds 0.12, and the clean rows that
return sit 0.36 from the coast while the tolerance has grown to 0.21. Two
constants are undersized for this IMU in these stretches, not the design:
ALPHA, which at 0.15 brings the standing lag to about 0.03, and DRIFT_RATE,
which at 0.010 covers the 90th percentile with margin; the floor can stay.
The variants on the records said otherwise: ALPHA 0.15 took the false
flags from 531 to 422 and lost one small episode, DRIFT_RATE 0.010 changed
nothing, since a tolerance growing at the coast's own drift rate never
catches a coast that starts 0.24 off, and a floor of 0.25 dropped the
portal's coverage to 40 percent, which settles the floor at 0.2. The
attribution of the false runs then gave the real shape: the coast error at
the runs' start had a median of 2.98 m/s, the twelve early runs all under
0.4, so the other 23 were the hour after the tunnel, where a coast carried
blind for 475 s was metres per second off and flagged the honest rows at
the exit for the rest of the run; the first form had hidden this, its rows'
velocity carrying the state and so looking like the nav itself after the
tunnel. Hence the two rules above, the blind limit and the re-anchor on
the filter's own velocity at the hold limit, both from quantities the
screen already has; on run 10's records, where the filter did use the
rows, the emulation of the second shows its best case, and run 11 shows
the real one. The same read, on the epochs without rows,
gives the real error growth as 1.47 times what the filter's model expects,
so the total-state coast, which removes the modelled part, is the better
of the two forms by a factor of three. The other coast form is
in the tool for the comparison on
the first run that logs the hand-over: the rows' velocity total-state in
both, the increment either from the prior and the step or from the
mechanisation's change plus the correction it applied, the two differing
by the error growth the filter's model expects between updates, with one
line giving each coast's error against the reference; the smaller decides
which the host function takes, and if it is the hand-over the host passes
the applied velocity correction in place of the update's step. Open: the
loss read on run 10; that comparison; the floor, left at 0.2 for the first
run with the screen; the reversing stretch at 304 to 344 s, a receiver
behaviour to look at in its raw output.
