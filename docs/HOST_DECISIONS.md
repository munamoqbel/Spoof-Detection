# Host integration decisions

Companion to `HOST_2HZ_WIRING.m`. Records what the host's 2 Hz function does
with the operational `KF`, the `trialKF` and the gate outputs in every case,
and what each gate mode compares against. `GAIN` is the host's constant
feedback gain (0.4). "Mode" is the gate mode **after** this epoch's call to
`spoofMonitor2hz`; "in probation" is the mode **before** it.

## 1. Definitions

| Symbol | Meaning |
|---|---|
| `xPrior` | `activeKF.states` captured before calling `kfUpdate` (extrapolated at 100 Hz) |
| `xPost` | `kfUpdate` output: `activeKF.states + K*y` |
| `PPrior`, `PPost` | `activeKF.covariance` before, `kfUpdate` covariance after |
| `y`, `H`, `R`, `numMeas` | innovation `z - h(xPrior)`, Jacobian, measurement noise, valid rows |
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
| COAST -> PROBATION | `kfCommand.reseedKF`, `info.eventProbationStarted` | drain, as COAST | drain, as COAST | input unchanged | discarded; then `trialKF = KF`, `trialKF.covariance = reseedCov` (equals `KF.covariance`) |
| PROBATION (stays) | none | zero except biases (frozen) | extrapolated, unchanged | input unchanged | belongs to the trial (see section 3) |
| VETO (PROBATION -> COAST) | `info.eventProbationVetoed` | drain, as COAST | drain, as COAST | input unchanged | trial discarded |
| COMMIT (PROBATION -> NOMINAL) | `nav.applyCorrection`, `info.eventHandback` | `GAIN * nav.state` | `nav.state - GAIN * nav.state` | `nav.covar` (trial `PPost`) | `nav.state` is the trial's `xPost` |
| No usable anchor at alarm | `info.anchorMissing`, no correction | as COAST | as COAST | input unchanged | discarded |
| Alignment | gate not called (`zeroTel`), `firstCall = true` | as today | as today | as today | as today |

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
| Probation opens (`reseedKF`) | `trialKF = KF` after this epoch's drain; `trialKF.covariance = kfCommand.reseedCov` |
| Start of each PROBATION epoch | `trialKF.states = accumPhi * trialKF.states`; `trialKF.covariance = accumPhi * P * accumPhi' + accumQ`, symmetrised (the trial has no 100 Hz extrapolation of its own) |
| Update | `kfUpdate(measurement, trialKF)`; `h(x)` must use the states of the struct passed in |
| After the gate call | `trialKF = kfUpdate result` (states and covariance). No `setKF`, no feedback, no drain: the trial has no mechanization of its own and its states must keep the full increments |
| Veto | trial dropped |
| Commit | trial `xPost`, `PPost` go through the normal `setKF` as `nav.state`, `nav.covar`; trial no longer used |

## 4. What each mode compares, and against what

| Mode | Active filter | Test | Compares | Against | Statistic and threshold | Result |
|---|---|---|---|---|---|---|
| NOMINAL | operational KF, updated normally | SS, every open window, every monitored axis, every epoch | window separation (increments since the window opened) | that window's INS-only coast | `abs(d_axis) > K_FALSE_ALERT * sqrt(P_C - P_KF)` | any alarm: LATCH to the anchor |
| NOMINAL | same | CPI, when a window reaches N = 10 | normalised innovation projections on the axis, `xi = f'S^-1 y / sqrt(f'S^-1 f)` | N(0,1) under no attack | `sum(xi^2) > CPI_THRESHOLD` (Gamma, 45.64) | any alarm: LATCH; a window alarm-free for its whole life refreshes the anchor |
| COAST | operational KF **not** updated (scratch update feeds `y, H, R` only) | re-validation, every epoch | raw GNSS innovation `y` | the INS-only coast (the operational KF's own prior) with its grown covariance | `y' (H P_C H' + R)^-1 y < chi2inv(1 - 1e-3, numMeas)`; `numMeas = 0` cannot pass | 10 consecutive passes: PROBATION; a fail resets the dwell. Monitors do not run |
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
| Probation opens | `eventProbationStarted` | none | none | `reseedKF = true`, `reseedCov = coast covariance` |
| Veto | `eventProbationVetoed` | none | none | none |
| Commit | `applyCorrection`, `eventHandback` | trial `xPost` | trial `PPost` | none |
| Alarm without usable anchor | `eventLatched`, `anchorMissing` | none | none | none; mode still COAST |
| Every epoch | `mode`, `ssAlarm`, `cpiAlarm`, `alarmPerAxis`, `maxProtectionLevel`, `qReval`, `dwellCount`, `coastEpochs` | | `coastCov` (diagnostic) | |

After a commit the anchor is the committed solution (separation zero,
covariance the trial's `PPost`), the probation windows carry over, and a new
attack goes through the same cycle with no limit on the number of cycles.

## 6. Deferred

- Baro-only aiding during COAST and PROBATION (vertical channel).
- Coast-time budget (Implementation Guide rule 4): each cycle costs a coast and
  widens the coarse re-validation gate.
- `MAX_ANCHOR_AGE` is a startup-only guard today.
