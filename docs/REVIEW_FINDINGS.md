# Runtime-path review (2026-09-07/08)

Four independent reviewers (paper fidelity, FSM logic, integration hazards,
Coder/types) audited the deployed path `spoof_monitor_2hz -> protectedNav ->
monitorPool -> {cpiMonitor, ssMonitor}, insCoast, revalidation, STRUCT_SPF,
CST_spfParam`. Every finding was then challenged by two independent
skeptics; 47 verdicts came back (43 confirmed, 4 refuted, all refuted ones
were "info"), 33 verdicts were lost to a session usage limit and are marked
*unverified* below.

## Confirmed correct (paper lens)
- cpiMonitor implements Eq. 17/20/29/33/35; per-epoch normalisation makes
  the H0 statistic exactly Gamma(N/2, 2).
- Omega = sigma2_gamma_u * sigma_t^2 is the paper's Eq. 37; Eq. 38 via
  gammainc(T_N/(2(1+Omega)), N/2) reproduces the reference P_MD values.
- Baked constants match the offline design to 1e-12; the per-epoch
  false-alert budget totals <= P_FA+ + P_FA- = 1e-5.
- ssMonitor implements Eq. 49-52 with correct time alignment (post-update
  KF at k vs coast propagated k-1 -> k); window semantics match Fig. 9.
- revalidation's r'(H P_C H' + R)^-1 r is a valid chi-square(m) test.

## Confirmed findings and status

| # | Finding | Sev | Status |
|---|---|---|---|
| 1 | Measurement count frozen at init; any satellite rise/set breaks the pool, CPI and re-validation | blocker | **fixed** (MAX_MEAS buffers + per-epoch numMeas; `tests/test_variable_numMeas.m`) |
| 2 | Re-validation gate is a single chi2(16) threshold; must be chi2(m) per satellite count | major | **fixed** (`REVAL_THRESHOLD_TABLE`) |
| 3 | Re-validation assumes linear z = H x; host h(x) is nonlinear | major | **fixed** (residual = innovation - H (coast - prior)) |
| 4 | Anchor bring-forward replays every past interval with the current Phi_acc/Q_acc | major | **fixed** (anchor kept live, propagated every epoch) |
| 5 | sigma_SS = 0 turns the SS gate into "alarm on any nonzero separation" | minor | **fixed** (undefined test -> no alarm) |
| 6 | alarmPerAxis indexed by state index in monitorPool but by bank position in ssMonitor | major | **fixed** (bank position everywhere) |
| 7 | Struct fields change class at runtime (info.mode, anchor.epoch, eventAnchorEpoch, freeSlot) — Coder rejects | major | **fixed** |
| 8 | anchorMissing branch (alarm before the first clean close) freezes on the spoofed posterior | major | **mitigated**: the initial state is now the startup anchor; flag still only reported on the latch epoch |
| 9 | MAX_ANCHOR_AGE can never trigger at a latch (anchor is 1 epoch old); the drift budget it was meant to enforce is a coast-time budget | major | **documented**; coast-time budget deferred by the user; `info.coastEpochs` telemetry added |
| 10 | No maximum coast time / dead-reckoning declaration (guide rule 4) | major | **deferred** by the user until the system is integrated and tested |
| 11 | run_recovery tests `isnan(t_anchor)` but the sentinel is 0 | minor | **fixed** |
| 12 | Coast = Phi_acc * state is the error-state time update, not an INS mechanization; meaning of kfState / reset_state depends on the host's feedback convention | blocker | **fixed**: engine rewritten in update-increment form (`x+ - x_bar`), host owns its filters, gate returns corrections applied through the host's own setKF (see `docs/HOST_2HZ_WIRING.m`) |
| 13 | kfState mirror + Eq. 49 collapse if the host zeroes its error state after feedback | blocker | **fixed**: the host reduces `KF.states` by `gain*states` after feedback (case A); separation is now the accumulated increment `d <- Phi d + (x+ - x_bar)`, invariant to that bookkeeping |
| 14 | reset_state / reseed / feedback semantics only defined for the linear total-state harness | major | **fixed**: `nav.applyCorrection/correction/covar` at latch and commit, `kfCommand.startTrial` = copy the coasting KF; harness (`recovery_nav_sim`, `scheduler_sim`) follows the same contract |
| 15 | PROBATION trusts the host to apply the reseed and never verifies it | minor | **documented** (contract in spoof_monitor_2hz header) |
| 16 | Host KF left free-running on spoofed GNSS after LATCH/VETO until the next probation reseed | minor | by design (fed architecture); documented |
| 17 | Persistent sys/epoch: a host `clear functions` silently re-arms feedback | major | **documented**; host must pass resetRequest deliberately |
| 18 | SS window runs N-1 tests / coast propagations while CPI uses N innovations | info | by design (CPI-aligned window); documented |
| 19 | Anchor is the clean window's INS-only coast rather than the window's final KF solution | info | by design (conservative; the coast never ingested the suspect GNSS) |
| 20 | altXCheck struct unused; no independent altitude source | info | left for a later baro cross-check |

## Unverified (skeptic pass lost to the usage limit)
- PROBATION/handback require the host EKF to run open-loop for 18 epochs and accept a reseed (integration).
- Eq. 50 variance can go negative/~0 if the host's own P propagation differs from spoofInfo.phiAcc/qAcc (numerical) — mitigated by #5.
- CPI assumes a joint innovation vector with the full S against the prior (sequential-scalar hosts would need a change).
- NaN/Inf from an ill-posed host S or P are swallowed silently.
- Pool struct (~1.3 MB) copied by value several times per epoch; obsMatrixBuffer stores all 60 columns although only the position columns are used.
- "Deterministic rounding" scalar loops in cpiMonitor follow two mldivide calls, so they do not buy bit-exactness.
- Coast covariance symmetrised in insCoast but not in ssMonitor (no observable effect).
