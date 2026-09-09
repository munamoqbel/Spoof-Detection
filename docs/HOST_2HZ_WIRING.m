%% HOST_2HZ_WIRING.m   (read-me file - do not run)
% How the host's 2 Hz function drives spoofMonitor2hz. Written against the
% host's own design: operational KF struct 'KF' (.states, .covariance,
% .stateFB); the 100 Hz function extrapolates states and covariance every
% tick (state = phi*state), applies KF.stateFB to the mechanization and
% resets it (biases excluded), and accumulates propTel.accumPhi / accumQ
% in spfAccumProp. The navigation output sees ONLY stateFB; KF.states is
% used only for estimation and extrapolation.
%
% HOST setKF AS IT EXISTS TODAY (feedback gain GAIN, one constant)
%   normal        stateFB = GAIN*xPlus;  states = xPlus - GAIN*xPlus;  cov = PPlus
%   no satellite  stateFB zeroed except biases; states = xPlus except the
%                 pressure-altitude state (1-GAIN)*xPlus; cov = PPlus
%   invalid       stateFB = 0, states = 0
%
% WHY INCREMENTS. The gate never uses an absolute state: it accumulates
% the update increments  dx = xPost - xPrior  (= K*y) of whichever filter
% is ACTIVE, propagated with the host's own interval matrices. The filter's
% estimate moves by exactly dx per update whatever GAIN does with it, so
% the monitors are invariant to the feedback split. Every decision the gate
% returns is either a mode or a complete (x+, P+) with the SAME definition
% as the kfUpdate outputs, to be passed through the normal setKF.
%
% HOST BEHAVIOUR ADDED WITH THE GATE
%   COAST         no setKF update: states and covariance stay as the 100 Hz
%                 side extrapolated them; position/velocity feedback zero;
%                 bias feedback kept. Keep the coast PURE (no baro update
%                 either): the gate's coast covariance is Phi*P*Phi' + Q of
%                 exactly this filter.
%   LATCH/COMMIT  nav.applyCorrection: (nav.state, nav.covar) through the
%                 NORMAL setKF case, nothing else that epoch.
%   invalid epoch present it to the gate as numMeas = 0 with xPost = xPrior,
%                 PPost = PPrior (no GNSS, no increment).
%
%   DRAIN IN COAST (decided: smooth transition). The normal setKF leaves
%   (1-GAIN)*nav.state in KF.states at the latch and nothing would feed it
%   while coasting, so the OUTPUT would keep the GAIN fraction of the
%   spoof offset. Therefore on every COAST epoch apply the normal split to
%   the EXTRAPOLATED state of the operational KF instead of zeroing feedback:
%       stateFB   = GAIN * KF.states;            (pos/vel/att; biases as today)
%       KF.states = KF.states - GAIN * KF.states;
%       KF.covariance unchanged
%   The estimate does not move (only where it is held), so the gate's coast
%   covariance stays valid. The leftover decays by (1-GAIN) per epoch:
%   with GAIN = 0.4, a 30 m latch correction shows 18 m after the latch
%   epoch, 3.9 m after 3 epochs, 0.1 m after 10 epochs (5 s), no step.
%   A smaller GAIN in this branch only slows the transition further.
%   PROBATION: no drain needed and none on the trial. The dwell guarantees
%   >= 10 COAST epochs before a probation opens, so the residual is already
%   < 1 % of its latch value; freeze the operational KF (pass-through, zero
%   pos/vel feedback, biases kept). The TRIAL never goes through setKF: it
%   has no mechanization of its own, its states must keep the full
%   increments, and at COMMIT its x+ goes through the normal setKF once.
%   The size of the latch correction is bounded by the protection level
%   (the monitors fire before the spoofer gets further), not by the attack.
%   Verified by tests/test_feedback_split.m (tests/hostSplitSim.m emulates
%   this setKF; estimate identical to the bare host to 1e-15 for GAIN 0.1,
%   0.4 and 1.0; residual < 2 % after 10 drained epochs).
% ----------------------------------------------------------------------
%  inputs the gate needs from kfUpdate (this epoch, first numMeas rows)
% ----------------------------------------------------------------------
%   kfUpdate must add four outputs it already computes internally:
%   y       = z - h(x)        the innovation actually used by the update,
%             recomputed at kfUpdate level with the PRIOR states. With the
%             host's definitions z = rho(p_nav) - rho_meas and
%             h(x) = rho(p_nav) - rho(p_nav - dp), y is the residual of the
%             KF-corrected position: the prior enters through h(x), so NO
%             extra H*xPrior term is added.                          [numMeas x 1]
%   H       = H(xPrior)       Jacobian                                [numMeas x 60]
%   R       = diag(measurementZ)^2                                    [numMeas x numMeas]
%   numMeas                   valid rows this epoch (0 = no GNSS)
%   NOT needed from kfUpdate: S. The gate forms S = H*P_bar*H' + R from
%   P_bar = activeKF.covariance as it goes INTO kfUpdate (extrapolated at
%   100 Hz), so the SVD-based inverse deep inside stays untouched.
%   xPrior  = activeKF.states going INTO the update
%   xPost, PPost = x+, P+ straight OUT of kfUpdate, BEFORE setKF's
%                  feedback / gain bookkeeping. The host's kfUpdate returns
%                  x+ = activeKF.states + stateUpdate with stateUpdate = K*y,
%                  so xPost - xPrior = K*y is the increment the gate uses;
%                  keep returning x+ itself (the gate hands back nav.state in
%                  the same absolute states space), not the update term alone.
%   propTel.accumPhi, .accumQ from spfAccumProp (published before its reset)
%
%   kfMeas = STRUCT_SPF.kfMeasFromUpdate(y, H, R, numMeas, xPrior, PPrior, xPost, PPost);
%   (arrays may be exact-size or padded to CST_spfParam.MAX_MEAS = 30;
%    a host that does have S can call STRUCT_SPF.setKfMeas instead)
%
% ----------------------------------------------------------------------
%  2 Hz function (complete skeleton; ONE setKF per epoch, never two)
% ----------------------------------------------------------------------
% function [KF, kfTel, spoofTel] = twoHzFunction(navOut, KF, satellite, propTel, inNavigationMode, ...)
%
% persistent trialKF spoofMode firstCall
% if isempty(trialKF)
%     trialKF   = KF;                        % set once as a KF structure
%     spoofMode = CST_spfMode.NOMINAL;
%     firstCall = true;
% end
%
% inProbation = (spoofMode == CST_spfMode.PROBATION);   % mode BEFORE this epoch's gate call
%
% % ---- 1. active filter: trial in PROBATION, operational KF otherwise ----
% if inProbation
%     % bring the trial to this epoch with the SAME interval matrices the
%     % 100 Hz side used for KF, then let it eat the GNSS. NO setKF
%     % bookkeeping on the trial: its states keep the full increments.
%     trialKF.states     = propTel.accumPhi * trialKF.states;
%     trialKF.covariance = propTel.accumPhi * trialKF.covariance * propTel.accumPhi' + propTel.accumQ;
%     trialKF.covariance = (trialKF.covariance + trialKF.covariance') / 2;
%     activeKF = trialKF;
% else
%     activeKF = KF;                         % NOMINAL: normal; COAST: scratch copy
% end
%
% % ---- 2. your pipeline on the active filter (+ 4 new outputs) ----
% xPrior = activeKF.states;  PPrior = activeKF.covariance;
% [kfUpdate, y, H, R, numMeas] = kfUpdate(measurement, activeKF, ...);
% kfMeas = STRUCT_SPF.kfMeasFromUpdate(y, H, R, numMeas, xPrior, PPrior, ...
%                                      kfUpdate.states, kfUpdate.covariance);
%
% % ---- 3. gate (navigation mode only, see next section) ----
% if inNavigationMode
%     spoofTel  = spoofMonitor2hz(kfMeas, propTel, firstCall);
%     firstCall = false;
% else
%     spoofTel  = STRUCT_SPF.zeroTel;        % NOMINAL, no commands
%     firstCall = true;                      % re-arm on the first navigation epoch
% end
% spoofMode = spoofTel.info.mode;            % mode AFTER the gate call
%
% % ---- 4. apply the epoch's result: exactly ONE setKF on the operational KF ----
% if spoofTel.nav.applyCorrection
%     % LATCH (x+ of the operational KF minus the anchor separation) or
%     % COMMIT (the trial's x+): a complete (x+, P+) for the OPERATIONAL KF.
%     % Hand it to your normal setKF so stateFB and states are set as usual.
%     % Do NOT also run the normal setKF this epoch, and do NOT write
%     % KF.states / KF.covariance directly (the mechanization would never
%     % receive the correction).
%     kfClean            = kfUpdate;
%     kfClean.states     = spoofTel.nav.state;
%     kfClean.covariance = spoofTel.nav.covar;
%     KF = setKF(kfClean, KF);
% elseif ~inProbation && (spoofMode == CST_spfMode.NOMINAL)
%     KF = setKF(kfUpdate, KF);              % NOMINAL -> NOMINAL: your existing path
% end
% elseif spoofMode == CST_spfMode.COAST
%     % the operational KF is the INS-only coast; drain the leftover of the
%     % latch correction smoothly (covariance untouched, as extrapolated):
%     KF.stateFB = GAIN * KF.states;          % pos/vel/att; biases as today
%     KF.states  = KF.states - GAIN * KF.states;
% end
% % PROBATION: operational KF frozen as extrapolated, pos/vel feedback zero.
%
% % ---- 5. trial bookkeeping ----
% if inProbation
%     trialKF = kfUpdate;                    % trial keeps its update (no setKF); on
% end                                        % VETO it is simply never used again
% if spoofTel.kfCommand.reseedKF             % probation opens THIS epoch
%     trialKF            = KF;               % trial starts ON the coast
%     trialKF.covariance = spoofTel.kfCommand.reseedCov;   % (equals KF.covariance)
% end
%
% Epoch-by-epoch this gives:
%   NOMINAL->NOMINAL   setKF(kfUpdate)
%   NOMINAL->COAST     setKF(kfClean) only            (latch)
%   COAST->COAST       drain                          (KF coasts)
%   COAST->PROBATION   frozen; trialKF <- KF          (reseed)
%   PROBATION->PROB.   frozen; trialKF <- kfUpdate    (trial: no setKF ever)
%   PROBATION->COAST   drain                          (veto, trial dropped)
%   PROBATION->NOMINAL setKF(kfClean) only            (commit)
%
% If you prefer to reseed at the START of the next probation epoch (as in
% an earlier sketch), copy KF there and skip step 1's propagation for that
% one epoch: KF.states / KF.covariance were already extrapolated by the
% 100 Hz side, so a second propagation would double it.
%
% ----------------------------------------------------------------------
%  alignment vs navigation mode
% ----------------------------------------------------------------------
%   propTel is only accumulated in NAVIGATION mode; in ALIGNMENT the host
%   publishes the defaults (accumPhi = I, accumQ = 0). The gate must NOT
%   run on those: with Q = 0 its coast covariance never grows, the SS
%   variance sigma_SS^2 = P_C - P_KF is underestimated, false alarms
%   follow, and a latch/COAST during alignment would stop the KF updates
%   the alignment needs. Step 3 of the skeleton above handles it: in
%   alignment the gate is skipped, spoofTel is zeroTel (NOMINAL, no
%   commands) and firstCall is re-armed.
%
%   On the first navigation epoch resetRequest = true makes that epoch's
%   (x+, P+) the startup anchor and the pool starts empty; the first
%   window closes clean 10 epochs (5 s) later and takes over as anchor.
%
% ----------------------------------------------------------------------
%  what the 100 Hz side must do
% ----------------------------------------------------------------------
%   - nothing new: spfAccumProp already provides propTel; KF.stateFB is
%     only written by setKF, so COAST/PROBATION produce no feedback and
%     the latch/commit corrections arrive through the same setKF path.
%   - optional annunciation: spoofTel.info.mode, .alarmPerAxis,
%     .maxProtectionLevel, .coastEpochs.
%
% ----------------------------------------------------------------------
%  sign convention check (do once)
% ----------------------------------------------------------------------
%   nav.state is built from this epoch's x+ of the active filter in your
%   own state space (x+ minus the anchor separation at LATCH; the trial's
%   x+ at COMMIT), so setKF sees exactly what a kfUpdate result looks like
%   and no sign is chosen by the gate. Verify with the QUICKSTART step-0 test
%   (add +10 m to z along H(:,idxD); the solution must move +10 m), then
%   with a scripted latch in shadow mode.
