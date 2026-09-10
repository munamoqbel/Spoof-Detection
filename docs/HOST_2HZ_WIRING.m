%% HOST_2HZ_WIRING.m   (read-me file - do not run)
% Decision tables (setKF per case, trial handling, what each mode compares
% against, gate outputs per event): docs/HOST_DECISIONS.md
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
%       input to setKF: states = extrapolated KF.states, covar = KF.covariance
%       setKF normal case then gives stateFB = GAIN*states,
%       states = states - GAIN*states, covar unchanged (pos/vel/att; biases as today)
%   The estimate does not move (only where it is held), so the gate's coast
%   covariance stays valid. The leftover decays by (1-GAIN) per epoch:
%   with GAIN = 0.4, a 30 m latch correction shows 18 m after the latch
%   epoch, 3.9 m after 3 epochs, 0.1 m after 10 epochs (5 s), no step.
%   A smaller GAIN in this branch only slows the transition further.
%   PROBATION: same drain branch (one branch is simpler; the dwell
%   guarantees >= 10 COAST epochs before a probation opens, so the residual
%   is already < 1 % of its latch value and what the drain moves during
%   probation is negligible). The TRIAL never goes through setKF: it
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
%   GNSS ROWS ONLY: build spfMeas from the range / range-rate rows and drop
%   the pressure row from y, H and R; numMeas = number of GNSS rows (0 on
%   an outage). The increment xPost - xPrior still contains everything the
%   filter did, baro included, which is what the monitors expect. Keeping
%   the baro row out matters in COAST: ten baro-only epochs would otherwise
%   pass re-validation with no GNSS evidence.
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
% function [KF, ...] = twoHzFunction(KF, navMode, propTel, satData, ...)
%
% persistent trialKF spfMode kfCommand
% if isempty(trialKF)
%     trialKF   = KF;                        % set once as a KF structure
%     spfMode   = CST_spfMode.NOMINAL;
%     kfCommand = STRUCT_SPF.zeroCommand;
%     firstCall = true;
% else
%     firstCall = false;
% end
%
% inProbation = (spfMode == CST_spfMode.PROBATION);   % mode BEFORE this epoch's gate call
%
% % ---- 1. active filter: trial in PROBATION, operational KF otherwise ----
% if inProbation
%     if kfCommand.startTrial
%         % probation opened last epoch: the trial starts as the coasting KF.
%         % KF has ALREADY been extrapolated to this epoch by the 100 Hz side,
%         % so copy it as is (states AND covariance) and do not propagate.
%         trialKF = KF;
%     else
%         % bring the trial to this epoch with the SAME interval matrices the
%         % 100 Hz side used for KF (the trial has no 100 Hz extrapolation)
%         trialKF.states     = propTel.accumPhi * trialKF.states;
%         trialKF.covariance = propTel.accumPhi * trialKF.covariance * transpose(propTel.accumPhi) + propTel.accumQ;
%         trialKF.covariance = (trialKF.covariance + transpose(trialKF.covariance)) / 2;
%     end
%     activeKF = trialKF;
% else
%     activeKF = KF;                         % NOMINAL: normal; COAST: scratch copy
% end
%
% % ---- 2. your pipeline on the active filter ----
% measurement       = calMeas(satData, ..., activeKF);     % z and h(x) with the ACTIVE states
% [kfPost, spfMeas] = kfUpdate(measurement, activeKF);
% %   inside kfUpdate: spfMeas = STRUCT_SPF.kfMeasFromUpdate(y, H, R, numMeas, ...
% %       activeKF.states, activeKF.covariance, kfPost.states, kfPost.covariance)
% %   with the INPUT states/covariance as the prior; on an invalid update pass
% %   numMeas = 0, kfPost.states = activeKF.states, kfPost.covariance = activeKF.covariance.
% %   (kfPost, not 'kfUpdate': a variable named like the function shadows it.)
%
% % ---- 3. gate: EVERY 2 Hz epoch, also with numMeas = 0 ----
% spfTel    = spoofMonitor2hz(spfMeas, propTel, navMode.navigation, firstCall);
% %   navActive = false (alignment): the gate returns zeroTel (NOMINAL, no
% %   commands), drops its state and epoch counter and re-initialises on the
% %   first active call. firstCall is your persistent-init flag only.
% spfMode   = spfTel.info.mode;              % mode AFTER the gate call
% kfCommand = spfTel.kfCommand;              % startTrial acted on at the START of the next epoch
%
% % ---- 4. apply the epoch's result: exactly ONE setKF on the operational KF ----
% %  Your setKF (normal / no-satellite / failed cases) is not changed. The
% %  INPUT is selected before the call:
% kfEstimatedUpdate = kfPost;
% if spfTel.nav.applyCorrection
%     % LATCH (x+ of the operational KF minus the anchor separation) or
%     % COMMIT (the trial's x+): a complete (x+, P+) for the OPERATIONAL KF,
%     % through the normal case. Never write KF.states / KF.covariance
%     % directly and then call setKF on kfPost: that overwrites the clean
%     % state with the scratch update.
%     kfEstimatedUpdate.states = spfTel.nav.state;
%     kfEstimatedUpdate.covar  = spfTel.nav.covar;
%     kfEstimatedUpdate.failed = false;      % force the normal case (kfPost may be failed / no-sat)
% elseif spfMode ~= CST_spfMode.NOMINAL
%     % COAST or PROBATION: the operational KF is the INS-only coast. Feed
%     % the EXTRAPOLATED prior (KF is untouched at this point), not the
%     % scratch update. setKF reads only its input struct, so with this
%     % input its normal case IS the drain:
%     %   stateFB = GAIN*kfEstimatedUpdate.states
%     %   states  = kfEstimatedUpdate.states - GAIN*kfEstimatedUpdate.states
%     %   covar   = kfEstimatedUpdate.covar (= the extrapolated covariance)
%     %   imuDrift integrates as usual.
%     kfEstimatedUpdate.states = KF.states;
%     kfEstimatedUpdate.covar  = KF.covariance;
%     kfEstimatedUpdate.failed = false;      % a failed scratch update must not reset the KF
% end
% % NOMINAL (and alignment): kfPost as today, your existing cases apply.
% KF = setKF(kfEstimatedUpdate, KF);
%
% % ---- 5. trial bookkeeping ----
% if inProbation
%     trialKF = kfPost;                      % trial keeps its update (no setKF, no drain);
% end                                        % on VETO it is simply never used again
%
% Epoch-by-epoch this gives:
%   NOMINAL->NOMINAL   setKF(kfPost)
%   NOMINAL->COAST     setKF(nav.state, nav.covar)    (latch)
%   COAST->COAST       setKF(KF prior) = drain        (KF coasts)
%   COAST->PROBATION   drain; next epoch trialKF <- KF (start trial)
%   PROBATION->PROB.   drain; trialKF <- kfPost       (trial: no setKF ever)
%   PROBATION->COAST   drain                          (veto, trial dropped)
%   PROBATION->NOMINAL setKF(nav.state, nav.covar)    (commit)
%
% Same-epoch alternative: copy the trial right after the gate call
% (trialKF = KF, states and covariance) and
% then step 1 propagates it on EVERY probation epoch. Do not mix the two.
%
% ----------------------------------------------------------------------
%  alignment vs navigation mode, and GNSS outages
% ----------------------------------------------------------------------
%   propTel is only accumulated in NAVIGATION mode; in ALIGNMENT the host
%   publishes the defaults (accumPhi = I, accumQ = 0). The gate must not
%   run on those: with Q = 0 its coast covariance never grows, the SS
%   variance sigma_SS^2 = P_C - P_KF is underestimated, false alarms
%   follow, and a latch/COAST during alignment would stop the KF updates
%   the alignment needs. Hence the navActive input: pass
%   navMode.navigation; the gate handles the rest (zeroTel, reset, re-init
%   on the first navigation epoch, whose (x+, P+) become the startup anchor).
%
%   GNSS outage (only the pressure row available): STILL call the gate,
%   with numMeas = 0. It contributes xi = 0 to the CPI windows, cannot pass
%   re-validation, and keeps windows, anchor and coast covariances
%   propagating. Skipping the call would lose one interval of Phi/Q and one
%   increment. The epoch counter counts gate calls (anchor age is a time),
%   and is reset only in alignment.
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
