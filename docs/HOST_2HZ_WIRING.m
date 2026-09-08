%% HOST_2HZ_WIRING.m   (read-me file - do not run)
% How the host's 2 Hz function drives spoofMonitor2hz. Written against the
% host's own design: operational KF struct 'KF' (.states, .covariance,
% .stateFB, .inValid); the 100 Hz function extrapolates KF every tick,
% applies KF.stateFB to the mechanization and resets it (biases excluded),
% and accumulates propTel.accumPhi / accumQ in spfAccumProp; the 2 Hz
% setKF loads stateFB from the updated states and reduces KF.states by
% gain*states.
%
% WHY INCREMENTS. Because setKF splits the estimate between the
% mechanization (via stateFB) and the residual KF.states, KF.states alone
% is not the navigation error. The gate therefore never uses an absolute
% state: it accumulates the update increments  dx = xPost - xPrior  (= K*y)
% of whichever filter is ACTIVE, propagated with the host's own interval
% matrices. Every decision it returns is either a mode or a complete
% (x+, P+) for the operational KF to be passed through the host's normal
% setKF, unchanged.
%
% ----------------------------------------------------------------------
%  inputs the gate needs from kfUpdate (this epoch, first numMeas rows)
% ----------------------------------------------------------------------
%   kfUpdate must add four outputs it already computes internally:
%   y       = z - h(xPrior)   the innovation actually used by the update,
%             recomputed at kfUpdate level (same expression as inside,
%             including any -H*xPrior term if your states persist)   [numMeas x 1]
%   H       = H(xPrior)       Jacobian                                [numMeas x 60]
%   R       = diag(measurementZ)^2                                    [numMeas x numMeas]
%   numMeas                   valid rows this epoch (0 = no GNSS)
%   NOT needed from kfUpdate: S. The gate forms S = H*P_bar*H' + R from
%   P_bar = activeKF.covariance as it goes INTO kfUpdate (extrapolated at
%   100 Hz), so the SVD-based inverse deep inside stays untouched.
%   xPrior  = activeKF.states going INTO the update
%   xPost, PPost = x+, P+ straight OUT of kfUpdate, BEFORE setKF's
%                  feedback / gain bookkeeping
%   propTel.accumPhi, .accumQ from spfAccumProp (published before its reset)
%
%   kfMeas = STRUCT_SPF.kfMeasFromUpdate(y, H, R, numMeas, xPrior, PPrior, xPost, PPost);
%   (arrays may be exact-size or padded to CST_spfParam.MAX_MEAS = 30;
%    a host that does have S can call STRUCT_SPF.setKfMeas instead)
%
% ----------------------------------------------------------------------
%  2 Hz function
% ----------------------------------------------------------------------
% function [KF, kfTel, spoofTel] = twoHzFunction(navOut, KF, satellite, propTel, ...)
%
% persistent trialKF spoofMode firstCall
% if isempty(trialKF)
%     trialKF   = KF;
%     spoofMode = CST_spfMode.NOMINAL;
%     firstCall = true;
% end
%
% inProbation = (spoofMode == CST_spfMode.PROBATION);
% kfUpdated   = (spoofMode == CST_spfMode.NOMINAL);
%
% if inProbation
%     % bring the trial to this epoch with the SAME interval matrices the
%     % 100 Hz side used for KF, then let it eat the GNSS. NO setKF
%     % bookkeeping on the trial: its states keep the full increments.
%     trialKF.states     = propTel.accumPhi * trialKF.states;
%     trialKF.covariance = propTel.accumPhi * trialKF.covariance * propTel.accumPhi' + propTel.accumQ;
%     trialKF.covariance = (trialKF.covariance + trialKF.covariance') / 2;
%     activeKF = trialKF;
% else
%     activeKF = KF;                      % NOMINAL: normal; COAST: scratch copy
% end
%
% xPrior = activeKF.states;  PPrior = activeKF.covariance;
% [kfUpdate, y, H, R, numMeas] = kfUpdate(measurement, activeKF, ...);      % your pipeline + 4 outputs
% kfMeas = STRUCT_SPF.kfMeasFromUpdate(y, H, R, numMeas, xPrior, PPrior, ...
%                                      kfUpdate.states, kfUpdate.covariance);
%
% spoofTel  = spoofMonitor2hz(kfMeas, propTel, firstCall);
% firstCall = false;
% spoofMode = spoofTel.info.mode;
%
% if inProbation
%     trialKF = kfUpdate;                                 % trial keeps the update (no setKF)
% elseif kfUpdated
%     KF = setKF(kfUpdate, KF);                           % NOMINAL: your existing path
% end
% % COAST: KF is left as extrapolated by the 100 Hz side (= the INS-only
% % coast); no stateFB is produced, so the mechanization coasts by itself.
%
% if spoofTel.nav.applyCorrection                          % LATCH or COMMIT
%     % (nav.state, nav.covar) is a complete update result for the
%     % OPERATIONAL KF: hand it to your normal setKF, unchanged
%     kfClean            = kfUpdate;
%     kfClean.states     = spoofTel.nav.state;
%     kfClean.covariance = spoofTel.nav.covar;
%     KF = setKF(kfClean, KF);                             % sets stateFB and states as usual
% end
%
% if spoofTel.kfCommand.reseedKF                           % probation opens
%     trialKF = KF;                                        % trial starts ON the coast
%     trialKF.covariance = spoofTel.kfCommand.reseedCov;   % (equals KF.covariance)
% end
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
