%% HOST_2HZ_WIRING.m   (read-me file - do not run)
% How the host's 2 Hz function drives spoof_monitor_2hz, written against
% the host's own sketch (operational KF struct 'KF' with .states,
% .covariance, .stateFB, .inValid; 100 Hz function extrapolates KF.states /
% KF.covariance every tick and applies KF.stateFB to the mechanization).
%
% Design choice made here (matches the host sketch): the OPERATIONAL KF is
% never updated while the gate is in COAST or PROBATION. The 100 Hz
% extrapolation then makes KF the INS-only coast by itself, and no feedback
% is generated, so the 100 Hz side needs no extra gating. The trial filter
% of PROBATION is a private persistent copy inside the 2 Hz function.
%
% ----------------------------------------------------------------------
%  inputs the gate needs from kfUpdate (this epoch, first numMeas rows)
% ----------------------------------------------------------------------
%   y      = z - h(x_bar)                 innovation           [numMeas x 1]
%   S      = H*P_bar*H' + R               innovation covariance [numMeas x numMeas]
%   H      = H(x_bar)                     Jacobian             [numMeas x 60]
%   R      = diag(measurementZ)^2                              [numMeas x numMeas]
%   xPrior = KF.states before the update  (the x_bar the innovation refers to)
%   xPost, PPost = x+, P+                 after the update
%   numMeas                               valid rows this epoch (0 = no GNSS)
%   spoofInfo.phiAcc, spoofInfo.qAcc      product / accumulation of the per-tick
%                                         phi and Q_step over the 50 ticks since
%                                         the last GNSS epoch (ins_100hz_template
%                                         step 2), restarted after each epoch
%
% ----------------------------------------------------------------------
%  2 Hz function (sketch -> exact names)
% ----------------------------------------------------------------------
% function [KF, kfTel, spoofTel] = twoHz_function(navOut, KF, satellite, spoofInfo, ...)
%
% persistent trialKF spoofMode kfCommand firstCall
% if isempty(trialKF)
%     trialKF   = KF;
%     spoofMode = CST_spfMode.NOMINAL;
%     kfCommand = STRUCT_SPF.zeroCommand;
%     firstCall = true;
% end
%
% inProbation = (spoofMode == CST_spfMode.PROBATION);
% if inProbation
%     if kfCommand.reseedKF                       % issued by the gate when probation opens
%         trialKF.states     = kfCommand.reseedState;
%         trialKF.covariance = kfCommand.reseedCov;
%     end
%     % bring the trial to this epoch with the SAME interval matrices the
%     % 100 Hz side used for KF (so trial and coast stay time-aligned)
%     trialKF.states     = spoofInfo.phiAcc * trialKF.states;
%     trialKF.covariance = spoofInfo.phiAcc * trialKF.covariance * spoofInfo.phiAcc' + spoofInfo.qAcc;
%     trialKF.covariance = (trialKF.covariance + trialKF.covariance') / 2;
%     activeKF = trialKF;
% else
%     activeKF = KF;                              % NOMINAL or COAST
% end
%
% xPrior = activeKF.states;
% [activeKF, y, S, H, R, numMeas] = kfUpdate(measurement, activeKF, ...);   % your pipeline
%
% [~, ~, ~, kfCommand, spoofTel] = spoof_monitor_2hz( ...
%     y, S, H, numMeas, xPrior, activeKF.states, activeKF.covariance, R, ...
%     spoofInfo, KF.states, KF.covariance, firstCall);
% firstCall = false;
% spoofMode = spoofTel.info.mode;
%
% if inProbation
%     trialKF = activeKF;                         % trial keeps eating GNSS
%     if spoofTel.info.eventHandback
%         KF = activeKF;                          % COMMIT: trial becomes operational
%     end
%     % VETO: nothing to do; next epoch is COAST and the trial is reseeded
%     % again when the next probation opens
% else
%     if spoofMode == CST_spfMode.NOMINAL
%         KF = activeKF;                          % normal closed-loop update
%     end
%     if spoofTel.info.eventLatched
%         KF.states     = spoofTel.nav.state;     % LATCH: operational KF <- clean anchor
%         KF.covariance = spoofTel.nav.covar;     % (see OPEN POINT below)
%     end
%     % COAST: KF is left as extrapolated by the 100 Hz side (= the coast);
%     % the update on activeKF was computed only to feed the gate
% end
%
% ----------------------------------------------------------------------
%  what the 100 Hz side must do (ins_100hz_template.m)
% ----------------------------------------------------------------------
%   - accumulate spoofInfo.phiAcc / .qAcc per tick, restart after each epoch
%   - nothing else if KF.stateFB is only ever written by kfUpdate on the
%     OPERATIONAL KF (COAST/PROBATION then generate no feedback by
%     construction); otherwise gate the feedback on
%     spoofTel.info.mode == CST_spfMode.NOMINAL
%
% ----------------------------------------------------------------------
%  OPEN POINT - the latch reset and the feedback convention
% ----------------------------------------------------------------------
% Loading KF.states <- spoofTel.nav.state moves the NAVIGATION OUTPUT to the
% clean coast only if navOut is derived from KF.states. If the pos/vel/att
% feedback is absorbed permanently into the mechanization and the
% corresponding KF.states are zeroed afterwards, the anchor must instead be
% applied as a one-shot correction (the accumulated separation
% KF.states - spoofTel.nav.state at the latch epoch, projected onto the
% fed-back components), and the SS separation must be tracked as an
% accumulated quantity. The gate supports both; the host's convention
% decides which form is switched on.
