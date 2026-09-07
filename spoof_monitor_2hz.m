%******************************************************************************************
% DESCRIPTION:
% 2 Hz spoofing-detection GATE for the host GNSS/INS EKF (persistent state).
%
% Call ONCE per GNSS epoch from the host's 2 Hz function, AFTER the host
% has computed its innovation and run its own KF update. The gate never
% touches the host filter directly: it watches it (CPI + SS monitors on
% overlapping windows), navigates on a certified-clean INS-only coast
% when an alarm latches, tests the raw GNSS against that coast, runs a
% probation, and hands back. It only ever COMMANDS the host through the
% outputs below.
%
%   100 Hz side ----> spoofInfo.phiAcc / .qAcc   (accumulated over the interval)
%   100 Hz side <---- mode_flag                  (feedback ON iff NOMINAL)
%               <---- reset_flag, reset_state    (consume once, at latch/commit)
%   host KF     <---- reseed.reseedKF/State/Cov  (apply when reseedKF is true)
%
% INPUTS (host error-state EKF, this epoch):
%   - innovation         [MAX_MEAS x 1]  y = z - h(x_prior)   (rows 1:numMeas valid)
%   - innovationCov      [MAX_MEAS x MAX_MEAS]  S = H P_prior H' + R
%   - obsMatrix          [MAX_MEAS x n]  H, host linearisation
%   - numMeas            scalar  valid measurement rows this epoch (0 = no GNSS)
%   - priorState         [n x 1]  predicted state the innovation refers to
%   - kfState            [n x 1]  post-update state
%   - kfCovariance       [n x n]  post-update covariance
%   - measNoiseCov       [MAX_MEAS x MAX_MEAS]  R
%   - spoofInfo          struct .phiAcc [n x n], .qAcc [n x n]: interval
%                        Phi / Q accumulated over the 100 Hz ticks since
%                        the previous GNSS epoch
%   - initialState       [n x 1]  used ONLY at (re)initialisation
%   - initialCovariance  [n x n]  used ONLY at (re)initialisation
%   - resetRequest       logical  true forces a full re-init (first call /
%                        commanded restart)
%   Inputs may be passed at their exact size (numMeas rows) or padded to
%   MAX_MEAS rows; only rows/cols 1:numMeas are read.
%
% OUTPUTS:
%   - mode_flag    CST_spfMode.NOMINAL / COAST / PROBATION
%                  (100 Hz feedback enabled iff mode_flag == NOMINAL)
%   - reset_flag   true exactly on latch and on commit epochs
%   - reset_state  [n x 1] protected navigation state the 100 Hz side must
%                  load when reset_flag is true (consume once, next tick)
%   - reseed       STRUCT_SPF.setCommand: when .reseedKF is true the host
%                  KF must load .reseedState / .reseedCov (probation open)
%   - spoofTel     STRUCT_SPF.setTel telemetry: .nav (protected output),
%                  .info (mode, alarms per axis, protection level, qReval,
%                  dwell, coastEpochs, events)
%
% ASSUMPTIONS AND LIMITATIONS:
% - Position error states are CST_spfParam.MONITORED_AXES (NED, 1:3).
% - The host error state is NOT zeroed at feedback (it persists and is
%   propagated at 100 Hz), so separations are computed in state space.
% - Persistents survive between runs: call 'clear functions' (or pass
%   resetRequest = true) to start fresh. Persistent-state variant; the
%   batch harness (recovery_nav_sim) keeps 'sys' explicitly instead.
%
% REQUIREMENT TRACEABILITY:
%
%******************************************************************************************
%#codegen
function [mode_flag, reset_flag, reset_state, reseed, spoofTel] = ...
    spoof_monitor_2hz(innovation, innovationCov, obsMatrix, numMeas, ...
    priorState, kfState, kfCovariance, measNoiseCov, spoofInfo, ...
    initialState, initialCovariance, resetRequest)

persistent sys epoch

if resetRequest || isempty(sys)
    mode = CST_spfMode.NOMINAL;
    filter = STRUCT_SPF.setFilter(initialState, initialCovariance);
    trial = STRUCT_SPF.setTrial(initialState, initialCovariance);
    pool = STRUCT_SPF.zeroMonitorPool;
    anchor = STRUCT_SPF.setAnchor(false, initialState, initialCovariance, uint32(0));

    sys = STRUCT_SPF.setSys(mode, filter, trial, pool, anchor, 0, 0, 0);
    epoch = uint32(0);
end

epoch = epoch + uint32(1);

[sys, spoofTel] = protectedNav(sys, innovation, innovationCov, ...
    obsMatrix, numMeas, priorState, kfState, kfCovariance, ...
    measNoiseCov, spoofInfo, epoch);

mode_flag = spoofTel.info.mode;
reset_flag = spoofTel.info.eventLatched || spoofTel.info.eventHandback;
reset_state = spoofTel.nav.state; % for the 100 Hz side
reseed = spoofTel.kfCommand; % for the host KF

end
% ------------------------------------------------------------------------
