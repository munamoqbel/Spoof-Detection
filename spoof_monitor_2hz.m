%******************************************************************************************
% DESCRIPTION:
% Deployment wrapper for YOUR 2 Hz GNSS function.
%
% Keeps sys / cfg / epoch alive between calls with PERSISTENT storage,
% per the two-function architecture:
%
%   100 Hz side ----> Phi_acc, Q_acc (accumulated over the interval)
%   100 Hz side <---- mode_flag  (feedback ON iff NOMINAL)
%               <---- reset_flag, reset_state (consume at latch/commit)
%
% Call this once per GNSS epoch from your 2 Hz function (or merge its
% body into it). In your real system, replace the internal
% kalman_update_step call chain by feeding YOUR EKF's innovation and
% post-update state into protected_nav_step's structure - see the
% integration guide.
%
% INPUTS:
%   - measurement, obs_matrix, meas_noise_cov   this epoch's z, H, V(R)
%   - Phi_acc, Q_acc        interval matrices from the 100 Hz accumulators
%   - initial_state, initial_covariance   used ONLY at (re)initialisation
%   - reset_request         true forces a full re-init (commanded restart)
%
% OUTPUTS:
%   - mode_flag     NavMode.NOMINAL / COAST / PROBATION  -> 100 Hz side
%                   (feedback enabled iff mode_flag == NavMode.NOMINAL)
%   - reset_flag    true exactly on latch and on commit epochs
%   - reset_state   [n x 1] state the 100 Hz mechanization must load when
%                   reset_flag is true (consume once, at its next cycle)
%   - nav, info     navigation output and telemetry (see NavInfo.m)
%
% ASSUMPTIONS AND LIMITATIONS:
% NOTE for scripts/tests: persistents survive between runs - call
% 'clear functions' (or pass reset_request = true) to start fresh.
%
% REQUIREMENT TRACEABILITY:
%
%******************************************************************************************
%#codegen
function [mode_flag, reset_flag, reset_state, reseed, spoofTel] = ...
    spoof_monitor_2hz(innovation, innovationCov, obsMatrix, ...
    kfState, kfCovariance, measurement, measNoiseCov, spoofInfo, ...
    initialState, initialCovariance, resetRequest, numMeas)

persistent sys epoch

if resetRequest || isempty(sys)
    mode = CST_spfMode.NOMINAL;
    filter = STRUCT_SPF.setFilter(initialState, initialCovariance);
    trial = STRUCT_SPF.setTrial(initialState, initialCovariance);
    pool = STRUCT_SPF.zeroMonitorPool(numMeas);
    anchor = STRUCT_SPF.setAnchor(false, initialState, initialCovariance, 0);

    sys = STRUCT_SPF.setSys(mode, filter, trial, pool, anchor, 0, 0);
    epoch = uint32(0);
end

epoch = epoch + uint32(1);

[sys, spoofTel] = protectedNav(sys, innovation, innovationCov, ...
    obsMatrix, kfState, kfCovariance, measurement, ...
    measNoiseCov, spoofInfo, epoch);

mode_flag = spoofTel.info.mode;
reset_flag = spoofTel.info.eventLatched || spoofTel.info.eventHandback;
reset_state = spoofTel.nav.state; % for the 100 Hz side
reseed = spoofTel.kfCommand; % for YOUR KF (2d)

end
% ------------------------------------------------------------------------------------------
