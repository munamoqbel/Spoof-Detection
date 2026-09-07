%******************************************************************************************
% DESCRIPTION:
% One epoch of the spoofing-protected navigator.
%
% MODES (constants in CST_spfMode)
%   NOMINAL    KF + monitors; clean-closed windows refresh the anchor;
%              any alarm => LATCH (sticky).
%   COAST      output = INS-only; revalidation_test counts consecutive
%              passes; full dwell => PROBATION.
%              Monitor silence NEVER exits this mode (capture).
%   PROBATION  trial filter takes GNSS under monitor watch; output stays
%              on the coast; alarm => VETO; full quiet probation =>
%              COMMIT (handback).
%
% INPUTS (all from the HOST's 2 Hz GNSS update of this epoch):
%   - sys            structure of type STRUCT_SPF.setSys (persistent)
%   - innovation     [MAX_MEAS x 1]  host innovation  y = z - h(x_prior)
%   - innovationCov  [MAX_MEAS x MAX_MEAS]  S = H P_prior H' + R
%   - obsMatrix      [MAX_MEAS x n]  H (host linearisation)
%   - numMeas        scalar  valid measurement rows this epoch (0 = none)
%   - priorState     [n x 1]  host predicted state the innovation refers to
%   - kfState        [n x 1]  host post-update state
%   - kfCovariance   [n x n]  host post-update covariance
%   - measNoiseCov   [MAX_MEAS x MAX_MEAS]  R
%   - spoofInfo      .phiAcc / .qAcc  interval Phi / Q since the last epoch
%   - epoch          scalar  epoch counter (event stamps, anchor age)
%
% OUTPUTS:
%   - sys          structure of type STRUCT_SPF.setSys
%   - spoofTel     structure of type STRUCT_SPF.setTel
%                  .kfCommand.reseedKF/State/Cov  -> host KF (probation open)
%                  .nav.state/covar/sigmaPosition -> protected navigation output
%                  .info.*                        -> telemetry / events
%
% ASSUMPTIONS AND LIMITATIONS:
%
% REQUIREMENT TRACEABILITY:
%
%******************************************************************************************
%#codegen
function [sys, spoofTel] = protectedNav(sys, innovation, ...
    innovationCov, obsMatrix, numMeas, priorState, kfState, kfCovariance, ...
    measNoiseCov, spoofInfo, epoch)

info = STRUCT_SPF.zeroInfo;
kfCommand = STRUCT_SPF.setCommand(false, sys.filter.state, ...
    sys.filter.covariance);

switch sys.mode
    % ==================================================================
    case CST_spfMode.NOMINAL
    % ==================================================================

        % ---- 1. normal Kalman filter epoch ----
        % ---- CHANGED: mirror the caller's updated KF (was: update) ----
        sys.filter.state = kfState;
        sys.filter.covariance = kfCovariance;

        % ---- 2. monitor bank on the filter ----
        [sys.pool, report] = monitorPool(sys.pool, ...
            innovation, innovationCov, obsMatrix, numMeas, ...
            sys.filter.state, sys.filter.covariance, spoofInfo);

        % ---- 2b. keep the anchor LIVE: propagate it INS-only to this
        %          epoch (E28/E31) so a latch can use it without a
        %          bring-forward loop (exact for time-varying Phi/Q) ----
        if (sys.anchor.valid)
            [sys.anchor.state, sys.anchor.covariance] = insCoast( ...
                sys.anchor.state, sys.anchor.covariance, spoofInfo);
        end % ELSE is trivial

        info.ssAlarm            = report.ssAlarm;
        info.cpiAlarm           = report.cpiAlarm;
        info.alarmPerAxis       = report.alarmPerAxis;
        info.maxProtectionLevel = report.maxProtectionLevel;

        % ---- 3. refresh anchor (only on alarm-free epochs, so the
        %         anchor always ends strictly before detection) ----
        if (report.cleanCloseFound) && (~report.anyAlarm)
            sys.anchor.valid      = true;
            sys.anchor.state      = report.cleanCloseState;   % already at 'now'
            sys.anchor.covariance = report.cleanCloseCovar;
            sys.anchor.epoch      = uint32(epoch);
        end

        % ---- 4. latch on any alarm, else open the next window ----
        if (report.anyAlarm)
            info.eventLatched = true;
            anchorAge = uint32(epoch) - sys.anchor.epoch;

            if (sys.anchor.valid) && (anchorAge <= CST_spfParam.MAX_ANCHOR_AGE)
                % fall back to the certified-clean coast (kept live at
                % 'now' by step 2b)
                sys.filter.state      = sys.anchor.state;
                sys.filter.covariance = sys.anchor.covariance;
                info.eventAnchorEpoch = sys.anchor.epoch;
            else
                % no clean window ever closed: fallback unavailable.
                info.eventAnchorEpoch = 0;
                info.anchorMissing    = true;
            end

            sys.pool       = STRUCT_SPF.closeAllWindows(sys.pool);
            sys.dwellCount = 0;
            sys.coastCount = 0;
            sys.mode       = CST_spfMode.COAST;

        else
            sys.pool = STRUCT_SPF.openWindow(sys.pool, ...
                sys.filter.state, sys.filter.covariance, ...
                innovation, innovationCov, obsMatrix, numMeas);
        end

    % ==================================================================
    case CST_spfMode.COAST
    % ==================================================================

        % ---- 1. output = INS-only propagation (GNSS severed) ----
        [sys.filter.state, sys.filter.covariance] = insCoast( ...
            sys.filter.state, sys.filter.covariance, spoofInfo);

        sys.coastCount = sys.coastCount + 1;

        % ---- 2. test the (untrusted) GNSS against the coast ----
        [passed, q_value] = revalidation(innovation, obsMatrix, ...
            measNoiseCov, numMeas, priorState, ...
            sys.filter.state, sys.filter.covariance);
        info.qReval = q_value;
        info.revalComputed = true;

        if (passed)
            sys.dwellCount = sys.dwellCount + 1;
        else
            sys.dwellCount = 0;
        end

        % ---- 3. dwell full: open PROBATION (never direct handback) ----
        if (sys.dwellCount >= CST_spfParam.REVAL_DWELL_REQUIRED)
            info.eventProbationStarted = true;
            sys.trial.state       = sys.filter.state;
            sys.trial.covariance  = sys.filter.covariance;
            sys.probationCount    = 0;
            sys.dwellCount        = 0;
            sys.pool              = STRUCT_SPF.closeAllWindows(sys.pool);
            sys.mode              = CST_spfMode.PROBATION;
            % ---- CHANGED: order the caller to reseed its KF ----
            kfCommand = STRUCT_SPF.setCommand(true, ...
                sys.filter.state, sys.filter.covariance);
        end

    % ==================================================================
    case CST_spfMode.PROBATION
    % ==================================================================

        % ---- 1. output STAYS on the coast ----
        [sys.filter.state, sys.filter.covariance] = insCoast( ...
            sys.filter.state, sys.filter.covariance, spoofInfo);

        % ---- 2. trial filter takes the GNSS update ----
        % ---- CHANGED: mirror the caller's updated KF (was: update) ----
        sys.trial.state = kfState;
        sys.trial.covariance = kfCovariance;

        sys.coastCount = sys.coastCount + 1;

        % ---- 3. diagnostic: GNSS-vs-coast, log only ----
        [~, q_value] = revalidation(innovation, obsMatrix, ...
            measNoiseCov, numMeas, priorState, ...
            sys.filter.state, sys.filter.covariance);
        info.qReval = q_value;
        info.revalComputed = true;

        % ---- 4. monitor bank on THE TRIAL filter ----
        [sys.pool, report] = monitorPool(sys.pool, ...
            innovation, innovationCov, obsMatrix, numMeas, ...
            sys.trial.state, sys.trial.covariance, spoofInfo);

        info.ssAlarm            = report.ssAlarm;
        info.cpiAlarm           = report.cpiAlarm;
        info.alarmPerAxis       = report.alarmPerAxis;
        info.maxProtectionLevel = report.maxProtectionLevel;
        % probation windows never refresh the anchor: the trial is not
        % yet trusted, so report.clean_close_* is deliberately ignored.

        if (report.anyAlarm)
            % ---- VETO: monitors caught the trial; discard it ----
            info.eventProbationVetoed = true;
            sys.pool       = STRUCT_SPF.closeAllWindows(sys.pool);
            sys.dwellCount = 0;
            sys.mode       = CST_spfMode.COAST;
        else
            sys.pool = STRUCT_SPF.openWindow(sys.pool, ...
                sys.trial.state, sys.trial.covariance, ...
                innovation, innovationCov, obsMatrix, numMeas);

            sys.probationCount = sys.probationCount + 1;

            if (sys.probationCount >= CST_spfParam.PROBATION_LENGTH)
                % ---- COMMIT (handback) ----
                info.eventHandback    = true;
                sys.filter.state      = sys.trial.state;
                sys.filter.covariance = sys.trial.covariance;

                % just certified by a full quiet probation: new anchor
                sys.anchor.valid      = true;
                sys.anchor.state      = sys.filter.state;
                sys.anchor.covariance = sys.filter.covariance;
                sys.anchor.epoch      = uint32(epoch);

                sys.coastCount = 0;
                sys.mode = CST_spfMode.NOMINAL;
                % pool carries over seamlessly
            end
        end
end

% ---- navigation output of this epoch ----
nav = STRUCT_SPF.setNav(sys.filter.state, sys.filter.covariance, zeros(3, 1));
for axisIdx = 1:3
    nav.sigmaPosition(axisIdx) = ...
        sqrt(max(sys.filter.covariance(axisIdx, axisIdx), 0.0));
end

info.mode        = sys.mode;
info.dwellCount  = sys.dwellCount;
info.coastEpochs = sys.coastCount;    % time-in-coast (COAST + PROBATION), epochs

spoofTel = STRUCT_SPF.setTel(info, kfCommand, nav);

end

%------------------------------------------------------------------------------------------
