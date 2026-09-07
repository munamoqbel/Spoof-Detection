%******************************************************************************************
% DESCRIPTION:
% One epoch of the overlapping-window dual monitor.
%
% For every OPEN window:
%   1. buffer this epoch's (gamma, S, H)
%   2. SS test via ss_monitor_step (propagates that window's coast
%      once, tests every monitored axis, Eq. 49-52)
%   3. at full length: CPI via cpi_monitor_window once per monitored
%      axis (Eq. 33/35), then free the slot. A window that lived its
%      whole life with NO alarm is reported as the anchor candidate.
%
% INPUTS:
%   - pool (PoolState)
%   - this epoch's innovation / innovation_cov / obs_matrix
%   - the watched filter's kf_state / kf_covariance (post-update)
%   - spoofInfo
%
% OUTPUTS:
%   - pool (updated)
%   - report (PoolReport - see that file)
%
% ASSUMPTIONS AND LIMITATIONS:
%
% REQUIREMENT TRACEABILITY:
%
%******************************************************************************************
%#codegen
function [poolOut, report] = monitorPool(poolIn, innovation, ...
    innovation_cov, obs_matrix, kf_state, kf_covariance, spoofInfo)

% Define variables
windowLength = CST_spfParam.WINDOW_LENGTH;
monitoredAxes = CST_spfParam.MONITORED_AXES;
numAxes = cast(numel(monitoredAxes), 'uint8');
poolOut = poolIn;

report = STRUCT_SPF.zeroMonitorReport;

for wIdx = 1:windowLength

    if (poolIn.active(wIdx))
        % ------------------------------------------------------------------
        %  1. age the window and buffer this epoch
        % ------------------------------------------------------------------
        age = poolIn.windowAge(wIdx) + 1;
        poolOut.windowAge(wIdx) = age;

        poolOut.innovationBuffer(:, age, wIdx)       = innovation;
        poolOut.innovationCovBuffer(:, :, age, wIdx) = innovation_cov;
        poolOut.obsMatrixBuffer(:, :, age, wIdx)     = obs_matrix;

        % ------------------------------------------------------------------
        %  2. Solution-Separation test (all monitored axes)
        % ------------------------------------------------------------------
        [newCoastState, newCoastCov, ssResult] = ssMonitor( ...
            poolIn.coastState(:, wIdx), poolIn.coastCovariance(:, :, wIdx), ...
            kf_state, kf_covariance, spoofInfo);

        poolOut.coastState(:, wIdx)         = newCoastState;
        poolOut.coastCovariance(:, :, wIdx) = newCoastCov;

        if (ssResult.anyAlarm)
            report.ssAlarm         = true;
            report.anyAlarm        = true;
            poolOut.hadAlarm(wIdx) = true;
            report.alarmPerAxis = report.alarmPerAxis | ssResult.alarmPerAxis;
        end % ELSE is trivial

        if (ssResult.maxProtectionLevel > report.maxProtectionLevel)
            report.maxProtectionLevel = ssResult.maxProtectionLevel;
        end % ELSE is trivial

        % ------------------------------------------------------------------
        %  3. window complete: CPI verdict per axis, then close
        % ------------------------------------------------------------------
        if (age >= windowLength)
            for idx = 1:numAxes
                axisIdx = monitoredAxes(idx);  % actual state index
                [cpiAlarm, ~, ~] = cpiMonitor( ...
                    poolOut.innovationBuffer(:, 1:windowLength, wIdx), ...
                    poolOut.innovationCovBuffer(:, :, 1:windowLength, wIdx), ...
                    poolOut.obsMatrixBuffer(:, :, 1:windowLength, wIdx), ...
                    axisIdx);

                if (cpiAlarm)
                    report.cpiAlarm                = true;
                    report.anyAlarm                = true;
                    poolOut.hadAlarm(wIdx)         = true;
                    report.alarmPerAxis(axisIdx)   = true;
                end % ELSE is trivial
            end

            if (~poolOut.hadAlarm(wIdx))
                % Whole life alarm-free: certified anchor candidate.
                % Caller commits it ONLY on a globally alarm-free epoch,
                % so the anchor always ends strictly before any detection.
                report.cleanCloseFound = true;
                report.cleanCloseState = poolOut.coastState(:, wIdx);
                report.cleanCloseCovar = poolOut.coastCovariance(:, :, wIdx);
            end % ELSE is trivial

            poolOut.active(wIdx) = false;        % free the slot
        end % ELSE is trivial
    end % ELSE is trivial
end

end

%------------------------------------------------------------------------------------------
