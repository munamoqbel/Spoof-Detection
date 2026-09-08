%******************************************************************************************
% DESCRIPTION:
% One epoch of the overlapping-window dual monitor.
%
% For every OPEN window:
%   1. buffer this epoch's (y, S, H, numMeas)
%   2. SS test via ssMonitor (accumulates this window's separation and
%      propagates its coast covariance once, tests every monitored axis,
%      Eq. 49-52)
%   3. at full length: CPI via cpiMonitor once per monitored axis
%      (Eq. 33/35), then free the slot. A window that lived its whole
%      life with NO alarm is reported as the anchor candidate.
%
% INPUTS:
%   - poolIn   (STRUCT_SPF.setMonitorPool)
%   - kfMeas   (STRUCT_SPF.setKfMeas) this epoch's host update quantities
%   - propTel  (STRUCT_SPF.setPropTel) .accumPhi / .accumQ
%
% OUTPUTS:
%   - poolOut (updated)
%   - report  (STRUCT_SPF.setMonitorReport)
%
% ASSUMPTIONS AND LIMITATIONS:
%
% REQUIREMENT TRACEABILITY:
%
%******************************************************************************************
%#codegen
function [poolOut, report] = monitorPool(poolIn, kfMeas, propTel)

% Define variables
windowLength = CST_spfParam.WINDOW_LENGTH;
monitoredAxes = CST_spfParam.MONITORED_AXES;
numAxes = cast(numel(monitoredAxes), 'uint8');
poolOut = poolIn;
kfIncrement = kfMeas.postState - kfMeas.priorState;     % K*y of the watched filter
m = kfMeas.numMeas;

report = STRUCT_SPF.zeroMonitorReport;

for wIdx = 1:windowLength

    if (poolIn.active(wIdx))
        % ------------------------------------------------------------------
        %  1. age the window and buffer this epoch
        % ------------------------------------------------------------------
        age = poolIn.windowAge(wIdx) + 1;
        poolOut.windowAge(wIdx) = age;

        poolOut.innovationBuffer(1:m, age, wIdx)          = kfMeas.innovation(1:m);
        poolOut.innovationCovBuffer(1:m, 1:m, age, wIdx)  = kfMeas.innovationCov(1:m, 1:m);
        poolOut.obsMatrixBuffer(1:m, :, age, wIdx)        = kfMeas.obsMatrix(1:m, :);
        poolOut.numMeasBuffer(age, wIdx)                  = uint8(m);

        % ------------------------------------------------------------------
        %  2. Solution-Separation test (all monitored axes)
        % ------------------------------------------------------------------
        [newSeparation, newCoastCov, ssResult] = ssMonitor( ...
            poolIn.separation(:, wIdx), poolIn.coastCovariance(:, :, wIdx), ...
            kfIncrement, kfMeas.postCov, propTel);

        poolOut.separation(:, wIdx)         = newSeparation;
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
                    poolOut.numMeasBuffer(1:windowLength, wIdx), ...
                    axisIdx);

                if (cpiAlarm)
                    report.cpiAlarm                = true;
                    report.anyAlarm                = true;
                    poolOut.hadAlarm(wIdx)         = true;
                    report.alarmPerAxis(idx)       = true;   % monitored-axis index (as ssMonitor)
                end % ELSE is trivial
            end

            if (~poolOut.hadAlarm(wIdx))
                % Whole life alarm-free: certified anchor candidate.
                % Caller commits it ONLY on a globally alarm-free epoch,
                % so the anchor always ends strictly before any detection.
                report.cleanCloseFound      = true;
                report.cleanCloseSeparation = poolOut.separation(:, wIdx);
                report.cleanCloseCovar      = poolOut.coastCovariance(:, :, wIdx);
            end % ELSE is trivial

            poolOut.active(wIdx) = false;        % free the slot
        end % ELSE is trivial
    end % ELSE is trivial
end

end

%------------------------------------------------------------------------
