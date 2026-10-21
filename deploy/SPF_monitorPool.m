%***************************************************************************
% DESCRIPTION:
% One epoch of the overlapping-window monitor bank. Every open window
% buffers this epoch's rows, runs the solution-separation test on all
% monitored axes, and at full length runs the CPI test per axis and
% closes. A window that lived its whole life without an alarm is
% reported as the anchor candidate, with its separation and coast
% covariance at the close.
%
% INPUTS:
%   - poolIn   STRUCT_SPF.setMonitorPool: the N window slots
%   - kfMeas   STRUCT_SPF.setKfMeas: this epoch's update quantities of
%              the watched filter
%   - propTel  STRUCT_SPF.setPropTel: interval Phi and Q
%
% OUTPUTS:
%   - poolOut  the slots after this epoch
%   - report   STRUCT_SPF.setMonitorReport: alarms (per axis, SS, CPI),
%              the largest protection level, the margins |d| / gate and
%              q / T_N, the clean-close candidate, solveFault
%
% ASSUMPTIONS AND LIMITATIONS:
% Windows are opened by the caller (STRUCT_SPF.openWindow), one per
% epoch, N slots in all: with N epochs per window every slot is in use
% once the bank is full. The nominal bias B is CST_spfParam.NOMINAL_BIAS.
%
% REQUIREMENT TRACEABILITY:
% Paper (Kujur 2024): the overlapping windows of the monitor; E24, a
% window opens on the post-update solution; the SS test Eq. 49 to 52 in
% SPF_ssMonitor; the CPI test Eq. 33 and 35 in SPF_cpiMonitor at the
% close.
%***************************************************************************
%#codegen
function [poolOut, report] = SPF_monitorPool(poolIn, kfMeas, propTel)

windowLength  = CST_spfParam.WINDOW_LENGTH;
monitoredAxes = CST_spfParam.MONITORED_AXES;
numAxes       = cast(numel(monitoredAxes), 'uint8');
poolOut       = poolIn;
kfIncrement   = kfMeas.postState - kfMeas.priorState;       % K y of the watched filter

report = STRUCT_SPF.zeroMonitorReport;

for wIdx = 1:windowLength

    if (poolIn.active(wIdx))
        % 1. age the window and buffer this epoch
        age = poolIn.windowAge(wIdx) + 1;
        poolOut.windowAge(wIdx) = age;

        poolOut.innovationBuffer(:, age, wIdx)       = kfMeas.innovation;
        poolOut.innovationCovBuffer(:, :, age, wIdx) = kfMeas.innovationCov;
        poolOut.obsMatrixBuffer(:, :, age, wIdx)     = kfMeas.obsMatrix;
        poolOut.numMeasBuffer(age, wIdx)             = kfMeas.numMeas;

        % 2. solution-separation test (Eq. 49 to 52)
        [newSeparation, newCoastCov, ssResult] = SPF_ssMonitor( ...
            poolIn.separation(:, wIdx), poolIn.coastCovariance(:, :, wIdx), ...
            kfIncrement, kfMeas.postCov, propTel, CST_spfParam.NOMINAL_BIAS);

        poolOut.separation(:, wIdx)         = newSeparation;
        poolOut.coastCovariance(:, :, wIdx) = newCoastCov;

        if (ssResult.anyAlarm)
            report.ssAlarm         = true;
            report.anyAlarm        = true;
            poolOut.hadAlarm(wIdx) = true;
            report.alarmPerAxis    = report.alarmPerAxis | ssResult.alarmPerAxis;
        end

        if (ssResult.maxProtectionLevel > report.maxProtectionLevel)
            report.maxProtectionLevel = ssResult.maxProtectionLevel;
        end

        for idx = 1:numAxes                                    % SS margin |d| / (k_FA sigma_SS + B)
            if (ssResult.sigmaSeparation(idx) > 0.0)
                ratio = abs(ssResult.separation(idx)) / (CST_spfParam.K_FALSE_ALERT * ssResult.sigmaSeparation(idx) ...
                    + CST_spfParam.NOMINAL_BIAS(monitoredAxes(idx)));
                if (ratio > report.ssRatio(idx))
                    report.ssRatio(idx) = ratio;
                end
            end
        end

        % 3. window complete: CPI verdict per axis (Eq. 33, 35), then close
        if (age >= windowLength)
            for idx = 1:numAxes
                axisIdx = monitoredAxes(idx);
                [cpiAlarm, qStatistic, ~, cpiFault] = SPF_cpiMonitor( ...
                    poolOut.innovationBuffer(:, 1:windowLength, wIdx), ...
                    poolOut.innovationCovBuffer(:, :, 1:windowLength, wIdx), ...
                    poolOut.obsMatrixBuffer(:, :, 1:windowLength, wIdx), ...
                    poolOut.numMeasBuffer(1:windowLength, wIdx), ...
                    axisIdx, CST_spfParam.NOMINAL_BIAS(axisIdx));

                if (cpiFault)
                    report.solveFault = true;
                end
                ratio = qStatistic / CST_spfParam.CPI_THRESHOLD;   % CPI margin q / T_N
                if (ratio > report.cpiRatio(idx))
                    report.cpiRatio(idx) = ratio;
                end

                if (cpiAlarm)
                    report.cpiAlarm          = true;
                    report.anyAlarm          = true;
                    poolOut.hadAlarm(wIdx)   = true;
                    report.alarmPerAxis(idx) = true;
                end
            end

            if (~poolOut.hadAlarm(wIdx))                       % alarm-free life: the anchor candidate
                report.cleanCloseFound      = true;
                report.cleanCloseSeparation = poolOut.separation(:, wIdx);
                report.cleanCloseCovar      = poolOut.coastCovariance(:, :, wIdx);
            end

            poolOut.active(wIdx) = false;
        end
    end
end

end
%------------------------------------------------------------------------
