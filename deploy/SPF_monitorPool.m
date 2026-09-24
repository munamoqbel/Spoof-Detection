%******************************************************************************************
% DESCRIPTION:
% One epoch of the overlapping-window monitor bank: buffer the epoch in every open
% window, SS test each, CPI test per axis on the windows that reach full length, close them.
%
% INPUTS:  poolIn (STRUCT_SPF.setMonitorPool), kfMeas, propTel
% OUTPUTS: poolOut, report (STRUCT_SPF.setMonitorReport)
%******************************************************************************************
%#codegen
function [poolOut, report] = SPF_monitorPool(poolIn, kfMeas, propTel)

windowLength  = CST_spfParam.WINDOW_LENGTH;
monitoredAxes = CST_spfParam.MONITORED_AXES;
numAxes       = cast(numel(monitoredAxes), 'uint8');
poolOut       = poolIn;
kfIncrement   = kfMeas.postState - kfMeas.priorState;       % K*y of the watched filter

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

        % 2. solution-separation test
        [newSeparation, newCoastCov, ssResult] = SPF_ssMonitor( ...
            poolIn.separation(:, wIdx), poolIn.coastCovariance(:, :, wIdx), ...
            kfIncrement, kfMeas.postCov, propTel);

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

        for idx = 1:numAxes                                    % SS margin telemetry
            if (ssResult.sigmaSeparation(idx) > 0.0)
                ratio = abs(ssResult.separation(idx)) / (CST_spfParam.K_FALSE_ALERT * ssResult.sigmaSeparation(idx));
                if (ratio > report.ssRatio(idx))
                    report.ssRatio(idx) = ratio;
                end
            end
        end

        % 3. window complete: CPI verdict per axis, then close
        if (age >= windowLength)
            for idx = 1:numAxes
                axisIdx = monitoredAxes(idx);
                [cpiAlarm, qStatistic, ~, cpiFault] = SPF_cpiMonitor( ...
                    poolOut.innovationBuffer(:, 1:windowLength, wIdx), ...
                    poolOut.innovationCovBuffer(:, :, 1:windowLength, wIdx), ...
                    poolOut.obsMatrixBuffer(:, :, 1:windowLength, wIdx), ...
                    poolOut.numMeasBuffer(1:windowLength, wIdx), ...
                    axisIdx);

                if (cpiFault)
                    report.solveFault = true;
                end
                ratio = qStatistic / CST_spfParam.CPI_THRESHOLD;   % CPI margin telemetry
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

            if (~poolOut.hadAlarm(wIdx))                       % alarm-free life: anchor candidate
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
