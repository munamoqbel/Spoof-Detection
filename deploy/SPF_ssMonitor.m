%******************************************************************************************
% DESCRIPTION:
% Solution-separation test of one window for one epoch, all monitored axes (increment
% form: d_k = Phi*d_(k-1) + (x+ - x_bar), P_C = Phi*P_C*Phi' + Q).
%
% INPUTS:  separation [n x 1], coastCovariance [n x n], kfIncrement [n x 1],
%          kfCovariance [n x n] (P+), propTel
% OUTPUTS: separation, coastCovariance (propagated), ssResult (STRUCT_SPF.setSSmonitor)
%******************************************************************************************
%#codegen
function [separation, coastCovariance, ssResult] = SPF_ssMonitor...
    (separation, coastCovariance, kfIncrement, kfCovariance, propTel)

monitoredAxes      = CST_spfParam.MONITORED_AXES;
kFalseAlert        = CST_spfParam.K_FALSE_ALERT;
kMissedDetection   = CST_spfParam.K_MISSED_DETECTION;
numAxes            = cast(numel(monitoredAxes), 'uint8');
sepAxis            = zeros(1, numAxes);
sigmaSeparation    = zeros(1, numAxes);
protectionLevel    = zeros(1, numAxes);
alarmPerAxis       = zeros(1, numAxes, 'logical');
anyAlarm           = false;
maxProtectionLevel = 0.0;

Phi = propTel.accumPhi;
Q   = propTel.accumQ;
separation      = Phi * separation + kfIncrement;
coastCovariance = Phi * coastCovariance * Phi' + Q;

for idx = 1:numAxes

    axisIdx = monitoredAxes(idx);

    sepAxis(idx)       = separation(axisIdx);
    varianceSeparation = coastCovariance(axisIdx, axisIdx) - kfCovariance(axisIdx, axisIdx);
    varianceCoast      = coastCovariance(axisIdx, axisIdx);

    testable = false;                                  % sigma_SS = 0: test undefined, no alarm
    if (varianceSeparation > 0.0)
        sigmaSeparation(idx) = sqrt(varianceSeparation);
        testable = true;
    end

    if (varianceCoast > 0.0)
        sigmaCoast = sqrt(varianceCoast);
    else
        sigmaCoast = 0.0;
    end

    protectionLevel(idx) = kFalseAlert * sigmaSeparation(idx) + kMissedDetection * sigmaCoast;

    if testable && (abs(sepAxis(idx)) > (kFalseAlert * sigmaSeparation(idx)))
        alarmPerAxis(idx) = true;
        anyAlarm = true;
    end

    if (protectionLevel(idx) > maxProtectionLevel)
        maxProtectionLevel = protectionLevel(idx);
    end
end

ssResult = STRUCT_SPF.setSSmonitor(alarmPerAxis, anyAlarm, ...
    sepAxis, sigmaSeparation, protectionLevel, maxProtectionLevel);

end
%------------------------------------------------------------------------
