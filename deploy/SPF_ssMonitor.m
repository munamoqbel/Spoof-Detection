%***************************************************************************
% DESCRIPTION:
% Solution-separation test of one window for one epoch on all monitored
% axes, in increment form. The window's coast is the solution at the
% window open propagated INS-only; the separation between the host
% solution and that coast is the accumulation of the host's update
% increments since the open,
%   d_k = Phi d_(k-1) + (x+ - x_bar)        P_C = Phi P_C Phi' + Q
% so no absolute state is needed.
%
% INPUTS:
%   - separation       [n x 1] d_(k-1) of this window
%   - coastCovariance  [n x n] P_C of the window's coast, previous epoch
%   - kfIncrement      [n x 1] x+ - x_bar of the watched filter
%   - kfCovariance     [n x n] P+ of the watched filter
%   - propTel          .accumPhi, .accumQ interval Phi and Q
%   - nominalBias      [1 x 3] B per NED position state (m)
%
% OUTPUTS:
%   - separation       d_k
%   - coastCovariance  P_C propagated to this epoch
%   - ssResult         STRUCT_SPF.setSSmonitor: alarm per axis, the
%                      separation, sigma_SS and the protection level per
%                      axis, and the largest protection level
%
% ASSUMPTIONS AND LIMITATIONS:
% An axis with sigma_SS = 0 (P_C not above P+) is untestable and reports
% no alarm. k_FA and k_MD are CST_spfParam.K_FALSE_ALERT and
% K_MISSED_DETECTION.
%
% REQUIREMENT TRACEABILITY:
% Paper (Kujur 2024):
%   E28, E31  coast propagation of the increment and of its covariance
%   Eq. 49    separation q_SS = u' (x_KF - x_C) = u' d
%   Eq. 50    sigma_SS^2 = u' (P_C - P_KF) u
%   Eq. 52    sigma_C^2 = u' P_C u
%   Eq. 51    PL = k_FA sigma_SS + B + k_MD sigma_C
%   alarm when |q_SS| > k_FA sigma_SS + B (B = 0 in the paper)
%***************************************************************************
%#codegen
function [separation, coastCovariance, ssResult] = SPF_ssMonitor...
    (separation, coastCovariance, kfIncrement, kfCovariance, propTel, nominalBias)

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
separation      = Phi * separation + kfIncrement;          % E28 in increment form
coastCovariance = Phi * coastCovariance * Phi' + Q;        % E31

for idx = 1:numAxes

    axisIdx = monitoredAxes(idx);

    sepAxis(idx)       = separation(axisIdx);                                                  % Eq. 49
    varianceSeparation = coastCovariance(axisIdx, axisIdx) - kfCovariance(axisIdx, axisIdx);   % Eq. 50
    varianceCoast      = coastCovariance(axisIdx, axisIdx);                                    % Eq. 52

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

    gate = kFalseAlert * sigmaSeparation(idx) + nominalBias(axisIdx);   % k_FA sigma_SS + B
    protectionLevel(idx) = gate + kMissedDetection * sigmaCoast;        % Eq. 51

    if testable && (abs(sepAxis(idx)) > gate)                           % alarm rule
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
