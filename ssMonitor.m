%******************************************************************************************
% DESCRIPTION:
% Solution-Separation test: one window, one epoch, ALL monitored axes.
%
% INCREMENT FORM. The window's coast is the solution at window open
% propagated INS-only. Because the host solution and the coast share the
% same propagation, their difference is the accumulation of the host KF's
% update increments since the window opened:
%     d_k = Phi_acc * d_(k-1) + (x+_k - x_bar_k)          (paper E28 / Eq. 49)
%     P_C = Phi_acc * P_C * Phi_acc' + Q_acc               (paper E31)
% so no absolute state is needed and the result is invariant to how the
% host splits its estimate between feedback and residual states.
%
% INPUTS:
%   - separation           [n x 1]  d_(k-1) for this window
%   - coastCovariance      [n x n]  P_C of this window's coast (previous epoch)
%   - kfIncrement          [n x 1]  x+ - x_bar of the watched filter this epoch
%   - kfCovariance         [n x n]  P+ of the watched filter this epoch
%   - propTel              .accumPhi / .accumQ  interval Phi / Q
%   - kFalseAlert, kMissedDetection (OPTIONAL, tests only) override
%                          the CST_spfParam constants
%
% OUTPUTS:
%   - separation           d_k
%   - coastCovariance      P_C propagated to this epoch
%   - ssResult             structure of type STRUCT_SPF.setSSmonitor
%
% ASSUMPTIONS AND LIMITATIONS:
% An undefined test (sigma_SS = 0) reports no alarm.
%
% REQUIREMENT TRACEABILITY:
% - PAPER MAPPING
%    coast propagation:     E28 (state, in increment form), E31 (covariance)
%    separation:            Eq. 49  q_SS = u'(x_KF - x_C) = u' d
%    separation variance:   Eq. 50  sigma_SS^2 = u'(P_C - P_KF)u
%    coast variance:        Eq. 52  sigma_C^2  = u' P_C u
%    protection level:      Eq. 51  PL = k_FA*sigma_SS + k_MD*sigma_C
%    alarm rule:            |q_SS| > k_FA * sigma_SS
%
%******************************************************************************************
%#codegen
function [separation, coastCovariance, ssResult] = ssMonitor...
    (separation, coastCovariance, kfIncrement, kfCovariance, propTel, ...
     kFalseAlert, kMissedDetection)

% Define variables
monitoredAxes      = CST_spfParam.MONITORED_AXES;
if (nargin < 7)
    kFalseAlert      = CST_spfParam.K_FALSE_ALERT;
    kMissedDetection = CST_spfParam.K_MISSED_DETECTION;
end
numAxes            = cast(numel(monitoredAxes), 'uint8');
sepAxis            = zeros(1, numAxes);
sigmaSeparation    = zeros(1, numAxes);
protectionLevel    = zeros(1, numAxes);
alarmPerAxis       = zeros(1, numAxes, 'logical');
anyAlarm           = false;
maxProtectionLevel = 0.0;

% Propagate the coast once (shared by all axes)
Phi = propTel.accumPhi;
Q   = propTel.accumQ;
separation      = Phi * separation + kfIncrement;                    % E28 (increment form)
coastCovariance = Phi * coastCovariance * Phi' + Q;                  % E31

% Scalar test per monitored axis
for idx = 1:numAxes

    axisIdx = monitoredAxes(idx);    % actual state index

    sepAxis(idx)       = separation(axisIdx);                            % Eq. 49
    varianceSeparation = coastCovariance(axisIdx, axisIdx) ...
        - kfCovariance(axisIdx, axisIdx);                                % Eq. 50
    varianceCoast      = coastCovariance(axisIdx, axisIdx);              % Eq. 52

    % Exception handler: sigma_SS = 0 (window just opened, or numerical)
    % means the test is undefined; report no alarm rather than |q| > 0.
    testable = false;
    if (varianceSeparation > 0.0)
        sigmaSeparation(idx) = sqrt(varianceSeparation);
        testable = true;
    end % ELSE is trivial

    if (varianceCoast > 0.0)
        sigmaCoast = sqrt(varianceCoast);
    else % Exception handler
        sigmaCoast = 0.0;
    end

    protectionLevel(idx) = kFalseAlert * sigmaSeparation(idx) ...
        + kMissedDetection * sigmaCoast;                                 % Eq. 51

    % alarm rule (section 5 under Eq. 52)
    if testable && (abs(sepAxis(idx)) > (kFalseAlert * sigmaSeparation(idx)))
        alarmPerAxis(idx) = true;
        anyAlarm = true;
    end % ELSE is trivial

    if (protectionLevel(idx) > maxProtectionLevel)
        maxProtectionLevel = protectionLevel(idx);
    end
end

ssResult = STRUCT_SPF.setSSmonitor(alarmPerAxis, anyAlarm, ...
    sepAxis, sigmaSeparation, protectionLevel, maxProtectionLevel);

end

%------------------------------------------------------------------------
