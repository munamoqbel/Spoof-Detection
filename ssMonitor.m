%******************************************************************************************
% DESCRIPTION:
% Solution-Separation test: one window, one epoch, ALL monitored axes.
%
% The coast is propagated ONCE, then the scalar test is evaluated for each
% axis.
%
% INPUTS:
%   - coast_state          [n x 1]  this window's coast (INS-only) state
%   - coast_covariance     [n x n]  its covariance
%   - kf_state             [n x 1]  current filter state (post-update)
%   - kf_covariance        [n x n]  current filter covariance (post-update)
%   - Phi, Q               [n x n]  state transition and process noise
%
% OUTPUTS:
%   - coast_state          propagated to this epoch
%   - coast_covariance     propagated to this epoch
%   - ssResult             structure of type STRUCT_SPF.setSSmonitor
%
% ASSUMPTIONS AND LIMITATIONS:
%
% REQUIREMENT TRACEABILITY:
% - PAPER MAPPING
%    coast propagation:     E28 (state), E31 (covariance)
%    separation:            Eq. 49  q_SS = u'(x_KF - x_C)
%    separation variance:   Eq. 50  sigma_SS^2 = u'(P_C - P_KF)u
%    coast variance:        Eq. 52  sigma_C^2  = u' P_C u
%    protection level:      Eq. 51  PL = k_FA*sigma_SS + k_MD*sigma_C
%    alarm rule:            |q_SS| > k_FA * sigma_SS
%
%******************************************************************************************
%#codegen
function [coastState, coastCovariance, ssResult] = ssMonitor...
    (coastState, coastCovariance, kfState, kfCovariance, spoofInfo)

% Define variables
monitoredAxes      = CST_spfParam.MONITORED_AXES;
kFalseAlert        = CST_spfParam.K_FALSE_ALERT;
kMissedDetection   = CST_spfParam.K_MISSED_DETECTION;
numAxes            = cast(numel(monitoredAxes), 'uint8');
separation         = zeros(1, numAxes);
sigmaSeparation    = zeros(1, numAxes);
protectionLevel    = zeros(1, numAxes);
alarmPerAxis       = zeros(1, numAxes, 'logical');
anyAlarm           = false;
maxProtectionLevel = 0.0;

% Propagate the coast once (shared by all axes)
Phi = spoofInfo.phiAcc;
Q   = spoofInfo.qAcc;
coastState      = Phi * coastState;                                  % E28
coastCovariance = Phi * coastCovariance * Phi' + Q;                  % E31

% Scalar test per monitored axis
for idx = 1:numAxes

    axisIdx = monitoredAxes(idx);    % actual state index

    separation(idx)     = kfState(axisIdx) - coastState(axisIdx);    % Eq. 49
    varianceSeparation = coastCovariance(axisIdx, axisIdx) ...
        - kfCovariance(axisIdx, axisIdx);                                % Eq. 50
    varianceCoast       = coastCovariance(axisIdx, axisIdx);             % Eq. 52

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
    if testable && (abs(separation(idx)) > (kFalseAlert * sigmaSeparation(idx)))
        alarmPerAxis(idx) = true;
        anyAlarm = true;
    end % ELSE is trivial

    if (protectionLevel(idx) > maxProtectionLevel)
        maxProtectionLevel = protectionLevel(idx);
    end
end

ssResult = STRUCT_SPF.setSSmonitor(alarmPerAxis, anyAlarm, ...
    separation, sigmaSeparation, protectionLevel, maxProtectionLevel);

end

%------------------------------------------------------------------------------------------
