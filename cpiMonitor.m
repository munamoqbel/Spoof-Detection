%******************************************************************************************
% DESCRIPTION:
% Cumulative Position-domain Innovation monitor, one complete window.
%
% INPUTS:
%   - innovationBuffer          [MAX_MEAS x N] gamma history of the window
%   - innovationCovBuffer       [MAX_MEAS x MAX_MEAS x N] S history
%   - obsMatrixBuffer           [MAX_MEAS x n x N] H history
%   - numMeasBuffer             [N x 1] valid measurement count per epoch
%                               (rows/cols beyond it are padding)
%   - axisIdx                   monitored position state index
%   - windowLength, cpiThreshold (OPTIONAL, tests only) override the
%                               CST_spfParam constants
%
% OUTPUTS:
%   - cpiAlarm                  logical
%   - qStatistic                double Eq. 33 statistic
%   - xiHistory                 [N x 1] normalised projections (diagnostics)
%
% ASSUMPTIONS AND LIMITATIONS:
% An epoch with numMeas = 0 (no GNSS) contributes xi = 0 to the window.
%
% REQUIREMENT TRACEABILITY:
% - PAPER MAPPING
%    Equations (Kujur 2024)
%    Eq. 17: gamma_u = f'*(S\gamma), f = H(:,axis_idx) = (H * u)
%    Eq. 20: sigma2 = f'*(S\f)
%    Eq. 29: xi = gamma_u / sqrt(sigma2) ~ N(0,1) under H0
%    Eq. 33: q = sum(xi.^2) ~ Gamma(N/2,2) under H0
%    Eq. 35: alarm iff q > cpi_threshold (T_N, precomputed offline)
%
%******************************************************************************************
%#codegen
function [cpiAlarm, qStatistic, xiHistory] = cpiMonitor...
    (innovationBuffer, innovationCovBuffer, obsMatrixBuffer, numMeasBuffer, axisIdx, ...
     windowLength, cpiThreshold)

% Define variables
if (nargin < 7)
    windowLength = CST_spfParam.WINDOW_LENGTH;
    cpiThreshold = CST_spfParam.CPI_THRESHOLD;
end
xiHistory    = zeros(windowLength, 1);
qStatistic   = 0.0;
cpiAlarm     = false;

for idx = 1:windowLength
    numMeas = numMeasBuffer(idx);
    xiNormalised = 0.0;

    if (numMeas > 0)
        innovation    = innovationBuffer(1:numMeas, idx);               % gamma. Eq.3
        innovationCov = innovationCovBuffer(1:numMeas, 1:numMeas, idx); % S. Defined under Eq. 4
        projection    = obsMatrixBuffer(1:numMeas, axisIdx, idx);       % f = H(:,axis). From Eq. 17

        % ------------------------------------------------------------------------
        % WARNING: Exception handler need to be added!
        % ------------------------------------------------------------------------
        sInvInnovation = innovationCov \ innovation;     % S^{-1} gamma
        sInvProjection = innovationCov \ projection;     % S^{-1} f
        % ------------------------------------------------------------------------

        gammaProjection  = 0.0;  % Eq. 17
        sigma2Projection = 0.0;  % Eq. 20
        for idxMeas = 1:numMeas
            gammaProjection  = gammaProjection + projection(idxMeas) * sInvInnovation(idxMeas);   % = projection' * sInvInnovation (Eq.17); loop used for deterministic rounding
            sigma2Projection = sigma2Projection + projection(idxMeas) * sInvProjection(idxMeas);  % = projection' * sInvProjection (Eq.20); loop used for deterministic rounding
        end

        % Exception handler
        if (sigma2Projection > 0.0)
            xiNormalised = gammaProjection / sqrt(sigma2Projection);   % Eq. 29
        end % ELSE: axis unobservable this epoch -> xi = 0
    end % ELSE: no GNSS this epoch -> xi = 0

    xiHistory(idx) = xiNormalised;
    qStatistic     = qStatistic + (xiNormalised * xiNormalised);   % Eq. 33
end

if (qStatistic > cpiThreshold)    % alarm rule (section 2 under Eq. 4)
    cpiAlarm = true;
end % ELSE is trivial

end

%------------------------------------------------------------------------
