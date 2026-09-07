%******************************************************************************************
% DESCRIPTION:
% Cumulative Position-domain Innovation monitor, one complete window.
%
% INPUTS:
%   - innovation_buffer          [m x N] gamma history of the window
%   - innovation_cov_buffer      [m x m x N] S history
%   - obs_matrix_buffer          [m x n x N] H history
%   - axis_idx                   monitored position state index
%
% OUTPUTS:
%   - is_alarm                   logical
%   - q_statistic                double Eq. 33 statistic
%   - xi_history                 [N x 1] normalised projections (diagnostics)
%
% ASSUMPTIONS AND LIMITATIONS:
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
    (innovationBuffer, innovationCovBuffer, obsMatrixBuffer, axisIdx)

% Define variables
windowLength = CST_spfParam.WINDOW_LENGTH;
cpiThreshold = CST_spfParam.CPI_THRESHOLD;
numMeas      = size(innovationBuffer, 1); % To be replaced by numMeas from KFL_IS22_stateUpdates
xiHistory    = zeros(windowLength, 1);
qStatistic   = 0.0;
cpiAlarm     = false;

for idx = 1:windowLength
    innovation    = innovationBuffer(:, idx);           % gamma. Eq.3
    innovationCov = innovationCovBuffer(:, :, idx);     % S. Defined under Eq. 4
    projection    = obsMatrixBuffer(:, axisIdx, idx);   % f = H(:,axis). From Eq. 17

    % ----------------------------------------------------------------------
    % WARNING: Exception handler need to be added!
    % ----------------------------------------------------------------------
    sInvInnovation = innovationCov \ innovation;     % S^{-1} gamma
    sInvProjection = innovationCov \ projection;     % S^{-1} f
    % ----------------------------------------------------------------------

    gammaProjection  = 0.0;  % Eq. 17
    sigma2Projection = 0.0;  % Eq. 20
    for idxMeas = 1:numMeas
        gammaProjection  = gammaProjection + projection(idxMeas) * sInvInnovation(idxMeas);   % = projection' * sInvInnovation (Eq.17); loop used for deterministic rounding
        sigma2Projection = sigma2Projection + projection(idxMeas) * sInvProjection(idxMeas);  % = projection' * sInvProjection (Eq.20); loop used for deterministic rounding
    end

    % Exception handler
    if (sigma2Projection > 0.0)
        xiNormalised = gammaProjection / sqrt(sigma2Projection);   % Eq. 29
    else
        xiNormalised = 0.0; % axis unobservable this epoch
    end

    xiHistory(idx) = xiNormalised;
    qStatistic     = qStatistic + (xiNormalised * xiNormalised);   % Eq. 33
end

if (qStatistic > cpiThreshold)    % alarm rule (section 2 under Eq. 4)
    cpiAlarm = true;
end % ELSE is trivial

end

%------------------------------------------------------------------------------------------
