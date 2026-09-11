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
%   Window length N and threshold T_N are CST_spfParam.WINDOW_LENGTH and
%   CST_spfParam.CPI_THRESHOLD; the buffers are N deep.
%
% OUTPUTS:
%   - cpiAlarm                  logical
%   - qStatistic                double Eq. 33 statistic
%   - xiHistory                 [N x 1] normalised projections (diagnostics)
%   - solveFault                logical, true if matrixInv flagged S of at
%                               least one epoch as unusable (that epoch
%                               counted xi = 0)
%
% ASSUMPTIONS AND LIMITATIONS:
% An epoch with numMeas = 0 (no GNSS) contributes xi = 0 to the window.
% Everything runs on the full MAX_MEAS layout (rows/cols beyond numMeas
% are zero): no variable-size expression, so it compiles with MATLAB
% Coder variable sizing off. S is inverted by the host's matrixInv (SVD
% pseudo-inverse: the padding stays zero, an unusable S is flagged), the
% same routine kfUpdate uses for its own S.
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
function [cpiAlarm, qStatistic, xiHistory, solveFault] = SPF_cpiMonitor...
    (innovationBuffer, innovationCovBuffer, obsMatrixBuffer, numMeasBuffer, axisIdx)

% Define variables
windowLength = CST_spfParam.WINDOW_LENGTH;
cpiThreshold = CST_spfParam.CPI_THRESHOLD;
maxMeas      = uint8(size(innovationBuffer, 1));
xiHistory    = zeros(windowLength, 1);
qStatistic   = 0.0;
cpiAlarm     = false;
solveFault   = false;

for idx = 1:windowLength
    numMeas = min(numMeasBuffer(idx), maxMeas);       % defensive: never past the buffer
    xiNormalised = 0.0;

    if (numMeas > 0)
        % full fixed-size layout; rows/cols beyond numMeas are zero padding
        innovation    = innovationBuffer(:, idx);                % gamma. Eq.3
        innovationCov = innovationCovBuffer(:, :, idx);          % S. Defined under Eq. 4
        projection    = obsMatrixBuffer(:, axisIdx, idx);        % f = H(:,axis). From Eq. 17

        % the host's matrixInv flags a singular input, and zero padding is
        % singular: fill the padding diagonal with the largest live variance
        % (blkdiag(S, p*I) inverts to blkdiag(S^-1, I/p); the padded rows of
        % gamma and f are zero, so the projections are unchanged)
        padValue = max(max(diag(innovationCov)), 1.0);
        for rowIdx = (double(numMeas) + 1):double(maxMeas)
            innovationCov(rowIdx, rowIdx) = padValue;
        end

        % Exception handler: a covariance with a non-positive variance on a
        % live row is unusable; the host's SVD inverse then handles the
        % padding (pseudo-inverse) and flags a non-finite S. Either case
        % drops the epoch (xi = 0) and reports solveFault.
        varianceOk = true;
        for rowIdx = 1:double(numMeas)
            if ~(innovationCov(rowIdx, rowIdx) > 0.0)
                varianceOk = false;
            end % ELSE is trivial
        end
        [sInverse, invInvalid] = matrixInv(innovationCov);
        if (varianceOk) && (~invInvalid)
            sInvInnovation = sInverse * innovation;              % S^-1 gamma
            sInvProjection = sInverse * projection;              % S^-1 f

            gammaProjection  = 0.0;  % Eq. 17
            sigma2Projection = 0.0;  % Eq. 20
            for idxMeas = 1:double(maxMeas)
                gammaProjection  = gammaProjection + projection(idxMeas) * sInvInnovation(idxMeas);   % = f' * S^-1 gamma (Eq.17); loop used for deterministic rounding
                sigma2Projection = sigma2Projection + projection(idxMeas) * sInvProjection(idxMeas);  % = f' * S^-1 f     (Eq.20); loop used for deterministic rounding
            end

            % Exception handler: axis unobservable this epoch (f = 0) -> xi = 0
            if (sigma2Projection > 0.0) && isfinite(sigma2Projection) && isfinite(gammaProjection)
                xiNormalised = gammaProjection / sqrt(sigma2Projection);   % Eq. 29
            end % ELSE is trivial
        else
            solveFault = true;
        end
    end % ELSE: no GNSS this epoch -> xi = 0

    xiHistory(idx) = xiNormalised;
    qStatistic     = qStatistic + (xiNormalised * xiNormalised);   % Eq. 33
end

if (qStatistic > cpiThreshold)    % alarm rule (section 2 under Eq. 4)
    cpiAlarm = true;
end % ELSE is trivial

end

%------------------------------------------------------------------------
