%***************************************************************************
% DESCRIPTION:
% Cumulative position-domain innovation (CPI) test of one complete
% window on one monitored axis: the innovation of each of the N epochs
% projected on the axis and normalised, the squares summed and compared
% with T_N.
%
% INPUTS:
%   - innovationBuffer     [MAX_MEAS x N] innovation per epoch
%   - innovationCovBuffer  [MAX_MEAS x MAX_MEAS x N] S per epoch
%   - obsMatrixBuffer      [MAX_MEAS x n x N] H per epoch
%   - numMeasBuffer        [N x 1] valid rows per epoch; rows and columns
%                          beyond it are zero padding
%   - axisIdx              position state index of the axis
%   - nominalBias          B of this axis (m): deadband on the
%                          position-domain innovation
%
% OUTPUTS:
%   - cpiAlarm     q > T_N
%   - qStatistic   q of Eq. 33
%   - xiHistory    [N x 1] the normalised projections (diagnostics)
%   - solveFault   true when matrixInv flagged the S of an epoch
%                  unusable; that epoch counts xi = 0
%
% ASSUMPTIONS AND LIMITATIONS:
% An epoch without rows contributes xi = 0. Everything runs on the full
% MAX_MEAS layout (padding diagonal filled before the inverse): no
% variable-size expression. S is inverted by the host's matrixInv, the
% routine its own update uses. N and T_N are CST_spfParam.WINDOW_LENGTH
% and CPI_THRESHOLD.
%
% REQUIREMENT TRACEABILITY:
% Paper (Kujur 2024):
%   Eq. 3   gamma, the innovation; S is defined under Eq. 4
%   Eq. 17  gamma_u = f' S^-1 gamma, f = H(:, axis)
%   Eq. 20  sigma_u^2 = f' S^-1 f
%   Eq. 29  xi = gamma_u / sigma_u ~ N(0, 1) under H0; with the
%           allowance B: xi = max(|gamma_u| / sigma_u - B sigma_u, 0)
%   Eq. 33  q = sum of xi^2 over the window
%   Eq. 35  alarm when q > T_N
%***************************************************************************
%#codegen
function [cpiAlarm, qStatistic, xiHistory, solveFault] = SPF_cpiMonitor...
    (innovationBuffer, innovationCovBuffer, obsMatrixBuffer, numMeasBuffer, axisIdx, nominalBias)

windowLength = CST_spfParam.WINDOW_LENGTH;
cpiThreshold = CST_spfParam.CPI_THRESHOLD;
maxMeas      = uint8(size(innovationBuffer, 1));
xiHistory    = zeros(windowLength, 1);
qStatistic   = 0.0;
cpiAlarm     = false;
solveFault   = false;

for idx = 1:windowLength
    numMeas = min(numMeasBuffer(idx), maxMeas);
    xiNormalised = 0.0;

    if (numMeas > 0)
        innovation    = innovationBuffer(:, idx);                % gamma (Eq. 3)
        innovationCov = innovationCovBuffer(:, :, idx);          % S (under Eq. 4)
        projection    = obsMatrixBuffer(:, axisIdx, idx);        % f = H(:, axis) (Eq. 17)

        % padding diagonal filled so that matrixInv sees a non-singular matrix
        padValue = max(max(diag(innovationCov)), 1.0);
        for rowIdx = (double(numMeas) + 1):double(maxMeas)
            innovationCov(rowIdx, rowIdx) = padValue;
        end

        varianceOk = true;
        for rowIdx = 1:double(numMeas)
            if ~(innovationCov(rowIdx, rowIdx) > 0.0)
                varianceOk = false;
            end
        end
        [sInverse, invInvalid] = matrixInv(innovationCov);
        if (varianceOk) && (~invInvalid)
            sInvInnovation = sInverse * innovation;
            sInvProjection = sInverse * projection;

            gammaProjection  = 0.0;                              % f' S^-1 gamma (Eq. 17)
            sigma2Projection = 0.0;                              % f' S^-1 f (Eq. 20)
            for idxMeas = 1:double(maxMeas)                      % element loop: deterministic rounding
                gammaProjection  = gammaProjection + projection(idxMeas) * sInvInnovation(idxMeas);
                sigma2Projection = sigma2Projection + projection(idxMeas) * sInvProjection(idxMeas);
            end

            % Eq. 29, the part of the projection beyond the deadband B
            if (sigma2Projection > 0.0) && isfinite(sigma2Projection) && isfinite(gammaProjection)
                xiNormalised = max(abs(gammaProjection) / sqrt(sigma2Projection) ...
                    - nominalBias * sqrt(sigma2Projection), 0.0);
            end
        else
            solveFault = true;
        end
    end

    xiHistory(idx) = xiNormalised;
    qStatistic     = qStatistic + (xiNormalised * xiNormalised);   % Eq. 33
end

if (qStatistic > cpiThreshold)                                     % Eq. 35
    cpiAlarm = true;
end

end
%------------------------------------------------------------------------
