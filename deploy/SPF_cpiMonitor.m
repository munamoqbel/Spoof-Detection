%******************************************************************************************
% DESCRIPTION:
% Cumulative position-domain innovation (CPI) test of one complete window, one axis.
% Fixed MAX_MEAS layout: rows/cols beyond numMeas are zero padding.
%
% INPUTS:  innovationBuffer [MAX_MEAS x N], innovationCovBuffer [MAX_MEAS x MAX_MEAS x N],
%          obsMatrixBuffer [MAX_MEAS x n x N], numMeasBuffer [N x 1], axisIdx
% OUTPUTS: cpiAlarm, qStatistic, xiHistory [N x 1], solveFault (an S was unusable)
%******************************************************************************************
%#codegen
function [cpiAlarm, qStatistic, xiHistory, solveFault] = SPF_cpiMonitor...
    (innovationBuffer, innovationCovBuffer, obsMatrixBuffer, numMeasBuffer, axisIdx)

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
        innovation    = innovationBuffer(:, idx);
        innovationCov = innovationCovBuffer(:, :, idx);
        projection    = obsMatrixBuffer(:, axisIdx, idx);        % f = H(:, axis)

        % fill the padding diagonal so matrixInv sees a non-singular matrix
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

            gammaProjection  = 0.0;                              % f' S^-1 gamma
            sigma2Projection = 0.0;                              % f' S^-1 f
            for idxMeas = 1:double(maxMeas)
                gammaProjection  = gammaProjection + projection(idxMeas) * sInvInnovation(idxMeas);
                sigma2Projection = sigma2Projection + projection(idxMeas) * sInvProjection(idxMeas);
            end

            if (sigma2Projection > 0.0) && isfinite(sigma2Projection) && isfinite(gammaProjection)
                xiNormalised = gammaProjection / sqrt(sigma2Projection);
            end
        else
            solveFault = true;
        end
    end

    xiHistory(idx) = xiNormalised;
    qStatistic     = qStatistic + (xiNormalised * xiNormalised);
end

if (qStatistic > cpiThreshold)
    cpiAlarm = true;
end

end
%------------------------------------------------------------------------
