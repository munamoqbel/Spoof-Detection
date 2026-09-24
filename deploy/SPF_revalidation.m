%******************************************************************************************
% DESCRIPTION:
% Chi-square consistency test of the raw GNSS rows against the clean coast:
% residual = innovation - H*coastMinusPrior, covariance H*P_C*H' + R.
%
% INPUTS:  innovation [MAX_MEAS x 1], obsMatrix [MAX_MEAS x n], measNoiseCov [MAX_MEAS x MAX_MEAS],
%          numMeas, coastMinusPrior [n x 1], coastCovariance [n x n]
% OUTPUTS: passed (q below the table threshold and numMeas >= REVAL_MIN_MEAS),
%          qValue, solveFault (residual covariance unusable)
%******************************************************************************************
%#codegen
function [passed, qValue, solveFault] = SPF_revalidation(innovation, obsMatrix, ...
    measNoiseCov, numMeas, coastMinusPrior, coastCovariance)

passed     = false;
qValue     = 0.0;
solveFault = false;
maxMeas    = double(size(measNoiseCov, 1));
numMeas    = min(double(numMeas), maxMeas);

if (numMeas > 0)
    tableIdx  = min(numMeas, numel(CST_spfParam.REVAL_THRESHOLD_TABLE));
    threshold = CST_spfParam.REVAL_THRESHOLD_TABLE(tableIdx);

    residual    = innovation - obsMatrix * coastMinusPrior;
    residualCov = obsMatrix * coastCovariance * obsMatrix' + measNoiseCov;
    residualCov = (residualCov + residualCov') / 2;

    % fill the padding diagonal so matrixInv sees a non-singular matrix
    padValue = max(max(diag(residualCov)), 1.0);
    for rowIdx = (numMeas + 1):maxMeas
        residualCov(rowIdx, rowIdx) = padValue;
    end

    varianceOk = true;
    for rowIdx = 1:numMeas
        if ~(residualCov(rowIdx, rowIdx) > 0.0)
            varianceOk = false;
        end
    end
    [sInverse, invInvalid] = matrixInv(residualCov);
    if (varianceOk) && (~invInvalid)
        qValue = residual' * (sInverse * residual);
        if isfinite(qValue) && (qValue >= 0.0) && (qValue < threshold) ...
                && (numMeas >= double(CST_spfParam.REVAL_MIN_MEAS))
            passed = true;
        end
    else
        solveFault = true;
    end
end

end
%------------------------------------------------------------------------
