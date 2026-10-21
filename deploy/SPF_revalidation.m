%***************************************************************************
% DESCRIPTION:
% Chi-square consistency test of the raw GNSS rows against the clean
% coast: the end-of-attack detector, a decision in COAST and a
% diagnostic in PROBATION. Formed on the host's own innovation so that a
% nonlinear h(x) is handled by the host's linearisation:
%   residual   = z - h(x_coast) ~ innovation - H (x_coast - x_prior)
%   covariance = H P_C H' + R
% One pass means nothing on its own: the caller needs a dwell of
% consecutive passes, and even then only opens PROBATION.
%
% INPUTS:
%   - innovation       [MAX_MEAS x 1] host innovation z - h(x_prior)
%   - obsMatrix        [MAX_MEAS x n] H of the host's linearisation
%   - measNoiseCov     [MAX_MEAS x MAX_MEAS] R
%   - numMeas          valid rows this epoch (0 = no measurement)
%   - coastMinusPrior  [n x 1] x_coast - x_prior in increment space, zero
%                      while the operational KF is not updated
%   - coastCovariance  [n x n] P_C of the coast
%
% OUTPUTS:
%   - passed      0 <= qValue < threshold(numMeas) and numMeas >=
%                 REVAL_MIN_MEAS; false with no row or an unusable
%                 covariance
%   - qValue      the chi-square statistic, 0 when it could not be formed
%   - solveFault  true when matrixInv flagged H P_C H' + R unusable
%
% ASSUMPTIONS AND LIMITATIONS:
% Threshold from CST_spfParam.REVAL_THRESHOLD_TABLE(numMeas), the index
% clamped to the table length. The rows may include non-GNSS rows (the
% pressure altitude); REVAL_MIN_MEAS keeps an epoch with too few rows
% from counting as a pass. Everything runs on the full MAX_MEAS layout
% (rows and columns beyond numMeas are zero, their diagonal padded
% before the inverse): no variable-size expression. A negative or
% non-finite q cannot pass.
%
% REQUIREMENT TRACEABILITY:
% Implementation Guide Sec. 5: one chi-square gate per possible row
% count. Paper Sec. 5 recovery: the re-validation that precedes the
% probation.
%***************************************************************************
%#codegen
function [passed, qValue, solveFault] = SPF_revalidation(innovation, obsMatrix, ...
    measNoiseCov, numMeas, coastMinusPrior, coastCovariance)

passed     = false;
qValue     = 0.0;
solveFault = false;
maxMeas    = double(size(measNoiseCov, 1));
numMeas    = min(double(numMeas), maxMeas);

if (numMeas > 0)
    tableIdx  = min(numMeas, numel(CST_spfParam.REVAL_THRESHOLD_TABLE));   % one gate per row count
    threshold = CST_spfParam.REVAL_THRESHOLD_TABLE(tableIdx);

    residual    = innovation - obsMatrix * coastMinusPrior;                 % z - h(x_coast)
    residualCov = obsMatrix * coastCovariance * obsMatrix' + measNoiseCov;  % H P_C H' + R
    residualCov = (residualCov + residualCov') / 2;

    % padding diagonal filled so that matrixInv sees a non-singular matrix
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
        qValue = residual' * (sInverse * residual);                         % the chi-square statistic
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
