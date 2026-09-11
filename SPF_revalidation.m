%******************************************************************************************
% DESCRIPTION:
% Check if GNSS is statistically consistent with the clean coast.
%
% Chi-square consistency test of the raw measurements against the
% INS-only coast. This is the end-of-attack detector used during COAST:
% while spoofing is active, the offset (vs mm-level carrier noise) makes
% qValue astronomically large; after the spoofer stops, qValue drops
% to ~chi2(m) and passes. NOTE: one pass means nothing on its own - the
% caller requires a dwell of consecutive passes, and even then only
% opens PROBATION, never a direct handback.
%
% Formulated on the HOST's innovation so that a nonlinear h(x) is handled
% by the host's own linearisation:
%   residual = z - h(x_coast)
%            ~ [z - h(x_prior)] - H (x_coast - x_prior)
%            =  innovation      - H * coastMinusPrior
% coastMinusPrior is supplied by the caller in increment form
% (= -Phi_acc * accumulated increments of the active filter since the
% protected solution was set; zero when the host has not updated it).
%
% INPUTS:
%   - innovation        [MAX_MEAS x 1] host innovation z - h(x_prior)
%   - obsMatrix         [MAX_MEAS x n] H (host linearisation)
%   - measNoiseCov      [MAX_MEAS x MAX_MEAS] R
%   - numMeas           scalar   valid rows this epoch (0 = no measurement)
%   - coastMinusPrior   [n x 1]  x_coast - x_prior (increment space)
%   - coastCovariance   [n x n]  P_C of the coast
%
% OUTPUTS:
%   - passed     logical, 0 <= qValue < threshold(numMeas) AND numMeas >=
%                CST_spfParam.REVAL_MIN_MEAS; false if numMeas = 0 or the
%                residual covariance is unusable
%   - qValue     the chi-square statistic (log it for diagnostics); 0 when
%                it could not be formed
%   - solveFault logical, true if matrixInv flagged H*P_C*H' + R as
%                unusable: test invalid, no pass
%
% ASSUMPTIONS AND LIMITATIONS:
% Threshold from CST_spfParam.REVAL_THRESHOLD_TABLE(numMeas). The rows may
% include non-GNSS measurements (pressure altitude); REVAL_MIN_MEAS keeps
% an epoch with too few rows from counting as a pass. Everything runs on
% the full MAX_MEAS layout (rows/cols beyond numMeas are zero): no
% variable-size expression. The residual covariance is inverted by the
% host's matrixInv (SVD pseudo-inverse, unusable matrix flagged); a
% negative or non-finite q cannot pass.
%
% REQUIREMENT TRACEABILITY:
%
%******************************************************************************************
%#codegen
function [passed, qValue, solveFault] = SPF_revalidation(innovation, obsMatrix, ...
    measNoiseCov, numMeas, coastMinusPrior, coastCovariance)

% Define variables
passed     = false;
qValue     = 0.0;
solveFault = false;
numMeas    = min(uint8(numMeas), CST_spfParam.MAX_MEAS);   % defensive: table bound

if (numMeas > 0)
    threshold = CST_spfParam.REVAL_THRESHOLD_TABLE(numMeas);

    % full fixed-size layout; rows/cols beyond numMeas are zero padding
    residual    = innovation - obsMatrix * coastMinusPrior;
    residualCov = obsMatrix * coastCovariance * obsMatrix' + measNoiseCov;
    residualCov = (residualCov + residualCov') / 2;

    % Exception handler: a residual covariance with a non-positive
    % variance on a live row is unusable; the host's SVD inverse then
    % handles the padding (pseudo-inverse) and flags a non-finite matrix.
    varianceOk = true;
    for rowIdx = 1:double(numMeas)
        if ~(residualCov(rowIdx, rowIdx) > 0.0)
            varianceOk = false;
        end % ELSE is trivial
    end
    [sInverse, invInvalid] = matrixInv(residualCov);
    if (varianceOk) && (~invInvalid)
        qValue = residual' * (sInverse * residual);      % r' S_r^-1 r
        if isfinite(qValue) && (qValue >= 0.0) && (qValue < threshold) ...
                && (numMeas >= CST_spfParam.REVAL_MIN_MEAS)
            passed = true;
        end % ELSE: inconsistent, too few rows to certify, or an indefinite S_r (q < 0)
    else
        solveFault = true;                               % cannot certify anything this epoch
    end
end % ELSE: no measurement this epoch -> cannot validate

end
%------------------------------------------------------------------------
