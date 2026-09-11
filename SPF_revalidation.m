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
%   - passed     logical, qValue < threshold(numMeas) AND numMeas >=
%                CST_spfParam.REVAL_MIN_MEAS; false if numMeas = 0 or the
%                residual covariance is not positive definite
%   - qValue     the chi-square statistic (log it for diagnostics); 0 when
%                it could not be formed
%   - solveFault logical, true if H*P_C*H' + R was not positive definite
%                (singular / indefinite / non-finite): test invalid, no pass
%
% ASSUMPTIONS AND LIMITATIONS:
% Threshold from CST_spfParam.REVAL_THRESHOLD_TABLE(numMeas). The rows may
% include non-GNSS measurements (pressure altitude); REVAL_MIN_MEAS keeps
% an epoch with too few rows from counting as a pass. The quadratic form
% is evaluated through the Cholesky factor of the residual covariance
% (S_r = L*L', q = |L\residual|^2 >= 0) with the fixed-size loops
% SPF_cholesky / SPF_forwardSubst: no explicit inverse, no library call,
% no division by a pivot that can be zero, and a non-PD S_r is flagged
% (plus a relative pivot floor, PIVOT_REL_TOL, for an S_r singular up to
% rounding) instead of producing a negative or NaN q.
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
maxMeas    = CST_spfParam.MAX_MEAS;
numMeas    = min(uint8(numMeas), maxMeas);   % defensive: table / buffer bound

if (numMeas > 0)
    threshold = CST_spfParam.REVAL_THRESHOLD_TABLE(numMeas);

    % Fixed-size formulation (no variable-size expressions, see
    % SPF_cpiMonitor): rows beyond numMeas of H, R and the innovation are
    % zeroed, so residualCov = H P_C H' + R is zero there, and its padding
    % diagonal is set to a positive value -> blkdiag(S_r, padValue*I).
    H        = obsMatrix;                                  % [MAX_MEAS x n]
    residual = innovation - H * coastMinusPrior;           % [MAX_MEAS x 1]
    R        = measNoiseCov;                               % [MAX_MEAS x MAX_MEAS]
    for rowIdx = (double(numMeas) + 1):double(maxMeas)
        H(rowIdx, :)  = 0.0;
        residual(rowIdx) = 0.0;
        R(rowIdx, :)  = 0.0;
        R(:, rowIdx)  = 0.0;
    end
    residualCov = H * coastCovariance * H' + R;
    residualCov = (residualCov + residualCov') / 2;
    padValue = max(max(diag(residualCov)), 1.0);
    for rowIdx = (double(numMeas) + 1):double(maxMeas)
        residualCov(rowIdx, rowIdx) = padValue;
    end

    % Exception handler: residual covariance must be positive definite
    [cholLower, pivotOk] = SPF_cholesky(residualCov, numMeas);   % S_r = L * L' (live block)
    if (pivotOk)
        minPivot = min(diag(cholLower)) ^ 2;           % rounding-level pivot = singular
        pivotOk  = (minPivot > CST_spfParam.PIVOT_REL_TOL * max(diag(residualCov)));
    end % ELSE is trivial
    if (pivotOk)
        whitened = SPF_forwardSubst(cholLower, residual, numMeas);   % L^-1 r  (zero in the padding)
        qValue   = whitened' * whitened;                 % r' S_r^-1 r  (>= 0)
        if (qValue < threshold) && (numMeas >= CST_spfParam.REVAL_MIN_MEAS)
            passed = true;
        end % ELSE: inconsistent, or too few rows to certify (e.g. pressure only)
    else
        solveFault = true;                               % cannot certify anything this epoch
    end
end % ELSE: no measurement this epoch -> cannot validate

end
%------------------------------------------------------------------------
