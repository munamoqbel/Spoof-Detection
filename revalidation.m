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
% (S_r = Rc'*Rc, q = |Rc'\residual|^2 >= 0): no explicit inverse, no
% division by a pivot that can be zero, and a non-PD S_r is detected by
% chol's second output (plus a relative pivot floor, PIVOT_REL_TOL, for an
% S_r singular up to rounding) instead of producing a negative or NaN q.
%
% REQUIREMENT TRACEABILITY:
%
%******************************************************************************************
%#codegen
function [passed, qValue, solveFault] = revalidation(innovation, obsMatrix, ...
    measNoiseCov, numMeas, coastMinusPrior, coastCovariance)

% Define variables
passed     = false;
qValue     = 0.0;
solveFault = false;
numMeas    = min(uint8(numMeas), CST_spfParam.MAX_MEAS);   % defensive: table / buffer bound

if (numMeas > 0)
    threshold = CST_spfParam.REVAL_THRESHOLD_TABLE(numMeas);
    H = obsMatrix(1:numMeas, :);

    residual    = innovation(1:numMeas) - H * coastMinusPrior;
    residualCov = H * coastCovariance * H' + measNoiseCov(1:numMeas, 1:numMeas);
    residualCov = (residualCov + residualCov') / 2;

    % Exception handler: residual covariance must be positive definite
    [cholFactor, cholFail] = chol(residualCov);
    pivotOk = (cholFail == 0);
    if (pivotOk)
        minPivot = min(diag(cholFactor)) ^ 2;          % rounding-level pivot = singular
        pivotOk  = (minPivot > CST_spfParam.PIVOT_REL_TOL * max(diag(residualCov)));
    end % ELSE is trivial
    if (pivotOk)
        whitened = cholFactor' \ residual;               % Rc'^-1 r
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
