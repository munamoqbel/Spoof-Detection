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
%            =  innovation      - H (coastState - priorState)
% For the linear harness (innovation = z - H x_prior) this is exact.
%
% INPUTS:
%   - innovation        [MAX_MEAS x 1] host innovation z - h(x_prior)
%   - obsMatrix         [MAX_MEAS x n] H (host linearisation)
%   - measNoiseCov      [MAX_MEAS x MAX_MEAS] R / V
%   - numMeas           scalar   valid rows this epoch (0 = no GNSS)
%   - priorState        [n x 1]  host predicted state the innovation refers to
%   - coastState        [n x 1]  quarantined INS-only navigation state
%   - coastCovariance   [n x n]  its covariance
%
% OUTPUTS:
%   - passed     logical, qValue < threshold(numMeas); false if numMeas = 0
%   - qValue     the chi-square statistic (log it for diagnostics)
%
% ASSUMPTIONS AND LIMITATIONS:
% Threshold from CST_spfParam.REVAL_THRESHOLD_TABLE(numMeas).
%
% REQUIREMENT TRACEABILITY:
%
%******************************************************************************************
%#codegen
function [passed, qValue] = revalidation(innovation, obsMatrix, ...
    measNoiseCov, numMeas, priorState, coastState, coastCovariance)

% Define variables
passed = false;
qValue = 0.0;

if (numMeas > 0)
    threshold = CST_spfParam.REVAL_THRESHOLD_TABLE(numMeas);
    H = obsMatrix(1:numMeas, :);

    residual    = innovation(1:numMeas) - H * (coastState - priorState);
    residualCov = H * coastCovariance * H' + measNoiseCov(1:numMeas, 1:numMeas);
    residualCov = (residualCov + residualCov') / 2;

    qValue = residual' * (residualCov \ residual);
    if (qValue < threshold)
        passed  = true;
    end % ELSE is trivial
end % ELSE: no GNSS this epoch -> cannot validate

end
%------------------------------------------------------------------------
