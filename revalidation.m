%******************************************************************************************
% DESCRIPTION:
% Ckech if GNSS statistically consistent with the clean coast.
%
% Chi-square consistency test of the raw measurements against the
% INS-only coast. This is the end-of-attack detector used during COAST:
% while spoofing is active, the offset (vs mm-level carrier noise) makes
% q_value astronomically large; after the spoofer stops, q_value drops
% to ~chi2(m) and passes. NOTE: one pass means nothing on its own - the
% caller requires a dwell of consecutive passes, and even then only
% opens PROBATION, never a direct handback.
%
% INPUTS:
%   - measurement       [m x 1]  raw GNSS measurement z (untrusted)
%   - obs_matrix        [m x n]  H (linearised at the coast in a real EKF)
%   - meas_noise_cov    [m x m]  V
%   - coast_state       [n x 1]  quarantined INS-only navigation state
%   - coast_covariance  [n x n]  its covariance
%   - threshold         scalar   precomputed gaminv(1-alpha, m/2, 2)
%
% OUTPUTS:
%   - passed     logical, q_value < threshold
%   - q_value    the chi-square statistic (log it for diagnostics)
%
% ASSUMPTIONS AND LIMITATIONS:
%
% REQUIREMENT TRACEABILITY:
%
%******************************************************************************************
%#codegen
function [passed, qValue] = revalidation(measurement, obsMatrix, ...
    measNoiseCov, coastState, coastCovariance)

% Define variables
threshold = CST_spfParam.REVAL_THRESHOLD;
passed = false;

residual    = measurement - obsMatrix * coastState;
residualCov = obsMatrix * coastCovariance * obsMatrix' + measNoiseCov;
residualCov = (residualCov + residualCov') / 2;

qValue = residual' * (residualCov \ residual);
if (qValue < threshold)
    passed  = true;
end % ELSE is trivial

end
%------------------------------------------------------------------------------------------
