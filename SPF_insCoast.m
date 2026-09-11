%******************************************************************************************
% DESCRIPTION:
% Propagate a state and covariance INS-only: no GNSS measurement update.
%
% INPUTS:
%   - state       [n x 1]   state (or increment vector) to propagate
%   - covariance  [n x n]   its covariance
%   - propTel     .accumPhi / .accumQ  interval Phi / Q from the 100 Hz side
%
% OUTPUTS:
%   - state, covariance propagated over the interval
%
% ASSUMPTIONS AND LIMITATIONS:
%
% REQUIREMENT TRACEABILITY:
% This is the "coasting" operation of the paper (E28 for the state,
% E31 for the covariance). Used by the gate for coast covariances and by
% the harness to extrapolate the operational KF while it is not updated.
%
%******************************************************************************************
%#codegen
function [state, covariance] = SPF_insCoast(state, covariance, propTel)

% Define variables
Phi = propTel.accumPhi;
Q   = propTel.accumQ;

% Output
state      = Phi * state;                                  % E28
covariance = Phi * covariance * Phi' + Q;                  % E31
covariance = (covariance + covariance') / 2;               % keep symmetric

end

%------------------------------------------------------------------------
