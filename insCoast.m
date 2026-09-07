%******************************************************************************************
% DESCRIPTION:
% Propagate a state INS-only: no GNSS measurement update.
%
% INPUTS:
%   - state       [n x 1]   state to propagate
%   - covariance  [n x n]   its covariance
%   - Phi, Q      [n x n]   state transition and process noise
%
% OUTPUTS:
%   - pool (updated)
%   - report (PoolReport - see that file)
%
% ASSUMPTIONS AND LIMITATIONS:
%
% REQUIREMENT TRACEABILITY:
% This is the "coasting" operation of the paper (E28 for the state,
% under E31 for the covariance). Used for:
%   - the navigation output while GNSS is severed (COAST / PROBATION)
%   - bringing the quarantine anchor forward to the current epoch
%
%******************************************************************************************
%#codegen
function [state, covariance] = insCoast(state, covariance, spoofInfo)

% Define variables
Phi = spoofInfo.phiAcc;
Q = spoofInfo.qAcc;

% Output
state      = Phi * state;                                  % E28
covariance = Phi * covariance * Phi' + Q;                  % covariane matrix. under E31
covariance = (covariance + covariance') / 2;               % keep symmetric to prevent negative diagonal terms.

end

%------------------------------------------------------------------------------------------
