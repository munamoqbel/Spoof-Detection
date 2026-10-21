%***************************************************************************
% DESCRIPTION:
% INS-only propagation of a state (or an increment) and its covariance
% over one epoch interval, with the interval matrices the host
% accumulated on the 100 Hz side: the coast of the paper.
%
% INPUTS:
%   - state       [n x 1] state or increment to propagate
%   - covariance  [n x n] its covariance
%   - propTel     .accumPhi [n x n], .accumQ [n x n] interval Phi and Q
%
% OUTPUTS:
%   - state, covariance  propagated over the interval
%
% ASSUMPTIONS AND LIMITATIONS:
% The covariance is symmetrised after the propagation. No measurement
% enters here.
%
% REQUIREMENT TRACEABILITY:
% Paper E28 (coast state) and E31 (coast covariance).
%***************************************************************************
%#codegen
function [state, covariance] = SPF_insCoast(state, covariance, propTel)

Phi = propTel.accumPhi;
Q   = propTel.accumQ;

state      = Phi * state;                                 % E28
covariance = Phi * covariance * Phi' + Q;                 % E31
covariance = (covariance + covariance') / 2;

end
%------------------------------------------------------------------------
