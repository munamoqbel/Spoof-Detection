%******************************************************************************************
% DESCRIPTION:
% INS-only propagation of a state (or increment) and its covariance over one interval.
%
% INPUTS:  state [n x 1], covariance [n x n], propTel (.accumPhi / .accumQ)
% OUTPUTS: state, covariance propagated
%******************************************************************************************
%#codegen
function [state, covariance] = SPF_insCoast(state, covariance, propTel)

Phi = propTel.accumPhi;
Q   = propTel.accumQ;

state      = Phi * state;
covariance = Phi * covariance * Phi' + Q;
covariance = (covariance + covariance') / 2;

end
%------------------------------------------------------------------------
