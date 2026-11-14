%******************************************************************************************
% DESCRIPTION:
% A covariance with sigma floors on named states: every listed diagonal below
% its floor is raised to it, the rest untouched. Adding to the diagonal keeps
% the matrix positive semi-definite.
%
% INPUTS:
%   covIn        [n x n]
%   floorStates  [1 x m]  state indices
%   floorSigmas  [1 x m]  the floors, in the states' units (sigma, not variance)
%
% OUTPUTS:
%   covOut       [n x n]
%******************************************************************************************
%#codegen
function [covOut] = SPF_floorCov(covIn, floorStates, floorSigmas)

covOut = covIn;
for k = 1:numel(floorStates)
    idx = floorStates(k);
    floorVar = floorSigmas(k) * floorSigmas(k);
    if covOut(idx, idx) < floorVar
        covOut(idx, idx) = floorVar;
    end
end
end
