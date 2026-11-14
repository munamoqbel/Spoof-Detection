%******************************************************************************************
% DESCRIPTION:
% The host's P as the gate reads it: floored on the states the host's own
% consistency reads say it underestimates (CST_spfParam.COAST_SIGMA_FLOOR_*).
% A host whose P is honest leaves the lists empty and the P passes unchanged.
% Applied once per epoch to kfMeas.postCov, so that every coast the gate
% starts (window, anchor, long coast) and the filter covariance the
% separation test subtracts carry the same floor: the floor itself cancels in
% the separation, its growth through the dynamics does not.
%
% INPUTS:
%   postCov  [n x n]  the host's P+ this epoch
%
% OUTPUTS:
%   covOut   [n x n]  the same with the floors
%******************************************************************************************
%#codegen
function [covOut] = SPF_coastStartCov(postCov)

covOut = SPF_floorCov(postCov, CST_spfParam.COAST_SIGMA_FLOOR_STATES, CST_spfParam.COAST_SIGMA_FLOOR_VALUES);
end
