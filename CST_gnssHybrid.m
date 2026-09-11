%******************************************************************************************
% DESCRIPTION:
% Stand-in for the simulation-side CST_gnssHybrid constant class. The
% target simulation defines its own copy; this file exists only so the
% standalone harness (run_recovery.m / run_scheduler_test.m) resolves
% CST_gnssHybrid.NO_STATES without the simulation on the path.
%
% The simulation class carries NO_STATES = 60 and MAX_MEASURES = uint8(51);
% only those two are used here. Keep this file OFF the path when running
% inside the simulation.
%
%******************************************************************************************
%#codegen
classdef CST_gnssHybrid
    properties (Constant = true)
        NO_STATES    = 60;          % total hybrid KF states
        MAX_MEASURES = uint8(51);   % max measurement rows per epoch (host arrays are sized to it)
    end
end
