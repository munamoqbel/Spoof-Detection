%******************************************************************************************
% DESCRIPTION:
% Stand-in for the simulation-side CST_gnssHybrid constant class. The
% target simulation defines its own copy; this file exists only so the
% standalone harness (run_recovery.m / run_scheduler_test.m) resolves
% CST_gnssHybrid.NO_STATES without the simulation on the path.
%
% Per the user, the simulation class carries only NO_STATES = 60.
% Keep this file OFF the path when running inside the simulation.
%
%******************************************************************************************
%#codegen
classdef CST_gnssHybrid
    properties (Constant = true)
        NO_STATES = 60;     % total hybrid KF states
    end
end
