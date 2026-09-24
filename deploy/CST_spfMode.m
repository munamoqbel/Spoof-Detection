%******************************************************************************************
% DESCRIPTION:
% FSM mode constants: NOMINAL (KF runs, monitors watch), COAST (alarm latched, INS-only
% output), PROBATION (trial filter on watch, output still coasting).
%******************************************************************************************
%#codegen
classdef CST_spfMode < uint8
    enumeration
        NOMINAL             (1)
        COAST               (2)
        PROBATION           (3)
    end
end
%------------------------------------------------------------------------------------------
