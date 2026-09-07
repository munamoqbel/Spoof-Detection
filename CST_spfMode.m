%******************************************************************************************
% DESCRIPTION:
% FSM mode constants. (Definition only.)
%
%   NavMode.NOMINAL   = 1   KF runs, monitors watch it
%   NavMode.COAST     = 2   alarm latched, output is INS-only
%   NavMode.PROBATION = 3   trial filter on watch, output still coasting
%
% ASSUMPTIONS AND LIMITATIONS:
%
% REQUIREMENT TRACEABILITY:
%
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
