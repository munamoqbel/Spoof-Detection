%***************************************************************************
% DESCRIPTION:
% Modes of the spoofing-protected navigator (SPF_protectedNav), as uint8
% so that the mode travels in the telemetry and compiles in a switch.
%   NOMINAL    host KF updated normally, the monitors watch its increments
%   COAST      alarm latched: INS-only output, raw GNSS re-validated
%   PROBATION  trial filter watched, output still on the coast
%
% INPUTS:
%   none (enumeration)
%
% OUTPUTS:
%   none
%
% ASSUMPTIONS AND LIMITATIONS:
% The numeric values are part of the telemetry format (spfLog.mode):
% 1, 2 and 3 are not to be renumbered.
%
% REQUIREMENT TRACEABILITY:
% Paper Sec. 5 recovery: detection (NOMINAL to COAST), re-validation and
% probation (COAST to PROBATION), hand-back (PROBATION to NOMINAL).
%***************************************************************************
%#codegen
classdef CST_spfMode < uint8
    enumeration
        NOMINAL             (1)
        COAST               (2)
        PROBATION           (3)
    end
end
%------------------------------------------------------------------------------------------
