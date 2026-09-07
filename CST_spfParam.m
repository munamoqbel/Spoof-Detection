%******************************************************************************************
% DESCRIPTION:
% Define all global constants specifically used for WSM1 project.
%
% ASSUMPTIONS AND LIMITATIONS:
% - MAX_ANCHOR_AGE:
% Max anchor age at latch; drift budget @ 10 m alert limit (harness Phi/Q,
% 02/09/2026); older anchor -> fallback refused + anchorMissing
%
% REQUIREMENT TRACEABILITY:
%
%******************************************************************************************
%#codegen
classdef CST_spfParam
    properties (Constant = true)

        %------------------------------------------------------------------
        % Kujur Configuration
        %------------------------------------------------------------------
        WINDOW_LENGTH = uint8(10);                        % N (epochs)
        CPI_THRESHOLD = 45.639677645740022;               % T_N  (Eq. 35)
        K_FALSE_ALERT = 5.233126417847868;                % SS gate
        K_MISSED_DETECTION = 4.753424308817088;           % PL term
        MONITORED_AXES = [1 2 3];
        REVAL_THRESHOLD = 39.252354790768472;     % chi-square gate
        REVAL_DWELL_REQUIRED = uint8(10);         % passes -> probation
        PROBATION_LENGTH = uint8(18);             % quiet epochs -> commit
        MAX_ANCHOR_AGE = uint32(234);             % epochs, 0.5 s each


    end
end
%------------------------------------------------------------------------
