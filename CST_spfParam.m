%******************************************************************************************
% DESCRIPTION:
% Define all global constants specifically used for WSM1 project.
%
% ASSUMPTIONS AND LIMITATIONS:
% - MAX_ANCHOR_AGE:
% Max anchor age at latch; drift budget @ 10 m alert limit (harness Phi/Q,
% 02/09/2026); older anchor -> fallback refused + anchorMissing.
% NOTE: the anchor is refreshed on every alarm-free clean close, so at a
% latch it is normally 1 epoch old; this guard only bites for the startup
% anchor (initial state) before the first clean close. The real coasting
% budget (max coast time, Implementation Guide rule 4) is not yet enforced.
% - MAX_MEAS:
% Upper bound on GNSS measurements per epoch. All monitor buffers are sized
% to it; the per-epoch valid count (numMeas) selects the live rows.
% - REVAL_THRESHOLD_TABLE(m) = chi2inv(1 - 1e-3, m), m = 1..MAX_MEAS
% (Implementation Guide Sec. 5: one gate per possible satellite count).
% Recompute offline if MAX_MEAS or the 1e-3 allocation changes.
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
        MONITORED_AXES = [1 2 3];                         % NED position state indices
        REVAL_THRESHOLD = 39.252354790768472;     % chi-square gate, m = 16 (kept for reference)
        REVAL_THRESHOLD_TABLE = [ ...                 % chi2inv(1 - 1e-3, m), m = 1..16
            10.8275661706627, 13.8155105579643, 16.2662361962381, 18.4668269529032, ...
            20.5150056524329, 22.4577444848253, 24.3218863478569, 26.1244815583761, ...
            27.8771648712566, 29.5882984450744, 31.2641336202400, 32.9094904073602, ...
            34.5281789748709, 36.1232736803981, 37.6972982183538, 39.2523547907685];
        REVAL_DWELL_REQUIRED = uint8(10);         % passes -> probation
        PROBATION_LENGTH = uint8(18);             % quiet epochs -> commit
        MAX_ANCHOR_AGE = uint32(234);             % epochs, 0.5 s each
        MAX_MEAS = uint8(16);                     % max GNSS measurements per epoch (8 SV x code/carrier)


    end
end
%------------------------------------------------------------------------
