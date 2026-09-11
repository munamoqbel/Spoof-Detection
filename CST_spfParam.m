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
% - REVAL_MIN_MEAS:
% Minimum valid rows for a re-validation pass. The host's measurement
% vector always carries the pressure-altitude row, so numMeas >= 1 even in
% a GNSS outage; without this floor ten baro-only epochs in COAST would
% open a probation with no GNSS evidence. qReval is still computed and
% logged whenever numMeas > 0.
% - PIVOT_REL_TOL:
% Every S / residual covariance is inverted through its Cholesky factor.
% chol's flag catches indefinite and non-finite matrices; an EXACTLY
% singular one (duplicated row with zero noise) can still factorise with a
% rounding-level pivot, so pivots below PIVOT_REL_TOL * max(diag) are
% rejected too (condition number above ~1e12). The epoch is then dropped
% (CPI: xi = 0; re-validation: no pass) and info.solveFault is raised.
% - REVAL_THRESHOLD_TABLE(m) = chi2inv(1 - 1e-3, m), m = 1..MAX_MEAS (30)
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
        REVAL_THRESHOLD_TABLE = [ ...                 % chi2inv(1 - 1e-3, m), m = 1..30
            10.82756617066, 13.81551055796, 16.26623619624, 18.4668269529, ...
            20.51500565243, 22.45774448483, 24.32188634786, 26.12448155838, ...
            27.87716487126, 29.58829844507, 31.26413362024, 32.90949040736, ...
            34.52817897487, 36.1232736804, 37.69729821835, 39.25235479077, ...
            40.7902167069, 42.31239633168, 43.82019596452, 45.31474661813, ...
            46.79703804156, 48.26794229084, 49.72823246643, 51.17859777738, ...
            52.61965577617, 54.05196238858, 55.47602020575, 56.89228539335, ...
            58.30117348979, 59.70306430443];
        REVAL_DWELL_REQUIRED = uint8(10);         % passes -> probation
        REVAL_MIN_MEAS = uint8(4);                % rows needed for a re-validation PASS (3-D fix + clock;
                                                  % a pressure-only or single-satellite epoch cannot certify)
        PROBATION_LENGTH = uint8(18);             % quiet epochs -> commit
        MAX_ANCHOR_AGE = uint32(234);             % epochs, 0.5 s each
        MAX_MEAS = uint8(30);                     % max measurements per epoch (host arrays are sized to 30)
        PIVOT_REL_TOL = 1.0e-12;                  % Cholesky pivot floor (relative to max diag): below it a
                                                  % covariance is treated as singular (SPF_cpiMonitor, SPF_revalidation)


    end
end
%------------------------------------------------------------------------
