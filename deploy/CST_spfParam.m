%***************************************************************************
% DESCRIPTION:
% Constants of the spoofing-detection gate: the window length and the
% three thresholds of the paper's design, the nominal bias allowance,
% the re-validation and probation rules of the recovery, the warm-up
% and the coast budget. The values are the published design; the host
% owns ARM_EPOCHS and NOMINAL_BIAS. The rationale of each constant is in
% the reference copy in the repository root.
%
% INPUTS:
%   none (constant properties)
%
% OUTPUTS:
%   none
%
% ASSUMPTIONS AND LIMITATIONS:
% - WINDOW_LENGTH, CPI_THRESHOLD, K_FALSE_ALERT and K_MISSED_DETECTION
%   belong together: a new window or a new allocation needs the offline
%   design (tools/designConstants, solve_N_min) for all four.
% - REVAL_THRESHOLD_TABLE must cover CST_gnssHybrid.MAX_MEASURES rows;
%   SPF_revalidation clamps the index to the table length.
% - NOMINAL_BIAS 0 is the design as published; the host sets per NED
%   axis the GNSS excursion it tolerates, and detection comes that much
%   later.
% - ARM_EPOCHS 10 is the value of the simulation references; the host
%   runs 120, the 60 s its biases take to converge after a reset.
%
% REQUIREMENT TRACEABILITY:
% Paper (Kujur 2024): N and T_N in Eq. 33 and 35; k_FA and k_MD in the
% SS gate and the protection level Eq. 51. Implementation Guide: Sec. 5,
% one chi-square gate per row count (REVAL_THRESHOLD_TABLE); rule 4,
% the coast budget on the protection level (ARM_PL_MAX and the time belt
% COAST_BUDGET_EPOCHS).
%***************************************************************************
%#codegen
classdef CST_spfParam
    properties (Constant = true)

        WINDOW_LENGTH = uint8(10);                        % N, epochs (Eq. 33)
        CPI_THRESHOLD = 45.639677645740022;               % T_N (Eq. 35)
        K_FALSE_ALERT = 5.233126417847868;                % k_FA, SS gate (Eq. 51)
        K_MISSED_DETECTION = 4.753424308817088;           % k_MD, PL term (Eq. 51)
        MONITORED_AXES = [1 2 3];                         % NED position state indices
        NOMINAL_BIAS = [0.0 0.0 0.0];                     % B per NED axis (m): SS gate + B, PL + B, CPI deadband B
        COAST_BUDGET_EPOCHS = uint32(1200);               % time belt behind the PL coast budget (rule 4)
        REVAL_THRESHOLD = 39.252354790768472;             % chi-square gate, m = 16 (reference only)
        REVAL_THRESHOLD_TABLE = [ ...                     % chi2inv(1 - 1e-3, m), m = 1..60 (Guide Sec. 5)
            10.82756617066, 13.81551055796, 16.26623619624, 18.46682695290, ...
            20.51500565243, 22.45774448483, 24.32188634786, 26.12448155838, ...
            27.87716487126, 29.58829844507, 31.26413362024, 32.90949040736, ...
            34.52817897487, 36.12327368040, 37.69729821835, 39.25235479077, ...
            40.79021670690, 42.31239633168, 43.82019596452, 45.31474661813, ...
            46.79703804156, 48.26794229084, 49.72823246643, 51.17859777738, ...
            52.61965577617, 54.05196238858, 55.47602020575, 56.89228539335, ...
            58.30117348979, 59.70306430443, 61.09830608106, 62.48721905709, ...
            63.87009852234, 65.24721746094, 66.61882884370, 67.98516762602, ...
            69.34645249624, 70.70288741151, 72.05466295199, 73.40195751899, ...
            74.74493839842, 76.08376270770, 77.41857824131, 78.74952422804, ...
            80.07673201082, 81.40032565871, 82.72042251912, 84.03713371722, ...
            85.35056460859, 86.66081519040, 87.96798047563, 89.27215083430, ...
            90.57341230530, 91.87184688166, 93.16753277223, 94.46054464188, ...
            95.75095383249, 97.03882856651, 98.32423413474, 99.60723306985];
        REVAL_DWELL_REQUIRED = uint8(10);                 % consecutive passes -> PROBATION
        REVAL_MIN_MEAS = uint8(4);                        % rows for a pass: a 3-D fix and the clock
        ARM_EPOCHS = uint8(10);                           % qualifying epochs before the monitors arm (host: 120)
        ARM_PL_MAX = 100.0;                               % m: qualifying when K_MISSED_DETECTION * sigma_pos < it
        PROBATION_LENGTH = uint8(18);                     % quiet epochs -> hand-back
        COAST_SIGMA_MAX = 3.0;                            % m: the coast sigma the re-admission tests may use at most, per bounded state
        COAST_BOUND_STATES = [1 2 3 16];                  % the states bounded: the NED position and the receiver clock bias
        RETURN_BANK = uint8(3);                           % long coasts kept beside the anchor, restarted in turn
        RETURN_RESTART = uint32(80);                      % epochs between restarts: the oldest spans 160 to 240 epochs
        RETURN_WANDER = 15.0;                             % m: the nav must have wandered this far from the long coast,
        RETURN_RADIUS = 15.0;                             % m: and the GNSS must land this close to it, for a return
        MAX_ANCHOR_AGE = uint32(234);                     % epochs: an older anchor is refused at a latch
        COAST_Q_DATASHEET = false;                        % the gate's coasts on the IMU datasheet Q (SPF_coastQ); false: the host's accumQ
        IMU_ARW = 0.05;                                   % deg/sqrt(h)
        IMU_VRW = 0.0147;                                 % m/s/sqrt(h)
        IMU_GYRO_BIAS_SIGMA = 0.03;                       % deg/h, in-run
        IMU_ACC_BIAS_SIGMA = 4.0;                         % micro-g, in-run
        IMU_BIAS_TAU = 300.0;                             % s
        IMU_GYRO_WALK = 0.016;                            % deg/h per sqrt(h)
        IMU_ACC_WALK = 10.0;                              % micro-g per sqrt(h)
        COAST_POS_STATES = [1 2 3];                       % the blocks in the host's state layout
        COAST_VEL_STATES = [4 5 6];
        COAST_ATT_STATES = [7 8 9];
        COAST_ACC_BIAS_STATES = [10 11 12];
        COAST_GYRO_BIAS_STATES = [13 14 15];
        COAST_Q_INTERVAL = 0.5;                           % s, when accumPhi carries no position-from-velocity entry
        COAST_SIGMA_FLOOR_STATES = zeros(1, 0);           % the host's P floored on these states
        COAST_SIGMA_FLOOR_VALUES = zeros(1, 0);           % at these sigmas
        MAX_MEAS = CST_gnssHybrid.MAX_MEASURES;           % rows per epoch, the host's bound

    end
end
%------------------------------------------------------------------------
