function prm = kujur_params(N_min_in, CPI_threshold_in, k_FA_in, k_MD_in)
%KUJUR_PARAMS  Parameter struct for the Kujur et al. (NAVIGATION 2024) dual monitor.
%
% CALL PATTERNS
% -------------
%   prm = kujur_params()
%       Bootstrap: structural fields only. Monitor constants are
%       placeholders (0). Use to get n_states/idx_z before design.
%
%   prm = kujur_params(N_min, CPI_threshold, k_FA, k_MD)
%       Full: monitor constants from solve_N_min() embedded.
%       This is the struct passed to dual_monitor().
%
% MATLAB CODER
% ------------
%   All fields are scalar doubles (N_min used as loop bound must be a
%   compile-time constant at codegen: wrap prm with coder.Constant in
%   the codegen script, see codegen_dual_monitor.m).
%   No varargin (unsupported); explicit named optional args + nargin.
%
% UNITS: SI (m, s, rad).

% ----------------------------------------------------------------------
%  A. Structural dimensions - EDIT to match your EKF
%     State layout (paper Appendix A, Eq. A1/A4/A12/A13):
%       1:3  pos | 4:6 vel | 7:9 att | 10:12 acc bias | 13:15 gyro bias
%       16:17 rcvr clock | 18:end per-SV GNSS error states
% ----------------------------------------------------------------------
prm.n_states = CST_gnssHybrid.NO_STATES;     % 60. total KF states
prm.m_meas   = 16;     % measurements/epoch = 2*n_sv (code+carrier)
prm.n_sv     = 8;
prm.idx_z    = 3;      % vertical position index (monitor direction u = e_idx_z)

% ----------------------------------------------------------------------
%  B. Timing
% ----------------------------------------------------------------------
prm.fs = 2.0;                          % GNSS rate (Hz), paper uses 2 Hz
prm.dt = 1.0 / prm.fs;

% ----------------------------------------------------------------------
%  C. Integrity budget (paper Section 5)
% ----------------------------------------------------------------------
prm.P_FA_plus  = 0.5e-5;        % CPI false-alert allocation (P_FA+)
prm.P_FA_minus = 0.5e-5;        % SS  false-alert allocation (P_FA-, total)
prm.P_MD_plus  = 1e-7;          % CPI missed-detection req.  (P_MD+)
prm.P_MD_minus = 1e-6;          % SS  missed-detection req.  (P_MD-)

% ----------------------------------------------------------------------
%  D. IMU process noise (Appendix F, Table F1, navigation grade)
% ----------------------------------------------------------------------
prm.VRW_sq       = (0.18 / 60.0)^2; %(1.43e-2 / 60.0)^2;                        % (m/s)^2/s
prm.ARW_sq       = (0.2 * pi/180.0 / 60.0)^2; %(1e-3 * pi/180.0 / 60.0)^2;      % rad^2/s
prm.acc_bias_sq  = (4e-2   * 9.81e-3)^2; %(1e-2   * 9.81e-3)^2;                 % (m/s^2)^2
prm.gyro_bias_sq = (7 * pi/180.0 / 3600.0)^2; %(3.5e-3 * pi/180.0 / 3600.0)^2;  % (rad/s)^2

% ----------------------------------------------------------------------
%  E. Solved monitor constants (from solve_N_min, offline)
%     gaminv/norminv are NOT Coder-supported: values are computed
%     offline and passed in here as plain doubles.
% ----------------------------------------------------------------------
if nargin >= 4
    prm.N_min         = N_min_in;
    prm.CPI_threshold = CPI_threshold_in;    % T_N  (Eq. 35)
    prm.k_FA          = k_FA_in;             % SS threshold multiplier
    prm.k_MD          = k_MD_in;             % SS PL multiplier (Eq. 51)
    prm.P_FA_ss_win   = prm.P_FA_minus / N_min_in;
else
    prm.N_min         = 20;                  % placeholder
    prm.CPI_threshold = 0.0;                 % placeholder - NOT valid
    prm.k_FA          = 0.0;
    prm.k_MD          = 0.0;
    prm.P_FA_ss_win   = prm.P_FA_minus / 20.0;
end

end
