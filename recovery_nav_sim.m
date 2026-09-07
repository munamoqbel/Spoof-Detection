function out = recovery_nav_sim(z_all, H_all, V, spoofInfo, P0, prm, use_fed)

if nargin<7
    use_fed = false;
end
%RECOVERY_NAV_SIM  Batch wrapper: run the protected navigator over a
%                  logged measurement set. THIN LOOP ONLY.
%
% State is managed EXPLICITLY here (no persistent) so batch runs are
% repeatable. The deployment wrapper with persistent state is
% spoof_monitor_2hz.m.
%
% CALL GRAPH
%   recovery_nav_sim
%     KujurConfig.fromPrmRec              (struct class - config)
%     NavState.create                     (struct class - sys; uses PoolState)
%     protected_nav_step                  (FSM: NOMINAL / COAST / PROBATION)
%       kalman_update_step
%       monitor_pool_step                 (windows; report = PoolReport.blank)
%         ss_monitor_step                 (SS per axis; result = SsResult.blank)
%         cpi_monitor_window              (CPI at close)
%       monitor_pool_open_window
%       monitor_pool_clear
%       ins_coast_step
%       revalidation_test
%
% Interface identical to previous versions: run_recovery.m unchanged.
%
% INPUTS
%   z_all [m x N], H_all [m x n x N], V [m x m], Phi/Q [n x n]
%   P0    [n x n]   Riccati-converged initial covariance
%   prm             kujur_params(N_min, T_N, k_FA, k_MD)
%                   optional prm.mon_axes (default: prm.idx_z only)
%   rec             .T_reval, .M_dwell, optional .M_prob
%
% OUTPUT struct 'out':
%   x_nav [n x N], sig_pos [3 x N], sig_z [1 x N], state [1 x N],
%   CPI_alarm/SS_alarm [N x 1], alarm_axis [N x 3], SS_PL [N x 1],
%   q_reval [N x 1], dwell [N x 1],
%   ev.t_detect / t_anchor / t_prob_start / t_prob_fail / t_handback

% ----------------------------------------------------------------------
%  configuration and initial state (struct classes)
% ----------------------------------------------------------------------
mode = CST_spfMode.NOMINAL;
filter = STRUCT_SPF.setFilter(zeros(prm.n_states, 1), P0);
trial = STRUCT_SPF.setTrial(zeros(prm.n_states, 1), P0);
pool = STRUCT_SPF.zeroMonitorPool(prm.m_meas);
anchor = STRUCT_SPF.setAnchor(false, zeros(prm.n_states, 1), P0, 0);
sys = STRUCT_SPF.setSys(mode, filter, trial, pool, anchor, 0, 0);

% ----------------------------------------------------------------------
%  allocate outputs
% ----------------------------------------------------------------------
num_epochs = size(z_all, 2);

out.x_nav          = zeros(prm.n_states, num_epochs);
out.sig_pos        = zeros(3, num_epochs);
out.sig_z          = zeros(1, num_epochs);
out.state          = zeros(1, num_epochs);
out.CPI_alarm      = false(num_epochs, 1);
out.SS_alarm       = false(num_epochs, 1);
out.alarm_axis     = false(num_epochs, 3);
out.SS_PL          = zeros(num_epochs, 1);
out.q_reval        = zeros(num_epochs, 1);
out.reval_computed = false(num_epochs, 1);
out.dwell          = zeros(num_epochs, 1);
out.anchor_missing  = false(num_epochs, 1);

event_detect     = [];
event_anchor     = [];
event_prob_start = [];
event_prob_fail  = [];
event_handback   = [];

kf_x = zeros(60, 1);
kf_P = P0;

% ----------------------------------------------------------------------
%  run
% ----------------------------------------------------------------------
for epoch = 1:num_epochs

    if (use_fed)
        % ---- your 2a/2b: the KF updates ALWAYS (free-running) ----
        [kf_x, kf_P, innov, innov_S, ~] = kalman_update_step(kf_x, kf_P, ...
            z_all(:, epoch), H_all(:, :, epoch), spoofInfo, V);

        % ---- your 2c: the fed FSM ----
        [sys, spoofTel] = protectedNav(sys, innov, innov_S, ...
            H_all(:, :, epoch), kf_x, kf_P, z_all(:, epoch), ...
            V, spoofInfo, epoch);

        % ---- your 2d: reseed on command ----
        if spoofTel.kfCommand.reseedKF
            kf_x = spoofTel.kfCommand.reseedState;
            kf_P = spoofTel.kfCommand.reseedCov;
        end
        out.reseed_log(epoch) = spoofTel.kfCommand.reseedKF;

    else
        [sys, spoofTel] = protected_nav_step(sys, z_all(:, epoch), ...
            H_all(:, :, epoch), V, spoofInfo, epoch);
    end

    % ---- per-epoch logs ----
    out.x_nav(:, epoch)         = spoofTel.nav.state;
    out.sig_pos(:, epoch)       = spoofTel.nav.sigmaPosition;
    out.sig_z(epoch)            = spoofTel.nav.sigmaPosition(prm.idx_z);
    out.state(epoch)            = spoofTel.info.mode;
    out.SS_alarm(epoch)         = spoofTel.info.ssAlarm;
    out.CPI_alarm(epoch)        = spoofTel.info.cpiAlarm;
    out.alarm_axis(epoch, :)    = spoofTel.info.alarmPerAxis;
    out.SS_PL(epoch)            = spoofTel.info.maxProtectionLevel;
    out.q_reval(epoch)          = spoofTel.info.qReval;
    out.reval_computed(epoch)   = spoofTel.info.revalComputed;
    out.dwell(epoch)            = spoofTel.info.dwellCount;
    out.anchor_missing(epoch)   = spoofTel.info.anchorMissing;

    % ---- event logs ----
    if spoofTel.info.eventLatched
        event_detect(end+1) = epoch;                                    %#ok<AGROW>
        event_anchor(end+1) = spoofTel.info.eventAnchorEpoch;           %#ok<AGROW>
    end
    if spoofTel.info.eventProbationStarted
        event_prob_start(end+1) = epoch;                                %#ok<AGROW>
    end
    if spoofTel.info.eventProbationVetoed
        event_prob_fail(end+1) = epoch;                                 %#ok<AGROW>
    end
    if spoofTel.info.eventHandback
        event_handback(end+1) = epoch;                                  %#ok<AGROW>
    end
end

out.ev.t_detect     = event_detect;
out.ev.t_anchor     = event_anchor;
out.ev.t_prob_start = event_prob_start;
out.ev.t_prob_fail  = event_prob_fail;
out.ev.t_handback   = event_handback;
end
