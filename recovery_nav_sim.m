function out = recovery_nav_sim(z_all, H_all, V, propTel, P0, prm)
%RECOVERY_NAV_SIM  Batch harness: a HOST that follows the gate contract
%                  (docs/HOST_2HZ_WIRING.m) over a logged measurement set.
%
% The host owns its filters:
%   - operational KF  kf_x / kf_P : updated in NOMINAL; only extrapolated in
%                     COAST / PROBATION (= the INS-only coast)
%   - trial KF        tr_x / tr_P : copy of the coasting KF when probation
%                     opens; updated in PROBATION
% and applies the gate's corrections at LATCH and COMMIT. The gate (the
% FSM in SPF_protectedNav) is called with the ACTIVE filter's update result.
% State is kept explicitly (no persistent) so batch runs are repeatable;
% the deployment wrapper with persistent state is SPF_spoofMonitor.m.
%
% INPUTS
%   z_all [m x N], H_all [m x n x N], V [m x m]
%   propTel         .accumPhi / .accumQ [n x n] interval matrices (constant here)
%   P0    [n x n]   Riccati-converged initial covariance
%   prm             kujur_params(...) (only n_states / idx_z are used)
%
% OUTPUT struct 'out':
%   x_nav [n x N] (host operational KF = protected solution), sig_pos [3 x N],
%   sig_z [1 x N], state [1 x N], CPI_alarm/SS_alarm [N x 1],
%   alarm_axis [N x 3], SS_PL [N x 1], q_reval [N x 1], reval_computed,
%   dwell, anchor_missing, coast_epochs, correction_log [N x 1],
%   ev.t_detect / t_anchor / t_prob_start / t_prob_fail / t_handback

n = prm.n_states;
num_epochs = size(z_all, 2);
num_meas   = size(z_all, 1);                 % constant in the harness

% ---- gate state (explicit) ----
sys = STRUCT_SPF.zeroSys;
sys.coastCov = P0;
sys.anchor   = STRUCT_SPF.setAnchor(true, zeros(n, 1), P0, uint32(0));   % startup anchor

% ---- host filters ----
kf_x = zeros(n, 1);  kf_P = P0;
tr_x = zeros(n, 1);  tr_P = P0;
mode = CST_spfMode.NOMINAL;

% ---- logs ----
out.x_nav          = zeros(n, num_epochs);
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
out.anchor_missing = false(num_epochs, 1);
out.coast_epochs   = zeros(num_epochs, 1);
out.correction_log = false(num_epochs, 1);

event_detect = []; event_anchor = []; event_prob_start = [];
event_prob_fail = []; event_handback = [];

for epoch = 1:num_epochs
    inProbation = (mode == CST_spfMode.PROBATION);
    kfUpdated   = (mode == CST_spfMode.NOMINAL);

    % ---- host: update the ACTIVE filter (your kfUpdate) ----
    if inProbation
        [tr_x, tr_P, y, ~, ~, prior, priorP] = kalman_update_step(tr_x, tr_P, ...
            z_all(:, epoch), H_all(:, :, epoch), propTel, V);
        post = tr_x; postP = tr_P;
    else
        [ax, aP, y, ~, ~, prior, priorP] = kalman_update_step(kf_x, kf_P, ...
            z_all(:, epoch), H_all(:, :, epoch), propTel, V);
        post = ax; postP = aP;
    end
    % as the host: kfUpdate exposes y, H, R, numMeas, x+, P+; S is formed by the gate
    kfMeas = STRUCT_SPF.kfMeasFromUpdate(y, H_all(:, :, epoch), V, num_meas, ...
        prior, priorP, post, postP);

    % ---- gate ----
    [sys, spoofTel] = SPF_protectedNav(sys, kfMeas, propTel, epoch);
    mode = spoofTel.info.mode;

    % ---- host: apply the gate's decisions to the OPERATIONAL KF ----
    if kfUpdated
        kf_x = post; kf_P = postP;                       % normal closed-loop update
    else
        [kf_x, kf_P] = SPF_insCoast(kf_x, kf_P, propTel);    % extrapolated only (coast)
    end
    if spoofTel.nav.applyCorrection                      % LATCH or COMMIT: setKF(nav.state, nav.covar)
        kf_x = spoofTel.nav.state;
        kf_P = spoofTel.nav.covar;
    end
    if spoofTel.kfCommand.startTrial                       % probation opens
        tr_x = kf_x; tr_P = kf_P;
    end

    % ---- per-epoch logs ----
    out.x_nav(:, epoch)         = kf_x;
    out.sig_pos(:, epoch)       = sqrt(max(diag(kf_P(1:3, 1:3)), 0));
    out.sig_z(epoch)            = out.sig_pos(prm.idx_z, epoch);
    out.state(epoch)            = spoofTel.info.mode;
    out.SS_alarm(epoch)         = spoofTel.info.ssAlarm;
    out.CPI_alarm(epoch)        = spoofTel.info.cpiAlarm;
    out.alarm_axis(epoch, :)    = spoofTel.info.alarmPerAxis;
    out.SS_PL(epoch)            = spoofTel.info.maxProtectionLevel;
    out.q_reval(epoch)          = spoofTel.info.qReval;
    out.reval_computed(epoch)   = spoofTel.info.revalComputed;
    out.dwell(epoch)            = spoofTel.info.dwellCount;
    out.anchor_missing(epoch)   = spoofTel.info.anchorMissing;
    out.coast_epochs(epoch)     = spoofTel.info.coastEpochs;
    out.correction_log(epoch)   = spoofTel.nav.applyCorrection;

    % ---- event logs ----
    if spoofTel.info.eventLatched
        event_detect(end+1) = epoch;                                    %#ok<AGROW>
        event_anchor(end+1) = double(spoofTel.info.eventAnchorEpoch);   %#ok<AGROW>
    end
    if spoofTel.info.eventProbationStarted, event_prob_start(end+1) = epoch; end %#ok<AGROW>
    if spoofTel.info.eventProbationVetoed,  event_prob_fail(end+1)  = epoch; end %#ok<AGROW>
    if spoofTel.info.eventHandback,         event_handback(end+1)   = epoch; end %#ok<AGROW>
end

out.ev.t_detect     = event_detect;
out.ev.t_anchor     = event_anchor;
out.ev.t_prob_start = event_prob_start;
out.ev.t_prob_fail  = event_prob_fail;
out.ev.t_handback   = event_handback;
end
