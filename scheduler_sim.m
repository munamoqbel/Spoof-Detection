function out = scheduler_sim(z_all, H_all, V, Phi_step, Q_step, ...
                             ticks_per_epoch, x0, P0)
%SCHEDULER_SIM  Test environment for the TWO-FUNCTION architecture.
%
% Emulates the real scheduler: a 100 Hz loop (the ins_100hz_template
% steps, implemented concretely for the linear harness) that calls
% spoof_monitor_2hz once every ticks_per_epoch ticks. This exercises
% everything the flat harness (recovery_nav_sim) cannot:
%
%   - spoof_monitor_2hz itself: persistent sys/cfg/epoch, reset_request
%   - Phi_acc / Q_acc accumulation and per-interval reset
%   - the feedback gate (ON only in NOMINAL)
%   - reset posting (latch / commit) and consume-at-next-tick timing
%   - the mechanization tracking the FSM output across all modes
%
% 100 Hz MODEL (linear-harness equivalents of the template):
%   mechanization        x_mech <- Phi_step * x_mech            (per tick)
%   accumulators         Phi_acc <- Phi_step*Phi_acc
%                        Q_acc   <- Phi_step*Q_acc*Phi_step' + Q_step
%   feedback (closed loop, linear equivalent): at each epoch boundary,
%                        if enabled, x_mech <- nav.state
%   reset                posted at boundary, CONSUMED AT THE NEXT TICK
%   (the 100 Hz P propagation is omitted: in this harness the 2 Hz KF
%    performs its own time update from Phi_acc/Q_acc, so a separate
%    100 Hz P would be redundant for the test)
%
% INPUTS
%   z_all [m x N], H_all [m x n x N], V [m x m]   the logged scenario
%   Phi_step, Q_step [n x n]   PER-TICK matrices (e.g. interval/50)
%   ticks_per_epoch            e.g. 50  (100 Hz / 2 Hz)
%   x0, P0                     initial state / covariance for init
%
% OUTPUT struct 'out' (recovery_nav_sim-compatible fields plus extras):
%   x_nav, sig_pos, state, CPI_alarm, SS_alarm, alarm_axis, SS_PL,
%   q_reval, dwell, ev.*                          (as recovery_nav_sim)
%   x_mech     [n x N]  mechanization state at each epoch boundary
%                       (AFTER feedback application)
%   feedback   [1 x N]  feedback flag applied at each boundary
%   reset_log  [1 x N]  reset_flag returned at each epoch
%   Phi_acc, Q_acc      the accumulated interval matrices (constant
%                       here since Phi_step is constant) - pass these
%                       to recovery_nav_sim for the equivalence check

num_states = size(Phi_step, 1);
num_epochs = size(z_all, 2);

% ---------------- 100 Hz side state (template persistents,
%                  local here so batch runs stay repeatable) ----------
x_mech           = x0;
Phi_acc          = eye(num_states);
Q_acc            = zeros(num_states);
feedback_enabled = true;
pending_reset    = false;
pending_state    = zeros(num_states, 1);

% ---------------- logs ----------------
out.x_nav      = zeros(num_states, num_epochs);
out.sig_pos    = zeros(3, num_epochs);
out.state      = zeros(1, num_epochs);
out.CPI_alarm  = false(num_epochs, 1);
out.SS_alarm   = false(num_epochs, 1);
out.alarm_axis = false(num_epochs, 3);
out.SS_PL      = zeros(num_epochs, 1);
out.q_reval    = zeros(num_epochs, 1);
out.reval_computed = false(num_epochs, 1);
out.dwell      = zeros(num_epochs, 1);
out.x_mech     = zeros(num_states, num_epochs);
out.feedback   = false(1, num_epochs);
out.reset_log  = false(1, num_epochs);

event_detect     = [];
event_anchor     = [];
event_prob_start = [];
event_prob_fail  = [];
event_handback   = [];

kf_x = x0;
kf_P = P0;
num_meas =  size(H_all(:, :, 1), 1);

% ======================================================================
%  the scheduler
% ======================================================================
for epoch = 1:num_epochs

    % ---------------- 100 Hz loop: ticks within this interval --------
    for tick = 1:ticks_per_epoch

        % template step 0: consume a pending reset FIRST
        if pending_reset
            x_mech        = pending_state;
            pending_reset = false;
        end

        % template step 1: mechanization (linear equivalent)
        x_mech = Phi_step * x_mech;

        % template step 2: accumulate the interval matrices
        Phi_acc = Phi_step * Phi_acc;
        Q_acc   = Phi_step * Q_acc * Phi_step' + Q_step;

        % template step 3: feedback is applied at the epoch boundary
        % in this linear model (see below), not per tick.
    end

    % ---------------- 2 Hz call (the GNSS epoch) ---------------------
    reset_request = (epoch == 1);              % fresh persistents at start
    spoofInfo.phiAcc = Phi_acc;
    spoofInfo.qAcc = Q_acc;

    [kf_x, kf_P, innov, innov_S, ~, kf_prior] = kalman_update_step(kf_x, kf_P, ...
        z_all(:, epoch), H_all(:, :, epoch), spoofInfo, V);

    [mode_flag, reset_flag, reset_state, reseed, spoofTel] = ...
        spoof_monitor_2hz(innov, innov_S, H_all(:, :, epoch), num_meas, ...
        kf_prior, kf_x, kf_P, V, spoofInfo, x0, P0, reset_request);
    if reseed.reseedKF % your step 2d
        kf_x = reseed.reseedState;
        kf_P = reseed.reseedCov;
    end



    % ---------------- template step 4: handoff back down -------------
    feedback_enabled = (mode_flag == CST_spfMode.NOMINAL);

    if reset_flag
        pending_reset = true;                    % consumed at the NEXT tick
        pending_state = reset_state;
    end

    if feedback_enabled
        x_mech = spoofTel.nav.state;                       % closed-loop correction
    end

    out.spoofInfo.phiAcc = Phi_acc;          % constant here; keep last
    out.spoofInfo.qAcc   = Q_acc;
    Phi_acc = eye(num_states);                 % restart the interval
    Q_acc   = zeros(num_states);

    % ---------------- logs ----------------
    out.x_nav(:, epoch)         = spoofTel.nav.state;
    out.sig_pos(:, epoch)       = spoofTel.nav.sigmaPosition;
    out.state(epoch)            = spoofTel.info.mode;
    out.SS_alarm(epoch)         = spoofTel.info.ssAlarm;
    out.CPI_alarm(epoch)        = spoofTel.info.cpiAlarm;
    out.alarm_axis(epoch, :)    = spoofTel.info.alarmPerAxis;
    out.SS_PL(epoch)            = spoofTel.info.maxProtectionLevel;
    out.q_reval(epoch)          = spoofTel.info.qReval;
    out.reval_computed(epoch)   = spoofTel.info.revalComputed;
    out.dwell(epoch)            = spoofTel.info.dwellCount;
    out.x_mech(:, epoch)        = x_mech;
    out.feedback(epoch)         = feedback_enabled;
    out.reset_log(epoch)        = reset_flag;

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
