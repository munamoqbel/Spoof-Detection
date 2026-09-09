function out = scheduler_sim(z_all, H_all, V, Phi_step, Q_step, ...
                             ticks_per_epoch, x0, P0)
%SCHEDULER_SIM  Test environment for the TWO-FUNCTION architecture.
%
% Emulates the real scheduler: a 100 Hz loop (the ins_100hz_template
% steps, implemented concretely for the linear harness) that calls the
% persistent 2 Hz gate spoofMonitor2hz once every ticks_per_epoch ticks,
% with a HOST 2 Hz side that follows docs/HOST_2HZ_WIRING.m. This
% exercises everything the flat harness (recovery_nav_sim) cannot:
%
%   - spoofMonitor2hz itself: persistent sys/epoch, resetRequest
%   - accumPhi / accumQ accumulation and per-interval reset (spfAccumProp)
%   - the feedback gate (ON only in NOMINAL)
%   - corrections at latch / commit consumed by the mechanization at the
%     NEXT tick
%   - the mechanization tracking the protected solution across all modes
%
% 100 Hz MODEL (linear-harness equivalents of the template):
%   mechanization        x_mech <- Phi_step * x_mech            (per tick)
%   accumulators         accumPhi <- Phi_step*accumPhi
%                        accumQ   <- Phi_step*accumQ*Phi_step' + Q_step
%   feedback (closed loop, linear equivalent): at each epoch boundary,
%                        if enabled, x_mech <- kf_x
%   reset                posted at boundary, CONSUMED AT THE NEXT TICK
%
% OUTPUT struct 'out' (recovery_nav_sim-compatible fields plus extras):
%   x_nav, sig_pos, state, CPI_alarm, SS_alarm, alarm_axis, SS_PL,
%   q_reval, reval_computed, dwell, ev.*         (as recovery_nav_sim)
%   x_mech     [n x N]  mechanization state at each epoch boundary
%   feedback   [1 x N]  feedback flag applied at each boundary
%   reset_log  [1 x N]  correction (latch/commit) returned at each epoch
%   propTel             the accumulated interval matrices (constant here)

num_states = size(Phi_step, 1);
num_epochs = size(z_all, 2);
num_meas   = size(z_all, 1);

% ---------------- 100 Hz side state (template persistents) ----------
x_mech           = x0;
accumPhi         = eye(num_states);
accumQ           = zeros(num_states);
feedback_enabled = true;
pending_reset    = false;
pending_state    = zeros(num_states, 1);

% ---------------- host 2 Hz side state ----------------
kf_x = x0;  kf_P = P0;
tr_x = x0;  tr_P = P0;
mode = CST_spfMode.NOMINAL;

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

event_detect = []; event_anchor = []; event_prob_start = [];
event_prob_fail = []; event_handback = [];

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

        % template step 2: accumulate the interval matrices (spfAccumProp)
        accumPhi = Phi_step * accumPhi;
        accumQ   = Phi_step * accumQ * Phi_step' + Q_step;
    end

    % ---------------- 2 Hz call (the GNSS epoch) ---------------------
    propTel = STRUCT_SPF.setPropTel(accumPhi, accumQ);
    inProbation = (mode == CST_spfMode.PROBATION);
    kfUpdated   = (mode == CST_spfMode.NOMINAL);

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

    spoofTel = spoofMonitor2hz(kfMeas, propTel, epoch == 1);
    mode = spoofTel.info.mode;

    if kfUpdated
        kf_x = post; kf_P = postP;
    else
        [kf_x, kf_P] = insCoast(kf_x, kf_P, propTel);
    end
    if spoofTel.nav.applyCorrection                      % LATCH or COMMIT: setKF(nav.state, nav.covar)
        kf_x = spoofTel.nav.state;
        kf_P = spoofTel.nav.covar;
    end
    if spoofTel.kfCommand.startTrial
        tr_x = kf_x; tr_P = kf_P;
    end

    % ---------------- template step 4: handoff back down -------------
    feedback_enabled = (mode == CST_spfMode.NOMINAL);
    reset_flag = spoofTel.nav.applyCorrection;

    if reset_flag
        pending_reset = true;                    % consumed at the NEXT tick
        pending_state = kf_x;
    end

    if feedback_enabled
        x_mech = kf_x;                           % closed-loop correction
    end

    out.propTel = propTel;                       % constant here; keep last
    accumPhi = eye(num_states);                  % restart the interval
    accumQ   = zeros(num_states);

    % ---------------- logs ----------------
    out.x_nav(:, epoch)         = kf_x;
    out.sig_pos(:, epoch)       = sqrt(max(diag(kf_P(1:3, 1:3)), 0));
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
