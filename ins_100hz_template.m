function ins_100hz_template()
%INS_100HZ_TEMPLATE  What YOUR 100 Hz INS function must add. TEMPLATE -
%                    merge these pieces into your mechanization; the
%                    TODO lines mark your existing code.
%
% Additions (three items only):
%   1. Phi_acc / Q_acc accumulators   (interval matrices for the 2 Hz side)
%   2. feedback ON/OFF                (OFF in COAST and PROBATION)
%   3. reset consumption              (load state at latch / commit)

persistent Phi_acc Q_acc feedback_enabled ...
           pending_reset pending_reset_state

if isempty(Phi_acc)
    n = CST_gnssHybrid.NO_STATES;                    % your state dimension
    Phi_acc             = eye(n);
    Q_acc               = zeros(n);
    feedback_enabled    = true;
    pending_reset       = false;
    pending_reset_state = zeros(n, 1);
end

% ======================================================================
%  0. consume a pending reset FIRST (posted by the 2 Hz side at
%     latch or commit) - overwrite the mechanization state once
% ======================================================================
if pending_reset
    % TODO: load pending_reset_state into your mechanization
    %       (position, velocity, attitude, bias states)
    pending_reset = false;
end

% ======================================================================
%  1. your normal mechanization step (unchanged)
% ======================================================================
% TODO: integrate IMU; build this step's Phi from the nav state;
%       propagate P as you already do:
%           Q_step = Phi*Q*Phi'*dt;                % YOUR discretization form
%           P      = Phi*P*Phi' + Q_step;

% ======================================================================
%  2. accumulate the interval matrices (three lines)
% ======================================================================
% Phi_acc = Phi * Phi_acc;
% Q_acc   = Phi * Q_acc * Phi' + Q_step;        % same Q_step as your P line

% ======================================================================
%  3. correction feedback (your existing code, now gated)
% ======================================================================
if feedback_enabled
    % TODO: apply the 2 Hz correction into the mechanization (existing)
end

% ======================================================================
%  4. at each GNSS epoch (every ~50th call): hand off and reset
% ======================================================================
% TODO: pass spoofInfo.phiAcc = Phi_acc, spoofInfo.qAcc = Q_acc to the
%       2 Hz function (spoof_monitor_2hz), receive
%       (mode_flag, reset_flag, reset_state), then:
%
%   feedback_enabled = (mode_flag == CST_spfMode.NOMINAL);
%   if reset_flag
%       pending_reset       = true;                 % consumed at next cycle
%       pending_reset_state = reset_state;
%   end
%   Phi_acc = eye(n);                               % restart the interval
%   Q_acc   = zeros(n);

end
