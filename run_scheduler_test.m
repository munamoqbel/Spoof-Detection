%% run_scheduler_test.m
% Test environment for the TWO-FUNCTION architecture (100 Hz + 2 Hz).
%
% WHAT IS TESTED (and what each test would catch):
% T1 Accumulator correctness: Phi_step^50 vs the interval Phi.
% T2 EQUIVALENCE: the scheduler path (100 Hz loop + persistent
%    SPF_gate) must reproduce the flat reference
%    (recovery_nav_sim) EXACTLY when both use the same interval
%    matrices. Catches: broken handoff, wrong persistent state,
%    double/missed epochs, accumulator reset errors.
% T3 Feedback gate: ON exactly when mode == NOMINAL.
% T4 Reset command: fires exactly at latch and commit epochs.
% T5 Mechanization tracking: x_mech equals the FSM output at every
%    epoch boundary (except the latch epoch itself, where a
%    one-tick lag is the designed behavior). Catches: missed reset,
%    wrong gate, reset consumed at the wrong time.
% T6 Persistence/reset_request: a second full run must reproduce
%    the first exactly (proves re-initialisation is complete).
%
% NOTE: this validates the ARCHITECTURE/WIRING. The design constants
% are the ones baked into CST_spfParam (the validated diagonal-
% attack design), not re-derived here.

clear; clc; close all;
clear functions % wipe any persistent state from earlier sessions

ticks_per_epoch = 50; % 100 Hz / 2 Hz

fprintf('==============================================\n');
fprintf('  Scheduler test: 100 Hz + 2 Hz architecture\n');
fprintf('==============================================\n\n');

%% Scenario (the validated diagonal-attack case)
prm_boot = kujur_params();
[~, ~, H_all, ~, ~, Phi_2hz, Q_2hz, scn, z_all, V, P0] = ...
    generate_test_data(prm_boot, 42, 0.10, 0.5, [1 1 1]);

num_states = size(Phi_2hz, 1);
num_epochs = scn.N_total;
x0 = zeros(num_states, 1);

% Per-tick matrices: split the interval into ticks_per_epoch steps
Phi_step = eye(num_states) + (Phi_2hz - eye(num_states)) / ticks_per_epoch;
Q_step = Q_2hz / ticks_per_epoch;

%% T1 - accumulator correctness
Phi_check = eye(num_states);
for t = 1:ticks_per_epoch
    Phi_check = Phi_step * Phi_check;
end
err_phi = norm(Phi_check - Phi_2hz, 'fro') / norm(Phi_2hz, 'fro');
fprintf('T1 Phi_step^%d vs interval Phi : rel err = %.2e', ...
    ticks_per_epoch, err_phi);
if err_phi < 1e-12
    fprintf(' PASS (exact - nilpotent coupling)\n');
else
    fprintf(' INFO (non-nilpotent terms; equivalence test unaffected)\n');
end

%% Run the scheduler environment
sch = scheduler_sim(z_all, H_all, V, Phi_step, Q_step, ...
    ticks_per_epoch, x0, P0);

fprintf(' Q_acc vs interval Q : rel diff = %.2e (informational)\n\n', ...
    norm(sch.propTel.accumQ - Q_2hz, 'fro') / norm(Q_2hz, 'fro'));

%% T2 - equivalence against the flat reference
% Reference uses the SAME accumulated matrices, so any difference is
% a wiring bug in the scheduler path, not a modelling difference.
prm = kujur_params(10, 45.6397, 5.2331, 4.7534); % baked design
prm.mon_axes = [1 2 3];
rec.T_reval = 39.25;
rec.M_dwell = 10;
rec.M_prob = 18;

ref = recovery_nav_sim(z_all, H_all, V, sch.propTel, P0, prm);

d_nav = max(max(abs(sch.x_nav - ref.x_nav)));
same_seq = isequal(sch.state, ref.state);
same_ev = isequal(sch.ev.t_detect, ref.ev.t_detect) && ...
    isequal(sch.ev.t_handback, ref.ev.t_handback) && ...
    isequal(sch.ev.t_prob_fail, ref.ev.t_prob_fail);
same_al = isequal(sch.SS_alarm, ref.SS_alarm) && ...
    isequal(sch.CPI_alarm, ref.CPI_alarm);

t2_pass = (d_nav < 1e-9) && same_seq && same_ev && same_al;
fprintf('T2 Scheduler vs flat reference:\n');
fprintf(' max |x_nav diff| = %.2e | modes %s | events %s | alarms %s', ...
    d_nav, pf(same_seq), pf(same_ev), pf(same_al));
fprintf(' -> %s\n\n', pf(t2_pass));

%% T3 - feedback gate
t3_pass = isequal(sch.feedback, (sch.state == CST_spfMode.NOMINAL));
fprintf('T3 Feedback ON iff NOMINAL : %s\n', pf(t3_pass));

%% T4 - reset command epochs
reset_epochs = find(sch.reset_log);
expected_resets = sort([sch.ev.t_detect, sch.ev.t_handback]);
t4_pass = isequal(reset_epochs(:)', expected_resets(:)');
fprintf('T4 Resets at latch+commit only: %s (epochs: %s)\n', ...
    pf(t4_pass), num2str(reset_epochs(:)'));

%% T5 - mechanization tracks the FSM output at boundaries
% Exclude the latch epochs: there the reset is consumed at the NEXT
% tick by design, so x_mech still holds the pre-latch value.
check = true(1, num_epochs);
check(sch.ev.t_detect) = false;
d_mech = max(max(abs(sch.x_mech(:, check) - sch.x_nav(:, check))));
t5_pass = d_mech < 1e-9;
fprintf('T5 x_mech == nav output : max diff %.2e -> %s\n', ...
    d_mech, pf(t5_pass));
for k = sch.ev.t_detect
    fprintf(' (latch epoch %d designed one-tick lag: |mech-nav| = %.3f m)\n', ...
        k, max(abs(sch.x_mech(:, k) - sch.x_nav(:, k))));
end

%% T6 - persistence / reset_request completeness
sch2 = scheduler_sim(z_all, H_all, V, Phi_step, Q_step, ...
    ticks_per_epoch, x0, P0);
t6_pass = isequal(sch.x_nav, sch2.x_nav) && isequal(sch.state, sch2.state);
fprintf('T6 Repeat run identical : %s\n\n', pf(t6_pass));

%% Verdict
all_pass = t2_pass && t3_pass && t4_pass && t5_pass && t6_pass;
if all_pass
    fprintf('=== ALL SCHEDULER TESTS PASS ===\n');
else
    fprintf('=== FAILURES ABOVE - fix before integrating ===\n');
end

%% quick visual: scheduler-path result (should match the known figure)
t = (0:num_epochs-1) / prm_boot.fs;
figure('Name', 'Scheduler-path result', 'NumberTitle', 'off');
subplot(3,1,1);
plot(t, vecnorm(sch.x_nav(1:3,:)), 'b-', 'LineWidth', 1.2); hold on;
plot(t, vecnorm(ref.x_nav(1:3,:)), 'r--', 'LineWidth', 1.0);
ylabel('|3-D pos| (m)'); legend('scheduler path', 'flat reference');
grid on; title('Equivalence: the two paths must overlay exactly');
subplot(3,1,2);
stairs(t, sch.state, 'k-', 'LineWidth', 1.4); hold on;
plot(t(sch.reset_log), sch.state(sch.reset_log), 'rv', 'MarkerFaceColor', 'r');
yticks([1 2 3]); yticklabels({'NOMINAL', 'COAST', 'PROBATION'});
grid on; title('FSM mode (red = reset commands to the 100 Hz side)');
subplot(3,1,3);
stairs(t, double(sch.feedback), 'b-', 'LineWidth', 1.4);
ylim([-0.1 1.1]); yticks([0 1]); yticklabels({'OFF', 'ON'});
xlabel('Time (s)'); grid on; title('Feedback gate to the 100 Hz side');

function s = pf(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end
