%% QUICKSTART_integration.m   (read-me file - do not run)
% How to put the spoofing monitor + recovery into YOUR EKF.
% Uses the modular files. Details: Implementation_Guide_v2.pdf
%
% Your setup: INS @ 100 Hz  |  GNSS EKF @ 2 Hz  |  NED  |  closed-loop.
%
% =======================================================================
%  STEP 0 - CHECK YOUR EKF (once)
% =======================================================================
% 1. Find where the correction feedback is applied to the INS.
% 2. Find your position indices (N, E, D). Example: 1, 2, 3.
% 3. Test: add +10 m to z along H(:,idxD). x_hat(idxD) must move +10.
%
% =======================================================================
%  STEP 1 - LOG (in your 2 Hz EKF, every update)
% =======================================================================
%   gamma = z - h(x_bar)        innovation
%   S     = H*P_bar*H' + V
%   H, x_hat, P_hat, Phi, Q, z, V
%
% =======================================================================
%  STEP 2 - DESIGN OFFLINE (once, on desktop, then HARD-CODE)
% =======================================================================
% Replay a clean run that includes your WORST sky (masking, few sats).
%   per axis a = 1,2,3:
%     s2(k)    = f'*(S\f),  f = H(:,a)
%     s2min(a) = min over all k
%     N_a(a)   = solve_N_min(s2min(a), 0.10, PFAp/3, PFAm/3, 1e-7, 1e-6)
%   N     = max(N_a) + 1                              (+1 = safety margin)
%   T_N   = gaminv(1 - PFAp/3, N/2, 2)
%   k_FA  = norminv(1 - (PFAm/3/N)/2)
%   k_MD  = norminv(1 - 1e-6)
%   T_reval(m) = gaminv(1 - 1e-3, m/2, 2)    for each possible sat count
% Hard-code all values. Re-run only if IMU / constellation / budget change.
%
% =======================================================================
%  STEP 3 - WIRE THE RUNTIME (in your 2 Hz function)
% =======================================================================
% Copy protected_nav_step.m and use it as your template. Per mode:
%
% NOMINAL:
%   your EKF update                  (replaces kalman_update_step)
%   monitor_pool_step                windows: SS every epoch, CPI at close
%   refresh anchor                   if clean window closed + no alarm now
%   any alarm -> LATCH:              reset nav to anchor (loop ins_coast_step
%                                    to bring it to now), clear pool, COAST
%   no alarm -> monitor_pool_open_window
%
% COAST (output = INS only):
%   ins_coast_step                   propagate output
%   revalidation_test                GNSS vs coast; 10 passes in a row
%                                    -> PROBATION (trial = coast copy)
%
% PROBATION (output STAYS on coast):
%   ins_coast_step                   output
%   your EKF update on TRIAL         trial eats GNSS
%   monitor_pool_step on TRIAL       monitors judge it
%   alarm  -> VETO, back to COAST
%   N+8 quiet epochs -> COMMIT: output = trial, anchor = trial, NOMINAL
%
% RULES (never break):
%   - alarm is sticky: silence never exits COAST
%   - during PROBATION the output never touches GNSS
%   - enforce a max coast time (INS drift budget)
%
% =======================================================================
%  STEP 4 - 100 Hz INS SIDE (only two additions)
% =======================================================================
%   feedback ON/OFF flag   (OFF during COAST and PROBATION)
%   state reset command    (used at latch and at commit)
%
% =======================================================================
%  STEP 5 - TEST (in this order)
% =======================================================================
% 1. Shadow mode: monitor only, no authority. Clean run = zero alarms.
% 2. Scripted attacks (dither / ramp / both / oblique):
%    latch < 1 s, error always < 3-sigma, every early probation vetoed.
% 3. Slow-ramp sweep: find and document the detection floor.
% 4. Give the FSM authority. Coder refactor last.
%
% Extras worth adding: advisory if s2(k) < s2min ("coverage degraded");
% annunciate mode, alarm axes, PL, time-in-coast.
