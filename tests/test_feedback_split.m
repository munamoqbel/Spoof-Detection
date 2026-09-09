%% test_feedback_split.m
% Host feedback split (setKF: stateFB = GAIN*x+, states = x+ - GAIN*x+) on
% the harness scenario, via tests/hostSplitSim.m:
%   1. INVARIANCE: the split host's ESTIMATE (mechanization + residual)
%      equals the bare-state host of recovery_nav_sim at every epoch, with
%      identical events, for GAIN = 0.4, 0.1 and 1.0 (full feedback).
%   2. DRAIN IN COAST: with the split applied to the extrapolated state on
%      COAST epochs, the residual left by the latch decays as (1-GAIN)^k and
%      the navigation OUTPUT converges onto the estimate (the clean coast).
%   3. NO DRAIN: the residual persists (and grows through Phi) for the whole
%      coast, so the output never receives that part of the correction.
%   4. COMMIT: the hand-back through the normal setKF drains in NOMINAL.
% Run from the repo root (MATLAB, or Octave with tools/octave_shim).

prm_boot = kujur_params();
prm = kujur_params(double(CST_spfParam.WINDOW_LENGTH), CST_spfParam.CPI_THRESHOLD, ...
    CST_spfParam.K_FALSE_ALERT, CST_spfParam.K_MISSED_DETECTION);
[~, ~, H_all, ~, ~, Phi, Q, scn, z_all, V, P0] = generate_test_data(prm_boot, 42, 0.10, 0.5, [1 1 1]);
propTel = STRUCT_SPF.setPropTel(Phi, Q);
ref = recovery_nav_sim(z_all, H_all, V, propTel, P0, prm);
GAIN = 0.4;

%% 1. invariance of the estimate to the feedback gain
fprintf('Split host vs bare host (estimate = mechanization + residual) ...\n');
ok1 = true;
for g = [GAIN 0.1 1.0]
    h = hostSplitSim(z_all, H_all, V, propTel, P0, prm, g, true, false);
    dmax = max(max(abs(h.estimate - ref.x_nav)));
    evOK = isequal(h.ev.t_detect, ref.ev.t_detect) && isequal(h.ev.t_prob_start, ref.ev.t_prob_start) ...
        && isequal(h.ev.t_prob_fail, ref.ev.t_prob_fail) && isequal(h.ev.t_handback, ref.ev.t_handback);
    okg = (dmax < 1e-9) && evOK;
    ok1 = ok1 && okg;
    fprintf('  GAIN = %.1f: max |estimate - x_nav| = %.2e, same events = %d -> %s\n', g, dmax, evOK, pf(okg));
end
fprintf('\n');

%% 2./3. drain vs no drain during COAST
drain   = hostSplitSim(z_all, H_all, V, propTel, P0, prm, GAIN, true,  false);
nodrain = hostSplitSim(z_all, H_all, V, propTel, P0, prm, GAIN, false, false);
td = drain.ev.t_detect(1);  th = drain.ev.t_handback(1);
rPos = @(h, k) norm(h.residual(1:3, k));
gap  = @(h, k) norm(h.output(1:3, k) - h.estimate(1:3, k));   % output - estimate = -residual

fprintf('Latch at epoch %d: |nav.state| pos = %.4f m, residual left by the normal setKF = %.4f m\n', ...
    td, norm(drain.latchState(1:3)), rPos(drain, td));
fprintf('  epoch   drain: |residual| output-estimate   no drain: |residual| output-estimate\n');
for kk = [0 1 3 5 10 20 40]
    fprintf('  +%3d    %.5f     %.5f            %.5f     %.5f\n', kk, ...
        rPos(drain, td+kk), gap(drain, td+kk), rPos(nodrain, td+kk), gap(nodrain, td+kk));
end
ok2 = (rPos(drain, td+10) < 0.02 * rPos(drain, td)) && (gap(drain, td+20) < 1e-3 * max(rPos(drain, td), 1e-6) + 1e-6);
ok3 = rPos(nodrain, td+10) > 0.5 * rPos(nodrain, td);
fprintf('  drain: residual at +10 < 2%% of latch value and output on the estimate by +20 -> %s\n', pf(ok2));
fprintf('  no drain: residual at +10 still > 50%% of latch value (never reaches the output) -> %s\n\n', pf(ok3));

%% 4. commit through the normal setKF drains in NOMINAL
fprintf('Commit at epoch %d: residual %.5f m -> +10: %.5f m\n', th, rPos(drain, th), rPos(drain, th+10));
ok4 = rPos(drain, th+10) < 0.05 * max(rPos(drain, th), 1e-9) + 1e-6;
fprintf('  hand-back leftover drained by the normal updates -> %s\n\n', pf(ok4));

if ok1 && ok2 && ok3 && ok4
    fprintf('=== test_feedback_split PASS ===\n');
else
    fprintf('=== test_feedback_split FAIL ===\n');
end

function s = pf(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end
