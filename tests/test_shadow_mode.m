%% test_shadow_mode.m
% Implementation Guide Step 5, items 1-2, on the harness scenario:
%   1. SHADOW MODE on clean data: the gate is called every epoch but its
%      outputs are ignored (host updates normally). Expect zero alarms and
%      a sub-metre protection level throughout.
%   2. SCRIPTED ATTACK with authority (recovery_nav_sim host contract):
%      expect a latch within 2 epochs of onset on every monitored axis,
%      every probation opened while the spoofer is on vetoed, exactly one
%      commit after it stops, and a bounded position error.
% Run from the repo root (MATLAB, or Octave with tools/octave_shim).

prm_boot = kujur_params();
n = prm_boot.n_states;
prm = kujur_params(double(CST_spfParam.WINDOW_LENGTH), CST_spfParam.CPI_THRESHOLD, ...
    CST_spfParam.K_FALSE_ALERT, CST_spfParam.K_MISSED_DETECTION);

%% 1. shadow mode, clean data (no attack injected)
fprintf('Shadow mode (clean data, gate has no authority) ...\n');
[~, ~, H_all, ~, ~, Phi, Q, scn, z_all, V, P0] = generate_test_data(prm_boot, 5, 0.0, 0.0, [0 0 1]);
propTel = STRUCT_SPF.setPropTel(Phi, Q);
N = scn.N_total; m = size(z_all, 1);
kf_x = zeros(n, 1); kf_P = P0;
nAlarm = 0; plMax = 0; modeOK = true;
for k = 1:N
    [kf_x, kf_P, y, ~, ~, xp, xpP] = kalman_update_step(kf_x, kf_P, z_all(:, k), H_all(:, :, k), propTel, V);
    kfMeas = STRUCT_SPF.kfMeasFromUpdate(y, H_all(:, :, k), V, m, xp, xpP, kf_x, kf_P);
    tel = spoofMonitor2hz(kfMeas, propTel, true, k == 1);  % outputs ignored (shadow)
    nAlarm = nAlarm + (tel.info.ssAlarm || tel.info.cpiAlarm);
    plMax  = max(plMax, tel.info.maxProtectionLevel);
    modeOK = modeOK && (tel.info.mode == CST_spfMode.NOMINAL);
end
ok1 = (nAlarm == 0) && modeOK && (plMax < 1.0);
fprintf('  alarms = %d (expect 0) | stayed NOMINAL = %d | max PL = %.3f m -> %s\n', ...
    nAlarm, modeOK, plMax, pf(ok1));

% alignment: navActive = false must return zeroTel and re-arm the gate; the
% first active call afterwards re-initialises (no alarm, coastEpochs 0)
alignOK = true;
for k = 1:5
    tel = spoofMonitor2hz(kfMeas, propTel, false, false);
    alignOK = alignOK && (tel.info.mode == CST_spfMode.NOMINAL) && ~tel.nav.applyCorrection ...
        && ~tel.kfCommand.startTrial && ~tel.info.ssAlarm && ~tel.info.cpiAlarm;
end
tel = spoofMonitor2hz(kfMeas, propTel, true, false);       % re-entry, no resetRequest
alignOK = alignOK && (tel.info.mode == CST_spfMode.NOMINAL) && ~tel.info.ssAlarm ...
    && ~tel.info.cpiAlarm && (tel.info.coastEpochs == 0);
fprintf('  alignment (navActive = false): zeroTel, then clean re-init on re-entry -> %s\n\n', pf(alignOK));
ok1 = ok1 && alignOK;

%% 2. scripted attack, gate has authority
fprintf('Scripted attack (oblique ramp + dither, gate has authority) ...\n');
[~, ~, H_all, xh_base, ~, Phi, Q, scn, z_all, V, P0] = generate_test_data(prm_boot, 42, 0.10, 0.5, [1 1 1]);
propTel = STRUCT_SPF.setPropTel(Phi, Q);
out = recovery_nav_sim(z_all, H_all, V, propTel, P0, prm);
err = vecnorm(out.x_nav(1:3, :));
Nw  = double(CST_spfParam.WINDOW_LENGTH);

latched   = ~isempty(out.ev.t_detect);
td        = out.ev.t_detect;
onsetOK   = latched && td(1) >= scn.spoof_start && td(1) <= scn.spoof_start + 2;
axesOK    = latched && all(out.alarm_axis(td(1), :));
earlyProb = out.ev.t_prob_start(out.ev.t_prob_start <= scn.spoof_end);
vetoOK    = numel(out.ev.t_prob_fail) >= numel(earlyProb) - 0 || isempty(earlyProb);   % every early probation vetoed
commitOK  = numel(out.ev.t_handback) == 1 && out.ev.t_handback(1) > scn.spoof_end;
errOK     = max(err(scn.spoof_start:end)) < 1.0;
unprot    = max(vecnorm(xh_base(1:3, scn.spoof_start:end)));
ok2 = onsetOK && axesOK && vetoOK && commitOK && errOK;
fprintf('  latch epoch %s (onset %d) -> %s | all 3 axes -> %s\n', mat2str(td), scn.spoof_start, pf(onsetOK), pf(axesOK));
fprintf('  probations while spoofed: %d started, %d vetoed -> %s\n', numel(earlyProb), numel(out.ev.t_prob_fail), pf(vetoOK));
fprintf('  commits: %s (spoof ends %d) -> %s\n', mat2str(out.ev.t_handback), scn.spoof_end, pf(commitOK));
fprintf('  max |3-D err| after onset: %.3f m protected vs %.3f m unprotected -> %s\n\n', ...
    max(err(scn.spoof_start:end)), unprot, pf(errOK));

if ok1 && ok2
    fprintf('=== test_shadow_mode PASS ===\n');
else
    fprintf('=== test_shadow_mode FAIL ===\n');
end

function s = pf(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end
