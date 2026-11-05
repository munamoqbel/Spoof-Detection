%% test_shadowArm.m
% The monitors in shadow through the warm-up: on a synthetic 60-state host (static truth, 14 range rows at 5 m) with
% the repository ARM_EPOCHS of 10,
%   A  a clean start arms at epoch ARM_EPOCHS and a 30 m step at 60 latches within two epochs;
%   B  a 40 m transient on epochs 4 to 7 raises no alarm flag and no latch while unarmed, shows in the logged margins,
%      restarts the count so the gate arms ARM_EPOCHS epochs after the last shadow alarm, and the 30 m step at 60 still
%      latches within two epochs;
%   C  the same transient and a 30 m step from epoch 30: the pool carried over the arming, the step latches within two
%      epochs of its onset.
% Run from the repo root (MATLAB, or Octave with tools/octave_shim).

n = CST_gnssHybrid.NO_STATES;  nArm = double(CST_spfParam.ARM_EPOCHS);
rng(5);
N = 160;  m = 14;  R = 25 * eye(m);
Phi = eye(n);  Q = zeros(n);  Q(1:3, 1:3) = 2.0 * eye(3);
propTel = STRUCT_SPF.setPropTel(Phi, Q);
el = 0.3 + 0.6 * rand(m, 1);  az = 2 * pi * rand(m, 1);
u = [cos(el) .* cos(az), cos(el) .* sin(az), -sin(el)];
H = zeros(m, n);  H(:, 1:3) = -u;
P0 = 0.01 * eye(n);  P0(1:3, 1:3) = eye(3);
function [res] = runSched(offN, n, m, N, H, R, P0, propTel)
    % offN(k): the north offset the ranges describe at epoch k
    sys = STRUCT_SPF.zeroSys;  sys.coastCov = P0;
    sys.anchor = STRUCT_SPF.setAnchor(true, zeros(n, 1), P0, uint32(0));
    kf_x = zeros(n, 1);  kf_P = P0;  tr_x = kf_x;  tr_P = P0;  mode = CST_spfMode.NOMINAL;
    armed = false(1, N);  alarm = false(1, N);  marg = zeros(1, N);  det = [];  hb = [];  st = zeros(1, N);
    for k = 1:N
        z = 5 * randn(m, 1) + H(:, 1:3) * [offN(k) 0 0]';
        inProb = (mode == CST_spfMode.PROBATION);  kfUpd = (mode == CST_spfMode.NOMINAL);
        if inProb
            [tr_x, tr_P, y, ~, ~, xp, xpP] = kalman_update_step(tr_x, tr_P, z, H, propTel, R);  post = tr_x;  postP = tr_P;
        else
            [ax, aP, y, ~, ~, xp, xpP] = kalman_update_step(kf_x, kf_P, z, H, propTel, R);  post = ax;  postP = aP;
        end
        kfMeas = STRUCT_SPF.kfMeasFromUpdate(y, H, R, m, xp, xpP, post, postP);
        [sys, tel] = SPF_protectedNav(sys, kfMeas, propTel, k);
        mode = tel.info.mode;
        if kfUpd, kf_x = post;  kf_P = postP; else, [kf_x, kf_P] = SPF_insCoast(kf_x, kf_P, propTel); end
        if tel.nav.applyCorrection, kf_x = tel.nav.state;  kf_P = tel.nav.covar; end
        if tel.kfCommand.startTrial, tr_x = kf_x;  tr_P = kf_P; end
        armed(k) = tel.info.armed;  alarm(k) = tel.info.ssAlarm || tel.info.cpiAlarm;  st(k) = double(mode);
        marg(k) = max([tel.info.ssRatio, tel.info.cpiRatio]);
        if tel.info.eventLatched, det(end+1) = k; end %#ok<AGROW>
        if tel.info.eventHandback, hb(end+1) = k; end %#ok<AGROW>
    end
    res = struct('armed', armed, 'alarm', alarm, 'marg', marg, 'det', det, 'hb', hb, 'st', st);
end
offA = zeros(1, N);  offA(60:120) = 30;
offB = offA;  offB(4:7) = 40;
offC = zeros(1, N);  offC(4:7) = 40;  offC(30:90) = 30;
rA = runSched(offA, n, m, N, H, R, P0, propTel);
rB = runSched(offB, n, m, N, H, R, P0, propTel);
rC = runSched(offC, n, m, N, H, R, P0, propTel);
kArmA = find(rA.armed, 1);  kArmB = find(rB.armed, 1);  kArmC = find(rC.armed, 1);
lastShadowB = find(rB.marg > 1 & ~rB.armed, 1, 'last');
fprintf('A clean: armed from epoch %d (ARM_EPOCHS %d), latch %s, hand-backs %s\n', kArmA, nArm, mat2str(rA.det), mat2str(rA.hb));
fprintf('B transient 4..7: shadow margins above 1 on epochs %s, alarm flags while unarmed %d, armed from epoch %d, latch %s\n', ...
    mat2str(find(rB.marg > 1 & ~rB.armed)), nnz(rB.alarm & ~rB.armed), kArmB, mat2str(rB.det));
fprintf('C transient 4..7, step from 30: armed from epoch %d, latch %s\n', kArmC, mat2str(rC.det));
okA = kArmA == nArm && ~isempty(rA.det) && rA.det(1) >= 60 && rA.det(1) <= 61 && all(rA.st(1:59) == 1) && ~any(rA.hb <= 120);
okB = ~isempty(lastShadowB) && lastShadowB >= 4 && kArmB >= lastShadowB + nArm && kArmB <= lastShadowB + nArm + 2 ...
    && ~any(rB.alarm & ~rB.armed) && ~any(rB.det < 60) && all(rB.st(1:59) == 1) && ~isempty(rB.det) && rB.det(1) >= 60 && rB.det(1) <= 61;
okC = kArmC < 30 && ~isempty(rC.det) && rC.det(1) >= 30 && rC.det(1) <= 31 && ~any(rC.hb <= 90);
fprintf('A %s | B %s | C %s\n', pf(okA), pf(okB), pf(okC));
if okA && okB && okC, fprintf('=== test_shadowArm PASS ===\n'); else, fprintf('=== test_shadowArm FAIL ===\n'); end

function s = pf(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end
