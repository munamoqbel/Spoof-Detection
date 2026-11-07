%% test_coastBound.m
% The re-admission path of the gate: the coast covariance bounded by COAST_SIGMA_MAX, the CPI of the innovations
% against the coast in the dwell, and the first probation window opened on the coast. (1) STRUCT_SPF.boundCoastCov: a
% matrix whose monitored position sigmas exceed the bound comes back with them at the bound, symmetric positive
% semi-definite, the other states and the cross terms scaled by congruence; one within the bound comes back as it is.
% (2) A synthetic 60-state host (static truth, 14 range rows at 5 m with a receiver clock column on state 16 and a
% pressure-altitude row at 5 m as the host carries them, a position and clock process noise that grows the coast to 6 m
% sigma within 20 epochs, an INS that drifts 0.02 m per coast epoch):
%   A  a 30 m step on every row from epoch 100 to 300 latches, stays refused through the spoof and is handed back after
%      it, while SPF_revalidation on the unbounded coast covariance would have passed mid-spoof;
%   B  a 16 m step, which the chi-square passes even on the bounded covariance, is held by the CPI against the coast:
%      no probation while it lasts, hand-back after it;
%   C  a 16 m step from 100 to 200 and then a 30 m spoof switched on at the epoch the probation opens: the trial jumps
%      onto it in its first update, the first window vetoes that epoch, the spoof stays refused until it ends at 400;
%   D  a 35 m down step, which every row sees with the same sign: with the position bounded alone the chi-square passes
%      mid-spoof because the clock variance absorbs the common part; with the clock bias bounded as well it fails, and
%      the step is held through the spoof;
%   E  the return path: a 0.15 m per epoch east ramp from 100 to 300 (30 m) is followed without a latch; the step back
%      at 301 latches and reads as a return (the nav had wandered 30 m from the long coast, the GNSS lands on it), the
%      nav is corrected back and handed back within 40 epochs, on the truth; without the return path the clean GNSS
%      would be refused for as long as it lasts;
%   G  a slow drag of 24 m (0.12 m per epoch, 100 to 300) and then a step 25 m further (301 to 400): the step latches
%      as an onset (the GNSS lands some 45 m from the long coast) and is held; at 401 the clean GNSS comes back 24 m from
%      the nav, which the strict bound would refuse, and the return check in COAST reads it (the nav had wandered some
%      20 m from the long coast, the GNSS lands within a few metres of it): corrected back and handed back.
% Run from the repo root (MATLAB, or Octave with tools/octave_shim).

n = CST_gnssHybrid.NO_STATES;  mMax = double(CST_spfParam.MAX_MEAS);
CAP = CST_spfParam.COAST_SIGMA_MAX;

%% 1. the bound itself
P = eye(n);  P(1, 1) = 36;  P(2, 2) = 16;  P(3, 3) = 4;  P(4, 4) = 0.25;  P(16, 16) = 100;
P(1, 4) = 2;  P(4, 1) = 2;  P(1, 2) = 10;  P(2, 1) = 10;                                % cross terms N-vN and N-E
Pb = STRUCT_SPF.boundCoastCov(P, zeros(n, 1));
fl = zeros(n, 1);  fl(1) = 5;  Pf = STRUCT_SPF.boundCoastCov(P, fl);                 % a floor above the bound on N
okBound = abs(Pb(1, 1) - CAP^2) < 1e-12 && abs(Pb(2, 2) - CAP^2) < 1e-12 && Pb(3, 3) == 4 && Pb(4, 4) == 0.25 ...
    && abs(Pb(16, 16) - CAP^2) < 1e-12 && isequal(CST_spfParam.COAST_BOUND_STATES, [1 2 3 16]) ...
    && abs(Pb(1, 4) - 2 * CAP / 6) < 1e-12 && abs(Pb(1, 2) - 10 * (CAP / 6) * (CAP / 4)) < 1e-12 ...
    && isequal(Pb, Pb') && min(eig(Pb(1:4, 1:4))) > -1e-12 && isequal(STRUCT_SPF.boundCoastCov(0.5 * eye(n), zeros(n, 1)), 0.5 * eye(n)) ...
    && abs(Pf(1, 1) - 25) < 1e-12 && abs(Pf(2, 2) - CAP^2) < 1e-12 && abs(Pf(1, 2) - 10 * (5 / 6) * (CAP / 4)) < 1e-12;

%% 2. the synthetic host
rng(3);
N = 500;  m = 15;  R = 25 * eye(m);                                                      % 14 range rows and the pressure row, 5 m each as the host carries them
Phi = eye(n);  Q = zeros(n);  Q(1:3, 1:3) = 2.0 * eye(3);  Q(16, 16) = 0.5;              % m^2 per epoch on the position and the clock
propTel = STRUCT_SPF.setPropTel(Phi, Q);
el = 0.3 + 0.6 * rand(14, 1);  az = 2 * pi * rand(14, 1);
u = [cos(el) .* cos(az), cos(el) .* sin(az), -sin(el)];
H = zeros(m, n);  H(1:14, 1:3) = -u;  H(1:14, 16) = 1;  H(15, 3) = 1;                     % range rows: line of sight and the clock; the pressure row
P0 = 0.01 * eye(n);  P0(1:3, 1:3) = eye(3);  P0(16, 16) = 1;
DRIFT = 0.02;                                                                             % m per coast epoch, north
function [res] = runCase(offsetAt, tOff, stepB, tOffB, n, mMax, m, N, H, R, P0, propTel, DRIFT)
    % offsetAt(k) = the [N E D] offset the ranges describe at epoch k; stepB (north) from the epoch after the first probation
    % that opens after tOff until tOffB
    sys = STRUCT_SPF.zeroSys;  sys.coastCov = P0;
    sys.anchor = STRUCT_SPF.setAnchor(true, zeros(n, 1), P0, uint32(0));
    kf_x = zeros(n, 1);  kf_P = P0;  tr_x = kf_x;  tr_P = P0;  mode = CST_spfMode.NOMINAL;
    st = zeros(1, N);  det = [];  hb = [];  veto = [];  prob = [];  ret = [];  qRaw = NaN;  qBound = NaN;  passRaw = false;  passBound = false;
    qPos = NaN;  passPos = false;  tOnB = Inf;  xNE = zeros(N, 2);  wander = NaN(N, 1);  landing = NaN(N, 1);
    for k = 1:N
        z = sqrt(diag(R)) .* randn(m, 1);
        off = offsetAt(k);
        z(1:14) = z(1:14) + H(1:14, 1:3) * off(:);                                             % the ranges describe a point off; the baro is not spoofed
        if k >= tOnB && k <= tOffB, z(1:14) = z(1:14) + H(1:14, 1:3) * [stepB 0 0]'; end
        inProb = (mode == CST_spfMode.PROBATION);  kfUpd = (mode == CST_spfMode.NOMINAL);
        if inProb
            [tr_x, tr_P, y, S, ~, xp, xpP] = kalman_update_step(tr_x, tr_P, z, H, propTel, R);
            post = tr_x;  postP = tr_P;
        else
            [ax, aP, y, S, ~, xp, xpP] = kalman_update_step(kf_x, kf_P, z, H, propTel, R);
            post = ax;  postP = aP;
        end
        kfMeas = STRUCT_SPF.kfMeasFromUpdate(y, H, R, m, xp, xpP, post, postP);
        if k == 200 && mode == CST_spfMode.COAST                                             % mid-spoof: the two covariances side by side
            yP = zeros(mMax, 1);  yP(1:m) = y;  HP = zeros(mMax, n);  HP(1:m, :) = H;  RP = zeros(mMax);  RP(1:m, 1:m) = R;
            [passRaw, qRaw] = SPF_revalidation(yP, HP, RP, m, zeros(n, 1), sys.coastCov);
            [passBound, qBound] = SPF_revalidation(yP, HP, RP, m, zeros(n, 1), STRUCT_SPF.boundCoastCov(sys.coastCov, zeros(n, 1)));
            sc = ones(n, 1);                                                                  % the bound on the position alone
            for a = 1:3, if sys.coastCov(a, a) > 9, sc(a) = 3 / sqrt(sys.coastCov(a, a)); end, end
            [passPos, qPos] = SPF_revalidation(yP, HP, RP, m, zeros(n, 1), (sc * sc') .* sys.coastCov);
        end
        [sys, tel] = SPF_protectedNav(sys, kfMeas, propTel, k);
        mode = tel.info.mode;
        if kfUpd, kf_x = post;  kf_P = postP;
        else, [kf_x, kf_P] = SPF_insCoast(kf_x, kf_P, propTel);  kf_x(1) = kf_x(1) + DRIFT; end
        if tel.nav.applyCorrection, kf_x = tel.nav.state;  kf_P = tel.nav.covar; end
        if tel.kfCommand.startTrial, tr_x = kf_x;  tr_P = kf_P; end
        st(k) = double(tel.info.mode);  xNE(k, :) = kf_x(1:2)';
        if tel.info.returnChecked, wander(k) = tel.info.returnWander;  landing(k) = tel.info.returnLanding; end
        if tel.info.eventReturn, ret(end+1) = k; end %#ok<AGROW>
        if tel.info.eventLatched, det(end+1) = k; end %#ok<AGROW>
        if tel.info.eventHandback, hb(end+1) = k; end %#ok<AGROW>
        if tel.info.eventProbationVetoed, veto(end+1) = k; end %#ok<AGROW>
        if tel.info.eventProbationStarted
            prob(end+1) = k; %#ok<AGROW>
            if stepB ~= 0 && k > tOff && ~isfinite(tOnB), tOnB = k + 1; end                  % the spoofer resumes as GNSS is re-admitted
        end
    end
    res = struct('st', st, 'det', det, 'hb', hb, 'veto', veto, 'prob', prob, 'ret', ret, 'tOnB', tOnB, 'qRaw', qRaw, 'qBound', qBound, ...
        'passRaw', passRaw, 'passBound', passBound, 'qPos', qPos, 'passPos', passPos, 'xNE', xNE, 'wander', wander, 'landing', landing);
end
function s = pf(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end
function v = latchVal(x, det)
% the value at the first latch epoch, NaN without one
if isempty(det), v = NaN; else, v = x(det(1)); end
end
stepOn = @(s, k, t1, t2) s * double(k >= t1 && k <= t2);
rA = runCase(@(k) stepOn([30 0 0], k, 100, 300), 300, 0, 0, n, mMax, m, N, H, R, P0, propTel, DRIFT);
rB = runCase(@(k) stepOn([16 0 0], k, 100, 300), 300, 0, 0, n, mMax, m, N, H, R, P0, propTel, DRIFT);
rC = runCase(@(k) stepOn([16 0 0], k, 100, 200), 200, 30, 400, n, mMax, m, N, H, R, P0, propTel, DRIFT);
rD = runCase(@(k) stepOn([0 0 35], k, 100, 300), 300, 0, 0, n, mMax, m, N, H, R, P0, propTel, DRIFT);
rE = runCase(@(k) stepOn([0 0.15 * (k - 100) 0], k, 100, 300), 300, 0, 0, n, mMax, m, N, H, R, P0, propTel, DRIFT);
rG = runCase(@(k) stepOn([0 0.12 * (k - 100) 0], k, 100, 300) + stepOn([0 49 0], k, 301, 400), 400, 0, 0, n, mMax, m, N, H, R, P0, propTel, DRIFT);
thr = CST_spfParam.REVAL_THRESHOLD_TABLE(m);
fprintf('A 30 m: latch %s, hand-backs %s, probations %s; mid-spoof q on the raw coast covariance %.1f (%s), on the bounded %.1f (%s), threshold %.1f\n', ...
    mat2str(rA.det), mat2str(rA.hb), mat2str(rA.prob), rA.qRaw, pf(rA.passRaw), rA.qBound, pf(~rA.passBound), thr);
fprintf('B 16 m: latch %s, hand-backs %s, probations %s; mid-spoof chi-square on the bounded covariance %.1f (%s), threshold %.1f\n', ...
    mat2str(rB.det), mat2str(rB.hb), mat2str(rB.prob), rB.qBound, pf(rB.passBound), thr);
fprintf('C 16 m then 30 m from epoch %g: latch %s, probations %s, vetoes %s, hand-backs %s\n', rC.tOnB, mat2str(rC.det), mat2str(rC.prob), mat2str(rC.veto), mat2str(rC.hb));
fprintf('D 35 m down: latch %s, hand-backs %s, probations %s; mid-spoof q with the position bounded alone %.1f (%s), with the clock too %.1f (%s), threshold %.1f\n', ...
    mat2str(rD.det), mat2str(rD.hb), mat2str(rD.prob), rD.qPos, pf(rD.passPos), rD.qBound, pf(~rD.passBound), thr);
okA = ~isempty(rA.det) && rA.det(1) >= 100 && rA.det(1) <= 110 && ~any(rA.hb <= 300) && any(rA.hb > 300 & rA.hb < 360) ...
    && all(rA.st(110:300) == 2) && rA.passRaw && ~rA.passBound && rA.qBound > thr && rA.qRaw < thr;
okB = ~isempty(rB.det) && rB.det(1) >= 100 && rB.det(1) <= 110 && ~any(rB.hb <= 300) && any(rB.hb > 300 & rB.hb < 360) ...
    && ~any(rB.prob > 100 & rB.prob <= 300) && rB.passBound && all(rB.st(110:300) == 2);
okC = isfinite(rC.tOnB) && rC.tOnB > 200 && rC.tOnB < 260 && any(rC.veto == rC.tOnB) && ~any(rC.hb <= 400) ...
    && any(rC.hb > 400 & rC.hb < 460) && all(rC.st(rC.tOnB + 1:400) ~= 1);
okD = ~isempty(rD.det) && rD.det(1) >= 100 && rD.det(1) <= 112 && ~any(rD.hb <= 300) && any(rD.hb > 300 & rD.hb < 370) ...
    && ~any(rD.prob > 100 & rD.prob <= 300) && rD.passPos && ~rD.passBound && rD.qBound > thr && rD.qPos < thr;
fprintf('E ramp to 30 m E, followed: latch %s, return %s, hand-backs %s, nav E %.1f m at 299, %.1f m at the end; at the latch wander %.1f m, landing %.1f m\n', ...
    mat2str(rE.det), mat2str(rE.ret), mat2str(rE.hb), rE.xNE(299, 2), rE.xNE(end, 2), latchVal(rE.wander, rE.det), latchVal(rE.landing, rE.det));
fprintf('G drag to 24 m E then 49 m: latch %s, return %s, hand-backs %s, nav E %.1f m at 399, %.1f m at the end; at the latch landing %.1f m\n', ...
    mat2str(rG.det), mat2str(rG.ret), mat2str(rG.hb), rG.xNE(399, 2), rG.xNE(end, 2), latchVal(rG.landing, rG.det));
okE = ~isempty(rE.det) && rE.det(1) >= 301 && rE.det(1) <= 312 && numel(rE.ret) == 1 && rE.ret(1) == rE.det(1) ...
    && ~any(rE.hb <= 300) && any(rE.hb > 300 & rE.hb < 350) && all(rE.st(1:300) == 1) && abs(rE.xNE(299, 2) - 30) < 2 ...
    && max(abs(rE.xNE(end, :))) < 3 && latchVal(rE.wander, rE.det) > 15 && latchVal(rE.landing, rE.det) < 15;
okG = ~isempty(rG.det) && rG.det(1) >= 301 && rG.det(1) <= 312 && ~any(rG.ret <= 400) && numel(rG.ret) == 1 && rG.ret(1) >= 401 && rG.ret(1) <= 403 ...
    && ~any(rG.hb <= 400) && any(rG.hb > 400 & rG.hb < 445) && all(rG.st(313:400) ~= 1) && abs(rG.xNE(399, 2) - 24) < 4 ...
    && max(abs(rG.xNE(end, :))) < 3 && latchVal(rG.landing, rG.det) > 15 && all(rG.st(1:300) == 1);
fprintf('bound: %s | A, the bound holds 30 m where the raw covariance admits it: %s | B, the CPI against the coast holds 16 m: %s | C, the first window vetoes a spoof resuming at the probation: %s | D, the clock bound holds a down step the position bound alone admits: %s\n', ...
    pf(okBound), pf(okA), pf(okB), pf(okC), pf(okD));
fprintf('E, the end of a followed ramp reads as a return and is handed back: %s | G, a drag then a step is held, the return read in COAST when the GNSS comes back: %s\n', pf(okE), pf(okG));
if okBound && okA && okB && okC && okD && okE && okG, fprintf('=== test_coastBound PASS ===\n'); else, fprintf('=== test_coastBound FAIL ===\n'); end
