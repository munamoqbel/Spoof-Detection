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
%   E  the return path: a 1 m per epoch east ramp from 100 to 300 (200 m) is followed without a latch; the step back
%      at 301 latches and reads as a return (the nav had wandered 200 m from the long coast, the GNSS lands on it), the
%      nav is corrected back and handed back within 40 epochs, on the truth; without the return path the clean GNSS
%      would be refused for as long as it lasts. The ramp is sized to the yardstick of the return check under this
%      host's Q: the long coast's sigma at its age (160 to 250 epochs of 2 m^2) is near 20 m per axis, so
%      K_MISSED_DETECTION * sigma_long is near 140 m, and only a followed spoof beyond that reads as a return. The
%      landing reads some 40 m, not the few metres of the geometry: the position the rows ask for is the CPI projection,
%      one axis at a time with the full S, and through these rows' lines of sight a 200 m east offset leaks about a
%      fifth of itself onto north (the leakage scales with the offset asked, so it stays inside a yardstick the wander
%      has to exceed anyway);
%   G  a slow drag of 200 m (1 m per epoch, 100 to 300) and then a step 50 m further (301 to 400): the step latches as
%      an onset (the GNSS lands 250 m from the long coast, beyond the yardstick) and is held; at 401 the clean GNSS comes
%      back 200 m from the nav, which the strict bound would refuse, and the return check in COAST reads it (the nav had
%      wandered 200 m from the long coast, beyond the yardstick, the GNSS lands within a few metres of it): corrected
%      back and handed back;
%   F  the 30 m ramp the return path was built on (0.15 m per epoch, 100 to 300), followed: the yardstick sits far above
%      the 30 m wander, the step back is not read as a return, and the clean GNSS stays refused through the strict bound
%      to the end of the record: the reach of the return path is the INS's own drift uncertainty at the long coast's age;
%   H  the false return measured on the vehicle: an INS that drifts 0.1 m per epoch north, corrected by the GNSS in
%      NOMINAL so that the long coast sits 16 to 24 m north of the nav at its age, and a 20 m north step at 300: the
%      step lands within 4 m of the long coast while the nav has wandered 16 m or more from it, which fixed 15 m radii
%      read as a return (the nav moved onto the drifted coast, the re-validation loosened to that coast's uncertainty,
%      the spoof handed back); against the yardstick it is an onset: no return, no hand-back while the spoof lasts;
%   I  the F ramp with the coast fast-forwarded by 1100 epochs at the latch (the hold counter, the coast covariance by
%      1100 Q, the INS drift: what 1100 more coast epochs would have left, Phi being the identity here), so that the hold
%      limit, COAST_BUDGET_EPOCHS = 1200 epochs, falls 100 epochs after the latch: the PL budget flag is up from the
%      latch on and changes nothing; at the limit the re-admission runs on the coast's own covariance, the paper's
%      recovery, the clean GNSS passes, the probation's windows open on that covariance so the trial's move is not
%      vetoed, and the hand-back follows, flagged handbackUncertified; the nav is back on the truth by the end;
%   J  a 30 m east step from 100 to 500, fast-forwarded the same way: held on the bound through the PL flag until the
%      limit, then admitted the same way, the hand-back into the spoof flagged; the step back at 500 latches again as an
%      onset and is held to the end (its own limit lies beyond the record): the hold on any spoof is bounded by the limit
%      alone (east, across the synthetic INS drift, which runs north and would otherwise carry the coast onto the spoof).
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
function [res] = runCase(offsetAt, tOff, stepB, tOffB, n, mMax, m, N, H, R, P0, propTel, DRIFT, driftAlways, fastForward)
    % offsetAt(k) = the [N E D] offset the ranges describe at epoch k; stepB (north) from the epoch after the first probation
    % that opens after tOff until tOffB; driftAlways = the INS drifts in NOMINAL too (the GNSS corrects it there, and the
    % long coast collects that drift); fastForward = coast epochs added at the first latch (the hold counter, the coast
    % covariance by that many Q, the INS drift), 0 for none
    sys = STRUCT_SPF.zeroSys;  sys.coastCov = P0;
    sys.anchor = STRUCT_SPF.setAnchor(true, zeros(n, 1), P0, uint32(0));
    kf_x = zeros(n, 1);  kf_P = P0;  tr_x = kf_x;  tr_P = P0;  mode = CST_spfMode.NOMINAL;
    st = zeros(1, N);  det = [];  hb = [];  veto = [];  prob = [];  ret = [];  qRaw = NaN;  qBound = NaN;  passRaw = false;  passBound = false;
    qPos = NaN;  passPos = false;  tOnB = Inf;  xNE = zeros(N, 2);  wander = NaN(N, 1);  landing = NaN(N, 1);  yard = NaN(N, 1);
    budget = false(N, 1);  hbUnc = false(N, 1);  forwarded = false;
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
        if kfUpd, kf_x = post;  kf_P = postP;  if driftAlways, kf_x(1) = kf_x(1) + DRIFT; end
        else, [kf_x, kf_P] = SPF_insCoast(kf_x, kf_P, propTel);  kf_x(1) = kf_x(1) + DRIFT; end
        if tel.nav.applyCorrection, kf_x = tel.nav.state;  kf_P = tel.nav.covar; end
        if fastForward > 0 && tel.info.eventLatched && ~forwarded
            sys.coastCount = sys.coastCount + fastForward;                                     % the hold counter
            sys.coastCov   = sys.coastCov + fastForward * propTel.accumQ;                       % the coast covariance, Phi = I
            kf_x(1) = kf_x(1) + fastForward * DRIFT;  kf_P = kf_P + fastForward * propTel.accumQ;   % the INS as it would have drifted
            forwarded = true;
        end
        if tel.kfCommand.startTrial, tr_x = kf_x;  tr_P = kf_P; end
        st(k) = double(tel.info.mode);  xNE(k, :) = kf_x(1:2)';
        budget(k) = tel.info.coastBudgetExceeded;  hbUnc(k) = tel.info.handbackUncertified;
        if tel.info.returnChecked, wander(k) = tel.info.returnWander;  landing(k) = tel.info.returnLanding;  yard(k) = tel.info.returnYardstick; end
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
        'passRaw', passRaw, 'passBound', passBound, 'qPos', qPos, 'passPos', passPos, 'xNE', xNE, 'wander', wander, 'landing', landing, 'yard', yard, ...
        'budget', budget, 'hbUnc', hbUnc);
end
function s = pf(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end
function v = latchVal(x, det)
% the value at the first latch epoch, NaN without one
if isempty(det), v = NaN; else, v = x(det(1)); end
end
stepOn = @(s, k, t1, t2) s * double(k >= t1 && k <= t2);
rA = runCase(@(k) stepOn([30 0 0], k, 100, 300), 300, 0, 0, n, mMax, m, N, H, R, P0, propTel, DRIFT, false, 0);
rB = runCase(@(k) stepOn([16 0 0], k, 100, 300), 300, 0, 0, n, mMax, m, N, H, R, P0, propTel, DRIFT, false, 0);
rC = runCase(@(k) stepOn([16 0 0], k, 100, 200), 200, 30, 400, n, mMax, m, N, H, R, P0, propTel, DRIFT, false, 0);
rD = runCase(@(k) stepOn([0 0 35], k, 100, 300), 300, 0, 0, n, mMax, m, N, H, R, P0, propTel, DRIFT, false, 0);
rE = runCase(@(k) stepOn([0 1.0 * (k - 100) 0], k, 100, 300), 300, 0, 0, n, mMax, m, N, H, R, P0, propTel, DRIFT, false, 0);
rG = runCase(@(k) stepOn([0 1.0 * (k - 100) 0], k, 100, 300) + stepOn([0 250 0], k, 301, 400), 400, 0, 0, n, mMax, m, N, H, R, P0, propTel, DRIFT, false, 0);
rF = runCase(@(k) stepOn([0 0.15 * (k - 100) 0], k, 100, 300), 300, 0, 0, n, mMax, m, N, H, R, P0, propTel, DRIFT, false, 0);
rH = runCase(@(k) stepOn([20 0 0], k, 300, 400), 400, 0, 0, n, mMax, m, N, H, R, P0, propTel, 0.1, true, 0);
FF = 1100;                                                                                % coast epochs fast-forwarded at the latch in I and J
rI = runCase(@(k) stepOn([0 0.15 * (k - 100) 0], k, 100, 300), 300, 0, 0, n, mMax, m, 600, H, R, P0, propTel, DRIFT, false, FF);
rJ = runCase(@(k) stepOn([0 30 0], k, 100, 500), 500, 0, 0, n, mMax, m, 600, H, R, P0, propTel, DRIFT, false, FF);
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
    && any(rC.hb > 400 & rC.hb < 460) && all(rC.st(rC.tOnB + 1:400) ~= 1) && ~any(rC.hbUnc);
okD = ~isempty(rD.det) && rD.det(1) >= 100 && rD.det(1) <= 112 && ~any(rD.hb <= 300) && any(rD.hb > 300 & rD.hb < 370) ...
    && ~any(rD.prob > 100 & rD.prob <= 300) && rD.passPos && ~rD.passBound && rD.qBound > thr && rD.qPos < thr;
fprintf('E ramp to 200 m E, followed: latch %s, return %s, hand-backs %s, nav E %.1f m at 299, %.1f m at the end; at the latch wander %.1f m, landing %.1f m, yardstick %.1f m\n', ...
    mat2str(rE.det), mat2str(rE.ret), mat2str(rE.hb), rE.xNE(299, 2), rE.xNE(end, 2), latchVal(rE.wander, rE.det), latchVal(rE.landing, rE.det), latchVal(rE.yard, rE.det));
fprintf('G drag to 200 m E then 250 m: latch %s, return %s, hand-backs %s, nav E %.1f m at 399, %.1f m at the end; at the latch landing %.1f m, yardstick %.1f m\n', ...
    mat2str(rG.det), mat2str(rG.ret), mat2str(rG.hb), rG.xNE(399, 2), rG.xNE(end, 2), latchVal(rG.landing, rG.det), latchVal(rG.yard, rG.det));
okE = ~isempty(rE.det) && rE.det(1) >= 301 && rE.det(1) <= 312 && numel(rE.ret) == 1 && rE.ret(1) == rE.det(1) ...
    && ~any(rE.hb <= 300) && any(rE.hb > 300 & rE.hb < 350) && all(rE.st(1:300) == 1) && abs(rE.xNE(299, 2) - 200) < 5 ...
    && max(abs(rE.xNE(end, :))) < 3 && latchVal(rE.wander, rE.det) > latchVal(rE.yard, rE.det) && latchVal(rE.yard, rE.det) > 100 ...
    && latchVal(rE.landing, rE.det) < 0.3 * latchVal(rE.wander, rE.det);                   % the projection's cross-axis leakage, see the header
okG = ~isempty(rG.det) && rG.det(1) >= 301 && rG.det(1) <= 312 && ~any(rG.ret <= 400) && numel(rG.ret) == 1 && rG.ret(1) >= 401 && rG.ret(1) <= 403 ...
    && ~any(rG.hb <= 400) && any(rG.hb > 400 & rG.hb < 445) && all(rG.st(313:400) ~= 1) && rG.xNE(399, 2) > 150 && rG.xNE(399, 2) < 205 ...
    && max(abs(rG.xNE(end, :))) < 3 && latchVal(rG.landing, rG.det) > latchVal(rG.yard, rG.det) && all(rG.st(1:300) == 1) ...
    && rG.yard(rG.ret(1)) < rG.wander(rG.ret(1)) && rG.landing(rG.ret(1)) < 0.3 * rG.wander(rG.ret(1));
fprintf('F ramp to 30 m E, followed: latch %s, return %s, hand-backs %s, nav E %.1f m at the end; at the latch wander %.1f m, yardstick %.1f m\n', ...
    mat2str(rF.det), mat2str(rF.ret), mat2str(rF.hb), rF.xNE(end, 2), latchVal(rF.wander, rF.det), latchVal(rF.yard, rF.det));
fprintf('H 20 m N onset on an INS drifting 0.1 m per epoch: latch %s, return %s, hand-backs %s, nav N %.1f m at 399; at the latch wander %.1f m, landing %.1f m, yardstick %.1f m\n', ...
    mat2str(rH.det), mat2str(rH.ret), mat2str(rH.hb), rH.xNE(399, 1), latchVal(rH.wander, rH.det), latchVal(rH.landing, rH.det), latchVal(rH.yard, rH.det));
okF = ~isempty(rF.det) && rF.det(1) >= 301 && rF.det(1) <= 312 && isempty(rF.ret) && isempty(rF.hb) && all(rF.st(320:end) == 2) ...
    && latchVal(rF.wander, rF.det) > 15 && latchVal(rF.yard, rF.det) > latchVal(rF.wander, rF.det);
okH = ~isempty(rH.det) && rH.det(1) >= 300 && rH.det(1) <= 312 && isempty(rH.ret) && ~any(rH.hb <= 400) && all(rH.st(313:400) == 2) ...
    && latchVal(rH.wander, rH.det) >= CST_spfParam.RETURN_WANDER && latchVal(rH.landing, rH.det) <= CST_spfParam.RETURN_RADIUS ...
    && latchVal(rH.yard, rH.det) > latchVal(rH.wander, rH.det) && rH.xNE(399, 1) < 15;
fprintf('I the F ramp, the coast fast-forwarded %d epochs at the latch: latch %s, PL flag from %s, probations %s, hand-backs %s (flagged %s), nav E %.1f m at the end\n', ...
    FF, mat2str(rI.det), mat2str(find(rI.budget, 1)), mat2str(rI.prob), mat2str(rI.hb), mat2str(find(rI.hbUnc)'), rI.xNE(end, 2));
fprintf('J 30 m E from 100 to 500, fast-forwarded the same: latches %s, PL flag from %s, hand-backs %s (flagged %s), nav E %.1f m at 490, %.1f m at the end\n', ...
    mat2str(rJ.det), mat2str(find(rJ.budget, 1)), mat2str(rJ.hb), mat2str(find(rJ.hbUnc)'), rJ.xNE(490, 2), rJ.xNE(end, 2));
lim = double(CST_spfParam.COAST_BUDGET_EPOCHS) - FF;                                      % the limit falls this many epochs after the latch
okI = ~isempty(rI.det) && rI.det(1) >= 301 && rI.det(1) <= 312 && isempty(rI.ret) && ~any(rI.hb < rI.det(1) + lim) ...
    && any(rI.hb >= rI.det(1) + lim & rI.hb < rI.det(1) + lim + 80) && ~isempty(rI.prob) && rI.hbUnc(rI.hb(end)) ...
    && max(abs(rI.xNE(end, :))) < 3 && rI.st(end) == 1 && all(rI.st(320:rI.det(1) + lim) == 2) && any(rI.budget(320:400));
okJ = ~isempty(rJ.det) && rJ.det(1) >= 100 && rJ.det(1) <= 110 && all(rJ.st(110:rJ.det(1) + lim) == 2) ...
    && any(rJ.hb > rJ.det(1) + lim & rJ.hb < 500) && rJ.hbUnc(rJ.hb(1)) && abs(rJ.xNE(490, 2) - 30) < 4 ...
    && any(rJ.det > 500 & rJ.det < 512) && ~any(rJ.hb > 500) && rJ.st(end) == 2 && any(rJ.budget(110:200));
fprintf('bound: %s | A, the bound holds 30 m where the raw covariance admits it: %s | B, the CPI against the coast holds 16 m: %s | C, the first window vetoes a spoof resuming at the probation: %s | D, the clock bound holds a down step the position bound alone admits: %s\n', ...
    pf(okBound), pf(okA), pf(okB), pf(okC), pf(okD));
fprintf('E, the end of a followed ramp reads as a return and is handed back: %s | G, a drag then a step is held, the return read in COAST when the GNSS comes back: %s\n', pf(okE), pf(okG));
fprintf('F, a 30 m followed ramp ends below the yardstick: no return, the clean GNSS held: %s | H, the false return on the INS''s own drift, which fixed radii read, is refused by the yardstick: %s\n', pf(okF), pf(okH));
fprintf('I, past the hold limit the clean GNSS is re-admitted on the coast''s own covariance, flagged: %s | J, a step is held through the PL flag to the limit, then admitted flagged, and its step back held: %s\n', pf(okI), pf(okJ));
if okBound && okA && okB && okC && okD && okE && okG && okF && okH && okI && okJ, fprintf('=== test_coastBound PASS ===\n'); else, fprintf('=== test_coastBound FAIL ===\n'); end
