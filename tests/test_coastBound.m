%% test_coastBound.m
% The re-admission path of the gate: the coast covariance bounded by COAST_SIGMA_MAX, the CPI of the innovations
% against the coast in the dwell, and the first probation window opened on the coast. (1) STRUCT_SPF.boundCoastCov: a
% matrix whose monitored position sigmas exceed the bound comes back with them at the bound, symmetric positive
% semi-definite, the other states and the cross terms scaled by congruence; one within the bound comes back as it is.
% (2) A synthetic 60-state host (static truth, 14 range rows at 5 m, a position process noise that grows the coast to
% 6 m sigma within 20 epochs, an INS that drifts 0.02 m per coast epoch):
%   A  a 30 m step on every row from epoch 100 to 300 latches, stays refused through the spoof and is handed back after
%      it, while SPF_revalidation on the unbounded coast covariance would have passed mid-spoof;
%   B  a 16 m step, which the chi-square passes even on the bounded covariance, is held by the CPI against the coast:
%      no probation while it lasts, hand-back after it;
%   C  a 16 m step from 100 to 200 and then a 30 m spoof switched on at the epoch the probation opens: the trial jumps
%      onto it in its first update, the first window vetoes that epoch, the spoof stays refused until it ends at 400.
% Run from the repo root (MATLAB, or Octave with tools/octave_shim).

n = CST_gnssHybrid.NO_STATES;  mMax = double(CST_spfParam.MAX_MEAS);
CAP = CST_spfParam.COAST_SIGMA_MAX;

%% 1. the bound itself
P = eye(n);  P(1, 1) = 36;  P(2, 2) = 16;  P(3, 3) = 4;  P(4, 4) = 0.25;
P(1, 4) = 2;  P(4, 1) = 2;  P(1, 2) = 10;  P(2, 1) = 10;                                % cross terms N-vN and N-E
Pb = STRUCT_SPF.boundCoastCov(P);
okBound = abs(Pb(1, 1) - CAP^2) < 1e-12 && abs(Pb(2, 2) - CAP^2) < 1e-12 && Pb(3, 3) == 4 && Pb(4, 4) == 0.25 ...
    && abs(Pb(1, 4) - 2 * CAP / 6) < 1e-12 && abs(Pb(1, 2) - 10 * (CAP / 6) * (CAP / 4)) < 1e-12 ...
    && isequal(Pb, Pb') && min(eig(Pb(1:4, 1:4))) > -1e-12 && isequal(STRUCT_SPF.boundCoastCov(0.5 * eye(n)), 0.5 * eye(n));

%% 2. the synthetic host
rng(3);
N = 500;  m = 14;  R = 25 * eye(m);
Phi = eye(n);  Q = zeros(n);  Q(1:3, 1:3) = 2.0 * eye(3);                                  % m^2 per epoch on the position
propTel = STRUCT_SPF.setPropTel(Phi, Q);
el = 0.3 + 0.6 * rand(m, 1);  az = 2 * pi * rand(m, 1);
u = [cos(el) .* cos(az), cos(el) .* sin(az), -sin(el)];
H = zeros(m, n);  H(:, 1:3) = -u;
P0 = 0.01 * eye(n);  P0(1:3, 1:3) = eye(3);
DRIFT = 0.02;                                                                             % m per coast epoch, north
function [res] = runCase(stepA, tOn, tOff, stepB, tOffB, n, mMax, m, N, H, R, P0, propTel, DRIFT)
    % stepA on every row from tOn to tOff; stepB from the epoch after the first probation that opens after tOff until tOffB
    sys = STRUCT_SPF.zeroSys;  sys.coastCov = P0;
    sys.anchor = STRUCT_SPF.setAnchor(true, zeros(n, 1), P0, uint32(0));
    kf_x = zeros(n, 1);  kf_P = P0;  tr_x = kf_x;  tr_P = P0;  mode = CST_spfMode.NOMINAL;
    st = zeros(1, N);  det = [];  hb = [];  veto = [];  prob = [];  qRaw = NaN;  qBound = NaN;  passRaw = false;  passBound = false;
    tOnB = Inf;
    for k = 1:N
        z = 5 * randn(m, 1);
        if k >= tOn && k <= tOff, z = z + H(:, 1:3) * [stepA 0 0]'; end                     % the ranges describe a point stepA north
        if k >= tOnB && k <= tOffB, z = z + H(:, 1:3) * [stepB 0 0]'; end
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
            [passBound, qBound] = SPF_revalidation(yP, HP, RP, m, zeros(n, 1), STRUCT_SPF.boundCoastCov(sys.coastCov));
        end
        [sys, tel] = SPF_protectedNav(sys, kfMeas, propTel, k);
        mode = tel.info.mode;
        if kfUpd, kf_x = post;  kf_P = postP;
        else, [kf_x, kf_P] = SPF_insCoast(kf_x, kf_P, propTel);  kf_x(1) = kf_x(1) + DRIFT; end
        if tel.nav.applyCorrection, kf_x = tel.nav.state;  kf_P = tel.nav.covar; end
        if tel.kfCommand.startTrial, tr_x = kf_x;  tr_P = kf_P; end
        st(k) = double(tel.info.mode);
        if tel.info.eventLatched, det(end+1) = k; end %#ok<AGROW>
        if tel.info.eventHandback, hb(end+1) = k; end %#ok<AGROW>
        if tel.info.eventProbationVetoed, veto(end+1) = k; end %#ok<AGROW>
        if tel.info.eventProbationStarted
            prob(end+1) = k; %#ok<AGROW>
            if stepB ~= 0 && k > tOff && ~isfinite(tOnB), tOnB = k + 1; end                  % the spoofer resumes as GNSS is re-admitted
        end
    end
    res = struct('st', st, 'det', det, 'hb', hb, 'veto', veto, 'prob', prob, 'tOnB', tOnB, 'qRaw', qRaw, 'qBound', qBound, ...
        'passRaw', passRaw, 'passBound', passBound);
end
rA = runCase(30, 100, 300, 0, 0, n, mMax, m, N, H, R, P0, propTel, DRIFT);
rB = runCase(16, 100, 300, 0, 0, n, mMax, m, N, H, R, P0, propTel, DRIFT);
rC = runCase(16, 100, 200, 30, 400, n, mMax, m, N, H, R, P0, propTel, DRIFT);
thr = CST_spfParam.REVAL_THRESHOLD_TABLE(m);
fprintf('A 30 m: latch %s, hand-backs %s, probations %s; mid-spoof q on the raw coast covariance %.1f (%s), on the bounded %.1f (%s), threshold %.1f\n', ...
    mat2str(rA.det), mat2str(rA.hb), mat2str(rA.prob), rA.qRaw, pf(rA.passRaw), rA.qBound, pf(~rA.passBound), thr);
fprintf('B 16 m: latch %s, hand-backs %s, probations %s; mid-spoof chi-square on the bounded covariance %.1f (%s), threshold %.1f\n', ...
    mat2str(rB.det), mat2str(rB.hb), mat2str(rB.prob), rB.qBound, pf(rB.passBound), thr);
fprintf('C 16 m then 30 m from epoch %g: latch %s, probations %s, vetoes %s, hand-backs %s\n', rC.tOnB, mat2str(rC.det), mat2str(rC.prob), mat2str(rC.veto), mat2str(rC.hb));
okA = ~isempty(rA.det) && rA.det(1) >= 100 && rA.det(1) <= 110 && ~any(rA.hb <= 300) && any(rA.hb > 300 & rA.hb < 360) ...
    && all(rA.st(110:300) == 2) && rA.passRaw && ~rA.passBound && rA.qBound > thr && rA.qRaw < thr;
okB = ~isempty(rB.det) && rB.det(1) >= 100 && rB.det(1) <= 110 && ~any(rB.hb <= 300) && any(rB.hb > 300 & rB.hb < 360) ...
    && ~any(rB.prob > 100 & rB.prob <= 300) && rB.passBound && all(rB.st(110:300) == 2);
okC = isfinite(rC.tOnB) && rC.tOnB > 200 && rC.tOnB < 260 && any(rC.veto == rC.tOnB) && ~any(rC.hb <= 400) ...
    && any(rC.hb > 400 & rC.hb < 460) && all(rC.st(rC.tOnB + 1:400) ~= 1);
fprintf('bound: %s | A, the bound holds 30 m where the raw covariance admits it: %s | B, the CPI against the coast holds 16 m: %s | C, the first window vetoes a spoof resuming at the probation: %s\n', ...
    pf(okBound), pf(okA), pf(okB), pf(okC));
if okBound && okA && okB && okC, fprintf('=== test_coastBound PASS ===\n'); else, fprintf('=== test_coastBound FAIL ===\n'); end

function s = pf(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end
