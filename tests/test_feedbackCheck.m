function test_feedbackCheck()
%% test_feedbackCheck.m
% tools/feedbackCheck on a synthetic 100 Hz nav log with 2 Hz updates. The mechanised solution integrates its
% velocity and receives the KF correction as a jump on the first sample at or after tKf: in full (A: D velocity
% and Z bias not applied; B: everything applied; C: nothing applied), with the records offset by one epoch
% from the nav log (D), slewed in with a 0.1 s time constant (E), through the host's split feedback with
% GAIN 0.4 applied in full and the residual propagated with Phi (F), and through a split that books 0.4 but
% applies 0.16 with the residual propagated with Phi transposed (G). The lag coefficients, the gain on the
% state, the step times, the bookkeeping shares and the verdicts must read each case. Run from anywhere:
% octave --no-gui --quiet --eval test_feedbackCheck
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'tools'));
randn('seed', 3);

dt = 0.01;  T = 400;  t = (0:dt:T)';  N = numel(t);
tKf = (1:0.5:T-1)';  K = numel(tKf);
accTrue = [0.3 * sin(2 * pi * t / 60), 0.2 * cos(2 * pi * t / 45), 0.05 * sin(2 * pi * t / 20)];
velTrue = cumsum(accTrue) * dt + repmat([20 5 0], N, 1);
near = @(x, v, tol) all(abs(x(:) - v) < tol);
cnt  = @(txt, s) numel(strfind(txt, s));

ok = true;
chk = struct();
% A: full feedback, D velocity and Z bias not applied
[runA, logA] = makeRun(t, tKf, velTrue, [1 1 0], 1, 0, 0, 1, 1, 'none', 0, 0.02, 1e-5);
txtA = evalc('resA = feedbackCheck(runA, logA);');
fprintf('%s', txtA);
chk.posA      = near(resA.coefPos(2, :), -1, 0.03) && near(resA.coefPos([1 3], :), 0, 0.03) && near(resA.gainPos(1, :), -1, 0.03);
chk.velA      = near(resA.coefVel(2, 1:2), -1, 0.03) && abs(resA.coefVel(2, 3)) < 0.05 && near(resA.coefVel([1 3], 1:2), 0, 0.03);
chk.stepA     = near(resA.t50(1, :), 0, 0.011) && near(resA.t90(1, :) - resA.t10(1, :), 0, 0.011) && abs(resA.tStep) < 0.011;
chk.keepA     = near([resA.keep.position, resA.keep.velocity], 0, 0.01);
chk.usedA     = resA.nUsed == K;
chk.biasA     = abs(resA.dBias(3) / resA.incBias(3)) < 0.05 && near(resA.dBias(1:2) ./ resA.incBias(1:2), -1, 1e-6);
chk.errA      = abs(resA.slopeE(3) / 9.80665e-3 - 2.04) < 0.05 && all(abs(resA.slopeE(1:2)) < 1e-4) && abs(resA.jumpE(3)) < 1e-3 && ...
                abs(resA.jumpE(1) + mean(resA.incVel(resA.used, 1))) < 1e-3;
chk.blocksA   = numel(resA.blocks.t0) == 4 && all(abs(resA.blocks.acc(:, 3) - 2.04) < 0.05) && ~isfinite(resA.blocks.slopeOnBias) && ...
                abs(resA.blocks.trendErr) < 0.1 && cnt(txtA, 'does not act on the vertical acceleration') == 0 && cnt(txtA, 'runs away from the truth') == 0;
chk.waveA     = all(abs(resA.slopeV(1:2)) < 0.005) && abs(resA.jumpV(3)) < 0.002 && ...
                abs(resA.jumpV(1) + mean(resA.incVel(resA.used, 1))) < 0.002 && abs(resA.meanWaveP(end, 3) + mean(resA.incPos(resA.used, 3))) < 0.003;
chk.verdictA  = cnt(txtA, 'the D velocity increment does not reach') == 1 && cnt(txtA, 'bias increments reach biasEst on some axes only') == 1 && cnt(txtA, 'VERDICT') == 2;
% B: everything applied, cell input, biasEst per nav sample, a window
[runB, logB] = makeRun(t, tKf, velTrue, [1 1 1], 1, 1, 0, 1, 1, 'none');
runB.biasEst = interp1(tKf, runB.biasEst, t, 'previous', 'extrap');
txtB = evalc('resB = feedbackCheck(runB, num2cell(logB), [100 300]);');
chk.coefB     = near(resB.coefPos(2, :), -1, 0.03) && near(resB.coefVel(2, :), -1, 0.03) && near(resB.gainVel, -1, 0.03);
chk.waveB     = all(abs(resB.slopeV(1:2)) < 0.003) && all(abs(resB.slopeE) < 1e-4) && all(abs(resB.meanWaveP(end, :) + mean(resB.incPos(resB.used, :), 1)) < 0.004);
chk.winB      = resB.nUsed == nnz(tKf >= 100 & tKf <= 300) && all(resB.t(resB.used) >= 100) && all(resB.t(resB.used) <= 300);
chk.biasB     = near(resB.dBias ./ resB.incBias, -1, 1e-6);
chk.quietB    = cnt(txtB, 'VERDICT') == 0;
% C: nothing applied
[runC, logC] = makeRun(t, tKf, velTrue, [0 0 0], 0, 0, 0, 1, 1, 'none');
txtC = evalc('resC = feedbackCheck(runC, logC);');
chk.coefC     = near(resC.coefPos, 0, 0.05) && near(resC.coefVel, 0, 0.05);
chk.verdictC  = cnt(txtC, 'does not reach the logged solution') == 6;
% D: the records one epoch ahead of the nav log (record k holds the increment applied at tKf(k+1))
[runD, logD] = makeRun(t, tKf, velTrue, [1 1 1], 1, 1, 0, 1, 1, 'none');
logD = [logD(2:end); logD(1)];  logD(end).postState(:) = 0;
txtD = evalc('resD = feedbackCheck(runD, logD);');
chk.coefD     = near(resD.coefPos(1, :), -1, 0.03) && near(resD.coefPos([2 3], :), 0, 0.03) && near(resD.coefVel(1, :), -1, 0.03);
chk.verdictD  = cnt(txtD, 'one epoch ahead of the nav log') == 6;
% E: the correction slewed in with a 0.1 s time constant
[runE, logE] = makeRun(t, tKf, velTrue, [1 1 1], 1, 1, 0.1, 1, 1, 'none');
txtE = evalc('resE = feedbackCheck(runE, logE);');
chk.totalE    = near(sum(resE.coefPos, 1), -1, 0.05) && near(resE.coefPos(2, :), -0.92, 0.05);
chk.stepE     = all(resE.t90(1, :) - resE.t10(1, :) > 0.15) && near(resE.t10(1, :), 0.01, 0.015);
chk.verdictE  = cnt(txtE, 'slewed feedback') == 3 && cnt(txtE, 'position step is spread') == 3;
% F: the host's split feedback, GAIN 0.4 booked and applied, residual propagated with Phi
[runF, logF] = makeRun(t, tKf, velTrue, [1 1 1], 1, 1, 0, 0.4, 0.4, 'right');
txtF = evalc('resF = feedbackCheck(runF, logF);');
fprintf('%s', txtF);
chk.lagF      = near(resF.coefPos(2, 1:2), -0.40, 0.03) && near(resF.coefPos(1, 1:2), -0.24, 0.04) && near(resF.coefPos(3, 1:2), 0, 0.03) && near(resF.coefPos(2, 3), -0.40, 0.1);
chk.gainF     = near(resF.gainPos(1, :), -0.40, 0.02) && all(resF.rxPos < -0.98) && near(resF.gainVel, -0.40, 0.06);
chk.keepF     = near([resF.keep.position, resF.keep.velocity], 0.6, 0.01) && near(resF.aPV, 0.30, 0.02) && near(resF.aVP, 0, 0.02);
chk.stepF     = near(resF.t50(1, :), 0, 0.011) && near(resF.t90(1, :) - resF.t10(1, :), 0, 0.011);
chk.verdictF  = cnt(txtF, 'VERDICT') == 0 && cnt(txtF, 'NOTE: ') == 6 && cnt(txtF, 'consistent') == 6;
% G: the split books 0.4 and applies 0.16, residual propagated with Phi transposed
[runG, logG] = makeRun(t, tKf, velTrue, [1 1 1], 1, 1, 0, 0.4, 0.16, 'transposed');
txtG = evalc('resG = feedbackCheck(runG, logG);');
chk.gainG     = near(resG.gainPos(1, :), -0.16, 0.02) && all(resG.rxPos < -0.98) && near(resG.coefPos(2, :), -0.16, 0.03);
chk.keepG     = near([resG.keep.position, resG.keep.velocity], 0.6, 0.01) && near(resG.aVP, 0.30, 0.02) && near(resG.aPV, 0, 0.02);
chk.verdictG  = cnt(txtG, 'lost between setKF and the mechanisation') == 6 && cnt(txtG, 'propagated with Phi transposed') == 1 && cnt(txtG, 'does not reach') == 0 && cnt(txtG, 'without the Phi coupling') == 0;
% H: the velocity correction applied at the step and taken back 0.25 s later
[runH, logH] = makeRun(t, tKf, velTrue, [1 1 1], 1, 1, 0, 1, 1, 'none', 0.25);
txtH = evalc('resH = feedbackCheck(runH, logH);');
chk.gainH     = near(resH.gainVel, -1, 0.1) && near(resH.gainVelLate, 0, 0.1) && near(resH.gainPos(1, :), -1, 0.03);
chk.verdictH  = cnt(txtH, 'something restores the velocity between the updates') == 3;
% I: a constant 2 mg vertical acceleration error, everything applied, the Z bias estimate drifting away from a constant truth
[runI, logI] = makeRun(t, tKf, velTrue, [1 1 1], 1, 1, 0, 1, 1, 'none', 0, 0.02, 1e-5);
txtI = evalc('resI = feedbackCheck(runI, logI);');
chk.blocksI   = all(abs(resI.blocks.acc(:, 3) - 2.04) < 0.05) && abs(resI.blocks.slopeOnBias) < 0.3 && resI.blocks.seOnBias < 0.3 && ...
                abs(resI.blocks.trendErr + 2.04) < 0.4 && abs(resI.slopeE(3) / 9.80665e-3 - 2.04) < 0.05;
chk.verdictI  = cnt(txtI, 'does not act on the vertical acceleration') == 1 && cnt(txtI, 'runs away from the truth') == 1 && cnt(txtI, 'VERDICT') == 2;
% the coupling block from the postCov of the records
chk.covA      = abs(resA.corrVelBias(3, 3) - 0.5) < 1e-9 && abs(resA.corrVelBias(1, 1) - 0.3) < 1e-9 && abs(resA.corrPosVel(1) - 0.7) < 1e-9 && ...
                resA.corrVelBias(2, 2) == 0 && cnt(txtA, 'coupling the covariance carries') == 1;
names = fieldnames(chk);
for i = 1:numel(names)
    if ~chk.(names{i}), fprintf('   FAILED sub-check: %s\n', names{i}); ok = false; end
end
fprintf('A at k: position %s velocity %s, step %+.2f s; D at k-1 %s; E spread %s s; F gain on x+ %s keep %.2f aPV %s; G gain %s aVP %s\n', ...
    sprintf('%+.3f ', resA.coefPos(2, :)), sprintf('%+.3f ', resA.coefVel(2, :)), resA.tStep, sprintf('%+.3f ', resD.coefPos(1, :)), ...
    sprintf('%.2f ', resE.t90(1, :) - resE.t10(1, :)), sprintf('%+.3f ', resF.gainPos(1, :)), resF.keep.position, sprintf('%+.3f ', resF.aPV), ...
    sprintf('%+.3f ', resG.gainPos(1, :)), sprintf('%+.3f ', resG.aVP));
if ok, fprintf('=== test_feedbackCheck PASS ===\n'); else, fprintf('=== test_feedbackCheck FAIL ===\n'); end
end

%------------------------------------------------------------------------
function [run, log] = makeRun(t, tKf, velTrue, applyVel, applyPos, applyBias, tauSlew, gKeep, gApplied, phiMode, undoAt, aErr, bDrift)
if nargin < 11, undoAt = 0; end
if nargin < 12, aErr = 0; end
if nargin < 13, bDrift = 0; end
% a mechanised solution with a velocity error e and a position error pErr. At each epoch the KF increment
% (random, plus a persistent +0.009 m/s on D and +0.05 m on D) is added to the carried state x- to give x+;
% the mechanisation receives gApplied x+ on the axes in applyVel (velocity) and, with applyPos, on every
% position axis, as a jump on the first sample at or after tKf or, with tauSlew > 0, as a first-order slew
% from that sample on; the state keeps (1 - gKeep) x+ and carries it to the next epoch with the Phi coupling
% ('right': position += 0.5 s x velocity; 'transposed': velocity += 0.5 s x position; 'none'). The Z bias
% increment reaches biasEst only with applyBias. The records hold x- and x+.
N = numel(t);  K = numel(tKf);  n = 22;  dt = t(2) - t(1);  DT = 0.5;
lat0 = 0.8;  lon0 = 0.1;  h0 = 100;
P = eye(n);  P(6, 12) = 0.5;  P(12, 6) = 0.5;  P(4, 10) = -0.3;  P(10, 4) = -0.3;  P(1, 4) = 0.7;  P(4, 1) = 0.7;
tUndo = -Inf;  undoV = zeros(1, 3);
e = [0.05 -0.03 0.21];  pErr = [1 -1 2];  biasFb = [1e-3 -2e-3 5e-3];
xm = zeros(1, 6);                                       % carried state: position 1:3, velocity 4:6
pendV = zeros(1, 3);  pendP = zeros(1, 3);
if tauSlew > 0, f = 1 - exp(-dt / tauSlew); else, f = 1; end
velNav = zeros(N, 3);  pErrArr = zeros(N, 3);  biasEst = zeros(K, 6);
log = repmat(struct('priorState', zeros(n, 1), 'postState', zeros(n, 1), 'numMeas', uint8(12), 'postCov', P), K, 1);
kNext = 1;
for i = 1:N
    if kNext <= K && t(i) >= tKf(kNext) - 1e-9
        incV = 0.02 * randn(1, 3) + [0 0 0.009];
        incP = 0.10 * randn(1, 3) + [0 0 0.05];
        incB = 2e-5 * randn(1, 3) + [0 0 bDrift];
        xp = xm + [incP, incV];
        pendV = pendV + gApplied * xp(4:6) .* applyVel;
        pendP = pendP + gApplied * xp(1:3) * applyPos;
        if undoAt > 0, tUndo = tKf(kNext) + undoAt;  undoV = gApplied * xp(4:6) .* applyVel; end
        biasFb = biasFb - incB .* [1 1 applyBias];
        x = zeros(n, 1);  x(1:6) = xp';  x(10:12) = incB';
        xmRec = zeros(n, 1);  xmRec(1:6) = xm';
        log(kNext).priorState = xmRec;
        log(kNext).postState = x;
        biasEst(kNext, :) = [biasFb, zeros(1, 3)];
        r = (1 - gKeep) * xp;                           % the residual kept, then propagated to the next epoch
        if strcmp(phiMode, 'right'),      r(1:3) = r(1:3) + DT * r(4:6);
        elseif strcmp(phiMode, 'transposed'), r(4:6) = r(4:6) + DT * r(1:3);
        end
        xm = r;
        kNext = kNext + 1;
    end
    e = e - f * pendV;  pErr = pErr - f * pendP;
    pendV = (1 - f) * pendV;  pendP = (1 - f) * pendP;
    if t(i) >= tUndo - 1e-9, e = e + undoV;  tUndo = Inf; end         % the velocity correction taken back
    e(3) = e(3) + aErr * dt;                                           % a constant vertical acceleration error
    velNav(i, :) = velTrue(i, :) + e;
    pErrArr(i, :) = pErr;
end
posNed = cumtrapz(t, velNav) + pErrArr;
% lat, lon, h integrated sample by sample with the radii at the sample, the way a mechanisation does
a = 6378137;  e2 = 6.69437999014e-3;
lat = zeros(N, 1);  lon = zeros(N, 1);  h = zeros(N, 1);
lat(1) = lat0;  lon(1) = lon0;  h(1) = h0;
dP = diff(posNed);
for i = 1:N-1
    RN = a * (1 - e2) / (1 - e2 * sin(lat(i))^2)^1.5;
    RE = a / sqrt(1 - e2 * sin(lat(i))^2);
    lat(i+1) = lat(i) + dP(i, 1) / (RN + h(i));
    lon(i+1) = lon(i) + dP(i, 2) / ((RE + h(i)) * cos(lat(i)));
    h(i+1)   = h(i) - dP(i, 3);
end
run = struct('t', t, 'navPos', [lat, lon, h], 'navVel', velNav, 'tKf', tKf, 'biasEst', biasEst, ...
    'tTrue', t, 'trueVel', velTrue, 'biasTrue', [1e-3 -2e-3 5e-3 0 0 0]);
end
