function test_feedbackCheck()
%% test_feedbackCheck.m
% tools/feedbackCheck on a synthetic 100 Hz nav log with 2 Hz updates. The mechanised solution integrates its
% velocity; at each epoch the KF increment is subtracted from the position and the velocity as a jump on the
% first sample at or after tKf (A: D velocity and Z bias not applied; B: everything applied; C: nothing
% applied), or the records are offset by one epoch from the nav log (D), or the correction is slewed into the
% solution with a 0.1 s time constant (E). The lag coefficients, the step times and the verdicts must read
% each case. Run from anywhere:
% octave --no-gui --quiet --eval test_feedbackCheck
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'tools'));
randn('seed', 3);

dt = 0.01;  T = 400;  t = (0:dt:T)';  N = numel(t);
tKf = (1:0.5:T-1)';  K = numel(tKf);
accTrue = [0.3 * sin(2 * pi * t / 60), 0.2 * cos(2 * pi * t / 45), 0.05 * sin(2 * pi * t / 20)];
velTrue = cumsum(accTrue) * dt + repmat([20 5 0], N, 1);
near = @(x, v, tol) all(abs(x(:) - v) < tol);

ok = true;
chk = struct();
% A: D velocity and Z bias not applied, sharp jumps
[runA, logA] = makeRun(t, tKf, velTrue, [1 1 0], 1, 0, 0);
txtA = evalc('resA = feedbackCheck(runA, logA);');
fprintf('%s', txtA);
chk.posA      = near(resA.coefPos(2, :), -1, 0.03) && near(resA.coefPos([1 3], :), 0, 0.03);
chk.velA      = near(resA.coefVel(2, 1:2), -1, 0.03) && abs(resA.coefVel(2, 3)) < 0.05 && near(resA.coefVel([1 3], 1:2), 0, 0.03);
chk.stepA     = near(resA.t50(1, :), 0, 0.011) && near(resA.t90(1, :) - resA.t10(1, :), 0, 0.011) && near(resA.t50(2, 1:2), 0, 0.011);
chk.usedA     = resA.nUsed == K;
chk.incA      = abs(mean(resA.incVel(resA.used, 3)) - 0.009) < 0.003;
chk.biasA     = abs(resA.dBias(3) / resA.incBias(3)) < 0.05 && near(resA.dBias(1:2) ./ resA.incBias(1:2), -1, 1e-6);
chk.verdictA  = ~isempty(strfind(txtA, 'the D velocity increment does not reach')) && ...
                ~isempty(strfind(txtA, 'bias increments reach biasEst on some axes only')) && ...
                numel(strfind(txtA, 'VERDICT')) == 2;
% B: everything applied, cell input, biasEst per nav sample, a window
[runB, logB] = makeRun(t, tKf, velTrue, [1 1 1], 1, 1, 0);
runB.biasEst = interp1(tKf, runB.biasEst, t, 'previous', 'extrap');
txtB = evalc('resB = feedbackCheck(runB, num2cell(logB), [100 300]);');
chk.coefB     = near(resB.coefPos(2, :), -1, 0.03) && near(resB.coefVel(2, :), -1, 0.03);
chk.winB      = resB.nUsed == nnz(tKf >= 100 & tKf <= 300) && all(resB.t(resB.used) >= 100) && all(resB.t(resB.used) <= 300);
chk.biasB     = near(resB.dBias ./ resB.incBias, -1, 1e-6);
chk.quietB    = isempty(strfind(txtB, 'VERDICT'));
% C: nothing applied
[runC, logC] = makeRun(t, tKf, velTrue, [0 0 0], 0, 0, 0);
txtC = evalc('resC = feedbackCheck(runC, logC);');
chk.coefC     = near(resC.coefPos, 0, 0.05) && near(resC.coefVel, 0, 0.05);
chk.verdictC  = numel(strfind(txtC, 'does not reach the logged solution')) == 6;
% D: the records one epoch ahead of the nav log (record k holds the increment applied at tKf(k+1))
[runD, logD] = makeRun(t, tKf, velTrue, [1 1 1], 1, 1, 0);
logD = [logD(2:end); logD(1)];  logD(end).postState(:) = 0;
txtD = evalc('resD = feedbackCheck(runD, logD);');
chk.coefD     = near(resD.coefPos(1, :), -1, 0.03) && near(resD.coefPos([2 3], :), 0, 0.03) && near(resD.coefVel(1, :), -1, 0.03);
chk.verdictD  = numel(strfind(txtD, 'one epoch ahead of the nav log')) == 6;
% E: the correction slewed in with a 0.1 s time constant
[runE, logE] = makeRun(t, tKf, velTrue, [1 1 1], 1, 1, 0.1);
txtE = evalc('resE = feedbackCheck(runE, logE);');
chk.totalE    = near(sum(resE.coefPos, 1), -1, 0.05) && near(resE.coefPos(2, :), -0.92, 0.05);
chk.stepE     = all(resE.t90(1, :) - resE.t10(1, :) > 0.15) && near(resE.t10(1, :), 0.01, 0.015);
chk.verdictE  = numel(strfind(txtE, 'slewed feedback')) >= 3;
names = fieldnames(chk);
for i = 1:numel(names)
    if ~chk.(names{i}), fprintf('   FAILED sub-check: %s\n', names{i}); ok = false; end
end
fprintf('A coefficients at k: position %s velocity %s, step %s s, bias ratio Z %+.3f; B %d updates; D at k-1 %s; E spread %s s\n', ...
    sprintf('%+.3f ', resA.coefPos(2, :)), sprintf('%+.3f ', resA.coefVel(2, :)), sprintf('%+.2f ', resA.t50(1, :)), ...
    resA.dBias(3) / resA.incBias(3), resB.nUsed, sprintf('%+.3f ', resD.coefPos(1, :)), sprintf('%.2f ', resE.t90(1, :) - resE.t10(1, :)));
if ok, fprintf('=== test_feedbackCheck PASS ===\n'); else, fprintf('=== test_feedbackCheck FAIL ===\n'); end
end

%------------------------------------------------------------------------
function [run, log] = makeRun(t, tKf, velTrue, applyVel, applyPos, applyBias, tauSlew)
% a mechanised solution with a velocity error e and a position error pErr; at each epoch the KF increment
% (random, plus a persistent +0.009 m/s on D and +0.05 m on D) is subtracted from the axes in applyVel
% (velocity) and, with applyPos, from every position axis, as a jump on the first sample at or after tKf or,
% with tauSlew > 0, as a first-order slew from that sample on. The Z bias increment reaches biasEst only
% with applyBias.
N = numel(t);  K = numel(tKf);  n = 22;  dt = t(2) - t(1);
lat0 = 0.8;  lon0 = 0.1;  h0 = 100;
e = [0.05 -0.03 0.21];  pErr = [1 -1 2];  biasFb = [1e-3 -2e-3 5e-3];
pendV = zeros(1, 3);  pendP = zeros(1, 3);
if tauSlew > 0, f = 1 - exp(-dt / tauSlew); else, f = 1; end
velNav = zeros(N, 3);  pErrArr = zeros(N, 3);  biasEst = zeros(K, 6);
log = repmat(struct('priorState', zeros(n, 1), 'postState', zeros(n, 1), 'numMeas', uint8(12)), K, 1);
kNext = 1;
for i = 1:N
    if kNext <= K && t(i) >= tKf(kNext) - 1e-9
        incV = 0.02 * randn(1, 3) + [0 0 0.009];
        incP = 0.10 * randn(1, 3) + [0 0 0.05];
        incB = 2e-5 * randn(1, 3) + [0 0 1e-5];
        pendV = pendV + incV .* applyVel;
        pendP = pendP + incP * applyPos;
        biasFb = biasFb - incB .* [1 1 applyBias];
        x = zeros(n, 1);  x(1:3) = incP';  x(4:6) = incV';  x(10:12) = incB';
        log(kNext).postState = x;
        biasEst(kNext, :) = [biasFb, zeros(1, 3)];
        kNext = kNext + 1;
    end
    e = e - f * pendV;  pErr = pErr - f * pendP;
    pendV = (1 - f) * pendV;  pendP = (1 - f) * pendP;
    velNav(i, :) = velTrue(i, :) + e;
    pErrArr(i, :) = pErr;
end
posNed = cumtrapz(t, velNav) + pErrArr;
a = 6378137;  e2 = 6.69437999014e-3;
RN = a * (1 - e2) / (1 - e2 * sin(lat0)^2)^1.5;
RE = a / sqrt(1 - e2 * sin(lat0)^2);
h   = h0 - posNed(:, 3);
lat = lat0 + posNed(:, 1) ./ (RN + h);
lon = lon0 + posNed(:, 2) ./ ((RE + h) * cos(lat0));
run = struct('t', t, 'navPos', [lat, lon, h], 'navVel', velNav, 'tKf', tKf, 'biasEst', biasEst);
end
