function test_feedbackCheck()
%% test_feedbackCheck.m
% tools/feedbackCheck on a synthetic 100 Hz nav log with 2 Hz updates. The mechanised solution integrates its
% velocity; at each epoch the position and the N and E velocity move by minus the KF increment while the D
% velocity keeps its value and the Z bias is never fed back (case A), or every increment and the bias are
% applied (case B), or nothing is applied at all (case C). The slopes must read -1 where the increment is
% applied and 0 where it is not, the verdicts must name the D velocity and the bias in case A and the missing
% jumps in case C, and the window must select the epochs. Run from anywhere:
% octave --no-gui --quiet --eval test_feedbackCheck
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'tools'));
randn('seed', 3);

dt = 0.01;  T = 400;  t = (0:dt:T)';  N = numel(t);
tKf = (1:0.5:T-1)';  K = numel(tKf);
accTrue = [0.3 * sin(2 * pi * t / 60), 0.2 * cos(2 * pi * t / 45), 0.05 * sin(2 * pi * t / 20)];
velTrue = cumsum(accTrue) * dt + repmat([20 5 0], N, 1);

ok = true;
chk = struct();
% A: D velocity and Z bias not applied
[runA, logA] = makeRun(t, tKf, velTrue, [1 1 0], 1, 0);
txtA = evalc('resA = feedbackCheck(runA, logA);');
fprintf('%s', txtA);
chk.slopeVelA = all(abs(resA.slopeVel(1:2) + 1) < 0.02) && abs(resA.slopeVel(3)) < 0.05;
chk.slopePosA = all(abs(resA.slopePos + 1) < 0.02);
chk.seA       = all(resA.seVel < 0.02) && all(resA.sePos < 0.02);
chk.usedA     = resA.nUsed == K - nnz(mod(1:K, 40) == 0);
chk.incA      = abs(mean(resA.incVel(resA.used, 3)) - 0.009) < 0.003;
chk.biasA     = abs(resA.dBias(3) / resA.incBias(3)) < 0.05 && all(abs(resA.dBias(1:2) ./ resA.incBias(1:2) + 1) < 1e-6);
chk.verdictA  = ~isempty(strfind(txtA, 'the D velocity increment does not reach')) && ...
                ~isempty(strfind(txtA, 'bias increments reach biasEst on some axes only')) && ...
                isempty(strfind(txtA, 'position increment does not reach'));
% B: everything applied, cell input, biasEst per nav sample, a window
[runB, logB] = makeRun(t, tKf, velTrue, [1 1 1], 1, 1);
runB.biasEst = interp1(tKf, runB.biasEst, t, 'previous', 'extrap');
txtB = evalc('resB = feedbackCheck(runB, num2cell(logB), [100 300]);');
chk.slopeB    = all(abs([resB.slopeVel, resB.slopePos] + 1) < 0.02);
chk.winB      = resB.nUsed == nnz(tKf >= 100 & tKf <= 300 & mod((1:K)', 40) ~= 0) && all(resB.t(resB.used) >= 100) && all(resB.t(resB.used) <= 300);
chk.biasB     = all(abs(resB.dBias ./ resB.incBias + 1) < 1e-6);
chk.quietB    = isempty(strfind(txtB, 'VERDICT'));
% C: nothing applied
[runC, logC] = makeRun(t, tKf, velTrue, [0 0 0], 0, 0);
txtC = evalc('resC = feedbackCheck(runC, logC);');
chk.slopeC    = all(abs([resC.slopeVel, resC.slopePos]) < 0.05);
chk.verdictC  = ~isempty(strfind(txtC, 'no jump on any axis'));
names = fieldnames(chk);
for i = 1:numel(names)
    if ~chk.(names{i}), fprintf('   FAILED sub-check: %s\n', names{i}); ok = false; end
end
fprintf('A velocity slopes %s, position slopes %s, bias ratio Z %+.3f; B slopes %s on %d updates; C slopes %s\n', ...
    sprintf('%+.3f ', resA.slopeVel), sprintf('%+.3f ', resA.slopePos), resA.dBias(3) / resA.incBias(3), ...
    sprintf('%+.3f ', [resB.slopeVel, resB.slopePos]), resB.nUsed, sprintf('%+.3f ', [resC.slopeVel, resC.slopePos]));
if ok, fprintf('=== test_feedbackCheck PASS ===\n'); else, fprintf('=== test_feedbackCheck FAIL ===\n'); end
end

%------------------------------------------------------------------------
function [run, log] = makeRun(t, tKf, velTrue, applyVel, applyPos, applyBias)
% a mechanised solution with a velocity error e and a position error pErr; at each epoch the KF increment
% (random, plus a persistent +0.009 m/s on D and +0.05 m on D) is subtracted from the axes in applyVel (velocity)
% and, with applyPos, from every position axis; the Z bias increment reaches biasEst only with applyBias. Every 40th epoch has
% no rows. The update lands on the first sample at or after tKf.
N = numel(t);  K = numel(tKf);  n = 22;
lat0 = 0.8;  lon0 = 0.1;  h0 = 100;
e = [0.05 -0.03 0.21];  pErr = [1 -1 2];  biasFb = [1e-3 -2e-3 5e-3];
velNav = zeros(N, 3);  pErrArr = zeros(N, 3);  biasEst = zeros(K, 6);
log = repmat(struct('priorState', zeros(n, 1), 'postState', zeros(n, 1), 'numMeas', uint8(0)), K, 1);
kNext = 1;
for i = 1:N
    if kNext <= K && t(i) >= tKf(kNext) - 1e-9
        if mod(kNext, 40) == 0
            incV = zeros(1, 3);  incP = zeros(1, 3);  incB = zeros(1, 3);  m = 0;
        else
            incV = 0.02 * randn(1, 3) + [0 0 0.009];
            incP = 0.10 * randn(1, 3) + [0 0 0.05];
            incB = 2e-5 * randn(1, 3) + [0 0 1e-5];
            m = 12;
        end
        e    = e - incV .* applyVel;
        pErr = pErr - incP * applyPos;
        biasFb = biasFb - incB .* [1 1 applyBias];
        x = zeros(n, 1);  x(1:3) = incP';  x(4:6) = incV';  x(10:12) = incB';
        log(kNext).postState = x;
        log(kNext).numMeas = uint8(m);
        biasEst(kNext, :) = [biasFb, zeros(1, 3)];
        kNext = kNext + 1;
    end
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
