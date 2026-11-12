function test_navConsistency()
%% test_navConsistency.m
% tools/navConsistency on a synthetic stationary run: errors drawn from known sigmas on position, velocity, attitude
% and biases, and a logged P that is either those sigmas (ratios near 1, about 0.3 % beyond 3 sigma), three times
% them (ratios near 3, nothing beyond 3 sigma) or a third of them (ratios near 0.33, many beyond 3 sigma); the
% window argument must restrict the epochs read; the three forms of P must read the same; a constant 3 m offset on
% the north position must read a low ratio and a mean-removed ratio near 1, with the mean reported. Run from anywhere:
% octave --no-gui --quiet --eval test_navConsistency
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'tools'));
rng(5);
K = 4000;  tKf = 0.5 * (1:K)';  n = 22;
lat0 = 0.43;  lon0 = 0.81;  h0 = 43;
a = 6378137;  e2 = 6.69437999014e-3;
Rm = a * (1 - e2) / (1 - e2 * sin(lat0)^2)^1.5;  Rn = a / sqrt(1 - e2 * sin(lat0)^2);
sPos = [2 3 4];  sVel = [0.03 0.03 0.05];  sAtt = [5e-4 5e-4 2e-2];  sAcc = 3e-4 * [1 1 1];  sGyr = 2e-6 * [1 1 1];
ep = randn(K, 3) .* sPos;  ev = randn(K, 3) .* sVel;  ea = randn(K, 3) .* sAtt;  eb = [randn(K, 3) .* sAcc, randn(K, 3) .* sGyr];
biasTrue = [1e-3 -2e-3 3e-3 1e-5 -1e-5 2e-5];  trueAtt = [0.01 0.2 3.1];
run = struct('t', tKf, 'navPos', [lat0 + ep(:, 1) / (Rm + h0), lon0 + ep(:, 2) / ((Rn + h0) * cos(lat0)), h0 - ep(:, 3)], ...
    'navVel', ev, 'navAtt', repmat(trueAtt, K, 1) + ea, 'biasEst', repmat(biasTrue, K, 1) + eb, ...
    'tTrue', [0; tKf(end) + 1], 'truePos', [lat0 lon0 h0; lat0 lon0 h0], 'trueVel', zeros(2, 3), 'trueAtt', [trueAtt; trueAtt], ...
    'biasTrue', biasTrue, 'pGyroFirst', true);
sig = zeros(1, n);  sig(1:3) = sPos;  sig(4:6) = sVel;  sig(7:9) = sAtt;  sig(10:12) = sGyr;  sig(13:15) = sAcc;  sig(16:22) = 1;
pDiag = repmat(sig .^ 2, K, 1);
fprintf('--- a consistent P ---\n');
r1 = navConsistency(run, tKf, pDiag);
% the tilt group compares a magnitude (Rayleigh) with a quadrature sigma: its ratio sits near 1 as well
ok1 = all(abs(r1.ratio - 1) < 0.1) && all(r1.beyond3 < 0.01) && numel(r1.group) == 10;
fprintf('--- a P three times the errors ---\n');
r3 = navConsistency(run, tKf, 9 * pDiag);
ok3 = all(abs(r3.ratio - 3) < 0.3) && all(r3.beyond3 == 0);
fprintf('--- an overconfident P, a third of the errors ---\n');
rT = navConsistency(run, tKf, pDiag / 9);
okT = all(abs(rT.ratio - 1 / 3) < 0.04) && all(rT.beyond3(1:3) > 0.2);
fprintf('--- a window, and the other forms of P ---\n');
rW = navConsistency(run, tKf, pDiag, [100 300; 1500 1600]);
okW = rW.K == nnz((tKf >= 100 & tKf <= 300) | (tKf >= 1500 & tKf <= 1600));
pCell = cell(K, 1);  for k = 1:K, pCell{k} = diag(sig .^ 2); end
evalc('rC = navConsistency(run, tKf, pCell);');
P3 = repmat(diag(sig .^ 2), [1 1 K]);
evalc('r3d = navConsistency(run, tKf, P3);');
okForm = max(abs(rC.ratio - r1.ratio)) < 1e-9 && max(abs(r3d.ratio - r1.ratio)) < 1e-9;
fprintf('--- a constant offset of 3 m on the north position ---\n');
runOff = run;  runOff.navPos(:, 1) = runOff.navPos(:, 1) + 3 / (Rm + h0);
rO = navConsistency(runOff, tKf, pDiag);
okOff = rO.ratio(1) < 0.7 && abs(rO.ratioC(1) - 1) < 0.1 && abs(rO.meanErr(1) - 3) < 0.2 && abs(rO.ratioC(2) - r1.ratioC(2)) < 1e-9 ...
    && all(abs(r1.ratioC - r1.ratio) < 0.05);
fprintf('consistent P reads near 1: %s | three times reads 3: %s | a third reads 0.33 with many beyond 3 sigma: %s | window: %s | forms of P: %s | offset: %s\n', ...
    pf(ok1), pf(ok3), pf(okT), pf(okW), pf(okForm), pf(okOff));
if ok1 && ok3 && okT && okW && okForm && okOff, fprintf('=== test_navConsistency PASS ===\n');
else, fprintf('=== test_navConsistency FAIL ===\n'); end
end

function s = pf(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end
