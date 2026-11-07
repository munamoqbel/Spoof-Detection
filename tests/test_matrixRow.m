function test_matrixRow()
%% test_matrixRow.m
% tools/matrixRow and tools/matrixTable on a synthetic gate log and nav run: a 20 m north step at 100 to 130 s that
% latches 2.5 s in on the SS N monitor, coasts, is vetoed once, hands back 9 s after tOff, with one false latch at
% 49.5 s and a nav error that follows the step; the same log read for a scenario at 150 s that nothing answered
% (NOMINAL, no latch) and without the run. Run from anywhere: octave --no-gui --quiet --eval test_matrixRow
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'tools'));
K = 400;  tKf = 0.5 * (0:K-1)';
L = SPF_logInit(K);  L.count = K;
L.mode(:) = 1;  L.armed(31:end) = true;
L.mode(100:105) = 2;  L.mode(106:109) = 3;                                   % a false latch at 49.5 s, handed back at 54.5 s
L.mode(206:260) = 2;  L.mode(261:278) = 3;                                   % the latch at 102.5 s, probation, hand-back at 139 s
L.evLatch = [100 206];  L.evProbStart = [106 261];  L.evVeto = 240;  L.evHandback = [110 279];
L.ssRatio = 0.3 + 0.1 * rand(K, 3);  L.cpiRatio = 0.2 + 0.1 * rand(K, 3);
L.ssRatio(203:206, 1) = [0.8 1.0 1.2 1.3];  L.ssAlarm(206) = true;  L.alarmAxis(206, 1) = true;
L.ssRatio(100, 2) = 1.1;  L.ssAlarm(100) = true;  L.alarmAxis(100, 2) = true;
L.cpiRatio(230, 1) = 0.9;
L.sigmaPos(:, :) = 3;  L.sigmaPos(279, :) = [12 11 20];  L.pl(:) = 10;  L.pl(279) = 25;
L.numMeas(:) = 14;  L.numMeas(206:250) = 9;                                 % the host rejects rows through the coast
L.revalDone(206:278) = true;  L.revalThr(206:278) = 36;  L.qReval(206:250) = 100;  L.qReval(251:278) = 10;   % 10 passes from 125 s
L.dwell(251:260) = 1:10;
L.inputFault(300) = true;
% the nav run: truth stationary, the nav 2 m north and 1 m east of it, pulled 20 m further north by the spoof
t = (0:0.01:199.99)';
lat0 = 0.43;  lon0 = 0.81;  h0 = 43;
a = 6378137;  e2 = 6.69437999014e-3;
Rm = a * (1 - e2) / (1 - e2 * sin(lat0)^2)^1.5;  Rn = a / sqrt(1 - e2 * sin(lat0)^2);
errN = 2 + 20 * min(max((t - 100) / 3, 0), 1) .* (t < 139);                 % ramps to 20 m over 3 s, back at the hand-back
errE = 1 * ones(size(t));
run = struct('t', t, 'navPos', [lat0 + errN / (Rm + h0), lon0 + errE / ((Rn + h0) * cos(lat0)), h0 * ones(size(t))], ...
    'tTrue', [0; 200], 'truePos', [lat0 lon0 h0; lat0 lon0 h0]);
scn = struct('kind', 'step', 'tOn', 100, 'tOff', 130, 'size', 20, 'dirNed', [1 0 0]);
fprintf('--- latched run ---\n');
r1 = matrixRow(L, tKf, run, scn, 'run 3');
scn2 = struct('kind', 'step', 'tOn', 150, 'tOff', 180, 'size', 5, 'dirNed', [1 0 0]);
fprintf('--- nothing answered ---\n');
r2 = matrixRow(L, tKf, run, scn2);
fprintf('--- without the run ---\n');
r3 = matrixRow(L, tKf, [], scn, 'run 3 gate');
fprintf('--- the table ---\n');
matrixTable([r1 r2 r3]);
matrixTable({r1, r2});
ok1 = abs(r1.latchDelay - 2.5) < 1e-9 && strcmp(r1.monitor, 'SS') && strcmp(r1.axis, 'N') && abs(r1.ratioAtLatch - 1.3) < 1e-9 ...
   && abs(r1.handbackDelay - 9) < 1e-9 && abs(r1.coastLen - 36.5) < 1e-9 && r1.nVeto == 1 && r1.falseAlarms == 1 && abs(r1.tFalse - 49.5) < 1e-9 ...
   && abs(r1.maxSsIn - 1.3) < 1e-9 && r1.maxCpiIn < 0.35 && abs(r1.maxSsClean - 1.1) < 1e-9 && r1.maxCpiClean < 0.35 ...   % the false latch is on the clean epochs; the 0.9 at epoch 230 is the coast's
   && abs(r1.errAtLatch(1) - 20 * (2.25 / 3)) < 1.5 && abs(r1.errAtLatch(2)) < 1e-6 && abs(r1.maxHorizIn - 20) < 1e-6 ...
   && abs(r1.errAtOff(1) - 20) < 1e-6 && abs(r1.meanErrIn(1) - 20) < 1e-6 && abs(r1.errAtHandback(1) - 20) < 0.5 ...
   && isequal(r1.sigmaAtHandback, [12 11 20]) && isequal(r1.faults, [1 0 0 0]) && ~isempty(strfind(r1.timeline, 'C 102.5')) ...
   && ~isempty(strfind(r1.timeline, 'P 130.0')) && ~isempty(strfind(r1.timeline, 'N 139.0'));
base2 = mean(errN(t >= 90 & t < 150));                                      % the minute before 150 s still holds the first spoof
ok2 = isnan(r2.latchDelay) && strcmp(r2.monitor, 'none') && isnan(r2.handbackDelay) && r2.maxSsIn < 0.5 && r2.falseAlarms == 2 ...
   && abs(r2.meanErrIn(1) - (2 - base2)) < 1e-6 && abs(r2.errAtOff(1) - (2 - base2)) < 1e-6 && strcmp(r2.name, 'step 5') ...
   && abs(r2.maxSsClean - 1.3) < 1e-9;
ok3 = abs(r3.latchDelay - 2.5) < 1e-9 && all(isnan(r3.errAtLatch)) && isnan(r3.maxHorizIn) && isequal(r3.sigmaAtHandback, [12 11 20]);
okC = r1.revalComputed == 73 && r1.revalPassed == 28 && r1.rowsCoast == 9 && r1.rowsBefore == 14 && abs(r1.tFirstPass - tKf(251)) < 1e-9 ...
   && ~r1.probationInSpoof && abs(r1.qRatioCoast - 100 / 36) < 1e-9 && isnan(r2.revalPassed);
% the same log read for a scenario that ends after the probation: a hand-back into the spoof, with fewer rows
scn3 = struct('kind', 'step', 'tOn', 100, 'tOff', 160, 'size', 20, 'dirNed', [1 0 0]);
fprintf('--- the probation opening while the spoof is on ---\n');
r4 = matrixRow(L, tKf, [], scn3, 'into the spoof');
% the latch read as an onset and a return in COAST at 112 s; the course table carries the landing column
L.returnWander(206) = 1;  L.returnLanding(206) = 40;  L.returnLanding(208) = 31.5;
L.evReturn = 225;  L.returnWander(225) = 22;  L.returnLanding(225) = 4;
fprintf('--- a return in COAST ---\n');
r5 = matrixRow(L, tKf, [], scn, 'run 3 return');
matrixTable([r1 r5]);
okC = okC && abs(r5.tReturn - tKf(225)) < 1e-9 && isnan(r1.tReturn) && isnan(r4.tReturn) ...
   && r4.probationInSpoof && abs(r4.handbackDelay - (139 - 160)) < 1e-9 && r4.plAtHandback == 25 && r1.plAtHandback == 25 ...
   && isfinite(r1.coastCpiMin) && r1.coastCpiMin < 0.35 && r1.coastCpiMed < 0.35 && isnan(r2.coastCpiMin);
if ok1 && ok2 && ok3 && okC, fprintf('=== test_matrixRow PASS ===\n');
else, fprintf('ok1 %d ok2 %d ok3 %d okC %d\n=== test_matrixRow FAIL ===\n', ok1, ok2, ok3, okC); end
end
