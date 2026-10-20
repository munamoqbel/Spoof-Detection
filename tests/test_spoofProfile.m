function test_spoofProfile()
%% test_spoofProfile.m
% tools/spoofProfile on synthetic spfMeas records: a 300 m offset on four of eight range rows between two epochs
% must be found with its onset, end, size and system; a majority spoofed group (five of eight) must still be the
% offset group; a clean log must report nothing. Run from anywhere: octave --no-gui --quiet --eval test_spoofProfile
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'tools'));
randn('seed', 9);  rand('seed', 9);
K = 200;  n = 22;  mMax = 51;  tKf = 0.5 * (0:K-1)';
function [log, sysId] = makeLog(K, n, mMax, offRows, kOn, kEnd, size)
    log = cell(K, 1);  sysId = zeros(K, mMax);
    for k = 1:K
        el = 0.2 + 0.8 * rand(8, 1);  az = 2 * pi * rand(8, 1);
        u  = [cos(el) .* cos(az), cos(el) .* sin(az), -sin(el)];
        H = zeros(mMax, n);  y = zeros(mMax, 1);
        for j = 1:8, H(j, 1:3) = -u(j, :);  H(j, 16) = 1;  y(j) = 2 * randn; end           % ranges
        for j = 1:8, H(8 + j, 4:6) = -u(j, :);  H(8 + j, 17) = 1;  y(8 + j) = 0.05 * randn; end   % range rates
        H(17, 3) = 1;  y(17) = 1.5 * randn;                                                % pressure
        if k >= kOn && k <= kEnd, y(offRows) = y(offRows) + size; end
        sysId(k, 1:8) = [1 1 1 1 2 2 2 2];  sysId(k, 9:16) = [1 1 1 1 2 2 2 2];  sysId(k, 17) = 0;
        log{k} = struct('innovation', y, 'innovationCov', eye(mMax), 'obsMatrix', H, 'measNoiseCov', eye(mMax), ...
            'numMeas', 17, 'priorState', zeros(n, 1), 'postState', zeros(n, 1), 'postCov', eye(n));
    end
end
[logA, sysA] = makeLog(K, n, mMax, 1:4, 80, 140, 300);
fprintf('--- A: 300 m on system 1 (four of eight rows), epochs 80 to 140 ---\n');
pA = spoofProfile(logA, tKf, sysA);
pA2 = spoofProfile(logA, tKf);
% the same through the host's pairing log: 4 used GPS slots and 4 used GLONASS slots per epoch, slot order
S = 12;
blank = struct('corrValid', true, 'geomRange', zeros(S, 1), 'geomRangeRate', zeros(S, 1));
blankSt = struct('position', zeros(S, 3), 'velocity', zeros(S, 3), 'LLA', zeros(S, 3), 'satValid', 1);
recL = struct('prLog', struct('gpsRange', blank, 'gloRange', blank), 'satPosLog', struct('navstarStatus', blankSt, 'glonassStatus', blankSt), 'sysIdLog', []);
loggedA = repmat(recL, K, 1);  loggedM = repmat(recL, K, 1);
for k = 1:K
    g = zeros(S, 1);  g([2 5 7 11]) = 2.2e7 + 1e5 * rand(4, 1);
    loggedA(k).prLog.gpsRange.geomRange = g;  loggedA(k).prLog.gloRange.geomRange = g([3 1 2 4 5 6 7 8 9 10 11 12]);
    loggedA(k).satPosLog.navstarStatus.position = 2.6e7 * (rand(S, 3) - 0.5);  loggedA(k).satPosLog.glonassStatus.position = 2.6e7 * (rand(S, 3) - 0.5);
    loggedA(k).sysIdLog = [1 1 1 1 2 2 2 2]';
    loggedM(k) = loggedA(k);  gM = g;  gM(9) = 2.2e7;  loggedM(k).prLog.gpsRange.geomRange = gM;   % five GPS slots: counts differ
end
fprintf('--- A through the pairing log ---\n');
pA3 = spoofProfile(logA, tKf, loggedA);
fprintf('--- A through a pairing log with one used slot too many ---\n');
pA4 = spoofProfile(logA, tKf, loggedM);
% the injection from the spoofed pairing log minus the clean one: +300 m on the GPS slots, epochs 80 to 140
loggedS = loggedA;
for k = 80:140, loggedS(k).prLog.gpsRange.geomRange([2 5 7 11]) = loggedA(k).prLog.gpsRange.geomRange([2 5 7 11]) + 300; end
fprintf('--- A with the clean twin: the injection itself ---\n');
pA5 = spoofProfile(logA, tKf, loggedS, loggedA);
[logB, sysB] = makeLog(K, n, mMax, 1:5, 50, 60, -40);
fprintf('--- B: -40 m on five of eight rows (a majority), epochs 50 to 60 ---\n');
pB = spoofProfile(logB, tKf, sysB);
[logC, ~] = makeLog(K, n, mMax, 1:4, 1, 0, 0);
fprintf('--- C: clean ---\n');
pC = spoofProfile(logC, tKf);
okA = abs(pA.onset - tKf(80)) < 1e-9 && abs(pA.tEnd - tKf(140)) < 1e-9 && abs(pA.sizeOnset - 300) < 5 && pA.sysOff == 1 ...
   && pA.nOff(100) == 4 && abs(pA.perSys(100, 1) - 300) < 5 && abs(pA.perSys(100, 2)) < 5 ...
   && abs(pA2.onset - tKf(80)) < 1e-9 && isnan(pA2.sysOff) && abs(pA2.sizeOnset - 300) < 5 ...
   && pA3.sysOff == 1 && pA3.nSysMismatch == 0 && abs(pA3.perSys(100, 1) - 300) < 5 && abs(pA3.perSys(100, 2)) < 5 ...
   && pA4.nSysMismatch == K ...
   && abs(pA5.injOnset - tKf(80)) < 1e-9 && abs(pA5.injEnd - tKf(140)) < 1e-9 && abs(pA5.injSize(1) - 300) < 1e-6 ...
   && abs(pA5.injSize(2)) < 1e-6 && pA5.nSlotChanged == 0;
okB = abs(pB.onset - tKf(50)) < 1e-9 && abs(pB.tEnd - tKf(60)) < 1e-9 && abs(pB.sizeOnset + 40) < 5 && pB.nOff(55) == 5 && pB.sysOff == 1;
okC = isnan(pC.onset) && ~any(pC.present);
if okA && okB && okC, fprintf('=== test_spoofProfile PASS ===\n'); else, fprintf('okA %d okB %d okC %d\n=== test_spoofProfile FAIL ===\n', okA, okB, okC); end
end
