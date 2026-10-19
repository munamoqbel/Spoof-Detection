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
[logB, sysB] = makeLog(K, n, mMax, 1:5, 50, 60, -40);
fprintf('--- B: -40 m on five of eight rows (a majority), epochs 50 to 60 ---\n');
pB = spoofProfile(logB, tKf, sysB);
[logC, ~] = makeLog(K, n, mMax, 1:4, 1, 0, 0);
fprintf('--- C: clean ---\n');
pC = spoofProfile(logC, tKf);
okA = abs(pA.onset - tKf(80)) < 1e-9 && abs(pA.tEnd - tKf(140)) < 1e-9 && abs(pA.sizeOnset - 300) < 5 && pA.sysOff == 1 ...
   && pA.nOff(100) == 4 && abs(pA.perSys(100, 1) - 300) < 5 && abs(pA.perSys(100, 2)) < 5 ...
   && abs(pA2.onset - tKf(80)) < 1e-9 && isnan(pA2.sysOff) && abs(pA2.sizeOnset - 300) < 5;
okB = abs(pB.onset - tKf(50)) < 1e-9 && abs(pB.tEnd - tKf(60)) < 1e-9 && abs(pB.sizeOnset + 40) < 5 && pB.nOff(55) == 5 && pB.sysOff == 1;
okC = isnan(pC.onset) && ~any(pC.present);
if okA && okB && okC, fprintf('=== test_spoofProfile PASS ===\n'); else, fprintf('okA %d okB %d okC %d\n=== test_spoofProfile FAIL ===\n', okA, okB, okC); end
end
