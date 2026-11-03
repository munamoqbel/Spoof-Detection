function test_spoofInjectCheck()
%% test_spoofInjectCheck.m
% tools/spoofInjectCheck on a synthetic pairing log of three systems (stationary receiver, clock ramp, inter-system
% biases, code noise) and its spoofed copies made with spoofInject on systems 1 and 2 only: a 5 m north step, a
% 1 m/s east ramp (against a clean twin shifted by two epochs, matched on the satellite positions) and a 30 m clock
% step must be read back exactly (the rates to the second-order term of the geometry), system 3 must report no injection, and a clean twin that differs on six epochs
% outside the interval must be reported as such. Run from anywhere: octave --no-gui --quiet --eval test_spoofInjectCheck
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'tools'));
randn('seed', 11);  rand('seed', 11);
K = 80;  tKf = 0.5 * (0:K-1)';
lat0 = 0.43;  lon0 = 0.81;  h0 = 43;
a = 6378137;  e2 = 6.69437999014e-3;
Nr = a / sqrt(1 - e2 * sin(lat0)^2);
refEcef = [(Nr + h0) * cos(lat0) * cos(lon0), (Nr + h0) * cos(lat0) * sin(lon0), (Nr * (1 - e2) + h0) * sin(lat0)];
S = 15;  used = {[1 2 3 5 9 11], [2 4 6 7 10 13], [3 8 10 12 14]};  isb = [0 40 -25];  drift0 = 0.3;
up = refEcef / norm(refEcef);
dirs = randn(S * 3, 3);  dirs = dirs - (dirs * up') * up;  dirs = dirs ./ repmat(sqrt(sum(dirs.^2, 2)), 1, 3);
elv = 0.2 + 0.6 * rand(S * 3, 1);
dirs = repmat(cos(elv), 1, 3) .* dirs + repmat(sin(elv), 1, 3) .* repmat(up, S * 3, 1);
satAll = repmat(refEcef, S * 3, 1) + 2.0e7 * dirs;
velAll = randn(S * 3, 3);  velAll = 3000 * velAll ./ repmat(sqrt(sum(velAll.^2, 2)), 1, 3);
blank = struct('corrValid', true, 'geomRange', zeros(S, 1), 'geomRangeRate', zeros(S, 1));
blankSt = struct('position', zeros(S, 3), 'velocity', zeros(S, 3), 'LLA', zeros(S, 3), 'satValid', 1);
rec = struct('prLog', struct('gpsRange', blank, 'gloRange', blank, 'galRange', blank), ...
    'satPosLog', struct('navstarStatus', blankSt, 'glonassStatus', blankSt, 'galileoStatus', blankSt), 'sysIdLog', []);
prN = {'gpsRange', 'gloRange', 'galRange'};  spN = {'navstarStatus', 'glonassStatus', 'galileoStatus'};
clean = repmat(rec, K, 1);
for k = 1:K
    clk = 1500 + 30 * tKf(k);
    for s = 1:3
        g = zeros(S, 1);  rr = zeros(S, 1);  pos = zeros(S, 3);  vel = zeros(S, 3);
        for j = used{s}
            i = (s - 1) * S + j;
            sat = satAll(i, :) + velAll(i, :) * tKf(k);                        % the sky moves: every epoch its own
            u = (sat - refEcef) / norm(sat - refEcef);
            g(j) = norm(sat - refEcef) + clk + isb(s) + 1.0 * randn;
            rr(j) = u * velAll(i, :)' + drift0 + 0.05 * randn;
            pos(j, :) = sat;  vel(j, :) = velAll(i, :);
        end
        clean(k).prLog.(prN{s}).geomRange = g;  clean(k).prLog.(prN{s}).geomRangeRate = rr;
        clean(k).satPosLog.(spN{s}).position = pos;  clean(k).satPosLog.(spN{s}).velocity = vel;
    end
end
function [sp] = spoofed(clean, scn)
    % the clean log with spoofInject on systems 1 and 2 at every epoch
    sp = clean;
    for k = 1:numel(sp)
        inj = spoofScenario(tKf(k), scn);
        for s = 1:2
            [g, rr] = spoofInject(sp(k).prLog.(prN{s}).geomRange, sp(k).prLog.(prN{s}).geomRangeRate, ...
                sp(k).satPosLog.(spN{s}).position, sp(k).satPosLog.(spN{s}).velocity, refEcef, [0 0 0], [lat0 lon0], inj);
            sp(k).prLog.(prN{s}).geomRange = g;  sp(k).prLog.(prN{s}).geomRangeRate = rr;
        end
    end
end
ref = [lat0 lon0 h0];
%% A: a 5 m north step, 15 to 30 s (epochs 31 to 60)
scnA = struct('kind', 'step', 'tOn', 15, 'tOff', 30, 'size', 5, 'dirNed', [1 0 0]);
fprintf('--- A: 5 m north step on systems 1 and 2, 15 to 30 s, against the clean twin ---\n');
rA = spoofInjectCheck(spoofed(clean, scnA), clean, tKf, [], ref, scnA);
%% B: a 1 m/s east ramp, against a clean twin shifted by two epochs
scnB = struct('kind', 'ramp', 'tOn', 15, 'tOff', 30, 'size', 1.0, 'dirNed', [0 1 0]);
fprintf('--- B: 1 m/s east ramp, the clean twin without its first two epochs ---\n');
rB = spoofInjectCheck(spoofed(clean, scnB), clean(3:end), tKf, [], ref, scnB);
%% C: a 30 m clock step, no scenario given
scnC = struct('kind', 'clock', 'tOn', 15, 'tOff', 30, 'size', 30, 'dirNed', [1 0 0]);
fprintf('--- C: 30 m clock step, read without the scenario ---\n');
rC = spoofInjectCheck(spoofed(clean, scnC), clean, tKf, [], ref);
%% D: the step of A against a clean twin that differs on six epochs outside the interval
cleanD = clean;
for k = 70:75, cleanD(k).prLog.gpsRange.geomRange(1) = cleanD(k).prLog.gpsRange.geomRange(1) + 0.5; end
fprintf('--- D: the step of A against a clean twin that differs on epochs 70 to 75 ---\n');
rD = spoofInjectCheck(spoofed(clean, scnA), cleanD, tKf, [], ref, scnA);

okA = all(abs(rA.tFirst(1:2) - tKf(31)) < 1e-9) && all(abs(rA.tLast(1:2) - tKf(60)) < 1e-9) && isnan(rA.tFirst(3)) ...
   && max(max(abs(rA.medFit(1:2, :) - repmat([5 0 0 0], 2, 1)))) < 1e-6 && max(abs(rA.fit(45, 1:4, 1) - [5 0 0 0])) < 1e-6 ...
   && rA.fit(45, 5, 1) < 1e-6 && rA.nChanged(45, 1) == 6 && rA.nChanged(45, 3) == 0 && rA.nUsed(45, 3) == 5 ...
   && isequal(rA.nMissing, [0 0 30]) && isequal(rA.nOutside, [0 0 0]) && rA.shift == 0 && max(abs(rA.fit(10, 1:4, 1))) < 1e-9 ...
   && max(max(abs(rA.medRate(1:2, :)))) < 5e-3 && isequal(size(rA.fit), [K 6 3]);
okB = rB.kClean(10) == 8 && isnan(rB.kClean(1)) && rB.shift == -2 ...
   && max(max(abs(rB.medRate(1:2, :) - repmat([0 1 0 0], 2, 1)))) < 1e-2 && max(abs(rB.fit(45, 1:4, 1) - [0 7 0 0])) < 1e-3 ...
   && max(abs(rB.fitRate(45, 1:4, 2) - [0 1 0 0])) < 1e-2 && isequal(rB.nMissing, [0 0 30]) && isequal(rB.nOutside, [0 0 0]);
okC = max(abs(rC.medFit(1, :) - [0 0 0 30])) < 1e-6 && max(abs(rC.medFit(2, :) - [0 0 0 30])) < 1e-6 && max(abs(rC.medRate(1, :))) < 1e-9 ...
   && all(isnan(rC.intended(:))) && all(isnan(rC.nMissing)) && isnan(rC.tFirst(3));
okD = isequal(rD.nOutside, [6 0 0]) && abs(rD.tLast(1) - tKf(75)) < 1e-9 && abs(rD.tLast(2) - tKf(60)) < 1e-9 ...
   && max(abs(rD.medFit(1, :) - [5 0 0 0])) < 1e-6 && rD.nChanged(72, 1) == 1;
if okA && okB && okC && okD, fprintf('=== test_spoofInjectCheck PASS ===\n');
else, fprintf('okA %d okB %d okC %d okD %d\n=== test_spoofInjectCheck FAIL ===\n', okA, okB, okC, okD); end
end
