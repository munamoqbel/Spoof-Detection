function test_spoofProfile()
%% test_spoofProfile.m
% tools/spoofProfile on synthetic spfMeas records: a 300 m offset on four of eight range rows between two epochs
% must be found with its onset, end, size and system; a majority spoofed group (five of eight) must still be the
% offset group; a clean log must report nothing. Then a run whose spoofed epochs alternate between GPS-only rows
% (all spoofed, no split possible) and all systems: the two-group reading finds no run of 3 epochs, the per-system
% course shows the row set, and the reference-based profile from a consistent pairing log (clock ramp, inter-system
% biases, 1 m code noise) reads the 600 m injection whole with its onset and end, the clock interpolated across the
% GPS-only epochs. Then a position spoof (G): GPS ranges that describe a point 800 m north, 300 m west and 100 m
% below the reference with a 50 m clock offset, one garbage slot-epoch (range and position of different
% satellites), and a gate log for the timeline: the fitted position and clock must come back, in the innovations
% and against the reference, with the bad slot dropped; the GPS rate rows carry a 0.5, -0.3, 0.1 m/s velocity and a 0.2 m/s
% drift on top of a 0.3 m/s receiver clock drift, and the Doppler reading must return them from the innovations and from the
% log's rates against the satellite velocities. Run from anywhere: octave --no-gui --quiet --eval test_spoofProfile
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'tools'));
randn('seed', 9);  rand('seed', 9);
K = 200;  n = 22;  mMax = 51;  tKf = 0.5 * (0:K-1)';
function [log, sysId] = makeLog(K, n, mMax, sysRows, offRows, kOn, kEnd, size, firstOnly)
    % sysRows: range rows per system; firstOnly(k): the epoch holds the first system's rows only
    log = cell(K, 1);  sysId = zeros(K, mMax);
    nR = sum(sysRows);  ids = [];
    for s = 1:numel(sysRows), ids = [ids, s * ones(1, sysRows(s))]; end %#ok<AGROW>
    for k = 1:K
        keep = true(1, nR);
        if firstOnly(k), keep = ids == 1; end
        el = 0.2 + 0.8 * rand(nR, 1);  az = 2 * pi * rand(nR, 1);
        u  = [cos(el) .* cos(az), cos(el) .* sin(az), -sin(el)];
        H = zeros(mMax, n);  y = zeros(mMax, 1);
        for j = 1:nR, H(j, 1:3) = -u(j, :);  H(j, 16) = 1;  y(j) = 2 * randn; end                 % ranges
        for j = 1:nR, H(nR + j, 4:6) = -u(j, :);  H(nR + j, 17) = 1;  y(nR + j) = 0.05 * randn; end  % range rates
        H(2 * nR + 1, 3) = 1;  y(2 * nR + 1) = 1.5 * randn;                                       % pressure
        if k >= kOn && k <= kEnd
            if numel(size) == 3, y(offRows) = y(offRows) + H(offRows, 1:3) * size(:);   % a position spoof: per-row offset along the line of sight
            else, y(offRows) = y(offRows) + size; end
        end
        rows = [find(keep), nR + find(keep), 2 * nR + 1];                                        % compacted as the host builds them
        m = numel(rows);
        Hk = zeros(mMax, n);  yk = zeros(mMax, 1);  Hk(1:m, :) = H(rows, :);  yk(1:m) = y(rows);
        sysId(k, 1:m) = [ids(keep), ids(keep), 0];
        log{k} = struct('innovation', yk, 'innovationCov', eye(mMax), 'obsMatrix', Hk, 'measNoiseCov', eye(mMax), ...
            'numMeas', m, 'priorState', zeros(n, 1), 'postState', zeros(n, 1), 'postCov', eye(n));
    end
end
none = false(K, 1);
[logA, sysA] = makeLog(K, n, mMax, [4 4], 1:4, 80, 140, 300, none);
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
    loggedM(k).sysIdLog = [1 1 1 1 1 4 4 4 4]';                                                  % and GLONASS coded as 4
end
fprintf('--- A through the pairing log ---\n');
pA3 = spoofProfile(logA, tKf, loggedA);
fprintf('--- A through a pairing log with one used slot too many ---\n');
pA4 = spoofProfile(logA, tKf, loggedM);
[logB, sysB] = makeLog(K, n, mMax, [4 4], 1:5, 50, 60, -40, none);
fprintf('--- B: -40 m on five of eight rows (a majority), epochs 50 to 60 ---\n');
pB = spoofProfile(logB, tKf, sysB);
[logC, ~] = makeLog(K, n, mMax, [4 4], 1:4, 1, 0, 0, none);
fprintf('--- C: clean ---\n');
pC = spoofProfile(logC, tKf);

%% F: 600 m on GPS, epochs 80 to 140, the odd epochs of the stretch with GPS rows only; a consistent pairing log
% with a stationary reference, a clock ramp of 30 m/s, inter-system biases and 1 m code noise
firstOnly = false(K, 1);  firstOnly(81:2:139) = true;
[logF, sysF] = makeLog(K, n, mMax, [4 4 4], 1:4, 80, 140, 600, firstOnly);
for k = 1:8, logF{k}.innovation(1:4) = logF{k}.innovation(1:4) - 12; end      % an early transient: its own stretch
lat0 = 0.90;  lon0 = 0.10;  h0 = 50;
a = 6378137;  e2 = 6.69437999014e-3;
Nr = a / sqrt(1 - e2 * sin(lat0)^2);
refEcef = [(Nr + h0) * cos(lat0) * cos(lon0), (Nr + h0) * cos(lat0) * sin(lon0), (Nr * (1 - e2) + h0) * sin(lat0)];
S = 15;  usedF = {[1 2 3 5], [2 4 6 7], [3 8 10 12]};  isbF = [0 40 -25];
up = refEcef / norm(refEcef);
dirs = randn(S * 3, 3);  dirs = dirs - (dirs * up') * up;  dirs = dirs ./ repmat(sqrt(sum(dirs.^2, 2)), 1, 3);
elv = 0.2 + 0.6 * rand(S * 3, 1);
dirs = repmat(cos(elv), 1, 3) .* dirs + repmat(sin(elv), 1, 3) .* repmat(up, S * 3, 1);
satAll = repmat(refEcef, S * 3, 1) + 2.0e7 * dirs;
blank = struct('corrValid', true, 'geomRange', zeros(S, 1), 'geomRangeRate', zeros(S, 1));
blankSt = struct('position', zeros(S, 3), 'velocity', zeros(S, 3), 'LLA', zeros(S, 3), 'satValid', 1);
recF = struct('prLog', struct('gpsRange', blank, 'gloRange', blank, 'galRange', blank), ...
    'satPosLog', struct('navstarStatus', blankSt, 'glonassStatus', blankSt, 'galileoStatus', blankSt), 'sysIdLog', []);
loggedF = repmat(recF, K, 1);
prN = {'gpsRange', 'gloRange', 'galRange'};  spN = {'navstarStatus', 'glonassStatus', 'galileoStatus'};
for k = 1:K
    clk = 1500 + 30 * tKf(k);
    ids = [];
    for s = 1:3
        g = zeros(S, 1);  pos = zeros(S, 3);
        if s == 1 || ~firstOnly(k)
            for j = usedF{s}
                sat = satAll((s - 1) * S + j, :);
                g(j) = norm(sat - refEcef) + clk + isbF(s) + 1.0 * randn;
                if s == 1 && k >= 80 && k <= 140, g(j) = g(j) + 600; end
                pos(j, :) = sat;
            end
            ids = [ids; s * ones(numel(usedF{s}), 1)]; %#ok<AGROW>
        end
        loggedF(k).prLog.(prN{s}).geomRange = g;
        loggedF(k).satPosLog.(spN{s}).position = pos;
    end
    loggedF(k).sysIdLog = ids;
end
fprintf('--- F: 600 m on GPS, epochs 80 to 140, GPS rows only on the odd epochs; against the reference position ---\n');
pF = spoofProfile(logF, tKf, loggedF, [], [lat0 lon0 h0]);

%% G: a position spoof on GPS, epochs 60 to 120, every epoch with all systems; a garbage slot at epoch 70; a gate log
dP = [800 -300 100];  dClk = 50;  dV = [0.5 -0.3 0.1];  dDrift = 0.2;  drift0 = 0.3;
[logG, ~] = makeLog(K, n, mMax, [4 4 4], 1:4, 60, 120, dP, none);
for k = 60:120
    logG{k}.innovation(1:4) = logG{k}.innovation(1:4) + dClk;
    Hk = logG{k}.obsMatrix;                                                   % the GPS rate rows sit at 13:16
    logG{k}.innovation(13:16) = logG{k}.innovation(13:16) + Hk(13:16, 4:6) * dV' + dDrift;
end
C = [-sin(lat0) * cos(lon0), -sin(lat0) * sin(lon0), cos(lat0); -sin(lon0), cos(lon0), 0; -cos(lat0) * cos(lon0), -cos(lat0) * sin(lon0), -sin(lat0)];
spoofEcef = refEcef + (C' * dP')';
satVelAll = randn(S * 3, 3);  satVelAll = 3000 * satVelAll ./ repmat(sqrt(sum(satVelAll.^2, 2)), 1, 3);   % 3 km/s
spoofVelEcef = (C' * dV')';
loggedG = loggedF;
for k = 1:K
    clk = 1500 + 30 * tKf(k);
    ids = [];
    for s = 1:3
        g = zeros(S, 1);  pos = zeros(S, 3);  rr = zeros(S, 1);  vel = zeros(S, 3);
        for j = usedF{s}
            sat = satAll((s - 1) * S + j, :);  vSat = satVelAll((s - 1) * S + j, :);  u = dirs((s - 1) * S + j, :);
            if s == 1 && k >= 60 && k <= 120
                g(j) = norm(sat - spoofEcef) + clk + dClk + isbF(s) + 1.0 * randn;
                rr(j) = u * (vSat - spoofVelEcef)' + drift0 + dDrift + 0.05 * randn;
            else
                g(j) = norm(sat - refEcef) + clk + isbF(s) + 1.0 * randn;
                rr(j) = u * vSat' + drift0 + 0.05 * randn;
            end
            pos(j, :) = sat;  vel(j, :) = vSat;
        end
        if s == 2 && k == 70, pos(2, :) = 1.1 * satAll(S + 9, :); end  % the range of slot 2 with another satellite's position
        loggedG(k).prLog.(prN{s}).geomRange = g;
        loggedG(k).prLog.(prN{s}).geomRangeRate = rr;
        loggedG(k).satPosLog.(spN{s}).position = pos;
        loggedG(k).satPosLog.(spN{s}).velocity = vel;
        if s == 2 && k >= 65 && k <= 68                                  % a range the host did not use
            v = false(S, 1);  v(usedF{2}) = true;  v(4) = false;
            loggedG(k).satPosLog.(spN{s}).satValid = v;
        end
        if s == 2 && k >= 75 && k <= 77, loggedG(k).prLog.(prN{s}).corrValid = false; end   % corrections flagged invalid
        ids = [ids; s * ones(numel(usedF{s}), 1)]; %#ok<AGROW>
    end
    loggedG(k).sysIdLog = ids;
end
spfLogG = struct('count', K, 'mode', ones(K, 1), 'evLatch', 62, 'evProbStart', 118, 'evVeto', [], 'evHandback', 124);
spfLogG.mode(62:117) = 2;  spfLogG.mode(118:123) = 3;
fprintf('--- G: a position spoof on GPS (800 N, -300 E, 100 D, clock 50 m), epochs 60 to 120; one garbage slot; gate log ---\n');
pG = spoofProfile(logG, tKf, loggedG, [], [lat0 lon0 h0], spfLogG);
[~, ~, ~, nRejG, corrG] = pairingLogUnpack(loggedG);

okA = abs(pA.onset - tKf(80)) < 1e-9 && abs(pA.tEnd - tKf(140)) < 1e-9 && abs(pA.sizeOnset - 300) < 5 && pA.sysOff == 1 ...
   && pA.nOff(100) == 4 && abs(pA.perSys(100, 1) - 300) < 5 && abs(pA.perSys(100, 2)) < 5 ...
   && abs(pA2.onset - tKf(80)) < 1e-9 && isnan(pA2.sysOff) && abs(pA2.sizeOnset - 300) < 5 ...
   && pA3.sysOff == 1 && pA3.nSysMismatch == 0 && abs(pA3.perSys(100, 1) - 300) < 5 && abs(pA3.perSys(100, 2)) < 5 ...
   && pA4.nSysMismatch == K;
okB = abs(pB.onset - tKf(50)) < 1e-9 && abs(pB.tEnd - tKf(60)) < 1e-9 && abs(pB.sizeOnset + 40) < 5 && pB.nOff(55) == 5 && pB.sysOff == 1;
okC = isnan(pC.onset) && ~any(pC.present);
okF = isnan(pF.onset) && pF.nSysMismatch == 0 && abs(pF.perSys(101, 1) - 600) < 5 && pF.nPerSys(101, 2) == 0 ...
   && abs(pF.injOnset - tKf(80)) < 1e-9 && abs(pF.injEnd - tKf(140)) < 1e-9 && pF.injSys == 1 && abs(pF.injSize - 600) < 5 ...
   && abs(pF.inj(100, 2)) < 5 && abs(pF.inj(100, 3)) < 5 && abs(pF.inj(101, 1) - 600) < 5 && pF.nSlots(101, 2) == 0 ...
   && abs(pF.clockGap(101) - 0.5) < 1e-9 && pF.clockGap(100) == 0 && max(max(abs(pF.inj(1:70, :)))) < 10 ...
   && abs(pF.isb(2) - 40) < 2 && abs(pF.isb(3) + 25) < 2;
fG = median(pG.fitRef(60:120, :, 1), 1);  fI = median(pG.fitInn(60:120, :, 1), 1);   % 4 rows per fit: single epochs scatter in D and clock
medFin = @(X) median(X(isfinite(X(:, 1)), :), 1);                                    % epochs without a fit are NaN rows
fRI = medFin(pG.fitRateInn(60:120, :, 1));  fRR = medFin(pG.fitRateRef(60:120, :, 1));  fR2 = medFin(pG.fitRateRef(60:120, :, 2));
okR = max(abs(fRI(1:2) - dV(1:2))) < 0.1 && abs(fRI(3) - dV(3)) < 0.2 && abs(fRI(4) - dDrift) < 0.2 ...
   && max(abs(fRR(1:2) - dV(1:2))) < 0.1 && abs(fRR(3) - dV(3)) < 0.2 && abs(fRR(4) - dDrift) < 0.2 ...
   && max(abs(fR2(1:2))) < 0.1 && abs(median(pG.rateDrift(isfinite(pG.rateDrift))) - drift0) < 0.05 ...
   && max(abs(pG.rateAfter(1, 1:2))) < 0.1 && max(abs(pG.rateBefore(1, 1:2))) < 0.1;
okG = pG.injSys == 1 && abs(pG.injOnset - tKf(60)) < 1e-9 && abs(pG.injEnd - tKf(120)) < 1e-9 && pG.nBadSlot == 1 ...
   && pG.nSlots(70, 2) == 3 && max(abs(fG(1:2) - dP(1:2))) < 5 && abs(fG(3) - dP(3)) < 10 && abs(fG(4) - dClk) < 10 && fG(5) < 3 && fG(6) == 4 ...
   && pG.offSysInn == 1 && max(abs(fI(1:2) - dP(1:2))) < 5 && abs(fI(3) - dP(3)) < 10 && abs(fI(4) - dClk) < 10 && fI(6) == 4 ...
   && max(abs(pG.fitRef(90, 1:2, 1) - dP(1:2))) < 15 && abs(pG.fitRef(90, 3, 1) - dP(3)) < 60 ...
   && max(abs(pG.fitRef(30, 1:2, 1))) < 15 && pG.fitRef(30, 6, 1) == 4 && abs(pG.isb(2) - 40) < 2 && pF.injSys == 1 ...
   && max(abs(pG.fixBefore(1, 1:2))) < 10 && max(abs(pG.fixAfter(1, 1:2))) < 10 && max(abs(pG.fixAfter(2, 1:2))) < 10 ...
   && nnz(corrG(:, 2)) == 3 && nnz(corrG(:, [1 3])) == 0 && size(pG.fitRef, 3) == 3 && nnz(pG.badEpochs) == 1 && pG.badEpochs(70) ...
   && pG.fitSign == 1 && pG.nRecNoRows == 0 && pG.nLogNoSlots == 0 && pF.nRecNoRows == 0 ...
   && size(pF.stretchesInn, 1) == 2 && pF.stretchesInn(1, 1) == 1 && pF.stretchesInn(2, 1) == 80 && pF.offSysInn == 1 ...
   && size(pF.stretchesRef, 1) == 1 && nnz(nRejG(:, 2)) == 4 && all(nRejG(65:68, 2) == 1) && nnz(nRejG(:, [1 3])) == 0;
if okA && okB && okC && okF && okG && okR, fprintf('=== test_spoofProfile PASS ===\n');
else, fprintf('okA %d okB %d okC %d okF %d okG %d okR %d\n=== test_spoofProfile FAIL ===\n', okA, okB, okC, okF, okG, okR); end
end
