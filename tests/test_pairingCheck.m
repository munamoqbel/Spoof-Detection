%% test_pairingCheck.m
% pairingCheck on the host's per-epoch log format (prLog / satPosLog / sysIdLog, one block of
% 15 slots per constellation) through pairingLogUnpack: a vehicle at 30 m/s, GPS and GLONASS in
% non-contiguous slots, Galileo switched off (all zeros), a receiver clock per system, 1 m code
% noise, one epoch without GNSS, satValid per slot on one system and scalar on the other. The
% pseudoranges are measured 0.09 s before the epoch they are tagged with. Checks: the unpacked
% arrays (slots, NaN padding, ids, positions), the best shift at -0.09 s with a residual at the
% noise level, the struct call equal to the array call, and the sysIdLog mismatch warning path.
% Run from the repo root (Octave with the shim).
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'tools'));
randn('seed', 3);  rand('seed', 3);

a = 6378137;  e2 = 6.69437999014e-3;
lla2ecef = @(p) [ (a ./ sqrt(1 - e2 * sin(p(:, 1)).^2) + p(:, 3)) .* cos(p(:, 1)) .* cos(p(:, 2)), ...
                  (a ./ sqrt(1 - e2 * sin(p(:, 1)).^2) + p(:, 3)) .* cos(p(:, 1)) .* sin(p(:, 2)), ...
                  (a ./ sqrt(1 - e2 * sin(p(:, 1)).^2) * (1 - e2) + p(:, 3)) .* sin(p(:, 1)) ];

% reference at 100 Hz: 30 m/s east with a slow north wobble
tRef = (0:0.01:300)';
lat0 = 0.90;  lon0 = 0.10;  h0 = 50;
vE = 30;  vN = 5 * sin(2 * pi * tRef / 60);
pN = cumsum([0; 0.5 * (vN(1:end-1) + vN(2:end)) * 0.01]);
pE = vE * tRef;
Rm = a * (1 - e2) / (1 - e2 * sin(lat0)^2)^1.5;  Rn = a / sqrt(1 - e2 * sin(lat0)^2);
posRef = [lat0 + pN / (Rm + h0), lon0 + pE / ((Rn + h0) * cos(lat0)), h0 * ones(size(tRef))];
ecefRef = lla2ecef(posRef);

% satellites: fixed directions above the horizon, 2.0e7 m away
S = 15;  nSat = [8 6 0];  usedSlots = {[1 2 3 5 6 9 11 14], [2 4 6 7 8 10], []};
up = ecefRef(1, :) / norm(ecefRef(1, :));
dirs = randn(S * 3, 3);  dirs = dirs - (dirs * up') * up;  dirs = dirs ./ repmat(sqrt(sum(dirs.^2, 2)), 1, 3);
elv = 0.2 + 0.6 * rand(S * 3, 1);
dirs = repmat(cos(elv), 1, 3) .* dirs + repmat(sin(elv), 1, 3) .* repmat(up, S * 3, 1);
satAll = repmat(ecefRef(1, :), S * 3, 1) + 2.0e7 * dirs;             % [45 x 3], slot (s-1)*15 + j
clk = [1500 1540 0];  tau = -0.09;

tKf = (10:0.5:290)';  K = numel(tKf);
blank = struct('corrValid', false, 'geomRange', zeros(S, 1), 'geomRangeRate', zeros(S, 1));
blankSt = struct('time', zeros(S, 1), 'satCbias', zeros(S, 1), 'satCdrift', zeros(S, 1), ...
    'position', zeros(S, 3), 'velocity', zeros(S, 3), 'LLA', zeros(S, 3), 'satValid', 1);
rec = struct('prLog', struct('gpsRange', blank, 'gloRange', blank, 'galRange', blank), ...
    'satPosLog', struct('navstarStatus', blankSt, 'glonassStatus', blankSt, 'galileoStatus', blankSt), 'sysIdLog', []);
loggedData = repmat(rec, K, 1);
prNames = {'gpsRange', 'gloRange', 'galRange'};  spNames = {'navstarStatus', 'glonassStatus', 'galileoStatus'};
for k = 2:K                                                            % epoch 1 stays without GNSS
    pMeas = interp1(tRef, ecefRef, tKf(k) + tau);                     % where the vehicle was when measured
    ids = [];
    for s = 1:3
        g = zeros(S, 1);  v = false(S, 1);  pos = zeros(S, 3);
        for j = usedSlots{s}
            sat = satAll((s - 1) * S + j, :);
            g(j) = norm(sat - pMeas) + clk(s) + 1.0 * randn;
            pos(j, :) = sat;  v(j) = true;
        end
        loggedData(k).prLog.(prNames{s}).geomRange = g;
        loggedData(k).prLog.(prNames{s}).corrValid = any(v);
        loggedData(k).satPosLog.(spNames{s}).position = pos;
        if s == 1, loggedData(k).satPosLog.(spNames{s}).satValid = v; end   % per slot on GPS, scalar elsewhere
        ids = [ids; s * ones(nSat(s), 1)]; %#ok<AGROW>
    end
    loggedData(k).sysIdLog = ids;
end

%% unpack
[pr, satPos, sysId] = pairingLogUnpack(loggedData);
okUnpack = isequal(size(pr), [K 3 * S]) && all(isnan(pr(1, :))) && nnz(isfinite(pr(2, :))) == 14 ...
    && isequal(find(isfinite(pr(2, :))), [usedSlots{1}, S + usedSlots{2}]) ...
    && isequal(sysId(2, isfinite(pr(2, :))), [ones(1, 8) 2 * ones(1, 6)]) ...
    && isequal(squeeze(satPos(2, 1, :))', satAll(1, :)) && isequal(squeeze(satPos(2, S + 4, :))', satAll(S + 4, :)) ...
    && all(all(isnan(pr(:, 2 * S + 1:end))));

%% the shift, struct call against array call
shifts = -0.3:0.01:0.3;
[costS, shS] = pairingCheck(tKf, loggedData, [], tRef, posRef, shifts);
[costA, shA] = pairingCheck(tKf, pr, satPos, tRef, posRef, shifts, sysId);
[cS, iS] = min(costS);
okShift = abs(shS(iS) - tau) < 0.011 && cS < 1.5 && isequal(costS, costA) && isequal(shS, shA);

%% sysIdLog that disagrees with the slots: a warning, the arrays unchanged
bad = loggedData;
bad(5).sysIdLog = [ones(7, 1); 2 * ones(7, 1)];
[prB, satB, idB] = pairingLogUnpack(bad);
okWarn = isequaln(prB, pr) && isequaln(satB, satPos) && isequaln(idB, sysId);

fprintf('unpack %d, shift %+.2f s (%.2f m) %d, warning path %d\n', okUnpack, shS(iS), cS, okShift, okWarn);
if okUnpack && okShift && okWarn
    fprintf('=== test_pairingCheck PASS ===\n');
else
    fprintf('=== test_pairingCheck FAIL ===\n');
end
