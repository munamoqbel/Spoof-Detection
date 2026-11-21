function test_dopplerScreen()
%% test_dopplerScreen.m
% SPF_dopplerScreen, the host-side Doppler episode screen, called epoch by epoch on a synthetic run and held
% against tools/dopplerEpisodeCheck on the same records with the host's correction logged: the two must flag the
% same epochs. The run: east at 20 m/s with a braking and a gentle turn under a 5 mrad heading error, the nav
% taking a tenth of the state per epoch; the receiver's Doppler wrong three times, 0.35 m/s for 5 s, a 0.4 m/s ramp
% for 20 s through a second braking, and 0.5 m/s for 35 s, longer than the hold; a gap without rows. Checks: the
% flags identical to the tool's; the first two episodes caught within two epochs and covered; no flag in the
% honest manoeuvres; the gap undecided with the coast kept; the long episode held through its hold limit and the
% rows' return afterwards accepted within the grown tolerance; the state and parameter constructors as documented.
% Run from the repo root (MATLAB, or Octave with tools/octave_shim).
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'tools'));
randn('seed', 5);  rand('seed', 5);

fs = 100;  T = 430;  tN = (0:1/fs:T)';  N = numel(tN);  n = CST_gnssHybrid.NO_STATES;  mMax = double(CST_spfParam.MAX_MEAS);
speed = 20 * ones(N, 1);
i1 = tN >= 100 & tN < 110;   speed(i1) = 20 - 1.2 * (tN(i1) - 100);   speed(tN >= 110 & tN < 130) = 8;
i2 = tN >= 130 & tN < 140;   speed(i2) = 8 + 1.2 * (tN(i2) - 130);
i5 = tN >= 298 & tN < 310;   speed(i5) = 20 - 1.17 * (tN(i5) - 298);  speed(tN >= 310 & tN < 330) = 6;
i6 = tN >= 330 & tN < 345;   speed(i6) = 6 + 0.93 * (tN(i6) - 330);
trap = @(t, t0, t1) min(min(max((t - t0) / 0.5, 0), 1), min(max((t1 - t) / 0.5, 0), 1));
rateTrue = -(10 * pi / 180) * trap(tN, 150, 159) + (10 * pi / 180) * trap(tN, 250, 259);
course = pi / 2 + cumsum(rateTrue) / fs;
vTrue = [speed .* cos(course), speed .* sin(course), zeros(N, 1)];
rowsErr = zeros(N, 3);
rowsErr(tN >= 200 & tN < 205, 1) = 0.35;
ramp = @(t, t0) min(max((t - t0) / 3, 0), 1);
rowsErr(:, 2) = -0.4 * (ramp(tN, 300) - ramp(tN, 320));
rowsErr(tN >= 360 & tN < 395, 1) = rowsErr(tN >= 360 & tN < 395, 1) + 0.5;
vRowsTrue = vTrue + rowsErr;
% the nav: a tenth of the state fed back at each update, the hand-over kept
psi = 5e-3;  G = 0.1;  KF = 0.5;
tKf = (0.5:0.5:T - 0.5)';  K = numel(tKf);
vNav = zeros(N, 3);  e = [0.02, -0.01, 0];  vNav(1, :) = vTrue(1, :) + e;
xPlus = zeros(K, 3);  kNext = 1;
for j = 2:N
    dvT = vTrue(j, :) - vTrue(j - 1, :);
    e = e + psi * [-dvT(2), dvT(1), 0];
    if kNext <= K && tN(j) >= tKf(kNext)
        jm = j - 1;
        dis = (vTrue(jm, :) + e) - vRowsTrue(jm, :);
        xPlus(kNext, :) = KF * dis;
        e = e - G * xPlus(kNext, :);
        kNext = kNext + 1;
    end
    vNav(j, :) = vTrue(j, :) + e + 0.003 * randn(1, 3);
end
fbVel = G * xPlus;
% the records: 10 range and 10 rate rows on two constellations, none in the gap
el = 0.2 + 1.1 * rand(10, 1);  az = 2 * pi * rand(10, 1);
u = [cos(el) .* cos(az), cos(el) .* sin(az), -sin(el)];
log = cell(K, 1);  jBefore = zeros(K, 1);
for k = 1:K
    j = find(tN < tKf(k), 1, 'last');  jBefore(k) = j;
    H = zeros(mMax, n);  y = zeros(mMax, 1);
    for s = 1:10
        clk = 16 + 2 * double(s > 5);
        H(s, 1:3) = -u(s, :);          H(s, clk) = 1;      y(s) = 3 * randn;
        H(10 + s, 4:6) = -u(s, :);     H(10 + s, clk + 1) = 1;
        y(10 + s) = H(10 + s, 4:6) * (vNav(j, :) - vRowsTrue(j, :))' + 0.05 * randn + 0.02 * double(s > 5);
    end
    numMeas = uint8(20);
    if tKf(k) >= 350 && tKf(k) < 356, numMeas = uint8(0); end
    log{k} = struct('innovation', y, 'obsMatrix', H, 'numMeas', numMeas, 'postState', [zeros(3, 1); xPlus(k, :)'; zeros(n - 6, 1)]);
end

%% the host's screen, epoch by epoch
prm = STRUCT_SPF.dopplerScreenParam;
screen = STRUCT_SPF.dopplerScreenState;
flagHost = false(K, 1);  validHost = false(K, 1);  holdHost = zeros(K, 1);
for k = 1:K
    velNavPre = vNav(jBefore(k), :)';
    if k > 1, velCorr = fbVel(k - 1, :)'; else, velCorr = zeros(3, 1); end
    [screen, drop, info] = SPF_dopplerScreen(screen, velNavPre, velCorr, log{k}.obsMatrix, log{k}.innovation, ...
        log{k}.numMeas, [17 19 22], prm);
    flagHost(k) = drop;  validHost(k) = info.valid;  holdHost(k) = info.hold;
    if k == 100, rateRows100 = info.rateRows; end
end

%% the offline tool on the same records, the correction logged and delivered in the next interval
run = struct('t', tN, 'navVel', vNav, 'tTrue', tN, 'trueVel', vTrue, 'navAtt', [zeros(N, 2), course + psi], 'fbVel', fbVel);
txt = evalc('res = dopplerEpisodeCheck(log, tKf, run, struct(''fbLagWeights'', [0 1]));');
fprintf('%s', txt);

ok = true;  chk = struct();
chk.same     = isequal(flagHost, res.flag);
ep1 = tKf >= 200 & tKf < 205;  ep2 = tKf >= 302 & tKf < 320;  ep3 = tKf >= 360 & tKf < 395;
chk.caught1  = find(flagHost & ep1, 1) - find(ep1, 1) <= 2 && mean(flagHost(ep1)) >= 0.7;
chk.caught2  = find(flagHost & ep2, 1) - find(ep2, 1) <= 2 && mean(flagHost(ep2)) >= 0.8;
chk.quiet    = ~any(flagHost(tKf >= 95 & tKf < 170)) && ~any(flagHost(tKf >= 245 & tKf < 270));
chk.gap      = ~any(validHost(tKf >= 350 & tKf < 356)) && ~any(flagHost(tKf >= 350 & tKf < 356)) && ...
    validHost(find(tKf >= 356, 1));
% the long episode: flagged from its start, the hold limit reached at 20 s, the episode reopened at once since the
% rows are still beyond the grown tolerance, and the rows' return accepted within it, so the return costs little
kHold = find(holdHost >= prm.HOLD - 1 & tKf > 360, 1);
chk.holdLimit = ~isempty(kHold) && tKf(kHold) > 378 && tKf(kHold) < 382 && mean(flagHost(tKf >= 361 & tKf < 395)) >= 0.9;
chk.return   = nnz(flagHost(tKf >= 396 & tKf < 430)) <= 8;
chk.structs  = islogical(screen.anchored) && isfield(screen, 'speedSince') && prm.FLOOR == 0.20 && prm.ROWS_SIGN == -1 && ...
    prm.RELEASE_GROW == 10;
chk.rateRows = nnz(rateRows100) == 10 && all(find(rateRows100) == (11:20)');
names = fieldnames(chk);
for i = 1:numel(names)
    if ~chk.(names{i}), fprintf('   FAILED sub-check: %s\n', names{i}); ok = false; end
end
fprintf('host flags %d, tool flags %d, differing epochs %d; long episode hold limit at %.1f s\n', nnz(flagHost), nnz(res.flag), ...
    nnz(flagHost ~= res.flag), tKf(max(kHold, 1)));
if ok
    fprintf('=== test_dopplerScreen PASS ===\n');
else
    fprintf('=== test_dopplerScreen FAIL ===\n');
end
end
