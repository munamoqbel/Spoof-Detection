function [run, log, tKf, aux] = dopplerSynthRun(opt)
%DOPPLERSYNTHRUN  A synthetic run for the Doppler screen tests, in the host's error-state form.
%   [run, log, tKf, aux] = dopplerSynthRun(opt)
%   A vehicle heading east at 20 m/s with two brakings, a left and a right turn and, by default, two sharp corners;
%   a receiver whose Doppler rows go wrong in given episodes and carry a lever-arm error; a nav whose increment
%   carries a heading error and an accelerometer bias; and an error-state filter in the host's form: at each epoch
%   the rows are measured against the nav velocity less the prior state x- (the innovation is predicted minus
%   measured minus H x-), the update's step is a share of that disagreement, the nav takes a tenth of the posterior
%   state and the state is reduced by the same, which the mechanisation may deliver only in part (opt.delivery).
%   No update and no hand-over on the epochs without rows. The records carry priorState and postState, the run the
%   nav velocity and attitude and the truth at a reference point offset from the IMU.
%   opt (defaults): T 400 s, seed 21, corners true, episodes [t0 t1 dN dE ramp] rows ([200 205 0.35 0 0; 300 323 0
%   -0.4 3]: the error dN, dE m/s between t0 and t1, ramped in and out over 'ramp' s), gap [350 360] s (no rows),
%   arm [-0.18 0.37] m per rad/s (the host's lever-arm error, along and across, in the rows), truthPoint [0.25 0.40]
%   m (the reference ahead and left of the IMU), psi 5e-3 rad, accBias 0.002 m/s2 north, G 0.1, KF 0.5, delivery 1,
%   noise 0.003 m/s, limiter Inf (m/s: above this horizontal disagreement the epoch is a limited one, the nav takes a
%   quarter of the posterior and the state keeps a quarter, LM 0.25), screen false (true: the host screen
%   SPF_dopplerScreen runs in the loop on each epoch's records before the update, and an epoch it flags leaves the
%   rows out of the update, as the host's R at its maximum does; its flags come back in aux). gap may have several
%   rows, one per stretch without rows; fewRows [] ([t0 t1] rows: stretches with 4 range and 4 rate rows on two
%   constellations, too few for the screen's fit against its five unknowns, so the epoch is not judged). The filter
%   also carries a position error state aided by the ranges every
%   epoch with rows (posNoise 0.5 m, Kp 0.3 on position, Kv 0.03 on velocity), which keeps working while the Doppler
%   rows are held out and pulls the velocity state back as the host's ranges do.
%   run: t, navVel, tTrue, trueVel, navAtt, fbVel ([K x 3], the hand-over per epoch). log: cell of records with
%   innovation, obsMatrix, numMeas, priorState, postState. aux: vTrue, vRowsTrue [N x 3], xPri, xPost, inc [K x 3],
%   jBefore [K x 1] (the nav sample before each epoch), rateTrue, course [N x 1], ARM, TS, limited [K x 1], and with
%   opt.screen hostFlag, hostValid, hostHold [K x 1] and rateRows100 (the screen's row mask at epoch 100).
if nargin < 1, opt = struct(); end
o = struct('T', 400, 'seed', 21, 'corners', true, 'episodes', [200 205 0.35 0 0; 300 323 0 -0.4 3], 'gap', [350 360], ...
    'arm', [-0.18, 0.37], 'truthPoint', [0.25, 0.40], 'psi', 5e-3, 'accBias', 0.002, 'G', 0.1, 'KF', 0.5, 'delivery', 1, ...
    'noise', 0.003, 'limiter', Inf, 'LM', 0.25, 'screen', false, 'posNoise', 0.5, 'Kp', 0.3, 'Kv', 0.03, 'fewRows', []);
f = fieldnames(opt);
for i = 1:numel(f), o.(f{i}) = opt.(f{i}); end
randn('seed', o.seed);  rand('seed', o.seed);
fs = 100;  tN = (0:1/fs:o.T)';  N = numel(tN);  n = CST_gnssHybrid.NO_STATES;  mMax = double(CST_spfParam.MAX_MEAS);

% truth: speed and course, the course integrated from a yaw-rate profile with every rate ramping over 0.5 s
speed = 20 * ones(N, 1);
i1 = tN >= 100 & tN < 110;   speed(i1) = 20 - 1.2 * (tN(i1) - 100);   speed(tN >= 110 & tN < 130) = 8;
i2 = tN >= 130 & tN < 140;   speed(i2) = 8 + 1.2 * (tN(i2) - 130);
i5 = tN >= 298 & tN < 310;   speed(i5) = 20 - 1.17 * (tN(i5) - 298);  speed(tN >= 310 & tN < 330) = 6;
i6 = tN >= 330 & tN < 345;   speed(i6) = 6 + 0.93 * (tN(i6) - 330);
trap = @(t, t0, t1) min(min(max((t - t0) / 0.5, 0), 1), min(max((t1 - t) / 0.5, 0), 1));
rateTrue = -(10 * pi / 180) * trap(tN, 150, 159) + (10 * pi / 180) * trap(tN, 250, 259);
if o.corners
    W = 45 * pi / 180;
    rateTrue = rateTrue + W * trap(tN, 118, 121) - W * trap(tN, 125, 128);
end
course = pi / 2 + cumsum(rateTrue) / fs;
vTrue = [speed .* cos(course), speed .* sin(course), zeros(N, 1)];
uAl = [cos(course), sin(course)];  uAc = [-sin(course), cos(course)];

% the receiver's episodes and the lever-arm error in the rows (nav minus yaw rate times arm, leverArmCheck's read)
rowsErr = zeros(N, 3);
for i = 1:size(o.episodes, 1)
    e = o.episodes(i, :);
    if e(5) > 0
        prof = min(max((tN - e(1)) / e(5), 0), 1) - min(max((tN - (e(2) - e(5))) / e(5), 0), 1);
    else
        prof = double(tN >= e(1) & tN < e(2));
    end
    rowsErr(:, 1) = rowsErr(:, 1) + e(3) * prof;
    rowsErr(:, 2) = rowsErr(:, 2) + e(4) * prof;
end
armErr = -rateTrue .* (o.arm(1) * uAl + o.arm(2) * uAc);
vRowsTrue = vTrue + rowsErr + [armErr, zeros(N, 1)];

% the nav and the filter: the nav error e grows with the heading error on the increment and the accelerometer bias;
% at each epoch the prior x- is what the rows are measured against, the step is KF times the disagreement, the
% hand-over is G times the posterior, taken off the state in full and off the nav with the given delivery. The
% records are built at the update instant from the nav sample before it; with opt.screen the host screen runs on
% them first and a flagged epoch leaves the rows out of the update, as the host's R at its maximum does
tKf = (0.5:0.5:o.T - 0.5)';  K = numel(tKf);
inGap = false(K, 1);
for g = 1:size(o.gap, 1)
    inGap = inGap | (tKf >= o.gap(g, 1) & tKf < o.gap(g, 2));
end
inFew = false(K, 1);
for g = 1:size(o.fewRows, 1)
    inFew = inFew | (tKf >= o.fewRows(g, 1) & tKf < o.fewRows(g, 2));
end
el = 0.2 + 1.1 * rand(10, 1);  az = 2 * pi * rand(10, 1);
u = [cos(el) .* cos(az), cos(el) .* sin(az), -sin(el)];
rowNoise = 0.05 * randn(K, 10);  rangeNoise = 3 * randn(K, 10);  posMeasNoise = o.posNoise * randn(K, 3);
vNav = zeros(N, 3);  e = [0.02, -0.01, 0];  vNav(1, :) = vTrue(1, :) + e;
x = zeros(1, 3);  pErr = zeros(1, 3);  xp = zeros(1, 3);                     % the nav position error and its state
xPri = zeros(K, 3);  xPost = zeros(K, 3);  inc = zeros(K, 3);  fbVel = zeros(K, 3);  jBefore = zeros(K, 1);
limited = false(K, 1);  hostFlag = false(K, 1);  hostValid = false(K, 1);  hostHold = zeros(K, 1);  rateRows100 = [];
log = cell(K, 1);
if o.screen
    prm = STRUCT_SPF.dopplerScreenParam;
    screen = STRUCT_SPF.dopplerScreenState;
end
kNext = 1;
for j = 2:N
    dvT = vTrue(j, :) - vTrue(j - 1, :);
    e = e + o.psi * [-dvT(2), dvT(1), 0] + [o.accBias / fs, 0, 0];
    pErr = pErr + e / fs;
    if kNext <= K && tN(j) >= tKf(kNext)
        k = kNext;  jm = j - 1;  jBefore(k) = jm;
        xPri(k, :) = x;  xpPri = xp;
        % the records of this epoch: 10 range rows and 10 rate rows on two constellations, the innovations
        % predicted minus measured minus H times the prior, none in a gap
        H = zeros(mMax, n);  y = zeros(mMax, 1);
        nSat = 10;  sats = 1:10;
        if inFew(k), nSat = 4;  sats = [1 2 6 7]; end                        % two per constellation: five unknowns, four rows
        for i = 1:nSat
            sIdx = sats(i);
            clk = 16 + 2 * double(sIdx > 5);
            H(i, 1:3) = -u(sIdx, :);          H(i, clk) = 1;      y(i) = rangeNoise(k, sIdx);
            H(nSat + i, 4:6) = -u(sIdx, :);   H(nSat + i, clk + 1) = 1;
            y(nSat + i) = H(nSat + i, 4:6) * (vNav(jm, :) - xPri(k, :) - vRowsTrue(jm, :))' + rowNoise(k, sIdx) + 0.02 * double(sIdx > 5);
        end
        numMeas = uint8(2 * nSat);
        if inGap(k), numMeas = uint8(0); end
        useRows = ~inGap(k);
        if o.screen
            [screen, drop, info] = SPF_dopplerScreen(screen, vNav(jm, :)', xPri(k, :)', H, y, numMeas, [17 19 22], prm);
            hostFlag(k) = drop;  hostValid(k) = info.valid;  hostHold(k) = info.hold;
            if k == 100, rateRows100 = info.rateRows; end
            if drop, useRows = false; end
        end
        dis = zeros(1, 3);
        if ~inGap(k)
            zp = (pErr - xpPri) + posMeasNoise(k, :);                        % the ranges: nav position less prior, against the receiver
            inc(k, :) = o.Kv * zp;                                            % the position aiding's share on the velocity
            xp = xpPri + o.Kp * zp;
        end
        if useRows
            dis = (vTrue(jm, :) + e - x) - vRowsTrue(jm, :);                % nav less prior, against the rows
            inc(k, :) = inc(k, :) + o.KF * dis;
        end
        xPost(k, :) = x + inc(k, :);
        if useRows && norm(dis(1:2)) > o.limiter
            h = [o.LM * xPost(k, 1:2), 0];                                    % a limited epoch: a quarter handed over,
            fbVel(k, :) = h;                                                  % the state keeps a quarter
            limited(k) = true;
            e = e - o.delivery * h;
            x = [o.LM * xPost(k, 1:2), xPost(k, 3)];
        elseif ~inGap(k)
            h = o.G * xPost(k, :);
            fbVel(k, :) = h;
            e = e - o.delivery * h;
            x = xPost(k, :) - h;
        else
            x = xPost(k, :);                                                  % no update: no hand-over either
        end
        if ~inGap(k)
            hp = o.G * xp;  pErr = pErr - hp;  xp = xp - hp;                  % the position hand-over, the same tenth
        end
        if o.screen, screen.incPrev = inc(k, :); end
        log{k} = struct('innovation', y, 'obsMatrix', H, 'numMeas', numMeas, ...
            'priorState', [xpPri'; xPri(k, :)'; zeros(n - 6, 1)], 'postState', [xp'; xPost(k, :)'; zeros(n - 6, 1)]);
        kNext = kNext + 1;
    end
    vNav(j, :) = vTrue(j, :) + e + o.noise * randn(1, 3);
end

% the reference point ahead and left of the IMU: its velocity carries yaw rate times that offset
trueVelRef = vTrue;  trueVelRef(:, 1:2) = vTrue(:, 1:2) + rateTrue .* (o.truthPoint(1) * uAl + o.truthPoint(2) * uAc);
run = struct('t', tN, 'navVel', vNav, 'tTrue', tN, 'trueVel', trueVelRef, 'navAtt', [zeros(N, 2), course + o.psi], 'fbVel', fbVel);

aux = struct('vTrue', vTrue, 'vRowsTrue', vRowsTrue, 'xPri', xPri, 'xPost', xPost, 'inc', inc, 'jBefore', jBefore, ...
    'rateTrue', rateTrue, 'course', course, 'ARM', o.arm, 'TS', o.truthPoint, 'limited', limited, ...
    'hostFlag', hostFlag, 'hostValid', hostValid, 'hostHold', hostHold, 'rateRows100', rateRows100);
end
