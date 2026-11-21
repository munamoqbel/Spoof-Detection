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
%   noise 0.003 m/s.
%   run: t, navVel, tTrue, trueVel, navAtt, fbVel ([K x 3], the hand-over per epoch). log: cell of records with
%   innovation, obsMatrix, numMeas, priorState, postState. aux: vTrue, vRowsTrue [N x 3], xPri, xPost, inc [K x 3],
%   jBefore [K x 1] (the nav sample before each epoch), rateTrue, course [N x 1], ARM, TS.
if nargin < 1, opt = struct(); end
o = struct('T', 400, 'seed', 21, 'corners', true, 'episodes', [200 205 0.35 0 0; 300 323 0 -0.4 3], 'gap', [350 360], ...
    'arm', [-0.18, 0.37], 'truthPoint', [0.25, 0.40], 'psi', 5e-3, 'accBias', 0.002, 'G', 0.1, 'KF', 0.5, 'delivery', 1, ...
    'noise', 0.003);
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
% hand-over is G times the posterior, taken off the state in full and off the nav with the given delivery
tKf = (0.5:0.5:o.T - 0.5)';  K = numel(tKf);
inGap = tKf >= o.gap(1) & tKf < o.gap(2);
vNav = zeros(N, 3);  e = [0.02, -0.01, 0];  vNav(1, :) = vTrue(1, :) + e;
x = zeros(1, 3);
xPri = zeros(K, 3);  xPost = zeros(K, 3);  inc = zeros(K, 3);  fbVel = zeros(K, 3);  jBefore = zeros(K, 1);
kNext = 1;
for j = 2:N
    dvT = vTrue(j, :) - vTrue(j - 1, :);
    e = e + o.psi * [-dvT(2), dvT(1), 0] + [o.accBias / fs, 0, 0];
    if kNext <= K && tN(j) >= tKf(kNext)
        jm = j - 1;  jBefore(kNext) = jm;
        xPri(kNext, :) = x;
        if ~inGap(kNext)
            dis = (vTrue(jm, :) + e - x) - vRowsTrue(jm, :);                % nav less prior, against the rows
            inc(kNext, :) = o.KF * dis;
        end
        xPost(kNext, :) = x + inc(kNext, :);
        if ~inGap(kNext)
            h = o.G * xPost(kNext, :);
            fbVel(kNext, :) = h;
            e = e - o.delivery * h;
            x = xPost(kNext, :) - h;
        else
            x = xPost(kNext, :);
        end
        kNext = kNext + 1;
    end
    vNav(j, :) = vTrue(j, :) + e + o.noise * randn(1, 3);
end

% the reference point ahead and left of the IMU: its velocity carries yaw rate times that offset
trueVelRef = vTrue;  trueVelRef(:, 1:2) = vTrue(:, 1:2) + rateTrue .* (o.truthPoint(1) * uAl + o.truthPoint(2) * uAc);
run = struct('t', tN, 'navVel', vNav, 'tTrue', tN, 'trueVel', trueVelRef, 'navAtt', [zeros(N, 2), course + o.psi], 'fbVel', fbVel);

% the records: 10 range rows and 10 rate rows on two constellations, the innovations predicted minus measured minus
% H times the prior, none in the gap
el = 0.2 + 1.1 * rand(10, 1);  az = 2 * pi * rand(10, 1);
u = [cos(el) .* cos(az), cos(el) .* sin(az), -sin(el)];
log = cell(K, 1);
for k = 1:K
    j = jBefore(k);
    H = zeros(mMax, n);  y = zeros(mMax, 1);
    for s = 1:10
        clk = 16 + 2 * double(s > 5);
        H(s, 1:3) = -u(s, :);          H(s, clk) = 1;      y(s) = 3 * randn;
        H(10 + s, 4:6) = -u(s, :);     H(10 + s, clk + 1) = 1;
        y(10 + s) = H(10 + s, 4:6) * (vNav(j, :) - xPri(k, :) - vRowsTrue(j, :))' + 0.05 * randn + 0.02 * double(s > 5);
    end
    numMeas = uint8(20);
    if inGap(k), numMeas = uint8(0); end
    log{k} = struct('innovation', y, 'obsMatrix', H, 'numMeas', numMeas, ...
        'priorState', [zeros(3, 1); xPri(k, :)'; zeros(n - 6, 1)], 'postState', [zeros(3, 1); xPost(k, :)'; zeros(n - 6, 1)]);
end
aux = struct('vTrue', vTrue, 'vRowsTrue', vRowsTrue, 'xPri', xPri, 'xPost', xPost, 'inc', inc, 'jBefore', jBefore, ...
    'rateTrue', rateTrue, 'course', course, 'ARM', o.arm, 'TS', o.truthPoint);
end
