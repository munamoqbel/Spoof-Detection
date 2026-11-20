function test_dopplerEpisodeCheck()
%% test_dopplerEpisodeCheck.m
% A synthetic hour-like run: a vehicle heading east at 20 m/s with a braking, a left turn and a right turn under a
% 5 mrad heading error in the nav, and a receiver whose Doppler rows go wrong twice, 0.3 m/s for 5 s on a straight
% leg and a 0.4 m/s ramp for 20 s through a second braking. The nav takes a tenth of the state per epoch; the records
% carry the rate-row innovations the host would have built. Checks that tools/dopplerEpisodeCheck catches both
% episodes within 4 s with most of their length covered, raises no false flag on the honest manoeuvres (the flags
% that outlast an episode by the release count are booked as release tails, not false flags), fits the feedback
% share, survives a gap without rows, and reports a miss when the floor is set above the episode.
% Run from the repo root (MATLAB, or Octave with tools/octave_shim).
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'tools'));
randn('seed', 21);  rand('seed', 21);

fs = 100;  T = 400;  tN = (0:1/fs:T)';  N = numel(tN);  n = 22;  mMax = 51;
% truth: speed and course
speed = 20 * ones(N, 1);  course = (pi / 2) * ones(N, 1);
i1 = tN >= 100 & tN < 110;   speed(i1) = 20 - 1.2 * (tN(i1) - 100);   speed(tN >= 110 & tN < 130) = 8;
i2 = tN >= 130 & tN < 140;   speed(i2) = 8 + 1.2 * (tN(i2) - 130);
i3 = tN >= 150 & tN < 159;   course(i3) = pi / 2 - (10 * pi / 180) * (tN(i3) - 150);   course(tN >= 159 & tN < 250) = pi / 2 - pi / 2;
i4 = tN >= 250 & tN < 259;   course(i4) = (10 * pi / 180) * (tN(i4) - 250);            course(tN >= 259) = pi / 2;
i5 = tN >= 298 & tN < 310;   speed(i5) = 20 - 1.17 * (tN(i5) - 298);  speed(tN >= 310 & tN < 330) = 6;
i6 = tN >= 330 & tN < 345;   speed(i6) = 6 + 0.93 * (tN(i6) - 330);
vTrue = [speed .* cos(course), speed .* sin(course), zeros(N, 1)];
% the receiver's episodes: 0.3 m/s north for 5 s at 200 s; a ramp to 0.4 m/s along east over 3 s from 300 s, back over 3 s from 320 s
rowsErr = zeros(N, 3);
rowsErr(tN >= 200 & tN < 205, 1) = 0.3;
ramp = @(t, t0) min(max((t - t0) / 3, 0), 1);
rowsErr(:, 2) = -0.4 * (ramp(tN, 300) - ramp(tN, 320));
vRowsTrue = vTrue + rowsErr;
% the nav: the INS increment carries a 5 mrad heading error, the feedback takes a tenth of the state at each epoch
psi = 5e-3;  G = 0.1;  KF = 0.5;
tKf = (0.5:0.5:T - 0.5)';  K = numel(tKf);
vNav = zeros(N, 3);  e = [0.02, -0.01, 0];  vNav(1, :) = vTrue(1, :) + e;
xPlus = zeros(K, 3);  kNext = 1;
for j = 2:N
    dvT = vTrue(j, :) - vTrue(j - 1, :);
    e = e + psi * [-dvT(2), dvT(1), 0];                              % the increment rotated by the heading error
    if kNext <= K && tN(j) >= tKf(kNext)
        % the update at this epoch: the state estimates half of the nav's disagreement with the rows, the nav takes a tenth
        jm = j - 1;
        dis = (vTrue(jm, :) + e) - vRowsTrue(jm, :);
        xPlus(kNext, :) = KF * dis;
        e = e - G * xPlus(kNext, :);
        kNext = kNext + 1;
    end
    vNav(j, :) = vTrue(j, :) + e;
end
run = struct('t', tN, 'navVel', vNav, 'tTrue', tN, 'trueVel', vTrue);
% the records: 10 range rows and 10 rate rows on two constellations, the innovations predicted minus measured
el = 0.2 + 1.1 * rand(10, 1);  az = 2 * pi * rand(10, 1);
u = [cos(el) .* cos(az), cos(el) .* sin(az), -sin(el)];
log = cell(K, 1);
for k = 1:K
    j = find(tN < tKf(k), 1, 'last');
    H = zeros(mMax, n);  y = zeros(mMax, 1);
    for s = 1:10
        clk = 16 + 2 * double(s > 5);
        H(s, 1:3) = -u(s, :);          H(s, clk) = 1;      y(s) = 3 * randn;
        H(10 + s, 4:6) = -u(s, :);     H(10 + s, clk + 1) = 1;
        y(10 + s) = H(10 + s, 4:6) * (vNav(j, :) - vRowsTrue(j, :))' + 0.05 * randn + 0.02 * double(s > 5);
    end
    numMeas = uint8(20);
    if tKf(k) >= 350 && tKf(k) < 360, numMeas = uint8(0); end        % a gap without rows
    log{k} = struct('innovation', y, 'obsMatrix', H, 'numMeas', numMeas, 'postState', [zeros(3, 1); xPlus(k, :)'; zeros(n - 6, 1)]);
end

ok = true;  chk = struct();
txt = evalc('res = dopplerEpisodeCheck(log, tKf, run);');
fprintf('%s', txt);
ep = res.episodes;
chk.two      = numel(ep) == 2;
chk.caught   = numel(ep) == 2 && all([ep.detected]) && all([ep.delay] <= 4.0);
chk.covered  = numel(ep) == 2 && ep(1).coverage > 0.6 && ep(2).coverage > 0.7;
chk.times    = numel(ep) == 2 && abs(ep(1).tStart - 200) < 1 && abs(ep(2).tStart - 301) < 2.5 && abs(ep(2).tEnd - 321) < 2.5;
chk.noFalse  = res.nFalse <= 2;
chk.tail     = res.nTail <= 2 * 4;                                    % the release latency: at most 'release' epochs per episode
chk.fbGain   = abs(res.fbGain - G) < 0.02;
chk.gap      = ~any(res.flag(tKf >= 350 & tKf < 360)) && all(isnan(res.onset(tKf >= 350.5 & tKf < 355)));
chk.verdict  = ~isempty(strfind(txt, 'every episode caught, no false flag'));
% the manoeuvres alone are honest: no flag in the first braking and the two turns
chk.quiet    = ~any(res.flag(tKf >= 95 & tKf < 170)) && ~any(res.flag(tKf >= 245 & tKf < 270));
% a floor above the small episode: it is missed and said so; the large one still caught
txt2 = evalc('res2 = dopplerEpisodeCheck(log, tKf, run, struct(''floor'', 0.5, ''quiet'', false));');
chk.missed   = ~isempty(strfind(txt2, 'MISSED')) && ~isempty(strfind(txt2, 'episodes missed')) && numel(res2.episodes) == 2 && ...
    ~res2.episodes(1).detected;
% a given feedback share is used as given
res3 = dopplerEpisodeCheck(log, tKf, run, struct('fbGain', 0.1, 'quiet', true));
chk.given    = res3.fbGain == 0.1 && numel(res3.episodes) == 2 && all([res3.episodes.detected]);
names = fieldnames(chk);
for i = 1:numel(names)
    if ~chk.(names{i}), fprintf('   FAILED sub-check: %s\n', names{i}); ok = false; end
end
if numel(ep) == 2
    fprintf('episodes: %.1f s caught after %.1f s (%.0f%%), %.1f s caught after %.1f s (%.0f%%); false flags %d, release tails %d; feedback share %.3f\n', ...
        ep(1).tEnd - ep(1).tStart, ep(1).delay, 100 * ep(1).coverage, ep(2).tEnd - ep(2).tStart, ep(2).delay, 100 * ep(2).coverage, ...
        res.nFalse, res.nTail, res.fbGain);
end
if ok
    fprintf('=== test_dopplerEpisodeCheck PASS ===\n');
else
    fprintf('=== test_dopplerEpisodeCheck FAIL ===\n');
end
end
