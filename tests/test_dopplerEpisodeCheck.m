function test_dopplerEpisodeCheck()
%% test_dopplerEpisodeCheck.m
% A synthetic hour-like run: a vehicle heading east at 20 m/s with a braking, a left turn and a right turn under a
% 5 mrad heading error in the nav, and a receiver whose Doppler rows go wrong twice, 0.35 m/s for 5 s on a straight
% leg and a 0.4 m/s ramp for 20 s through a second braking. The nav takes a tenth of the state per epoch; the records
% carry the rate-row innovations the host would have built. Checks that tools/dopplerEpisodeCheck catches both
% episodes within 4 s with most of their length covered, raises no false flag on the honest manoeuvres (the flags
% that outlast an episode by the release count are booked as release tails, not false flags), fits the feedback
% share at the lag the logged velocity shows it, survives a gap without rows, reports a miss when the floor is set
% above the episode, reads a lever-arm error back from the truth in two sharp corners and takes it out of the rows
% with opts.armSlope, takes the reference point's own offset out of the scoring with opts.truthSlope, finds a run
% whose clock is 3 s ahead and shifts it back, reads rows paired with the wrong pulse as a time offset, finds a
% state vector recorded in ENU by the share appearing only with north and east swapped, a nav that took a share
% of the update's increment rather than of the state, and the host's correction logged per epoch or per sample.
% Then the same vehicle in the host's error-state form (tests/dopplerSynthRun, the records carrying the prior): the
% total-state model is the default, nothing is fitted, both episodes are caught with no false flag, the share of the
% state left in the increment reads zero, the flags are the same without the truth, a mechanisation that delivers
% only three quarters of the hand-over is read as such, a feedback branch in which the state keeps a quarter while the
% nav takes a quarter is read from the state's loss between an update and the next prior, the hand-over form of the
% coast gives the same flags where the loss is the hand-over and the exact coast where it is not, and the fallback
% models still run when forced.
% Run from the repo root (MATLAB, or Octave with tools/octave_shim).
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'tools'));
randn('seed', 21);  rand('seed', 21);

fs = 100;  T = 400;  tN = (0:1/fs:T)';  N = numel(tN);  n = 22;  mMax = 51;
% truth: speed and course, the course integrated from a yaw-rate profile
speed = 20 * ones(N, 1);
i1 = tN >= 100 & tN < 110;   speed(i1) = 20 - 1.2 * (tN(i1) - 100);   speed(tN >= 110 & tN < 130) = 8;
i2 = tN >= 130 & tN < 140;   speed(i2) = 8 + 1.2 * (tN(i2) - 130);
i5 = tN >= 298 & tN < 310;   speed(i5) = 20 - 1.17 * (tN(i5) - 298);  speed(tN >= 310 & tN < 330) = 6;
i6 = tN >= 330 & tN < 345;   speed(i6) = 6 + 0.93 * (tN(i6) - 330);
% a left turn and a right turn at 10 deg/s for 9 s, and two sharp corners at 8 m/s, 45 deg/s, right then left,
% every yaw rate ramping over 0.5 s at each end
W = 45 * pi / 180;  trap = @(t, t0, t1) min(min(max((t - t0) / 0.5, 0), 1), min(max((t1 - t) / 0.5, 0), 1));
rateTrue = -(10 * pi / 180) * trap(tN, 150, 159) + (10 * pi / 180) * trap(tN, 250, 259) + W * trap(tN, 118, 121) - W * trap(tN, 125, 128);
course = pi / 2 + cumsum(rateTrue) / fs;
vTrue = [speed .* cos(course), speed .* sin(course), zeros(N, 1)];
% the receiver's episodes: 0.3 m/s north for 5 s at 200 s; a ramp to 0.4 m/s along east over 3 s from 300 s, back over 3 s from 320 s
rowsErr = zeros(N, 3);
rowsErr(tN >= 200 & tN < 205, 1) = 0.35;
ramp = @(t, t0) min(max((t - t0) / 3, 0), 1);
rowsErr(:, 2) = -0.4 * (ramp(tN, 300) - ramp(tN, 320));
% the host's lever arm wrong by leverArmCheck's slopes: the rows read nav minus (yaw rate x arm), along -0.18 m, across +0.37 m
ARM = [-0.18, 0.37];
uAl = [cos(course), sin(course)];  uAc = [-sin(course), cos(course)];
armErr = -rateTrue .* (ARM(1) * uAl + ARM(2) * uAc);
vRowsTrue = vTrue + rowsErr + [armErr, zeros(N, 1)];
% the nav: the INS increment carries a 5 mrad heading error, the feedback takes a tenth of the state at each epoch
psi = 5e-3;  G = 0.1;  KF = 0.5;
tKf = (0.5:0.5:T - 0.5)';  K = numel(tKf);
vNav = zeros(N, 3);  e = [0.02, -0.01, 0];  vNav(1, :) = vTrue(1, :) + e;
xPlus = zeros(K, 3);  kNext = 1;  fbSample = zeros(N, 3);
for j = 2:N
    dvT = vTrue(j, :) - vTrue(j - 1, :);
    e = e + psi * [-dvT(2), dvT(1), 0];                              % the increment rotated by the heading error
    if kNext <= K && tN(j) >= tKf(kNext)
        % the update at this epoch: the state estimates half of the nav's disagreement with the rows, the nav takes a tenth
        jm = j - 1;
        dis = (vTrue(jm, :) + e) - vRowsTrue(jm, :);
        xPlus(kNext, :) = KF * dis;
        e = e - G * xPlus(kNext, :);
        fbSample(j, :) = G * xPlus(kNext, :);                           % the correction, logged at the sample it was applied
        kNext = kNext + 1;
    end
    vNav(j, :) = vTrue(j, :) + e + 0.003 * randn(1, 3);             % the logged velocity carries a little noise
end
% the reference sits 0.40 m ahead and 0.25 m left of the IMU: its velocity carries yaw rate x that offset,
% truth minus rows per rad/s along +0.25 m and across +0.40 m
TS = [0.25, 0.40];
trueVelRef = vTrue;  trueVelRef(:, 1:2) = vTrue(:, 1:2) + rateTrue .* (TS(1) * uAl + TS(2) * uAc);
run = struct('t', tN, 'navVel', vNav, 'tTrue', tN, 'trueVel', trueVelRef, 'navAtt', [zeros(N, 2), course + psi]);
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
% the rows as a correct arm would give them, scored against the truth at the rows' point: the main run
CORR = struct('armSlope', ARM, 'truthSlope', TS);
txt = evalc('res = dopplerEpisodeCheck(log, tKf, run, CORR);');
fprintf('%s', txt);
ep = res.episodes;
% the rows as recorded: the corners fire the screen, and the read-back is the arm plus the truth point's offset
txtRaw = evalc('resRaw = dopplerEpisodeCheck(log, tKf, run);');
% (nav minus rows also carries the heading error's work in the corners, so the host's arm reads to within 0.08 m)
chk.armRaw   = all(abs(resRaw.armFit - (ARM + TS)) < 0.06) && any(resRaw.flag(tKf >= 118 & tKf <= 130)) && ...
    all(abs(resRaw.armHost - ARM) < 0.08) && all(abs((resRaw.armFit - resRaw.armHost) - TS) < 0.1) && ...
    ~isempty(strfind(txtRaw, 'pass opts.armSlope = [')) && ~isempty(strfind(txtRaw, 'and opts.truthSlope = ['));
% the arm taken out alone: the screen is quiet in the corners but the scoring still sees the truth point's offset there
txtArm = evalc('resArm = dopplerEpisodeCheck(log, tKf, run, struct(''armSlope'', ARM));');
chk.armOnly  = all(abs(resArm.armFit - TS) < 0.05) && numel(resArm.episodes) == 4 && ~any(resArm.flag(tKf >= 95 & tKf < 170)) && ...
    ~isempty(strfind(txtArm, 'opts.truthSlope = these two numbers'));
chk.armFixed = all(abs(res.armFit) < 0.05) && all(abs(res.armHost - ARM) < 0.08) && ~isempty(strfind(txt, 'after the arm correction')) && ...
    ~isempty(strfind(txt, 'truth point moved'));
% the records and the run on the same clock, and a run whose clock runs 1 s ahead, found and shifted back
% (the state is persistent, so the lag of the feedback is read to within one epoch)
chk.aligned  = res.alignLag >= 0 && res.alignLag <= 2 && isempty(strfind(txt, 'NOT ALIGNED'));
runAhead = run;  runAhead.t = tN + 3.0;  runAhead.tTrue = tN + 3.0;
txtS = evalc('resS = dopplerEpisodeCheck(log, tKf, runAhead, CORR);');
resS2 = dopplerEpisodeCheck(log, tKf, runAhead, setfield(setfield(CORR, 'tShift', -3.0), 'quiet', true)); %#ok<SFLD>
chk.shifted  = abs(resS.alignLag - 7) <= 1 && ~isempty(strfind(txtS, 'NOT ALIGNED')) && resS2.alignLag >= 0 && resS2.alignLag <= 2 && ...
    numel(resS2.episodes) == 2 && all([resS2.episodes.detected]) && resS2.nFalse <= 2;
% every third epoch's rows paired with the pulse half an epoch before: read as a time offset, not screened away
log3 = log;
for k = 3:3:K
    j = find(tN < tKf(k), 1, 'last') - 50;
    for s = 1:10
        log3{k}.innovation(10 + s) = log3{k}.obsMatrix(10 + s, 4:6) * (vNav(j + 50, :) - vRowsTrue(j, :))' + 0.05 * randn + 0.02 * double(s > 5);
    end
end
txt3 = evalc('res3p = dopplerEpisodeCheck(log3, tKf, run, CORR);');
chk.pairing  = res.pairJump < 0.05 && isempty(strfind(txt, 'PAIRING')) && res3p.pairJump > 0.4 && ~isempty(strfind(txt3, 'PAIRING')) && ...
    ~isempty(res3p.jumpClusters) && ~isempty(strfind(txt3, 'the jumps cluster at'));
% the same nav with records that carry a prior: the state accumulates a remainder of nine tenths, the nav took a
% tenth of each update's increment at lag 1, so the share of the state is diluted and the increment's delivery
% profile, a tenth at lag 1 and nothing elsewhere, is the one used
logI = log;  xAcc = zeros(1, 3);
for k = 1:K
    xPri = 0.9 * xAcc;  xAcc = xPri + xPlus(k, :);
    logI{k}.priorState = [zeros(3, 1); xPri'; zeros(n - 6, 1)];
    logI{k}.postState  = [zeros(3, 1); xAcc'; zeros(n - 6, 1)];
end
txtI = evalc('resI = dopplerEpisodeCheck(logI, tKf, run, setfield(CORR, ''fbModel'', ''profile''));'); %#ok<SFLD>
% (the synthetic increments are persistent, so the lags share the total between them; the total is what is read;
% these records are not in the error-state form, so the profile model is forced over the default total-state one)
chk.incModel = strcmp(resI.fbModel, 'profile') && abs(resI.firTotal - G) < 0.03 && resI.firTotal / resI.firTotalSe > 5 && ...
    ~isempty(strfind(txtI, 'slews the Kalman step in')) && abs(resI.fbScan(2, 2)) < 0.06 && ...
    numel(resI.episodes) == 2 && all([resI.episodes.detected]) && resI.nFalse <= 2;
% the host's velocity correction logged, per epoch and per sample: nothing fitted, the flags as with the fitted share
runL = run;  runL.fbVel = G * xPlus;
txtL = evalc('resL = dopplerEpisodeCheck(log, tKf, runL, CORR);');
runS = run;  runS.fbVel = fbSample;
txtS2 = evalc('resS2b = dopplerEpisodeCheck(log, tKf, runS, CORR);');
resW = dopplerEpisodeCheck(log, tKf, runL, setfield(setfield(CORR, 'fbLagWeights', [0 1]), 'quiet', true)); %#ok<SFLD>
chk.given2   = strcmp(resW.fbModel, 'logged') && numel(resW.episodes) == 2 && all([resW.episodes.detected]) && resW.nFalse <= 2 && ...
    nnz(resW.flag ~= res.flag) <= 2;
chk.logged   = strcmp(resL.fbModel, 'logged') && abs(resL.lgTotal - 1) < 0.1 && ~isempty(strfind(txtL, 'USED: the logged velocity correction, per epoch')) && ...
    strcmp(resS2b.fbModel, 'logged') && abs(resS2b.lgTotal - 1) < 0.1 && ~isempty(strfind(txtS2, 'USED: the logged velocity correction, per sample')) && ...
    numel(resL.episodes) == 2 && all([resL.episodes.detected]) && resL.nFalse <= 2 && ...
    numel(resS2b.episodes) == 2 && all([resS2b.episodes.detected]) && resS2b.nFalse <= 2 && nnz(resS2b.flag ~= res.flag) <= 2;
% the same records with the state vector and the H columns in ENU: read as NED the share is lost and the swap
% finds it; read as ENU everything is as before
logE = log;
for k = 1:K
    xk = log{k}.postState;  Hk = log{k}.obsMatrix;
    logE{k}.postState = [xk(2); xk(1); -xk(3); xk(5); xk(4); -xk(6); xk(7:end)];
    logE{k}.obsMatrix = [Hk(:, 2), Hk(:, 1), -Hk(:, 3), Hk(:, 5), Hk(:, 4), -Hk(:, 6), Hk(:, 7:end)];
end
txtEn = evalc('resEn = dopplerEpisodeCheck(logE, tKf, run, CORR);');
resEe = dopplerEpisodeCheck(logE, tKf, run, setfield(setfield(CORR, 'frame', 'ENU'), 'quiet', true)); %#ok<SFLD>
[~, iSw] = max(abs(resEn.fbScanSwap(:, 3)));
chk.frame    = abs(resEn.fbGain) < 0.03 && resEn.frameOff && ~isempty(strfind(txtEn, 'pass opts.frame = ''ENU''')) && ...
    abs(resEn.fbScanSwap(iSw, 2) - G) < 0.02 && abs(resEe.fbGain - G) < 0.02 && ~resEe.frameOff && ...
    numel(resEe.episodes) == 2 && all([resEe.episodes.detected]) && resEe.nFalse <= 2 && isequal(resEe.flag, res.flag);
% the feedback found at lag 1, and at lag 0 when the logged sample before an epoch already carries that epoch's feedback
chk.fbLag    = res.fbLag <= 1 && res.fbScan(2, 3) > 10;             % a persistent state shows at every lag; the share is what matters
res0 = dopplerEpisodeCheck(log, tKf + 0.011, run, setfield(CORR, 'quiet', true)); %#ok<SFLD>
chk.lag0     = res0.fbLag <= 1 && abs(res0.fbGain - G) < 0.02 && numel(res0.episodes) == 2 && all([res0.episodes.detected]) && res0.nFalse <= 2;
chk.two      = numel(ep) == 2;
chk.caught   = numel(ep) == 2 && all([ep.detected]) && all([ep.delay] <= 4.0);
chk.covered  = numel(ep) == 2 && ep(1).coverage > 0.6 && ep(2).coverage > 0.7;
chk.times    = numel(ep) == 2 && abs(ep(1).tStart - 200) < 1 && abs(ep(2).tStart - 301) < 2.5 && abs(ep(2).tEnd - 321) < 2.5;
chk.noFalse  = res.nFalse <= 2;
chk.tail     = res.nTail <= 4 + 7;                                    % the release latency: 4 after the short episode, 4 + 3 after the 19 s one
chk.fbGain   = abs(res.fbGain - G) < 0.02;
chk.gap      = ~any(res.flag(tKf >= 350 & tKf < 360)) && all(isnan(res.onset(tKf >= 350.5 & tKf < 355)));
chk.verdict  = ~isempty(strfind(txt, 'every episode caught, no false flag'));
% the manoeuvres alone are honest: no flag in the first braking and the two turns
chk.quiet    = ~any(res.flag(tKf >= 95 & tKf < 170)) && ~any(res.flag(tKf >= 245 & tKf < 270));
% a floor above the small episode: it is missed and said so; the large one still caught
txt2 = evalc('res2 = dopplerEpisodeCheck(log, tKf, run, setfield(CORR, ''floor'', 0.5));'); %#ok<SFLD>
chk.missed   = ~isempty(strfind(txt2, 'MISSED')) && ~isempty(strfind(txt2, 'episodes missed')) && numel(res2.episodes) == 2 && ...
    ~res2.episodes(1).detected;
% a given feedback share is used as given
res3 = dopplerEpisodeCheck(log, tKf, run, setfield(setfield(CORR, 'fbGain', 0.1), 'quiet', true)); %#ok<SFLD>
chk.given    = res3.fbGain == 0.1 && res3.fbLag <= 1 && numel(res3.episodes) == 2 && all([res3.episodes.detected]);
%% the host's error-state form: the rows measured against the nav less the prior, the prior and the step in the records
[runT, logT, tKfT, auxT] = dopplerSynthRun(struct('seed', 21));
txtT = evalc('resT = dopplerEpisodeCheck(logT, tKfT, runT, CORR);');
fprintf('%s', txtT);
epT = resT.episodes;
chk.totalUsed    = strcmp(resT.fbModel, 'total') && isnan(resT.fbGain) && ...
    ~isempty(strfind(txtT, 'USED: the nav velocity less the prior state'));
chk.totalCaught  = numel(epT) == 2 && all([epT.detected]) && all([epT.delay] <= 4.0) && epT(1).coverage > 0.6 && epT(2).coverage > 0.7;
chk.totalClean   = resT.nFalse <= 2 && resT.nTail <= 4 + 7 && ~any(resT.flag(tKfT >= 95 & tKfT < 170)) && ~any(resT.flag(tKfT >= 245 & tKfT < 270));
% the hand-over cancels in nav less prior: the increment's error against the truth is the logged noise, and no share
% of the state is left in it
chk.totalShare   = abs(resT.totShare) < 0.015 && abs(resT.totShareT) < 3 && resT.totRms < 0.01 && isempty(strfind(txtT, 'a share left'));
% the state's loss between an update and the next prior: a tenth of the posterior on every sizeable epoch, and the nav
% took all of it
tenth = resT.lossClass(strcmp({resT.lossClass.name}, 'a tenth'));
chk.totalLoss    = ~resT.lossNotTaken && tenth.n > 100 && abs(tenth.slope) < 0.1 && ...
    sum([resT.lossClass.n]) - tenth.n <= 2 && isempty(strfind(txtT, 'LOSS NOT TAKEN'));
% the host's arm error read from the nav and the rows themselves, not diluted by what the state absorbed in the corner
chk.totalArm     = all(abs(resT.armHost - ARM) < 0.08) && all(abs(resT.armFit) < 0.05);
% without the truth and without a logged correction: nothing to fit, the same flags
runNT = rmfield(rmfield(rmfield(runT, 'tTrue'), 'trueVel'), 'fbVel');
resNT = dopplerEpisodeCheck(logT, tKfT, runNT, struct('armSlope', ARM, 'quiet', true));
chk.totalNoTruth = strcmp(resNT.fbModel, 'total') && isequal(resNT.flag, resT.flag);
% a mechanisation that takes only three quarters of the hand-over: the quarter left reads as a share of the state
[runP, logP, tKfP] = dopplerSynthRun(struct('seed', 21, 'delivery', 0.75));
txtP = evalc('resP = dopplerEpisodeCheck(logP, tKfP, runP, CORR);');
chk.partial      = resP.totShare > 0.015 && resP.totShare < 0.045 && resP.totShareT >= 3 && ~isempty(strfind(txtP, 'a share left'));
% a feedback branch that fires on a large disagreement: the nav takes a quarter of the posterior and the state keeps a
% quarter, so half of it leaves nav minus prior at once; read as a loss of three quarters of which the nav took a third,
% the increment errors above 0.10 m/s sit on those epochs, and the verdict names it
% (the first episode is made 12 s long here: the state keeps a quarter, so the disagreement stays and the branch fires on
% most of its epochs; the corners fire it too, but a sharp corner is outside the read's 3 m/s2 gate)
[runQ, logQ, tKfQ, auxQ] = dopplerSynthRun(struct('seed', 21, 'limiter', 0.15, 'episodes', [200 212 0.35 0 0; 300 323 0 -0.4 3]));
txtQ = evalc('resQ = dopplerEpisodeCheck(logQ, tKfQ, runQ, CORR);');
threeQ = resQ.lossClass(strcmp({resQ.lossClass.name}, 'three quarters'));
chk.limiter      = nnz(auxQ.limited) >= 10 && threeQ.n >= 10 && threeQ.nUsed >= 10 && abs(threeQ.slope - 2 / 3) < 0.15 && ...
    resQ.lossNotTaken && ~isempty(strfind(txtQ, 'LOSS NOT TAKEN')) && ~isempty(strfind(txtQ, 'three quarters on'));
% the hand-over form of the coast: the same rows' velocity, the increment from the nav's change plus the logged
% correction; without a model term between updates the two forms are one and the same, so the flags are identical,
% and on the limited run the hand-over form has the exact coast where the total-state one has the loss
HAND = setfield(setfield(CORR, 'fbModel', 'handover'), 'fbLagWeights', [0 1]); %#ok<SFLD>
txtH = evalc('resH = dopplerEpisodeCheck(logT, tKfT, runT, HAND);');
resHQ = dopplerEpisodeCheck(logQ, tKfQ, runQ, setfield(HAND, 'quiet', true)); %#ok<SFLD>
chk.handover     = strcmp(resH.fbModel, 'handover') && isequal(resH.flag, resT.flag) && abs(resH.hoRms - resT.totRms) < 1e-3 && ...
    ~isempty(strfind(txtH, 'USED: the nav velocity less the prior state for the rows')) && ~isempty(strfind(txtH, 'the two coasts against the truth')) && ...
    resHQ.hoRms < resQ.totRms && resHQ.nFalse < resQ.nFalse && numel(resHQ.episodes) == 2 && all([resHQ.episodes.detected]);
% the fallback models still run on these records when forced; their rows' velocity is off by the prior state, which
% absorbs half of the small episode within two epochs, so only the large one is read and caught
resTL = dopplerEpisodeCheck(logT, tKfT, runT, setfield(setfield(CORR, 'fbModel', 'logged'), 'quiet', true)); %#ok<SFLD>
chk.totalForced  = strcmp(resTL.fbModel, 'logged') && numel(resTL.episodes) >= 1 && all([resTL.episodes.detected]) && resTL.nFalse <= 2;
names = fieldnames(chk);
for i = 1:numel(names)
    if ~chk.(names{i}), fprintf('   FAILED sub-check: %s\n', names{i}); ok = false; end
end
fprintf('error-state form: flags %d, false %d, tails %d, share of the state left %+.3f (t %.1f), increment rms %.4f m/s; three-quarter delivery read %+.3f (t %.1f)\n', ...
    resT.nFlag, resT.nFalse, resT.nTail, resT.totShare, resT.totShareT, resT.totRms, resP.totShare, resP.totShareT);
fprintf('limited feedback: %d limited epochs, loss of three quarters on %d, slope %+.2f on %d (the nav took %.2f); flags %d, false %d\n', ...
    nnz(auxQ.limited), threeQ.n, threeQ.slope, threeQ.nUsed, 1 - threeQ.slope, resQ.nFlag, resQ.nFalse);
fprintf('hand-over form: flags differing from total-state %d; on the limited run the coasts'' rms %.4f against %.4f m/s, false %d against %d\n', ...
    nnz(resH.flag ~= resT.flag), resHQ.hoRms, resQ.totRms, resHQ.nFalse, resQ.nFalse);
if numel(ep) == 2
    fprintf('episodes: %.1f s caught after %.1f s (%.0f%%), %.1f s caught after %.1f s (%.0f%%); false flags %d, release tails %d; feedback share %.3f at lag %d\n', ...
        ep(1).tEnd - ep(1).tStart, ep(1).delay, 100 * ep(1).coverage, ep(2).tEnd - ep(2).tStart, ep(2).delay, 100 * ep(2).coverage, ...
        res.nFalse, res.nTail, res.fbGain, res.fbLag);
    fprintf('wide scan lag %d; run 3 s ahead read at lag %d, shifted back %d; arm read-back raw [%+.2f %+.2f], arm only [%+.2f %+.2f]; pairing jumps %.0f%% / %.0f%%\n', ...
        res.alignLag, resS.alignLag, resS2.alignLag, resRaw.armFit(1), resRaw.armFit(2), resArm.armFit(1), resArm.armFit(2), ...
        100 * res.pairJump, 100 * res3p.pairJump);
end
if ok
    fprintf('=== test_dopplerEpisodeCheck PASS ===\n');
else
    fprintf('=== test_dopplerEpisodeCheck FAIL ===\n');
end
end
