function [res] = dopplerEpisodeCheck(spfMeasLog, tKf, run, opts)
%DOPPLEREPISODECHECK  Would a screen of the Doppler rows against the INS velocity catch a receiver episode,
%   and what would it cost on the clean hour? An offline emulation on the records and the run.
%   res = dopplerEpisodeCheck(spfMeasLog, tKf, run)
%   res = dopplerEpisodeCheck(spfMeasLog, tKf, run, opts)
%
%   The receiver's Doppler goes wrong in some episodes of strong transient dynamics (the vehicle data: 0.4 m/s
%   for 20 s at a tunnel portal, all rows together, so no consensus screen sees it) and the filter, following
%   its rows, carries the error into the states. The screen emulated here compares the rows' own velocity with
%   the INS over a short span and, once they part, holds the Doppler rows out until they agree with an INS
%   velocity coast again or a hold limit runs out:
%     rows' velocity  v_rows(k) = v_nav-(k) - dv(k): the nav velocity before the update minus the least-squares
%                     fit of the rate-row innovations to a velocity and one clock drift per constellation (the
%                     host's innovation is predicted minus measured; opts.sign = +1 for the other convention)
%     INS increment   dvIns(k) = v_nav-(k) - v_nav-(k-1) minus the feedback the nav took in between (fbGain
%                     times the state's velocity, the gain fitted against the truth when it is given, at the lag
%                     where the logged velocity shows it: lag 1 when the sample before an epoch precedes that
%                     epoch's feedback, lag 0 when it already carries it)
%     lever arm       with opts.armSlope = [along across] (m per rad/s of yaw rate, leverArmCheck's slopes) and
%                     run.navAtt, the lever-arm error the rows carry is taken out of their velocity, as a
%                     corrected arm input would: in the corners that error is yaw rate times arm and reads as
%                     an episode otherwise
%     onset test      d(k) = [v_rows(k) - v_rows(k-span)] - sum(dvIns) over the span, horizontal, averaged over
%                     'avg' consecutive epochs, against T = sqrt(floor^2 + (headingSigma * sum|dvIns|)^2): the
%                     rows changed by more than the INS says, beyond what the row noise and a heading error on
%                     that much speed change allow; 'confirm' consecutive epochs above T fire the onset
%     episode         from the onset, a coast v_c anchored at the rows' velocity before the span and propagated
%                     with dvIns; the rows flagged while |v_rows - v_c| stays above the floor, released after
%                     'release' consecutive epochs below it, or at 'hold' epochs at the latest; after a release
%                     the onset test rests for span + avg epochs, since its reference epochs still sit in the
%                     episode and the rows' return would otherwise fire it again
%   Scoring needs the truth (run.tTrue, run.trueVel): the rows' error v_rows - v_truth labels the epochs the
%   receiver was actually wrong (above 'episodeLevel' for at least 'episodeMin' epochs); the printout gives
%   each true episode with its detection delay and coverage, the release tails (flags within 'release' epochs
%   after an episode, the rows already back within 'cleanLevel' of the truth: the latency of the release, not
%   an error), the false flags (flagged while the rows were within 'cleanLevel' of the truth, not a tail)
%   grouped into runs, and a per-epoch table around the largest episode. With run.navAtt it also regresses
%   the rows' error against the truth on the yaw rate over the turning epochs: the lever-arm error the rows
%   carry against the truth point, in leverArmCheck's convention, before or after the arm correction.
%
%   run: t, navVel (100 Hz), tTrue, trueVel (NED, m/s), navAtt ([roll pitch yaw], rad, for the arm); opts
%   (defaults): span 10 epochs, avg 2, confirm 2, floor 0.20 m/s, headingSigma 5e-3 rad, hold 40 epochs,
%   release 4 epochs, sign -1, fbGain [] (fitted; a number is used as given, with its sign), fbLag [] (0, 1 or 2
%   epochs; default the lag the fit finds, 1 without truth), armSlope [] (no arm correction),
%   episodeLevel 0.25 m/s, episodeMin 4 epochs, cleanLevel 0.15 m/s, tWin [] (the
%   epochs to run on; default all), tTable [] (the window of the per-epoch table; default the largest episode
%   with 10 s either side), quiet false.
%   res: t, vRows [K x 3], dvIns [K x 3], onset [K x 1], thresh [K x 1], flag [K x 1] logical, sep [K x 1],
%   rowsErr [K x 3] (NaN without truth), episodes (struct array: tStart, tEnd, detected, delay, coverage),
%   falseRuns (struct array: tStart, tEnd, epochs), nFlag, nFalse, nTail, fbGain, fbLag, fbScan ([3 x 3]: lag,
%   share, t-statistic), yawRate [K x 1], armFit [1 x 2] (along, across, m per rad/s, NaN without navAtt or
%   truth), armFitSig, opts.
%   MATLAB desktop tool, not codegen. Runs under Octave.

if iscell(spfMeasLog), spfMeasLog = [spfMeasLog{:}]; end
if nargin < 4 || isempty(opts), opts = struct(); end
o = struct('span', 10, 'avg', 2, 'confirm', 2, 'floor', 0.20, 'headingSigma', 5e-3, 'hold', 40, 'release', 4, 'sign', -1, 'fbGain', [], ...
    'fbLag', [], 'armSlope', [], 'episodeLevel', 0.25, 'episodeMin', 4, 'cleanLevel', 0.15, 'tWin', [], 'tTable', [], 'quiet', false);
f = fieldnames(opts);
for i = 1:numel(f), o.(f{i}) = opts.(f{i}); end
K = numel(spfMeasLog);
t = tKf(:);  t = t(1:K);
tNav = run.t(:);  navVel = run.navVel;
hasTruth = isfield(run, 'tTrue') && isfield(run, 'trueVel') && ~isempty(run.trueVel);
hasAtt = isfield(run, 'navAtt') && ~isempty(run.navAtt);

% the nav velocity just before each update, the rows' velocity, the state's velocity after the update
vNavM = NaN(K, 3);  dv = NaN(K, 3);  xPlus = NaN(K, 3);  nRate = zeros(K, 1);  tSample = t;
for k = 1:K
    j = find(tNav < t(k), 1, 'last');
    if isempty(j), continue; end
    vNavM(k, :) = navVel(j, :);  tSample(k) = tNav(j);               % the rows and the truth are read at this sample
    s = spfMeasLog(k);  m = double(s.numMeas);
    if m == 0, continue; end
    xPlus(k, :) = s.postState(4:6)';
    H = s.obsMatrix(1:m, :);  y = s.innovation(1:m);  y = y(:);
    rate = sqrt(sum(H(:, 4:6) .^ 2, 2)) > 0.5;
    if nnz(rate) < 4, continue; end
    cols = [4 5 6, 6 + find(any(H(rate, 7:end) ~= 0, 1))];
    sol = H(rate, cols) \ y(rate);
    dv(k, :) = sol(1:3)';  nRate(k) = nnz(rate);
end
inWin = true(K, 1);
if ~isempty(o.tWin), inWin = t >= o.tWin(1) & t <= o.tWin(2); end

% the INS increment between updates: the nav's change minus the feedback it took
dNav = [NaN(1, 3); diff(vNavM)];
if hasTruth
    vTrue = interp1(run.tTrue(:), run.trueVel, tSample);              % at the nav sample the rows are compared with
    dTrue = [NaN(1, 3); diff(vTrue)];
else
    vTrue = NaN(K, 3);  dTrue = NaN(K, 3);
end
% the share of the state the nav took, and the lag at which the logged velocity shows it:
% (nav change - truth change)(k) = -fb * x+(k - lag) + INS drift over the quiet epochs, both axes stacked
LAGS = 0:2;  T_MIN = 3;
fbScan = [LAGS(:), NaN(numel(LAGS), 2)];
if hasTruth
    for iL = 1:numel(LAGS)
        L = LAGS(iL);
        kk = (L + 2:K)';
        a = xPlus(kk - L, 1:2);  b = dNav(kk, 1:2) - dTrue(kk, 1:2);
        ok = all(isfinite([a b]), 2) & inWin(kk) & sqrt(sum(dTrue(kk, 1:2) .^ 2, 2)) < 0.5;
        av = a(ok, :);  bv = b(ok, :);  av = av(:);  bv = bv(:);
        if numel(av) >= 10 && av' * av > 0
            c = -(av' * bv) / (av' * av);
            resid = bv + c * av;
            se = sqrt((resid' * resid) / max(numel(av) - 1, 1) / (av' * av));
            fbScan(iL, 2:3) = [c, c / max(se, eps)];
        end
    end
end
scanOk = any(isfinite(fbScan(:, 3)));
if scanOk, [tBest, iBest] = max(abs(fbScan(:, 3))); else, tBest = NaN;  iBest = 2; end
fbFound = isfinite(tBest) && tBest >= T_MIN;
if isempty(o.fbGain)
    if ~hasTruth
        error('dopplerEpisodeCheck:fbGain', 'give opts.fbGain (the share of the state the nav takes per epoch, signed) when the run has no truth');
    end
    if fbFound
        fbGain = fbScan(iBest, 2);  fbLag = fbScan(iBest, 1);
    else
        fbGain = 0;  fbLag = 1;                                       % no feedback of the velocity state in the logged velocity
    end
else
    fbGain = o.fbGain;
    if ~isempty(o.fbLag), fbLag = o.fbLag;
    elseif fbFound, fbLag = fbScan(iBest, 1);
    else, fbLag = 1; end
end
xLag = NaN(K, 3);  xLag(fbLag + 1:K, :) = xPlus(1:K - fbLag, :);
dvIns = dNav + fbGain * xLag;
% the rows' velocity: the nav velocity the host predicted with, before that epoch's feedback, minus the fitted
% innovation velocity (predicted minus measured: the rows sit at nav - dv)
vNavPre = vNavM;
if fbLag == 0, vNavPre = vNavM + fbGain * xPlus; end
vRows = vNavPre + o.sign * dv;
% the yaw rate and the along / across unit vectors, for the lever arm
yawRate = NaN(K, 1);  uAlong = NaN(K, 2);  uAcross = NaN(K, 2);
if hasAtt
    yaw = unwrap(run.navAtt(:, 3));
    dtN = median(diff(tNav));  w = max(1, round(0.5 / dtN));
    rate = conv([diff(yaw); 0] / dtN, ones(w, 1) / w, 'same');              % centred over 0.5 s
    yawRate = interp1(tNav, rate, tSample, 'linear', NaN);
    yawK = interp1(tNav, yaw, tSample, 'linear', NaN);
    uAlong = [cos(yawK), sin(yawK)];  uAcross = [-sin(yawK), cos(yawK)];
end
if ~isempty(o.armSlope)
    if ~hasAtt, error('dopplerEpisodeCheck:armSlope', 'opts.armSlope needs run.navAtt for the yaw rate'); end
    armVel = yawRate .* (o.armSlope(1) * uAlong + o.armSlope(2) * uAcross);
    vRows(:, 1:2) = vRows(:, 1:2) + armVel;                                 % the rows as a correct arm would give them
end

% the emulated screen
% the single-epoch statistic: the rows' change over the span against the INS's
dSingle = NaN(K, 2);  speedChange = NaN(K, 1);
for k = o.span + 1:K
    if ~inWin(k), continue; end
    if all(isfinite(vRows(k, :))) && all(isfinite(vRows(k - o.span, :))) && all(all(isfinite(dvIns(k - o.span + 1:k, :))))
        dSingle(k, :) = (vRows(k, 1:2) - vRows(k - o.span, 1:2)) - sum(dvIns(k - o.span + 1:k, 1:2), 1);
        speedChange(k) = sum(sqrt(sum(dvIns(k - o.span + 1:k, 1:2) .^ 2, 2)));
    end
end
onset = NaN(K, 1);  thresh = NaN(K, 1);  flag = false(K, 1);  sep = NaN(K, 1);
inEp = false;  vCoast = zeros(1, 3);  hold = 0;  rel = 0;  above = 0;  refractory = 0;
for k = o.span + o.avg:K
    if ~inWin(k), continue; end
    dWin = dSingle(k - o.avg + 1:k, :);
    have = all(isfinite(dWin(:)));
    if have
        onset(k) = norm(mean(dWin, 1));                          % the statistic averaged over 'avg' epochs
        thresh(k) = sqrt(o.floor ^ 2 + (o.headingSigma * speedChange(k)) ^ 2);
    end
    if ~inEp
        if refractory > 0
            refractory = refractory - 1;                          % after a release the reference epochs still sit in the episode
            above = 0;
        elseif have && onset(k) > thresh(k)
            above = above + 1;
            if above >= o.confirm
                inEp = true;  hold = 0;  rel = 0;  above = 0;
                k0 = k - o.span - o.avg + 1;                      % the last epoch the rows and the INS agreed
                vCoast = vRows(k0, :) + sum(dvIns(k0 + 1:k, :), 1);
                flag(k) = true;  sep(k) = norm(vRows(k, 1:2) - vCoast(1:2));
            end
        else
            above = 0;
        end
    else
        if all(isfinite(dvIns(k, :))), vCoast = vCoast + dvIns(k, :); end
        hold = hold + 1;  flag(k) = true;
        if all(isfinite(vRows(k, :)))
            sep(k) = norm(vRows(k, 1:2) - vCoast(1:2));
            if sep(k) < o.floor, rel = rel + 1; else, rel = 0; end
        end
        if rel >= o.release || hold >= o.hold
            inEp = false;  refractory = o.span + o.avg;
        end
    end
end

% scoring against the truth
rowsErr = vRows - vTrue;
% the lever-arm error the rows carry against the truth point: the rows' error along and across the heading
% regressed on the yaw rate over the turning epochs, with an intercept; reported as innovation per rad/s
% (nav minus rows), leverArmCheck's convention
armFit = NaN(1, 2);  armFitSig = NaN(1, 2);  nTurn = 0;
if hasTruth && hasAtt
    turn = abs(yawRate) > 2 * pi / 180 & all(isfinite(rowsErr(:, 1:2)), 2) & inWin;
    nTurn = nnz(turn);
    if nTurn >= 10
        eA = sum(rowsErr(turn, 1:2) .* uAlong(turn, :), 2);  eX = sum(rowsErr(turn, 1:2) .* uAcross(turn, :), 2);
        A = [yawRate(turn), ones(nTurn, 1)];
        for p = 1:2
            if p == 1, yv = eA; else, yv = eX; end
            c = A \ yv;  r = yv - A * c;
            covC = (r' * r) / max(nTurn - 2, 1) * inv(A' * A);
            armFit(p) = -c(1);  armFitSig(p) = sqrt(covC(1, 1));
        end
    end
end
episodes = struct('tStart', {}, 'tEnd', {}, 'detected', {}, 'delay', {}, 'coverage', {});
falseRuns = struct('tStart', {}, 'tEnd', {}, 'epochs', {});
nFalse = 0;  nTail = 0;
if hasTruth
    errH = sqrt(sum(rowsErr(:, 1:2) .^ 2, 2));
    bad = errH > o.episodeLevel & inWin;
    runsB = runsOf(bad);
    for i = 1:size(runsB, 1)
        a = runsB(i, 1);  b = runsB(i, 2);
        if b - a + 1 < o.episodeMin, continue; end
        fl = flag(a:b);
        e = struct('tStart', t(a), 'tEnd', t(b), 'detected', any(fl), 'delay', NaN, 'coverage', mean(fl));
        if any(fl), e.delay = t(a - 1 + find(fl, 1)) - t(a); end
        episodes(end + 1) = e; %#ok<AGROW>
    end
    % a flag on a clean epoch within 'release' epochs of a flagged epoch that was not clean is the release latency
    % after an episode, booked apart from the false flags
    cleanFlag = flag & errH < o.cleanLevel & inWin;
    tailFlag = false(K, 1);
    for k = find(cleanFlag)'
        a = max(k - o.release, 1);
        tailFlag(k) = any(flag(a:k - 1) & errH(a:k - 1) >= o.cleanLevel);
    end
    falseFlag = cleanFlag & ~tailFlag;
    nFalse = nnz(falseFlag);
    nTail = nnz(tailFlag);
    runsF = runsOf(falseFlag);
    for i = 1:size(runsF, 1)
        falseRuns(end + 1) = struct('tStart', t(runsF(i, 1)), 'tEnd', t(runsF(i, 2)), 'epochs', runsF(i, 2) - runsF(i, 1) + 1); %#ok<AGROW>
    end
end

res = struct('t', t, 'vRows', vRows, 'dvIns', dvIns, 'onset', onset, 'thresh', thresh, 'flag', flag, 'sep', sep, ...
    'rowsErr', rowsErr, 'episodes', episodes, 'falseRuns', falseRuns, 'nFlag', nnz(flag), 'nFalse', nFalse, 'nTail', nTail, ...
    'fbGain', fbGain, 'fbLag', fbLag, 'fbScan', fbScan, 'yawRate', yawRate, 'armFit', armFit, 'armFitSig', armFitSig, 'opts', o);
if o.quiet, return; end

% print
fprintf('dopplerEpisodeCheck: %d epochs, %d with rate rows; span %d epochs, averaged over %d, confirmed over %d, floor %.2f m/s, heading sigma %.1f mrad, hold %d, release %d\n', ...
    nnz(inWin), nnz(nRate > 0 & inWin), o.span, o.avg, o.confirm, o.floor, o.headingSigma * 1e3, o.hold, o.release);
scanTxt = '';
if scanOk
    scanTxt = sprintf('lag %d %+.3f t %.1f, lag %d %+.3f t %.1f, lag %d %+.3f t %.1f', fbScan');
end
if isempty(o.fbGain) && ~fbFound
    fprintf('   feedback share of the state per epoch in the nav velocity: NONE FOUND (%s): the logged velocity does not take the\n', scanTxt);
    fprintf('      velocity state at any lag, or the records and the run are not on the same time; the INS increment is used as logged\n');
else
    fprintf('   feedback share of the state per epoch in the nav velocity: %+.3f at lag %d (%s; %s)\n', fbGain, fbLag, ...
        pick(isempty(o.fbGain), 'fitted against the truth', 'given'), pick(scanOk, scanTxt, 'no scan without truth'));
end
if hasAtt && hasTruth
    if nTurn >= 10
        armTxt = '';
        if ~isempty(o.armSlope), armTxt = sprintf(', after the arm correction [%+.2f %+.2f]', o.armSlope(1), o.armSlope(2)); end
        fprintf('   rows minus truth per rad/s of yaw rate over %d turning epochs%s: along %+.2f m, across %+.2f m (sigma %.2f / %.2f)\n', ...
            nTurn, armTxt, armFit(1), armFit(2), armFitSig(1), armFitSig(2));
        if isempty(o.armSlope) && any(abs(armFit) > 2 * armFitSig & abs(armFit) > 0.1)
            fprintf('      (the rows carry a lever-arm error against the truth point; opts.armSlope = these two numbers emulates a corrected arm)\n');
        end
    else
        fprintf('   rows minus truth per rad/s of yaw rate: fewer than 10 turning epochs, no arm read\n');
    end
end
fprintf('   flagged: %d epochs (%.1f s), %.2f%% of the epochs with rows\n', nnz(flag), nnz(flag) * median(diff(t)), 100 * nnz(flag) / max(nnz(nRate > 0 & inWin), 1));
if hasTruth
    fprintf('   true episodes (rows more than %.2f m/s from the truth for %d epochs or more): %d\n', o.episodeLevel, o.episodeMin, numel(episodes));
    for i = 1:numel(episodes)
        e = episodes(i);
        fprintf('      %8.2f to %8.2f s (%5.1f s)  %s   delay %5.1f s   coverage %3.0f%%\n', e.tStart, e.tEnd, e.tEnd - e.tStart, ...
            pick(e.detected, 'caught ', 'MISSED '), e.delay, 100 * e.coverage);
    end
    fprintf('   release tails (flagged within %d epochs after an episode, with the rows back within %.2f m/s of the truth): %d epochs\n', ...
        o.release, o.cleanLevel, nTail);
    fprintf('   false flags (flagged while the rows were within %.2f m/s of the truth, not the tail of an episode): %d epochs in %d runs\n', ...
        o.cleanLevel, nFalse, numel(falseRuns));
    for i = 1:min(numel(falseRuns), 12)
        fprintf('      %8.2f to %8.2f s  %3d epochs\n', falseRuns(i).tStart, falseRuns(i).tEnd, falseRuns(i).epochs);
    end
    if numel(falseRuns) > 12, fprintf('      ... %d more runs\n', numel(falseRuns) - 12); end
    % the per-epoch table around the largest episode
    tt = o.tTable;
    if isempty(tt) && ~isempty(episodes)
        [~, iBig] = max([episodes.tEnd] - [episodes.tStart]);
        tt = [episodes(iBig).tStart - 10, episodes(iBig).tEnd + 10];
    end
    if ~isempty(tt)
        fprintf('   per epoch, %.0f to %.0f s: rows - truth, the onset statistic against its threshold, the flag and the separation from the coast\n', tt(1), tt(2));
        fprintf('        t      rows-truth N  E   |.|     onset  thresh   flag    sep\n');
        for k = find(t >= tt(1) & t <= tt(2))'
            fprintf('   %8.2f    %+6.3f %+6.3f  %5.3f    %5.3f   %5.3f    %s   %5.3f\n', t(k), rowsErr(k, 1), rowsErr(k, 2), errH(k), ...
                onset(k), thresh(k), pick(flag(k), 'FLAG', '    '), sep(k));
        end
    end
    if isempty(episodes)
        fprintf('   VERDICT: no episode in this window; %d false flags\n', nFalse);
    elseif all([episodes.detected]) && nFalse == 0
        fprintf('   VERDICT: every episode caught, no false flag; the release latency cost %.1f s of Doppler after the episodes\n', nTail * median(diff(t)));
    elseif all([episodes.detected])
        fprintf('   VERDICT: every episode caught; the false flags cost %.1f s of Doppler over the window, the release latency %.1f s\n', ...
            nFalse * median(diff(t)), nTail * median(diff(t)));
    else
        fprintf('   VERDICT: %d of %d episodes missed: the floor or the span is too generous for this receiver, or the heading sigma hides it\n', ...
            nnz(~[episodes.detected]), numel(episodes));
    end
else
    fprintf('   (no truth in the run: flags only, no scoring)\n');
end
end

%------------------------------------------------------------------------
function [r] = runsOf(mask)
% start and end indices of the runs of true in a logical vector, [n x 2]
mask = mask(:)';
d = diff([false, mask, false]);
r = [find(d == 1)', find(d == -1)' - 1];
end

function [out] = pick(cond, a, b)
if cond, out = a; else, out = b; end
end
