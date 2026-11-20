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
%   Three reads on the inputs come first, since the screen is only as good as what it is fed:
%     feedback share  fitted at lags 0 to 2 with a t-statistic each, NONE FOUND when no lag reaches 3; the same
%                     fit with north and east swapped: a share that appears only then is a state vector in the
%                     other frame (opts.frame = 'ENU' reads the records' states and H columns as east, north, up)
%     alignment       the same fit over lags -12 to 12: the share showing up sharpest well outside 0 to 2 is
%                     the feedback where the run's clock put it, the records and the run not on the same
%                     clock (opts.tShift moves the run)
%     time pairing    on the epochs faster than 1 m/s2 the rows' error against the truth projected on the
%                     acceleration is a time offset; an offset that jumps between consecutive epochs is rows
%                     paired with the wrong pulse on some of them, which no screen can fix
%   Scoring needs the truth (run.tTrue, run.trueVel): the rows' error v_rows - v_truth labels the epochs the
%   receiver was actually wrong (above 'episodeLevel' for at least 'episodeMin' epochs); the printout gives
%   each true episode with its detection delay and coverage, the release tails (flags within 'release' epochs
%   after an episode, the rows already back within 'cleanLevel' of the truth: the latency of the release, not
%   an error), the false flags (flagged while the rows were within 'cleanLevel' of the truth, not a tail)
%   grouped into runs, and a per-epoch table around the largest episode. With run.navAtt it also regresses
%   the rows' error against the truth on the yaw rate over the turning epochs: the lever-arm error the rows
%   carry against the truth point, in leverArmCheck's convention, before or after the arm correction; what is
%   left after it is the truth point's own offset from the rows' point, taken out of the scoring with
%   opts.truthSlope (the printed read-back).
%
%   run: t, navVel (100 Hz), tTrue, trueVel (NED, m/s), navAtt ([roll pitch yaw], rad, for the arm); opts
%   (defaults): span 10 epochs, avg 2, confirm 2, floor 0.20 m/s, headingSigma 5e-3 rad, hold 40 epochs,
%   release 4 epochs, sign -1, fbGain [] (fitted; a number is used as given, with its sign), fbLag [] (0, 1 or 2
%   epochs; default the lag the fit finds, 1 without truth), armSlope [] (no arm correction), truthSlope []
%   (no truth-point correction), tShift 0 s (added to the run's time stamps to put them on the records'
%   clock), frame 'NED' (the records' position and velocity components; 'ENU' maps them to NED), episodeLevel
%   0.25 m/s, episodeMin 4 epochs, cleanLevel 0.15 m/s, tWin [] (the
%   epochs to run on; default all), tTable [] (the window of the per-epoch table; default the largest episode
%   with 10 s either side), quiet false.
%   res: t, vRows [K x 3], dvIns [K x 3], onset [K x 1], thresh [K x 1], flag [K x 1] logical, sep [K x 1],
%   rowsErr [K x 3] (NaN without truth), episodes (struct array: tStart, tEnd, detected, delay, coverage),
%   falseRuns (struct array: tStart, tEnd, epochs), nFlag, nFalse, nTail, fbGain, fbLag, fbScan ([3 x 3]: lag,
%   share, t-statistic), yawRate [K x 1], armFit [1 x 2] (along, across, m per rad/s, NaN without navAtt or
%   truth), armFitSig, fbScanWide ([25 x 3], lags -12 to 12), alignLag (its sharpest lag), fbScanSwap ([3 x 3],
%   north and east swapped), frameOff (true when the share appears only swapped), jumpClusters ([n x 3]: start,
%   end, count of the offset jumps 5 s apart or less), tauPair [K x 1] (the
%   time offset per dynamic epoch, NaN elsewhere), pairShare (the share of them beyond 0.35 s), pairJump (the
%   share of consecutive dynamic epochs whose offset jumps by more than 0.35 s), opts.
%   MATLAB desktop tool, not codegen. Runs under Octave.

if iscell(spfMeasLog), spfMeasLog = [spfMeasLog{:}]; end
if nargin < 4 || isempty(opts), opts = struct(); end
o = struct('span', 10, 'avg', 2, 'confirm', 2, 'floor', 0.20, 'headingSigma', 5e-3, 'hold', 40, 'release', 4, 'sign', -1, 'fbGain', [], ...
    'fbLag', [], 'armSlope', [], 'truthSlope', [], 'tShift', 0, 'episodeLevel', 0.25, 'episodeMin', 4, 'cleanLevel', 0.15, ...
    'tWin', [], 'tTable', [], 'quiet', false, 'frame', 'NED');
f = fieldnames(opts);
for i = 1:numel(f), o.(f{i}) = opts.(f{i}); end
if ~any(strcmp(o.frame, {'NED', 'ENU'})), error('dopplerEpisodeCheck:frame', 'opts.frame is ''NED'' or ''ENU'''); end
isENU = strcmp(o.frame, 'ENU');
K = numel(spfMeasLog);
t = tKf(:);  t = t(1:K);
tNav = run.t(:) + o.tShift;  navVel = run.navVel;                       % the run on the records' clock
tTrueV = [];
if isfield(run, 'tTrue'), tTrueV = run.tTrue(:) + o.tShift; end
hasTruth = isfield(run, 'tTrue') && isfield(run, 'trueVel') && ~isempty(run.trueVel);
hasAtt = isfield(run, 'navAtt') && ~isempty(run.navAtt);

% the nav velocity just before each update, the rows' innovation velocity, the state's velocity after the update
vNavM = NaN(K, 3);  dv = NaN(K, 3);  xPlus = NaN(K, 3);  nRate = zeros(K, 1);  tSample = t;
for k = 1:K
    j = find(tNav < t(k), 1, 'last');
    if isempty(j), continue; end
    vNavM(k, :) = navVel(j, :);  tSample(k) = tNav(j);               % the rows and the truth are read at this sample
    s = spfMeasLog(k);  m = double(s.numMeas);
    if m == 0, continue; end
    xPlus(k, :) = toNED(s.postState(4:6)', isENU);
    H = s.obsMatrix(1:m, :);  y = s.innovation(1:m);  y = y(:);
    rate = sqrt(sum(H(:, 4:6) .^ 2, 2)) > 0.5;
    if nnz(rate) < 4, continue; end
    cols = [4 5 6, 6 + find(any(H(rate, 7:end) ~= 0, 1))];
    sol = H(rate, cols) \ y(rate);
    dv(k, :) = toNED(sol(1:3)', isENU);  nRate(k) = nnz(rate);
end
inWin = true(K, 1);
if ~isempty(o.tWin), inWin = t >= o.tWin(1) & t <= o.tWin(2); end

% the yaw rate and the along / across unit vectors at the samples, for the lever arm and the truth point
yawRate = NaN(K, 1);  uAlong = NaN(K, 2);  uAcross = NaN(K, 2);
if hasAtt
    yaw = unwrap(run.navAtt(:, 3));
    dtN = median(diff(tNav));  w = max(1, round(0.5 / dtN));
    rate = conv([diff(yaw); 0] / dtN, ones(w, 1) / w, 'same');              % centred over 0.5 s
    yawRate = interp1(tNav, rate, tSample, 'linear', NaN);
    yawK = interp1(tNav, yaw, tSample, 'linear', NaN);
    uAlong = [cos(yawK), sin(yawK)];  uAcross = [-sin(yawK), cos(yawK)];
end
if ~isempty(o.armSlope) && ~hasAtt, error('dopplerEpisodeCheck:armSlope', 'opts.armSlope needs run.navAtt for the yaw rate'); end
if ~isempty(o.truthSlope) && ~hasAtt, error('dopplerEpisodeCheck:truthSlope', 'opts.truthSlope needs run.navAtt for the yaw rate'); end

% the truth at the samples, moved to the nav's point when the truth point's offset is given: the arm correction
% below puts the rows at the nav's point too, so nav, rows and truth are then compared at one point
if hasTruth
    vTrue = interp1(tTrueV, run.trueVel, tSample);
    if ~isempty(o.truthSlope)
        vTrue(:, 1:2) = vTrue(:, 1:2) - yawRate .* (o.truthSlope(1) * uAlong + o.truthSlope(2) * uAcross);
    end
    dTrue = [NaN(1, 3); diff(vTrue)];
else
    vTrue = NaN(K, 3);  dTrue = NaN(K, 3);
end

% the INS increment between updates: the nav's change minus the feedback it took
dNav = [NaN(1, 3); diff(vNavM)];
% the share of the state the nav took, and the lag at which the logged velocity shows it:
% (nav change - truth change)(k) = -fb * x+(k - lag) + INS drift over the quiet epochs, both axes stacked
LAGS = 0:2;  T_MIN = 3;  WIDE = -12:12;
fbScan = [LAGS(:), NaN(numel(LAGS), 2)];  fbScanWide = [WIDE(:), NaN(numel(WIDE), 2)];
if hasTruth
    for iL = 1:numel(WIDE)
        fbScanWide(iL, 2:3) = shareAtLag(xPlus, dNav, dTrue, inWin, WIDE(iL), K);
    end
    for iL = 1:numel(LAGS)
        fbScan(iL, 2:3) = fbScanWide(WIDE == LAGS(iL), 2:3);
    end
end
% the same share with north and east swapped: a share that appears only then is a state vector in the other
% frame (ENU against NED), the first two components exchanged
fbScanSwap = [LAGS(:), NaN(numel(LAGS), 2)];
if hasTruth
    xSwap = xPlus(:, [2 1 3]);
    for iL = 1:numel(LAGS)
        fbScanSwap(iL, 2:3) = shareAtLag(xSwap, dNav, dTrue, inWin, LAGS(iL), K);
    end
end
frameOff = false;
if any(isfinite(fbScanSwap(:, 3))) && any(isfinite(fbScan(:, 3)))
    frameOff = max(abs(fbScanSwap(:, 3))) >= T_MIN && max(abs(fbScanSwap(:, 3))) > 1.5 * max(abs(fbScan(:, 3)));
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
% innovation velocity (predicted minus measured: the rows sit at nav - dv), then the lever-arm error taken out,
% which puts the rows at the nav's point
vNavPre = vNavM;
if fbLag == 0, vNavPre = vNavM + fbGain * xPlus; end
vRows = vNavPre + o.sign * dv;
if ~isempty(o.armSlope)
    armVel = yawRate .* (o.armSlope(1) * uAlong + o.armSlope(2) * uAcross);
    vRows(:, 1:2) = vRows(:, 1:2) + armVel;                                 % the rows as a correct arm would give them
end

% alignment of the records with the run: the share found at a lag well outside 0 to 2 epochs, sharper than
% there, is the feedback showing up where the run's clock put it
alignLag = NaN;  alignOff = false;
if any(isfinite(fbScanWide(:, 3)))
    [tWide, iW] = max(abs(fbScanWide(:, 3)));
    alignLag = fbScanWide(iW, 1);
    tNear = max(abs(fbScan(:, 3)));
    alignOff = tWide >= T_MIN && (alignLag < 0 || alignLag > 2) && (~isfinite(tNear) || tWide > 1.5 * tNear);
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
% rows paired with the wrong time: on the epochs faster than 1 m/s2 the rows' error projected on the truth's
% acceleration is the time offset that would explain it; a receiver episode reads as a steady offset too, so the
% mark of a pairing fault is the offset jumping between consecutive epochs
tauPair = NaN(K, 1);  nDyn = 0;  pairShare = NaN;  pairHalf = NaN;  pairBig = NaN;  pairOk = NaN;  pairJump = NaN;  tauMed = NaN;
jumpClusters = zeros(0, 3);
if hasTruth
    DT_A = 0.25;
    vA = interp1(tTrueV, run.trueVel(:, 1:2), tSample + DT_A);  vB = interp1(tTrueV, run.trueVel(:, 1:2), tSample - DT_A);
    accH = (vA - vB) / (2 * DT_A);
    aMag = sqrt(sum(accH .^ 2, 2));
    dyn = aMag > 1.0 & all(isfinite(rowsErr(:, 1:2)), 2) & inWin;
    nDyn = nnz(dyn);
    tauPair(dyn) = -sum(rowsErr(dyn, 1:2) .* accH(dyn, :), 2) ./ aMag(dyn) .^ 2;
    if nDyn >= 10
        ta = abs(tauPair(dyn));
        pairOk = mean(ta < 0.10);  pairHalf = mean(ta >= 0.35 & ta < 0.65);  pairBig = mean(ta >= 0.65);
        pairShare = pairHalf + pairBig;  tauMed = median(tauPair(dyn));
        kd = find(dyn);  both = kd(1:end-1) + 1 == kd(2:end);
        if any(both)
            dTau = abs(tauPair(kd([both; false])) - tauPair(kd([false; both])));
            pairJump = mean(dTau > 0.35);
            tJump = t(kd([false; both]));  tJump = tJump(dTau > 0.35);
            if ~isempty(tJump)
                gaps = find(diff(tJump) > 5);                                % clusters of jumps, 5 s apart or less
                cStart = [1; gaps + 1];  cEnd = [gaps; numel(tJump)];
                jumpClusters = [tJump(cStart), tJump(cEnd), cEnd - cStart + 1];
            end
        end
    end
end
% the lever-arm error the rows carry against the truth point: the rows' error along and across the heading
% regressed on the yaw rate over the turning epochs, with an intercept; reported as truth minus rows per rad/s,
% which is leverArmCheck's convention (nav minus rows) when the truth is at the nav's point
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
    'fbGain', fbGain, 'fbLag', fbLag, 'fbScan', fbScan, 'yawRate', yawRate, 'armFit', armFit, 'armFitSig', armFitSig, ...
    'fbScanWide', fbScanWide, 'alignLag', alignLag, 'fbScanSwap', fbScanSwap, 'frameOff', frameOff, 'tauPair', tauPair, ...
    'pairShare', pairShare, 'pairJump', pairJump, 'jumpClusters', jumpClusters, 'opts', o);
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
if hasTruth && isfinite(alignLag)
    iW = find(fbScanWide(:, 1) == alignLag, 1);
    fprintf('   the same share over lags %d to %d epochs: sharpest %+.3f at lag %d (t %.1f)%s\n', WIDE(1), WIDE(end), fbScanWide(iW, 2), alignLag, ...
        fbScanWide(iW, 3), pick(o.tShift ~= 0, sprintf('; run shifted by %+.2f s', o.tShift), ''));
    if alignOff
        fprintf('      NOT ALIGNED: the feedback shows up %d epochs from where it belongs; the records and the run are not on the same clock,\n', alignLag - 1);
        fprintf('      give opts.tShift = %+.2f s (added to the run''s time)\n', o.tShift - (alignLag - 1) * median(diff(t)));
    end
end
if hasTruth && any(isfinite(fbScanSwap(:, 3)))
    [tSw, iSw] = max(abs(fbScanSwap(:, 3)));
    fprintf('   the share with north and east swapped: %+.3f at lag %d (t %.1f); records read as %s\n', fbScanSwap(iSw, 2), fbScanSwap(iSw, 1), tSw, o.frame);
    if frameOff
        fprintf('      FRAME: the feedback appears only with the axes swapped: the records'' state vector is in the other frame, pass opts.frame = ''%s''\n', ...
            pick(isENU, 'NED', 'ENU'));
    end
end
if hasAtt && hasTruth
    if nTurn >= 10
        armTxt = '';
        if ~isempty(o.armSlope), armTxt = sprintf(', after the arm correction [%+.2f %+.2f]', o.armSlope(1), o.armSlope(2)); end
        fprintf('   rows minus truth per rad/s of yaw rate over %d turning epochs%s: along %+.2f m, across %+.2f m (sigma %.2f / %.2f)\n', ...
            nTurn, armTxt, armFit(1), armFit(2), armFitSig(1), armFitSig(2));
        if ~isempty(o.truthSlope)
            fprintf('      (truth point moved by [%+.2f %+.2f] m per rad/s to the rows'' point)\n', o.truthSlope(1), o.truthSlope(2));
        end
        if any(abs(armFit) > 2 * armFitSig & abs(armFit) > 0.1)
            if isempty(o.armSlope)
                fprintf('      (the rows carry a lever-arm error against the truth point; opts.armSlope = leverArmCheck''s slopes emulates a corrected arm)\n');
            else
                fprintf('      (what is left is the truth point''s offset from the rows'' point; opts.truthSlope = these two numbers takes it out of the scoring)\n');
            end
        end
    else
        fprintf('   rows minus truth per rad/s of yaw rate: fewer than 10 turning epochs, no arm read\n');
    end
end
if hasTruth && nDyn >= 10
    fprintf('   time pairing on %d epochs faster than 1 m/s2, the rows'' error as a time offset: under 0.10 s %.0f%%, 0.35 to 0.65 s %.0f%%, over 0.65 s %.0f%%; median %+.2f s; jumps over 0.35 s between consecutive epochs %.0f%%\n', ...
        nDyn, 100 * pairOk, 100 * pairHalf, 100 * pairBig, tauMed, 100 * pairJump);
    if pairJump > 0.2
        fprintf('      PAIRING: the offset jumps between consecutive epochs on %.0f%% of them: rows paired with the wrong pulse on some epochs, in the\n', 100 * pairJump);
        fprintf('      receiver''s time tag or in the host; a screen cannot fix a time tag, and these epochs read as episodes and flags below\n');
    end
    if ~isempty(jumpClusters)
        [~, order] = sort(jumpClusters(:, 3), 'descend');
        nShow = min(numel(order), 6);
        fprintf('      the jumps cluster at:');
        for i = 1:nShow
            fprintf(' %.0f-%.0f s (%d)', jumpClusters(order(i), 1), jumpClusters(order(i), 2), jumpClusters(order(i), 3));
            if i < nShow, fprintf(','); end
        end
        fprintf('\n');
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
        fprintf('   per epoch, %.0f to %.0f s: rows - truth, the onset statistic against its threshold, the flag, the separation from the coast,\n', tt(1), tt(2));
        fprintf('   and on the epochs faster than 1 m/s2 the rows'' error as a time offset\n');
        fprintf('        t      rows-truth N  E   |.|     onset  thresh   flag    sep    offset\n');
        for k = find(t >= tt(1) & t <= tt(2))'
            tauTxt = '      ';
            if isfinite(tauPair(k)), tauTxt = sprintf('%+6.2f', tauPair(k)); end
            fprintf('   %8.2f    %+6.3f %+6.3f  %5.3f    %5.3f   %5.3f    %s   %5.3f   %s\n', t(k), rowsErr(k, 1), rowsErr(k, 2), errH(k), ...
                onset(k), thresh(k), pick(flag(k), 'FLAG', '    '), sep(k), tauTxt);
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
function [v] = toNED(v, isENU)
% a [n x 3] vector set in the records' frame, returned in NED: ENU exchanges the first two and flips the third
if isENU, v = [v(:, 2), v(:, 1), -v(:, 3)]; end
end

function [out] = shareAtLag(xPlus, dNav, dTrue, inWin, L, K)
% the share of the state x+(k - L) in (nav change - truth change)(k) over the quiet epochs, both axes stacked,
% through the origin, with one trim of the outliers: [share, t-statistic]; NaN with too few epochs
out = NaN(1, 2);
kk = (max(1, 1 + L):min(K, K + L))';  kk = kk(kk >= 2);
a = xPlus(kk - L, 1:2);  b = dNav(kk, 1:2) - dTrue(kk, 1:2);
ok = all(isfinite([a b]), 2) & inWin(kk) & sqrt(sum(dTrue(kk, 1:2) .^ 2, 2)) < 0.5;
av = a(ok, :);  bv = b(ok, :);  av = av(:);  bv = bv(:);
if numel(av) < 10 || av' * av <= 0, return; end
c = -(av' * bv) / (av' * av);
resid = bv + c * av;
% the scale floored on the rms so that a residual that is mostly zero does not trim everything but the zeros
madS = max(1.4826 * median(abs(resid - median(resid))), 0.2 * sqrt(mean(resid .^ 2)));
if madS > 0
    keep = abs(resid) < 4 * madS;
    if nnz(keep) >= 10 && av(keep)' * av(keep) > 0
        av = av(keep);  bv = bv(keep);
        c = -(av' * bv) / (av' * av);
        resid = bv + c * av;
    end
end
se = sqrt((resid' * resid) / max(numel(av) - 1, 1) / (av' * av));
out = [c, c / max(se, eps)];
end

function [r] = runsOf(mask)
% start and end indices of the runs of true in a logical vector, [n x 2]
mask = mask(:)';
d = diff([false, mask, false]);
r = [find(d == 1)', find(d == -1)' - 1];
end

function [out] = pick(cond, a, b)
if cond, out = a; else, out = b; end
end
