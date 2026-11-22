function [res] = dopplerEpisodeCheck(spfMeasLog, tKf, run, opts)
%DOPPLEREPISODECHECK  The Doppler episode screen emulated offline on the records and the run, read and scored.
%   res = dopplerEpisodeCheck(spfMeasLog, tKf, run)
%   res = dopplerEpisodeCheck(spfMeasLog, tKf, run, opts)
%
%   WHAT IT ANSWERS
%   The receiver's Doppler goes wrong in episodes of strong transient dynamics (the vehicle data: 0.4 m/s for
%   20 s at a tunnel portal, all rows together, so no consensus screen sees it), and a filter following its rows
%   carries the error into the states. SPF_dopplerScreen holds such rows out of the update; this tool runs the
%   same logic on the spfMeas records and the run, reads the inputs it is fed, and, with the truth, scores it:
%   episodes caught, release tails, false flags. It also scores the host screen's own log (opts.hostLog).
%
%   THE SCREEN (as SPF_dopplerScreen, emulated epoch by epoch)
%     rows' velocity  v_rows(k) = v_hat(k) - dv(k): the velocity the host measured the rows against minus the
%                     least-squares fit of the rate-row innovations to a velocity and one clock drift per
%                     constellation, fitted only when the rows outnumber the unknowns as the host's fit does
%                     (predicted minus measured; opts.sign = +1 for the other convention). The host's innovation
%                     is total-state, predicted minus measured minus H times the prior error state, so with the
%                     prior in the records (priorState) v_hat is the nav velocity before the update less the
%                     prior state's velocity; without it, the nav velocity alone
%     increment       with the prior: dvIns(k) = v_hat(k) - v_hat(k-1) + (x+ - x-)(k-1), the change of v_hat with
%                     the previous update's own step added back, the IMU's increment less the error growth the
%                     model expects; a hand-over moves the nav and the state alike and leaves v_hat where it was,
%                     so nothing is fitted (fbModel 'total', the default with a prior)
%     coast           v_c follows the increment every epoch; each accepted epoch re-anchors it on the rows through
%                     a first-order filter ('alpha' of the difference); the tolerance T = sqrt(floor^2 +
%                     (driftRate * t_anchor)^2 + (headingSigma * sum|dvIns| since the anchor)^2) grows with the
%                     time the coast has run on the IMU alone
%     flag            'confirm' consecutive epochs with |v_rows - v_c| > T start an episode; 'release' consecutive
%                     epochs back within T end it, plus one per 'releaseGrow' epochs held, the coast kept; the hold
%                     is time, judged epochs or not; at 'hold' epochs the episode is closed and the coast re-anchored
%                     on v_hat, the filter's own velocity, which in the host has run without the range-rate rows
%                     for the whole hold and carries the ranges' aiding (on records from a host that used the rows,
%                     v_hat follows them, so the emulation shows that rule's best case); without the prior the
%                     coast is kept; an epoch with too few rows to judge neither confirms nor releases, and stays
%                     flagged while an episode is on
%     blind           after 'hold' consecutive epochs without rows the coast is lost and the next rows anchor it
%                     afresh, so a tunnel does not leave a coast metres per second off at its exit
%     lever arm       opts.armSlope = [along across] (m per rad/s of yaw rate) takes a lever-arm error out of the
%                     rows' velocity, as a corrected arm input would; opts.truthSlope moves the truth to the rows'
%                     point for the scoring; a call without them names both vectors
%
%   THE READS ON THE INPUTS (printed first, since the screen is only as good as what it is fed)
%     increment       with the prior and the truth: the total-state increment against the truth over the quiet
%                     epochs, its rms and the share of the state left in it, 0 when the mechanisation takes the
%                     whole hand-over and the prior is the state the rows were predicted with; trimmed, so sparse
%                     jumps escape it, and the next read is for them
%     state loss      what the state lost between an update and the next prior, as a share of that update's
%                     posterior (horizontal), classed per epoch (negative, none, a tenth, a quarter, a half, three
%                     quarters, all); the increment error on the loss per class says how much of it the nav took
%                     (slope 0: all; slope s: a share 1 - s); a positive share lost and not taken is a step in nav
%                     minus prior the increment reads as motion, flagged LOSS NOT TAKEN
%     two coasts      with the prior and run.fbVel: the total-state increment and the nav's change plus the logged
%                     correction, each against the truth; they differ by the error growth the model expects
%     lever arm       nav minus rows per rad/s of yaw rate over the turning epochs, along and across, the host's
%                     own arm error (leverArmCheck's read); and rows minus truth, the rows' error against the truth
%                     point, before or after the corrections
%     arm or latency  an arm and a pairing latency read the same per rad/s at one speed, since the latency's error
%                     in a turn is the centripetal acceleration, speed times yaw rate, times the latency; they
%                     part with speed. Nav minus rows is regressed jointly on the yaw rate (an arm, m) and on
%                     speed times yaw rate (a latency, s), and the across slope per rad/s is given in a slow and
%                     a fast band; the filter dilutes a consistent row error within a turn, so the host's view is
%                     a lower bound, and with the truth the same fit is made on truth minus rows and decides;
%                     turns all at one speed cannot tell the two apart and the read says so
%     time pairing    on the epochs faster than 1 m/s2 the rows' error against the truth projected on the
%                     acceleration is a time offset; an offset jumping between consecutive epochs is rows paired
%                     with the wrong pulse on some of them, which no screen can fix
%     coast           the coast's own error against the truth and its drift rate over every flagged stretch of
%                     10 epochs or more, the IMU drift the tolerance has to cover, against 'driftRate'; and the
%                     times the coast was lost
%     legacy          under the fallback models, or with opts.verbose: the feedback share of the state at lags 0
%                     to 2 and over lags -12 to 12 (a share far from its place is the run's clock off the records',
%                     opts.tShift), the same with north and east swapped (a share only there is a state vector in
%                     the other frame, opts.frame), and the delivery profile of the update's increment or of the
%                     logged correction; under the total-state models these are notes at most
%
%   THE FALLBACK MODELS (records without the prior)
%     fbModel 'state': the increment is the nav's change minus fbGain times the state at fbLag (fitted against
%     the truth, or given); 'profile': minus the update's increment delivered with a fitted profile; 'logged':
%     minus run.fbVel as delivered (fitted, or fbLagWeights); 'handover': the total-state rows' velocity with the
%     logged correction's increment, for the two-coasts comparison. The default is 'total' with the prior, else
%     'logged' with run.fbVel, else the one of 'state' and 'profile' implying the larger feedback.
%
%   SCORING (needs the truth)
%   The rows' error v_rows - v_truth labels the epochs the receiver was wrong (above 'episodeLevel' for at least
%   'episodeMin' epochs); each episode is listed with its detection delay and coverage; release tails are flags
%   within the longest release after an episode with the rows back within 'cleanLevel' of the truth; false flags
%   are flags on clean rows that are not a tail, grouped into runs, each with the epochs without rows in the 30 s
%   before it and the coast's error at its start; a per-epoch table around the largest episode (or opts.tTable)
%   gives rows minus truth, the separation against the tolerance, the flag, the time offset, the increment's
%   error, the state's loss and the coast's error.
%
%   INPUTS
%     spfMeasLog  the records, one per epoch: innovation, obsMatrix, numMeas, priorState, postState
%     tKf         [K x 1] the epochs' host time
%     run         t, navVel (100 Hz, NED m/s); tTrue, trueVel (NED, m/s) for the scoring; navAtt ([roll pitch yaw],
%                 rad) for the arm reads; fbVel (the velocity correction the host took off the nav, NED m/s,
%                 [K x 3] per epoch or [N x 3] per sample) for the fallback and hand-over models
%
%   OPTIONS (defaults)
%     the screen   floor 0.20 m/s, driftRate 0.005 m/s2, headingSigma 5e-3 rad, alpha 0.10, confirm 2, release 4,
%                  releaseGrow 10, hold 40 (epochs), sign -1
%     the model    fbModel '' (see above), fbGain [] (fitted; a number is used as given), fbLag [] (0, 1 or 2;
%                  default the lag the fit finds, 1 without truth), fbLagWeights [] (a per-epoch run.fbVel's
%                  delivery given instead of fitted, weights by lag from 0, e.g. [0 1])
%     corrections  armSlope [] (m per rad/s, along and across), truthSlope [], tShift 0 s (added to the run's time
%                  to put it on the records' clock), frame 'NED' ('ENU' maps the records' components to NED)
%     scoring      episodeLevel 0.25 m/s, episodeMin 4 epochs, cleanLevel 0.15 m/s, tWin [] (the epochs to run
%                  on), tTable [] (the table's window), hostLog [] ([K x 5]: flag, separation, tolerance, hold,
%                  rate rows per epoch from the host screen; given, the host's flags are scored in place of the
%                  emulation's and its separation and tolerance fill the table)
%     output       quiet false, verbose false
%
%   OUTPUT res
%     the screen   t, vRows, dvIns [K x 3], sep, onset, thresh [K x 1], flag [K x 1] logical, vCoast [K x 3],
%                  nFlag, coastLost ([n x 1], the times the coast was lost)
%     the scoring  rowsErr [K x 3], episodes (tStart, tEnd, detected, delay, coverage), falseRuns (tStart, tEnd,
%                  epochs, blindBefore, coastErrStart), nFalse, nTail, coastErr [K x 1], coastDrift ([n x 3]: start,
%                  end, m/s per s)
%     the reads    fbModel, fbGain, fbLag (NaN under 'total'); totShare, totShareT, totRms, incErr [K x 1];
%                  lossRatio [K x 1], lossClass (name, range, n, nUsed, slope, rms), lossNotTaken; hoRms; yawRate
%                  [K x 1], armHost, armHostSig, armFit, armFitSig ([1 x 2], along and across, m per rad/s);
%                  armJoint, armJointSig, latFit, latFitSig, bandSlope ([2 x 2]), bandSpeed ([1 x 2]) and the same
%                  against the truth as armTrueJoint, armTrueJointSig, latTrue, latTrueSig, bandSlopeTrue,
%                  bandSpeedTrue; tauPair [K x 1], pairShare, pairJump, jumpClusters ([n x 3]); the legacy reads
%                  fbScan ([3 x 4]), fbScanWide ([25 x 4]), alignLag, fbScanSwap, frameOff, firC, firSe, firTotal,
%                  firTotalSe, lgTotal, lgC; opts
%   MATLAB desktop tool, not codegen. Runs under Octave.

if iscell(spfMeasLog), spfMeasLog = [spfMeasLog{:}]; end
if nargin < 4 || isempty(opts), opts = struct(); end
o = struct('floor', 0.20, 'driftRate', 0.005, 'headingSigma', 5e-3, 'alpha', 0.10, 'confirm', 2, 'release', 4, 'releaseGrow', 10, 'hold', 40, 'sign', -1, 'fbGain', [], ...
    'fbLag', [], 'fbModel', '', 'fbLagWeights', [], 'armSlope', [], 'truthSlope', [], 'tShift', 0, 'episodeLevel', 0.25, 'episodeMin', 4, 'cleanLevel', 0.15, ...
    'tWin', [], 'tTable', [], 'quiet', false, 'frame', 'NED', 'hostLog', [], 'verbose', false);
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
vNavM = NaN(K, 3);  dv = NaN(K, 3);  xPlus = NaN(K, 3);  xPri = NaN(K, 3);  nRate = zeros(K, 1);  tSample = t;
hasPrior = isfield(spfMeasLog, 'priorState');
for k = 1:K
    j = find(tNav < t(k), 1, 'last');
    if isempty(j), continue; end
    vNavM(k, :) = navVel(j, :);  tSample(k) = tNav(j);               % the rows and the truth are read at this sample
    s = spfMeasLog(k);  m = double(s.numMeas);
    if hasPrior && ~isempty(s.priorState), xPri(k, :) = toNED(s.priorState(4:6)', isENU); end   % the prior, update or not
    if m == 0, continue; end
    xPlus(k, :) = toNED(s.postState(4:6)', isENU);
    H = s.obsMatrix(1:m, :);  y = s.innovation(1:m);  y = y(:);
    rate = sqrt(sum(H(:, 4:6) .^ 2, 2)) > 0.5;
    if nnz(rate) < 4, continue; end
    cols = [4 5 6, 6 + find(any(H(rate, 7:end) ~= 0, 1))];
    nRate(k) = nnz(rate);                                                 % rate rows present, judged or not
    if nnz(rate) <= numel(cols), continue; end                            % the host's rule: more rows than unknowns, else unjudged
    sol = H(rate, cols) \ y(rate);
    dv(k, :) = toNED(sol(1:3)', isENU);
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
fbScan = [LAGS(:), NaN(numel(LAGS), 3)];  fbScanWide = [WIDE(:), NaN(numel(WIDE), 3)];
if hasTruth
    for iL = 1:numel(WIDE)
        fbScanWide(iL, 2:4) = shareAtLag(xPlus, dNav, dTrue, inWin, WIDE(iL), K);
    end
    for iL = 1:numel(LAGS)
        fbScan(iL, 2:4) = fbScanWide(WIDE == LAGS(iL), 2:4);
    end
end
% the same share with north and east swapped: a share that appears only then is a state vector in the other
% frame (ENU against NED), the first two components exchanged
fbScanSwap = [LAGS(:), NaN(numel(LAGS), 3)];
if hasTruth
    xSwap = xPlus(:, [2 1 3]);
    for iL = 1:numel(LAGS)
        fbScanSwap(iL, 2:4) = shareAtLag(xSwap, dNav, dTrue, inWin, LAGS(iL), K);
    end
end
% and the update's own increment, x+ minus x-, when the records carry the prior, delivered to the nav over the
% following epochs: a joint fit of (nav change - truth change)(k) on the increments at lags 0 to FIR_LAGS, both
% axes stacked, over the quiet epochs, with one trim; a mechanisation that takes the whole Kalman step at once
% reads 1 at one lag, one that slews it in reads a profile whose total is the share delivered
FIR_LAGS = 24;
hasInc = hasPrior && any(isfinite(xPri(:, 1)));
xInc = xPlus - xPri;                                                       % NaN without a prior or without an update
incZ = xInc;  incZ(~isfinite(incZ)) = 0;                                  % no update, no increment
% the total-state quantities: v_hat, the velocity the host measured the rows against (the innovation is predicted
% minus measured minus H x-), and its change with the previous update's own step added back, which is the IMU's
% increment less the error growth the model expects; a hand-over moves the nav and the state alike and leaves
% v_hat where it was, so the feedback never enters
vHat = vNavM - xPri;
dHat = [NaN(1, 3); diff(vHat)];
dvTotal = dHat + [zeros(1, 3); incZ(1:K - 1, :)];
% against the truth: the rms of that increment's error over the quiet epochs, and the share of the state left in it
% (dvTotal - dTrue = share * x+(k - 1)), 0 when the mechanisation takes the whole hand-over and the prior is the
% state the rows were predicted with
totShare = NaN;  totShareT = NaN;  totRms = NaN;
if hasInc && hasTruth
    sc = shareAtLag(xPlus, dvTotal, dTrue, inWin, 1, K);
    totShare = -sc(1);  totShareT = -sc(2);
    rT = dvTotal(:, 1:2) - dTrue(:, 1:2);
    qT = all(isfinite(rT), 2) & inWin & sqrt(sum(dTrue(:, 1:2) .^ 2, 2)) < 0.5;
    if nnz(qT) >= 10, rr = rT(qT, :);  totRms = sqrt(mean(rr(:) .^ 2)); end
end
% the state's loss between an update and the next prior, x+(k - 1) - x-(k), as a share of x+(k - 1) (horizontal, for a
% posterior above 0.03 m/s), classed per epoch; a feedback that hands over a share of the state and takes the same off
% it reads that share, and the nav took what the state lost; the increment error dvTotal - dTrue regressed on the loss
% per class says how much of the loss the nav took: slope 0 all of it, slope s a share 1 - s, and the rest is a step
% in nav minus prior that the total-state increment cannot see
incErr = sqrt(sum((dvTotal(:, 1:2) - dTrue(:, 1:2)) .^ 2, 2));
lossH = NaN(K, 2);  lossRatio = NaN(K, 1);
for k = 2:K
    xp = xPlus(k - 1, 1:2);
    if all(isfinite(xp)) && all(isfinite(xPri(k, 1:2)))
        lossH(k, :) = xp - xPri(k, 1:2);
        if norm(xp) > 0.03, lossRatio(k) = (lossH(k, :) * xp') / (xp * xp'); end
    end
end
LOSS_CLASSES = {'negative', -Inf, -0.03; 'none', -0.03, 0.03; 'a tenth', 0.03, 0.17; 'a quarter', 0.17, 0.4; 'a half', 0.4, 0.6; ...
    'three quarters', 0.6, 0.9; 'all', 0.9, 1.1; 'other', 1.1, Inf};
lossClass = struct('name', {}, 'range', {}, 'n', {}, 'nUsed', {}, 'slope', {}, 'rms', {});
lossNotTaken = false;
if hasInc
    % the loss at such an epoch is large against the truth's own error in dynamics, so the read takes every epoch
    % under 3 m/s2, not only the quiet ones
    usableL = isfinite(incErr) & inWin & sqrt(sum(dTrue(:, 1:2) .^ 2, 2)) < 1.5;
    for i = 1:size(LOSS_CLASSES, 1)
        inC = inWin & isfinite(lossRatio) & lossRatio >= LOSS_CLASSES{i, 2} & lossRatio < LOSS_CLASSES{i, 3};
        c = struct('name', LOSS_CLASSES{i, 1}, 'range', [LOSS_CLASSES{i, 2}, LOSS_CLASSES{i, 3}], 'n', nnz(inC), 'nUsed', 0, 'slope', NaN, 'rms', NaN);
        use = inC & usableL;
        c.nUsed = nnz(use);
        if hasTruth && c.nUsed >= 5
            eL = dvTotal(use, 1:2) - dTrue(use, 1:2);  jL = lossH(use, :);
            if jL(:)' * jL(:) > 0, c.slope = (eL(:)' * jL(:)) / (jL(:)' * jL(:)); end
            c.rms = sqrt(mean(eL(:) .^ 2));
            if c.nUsed >= 10 && c.slope > 0.2 && LOSS_CLASSES{i, 2} >= 0.03, lossNotTaken = true; end   % a positive share lost and not taken
        end
        lossClass(end + 1) = c; %#ok<AGROW>
    end
end
[firC, firSe, firTotal, firTotalSe, firOk, firApplied, firLag0Dropped] = firFit(incZ, dNav, dTrue, inWin, K, FIR_LAGS, T_MIN, hasTruth && hasInc);
% the velocity correction the host logged per epoch (run.fbVel, [K x 3] NED m/s, the amount taken off the nav
% velocity after that epoch's update) or per nav sample ([N x 3], at the sample it was applied): the first is
% fitted for its delivery profile like the increment, the second is summed into the intervals and needs no fit
hasFb = isfield(run, 'fbVel') && ~isempty(run.fbVel);
fbLogged = zeros(K, 3);  fbPerSample = false;
lgC = NaN(FIR_LAGS + 1, 1);  lgSe = NaN(FIR_LAGS + 1, 1);  lgTotal = NaN;  lgTotalSe = NaN;  lgOk = false;  lgApplied = zeros(K, 3);
lgLag0Dropped = false;
if hasFb
    fbV = toNED(run.fbVel, isENU);
    if size(fbV, 1) == numel(tNav)
        fbPerSample = true;
        for k = 2:K
            inInt = tNav > tSample(k - 1) & tNav <= tSample(k);
            fbLogged(k, :) = sum(fbV(inInt, :), 1);
        end
        lgApplied = fbLogged;  lgOk = true;
        if hasTruth
            lgScan = shareAtLag(fbLogged, dNav, dTrue, inWin, 0, K);          % must read 1 at lag 0
            lgTotal = lgScan(1);  lgTotalSe = abs(lgScan(1)) / max(abs(lgScan(2)), eps);
        end
    elseif size(fbV, 1) == K
        fbLogged = fbV;  fbLogged(~isfinite(fbLogged)) = 0;
        [lgC, lgSe, lgTotal, lgTotalSe, lgOk, lgApplied, lgLag0Dropped] = firFit(fbLogged, dNav, dTrue, inWin, K, FIR_LAGS, T_MIN, hasTruth);
        if ~isempty(o.fbLagWeights)
            % the delivery given, not fitted: weights per lag from 0, e.g. [0 1] for the whole correction in the
            % next interval; the fitted profile stays in the printout as the check
            w = o.fbLagWeights(:);
            lgApplied = zeros(K, 3);
            for j = 0:numel(w) - 1
                lgApplied(j + 1:K, :) = lgApplied(j + 1:K, :) + w(j + 1) * fbLogged(1:K - j, :);
            end
            lgOk = true;
        elseif ~hasTruth
            lgApplied = [zeros(1, 3); fbLogged(1:K - 1, :)];  lgOk = true;    % without truth: applied in the next interval
        end
    else
        error('dopplerEpisodeCheck:fbVel', 'run.fbVel must have one row per epoch or one per nav sample');
    end
end
frameOff = false;
if any(isfinite(fbScanSwap(:, 3))) && any(isfinite(fbScan(:, 3)))
    frameOff = max(abs(fbScanSwap(:, 3))) >= T_MIN && max(abs(fbScanSwap(:, 3))) > 1.5 * max(abs(fbScan(:, 3)));
end

% the model of the feedback: a share of the state at one lag, or the update's increment delivered over the
% following epochs, whichever implies the larger feedback in the logged velocity; opts.fbModel forces one
scanOk = any(isfinite(fbScan(:, 3)));
if scanOk, [tBest, iBest] = max(abs(fbScan(:, 3))); else, tBest = NaN;  iBest = 2; end
fbFound = isfinite(tBest) && tBest >= T_MIN;
impliedState = 0;  impliedFir = 0;                                        % the feedback each model implies, m/s per epoch rms
if fbFound, impliedState = abs(fbScan(iBest, 2)) * fbScan(iBest, 4); end
if firOk, fa = firApplied(isfinite(xInc(:, 1)) & inWin, 1:2);  impliedFir = sqrt(mean(fa(:) .^ 2)); end
if ~isempty(o.fbModel)
    fbModel = o.fbModel;
elseif hasInc
    fbModel = 'total';                                                     % the prior in the records: nothing to fit
elseif hasFb && lgOk
    fbModel = 'logged';
elseif firOk && (~fbFound || impliedFir > impliedState)
    fbModel = 'profile';
else
    fbModel = 'state';
end
if ~any(strcmp(fbModel, {'total', 'handover', 'state', 'profile', 'logged'})), error('dopplerEpisodeCheck:fbModel', 'opts.fbModel is ''total'', ''handover'', ''state'', ''profile'' or ''logged'''); end
if any(strcmp(fbModel, {'logged', 'handover'})) && ~hasFb, error('dopplerEpisodeCheck:fbModel', 'opts.fbModel ''%s'' needs run.fbVel', fbModel); end
if any(strcmp(fbModel, {'total', 'handover'})) && ~hasInc, error('dopplerEpisodeCheck:fbModel', 'opts.fbModel ''%s'' needs priorState in the records', fbModel); end
useTotal = strcmp(fbModel, 'total');  useHand = strcmp(fbModel, 'handover');  useFir = strcmp(fbModel, 'profile');  useLog = strcmp(fbModel, 'logged');
found = fbFound;
if useFir, found = firOk; end
if useLog || useHand, found = lgOk; end
if useTotal, found = true; end
if isempty(o.fbGain)
    if ~hasTruth && ~useLog && ~useTotal && ~useHand
        error('dopplerEpisodeCheck:fbGain', 'give opts.fbGain (the share of the state the nav takes per epoch, signed) when the run has no truth');
    end
    if useTotal
        fbGain = NaN;  fbLag = NaN;                                    % nothing fitted: the prior carries the feedback
    elseif useHand
        fbGain = lgTotal;  fbLag = 1;                                  % the logged correction, as delivered
    elseif useLog
        fbGain = lgTotal;  fbLag = 1;                                  % the logged correction, as delivered
    elseif useFir
        fbGain = firTotal;  fbLag = 1;                                 % the total delivered; the profile itself is applied
    elseif fbFound
        fbGain = fbScan(iBest, 2);  fbLag = fbScan(iBest, 1);
    else
        fbGain = 0;  fbLag = 1;                                       % no feedback in the logged velocity under either model
    end
else
    fbGain = o.fbGain;  useTotal = false;  useHand = false;  useFir = false;  useLog = false;  fbModel = 'state';
    if ~isempty(o.fbLag), fbLag = o.fbLag;
    elseif fbFound, fbLag = fbScan(iBest, 1);
    else, fbLag = 1; end
end
if useTotal
    dvIns = dvTotal;
elseif useLog || useHand
    dvIns = dNav + lgApplied;
elseif useFir
    dvIns = dNav + firApplied;
else
    xLag = NaN(K, 3);  xLag(fbLag + 1:K, :) = xPlus(1:K - fbLag, :);
    xLag(~isfinite(xLag)) = 0;                                            % no update, no feedback
    dvIns = dNav + fbGain * xLag;
end

% the hand-over increment, the nav's change plus the logged correction as delivered, against the truth over the quiet
% epochs: beside the total-state one it says which coast drifts less, the two differing by the error growth the
% filter's model expects between updates
hoRms = NaN;
if hasInc && hasFb && lgOk && hasTruth
    rH = (dNav(:, 1:2) + lgApplied(:, 1:2)) - dTrue(:, 1:2);
    qH = all(isfinite(rH), 2) & inWin & sqrt(sum(dTrue(:, 1:2) .^ 2, 2)) < 0.5;
    if nnz(qH) >= 10, rr = rH(qH, :);  hoRms = sqrt(mean(rr(:) .^ 2)); end
end

% the rows' velocity: the nav velocity the host predicted with, before that epoch's feedback, minus the fitted
% innovation velocity (predicted minus measured: the rows sit at nav - dv), then the lever-arm error taken out,
% which puts the rows at the nav's point
vNavPre = vNavM;
if useTotal || useHand
    vNavPre = vHat;                                                      % the velocity the host measured the rows against
elseif useLog && ~fbPerSample && ~isempty(o.fbLagWeights)
    vNavPre = vNavM + o.fbLagWeights(1) * fbLogged;                      % the part of this epoch's correction already in the sample
elseif useLog && ~fbPerSample && isfinite(lgC(1))
    vNavPre = vNavM + lgC(1) * fbLogged;
elseif useFir
    vNavPre = vNavM + firC(1) * incZ;                                    % the part of this epoch's step already in the sample
elseif fbLag == 0
    xNow = xPlus;  xNow(~isfinite(xNow)) = 0;
    vNavPre = vNavM + fbGain * xNow;
end
vRows = vNavPre + o.sign * dv;
vRowsRaw = vRows;
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

% the emulated screen: an INS coast anchored on the accepted rows and propagated with the INS increment; the
% rows flagged while they sit further from it than a tolerance that grows with the time the coast has run on
% the INS alone, released when they agree with it again
onset = NaN(K, 1);  thresh = NaN(K, 1);  flag = false(K, 1);  sep = NaN(K, 1);  vCoastK = NaN(K, 3);  coastLost = zeros(0, 1);
inEp = false;  anchored = false;  vCoast = zeros(1, 3);  hold = 0;  rel = 0;  above = 0;  blind = 0;
tAnchor = 0;  speedSince = 0;  dtEp = median(diff(t));
for k = 1:K
    if ~inWin(k), continue; end
    haveIns = all(isfinite(dvIns(k, :)));
    haveRows = all(isfinite(vRows(k, :)));
    if anchored && haveIns
        vCoast = vCoast + dvIns(k, :);                                        % the coast follows the INS
        tAnchor = tAnchor + dtEp;  speedSince = speedSince + norm(dvIns(k, 1:2));
    elseif anchored
        anchored = false;  inEp = false;  above = 0;                          % no INS increment: the coast is lost
    end
    if anchored, vCoastK(k, :) = vCoast; end                                  % the coast as propagated, rows or not
    if ~haveRows
        if inEp && nRate(k) > 0, flag(k) = true; end                          % too few rate rows to judge: an episode in progress stands
        if inEp
            hold = hold + 1;                                                  % the hold is time, judged or not
            if hold >= o.hold
                inEp = false;  rel = 0;  hold = 0;  flag(k) = false;          % the hold up on an unjudged epoch: the same re-anchor
                if (useTotal || useHand) && all(isfinite(vHat(k, :))), vCoast = vHat(k, :); end
            end
        end
        if anchored
            blind = blind + 1;
            if blind > o.hold                                                 % blind for longer than the IMU alone can carry it
                anchored = false;  inEp = false;  above = 0;  blind = 0;
                coastLost(end + 1, 1) = t(k); %#ok<AGROW>
            end
        end
        continue;
    end
    blind = 0;
    if ~anchored
        vCoast = vRows(k, :);  anchored = true;  tAnchor = 0;  speedSince = 0;   % the first rows anchor the coast
        vCoastK(k, :) = vCoast;
        continue;
    end
    sep(k) = norm(vRows(k, 1:2) - vCoast(1:2));
    thresh(k) = sqrt(o.floor ^ 2 + (o.driftRate * tAnchor) ^ 2 + (o.headingSigma * speedSince) ^ 2);
    onset(k) = sep(k);
    if ~inEp
        if sep(k) > thresh(k)
            above = above + 1;                                                % not accepted: the coast keeps to the INS
            if above >= o.confirm
                inEp = true;  hold = 0;  rel = 0;  above = 0;  flag(k) = true;
            end
        else
            above = 0;
            vCoast = vCoast + o.alpha * (vRows(k, :) - vCoast);               % accepted: the coast re-anchored, filtered
            tAnchor = 0;  speedSince = 0;
        end
    else
        hold = hold + 1;  flag(k) = true;
        if sep(k) <= thresh(k), rel = rel + 1; else, rel = 0; end
        if rel >= o.release + floor(hold / o.releaseGrow)                     % the longer the rows were wrong, the longer they must agree
            inEp = false;  rel = 0;  hold = 0;                                % released: the coast stays, the rows must keep agreeing
            vCoast = vCoast + o.alpha * (vRows(k, :) - vCoast);
            tAnchor = 0;  speedSince = 0;
        elseif hold >= o.hold
            inEp = false;  rel = 0;  hold = 0;                                % the hold is up: the coast re-anchored on the filter's own
            if (useTotal || useHand) && all(isfinite(vHat(k, :)))            % velocity, its tolerance kept as grown (the total-state
                vCoast = vHat(k, :);                                          % forms; the others keep the coast)
            end
        end
    end
end
% the host's own log scored in place of the emulation: its flags, separation and tolerance per epoch
hostScored = ~isempty(o.hostLog);
if hostScored
    if size(o.hostLog, 1) ~= K || size(o.hostLog, 2) < 3
        error('dopplerEpisodeCheck:hostLog', 'opts.hostLog must have one row per epoch (%d) and at least 3 columns, has %d x %d', K, size(o.hostLog, 1), size(o.hostLog, 2));
    end
    flag = logical(o.hostLog(:, 1)) & inWin;
    sep = o.hostLog(:, 2);  thresh = o.hostLog(:, 3);
    sep(~inWin | nRate == 0) = NaN;  thresh(~inWin | nRate == 0) = NaN;
end

% the coast's own error against the truth, and its growth over every flagged stretch of 10 epochs or more: the IMU
% drift the tolerance has to cover, to set against DRIFT_RATE
coastErr = sqrt(sum((vCoastK(:, 1:2) - vTrue(:, 1:2)) .^ 2, 2));
coastDrift = zeros(0, 3);
if hasTruth
    runsFl = runsOf(flag & isfinite(coastErr));
    for i = 1:size(runsFl, 1)
        a = runsFl(i, 1);  b = runsFl(i, 2);
        if b - a + 1 >= 10 && t(b) > t(a)
            coastDrift(end + 1, :) = [t(a), t(b), (coastErr(b) - coastErr(a)) / (t(b) - t(a))]; %#ok<AGROW>
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
% and the host's own lever-arm error, nav minus the rows as recorded per rad/s, leverArmCheck's read in the
% chosen frame; with neither correction given, the two vectors to pass follow: armSlope the host's arm error,
% truthSlope the truth's read-back less it
armFit = NaN(1, 2);  armFitSig = NaN(1, 2);  armHost = NaN(1, 2);  armHostSig = NaN(1, 2);  nTurn = 0;
if hasAtt
    navMinusRows = vNavPre(:, 1:2) - vRowsRaw(:, 1:2);
    if useTotal || useHand, navMinusRows = vNavM(:, 1:2) - vRowsRaw(:, 1:2); end   % the nav itself, not nav less the prior: the state absorbs part of the arm in a corner
    turnH = abs(yawRate) > 2 * pi / 180 & all(isfinite(navMinusRows), 2) & inWin;
    if nnz(turnH) >= 10
        [armHost, armHostSig] = slopeOnYaw(navMinusRows, yawRate, uAlong, uAcross, turnH);
    end
end
% an arm and a latency part with speed: nav minus rows on the yaw rate and on speed times yaw rate jointly, and the
% slope per rad/s in a slow and a fast band
armJoint = NaN(1, 2);  armJointSig = NaN(1, 2);  latFit = NaN(1, 2);  latFitSig = NaN(1, 2);  bandSlope = NaN(2, 2);  bandSpeed = NaN(1, 2);
speedK = sqrt(sum(vNavM(:, 1:2) .^ 2, 2));
SPEED_SPREAD_MIN = 0.15;                                                  % the turns must span speeds for the two to part
speedSpreadH = NaN;
if hasAtt && nnz(turnH) >= 20, speedSpreadH = std(speedK(turnH)) / max(mean(speedK(turnH)), eps); end
if hasAtt && nnz(turnH) >= 20 && speedSpreadH >= SPEED_SPREAD_MIN
    [armJoint, armJointSig, latFit, latFitSig] = armLatencyFit(navMinusRows, yawRate, speedK, uAlong, uAcross, turnH);
    medSpeed = median(speedK(turnH));
    slow = turnH & speedK <= medSpeed;  fast = turnH & speedK > medSpeed;
    if nnz(slow) >= 10 && nnz(fast) >= 10
        bandSlope(1, :) = slopeOnYaw(navMinusRows, yawRate, uAlong, uAcross, slow);
        bandSlope(2, :) = slopeOnYaw(navMinusRows, yawRate, uAlong, uAcross, fast);
        bandSpeed = [mean(speedK(slow)), mean(speedK(fast))];
    end
end
if hasTruth && hasAtt
    turn = abs(yawRate) > 2 * pi / 180 & all(isfinite(rowsErr(:, 1:2)), 2) & inWin;
    nTurn = nnz(turn);
    if nTurn >= 10
        [armFit, armFitSig] = slopeOnYaw(-rowsErr(:, 1:2), yawRate, uAlong, uAcross, turn);
    end
end
% the same against the truth, truth minus rows, which the filter does not dilute
armTrueJoint = NaN(1, 2);  armTrueJointSig = NaN(1, 2);  latTrue = NaN(1, 2);  latTrueSig = NaN(1, 2);  bandSlopeTrue = NaN(2, 2);  bandSpeedTrue = NaN(1, 2);
if hasTruth && hasAtt && nTurn >= 20 && std(speedK(turn)) / max(mean(speedK(turn)), eps) >= SPEED_SPREAD_MIN
    [armTrueJoint, armTrueJointSig, latTrue, latTrueSig] = armLatencyFit(-rowsErr(:, 1:2), yawRate, speedK, uAlong, uAcross, turn);
    medSpeed = median(speedK(turn));
    slow = turn & speedK <= medSpeed;  fast = turn & speedK > medSpeed;
    if nnz(slow) >= 10 && nnz(fast) >= 10
        bandSlopeTrue(1, :) = slopeOnYaw(-rowsErr(:, 1:2), yawRate, uAlong, uAcross, slow);
        bandSlopeTrue(2, :) = slopeOnYaw(-rowsErr(:, 1:2), yawRate, uAlong, uAcross, fast);
        bandSpeedTrue = [mean(speedK(slow)), mean(speedK(fast))];
    end
end
episodes = struct('tStart', {}, 'tEnd', {}, 'detected', {}, 'delay', {}, 'coverage', {});
falseRuns = struct('tStart', {}, 'tEnd', {}, 'epochs', {}, 'blindBefore', {}, 'coastErrStart', {});
BLIND_WIN = 30;                                                           % s before a false run searched for epochs without rows
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
    tailWin = o.release + floor(o.hold / o.releaseGrow);                 % the longest release the rule can ask for
    for k = find(cleanFlag)'
        a = max(k - tailWin, 1);
        tailFlag(k) = any(flag(a:k - 1) & errH(a:k - 1) >= o.cleanLevel);
    end
    falseFlag = cleanFlag & ~tailFlag;
    nFalse = nnz(falseFlag);
    nTail = nnz(tailFlag);
    runsF = runsOf(falseFlag);
    for i = 1:size(runsF, 1)
        a = runsF(i, 1);
        before = t >= t(a) - BLIND_WIN & t < t(a) & inWin;
        falseRuns(end + 1) = struct('tStart', t(a), 'tEnd', t(runsF(i, 2)), 'epochs', runsF(i, 2) - a + 1, ...
            'blindBefore', nnz(before & nRate == 0), 'coastErrStart', coastErr(a)); %#ok<AGROW>
    end
end

res = struct('t', t, 'vRows', vRows, 'dvIns', dvIns, 'onset', onset, 'thresh', thresh, 'flag', flag, 'sep', sep, ...
    'rowsErr', rowsErr, 'episodes', episodes, 'falseRuns', falseRuns, 'nFlag', nnz(flag), 'nFalse', nFalse, 'nTail', nTail, ...
    'fbGain', fbGain, 'fbLag', fbLag, 'fbModel', fbModel, 'totShare', totShare, 'totShareT', totShareT, 'totRms', totRms, ...
    'incErr', incErr, 'lossRatio', lossRatio, 'lossClass', lossClass, 'lossNotTaken', lossNotTaken, 'hoRms', hoRms, ...
    'vCoast', vCoastK, 'coastErr', coastErr, 'coastDrift', coastDrift, 'coastLost', coastLost, ...
    'fbScan', fbScan, 'yawRate', yawRate, 'armFit', armFit, 'armFitSig', armFitSig, ...
    'armHost', armHost, 'armHostSig', armHostSig, 'armJoint', armJoint, 'armJointSig', armJointSig, 'latFit', latFit, ...
    'latFitSig', latFitSig, 'bandSlope', bandSlope, 'bandSpeed', bandSpeed, 'armTrueJoint', armTrueJoint, ...
    'armTrueJointSig', armTrueJointSig, 'latTrue', latTrue, 'latTrueSig', latTrueSig, 'bandSlopeTrue', bandSlopeTrue, ...
    'bandSpeedTrue', bandSpeedTrue, ...
    'fbScanWide', fbScanWide, 'alignLag', alignLag, 'fbScanSwap', fbScanSwap, 'frameOff', frameOff, 'firC', firC, 'firSe', firSe, ...
    'firTotal', firTotal, 'firTotalSe', firTotalSe, 'lgTotal', lgTotal, 'lgC', lgC, 'tauPair', tauPair, ...
    'pairShare', pairShare, 'pairJump', pairJump, 'jumpClusters', jumpClusters, 'opts', o);
if o.quiet, return; end

% print
fprintf('dopplerEpisodeCheck: %d epochs, %d with rate rows; floor %.2f m/s, drift %.3f m/s2, heading sigma %.1f mrad, anchor filter %.2f, confirm %d, release %d + 1 per %d held, hold %d\n', ...
    nnz(inWin), nnz(nRate > 0 & inWin), o.floor, o.driftRate, o.headingSigma * 1e3, o.alpha, o.confirm, o.release, o.releaseGrow, o.hold);
legacy = o.verbose || ~(useTotal || useHand);                         % the feedback reads matter under a fallback model, or on request
fprintf('   inputs:\n');
if isempty(o.fbGain) && ~found
    fprintf('   USED: no feedback found under either model; the INS increment is taken as logged\n');
elseif useTotal
    fprintf('   USED: the nav velocity less the prior state, the previous update''s step added back; nothing fitted, no feedback model\n');
elseif useHand
    fprintf('   USED: the nav velocity less the prior state for the rows; the nav''s change plus the logged correction for the increment%s\n', ...
        pick(~fbPerSample && ~isempty(o.fbLagWeights), sprintf(', delivered as given: weights by lag %s', sprintf('%.2f ', o.fbLagWeights)), ...
        pick(fbPerSample, ', per sample', sprintf(', per epoch, delivered with total %+.2f over %d epochs', lgTotal, FIR_LAGS))));
elseif useLog && ~fbPerSample && ~isempty(o.fbLagWeights)
    fprintf('   USED: the logged velocity correction, per epoch, delivered as given: weights by lag %s\n', sprintf('%.2f ', o.fbLagWeights));
elseif useLog
    fprintf('   USED: the logged velocity correction%s\n', pick(fbPerSample, ', per sample', sprintf(', per epoch, delivered with total %+.2f over %d epochs', lgTotal, FIR_LAGS)));
elseif useFir
    fprintf('   USED: the increment''s delivery profile, total %+.2f over %d epochs (fitted; the mechanisation slews the Kalman step in)\n', firTotal, FIR_LAGS);
else
    fprintf('   USED: %+.3f of the state at lag %d (%s)\n', fbGain, fbLag, pick(isempty(o.fbGain), 'fitted', 'given'));
end
if hasInc && hasTruth && isfinite(totRms)
    fprintf('   the total-state increment, the change of nav minus prior with the previous update''s step added back, against the truth: rms %.4f m/s per epoch over the quiet epochs;\n', totRms);
    fprintf('      share of the state left in it %+.3f (t %.1f; 0 when the mechanisation takes the whole hand-over and the prior is the state the rows were predicted with)\n', ...
        totShare, totShareT);
    if abs(totShareT) >= T_MIN && abs(totShare) > 0.01
        fprintf('      (a share left: part of the hand-over does not reach the velocity, or the prior is not the state the rows were predicted with)\n');
    end
end
if hasInc && ~isempty(lossClass)
    nL = [lossClass.n];
    txtL = '';
    for i = find(nL > 0)
        txtL = [txtL, sprintf(' %s on %d,', lossClass(i).name, lossClass(i).n)]; %#ok<AGROW>
    end
    fprintf('   the state''s loss between an update and the next prior, as a share of that update''s posterior (horizontal, posterior above 0.03 m/s):%s\n', txtL(1:end - 1));
    if hasTruth
        txtS = '';
        for i = find(isfinite([lossClass.slope]))
            txtS = [txtS, sprintf(' %s: slope %+.2f on %d epochs (the nav took %.2f of it), rms %.3f;', lossClass(i).name, lossClass(i).slope, ...
                lossClass(i).nUsed, 1 - lossClass(i).slope, lossClass(i).rms)]; %#ok<AGROW>
        end
        if ~isempty(txtS), fprintf('      increment error on the loss per class, epochs under 3 m/s2:%s\n', txtS(1:end - 1)); end
        big = incErr > 0.10 & inWin & isfinite(lossRatio);
        odd = big & ~(lossRatio >= 0.03 & lossRatio < 0.17);
        fprintf('      increment errors above 0.10 m/s: %d epochs, %d of them at a loss other than a tenth\n', nnz(big), nnz(odd));
        if lossNotTaken
            fprintf('      LOSS NOT TAKEN: on some epochs the state lost more than the nav took; nav minus prior steps there by the difference, which\n');
            fprintf('      the total-state increment reads as motion; the screen needs the hand-over itself on those epochs (the per-epoch table shows the loss)\n');
        end
    end
end
if isfinite(hoRms) && isfinite(totRms)
    fprintf('   the two coasts against the truth, rms per epoch over the quiet epochs: total-state %.4f m/s, hand-over %.4f m/s (the nav''s change plus the\n', totRms, hoRms);
    fprintf('      logged correction); they differ by the error growth the filter''s model expects between updates, and the smaller is the better coast\n');
end
if legacy || useLog || useHand
    if hasFb && hasTruth
        if fbPerSample
            fprintf('   the logged velocity correction, summed per interval: the nav''s change follows it with share %+.2f (t %.1f; 1 means the log and the sign are right)\n', ...
                lgTotal, lgTotal / max(lgTotalSe, eps));
        else
            fprintf('   the logged velocity correction per epoch, delivered to the nav over the following epochs: total %+.2f (t %.1f); by lag 0 to 7: %s%s\n', ...
                lgTotal, lgTotal / max(lgTotalSe, eps), sprintf('%+.3f ', lgC(1:8)), pick(lgLag0Dropped, '(lag 0 read negative, cause and effect reversed, dropped)', ''));
        end
        if isfinite(lgTotal) && abs(lgTotal - 1) > 0.2
            fprintf('      (a total away from 1: the logged correction is not what reaches the velocity, or its sign or frame differ)\n');
        end
    end
end
if hasAtt && all(isfinite(armHost))
    fprintf('   nav minus rows per rad/s of yaw rate, the host''s lever-arm error in this frame (leverArmCheck''s read): along %+.2f m, across %+.2f m (sigma %.2f / %.2f)\n', ...
        armHost(1), armHost(2), armHostSig(1), armHostSig(2));
end
if hasAtt && isfinite(speedSpreadH) && speedSpreadH < SPEED_SPREAD_MIN
    fprintf('   arm or latency: the turning epochs sit at one speed (%.0f m/s, spread %.0f%%), so an arm and a latency cannot be told apart here\n', ...
        mean(speedK(turnH)), 100 * speedSpreadH);
end
if hasAtt && all(isfinite(latFit))
    fprintf('   arm or latency, nav minus rows on the yaw rate and on speed times yaw rate jointly: across, an arm of %+.2f m (sigma %.2f) and a latency of %+.3f s (sigma %.3f);\n', ...
        armJoint(2), armJointSig(2), latFit(2), latFitSig(2));
    fprintf('      along, %+.2f m (sigma %.2f) and %+.3f s (sigma %.3f); the across latency is the read, the along one a check\n', ...
        armJoint(1), armJointSig(1), latFit(1), latFitSig(1));
    if all(isfinite(bandSpeed))
        fprintf('      the across slope per rad/s at %.0f m/s %+.2f m and at %.0f m/s %+.2f m (an arm reads the same, a latency in proportion to the speed)\n', ...
            bandSpeed(1), bandSlope(1, 2), bandSpeed(2), bandSlope(2, 2));
    end
    if all(isfinite(latTrue))
        fprintf('      against the truth, undiluted by the filter: across, an arm of %+.2f m (sigma %.2f) and a latency of %+.3f s (sigma %.3f); along %+.2f m and %+.3f s\n', ...
            armTrueJoint(2), armTrueJointSig(2), latTrue(2), latTrueSig(2), armTrueJoint(1), latTrue(1));
        if all(isfinite(bandSpeedTrue))
            fprintf('      its across slope per rad/s at %.0f m/s %+.2f m and at %.0f m/s %+.2f m\n', bandSpeedTrue(1), bandSlopeTrue(1, 2), bandSpeedTrue(2), bandSlopeTrue(2, 2));
        end
        vLat = latTrue(2);  vLatSig = latTrueSig(2);  vArm = armTrueJoint(2);  vArmSig = armTrueJointSig(2);
    else
        vLat = latFit(2);  vLatSig = latFitSig(2);  vArm = armJoint(2);  vArmSig = armJointSig(2);
    end
    if abs(vLat) > 2 * vLatSig && abs(vLat) > 0.01 && abs(vArm) < 2 * vArmSig
        fprintf('      LATENCY: the across slope is a pairing latency of %.0f ms, not an arm; the rows are paired with a nav sample that much apart from their time\n', ...
            1e3 * vLat);
    elseif abs(vLat) > 2 * vLatSig && abs(vLat) > 0.01
        fprintf('      (both an arm and a latency read above their sigma across track; a latency of %.0f ms and an arm of %.2f m)\n', 1e3 * vLat, vArm);
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
        if isempty(o.armSlope) && isempty(o.truthSlope) && all(isfinite(armHost))
            fprintf('      pass opts.armSlope = [%+.2f %+.2f] (the host''s arm error, out of the rows) and opts.truthSlope = [%+.2f %+.2f] (what is left, the\n', ...
                armHost(1), armHost(2), armFit(1) - armHost(1), armFit(2) - armHost(2));
            fprintf('      truth point''s offset from the nav''s, out of the scoring)\n');
        elseif any(abs(armFit) > 2 * armFitSig & abs(armFit) > 0.1)
            if isempty(o.armSlope)
                fprintf('      (the rows carry a lever-arm error against the truth point; opts.armSlope = the host''s arm error above emulates a corrected arm)\n');
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
if legacy
    scanTxt = '';
    if scanOk
        scanTxt = sprintf('lag %d %+.3f t %.1f, lag %d %+.3f t %.1f, lag %d %+.3f t %.1f', fbScan(:, 1:3)');
    end
    if scanOk
        fprintf('   feedback share of the state per epoch in the nav velocity: %+.3f at lag %d (t %.1f; %s)\n', fbScan(iBest, 2), fbScan(iBest, 1), tBest, scanTxt);
    else
        fprintf('   feedback share of the state per epoch in the nav velocity: no scan without truth\n');
    end
    if hasTruth && isfinite(alignLag)
        iW = find(fbScanWide(:, 1) == alignLag, 1);
        fprintf('   the same share over lags %d to %d epochs: sharpest %+.3f at lag %d (t %.1f)%s\n', WIDE(1), WIDE(end), fbScanWide(iW, 2), alignLag, ...
            fbScanWide(iW, 3), pick(o.tShift ~= 0, sprintf('; run shifted by %+.2f s', o.tShift), ''));
        if alignOff && ~useTotal && ~useHand
            fprintf('      NOT ALIGNED: the feedback shows up %d epochs from where it belongs; the records and the run are not on the same clock,\n', alignLag - 1);
            fprintf('      give opts.tShift = %+.2f s (added to the run''s time)\n', o.tShift - (alignLag - 1) * median(diff(t)));
        elseif alignOff
            fprintf('      (the share sits %d epochs from its place; this read assumes the nav follows the state, which a host holding rows out breaks, so it is a note here)\n', alignLag - 1);
        end
    end
    if hasTruth && any(isfinite(fbScanSwap(:, 3)))
        [tSw, iSw] = max(abs(fbScanSwap(:, 3)));
        fprintf('   the share with north and east swapped: %+.3f at lag %d (t %.1f); records read as %s\n', fbScanSwap(iSw, 2), fbScanSwap(iSw, 1), tSw, o.frame);
        if frameOff && ~useTotal && ~useHand
            fprintf('      FRAME: the feedback appears only with the axes swapped: the records'' state vector is in the other frame, pass opts.frame = ''%s''\n', ...
                pick(isENU, 'NED', 'ENU'));
        elseif frameOff
            fprintf('      (the share appears only swapped; the same note: with rows held out this read no longer tells a frame)\n');
        end
    end
    if hasTruth && isfinite(firTotal)
        fprintf('   the update''s increment x+ minus x- delivered to the nav over the following epochs: total %+.2f (t %.1f) over %d epochs;\n', ...
            firTotal, firTotal / max(firTotalSe, eps), FIR_LAGS);
        fprintf('      by lag 0 to 7: %s; implied feedback %.4f m/s per epoch against %.4f from the state''s share%s\n', ...
            sprintf('%+.3f ', firC(1:8)), impliedFir, impliedState, pick(firLag0Dropped, ' (lag 0 read negative, cause and effect reversed, dropped)', ''));
    end
end
fprintf('   coast:\n');
if hasTruth && ~isempty(coastDrift)
    dr = coastDrift(:, 3);  drS = sort(dr);
    fprintf('   the coast''s drift against the truth over the %d flagged stretches of 10 epochs or more: median %+.4f m/s per s, 90th percentile %+.4f, largest %+.4f;\n', ...
        numel(dr), median(dr), drS(max(1, ceil(0.9 * numel(dr)))), max(dr));
    fprintf('      the tolerance grows at %.4f per s; a coast that drifts faster is flagged on clean rows until the tolerance catches it\n', o.driftRate);
    if median(dr) > o.driftRate
        fprintf('      DRIFT: the coast drifts faster than the tolerance grows on most long stretches: the IMU''s coast is worse than DRIFT_RATE says\n');
    end
end
if ~isempty(coastLost)
    fprintf('   the coast lost after %d epochs without rows: %d times, at%s s; the next rows anchored it afresh\n', o.hold, numel(coastLost), ...
        sprintf(' %.0f', coastLost(1:min(end, 8))));
end
fprintf('   screen:\n');
if hostScored
    fprintf('   SCORING THE HOST''S OWN FLAGS (opts.hostLog): the separation and tolerance below are the host screen''s, the coast columns the emulation''s\n');
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
        o.release + floor(o.hold / o.releaseGrow), o.cleanLevel, nTail);
    fprintf('   false flags (flagged while the rows were within %.2f m/s of the truth, not the tail of an episode): %d epochs in %d runs\n', ...
        o.cleanLevel, nFalse, numel(falseRuns));
    if ~isempty(falseRuns)
        postBlind = [falseRuns.blindBefore] >= 10;
        fprintf('      after a blind stretch of 10 epochs or more within %d s: %d runs, %d epochs; coast error at the runs'' start: median %.3f m/s\n', ...
            BLIND_WIN, nnz(postBlind), sum([falseRuns(postBlind).epochs]), median([falseRuns.coastErrStart]));
        fprintf('      run                       epochs   blind epochs before   coast error at start\n');
    end
    for i = 1:min(numel(falseRuns), 12)
        fprintf('      %8.2f to %8.2f s  %3d        %3d                 %.3f\n', falseRuns(i).tStart, falseRuns(i).tEnd, falseRuns(i).epochs, ...
            falseRuns(i).blindBefore, falseRuns(i).coastErrStart);
    end
    if numel(falseRuns) > 12, fprintf('      ... %d more runs\n', numel(falseRuns) - 12); end
    % the per-epoch table around the largest episode
    tt = o.tTable;
    if isempty(tt) && ~isempty(episodes)
        [~, iBig] = max([episodes.tEnd] - [episodes.tStart]);
        tt = [episodes(iBig).tStart - 10, episodes(iBig).tEnd + 10];
    end
    if ~isempty(tt)
        fprintf('   per epoch, %.0f to %.0f s: rows - truth, the rows'' separation from the INS coast against the tolerance, the flag,\n', tt(1), tt(2));
        fprintf('   on the epochs faster than 1 m/s2 the rows'' error as a time offset, the increment''s error against the truth, the\n');
        fprintf('   state''s loss since the previous update as a share of its posterior, and the coast''s own error against the truth\n');
        fprintf('        t      rows-truth N  E   |.|      sep  tolerance   flag   offset  incErr   loss   coastErr\n');
        for k = find(t >= tt(1) & t <= tt(2))'
            tauTxt = '      ';
            if isfinite(tauPair(k)), tauTxt = sprintf('%+6.2f', tauPair(k)); end
            incTxt = '     ';
            if isfinite(incErr(k)), incTxt = sprintf('%5.3f', incErr(k)); end
            lossTxt = '     ';
            if isfinite(lossRatio(k)), lossTxt = sprintf('%+5.2f', lossRatio(k)); end
            coastTxt = '     ';
            if isfinite(coastErr(k)), coastTxt = sprintf('%5.3f', coastErr(k)); end
            fprintf('   %8.2f    %+6.3f %+6.3f  %5.3f    %5.3f   %5.3f     %s   %s  %s  %s   %s\n', t(k), rowsErr(k, 1), rowsErr(k, 2), errH(k), ...
                sep(k), thresh(k), pick(flag(k), 'FLAG', '    '), tauTxt, incTxt, lossTxt, coastTxt);
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
        fprintf('   VERDICT: %d of %d episodes missed: the floor or the tolerance is too generous for this receiver\n', ...
            nnz(~[episodes.detected]), numel(episodes));
    end
else
    fprintf('   (no truth in the run: flags only, no scoring)\n');
end
end

%------------------------------------------------------------------------
function [c, se, total, totalSe, ok, applied, lag0Dropped] = firFit(series, dNav, dTrue, inWin, K, nLags, tMin, doFit)
% the delivery profile of a per-epoch series into the nav's velocity: (nav change - truth change)(k) regressed
% jointly on series(k - j), j = 0..nLags, both axes stacked, over the quiet epochs, with one trim; applied is the
% series delivered with that profile, [K x 3]. A negative coefficient at lag 0 is the reverse of cause and
% effect, the correction of an epoch answering the error that grew in the interval before it, not a delivery:
% the fit is then repeated without lag 0
c = NaN(nLags + 1, 1);  se = NaN(nLags + 1, 1);  total = NaN;  totalSe = NaN;  ok = false;  applied = zeros(K, 3);
lag0Dropped = false;
if ~doFit, return; end
kk = (nLags + 2:K)';
bb = dNav(kk, 1:2) - dTrue(kk, 1:2);
use = all(isfinite(bb), 2) & inWin(kk) & sqrt(sum(dTrue(kk, 1:2) .^ 2, 2)) < 0.5;
nUse = nnz(use);
if nUse < 10 * (nLags + 1), return; end
kUse = kk(use);
A = zeros(2 * nUse, nLags + 1);
for j = 0:nLags
    A(:, j + 1) = [series(kUse - j, 1); series(kUse - j, 2)];
end
b = [bb(use, 1); bb(use, 2)];
cc = -(A \ b);
r = b + A * cc;
madS = max(1.4826 * median(abs(r - median(r))), 0.2 * sqrt(mean(r .^ 2)));
keep = abs(r) < 4 * madS;
if nnz(keep) >= 10 * (nLags + 1)
    A = A(keep, :);  b = b(keep);
    cc = -(A \ b);  r = b + A * cc;
end
if cc(1) < 0
    lag0Dropped = true;
    A = A(:, 2:end);
    cc = [0; -(A \ b)];  r = b + A * cc(2:end);
    s2 = (r' * r) / max(numel(r) - nLags, 1);
    covC = zeros(nLags + 1);  covC(2:end, 2:end) = s2 * inv(A' * A);
else
    s2 = (r' * r) / max(numel(r) - nLags - 1, 1);
    covC = s2 * inv(A' * A);
end
c = cc;  se = sqrt(diag(covC));
total = sum(cc);  totalSe = sqrt(max(sum(covC(:)), 0));
ok = totalSe > 0 && abs(total) / totalSe >= tMin;
for j = 0:nLags
    applied(j + 1:K, :) = applied(j + 1:K, :) + cc(j + 1) * series(1:K - j, :);
end
end

function [arm, armSig, lat, latSig] = armLatencyFit(vec, yawRate, speed, uAlong, uAcross, use)
% a [K x 2] horizontal vector projected along and across the heading and regressed on the yaw rate and on speed
% times yaw rate over the epochs 'use', with an intercept: the arm part (m per rad/s) and the latency part (s), each
% with its sigma, along and across
arm = NaN(1, 2);  armSig = NaN(1, 2);  lat = NaN(1, 2);  latSig = NaN(1, 2);
n = nnz(use);
eA = sum(vec(use, :) .* uAlong(use, :), 2);  eX = sum(vec(use, :) .* uAcross(use, :), 2);
A = [yawRate(use), speed(use) .* yawRate(use), ones(n, 1)];
if rank(A) < 3, return; end
for p = 1:2
    if p == 1, yv = eA; else, yv = eX; end
    c = A \ yv;  r = yv - A * c;
    covC = (r' * r) / max(n - 3, 1) * inv(A' * A);
    arm(p) = c(1);  armSig(p) = sqrt(covC(1, 1));  lat(p) = c(2);  latSig(p) = sqrt(covC(2, 2));
end
end

function [slope, sig] = slopeOnYaw(vec, yawRate, uAlong, uAcross, use)
% a [K x 2] horizontal vector projected along and across the heading and regressed on the yaw rate over the
% epochs 'use', with an intercept: [along across] slopes in m per rad/s and their sigma
slope = NaN(1, 2);  sig = NaN(1, 2);
n = nnz(use);
eA = sum(vec(use, :) .* uAlong(use, :), 2);  eX = sum(vec(use, :) .* uAcross(use, :), 2);
A = [yawRate(use), ones(n, 1)];
for p = 1:2
    if p == 1, yv = eA; else, yv = eX; end
    c = A \ yv;  r = yv - A * c;
    covC = (r' * r) / max(n - 2, 1) * inv(A' * A);
    slope(p) = c(1);  sig(p) = sqrt(covC(1, 1));
end
end

function [v] = toNED(v, isENU)
% a [n x 3] vector set in the records' frame, returned in NED: ENU exchanges the first two and flips the third
if isENU, v = [v(:, 2), v(:, 1), -v(:, 3)]; end
end

function [out] = shareAtLag(xPlus, dNav, dTrue, inWin, L, K)
% the share of the state x+(k - L) in (nav change - truth change)(k) over the quiet epochs, both axes stacked,
% through the origin, with one trim of the outliers: [share, t-statistic, rms of the regressor]; NaN with too
% few epochs
out = NaN(1, 3);
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
out = [c, c / max(se, eps), sqrt(mean(av .^ 2))];
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
