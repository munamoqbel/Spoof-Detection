function [row] = matrixRow(spfLog, tKf, run, scn, name)
%MATRIXROW  One row of the injection-matrix table from the gate log and the nav errors of one spoofed run.
%   row = matrixRow(spfLog, tKf, run, scn)        spfLog from SPF_logAppend; tKf the epoch times (s); run the struct
%                                                 runCompare takes (t, navPos, tTrue, truePos at least), or [] for the
%                                                 gate numbers alone; scn the scenario struct of spoofScenario
%   row = matrixRow(spfLog, tKf, run, scn, name)  a name for the row (default the scenario)
%   Prints the reading of the run and returns its numbers; matrixTable(rows) prints several rows as one table.
%   The nav position error is nav minus truth in NED at the nav samples, with its mean over the minute before tOn
%   taken out, so what is read is what the spoof did, not where the nav sat already.
%   What is read: the gate mode through the window; the first latch at or after tOn (its delay from tOn, the
%   monitor and axis that fired, its ratio); the largest SS and CPI ratio from tOn to the hand-back (or the end), and
%   the largest on the clean armed epochs (before tOn, and from a minute after the hand-back); the nav error at the
%   latch, its largest horizontal excursion from tOn to the hand-back, the error at tOff (how far the nav followed)
%   and its mean over the window after 10 s; the hand-back (its delay from tOff, the coast length, the vetoes, the
%   error and the sigma carried at the hand-back); false alarms (latches outside the window and the hand-back);
%   faults (inputFault, solveFault, clamped, budget exceeded). After a latch, the re-validation through the coast: per
%   epoch the rows the host handed over, the chi-square q against its threshold, pass or not and the dwell, with the
%   first pass and whether the probation opened while the spoof was still on (a hand-back into the spoof: the
%   re-validation passed on the rows and S the host handed over, fewer rows than before the latch when the host
%   rejected the rows carrying the spoof, or q far below its threshold on a full row set when S carries an R the
%   host scaled up on large innovations).
%   Result: name, kind, size, tOn, tOff, tEnd, latchDelay, monitor, axis, ratioAtLatch, maxSsIn, maxCpiIn,
%   maxSsClean, maxCpiClean, errAtLatch [1 x 3], maxHorizIn, errAtOff [1 x 3], meanErrIn [1 x 3], handbackDelay,
%   coastLen, nVeto, errAtHandback [1 x 3], sigmaAtHandback [1 x 3], falseAlarms, tFalse, faults [1 x 4], timeline,
%   revalComputed, revalPassed, rowsCoast, rowsBefore (medians), qRatioCoast (median q over its threshold on the
%   computed coast epochs), tFirstPass, probationInSpoof.
%   MATLAB desktop tool, not codegen. Runs under Octave.

if nargin < 5 || isempty(name), name = sprintf('%s %g', scn.kind, scn.size); end
BASE = 60;                                                   % s: the baseline of the nav error is its mean over this before tOn
SETTLE = 10;                                                 % s: the mean error over the window starts this long after tOn
CLEAN_AFTER = 60;                                            % s: clean epochs start this long after the hand-back
tKf = tKf(:);
n = min([double(spfLog.count), double(spfLog.numEpochs), numel(tKf)]);
tKf = tKf(1:n);
tOn = scn.tOn;  tOff = scn.tOff;  tEnd = tKf(n);
mode = spfLog.mode(1:n);  armed = logical(spfLog.armed(1:n));
ssR = spfLog.ssRatio(1:n, :);  cpiR = spfLog.cpiRatio(1:n, :);
ev = @(f) double(spfLog.(f)(spfLog.(f) >= 1 & spfLog.(f) <= n));
latches = ev('evLatch');  handbacks = ev('evHandback');  vetoes = ev('evVeto');  probs = ev('evProbStart');
fprintf('matrixRow: %s, %s %g %s on %.2f to %.2f s (scenario), records to %.2f s', name, scn.kind, scn.size, unitOf(scn.kind), tOn, tOff, tEnd);
if any(armed), fprintf('; armed from %.2f s\n', tKf(find(armed, 1))); else, fprintf('; never armed\n'); end
%% the latch and the hand-back
kL = latches(find(tKf(latches) >= tOn, 1));
latchDelay = NaN;  monitor = 'none';  axisName = '-';  ratioAtLatch = NaN;  kH = [];  handbackDelay = NaN;  coastLen = NaN;  nVeto = 0;
if ~isempty(kL)
    latchDelay = tKf(kL) - tOn;
    look = max(1, kL - 1):kL;
    ss = any(spfLog.ssAlarm(look));  cpi = any(spfLog.cpiAlarm(look));
    if ss && cpi, monitor = 'SS+CPI';  elseif ss, monitor = 'SS';  elseif cpi, monitor = 'CPI'; end
    ax = any(spfLog.alarmAxis(look, :), 1);
    letters = 'NED';  axisName = letters(ax);
    if isempty(axisName), axisName = '-'; end
    if ss, ratioAtLatch = max(max(ssR(look, :))); else, ratioAtLatch = max(max(cpiR(look, :))); end
    kH = handbacks(find(handbacks > kL, 1));
    if ~isempty(kH)
        handbackDelay = tKf(kH) - tOff;  coastLen = tKf(kH) - tKf(kL);
        nVeto = nnz(vetoes > kL & vetoes < kH);
    else
        coastLen = tEnd - tKf(kL);  nVeto = nnz(vetoes > kL);
    end
end
if isempty(kH), kStop = n; else, kStop = kH; end
inSpan = (1:n)' >= find(tKf >= tOn, 1) & (1:n)' <= kStop;            % tOn to the hand-back (or the end)
maxSsIn = NaN;  maxCpiIn = NaN;  axSs = '-';  axCpi = '-';
if any(inSpan)
    [maxSsIn, i] = max(max(ssR(inSpan, :), [], 1));  axSs = 'NED';  axSs = axSs(i);
    [maxCpiIn, i] = max(max(cpiR(inSpan, :), [], 1));  axCpi = 'NED';  axCpi = axCpi(i);
end
clean = armed & tKf < tOn;
if ~isempty(kH), clean = clean | (armed & tKf > tKf(kH) + CLEAN_AFTER);
elseif isempty(kL), clean = clean | (armed & tKf > tOff + CLEAN_AFTER); end
maxSsClean = NaN;  maxCpiClean = NaN;
if any(clean), maxSsClean = max(max(ssR(clean, :)));  maxCpiClean = max(max(cpiR(clean, :))); end
falseL = latches(tKf(latches) < tOn | (~isempty(kH) & tKf(latches) > tKf(max([kH 1]))));
if isempty(kL), falseL = latches(tKf(latches) < tOn | tKf(latches) > tOff + CLEAN_AFTER); end
falseL = setdiff(falseL, kL);
faults = [nnz(spfLog.inputFault(1:n)), nnz(spfLog.solveFault(1:n)), nnz(spfLog.clamped(1:n)), nnz(spfLog.budgetOut(1:n))];
%% the mode timeline from 10 s before tOn to 10 s after the hand-back (or tOff)
if isempty(kH), tStop = tOff; else, tStop = tKf(kH); end
w = find(tKf >= tOn - 10 & tKf <= tStop + 10);
timeline = '';
letters = '?NCP';
for i = w'
    if i == w(1) || mode(i) ~= mode(i - 1)
        timeline = [timeline, sprintf(' %s %.1f', letters(min(max(round(mode(i)) + 1, 1), 4)), tKf(i))]; %#ok<AGROW>
    end
end
timeline = strtrim(timeline);
%% the gate lines
if isempty(kL)
    fprintf('   gate: no latch at or after tOn; mode through the window%s: %s\n', ifelse(all(mode(inSpan & tKf < tOff) == 1), ' NOMINAL throughout', ''), timeline);
else
    fprintf('   gate: latch at %.2f s, %.1f s after tOn, %s on %s at ratio %.2f', tKf(kL), latchDelay, monitor, axisName, ratioAtLatch);
    if isempty(kH)
        fprintf('; no hand-back by the end of the run (%.1f s after tOff), %d vetoes\n', tEnd - tOff, nVeto);
    else
        kP = probs(find(probs > kL & probs <= kH, 1));
        if isempty(kP), sP = ''; else, sP = sprintf('probation at %.2f s, ', tKf(kP)); end
        fprintf('; coast %.1f s; %s%d vetoes; hand-back at %.2f s, %.1f s after tOff\n', coastLen, sP, nVeto, tKf(kH), handbackDelay);
    end
    fprintf('   timeline (mode and the time it starts): %s\n', timeline);
end
fprintf('   largest ratio from tOn to the hand-back or the end: SS %.2f (%s), CPI %.2f (%s); on the clean armed epochs SS %.2f, CPI %.2f\n', ...
    maxSsIn, axSs, maxCpiIn, axCpi, maxSsClean, maxCpiClean);
%% the re-validation through the coast: what let the gate hand back
revalComputed = NaN;  revalPassed = NaN;  rowsCoast = NaN;  rowsBefore = NaN;  qRatioCoast = NaN;  tFirstPass = NaN;  probationInSpoof = false;
if ~isempty(kL) && isfield(spfLog, 'qReval')
    minMeas = 4;
    try, minMeas = double(CST_spfParam.REVAL_MIN_MEAS); catch, minMeas = 4; end %#ok<CTCH>
    segC = kL:kStop;
    q = spfLog.qReval(segC);  thr = spfLog.revalThr(segC);  rows = spfLog.numMeas(segC);  done = logical(spfLog.revalDone(segC));
    pass = done & q >= 0 & q < thr & rows >= minMeas;
    revalComputed = nnz(done);  revalPassed = nnz(pass);
    rowsCoast = median(rows);
    bef = max(1, kL - 120):kL - 1;
    if ~isempty(bef), rowsBefore = median(spfLog.numMeas(bef)); end
    okQ = done & thr > 0;
    if any(okQ), qRatioCoast = median(q(okQ) ./ thr(okQ)); end
    fprintf('   re-validation through the coast, %d epochs from the latch to the hand-back or the end: computed on %d, passed on %d; rows the host handed\n', ...
        numel(segC), revalComputed, revalPassed);
    fprintf('   over median %d (the minute before the latch %d); q at %.2f of its threshold (median of the computed epochs)\n', round(rowsCoast), round(rowsBefore), qRatioCoast);
    iP = find(pass, 1);
    if ~isempty(iP)
        tFirstPass = tKf(segC(iP));
        fprintf('   first pass at %.2f s, %.1f s after the latch', tFirstPass, tFirstPass - tKf(kL));
        if tFirstPass < tOff, fprintf(', while the spoof was on'); end
        fprintf('\n');
    end
    kP = probs(find(probs > kL & probs <= kStop, 1));
    probationInSpoof = ~isempty(kP) && tKf(kP) < tOff;
    pick = unique([segC(1:min(8, numel(segC))), segC(round(linspace(1, numel(segC), min(10, numel(segC)))))]);
    fprintf('   course (the first 8 epochs in a row, then spread to the end):\n');
    fprintf('   %10s %5s %5s %9s %8s %5s %5s\n', 't [s]', 'mode', 'rows', 'q', 'thr', 'pass', 'dwell');
    letters2 = '?NCP';
    for k = pick
        i = k - kL + 1;
        fprintf('   %10.2f %5s %5d %9.1f %8.1f %5s %5d\n', tKf(k), letters2(min(max(round(mode(k)) + 1, 1), 4)), rows(i), q(i), thr(i), ifelse(pass(i), 'yes', '-'), spfLog.dwell(k));
    end
    if probationInSpoof
        fprintf('   VERDICT: the probation opened at %.2f s while the spoof was on (tOff %.2f s): the re-validation passed on the rows and S the host\n', tKf(kP), tOff);
        fprintf('            handed over, a hand-back into the spoof');
        if rowsCoast < 0.8 * rowsBefore
            fprintf('; the host handed fewer rows during the coast (%d against %d): the rows carrying the spoof were\n            rejected before the gate saw them\n', round(rowsCoast), round(rowsBefore));
        elseif qRatioCoast < 0.5
            fprintf('; q sat at %.2f of its threshold on a full row set: S = H P_C H'' + R carries an R the host scaled up on the\n            large innovations\n', qRatioCoast);
        else
            fprintf('\n');
        end
    elseif ~isempty(kP)
        fprintf('   the probation opened at %.2f s, %.1f s after tOff: the re-validation held through the spoof\n', tKf(kP), tKf(kP) - tOff);
    elseif ~isempty(kH)
        fprintf('   NOTE: a hand-back at %.2f s without a probation start in the log\n', tKf(kH));
    else
        fprintf('   no probation and no hand-back by the end of the records\n');
    end
end
%% the nav error
errAtLatch = NaN(1, 3);  maxHorizIn = NaN;  errAtOff = NaN(1, 3);  meanErrIn = NaN(1, 3);  errAtHandback = NaN(1, 3);  sigmaAtHandback = NaN(1, 3);
if ~isempty(kH), sigmaAtHandback = spfLog.sigmaPos(kH, :); end
if ~isempty(run)
    t = run.t(:);
    err = navErrNed(run, t);
    base = t >= tOn - BASE & t < tOn & all(isfinite(err), 2);
    if any(base), err = err - repmat(mean(err(base, :), 1), numel(t), 1); end
    if isempty(kH), tH = tEnd; else, tH = tKf(kH); end
    span = t >= tOn & t <= tH & all(isfinite(err), 2);
    if any(span), maxHorizIn = max(hypot(err(span, 1), err(span, 2))); end
    if ~isempty(kL), errAtLatch = meanAround(t, err, tKf(kL), 0.5); end
    errAtOff = meanAround(t, err, tOff, 5);
    win = t >= tOn + SETTLE & t < tOff & all(isfinite(err), 2);
    if any(win), meanErrIn = mean(err(win, :), 1); end
    if ~isempty(kH), errAtHandback = meanAround(t, err, tKf(kH), 1); end
    fprintf('   nav error (nav minus truth, NED, its mean over %.0f to %.0f s taken out): at tOff N %+.1f E %+.1f D %+.1f m (mean of the 5 s before);\n', ...
        tOn - BASE, tOn, errAtOff);
    fprintf('   mean over %.0f to %.0f s N %+.1f E %+.1f D %+.1f m; largest horizontal from tOn to the hand-back or the end %.1f m\n', ...
        tOn + SETTLE, tOff, meanErrIn, maxHorizIn);
    if ~isempty(kL), fprintf('   at the latch N %+.1f E %+.1f D %+.1f m', errAtLatch); end
    if ~isempty(kH), fprintf('; at the hand-back N %+.1f E %+.1f D %+.1f m, sigma carried N %.1f E %.1f D %.1f m', errAtHandback, sigmaAtHandback); end
    if ~isempty(kL), fprintf('\n'); end
else
    fprintf('   no run given: no nav error\n');
end
%% false alarms and faults
if isempty(falseL), fprintf('   false alarms: none');
else, fprintf('   false alarms: %d latch(es) outside the window and the hand-back, at%s', numel(falseL), sprintf(' %.2f s', tKf(falseL))); end
if any(faults), fprintf('; faults: inputFault %d, solveFault %d, clamped %d, budget exceeded %d epochs\n', faults);
else, fprintf('; faults: none\n'); end
row = struct('name', name, 'kind', scn.kind, 'size', scn.size, 'tOn', tOn, 'tOff', tOff, 'tEnd', tEnd, ...
    'latchDelay', latchDelay, 'monitor', monitor, 'axis', axisName, 'ratioAtLatch', ratioAtLatch, ...
    'maxSsIn', maxSsIn, 'maxCpiIn', maxCpiIn, 'maxSsClean', maxSsClean, 'maxCpiClean', maxCpiClean, ...
    'errAtLatch', errAtLatch, 'maxHorizIn', maxHorizIn, 'errAtOff', errAtOff, 'meanErrIn', meanErrIn, ...
    'handbackDelay', handbackDelay, 'coastLen', coastLen, 'nVeto', nVeto, 'errAtHandback', errAtHandback, ...
    'sigmaAtHandback', sigmaAtHandback, 'falseAlarms', numel(falseL), 'tFalse', tKf(falseL)', 'faults', faults, 'timeline', timeline, ...
    'revalComputed', revalComputed, 'revalPassed', revalPassed, 'rowsCoast', rowsCoast, 'rowsBefore', rowsBefore, ...
    'qRatioCoast', qRatioCoast, 'tFirstPass', tFirstPass, 'probationInSpoof', probationInSpoof);
end

%------------------------------------------------------------------------
function [e] = navErrNed(run, t)
% nav minus truth in NED (m) at the nav samples t: truth interpolated on its own time base
a = 6378137;  e2 = 6.69437999014e-3;
tr = NaN(numel(t), 3);
for c = 1:3, tr(:, c) = interp1(run.tTrue(:), run.truePos(:, c), t, 'linear', NaN); end
L = tr(:, 1);  s2 = sin(L).^2;
Rm = a * (1 - e2) ./ (1 - e2 * s2).^1.5;  Rn = a ./ sqrt(1 - e2 * s2);
e = [(run.navPos(:, 1) - tr(:, 1)) .* (Rm + tr(:, 3)), (run.navPos(:, 2) - tr(:, 2)) .* (Rn + tr(:, 3)) .* cos(L), -(run.navPos(:, 3) - tr(:, 3))];
end

function [m] = meanAround(t, err, t0, before)
% the mean of err over the 'before' seconds up to t0 (the nearest sample when none falls inside)
sel = t > t0 - before & t <= t0 & all(isfinite(err), 2);
if any(sel), m = mean(err(sel, :), 1);
else, [~, i] = min(abs(t - t0));  m = err(i, :); end
end

function [u] = unitOf(kind)
if any(strcmp(kind, {'ramp', 'clockRamp'})), u = 'm/s'; else, u = 'm'; end
end

function [s] = ifelse(c, a, b)
if c, s = a; else, s = b; end
end
