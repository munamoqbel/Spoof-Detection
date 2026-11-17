function [res] = cpiReplay(spfMeasLog, tKf, spfLog, tWin, nominalBias, rMode, minRows, run, floors)
%CPIREPLAY  The CPI monitor replayed on the saved spfMeas records: what its projection reads, epoch by epoch, against the gate log.
%   res = cpiReplay(spfMeasLog, tKf, spfLog, tWin)               B from CST_spfParam on the path (the host's values)
%   res = cpiReplay(spfMeasLog, tKf, spfLog, tWin, nominalBias)  B per NED position axis (m)
%   res = cpiReplay(..., nominalBias, 'nominal')                 the monitor replayed with the nominal R instead of the R
%                                                 the host logged: per row type (range, rate, other) the run median of
%                                                 the logged R diagonal is taken as nominal and S is rebuilt as
%                                                 S - R_logged + R_nominal (exact, S - R is H P H'); what the gate would
%                                                 have read had the host handed it S with the nominal R, while the
%                                                 host's own update kept its scaled R. Rows the host rejected are not
%                                                 in the records, so this replay is conservative for a spoofed run
%   res = cpiReplay(..., nominalBias, 'host')     records from setKfMeasBoth / kfMeasFromUpdateBoth: the monitor replayed
%                                                 on the S and R the host applied (innovationCovHost, measNoiseCovHost)
%                                                 instead of the gate's: what the gate would have read on the old interface
%   res = cpiReplay(..., rMode, minRows)          an epoch with fewer range rows than minRows does not count (xi = 0,
%                                                 the window sum from SPF_cpiMonitor sees it as an epoch without rows):
%                                                 the few satellites left in a partial outage carry reflected signals
%                                                 or the gap's drift, and a one-axis projection on three rows turns one
%                                                 bad row into tens of metres; this measures the rule offline
%   res = cpiReplay(..., rMode, minRows, run)     the detailed windows read against the truth: run holds the arrays
%                                                 runCompare reads (t, navPos [lat lon h], tTrue, truePos), and each
%                                                 detailed window prints the logged nav's error from the truth beside the
%                                                 reading, so the printout says whether the ranges or the nav carry it
%   res = cpiReplay(..., minRows, run, floors)   the replay with the gate's floors carried into S: floors is a struct with
%                                                 states and sigmas (the host's COAST_SIGMA_FLOOR lists, or any) and, or
%                                                 only, clockSigma, a sigma for the clock columns the range rows use; S
%                                                 becomes S + H dP H' with dP the floor variance above the record's P on
%                                                 those states, and the sums are read again (dHatFloor, sigFloor,
%                                                 xiFloor, qFloor): what the CPI would read were its S to carry the
%                                                 gate's covariance, with the step one window then latches on
%   spfLog may be [] (no comparison with the log); tWin [tOn tOff] in s on tKf, the minute before it is the baseline.
%   Per epoch and monitored axis, from the record's innovation gamma, S and H: dHat = f' S^-1 gamma / (f' S^-1 f), the
%   offset the CPI projection reads on the axis (m); sig = 1 / sqrt(f' S^-1 f), the sigma of that reading (m); sigR =
%   1 / sqrt(f' R^-1 f), the same with R alone (what the sigma would be if the host carried no position uncertainty in
%   S); xi = (|dHat| - B) / sig floored at zero, as SPF_cpiMonitor forms it; q/T_N of the ten epochs ending at the
%   epoch from SPF_cpiMonitor itself; the gate log's cpiRatio at the same epoch; and sigP = sqrt(P+) on the axis, the
%   position sigma the host carries. Prints the table over the window, the medians over the minute before and over
%   the window, whether the replay matches the log (then the gate was fed these records), and what sets the reading:
%   how large the offset read is, how many epochs it stays above B, and whether S is R or the host's P.
%   Row screen, a diagnostic read beside the replay, not in the gate: per epoch the range rows are fitted by one position
%   shift [N E D] and one clock per system (every column beyond the velocity that a range row uses) with the replay's R;
%   a row whose residual exceeds ROW_K sigma of the fit is dropped, the worst first, while the fit keeps one degree of
%   redundancy, and the projection is read again on the rows left (dHatScreen, xiScreen, qScreen, nDropped). A spoof
%   that shifts the position moves every row together and drops nothing; a reflected or re-acquired satellite disagrees
%   with the others and is dropped. Every ten-epoch sum that crosses 1 in the window (the largest one when none does) is
%   printed epoch by epoch: the rows (the range rows per system), the reading, the shift [N E D] the rows imply with
%   the clocks left free and the clock shift per system (the same fit: what the rows say about the position once a
%   common mode per system is allowed, where the CPI's projection holds the clocks at the host's P), the fit's worst
%   row (its index in the record, its elevation from H, its residual over the sigma of the fit), the rows dropped, the
%   reading on the rows left and, with run, the logged nav's error from the truth (nav minus truth, NED) on the axis.
%   Each detailed window ends with what carries the reading over the epochs that count: the rows dropped, a common
%   mode of the set (a clock jump, the inter-channel biases of a new set) or the position the rows agree on; and, with
%   run, whether the ranges or the nav carry it, or both. The clock shifts of the fit, epoch to epoch, also give the
%   receiver's common-mode steps over the window (rms and the largest, with its time): the size a clock floor must cover.
%   Result: t, dHat, sig, sigR, xi, qReplay, qLog [K x nAxes], sigP [K x nAxes], matches, window, axes, rMode,
%   rNominal, minRows, nRangeRows, peaks (a cell per axis of [t q] for every run of the ten-epoch sum above 0.25 in
%   the window, one peak per run: the near misses, not only the largest), dHatScreen, xiScreen, qScreen [K x nAxes],
%   nDropped, worstRow, worstNorm [K x 1], shiftFree [K x 3] (the shift of the free-clock fit), clockFree {K x 1} (the
%   clock shift per system), rowK, detail (a struct per detailed window: tEnd, axis, q, qScreen, nDropped, shiftFree
%   and carrier, the free-clock shift on the axis over the counting epochs and what carries the reading, 'rows',
%   'clock' or 'position'; navErr and excess, the medians over the counting epochs of nav minus truth on the axis and
%   of |dHat| - |nav - truth|, NaN without run, and verdict, 'ranges' when the reading exceeds the nav's error by three
%   sigma, 'nav' within 1.5 sigma of it, 'both' between), sigFree [K x 3] (the sigma of the free-clock shift), clockStep
%   (rms, max, tMax, col: the common-mode steps between consecutive epochs), and with floors dHatFloor, sigFloor, xiFloor,
%   qFloor [K x nAxes] and detail.qFloor.
%   Needs SPF_cpiMonitor, CST_spfParam and the host's matrixInv on the path. MATLAB desktop tool, not codegen. Runs under Octave.

if iscell(spfMeasLog), spfMeasLog = [spfMeasLog{:}]; end
K = numel(spfMeasLog);
if nargin < 2 || isempty(tKf), tKf = 0.5 * (0:K-1)'; end
tKf = tKf(:);
if nargin < 3, spfLog = []; end
if nargin < 4 || isempty(tWin), tWin = [tKf(1) tKf(end)]; end
if nargin < 5 || isempty(nominalBias), nominalBias = CST_spfParam.NOMINAL_BIAS; end
if nargin < 6 || isempty(rMode), rMode = 'record'; end
if nargin < 7 || isempty(minRows), minRows = 0; end
if nargin < 8, run = []; end
if nargin < 9, floors = []; end
useNom = strcmpi(rMode, 'nominal');
useHost = strcmpi(rMode, 'host');
if useHost && ~isfield(spfMeasLog, 'measNoiseCovHost')
    error('cpiReplay:host', 'the records carry no measNoiseCovHost: the mode ''host'' needs records from setKfMeasBoth / kfMeasFromUpdateBoth');
end
BEFORE = 60;                                                 % s: the baseline minute before the window
ROW_K = 3;                                                   % the row screen: a range row beyond this many sigma of the fit is dropped
DETAIL_MAX = 6;                                              % detailed windows printed at most
MATCH_TOL = 0.02;                                            % in q/T_N: replay and log closer than this agree
N = double(CST_spfParam.WINDOW_LENGTH);  T = CST_spfParam.CPI_THRESHOLD;
axes = double(CST_spfParam.MONITORED_AXES);  nA = numel(axes);  letters = 'NED';
kW = find(tKf >= min(tWin) & tKf <= max(tWin));
kB = find(tKf >= min(tWin) - BEFORE & tKf < min(tWin));
if isempty(kW), error('cpiReplay:window', 'the window %g to %g s holds no epoch of the records (%.2f to %.2f s)', tWin, tKf(1), tKf(end)); end
kAll = [kB; kW];
dHat = NaN(K, nA);  sig = NaN(K, nA);  sigR = NaN(K, nA);  xi = zeros(K, nA);  sigP = NaN(K, nA);  qReplay = NaN(K, nA);  qLog = NaN(K, nA);
dHatS = NaN(K, nA);  sigS = NaN(K, nA);  xiS = zeros(K, nA);  qScreen = NaN(K, nA);                 % the row screen
nDropped = zeros(K, 1);  worstRow = zeros(K, 1);  worstNorm = NaN(K, 1);  worstEl = NaN(K, 1);  dropList = cell(K, 1);
shiftFit = NaN(K, 3);  clkFit = cell(K, 1);  sysRows = cell(K, 1);                                    % the free-clock fit
haveRun = ~isempty(run);
haveFloors = ~isempty(floors);
if haveFloors
    if ~isstruct(floors), error('cpiReplay:floors', 'floors is a struct with states and sigmas, or clockSigma, or both'); end
    if ~isfield(floors, 'states'), floors.states = zeros(1, 0); end
    if ~isfield(floors, 'sigmas'), floors.sigmas = zeros(1, 0); end
    if ~isfield(floors, 'clockSigma'), floors.clockSigma = 0; end
    if numel(floors.states) ~= numel(floors.sigmas), error('cpiReplay:floors', 'floors.states and floors.sigmas differ in length'); end
end
dHatF = NaN(K, nA);  sigF = NaN(K, nA);  xiF = zeros(K, nA);  qFloor = NaN(K, nA);                 % the floors in S
sigFree = NaN(K, 3);  clkColsList = cell(K, 1);
if haveRun && ~(isstruct(run) && all(isfield(run, {'t', 'navPos', 'tTrue', 'truePos'})))
    error('cpiReplay:run', 'run needs the fields t, navPos [lat lon h], tTrue and truePos, as runCompare reads them');
end
maxMeas = size(spfMeasLog(1).innovation, 1);  nState = size(spfMeasLog(1).obsMatrix, 2);
rNom = NaN(1, 3);
if useNom
    acc = {[], [], []};                                      % the logged R diagonal by row type over the run
    for k = 1:K
        s = spfMeasLog(k);  m = double(s.numMeas);
        if m < 1, continue; end
        typ = rowTypes(s.obsMatrix(1:m, :));  r = diag(s.measNoiseCov(1:m, 1:m));
        for c = 1:3, acc{c} = [acc{c}; r(typ == c)]; end
    end
    for c = 1:3, rNom(c) = medF(acc{c}); end
end
nRangeRows = zeros(K, 1);
for k = 1:K
    m = double(spfMeasLog(k).numMeas);
    if m >= 1, nRangeRows(k) = nnz(rowTypes(spfMeasLog(k).obsMatrix(1:m, :)) == 1); end
end
counts = nRangeRows >= minRows;                              % epochs the CPI may count
for k = kAll'
    s = spfMeasLog(k);  m = double(s.numMeas);
    for a = 1:nA, sigP(k, a) = sqrt(max(s.postCov(axes(a), axes(a)), 0)); end
    if m < 1 || ~counts(k), continue; end
    H = s.obsMatrix(1:m, :);  y = s.innovation(1:m);
    [S, R] = covOf(s, m, useNom, rNom, useHost);
    if ~all(isfinite(S(:))) || ~all(isfinite(y)), continue; end
    Sy = S \ y;
    for a = 1:nA
        f = H(:, axes(a));
        if ~any(f), continue; end
        s2 = f' * (S \ f);
        if s2 > 0 && isfinite(s2)
            g = f' * Sy;
            dHat(k, a) = g / s2;  sig(k, a) = 1 / sqrt(s2);
            xi(k, a) = max(abs(g) / sqrt(s2) - nominalBias(axes(a)) * sqrt(s2), 0);
        end
        r2 = f' * (R \ f);
        if r2 > 0 && isfinite(r2), sigR(k, a) = 1 / sqrt(r2); end
    end
    % the row screen: the range rows against one shift and one clock per system, the rows beyond ROW_K sigma dropped
    [drop, wRow, wNorm, wEl, shiftK, clkK, sysK, clkColsK, sigFreeK] = rowScreen(H, y, R, ROW_K);
    dropList{k} = drop;  nDropped(k) = numel(drop);  worstRow(k) = wRow;  worstNorm(k) = wNorm;  worstEl(k) = wEl;
    shiftFit(k, :) = shiftK;  clkFit{k} = clkK;  sysRows{k} = sysK;  clkColsList{k} = clkColsK;  sigFree(k, :) = sigFreeK;
    if haveFloors
        SF = floorOf(S, H, s.postCov, floors, clkColsK);  SFy = SF \ y;
        for a = 1:nA
            f = H(:, axes(a));
            if ~any(f), continue; end
            s2 = f' * (SF \ f);
            if s2 > 0 && isfinite(s2)
                g = f' * SFy;
                dHatF(k, a) = g / s2;  sigF(k, a) = 1 / sqrt(s2);
                xiF(k, a) = max(abs(g) / sqrt(s2) - nominalBias(axes(a)) * sqrt(s2), 0);
            end
        end
    end
    if isempty(drop)
        dHatS(k, :) = dHat(k, :);  sigS(k, :) = sig(k, :);  xiS(k, :) = xi(k, :);
    else
        keep = setdiff(1:m, drop);
        Sk = S(keep, keep);  yk = y(keep);  Syk = Sk \ yk;
        for a = 1:nA
            f = H(keep, axes(a));
            if ~any(f), continue; end
            s2 = f' * (Sk \ f);
            if s2 > 0 && isfinite(s2)
                g = f' * Syk;
                dHatS(k, a) = g / s2;  sigS(k, a) = 1 / sqrt(s2);
                xiS(k, a) = max(abs(g) / sqrt(s2) - nominalBias(axes(a)) * sqrt(s2), 0);
            end
        end
    end
end
%% the window sums from SPF_cpiMonitor itself: the ten epochs ending at each epoch
for k = kAll(kAll >= N)'
    idx = k - N + 1:k;
    innB = zeros(maxMeas, N);  covB = zeros(maxMeas, maxMeas, N);  obsB = zeros(maxMeas, nState, N);  numB = zeros(N, 1);
    for i = 1:N
        s = spfMeasLog(idx(i));  mi = double(s.numMeas);
        innB(:, i) = s.innovation(:);  covB(:, :, i) = s.innovationCov;  obsB(:, :, i) = s.obsMatrix;  numB(i) = mi;
        if ~counts(idx(i)), numB(i) = 0; end                 % an epoch below the row floor counts as one without rows
        if (useNom || useHost) && mi >= 1, covB(1:mi, 1:mi, i) = covOf(s, mi, useNom, rNom, useHost); end
    end
    for a = 1:nA
        [~, q] = SPF_cpiMonitor(innB, covB, obsB, numB, axes(a), nominalBias(axes(a)));
        qReplay(k, a) = q / T;
    end
    % the same sum with the gate's floors carried into S
    if haveFloors
        covF = covB;
        for i = 1:N
            s = spfMeasLog(idx(i));  mi = numB(i);
            if mi < 1, continue; end
            covF(1:mi, 1:mi, i) = floorOf(covB(1:mi, 1:mi, i), s.obsMatrix(1:mi, :), s.postCov, floors, clkColsList{idx(i)});
        end
        for a = 1:nA
            [~, q] = SPF_cpiMonitor(innB, covF, obsB, numB, axes(a), nominalBias(axes(a)));
            qFloor(k, a) = q / T;
        end
    end
    % the same sum on the rows the screen leaves: each epoch's rows compacted without the dropped ones
    if any(nDropped(idx) > 0)
        for i = reshape(find(nDropped(idx) > 0), 1, [])
            mi = numB(i);
            if mi < 1, continue; end
            keep = setdiff(1:mi, dropList{idx(i)});  nk = numel(keep);
            yi = innB(:, i);  Si = covB(:, :, i);  Hi = obsB(:, :, i);
            innB(:, i) = 0;  covB(:, :, i) = 0;  obsB(:, :, i) = 0;
            innB(1:nk, i) = yi(keep);  covB(1:nk, 1:nk, i) = Si(keep, keep);  obsB(1:nk, :, i) = Hi(keep, :);  numB(i) = nk;
        end
        for a = 1:nA
            [~, q] = SPF_cpiMonitor(innB, covB, obsB, numB, axes(a), nominalBias(axes(a)));
            qScreen(k, a) = q / T;
        end
    else
        qScreen(k, :) = qReplay(k, :);
    end
end
haveLog = ~isempty(spfLog) && isstruct(spfLog) && isfield(spfLog, 'cpiRatio');
if haveLog
    nL = min([size(spfLog.cpiRatio, 1), double(spfLog.count), K]);
    kk = kAll(kAll <= nL);
    qLog(kk, 1:min(nA, size(spfLog.cpiRatio, 2))) = spfLog.cpiRatio(kk, 1:min(nA, size(spfLog.cpiRatio, 2)));
end
%% print
fprintf('cpiReplay: %d records; window %.2f to %.2f s (epochs %d to %d), the minute before as the baseline; N = %d epochs, T_N = %.2f, B = [%s] m\n', ...
    K, tKf(kW(1)), tKf(kW(end)), kW(1), kW(end), N, T, strtrim(sprintf('%g ', nominalBias(axes))));
if useHost
    fprintf('   R HOST: the monitor replayed on the S and R the host applied (innovationCovHost, measNoiseCovHost): what the gate would read on the old interface\n');
end
if useNom
    fprintf('   R NOMINAL: per row type the run median of the logged R, range rows %.2f m, rate rows %.3f m/s, other rows %.2f; S rebuilt as\n', ...
        sqrt(rNom(1)), sqrt(rNom(2)), sqrt(rNom(3)));
    fprintf('   S - R_logged + R_nominal: what the gate reads when the host hands it S with the nominal R (the log column keeps the host''s R)\n');
end
if minRows > 0
    fprintf('   ROW FLOOR: an epoch with fewer than %d range rows does not count; %d of the %d epochs of the window and the minute before are below it\n', ...
        minRows, nnz(~counts(kAll)), numel(kAll));
end
fprintf('   per epoch and axis: dHat = the offset the CPI projection reads (m), sig = its sigma 1/sqrt(f'' S^-1 f) (m), xi = (|dHat| - B)/sig floored\n');
fprintf('   at 0; q/T_N = the ten epochs ending here, from SPF_cpiMonitor (alarm at 1); log = the gate log''s cpiRatio there; sigP = the host''s sqrt(P+)\n');
pick = unique([kW(1:min(8, numel(kW)))', kW(round(linspace(1, numel(kW), min(10, numel(kW)))))']);
fprintf('   %8s', 't [s]');
for a = 1:nA, fprintf(' | %s %6s %5s %5s %6s %6s %5s', letters(a), 'dHat', 'sig', 'xi', 'q/T_N', 'log', 'sigP'); end
fprintf('\n');
for k = pick
    fprintf('   %8.2f', tKf(k));
    for a = 1:nA
        fprintf(' | %s %6s %5s %5s %6s %6s %5s', ' ', num(dHat(k, a), '%+6.1f'), num(sig(k, a), '%5.2f'), num(xi(k, a), '%5.2f'), ...
            num(qReplay(k, a), '%6.3f'), num(qLog(k, a), '%6.3f'), num(sigP(k, a), '%5.2f'));
    end
    fprintf('\n');
end
fprintf('   medians, the minute before (%d epochs) then the window (%d epochs):\n', numel(kB), numel(kW));
fprintf('   %-9s', '');
for a = 1:nA, fprintf(' | %s %6s %5s %5s %6s %5s %5s', letters(a), 'dHat', 'sig', 'sigR', 'xi', 'sigP', ''); end
fprintf('\n');
for side = 1:2
    if side == 1, idx = kB;  lab = 'before'; else, idx = kW;  lab = 'window'; end
    fprintf('   %-9s', lab);
    for a = 1:nA
        fprintf(' | %s %6s %5s %5s %6s %5s %5s', ' ', num(medF(dHat(idx, a)), '%+6.1f'), num(medF(sig(idx, a)), '%5.2f'), num(medF(sigR(idx, a)), '%5.2f'), ...
            num(medF(xi(idx, a)), '%6.2f'), num(medF(sigP(idx, a)), '%5.2f'), '');
    end
    fprintf('\n');
end
%% the verdicts
matches = NaN;
% the log carries the CPI on the NOMINAL armed epochs only (the monitors stand down in COAST and PROBATION)
% the column the log is compared with: a gate that carries the floors in S logs the floored sums
qCmp = qReplay;  cmpName = 'the replay';
if haveFloors, qCmp = qFloor;  cmpName = 'the replay with the floors in S'; end
cmpEp = kW(all(isfinite(qCmp(kW, :)) & isfinite(qLog(kW, :)), 2));
if haveLog && isfield(spfLog, 'mode'), cmpEp = cmpEp(round(spfLog.mode(cmpEp)) == 1); end
if haveLog && isfield(spfLog, 'armed'), cmpEp = cmpEp(logical(spfLog.armed(cmpEp))); end
if haveLog && useNom
    both = cmpEp;
    if ~isempty(both)
        fprintf('   with the nominal R the largest ten-epoch sum over the window is %.3f (N %.3f, E %.3f, D %.3f); the log, with the host''s R, %.3f\n', ...
            max(max(qReplay(both, :))), max(qReplay(both, :)), max(max(qLog(both, :))));
        if haveFloors
            fprintf('   with the floors in S it is %.3f (N %.3f, E %.3f, D %.3f); a gate that carries them logs this column\n', max(max(qFloor(both, :))), max(qFloor(both, :)));
            matches = matchVerdict(qCmp, qLog, both, kW, tKf, cmpName, MATCH_TOL);
        end
    end
elseif haveLog
    both = cmpEp;
    if ~isempty(both)
        matches = matchVerdict(qCmp, qLog, both, kW, tKf, cmpName, MATCH_TOL);
    else
        fprintf('   no NOMINAL armed epoch of the window has both a replayed and a logged q/T_N: no comparison\n');
    end
end
for a = 1:nA
    [xiMax, i] = max(xi(kW, a));                                              % the epoch that counted most
    kM = kW(i);
    above = nnz(abs(dHat(kW, a)) > nominalBias(axes(a)));
    [qMax, j] = max(qReplay(kW, a));
    fprintf('   %s: the reading counts most at %.2f s, %+.1f m with sigma %.2f m (%.2f m from R alone, the host''s position sigma %.2f m), xi %.2f;\n', ...
        letters(a), tKf(kM), dHat(kM, a), sig(kM, a), sigR(kM, a), sigP(kM, a), xiMax);
    fprintf('      the reading stays beyond B = %g m on %d of %d epochs; the ten-epoch sum peaks at q/T_N = %.3f', nominalBias(axes(a)), above, numel(kW), qMax);
    if isfinite(qMax) && ~isempty(j), fprintf(' (window ending %.2f s)', tKf(kW(j))); end
    fprintf('\n');
end
% every run of the ten-epoch sum above PEAK_MIN, one peak per run: the near misses, not only the largest
PEAK_MIN = 0.25;
peaks = cell(1, nA);
for a = 1:nA
    q = qReplay(kW, a);  q(~isfinite(q)) = 0;
    pk = zeros(0, 2);  i = 1;
    while i <= numel(q)
        if q(i) >= PEAK_MIN
            j = i;
            while j < numel(q) && q(j + 1) >= PEAK_MIN, j = j + 1; end
            [v, m] = max(q(i:j));  pk(end+1, :) = [tKf(kW(i + m - 1)), v]; %#ok<AGROW>
            i = j + 1;
        else
            i = i + 1;
        end
    end
    peaks{a} = pk;
    if isempty(pk)
        fprintf('   %s: no run of the ten-epoch sum above %.2f in the window\n', letters(a), PEAK_MIN);
    else
        fprintf('   %s: runs of the ten-epoch sum above %.2f (time of the peak, value):%s\n', letters(a), PEAK_MIN, sprintf(' %.2f s %.3f;', pk'));
    end
end
%% the row screen over the window, and the windows that cross 1 read epoch by epoch
nRangeTot = sum(nRangeRows(kW));  nDropTot = sum(nDropped(kW));
qS = qScreen(kW, :);  qS(~isfinite(qS)) = 0;  qR = qReplay(kW, :);  qR(~isfinite(qR)) = 0;
fprintf('   ROW SCREEN (a diagnostic read, not in the gate): the range rows of each epoch fitted by one shift [N E D] and one clock per system with the\n');
fprintf('   replay''s R; a row beyond %g sigma of the fit is dropped, the worst first, while one degree of redundancy is left. Over the window %d of %d range\n', ROW_K, nDropTot, nRangeTot);
fprintf('   rows are dropped, on %d of %d epochs; the largest ten-epoch sum on the rows left is %.3f (N %.3f, E %.3f, D %.3f) against %.3f without the screen\n', ...
    nnz(nDropped(kW) > 0), numel(kW), max(max(qS)), max(qS), max(max(qR)));
if nDropTot > 0
    fprintf('   (a spoof that shifts the position moves every row together and drops nothing; a row dropped is one that disagrees with the others)\n');
end
% the receiver's common-mode steps: the clock shift of the fit, per clock column, between consecutive epochs
allCols = zeros(1, 0);
for k = kW', allCols = union(allCols, clkColsList{k}); end
clockStep = struct('rms', NaN, 'max', NaN, 'tMax', NaN, 'col', NaN);
steps = zeros(0, 3);                                         % [step t col]
for c = allCols
    prev = NaN;
    for k = kW'
        j = find(clkColsList{k} == c, 1);
        if isempty(j) || isempty(clkFit{k}) || ~isfinite(clkFit{k}(j)), continue; end
        v = clkFit{k}(j);
        if isfinite(prev), steps(end+1, :) = [v - prev, tKf(k), c]; end %#ok<AGROW>
        prev = v;
    end
end
if ~isempty(steps)
    [mx, i] = max(abs(steps(:, 1)));
    clockStep = struct('rms', sqrt(mean(steps(:, 1).^2)), 'max', mx, 'tMax', steps(i, 2), 'col', steps(i, 3));
    fprintf('   the clock common mode per system (the free-clock fit) steps between consecutive epochs by %.1f m rms over the window; the largest step\n', clockStep.rms);
    fprintf('   is %.1f m at %.2f s (state %d): a receiver clock or set-change step a clock floor in S has to cover, and what the host''s clock P must absorb\n', ...
        clockStep.max, clockStep.tMax, clockStep.col);
end
if haveFloors
    qF = qFloor(kW, :);  qF(~isfinite(qF)) = 0;
    nomF = kW;
    if haveLog && isfield(spfLog, 'mode'), nomF = kW(round(spfLog.mode(kW)) == 1); end
    if isempty(nomF), nomF = kW; end
    sig0 = medF3(sig(nomF, :));  sig1 = medF3(sigF(nomF, :));  stepFac = sqrt(T / N);
    fprintf('   FLOORS IN S: S + H dP H'' with the floors on states [%s] at sigma [%s] and the clock columns at %g m, where they exceed the record''s P;\n', ...
        strtrim(sprintf('%d ', floors.states)), strtrim(sprintf('%g ', floors.sigmas)), floors.clockSigma);
    fprintf('   the largest ten-epoch sum with the floors in S is %.3f (N %.3f, E %.3f, D %.3f) against %.3f without them\n', max(max(qF)), max(qF), max(max(qR)));
    fprintf('   the sigma of the reading, medians over the NOMINAL epochs: N %.2f E %.2f D %.2f m with the floors against %.2f %.2f %.2f m without;\n', sig1, sig0);
    fprintf('   the offset one window latches on, B + %.2f sigma: N %.1f E %.1f D %.1f m with the floors against %.1f %.1f %.1f m without\n', stepFac, ...
        nominalBias(axes) + stepFac * sig1, nominalBias(axes) + stepFac * sig0);
end
sel = zeros(0, 2);                                           % [k axis] of the windows to detail
for a = 1:nA
    pk = peaks{a};
    for i = 1:size(pk, 1)
        if pk(i, 2) >= 1, sel(end+1, :) = [kW(find(tKf(kW) >= pk(i, 1) - 1e-9, 1)), a]; end %#ok<AGROW>
    end
end
if isempty(sel)
    [~, iM] = max(qR(:));  [iK, aM] = ind2sub(size(qR), iM);
    if qR(iM) > 0, sel = [kW(iK), aM]; end
end
if isempty(sel)
    fprintf('   no ten-epoch sum above zero in the window: nothing to detail\n');
elseif size(sel, 1) > DETAIL_MAX
    fprintf('   %d windows cross 1; the first %d are detailed\n', size(sel, 1), DETAIL_MAX);
    sel = sel(1:DETAIL_MAX, :);
end
detail = struct('tEnd', {}, 'axis', {}, 'q', {}, 'qScreen', {}, 'qFloor', {}, 'nDropped', {}, 'shiftFree', {}, 'carrier', {}, 'navErr', {}, 'excess', {}, 'verdict', {});
for i = 1:size(sel, 1)
    kE = sel(i, 1);  a = sel(i, 2);  idx = max(kE - N + 1, 1):kE;
    fl = '';
    if haveFloors, fl = sprintf(', with the floors in S %.3f', qFloor(kE, a)); end
    fprintf('   window ending %.2f s on %s: q/T_N %.3f, on the rows left %.3f%s; B = %g m\n', tKf(kE), letters(a), qReplay(kE, a), qScreen(kE, a), fl, nominalBias(axes(a)));
    fprintf('      %8s %4s %-9s | %6s %5s %5s | %s %6s %6s %6s  %-17s | %s %4s %5s %7s %7s | %6s %5s | %s\n', 't [s]', 'rows', 'rng/system', 'dHat', 'sig', 'xi', ...
        'clocks free:', 'N', 'E', 'D', 'clock per system', 'worst', 'row', 'elev', 'res/sig', 'dropped', 'dHat''', 'xi''', 'nav - truth N, E, D (m)');
    navErrW = NaN(numel(idx), 3);
    for j = 1:numel(idx)
        k = idx(j);  m = double(spfMeasLog(k).numMeas);
        if haveRun
            navErrW(j, :) = nedDelta(interpRows(run.t(:), run.navPos, tKf(k)), interpRows(run.tTrue(:), run.truePos, tKf(k)));
        end
        if m < 1
            fprintf('      %8.2f %4d %-9s | no rows\n', tKf(k), 0, '0');  continue;
        elseif ~counts(k)
            fprintf('      %8.2f %4d %-9s | below the row floor: not counted\n', tKf(k), m, sprintf('%d', nRangeRows(k)));  continue;
        end
        if worstRow(k) > 0, wr = sprintf('%4d %5.0f %7.1f', worstRow(k), worstEl(k), worstNorm(k)); else, wr = sprintf('%4s %5s %7s', '-', '-', 'no fit'); end
        if haveRun, te = sprintf('%+6.1f %+6.1f %+6.1f', navErrW(j, :)); else, te = '(give run)'; end
        rs = sprintf('%d', nRangeRows(k));
        if numel(sysRows{k}) > 1, rs = sprintf('%d=%s', nRangeRows(k), strjoin(arrayfun(@(v) sprintf('%d', v), sysRows{k}, 'UniformOutput', false), '+')); end
        if isempty(clkFit{k}), cs = '-'; else, cs = strjoin(arrayfun(@(v) sprintf('%+.1f', v), clkFit{k}, 'UniformOutput', false), ' '); end
        fprintf('      %8.2f %4d %-9s | %6s %5s %5s | %s %6s %6s %6s  %-17s | %s %s %7d | %6s %5s | %s\n', tKf(k), m, rs, num(dHat(k, a), '%+6.1f'), num(sig(k, a), '%5.2f'), ...
            num(xi(k, a), '%5.2f'), '            ', num(shiftFit(k, 1), '%+6.1f'), num(shiftFit(k, 2), '%+6.1f'), num(shiftFit(k, 3), '%+6.1f'), cs, ...
            '     ', wr, nDropped(k), num(dHatS(k, a), '%+6.1f'), num(xiS(k, a), '%5.2f'), te);
    end
    d = struct('tEnd', tKf(kE), 'axis', a, 'q', qReplay(kE, a), 'qScreen', qScreen(kE, a), 'qFloor', qFloor(kE, a), 'nDropped', sum(nDropped(idx)), ...
        'shiftFree', NaN, 'carrier', '', 'navErr', NaN, 'excess', NaN, 'verdict', '');
    % the epochs that count: those whose xi is above zero, else every epoch with a reading
    cnt = idx(xi(idx, a) > 0 & isfinite(dHat(idx, a)));
    if isempty(cnt), cnt = idx(isfinite(dHat(idx, a))); end
    if ~isempty(cnt)
        medD = median(abs(dHat(cnt, a)));  medS = medF(abs(dHatS(cnt, a)));  medFree = medF(shiftFit(cnt, a));  d.shiftFree = medFree;
        sgFree = medF(sigFree(cnt, a));  sg0 = medF(sig(cnt, a));
        if isfinite(medS) && medS < 0.5 * medD
            d.carrier = 'rows';
            fprintf('      over the %d epochs that count the reading is %.1f m and %.1f m on the rows left: the rows dropped carry it\n', numel(cnt), medD, medS);
        elseif isfinite(medFree) && abs(medFree) < 0.5 * medD
            d.carrier = 'clock';
            fprintf('      over the %d epochs that count the reading is %.1f m while the rows with their clocks free put the shift on %s at %+.1f m: a common mode\n', numel(cnt), medD, letters(a), medFree);
            fprintf('      of the set carries it (the receiver clock, or the biases of a changed satellite set), read as position because the projection holds the clocks at the host''s P\n');
            fprintf('      (the sigma of the shift with the clocks free is %.1f m on %s against the CPI''s %.1f m: what a clock floor in S costs on this axis)\n', sgFree, letters(a), sg0);
        elseif isfinite(medFree)
            d.carrier = 'position';
            fprintf('      over the %d epochs that count the reading is %.1f m and the rows with their clocks free put the shift on %s at %+.1f m: the rows agree on a position shift\n', numel(cnt), medD, letters(a), medFree);
        else
            fprintf('      over the %d epochs that count the reading is %.1f m; no free-clock fit (too few range rows for a redundancy)\n', numel(cnt), medD);
        end
    end
    if haveRun
        ok = cnt(isfinite(navErrW(cnt - idx(1) + 1, a)));
        if ~isempty(ok)
            jj = ok - idx(1) + 1;
            d.navErr = median(navErrW(jj, a));  d.excess = median(abs(dHat(ok, a)) - abs(navErrW(jj, a)));
            sg = medF(sig(ok, a));
            if d.excess > 3 * sg
                d.verdict = 'ranges';
                fprintf('      the reading exceeds the logged nav''s error from the truth by %.1f m (median of |dHat| - |nav - truth| over those epochs, %.1f sigma): the ranges carry it\n', d.excess, d.excess / sg);
            elseif abs(d.excess) <= 1.5 * sg
                d.verdict = 'nav';
                fprintf('      the reading is within %.1f m of the logged nav''s error from the truth (median of |dHat| - |nav - truth| over those epochs, %.1f sigma): the nav carries it\n', abs(d.excess), abs(d.excess) / sg);
            else
                d.verdict = 'both';
                fprintf('      the reading exceeds the logged nav''s error from the truth by %.1f m (median of |dHat| - |nav - truth| over those epochs, %.1f sigma): the nav carries\n', d.excess, d.excess / sg);
                fprintf('      %.1f m of it and the ranges the rest\n', abs(d.navErr));
            end
        else
            fprintf('      no epoch of the window has both a reading and a truth sample\n');
        end
    end
    detail(end+1) = d; %#ok<AGROW>
end
nomW = kW;
if haveLog && isfield(spfLog, 'mode'), nomW = kW(round(spfLog.mode(kW)) == 1); end
if isempty(nomW), nomW = kW; end
sigRatio = medF(sig(nomW, :)) ./ medF(sigR(nomW, :));
if any(sigRatio > 1.5)
    fprintf('   NOTE: the sigma of the reading is %.1f x what R alone gives (medians over the NOMINAL epochs of the window): S = H P H'' + R carries the\n', max(sigRatio));
    fprintf('         host''s position sigma, and the CPI is normalised by S, so a P the host carries above its true error desensitises the monitor by that factor\n');
else
    fprintf('   NOTE: S is R-dominated on the NOMINAL epochs of the window (sigma %.2f x what R alone gives): the host''s P does not desensitise the monitor here\n', max(sigRatio));
end
if haveLog && numel(nomW) < numel(kW)
    cw = setdiff(kW, nomW);
    fprintf('   on the %d COAST / PROBATION epochs of the window the record''s S carries the host''s coast covariance unbounded (position sigma up to\n', numel(cw));
    fprintf('   %.1f m): the replay''s sums there are what the gate would read without the bound; the gate''s CPI against the bounded coast is the log column\n', ...
        max(max(sigP(cw, :))));
end
res = struct('t', tKf, 'dHat', dHat, 'sig', sig, 'sigR', sigR, 'xi', xi, 'qReplay', qReplay, 'qLog', qLog, 'sigP', sigP, ...
    'matches', matches, 'window', [tKf(kW(1)) tKf(kW(end))], 'axes', axes, 'rMode', rMode, 'rNominal', rNom, ...
    'minRows', minRows, 'nRangeRows', nRangeRows, 'peaks', {peaks}, 'dHatScreen', dHatS, 'xiScreen', xiS, 'qScreen', qScreen, ...
    'nDropped', nDropped, 'worstRow', worstRow, 'worstNorm', worstNorm, 'shiftFree', shiftFit, 'clockFree', {clkFit}, 'sigFree', sigFree, ...
    'clockStep', clockStep, 'dHatFloor', dHatF, 'sigFloor', sigF, 'xiFloor', xiF, 'qFloor', qFloor, 'rowK', ROW_K, 'detail', detail);
end

%------------------------------------------------------------------------
function [S, R] = covOf(s, m, useNom, rNom, useHost)
% the live blocks of S and R of one record; with useHost the pair the host applied; with useNom, R replaced by the
% nominal R of each row's type and S rebuilt
S = s.innovationCov(1:m, 1:m);  R = s.measNoiseCov(1:m, 1:m);
if useHost, S = s.innovationCovHost(1:m, 1:m);  R = s.measNoiseCovHost(1:m, 1:m); end
if useNom
    Rn = diag(rNom(rowTypes(s.obsMatrix(1:m, :))));
    S = S - R + Rn;  R = Rn;
end
end

function [typ] = rowTypes(H)
% 1 = range row (unit line of sight in the position columns, two entries or more), 2 = rate row, 3 = anything else
normP = sqrt(sum(H(:, 1:3).^2, 2));  normV = sqrt(sum(H(:, 4:6).^2, 2));
typ = 3 * ones(size(H, 1), 1);
typ(normP > 0.5 & normV < 0.1 & sum(abs(H(:, 1:3)) > 1e-9, 2) >= 2) = 1;
typ(normV > 0.5 & normP < 0.1) = 2;
end

%------------------------------------------------------------------------
function [drop, worst, worstNorm, worstEl, shift, clk, sysK, clkCols, sigFree] = rowScreen(H, y, R, rowK)
% the range rows of one epoch fitted by one shift [N E D] and one clock per system (every column beyond the velocity a
% range row uses), weighted by R; the row with the largest residual over its sigma is dropped while it exceeds rowK and
% one degree of redundancy is left after the drop. worst is the row with the largest residual of the first fit (its
% index in the record), worstNorm that residual over its sigma, worstEl the elevation its line of sight gives (deg);
% shift [N E D] and clk (one per clock column) are the last fit, on the rows left, NaN without one; sysK the range rows
% per clock column; clkCols the clock columns (state indices); sigFree the sigma of the fitted shift per axis
drop = zeros(1, 0);  worst = 0;  worstNorm = NaN;  worstEl = NaN;  shift = NaN(1, 3);  clk = zeros(1, 0);  sysK = zeros(1, 0);
clkCols = zeros(1, 0);  sigFree = NaN(1, 3);
typ = rowTypes(H);  keep = find(typ == 1)';
if numel(keep) < 2, return; end
clkCols = 6 + find(any(abs(H(keep, 7:end)) > 1e-9, 1));
cols = [1 2 3, clkCols];
sysK = zeros(1, numel(clkCols));
for c = 1:numel(clkCols), sysK(c) = nnz(abs(H(keep, clkCols(c))) > 1e-9); end
clk = NaN(1, numel(clkCols));
first = true;
while true
    nK = numel(keep);
    w = 1 ./ diag(R(keep, keep));
    if any(~(w > 0) | ~isfinite(w)), break; end
    Aw = H(keep, cols) .* repmat(sqrt(w), 1, numel(cols));  yw = y(keep) .* sqrt(w);
    rk = rank(Aw);
    if nK - rk < 1, break; end                               % no redundancy: nothing can be judged
    pA = pinv(Aw);  hat = Aw * pA;
    delta = pA * yw;  shift = delta(1:3)';  clk = delta(4:end)';
    cv = pinv(Aw' * Aw);  sigFree = sqrt(max(diag(cv(1:3, 1:3)), 0))';
    rw = yw - hat * yw;
    rn = rw ./ sqrt(max(1 - diag(hat), 1e-6));
    [mx, i] = max(abs(rn));
    if first
        worst = keep(i);  worstNorm = mx;  worstEl = asin(min(abs(H(worst, 3)), 1)) * 180 / pi;  first = false;
    end
    if mx > rowK && (nK - 1) - rk >= 1
        drop(end+1) = keep(i);  keep(i) = []; %#ok<AGROW>
    else
        break;
    end
end
end

function [matches] = matchVerdict(qCmp, qLog, both, kW, tKf, cmpName, tol)
% the replayed column against the gate log's cpiRatio on the NOMINAL armed epochs
dq = max(abs(qCmp(both, :) - qLog(both, :)), [], 2);
matches = max(dq) < tol;
if matches
    fprintf('   VERDICT: %s matches the log on the %d NOMINAL armed epochs of the window (largest difference %.3f): the gate was fed\n', cmpName, numel(both), max(dq));
    fprintf('            these records');
    if numel(both) < numel(kW), fprintf('; on the other %d epochs the gate was in COAST or PROBATION, where the log carries no CPI', numel(kW) - numel(both)); end
    fprintf('\n');
else
    fprintf('   VERDICT: %s differs from the log on %d of %d epochs (largest difference %.3f, first at %.2f s): the gate''s buffers do not\n', ...
        cmpName, nnz(dq >= tol), numel(both), max(dq), tKf(both(find(dq >= tol, 1))));
    fprintf('            hold these records'' innovations, or the gate reads another S (other floors in S, another R, or windows that skip epochs)\n');
end
end

function [SF] = floorOf(S, H, P, floors, clkCols)
% S with the gate's floors carried in: S + H dP H', dP the floor variance above the record's P on the floored states and
% on the clock columns the range rows use (the live rows only: S and H are m x m and m x n)
n = size(H, 2);  dP = zeros(n, 1);
for i = 1:numel(floors.states)
    st = floors.states(i);
    if st >= 1 && st <= n, dP(st) = max(dP(st), floors.sigmas(i)^2 - P(st, st)); end
end
if floors.clockSigma > 0
    for c = clkCols, dP(c) = max(dP(c), floors.clockSigma^2 - P(c, c)); end
end
dP = max(dP, 0);
SF = S + H * diag(dP) * H';
SF = (SF + SF') / 2;
end

function [m] = medF3(X)
% the median per column over the finite rows
m = NaN(1, size(X, 2));
for c = 1:size(X, 2), m(c) = medF(X(:, c)); end
end

function [e] = nedDelta(p, q)
% p minus q, both [lat lon h] rows, in NED metres (WGS84 radii at p), as runCompare forms it
L  = p(:, 1);  h = p(:, 3);
a  = 6378137;  e2 = 6.69437999014e-3;
RN = a * (1 - e2) ./ (1 - e2 * sin(L).^2).^1.5;
RE = a ./ sqrt(1 - e2 * sin(L).^2);
e = [(p(:, 1) - q(:, 1)) .* (RN + h), ...
     (p(:, 2) - q(:, 2)) .* (RE + h) .* cos(L), ...
    -(p(:, 3) - q(:, 3))];
end

function [y] = interpRows(tt, Y, tq)
y = NaN(numel(tq), size(Y, 2));
for c = 1:size(Y, 2)
    y(:, c) = interp1(tt(:), Y(:, c), tq(:), 'linear', NaN);
end
end

function [s] = num(x, fmt)
if isnan(x), s = '-'; else, s = sprintf(fmt, x); end
end

function [m] = medF(x)
x = x(isfinite(x));
if isempty(x), m = NaN; else, m = median(x); end
end
