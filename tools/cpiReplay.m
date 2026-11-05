function [res] = cpiReplay(spfMeasLog, tKf, spfLog, tWin, nominalBias, rMode, minRows)
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
%   spfLog may be [] (no comparison with the log); tWin [tOn tOff] in s on tKf, the minute before it is the baseline.
%   Per epoch and monitored axis, from the record's innovation gamma, S and H: dHat = f' S^-1 gamma / (f' S^-1 f), the
%   offset the CPI projection reads on the axis (m); sig = 1 / sqrt(f' S^-1 f), the sigma of that reading (m); sigR =
%   1 / sqrt(f' R^-1 f), the same with R alone (what the sigma would be if the host carried no position uncertainty in
%   S); xi = (|dHat| - B) / sig floored at zero, as SPF_cpiMonitor forms it; q/T_N of the ten epochs ending at the
%   epoch from SPF_cpiMonitor itself; the gate log's cpiRatio at the same epoch; and sigP = sqrt(P+) on the axis, the
%   position sigma the host carries. Prints the table over the window, the medians over the minute before and over
%   the window, whether the replay matches the log (then the gate was fed these records), and what sets the reading:
%   how large the offset read is, how many epochs it stays above B, and whether S is R or the host's P.
%   Result: t, dHat, sig, sigR, xi, qReplay, qLog [K x nAxes], sigP [K x nAxes], matches, window, axes, rMode,
%   rNominal, minRows, nRangeRows, peaks (a cell per axis of [t q] for every run of the ten-epoch sum above 0.25 in
%   the window, one peak per run: the near misses, not only the largest).
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
useNom = strcmpi(rMode, 'nominal');
useHost = strcmpi(rMode, 'host');
if useHost && ~isfield(spfMeasLog, 'measNoiseCovHost')
    error('cpiReplay:host', 'the records carry no measNoiseCovHost: the mode ''host'' needs records from setKfMeasBoth / kfMeasFromUpdateBoth');
end
BEFORE = 60;                                                 % s: the baseline minute before the window
MATCH_TOL = 0.02;                                            % in q/T_N: replay and log closer than this agree
N = double(CST_spfParam.WINDOW_LENGTH);  T = CST_spfParam.CPI_THRESHOLD;
axes = double(CST_spfParam.MONITORED_AXES);  nA = numel(axes);  letters = 'NED';
kW = find(tKf >= min(tWin) & tKf <= max(tWin));
kB = find(tKf >= min(tWin) - BEFORE & tKf < min(tWin));
if isempty(kW), error('cpiReplay:window', 'the window %g to %g s holds no epoch of the records (%.2f to %.2f s)', tWin, tKf(1), tKf(end)); end
kAll = [kB; kW];
dHat = NaN(K, nA);  sig = NaN(K, nA);  sigR = NaN(K, nA);  xi = zeros(K, nA);  sigP = NaN(K, nA);  qReplay = NaN(K, nA);  qLog = NaN(K, nA);
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
cmpEp = kW(all(isfinite(qReplay(kW, :)) & isfinite(qLog(kW, :)), 2));
if haveLog && isfield(spfLog, 'mode'), cmpEp = cmpEp(round(spfLog.mode(cmpEp)) == 1); end
if haveLog && isfield(spfLog, 'armed'), cmpEp = cmpEp(logical(spfLog.armed(cmpEp))); end
if haveLog && useNom
    both = cmpEp;
    if ~isempty(both)
        fprintf('   with the nominal R the largest ten-epoch sum over the window is %.3f (N %.3f, E %.3f, D %.3f); the log, with the host''s R, %.3f\n', ...
            max(max(qReplay(both, :))), max(qReplay(both, :)), max(max(qLog(both, :))));
    end
elseif haveLog
    both = cmpEp;
    if ~isempty(both)
        dq = max(abs(qReplay(both, :) - qLog(both, :)), [], 2);
        matches = max(dq) < MATCH_TOL;
        if matches
            fprintf('   VERDICT: the replay matches the log on the %d NOMINAL armed epochs of the window (largest difference %.3f): the gate was fed\n', numel(both), max(dq));
            fprintf('            these records');
            if numel(both) < numel(kW), fprintf('; on the other %d epochs the gate was in COAST or PROBATION, where the log carries no CPI', numel(kW) - numel(both)); end
            fprintf('\n');
        else
            fprintf('   VERDICT: the replay differs from the log on %d of %d epochs (largest difference %.3f, first at %.2f s): the gate''s buffers do not\n', ...
                nnz(dq >= MATCH_TOL), numel(both), max(dq), tKf(both(find(dq >= MATCH_TOL, 1))));
            fprintf('            hold these records'' innovations (another vector, another S, or windows that skip epochs)\n');
        end
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
sigRatio = medF(sig(kW, :)) ./ medF(sigR(kW, :));
if any(sigRatio > 1.5)
    fprintf('   NOTE: the sigma of the reading is %.1f x what R alone gives (medians over the window): S = H P H'' + R carries the host''s position\n', max(sigRatio));
    fprintf('         sigma, and the CPI is normalised by S, so a P the host carries above its true error desensitises the monitor by that factor\n');
else
    fprintf('   NOTE: S is R-dominated on the window (sigma %.2f x what R alone gives): the host''s P does not desensitise the monitor here\n', max(sigRatio));
end
res = struct('t', tKf, 'dHat', dHat, 'sig', sig, 'sigR', sigR, 'xi', xi, 'qReplay', qReplay, 'qLog', qLog, 'sigP', sigP, ...
    'matches', matches, 'window', [tKf(kW(1)) tKf(kW(end))], 'axes', axes, 'rMode', rMode, 'rNominal', rNom, ...
    'minRows', minRows, 'nRangeRows', nRangeRows, 'peaks', {peaks});
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
function [s] = num(x, fmt)
if isnan(x), s = '-'; else, s = sprintf(fmt, x); end
end

function [m] = medF(x)
x = x(isfinite(x));
if isempty(x), m = NaN; else, m = median(x); end
end
