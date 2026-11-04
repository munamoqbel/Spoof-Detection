function [res] = measNoiseCheck(spfMeasLog, tKf, tWin)
%MEASNOISECHECK  The measurement noise the gate was handed, epoch by epoch: is R scaled on large innovations, are rows missing.
%   res = measNoiseCheck(spfMeasLog, tKf, tWin)   tWin [tOn tOff] in s on tKf; the minute before it is the baseline
%   The nominal R is the run median of the logged R diagonal per row type (range, rate, other). Per epoch of the window:
%   the rows in the record, the range rows, the largest sqrt(R) on the range rows, the range rows with R above 1.5 x the
%   nominal, and the innovation of the row with the largest R; with records from setKfMeasBoth / kfMeasFromUpdateBoth
%   the same for the R the host applied (measNoiseCovHost) and the rows the host did not apply (a host variance of
%   1e6 or more). The first 8 epochs of the window print every range row's sqrt(R), the host's beside when present.
%   Verdicts: the gate's R scaled in the window (the host hands the gate the R it applies: the interface change is not
%   in) or at the nominal while the host's R is scaled (the change is in), and range rows missing against the minute
%   before (rows the host rejected are not in the record).
%   Result: t, nRows, nRange, maxSqrtR, nScaled, innAtMaxR, maxSqrtRHost, nScaledHost, nRejectedHost [K x 1],
%   rNominal [1 x 3] (variances), window, hasHost, gateScaled, hostScaled, nMissing.
%   MATLAB desktop tool, not codegen. Runs under Octave.

if iscell(spfMeasLog), spfMeasLog = [spfMeasLog{:}]; end
K = numel(spfMeasLog);
if nargin < 2 || isempty(tKf), tKf = 0.5 * (0:K-1)'; end
tKf = tKf(:);
if nargin < 3 || isempty(tWin), tWin = [tKf(1) tKf(end)]; end
BEFORE = 60;                                                 % s: the baseline minute before the window
SCALED = 1.5;                                                % a variance above this times the nominal is scaled
REJECTED = 1e6;                                              % m^2: a host variance at or above this is a row the host did not apply
kW = find(tKf >= min(tWin) & tKf <= max(tWin));
if isempty(kW), error('measNoiseCheck:window', 'the window %g to %g s holds no epoch of the records (%.2f to %.2f s)', tWin, tKf(1), tKf(end)); end
kB = find(tKf >= min(tWin) - BEFORE & tKf < min(tWin));
hasHost = isfield(spfMeasLog, 'measNoiseCovHost');
%% the nominal R per row type over the run
acc = {[], [], []};
for k = 1:K
    s = spfMeasLog(k);  m = double(s.numMeas);
    if m < 1, continue; end
    typ = rowTypes(s.obsMatrix(1:m, :));  r = diag(s.measNoiseCov(1:m, 1:m));
    for c = 1:3, acc{c} = [acc{c}; r(typ == c)]; end
end
rNom = NaN(1, 3);
for c = 1:3, rNom(c) = medF(acc{c}); end
%% per epoch
nRows = NaN(K, 1);  nRange = NaN(K, 1);  maxSqrtR = NaN(K, 1);  nScaled = NaN(K, 1);  innAtMaxR = NaN(K, 1);
maxSqrtRHost = NaN(K, 1);  nScaledHost = NaN(K, 1);  nRejectedHost = NaN(K, 1);
rowR = cell(K, 1);  rowRH = cell(K, 1);
for k = [kB; kW]'
    s = spfMeasLog(k);  m = double(s.numMeas);
    nRows(k) = m;
    if m < 1, nRange(k) = 0;  continue; end
    typ = rowTypes(s.obsMatrix(1:m, :));  rg = find(typ == 1);
    nRange(k) = numel(rg);
    if isempty(rg), continue; end
    r = diag(s.measNoiseCov(1:m, 1:m));  y = s.innovation(1:m);
    [rMax, i] = max(r(rg));
    maxSqrtR(k) = sqrt(rMax);  nScaled(k) = nnz(r(rg) > SCALED * rNom(1));  innAtMaxR(k) = y(rg(i));
    rowR{k} = sqrt(r(rg))';
    if hasHost
        rh = diag(s.measNoiseCovHost(1:m, 1:m));
        live = rh(rg) < REJECTED;
        nRejectedHost(k) = nnz(~live);
        nScaledHost(k) = nnz(rh(rg) > SCALED * rNom(1) & live);
        if any(live), maxSqrtRHost(k) = sqrt(max(rh(rg(live)))); else, maxSqrtRHost(k) = NaN; end
        rowRH{k} = sqrt(rh(rg))';
    end
end
%% print
fprintf('measNoiseCheck: %d records; window %.2f to %.2f s (epochs %d to %d), the minute before as the baseline\n', K, tKf(kW(1)), tKf(kW(end)), kW(1), kW(end));
fprintf('   nominal R (run median per row type): range rows %.2f m, rate rows %.3f m/s, other rows %.2f; a range row counts as scaled above %.1f x the\n', ...
    sqrt(rNom(1)), sqrt(rNom(2)), sqrt(rNom(3)), SCALED);
fprintf('   nominal variance');
if hasHost, fprintf('; the records carry the R the host applied beside the gate''s (a host variance of %.0e or more: a row the host did not apply)', REJECTED); end
fprintf('\n');
pick = unique([kW(1:min(8, numel(kW)))', kW(round(linspace(1, numel(kW), min(12, numel(kW)))))']);
fprintf('   %10s %5s %5s %12s %7s %10s', 't [s]', 'rows', 'range', 'max sqrt(R)', 'scaled', 'inn there');
if hasHost, fprintf(' | %12s %7s %8s', 'host sqrt(R)', 'scaled', 'left out'); end
fprintf('\n');
for k = pick
    fprintf('   %10.2f %5s %5s %12s %7s %10s', tKf(k), num(nRows(k), '%d'), num(nRange(k), '%d'), num(maxSqrtR(k), '%.1f'), num(nScaled(k), '%d'), num(innAtMaxR(k), '%+.1f'));
    if hasHost, fprintf(' | %12s %7s %8s', num(maxSqrtRHost(k), '%.1f'), num(nScaledHost(k), '%d'), num(nRejectedHost(k), '%d')); end
    fprintf('\n');
end
fprintf('   sqrt(R) of every range row on the first %d epochs of the window:\n', min(8, numel(kW)));
for k = kW(1:min(8, numel(kW)))'
    fprintf('   %10.2f gate:%s\n', tKf(k), sprintf(' %5.1f', rowR{k}));
    if hasHost, fprintf('   %10s host:%s\n', '', sprintf(' %5.1f', min(rowRH{k}, 999.9))); end
end
%% verdicts
inW = kW(isfinite(nScaled(kW)));
gateScaled = any(nScaled(inW) > 0);
hostScaled = hasHost && any(nScaledHost(inW) > 0 | nRejectedHost(inW) > 0);
rowsBefore = medF(nRange(kB));
nMissing = nnz(nRange(kW) < rowsBefore - 1 & nRange(kW) > 0);
if gateScaled
    kF = inW(find(nScaled(inW) > 0, 1));
    fprintf('   VERDICT: the gate''s R is scaled above the nominal on %d of %d window epochs (first at %.2f s, up to %.1f m against %.1f m nominal)', ...
        nnz(nScaled(inW) > 0), numel(inW), tKf(kF), max(maxSqrtR(inW)), sqrt(rNom(1)));
    if hasHost
        fprintf(':\n            the records carry a host R, yet the gate''s R is not the nominal one: the host passed its scaled R as the gate''s\n');
    else
        fprintf(':\n            the host hands the gate the R it applies (setKfMeas / kfMeasFromUpdate); the gate''s S must carry the nominal R on all\n');
        fprintf('            the rows the host built, the host''s own R beside it (setKfMeasBoth / kfMeasFromUpdateBoth)\n');
    end
elseif hasHost
    fprintf('   VERDICT: the gate''s R stays at the nominal on every window epoch while the host''s R is scaled on %d epochs (up to %.1f m) and the host\n', ...
        nnz(nScaledHost(inW) > 0), max([maxSqrtRHost(inW); 0]));
    fprintf('            left out %d row-epochs: the interface change is in\n', sum(nRejectedHost(inW)));
else
    fprintf('   VERDICT: R stays at the nominal on every window epoch: no scaling in this window\n');
end
if nMissing > 0 && ~hasHost
    kM = kW(find(nRange(kW) < rowsBefore - 1 & nRange(kW) > 0, 1));
    fprintf('   fewer range rows than the minute before (median %d) on %d window epochs, first at %.2f s: rows the host rejected are not in the record\n', ...
        round(rowsBefore), nMissing, tKf(kM));
end
res = struct('t', tKf, 'nRows', nRows, 'nRange', nRange, 'maxSqrtR', maxSqrtR, 'nScaled', nScaled, 'innAtMaxR', innAtMaxR, ...
    'maxSqrtRHost', maxSqrtRHost, 'nScaledHost', nScaledHost, 'nRejectedHost', nRejectedHost, 'rNominal', rNom, ...
    'window', [tKf(kW(1)) tKf(kW(end))], 'hasHost', hasHost, 'gateScaled', gateScaled, 'hostScaled', hostScaled, 'nMissing', nMissing);
end

%------------------------------------------------------------------------
function [typ] = rowTypes(H)
% 1 = range row (unit line of sight in the position columns, two entries or more), 2 = rate row, 3 = anything else
normP = sqrt(sum(H(:, 1:3).^2, 2));  normV = sqrt(sum(H(:, 4:6).^2, 2));
typ = 3 * ones(size(H, 1), 1);
typ(normP > 0.5 & normV < 0.1 & sum(abs(H(:, 1:3)) > 1e-9, 2) >= 2) = 1;
typ(normV > 0.5 & normP < 0.1) = 2;
end

function [s] = num(x, fmt)
if isnan(x), s = '-'; else, s = sprintf(fmt, x); end
end

function [m] = medF(x)
x = x(isfinite(x));
if isempty(x), m = NaN; else, m = median(x); end
end
