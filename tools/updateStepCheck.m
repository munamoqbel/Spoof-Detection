function [res] = updateStepCheck(spfMeasLog, tKf, tWin)
%UPDATESTEPCHECK  Is the state increment at each update the Kalman step that the logged P and R imply?
%   res = updateStepCheck(spfMeasLog)             all epochs
%   res = updateStepCheck(spfMeasLog, tKf, tWin)  epochs with tKf in [tWin(1) tWin(2)]
%
%   For a standard Kalman update K y = P+ H' R^-1 y, whatever form the host computes it in
%   (P+ the posterior covariance, R the measurement noise the gain was built with). The spfMeas
%   records carry P+ (postCov), H, y, the R the host handed over (measNoiseCovHost when present,
%   else measNoiseCov) and the state before and after the update, so per epoch they give the
%   step the update should have added to the state with that P and that R, and the increment it
%   did add. The host's own S and K are not in the records; this is the test that does not need
%   them. With the single-R record call (STRUCT_SPF.kfMeasFromUpdate) measNoiseCovHost is a
%   copy of the gate's R, so a scaling the host applies inside its update is not on record and
%   shows here only through the row-type regression below; the Both call puts the applied R on
%   record and makes the slope the direct test.
%
%   Row types from H with the host's state layout (position 1:3, velocity 4:6, attitude 7:9,
%   gyro bias 10:12, accel bias 13:15, the rest clocks and others, NED error states in m, m/s,
%   rad): a row with a unit velocity partial is a range-rate row, a unit position partial a range
%   row, a single position entry the pressure row.
%
%   Printed per state group (the six groups above, then one line per state from 16 on): rms of
%   the Kalman step and of the recorded increment, the slope of
%   the recorded increment on the Kalman step (1 when the update is whole) and the fit of that
%   line (fraction of the increment's variance it explains). Then the split by row type: the
%   step is linear in R^-1, so an update that ran with R_used = c R on a row type adds a = 1/c of
%   the step that row type alone would give. Per group the recorded increment is regressed on the
%   three steps (range rows, range-rate rows, pressure row); a coefficient is reported only when
%   its row type carries a material share of the group's step. The same a across the groups is an
%   update that ran with another R than recorded, by c = 1/a per row type; coefficients that
%   agree within each group but differ between the groups are a factor per state; the rest is a
%   P in the gain that is not the P+ logged, or something between K y and the state.
%
%   Records written with STRUCT_SPF.kfMeasFromUpdateBoth carry the host's S and R beside the
%   gate's, and then a second block checks the covariance step: H P- H' is the host's S minus
%   its R, the logged H P+ H' comes from postCov, and a standard update would leave
%   H P- H' - H P- H' S^-1 H P- H' with S built from the applied R or from the nominal R. The
%   share of the predicted reduction that the logged P+ took, per row type, says which R the
%   covariance was updated with: the applied one (covariance and state consistent) or the
%   nominal one (the covariance books information the state did not use).
%
%   Verdict lines:
%     whole      every group with signal has a slope within 0.9..1.1 and a fit above 0.9
%     R          the row-type coefficients agree across the groups within 25 % and fit above 0.8:
%                the gain ran with R scaled by c per row type
%     per state  the coefficients agree within each group but not across them
%     otherwise  not a per-row-type scaling of the step
%
%   res: t, step [K x n] (Kalman step per epoch and state), incr [K x n] (recorded increment),
%   groups {name, idx, unit}, slope, fit, rmsStep, rmsIncr (per group), aRow [G x 3] (range,
%   range-rate, pressure; NaN where immaterial), fitRow, cRow (1 ./ median material a per row
%   type), verdict ('whole' | 'R' | 'per state' | 'other' | 'no signal'), covShare [2 x 3] (share of
%   the predicted reduction taken, rows applied / nominal R, columns range / range-rate / pressure),
%   covVerdict ('applied' | 'nominal' | 'neither' | 'same R' | 'no host S on record').
%   MATLAB desktop tool, not codegen. Runs under Octave.

if iscell(spfMeasLog), spfMeasLog = [spfMeasLog{:}]; end
K = numel(spfMeasLog);
if nargin < 2 || isempty(tKf), tKf = (1:K)'; end
if nargin < 3 || isempty(tWin), tWin = [-Inf Inf]; end
t = tKf(:);
n = numel(spfMeasLog(1).postState);
useHostR = isfield(spfMeasLog, 'measNoiseCovHost');

SHARE_MIN = 0.25;          % a row type is material for a group when its step is at least this share of the group's total step
AGREE = 0.25;              % coefficients agree when their spread is within this fraction of their median
MG = 9.80665e-3;

step = NaN(K, n);  incr = NaN(K, n);
bRange = NaN(K, n);  bRate = NaN(K, n);  bBaro = NaN(K, n);
nRows = zeros(K, 3);
hasBoth = useHostR && isfield(spfMeasLog, 'innovationCovHost');
cvPrior = NaN(K, 3);  cvLogged = NaN(K, 3);  cvApplied = NaN(K, 3);  cvNominal = NaN(K, 3);
shApplied = NaN(K, 3);  shNominal = NaN(K, 3);  rDiffers = false(K, 1);
for k = 1:K
    if t(k) < tWin(1) || t(k) > tWin(2), continue; end
    s = spfMeasLog(k);
    m = double(s.numMeas);
    if m == 0, continue; end
    H = s.obsMatrix(1:m, :);
    y = s.innovation(1:m);  y = y(:);
    if useHostR, R = s.measNoiseCovHost(1:m, 1:m); else, R = s.measNoiseCov(1:m, 1:m); end
    P = s.postCov;
    normP = sqrt(sum(H(:, 1:3).^2, 2));
    normV = sqrt(sum(H(:, 4:6).^2, 2));
    isRr = normV > 0.5 & normV > normP;
    isPr = normP > 0.5 & normV < 0.1;
    isBaro = isPr & sum(abs(H(:, 1:3)) > 1e-9, 2) == 1;
    isPr = isPr & ~isBaro;
    nRows(k, :) = [nnz(isPr), nnz(isRr), nnz(isBaro)];
    step(k, :) = (P * H' * (R \ y))';
    incr(k, :) = (s.postState(:) - s.priorState(:))';
    bRange(k, :) = typeStep(P, H, R, y, isPr);
    bRate(k, :)  = typeStep(P, H, R, y, isRr);
    bBaro(k, :)  = typeStep(P, H, R, y, isBaro);
    if hasBoth
        % the covariance step: H P- H' from the host's S minus its R, the logged H P+ H', and the H P+ H' a
        % standard update would leave with the applied R and with the nominal R
        Rn = s.measNoiseCov(1:m, 1:m);
        Sh = s.innovationCovHost(1:m, 1:m);  Sn = s.innovationCov(1:m, 1:m);
        rDiffers(k) = max(abs(diag(R) - diag(Rn))) > 1e-9 * max(1, max(abs(diag(Rn))));
        HPmH = Sh - R;  HPmH = (HPmH + HPmH') / 2;
        A = HPmH - HPmH * (Sh \ HPmH);  B = HPmH - HPmH * (Sn \ HPmH);  L = H * P * H';
        dPrior = diag(HPmH);  dL = diag(L);  dA = diag(A);  dB = diag(B);
        sel = {isPr, isRr, isBaro};
        for c = 1:3
            i = sel{c};
            if ~any(i), continue; end
            cvPrior(k, c) = median(dPrior(i));  cvLogged(k, c) = median(dL(i));
            cvApplied(k, c) = median(dA(i));    cvNominal(k, c) = median(dB(i));
            redA = dPrior(i) - dA(i);  redB = dPrior(i) - dB(i);  redL = dPrior(i) - dL(i);
            okA = redA > 1e-12 * max(dPrior(i));  okB = redB > 1e-12 * max(dPrior(i));
            if any(okA), shApplied(k, c) = median(redL(okA) ./ redA(okA)); end
            if any(okB), shNominal(k, c) = median(redL(okB) ./ redB(okB)); end
        end
    end
end
valid = all(isfinite(step), 2) & all(isfinite(incr), 2);

groups = {'position 1:3',   1:3,   1,               'm'; ...
          'velocity 4:6',   4:6,   1,               'm/s'; ...
          'tilt 7:8',       7:8,   1e3,             'mrad'; ...
          'heading 9',      9,     1e3,             'mrad'; ...
          'gyro bias 10:12', 10:12, 180 / pi * 3600, 'deg/h'; ...
          'accel bias 13:15', 13:15, 1 / MG,         'mg'};
keep = false(size(groups, 1), 1);
for g = 1:size(groups, 1)
    idx = groups{g, 2};
    keep(g) = ~isempty(idx) && max(idx) <= n;
end
groups = groups(keep, :);
for j = 16:n                                                 % clocks and the rest: one line each, so no units are mixed
    groups(end + 1, :) = {sprintf('state %d', j), j, 1, ''};
end
G = size(groups, 1);

slope = NaN(G, 1);  fit = NaN(G, 1);  rmsStep = NaN(G, 1);  rmsIncr = NaN(G, 1);
aRow = NaN(G, 3);  fitRow = NaN(G, 1);  share = NaN(G, 3);  signal = false(G, 1);
for g = 1:G
    idx = groups{g, 2};
    X = step(valid, idx);  D = incr(valid, idx);
    x = X(:);  d = D(:);
    if numel(x) < 10, continue; end
    rmsStep(g) = sqrt(mean(x.^2));  rmsIncr(g) = sqrt(mean(d.^2));
    if rmsStep(g) <= 0 || rmsIncr(g) <= 0, continue; end
    signal(g) = true;
    slope(g) = (x' * d) / (x' * x);
    fit(g) = 1 - sum((d - slope(g) * x).^2) / sum(d.^2);
    B = [reshape(bRange(valid, idx), [], 1), reshape(bRate(valid, idx), [], 1), reshape(bBaro(valid, idx), [], 1)];
    B(~isfinite(B)) = 0;
    for j = 1:3, share(g, j) = sqrt(mean(B(:, j).^2)) / rmsStep(g); end
    use = share(g, :) > SHARE_MIN;
    if any(use)
        a = B(:, use) \ d;
        aRow(g, use) = a';
        fitRow(g) = 1 - sum((d - B(:, use) * a).^2) / sum(d.^2);
    end
end

% verdict
sig = find(signal);
verdict = 'no signal';
cRow = NaN(1, 3);
if ~isempty(sig)
    whole = all(slope(sig) > 0.9 & slope(sig) < 1.1 & fit(sig) > 0.9);
    % across the groups, per row type (material coefficients only)
    acrossOk = true;  anyAcross = false;
    for j = 1:3
        a = aRow(sig, j);  a = a(isfinite(a));
        if numel(a) >= 2
            anyAcross = true;
            if (max(a) - min(a)) > AGREE * abs(median(a)), acrossOk = false; end
        end
        if ~isempty(a), cRow(j) = 1 / median(a); end
    end
    % within each group, across its material row types
    withinOk = true;  anyWithin = false;
    for g = sig'
        a = aRow(g, :);  a = a(isfinite(a));
        if numel(a) >= 2
            anyWithin = true;
            if (max(a) - min(a)) > AGREE * abs(median(a)), withinOk = false; end
        end
    end
    fitsOk = all(fitRow(sig) > 0.8);
    if whole
        verdict = 'whole';
    elseif anyAcross && acrossOk && fitsOk
        verdict = 'R';
    elseif anyWithin && withinOk && fitsOk && ...
            (max(slope(sig)) - min(slope(sig))) > AGREE * abs(median(slope(sig)))
        verdict = 'per state';
    else
        verdict = 'other';
    end
end

% print
nv = nnz(valid);
medRows = [0 0 0];
if nv > 0, medRows = round(median(nRows(valid, :), 1)); end
fprintf('updateStepCheck: %d epochs with rows in the window, %d states; rows per epoch (median) %d range, %d range-rate, %d pressure\n', ...
    nv, n, medRows(1), medRows(2), medRows(3));
fprintf('   step = P+ H'' R^-1 y from the record, R = %s; increment = postState - priorState\n', ...
    pick(useHostR, 'measNoiseCovHost (the host''s R as handed over)', 'measNoiseCov'));
fprintf('                         Kalman step   recorded incr   slope    fit  |  by row type: a_range  a_rate  a_baro    fit\n');
for g = 1:G
    sc = groups{g, 3};  un = groups{g, 4};
    if ~signal(g)
        fprintf('   %-18s %s\n', groups{g, 1}, pick(isfinite(rmsStep(g)) && rmsStep(g) == 0, 'no step (rows do not reach these states)', 'no signal'));
        continue
    end
    fprintf('   %-18s %9.4f %-5s %9.4f %-5s %6.2f  %5.2f  |              %s  %s  %s   %5.2f\n', groups{g, 1}, ...
        rmsStep(g) * sc, un, rmsIncr(g) * sc, un, slope(g), fit(g), ...
        fmtA(aRow(g, 1)), fmtA(aRow(g, 2)), fmtA(aRow(g, 3)), fitRow(g));
end
fprintf('   (slope 1 and fit near 1: the update is whole; a is 1/c for an update that ran with R scaled by c on that row type;\n');
fprintf('    a coefficient is left out when its row type carries less than %.0f%% of the group''s step)\n', 100 * SHARE_MIN);
switch verdict
    case 'whole'
        fprintf('   VERDICT: the state receives the Kalman step the logged P and R imply, on every group with signal\n');
    case 'R'
        fprintf('   VERDICT: the increments are the Kalman step for an R scaled per row type: range rows x %s, range-rate rows x %s, pressure x %s.\n', ...
            fmtC(cRow(1)), fmtC(cRow(2)), fmtC(cRow(3)));
        fprintf('            The gain was built with that R, not with the R handed to the record\n');
    case 'per state'
        fprintf('   VERDICT: each group takes one fraction of its whole Kalman step (the slopes above), the same for every row type:\n');
        fprintf('            a factor per state sits between K y and the state, not R\n');
    case 'other'
        fprintf('   VERDICT: the increments are not the Kalman step, nor a per-row-type scaling of it, nor one factor per state:\n');
        fprintf('            the P the gain is built with is not the P+ logged, or something between K y and the state is not a constant factor\n');
    otherwise
        fprintf('   VERDICT: no epoch with rows and a finite step in the window\n');
end

% the covariance step, when the records carry the host's S and R beside the gate's
covVerdict = 'no host S on record';
covShare = NaN(2, 3);
if hasBoth && nv > 0
    if ~any(rDiffers(valid))
        covVerdict = 'same R';
        fprintf('   covariance step: the host''s R on record equals the gate''s on every epoch, so the P''s R cannot be told apart here\n');
    else
        typeName = {'range rows', 'range-rate rows', 'pressure row'};
        typeUnit = {'m', 'm/s', 'm'};
        fprintf('   covariance step: H P- H'' from the host''s S minus its R, the logged H P+ H'', and the H P+ H'' a standard update leaves\n');
        fprintf('   with the applied R and with the nominal R (sqrt of the median diagonal; share = the logged reduction over the predicted one)\n');
        fprintf('                         prior        logged       applied R: predicted  share    nominal R: predicted  share\n');
        for c = 1:3
            v = valid & isfinite(cvPrior(:, c));
            if ~any(v), continue; end
            covShare(1, c) = median(shApplied(v, c));  covShare(2, c) = median(shNominal(v, c));
            fprintf('      %-16s %8.4f %-4s %8.4f %-4s %8.4f %-4s          %5.2f    %8.4f %-4s          %5.2f\n', typeName{c}, ...
                sqrt(median(cvPrior(v, c))), typeUnit{c}, sqrt(median(cvLogged(v, c))), typeUnit{c}, ...
                sqrt(median(cvApplied(v, c))), typeUnit{c}, covShare(1, c), sqrt(median(cvNominal(v, c))), typeUnit{c}, covShare(2, c));
        end
        sA = covShare(1, isfinite(covShare(1, :)));  sB = covShare(2, isfinite(covShare(2, :)));
        if ~isempty(sA) && all(abs(sA - 1) < 0.15)
            covVerdict = 'applied';
            fprintf('   VERDICT: the logged P+ is the posterior of an update with the applied R: the covariance and the state used the same R\n');
        elseif ~isempty(sB) && all(abs(sB - 1) < 0.15)
            covVerdict = 'nominal';
            fprintf('   VERDICT: the logged P+ is the posterior of an update with the NOMINAL R: the covariance books the information of an R\n');
            fprintf('            the state did not use, and shrinks faster than the state moves\n');
        else
            covVerdict = 'neither';
            fprintf('   VERDICT: the logged P+ matches neither posterior: the P update does not follow the standard form with either R\n');
        end
    end
end

res = struct('t', t, 'step', step, 'incr', incr, 'groups', {groups}, 'slope', slope, 'fit', fit, ...
    'rmsStep', rmsStep, 'rmsIncr', rmsIncr, 'aRow', aRow, 'fitRow', fitRow, 'share', share, 'cRow', cRow, ...
    'signal', signal, 'valid', valid, 'nRows', nRows, 'verdict', verdict, ...
    'covShare', covShare, 'covVerdict', covVerdict, 'shApplied', shApplied, 'shNominal', shNominal);
end

%------------------------------------------------------------------------
function [b] = typeStep(P, H, R, y, rows)
% the Kalman step the rows of one type alone would give: P+ H_type' R_type^-1 y_type
n = size(H, 2);
if ~any(rows)
    b = zeros(1, n);
    return
end
b = (P * H(rows, :)' * (R(rows, rows) \ y(rows)))';
end

function [s] = fmtA(a)
if isfinite(a), s = sprintf('%7.3f', a); else, s = '      -'; end
end

function [s] = fmtC(c)
if isfinite(c), s = sprintf('%.1f', c); else, s = '-'; end
end

function [out] = pick(cond, a, b)
if cond, out = a; else, out = b; end
end
