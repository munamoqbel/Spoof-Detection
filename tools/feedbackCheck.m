function [res] = feedbackCheck(run, spfMeasLog, win)
%FEEDBACKCHECK  Does the KF update reach the mechanisation, when, and in full? The change of the
%   logged nav solution around each 2 Hz epoch against the increment the filter computed, per NED
%   axis, with the instant and the profile of the change read from the log itself.
%   res = feedbackCheck(run, spfMeasLog)         every epoch
%   res = feedbackCheck(run, spfMeasLog, win)    epochs with win(1) <= tKf <= win(2)
%
%   run: the struct runCompare takes: t [N x 1] host time of the nav samples (100 Hz), navPos
%   [N x 3] [lat lon h] (rad, rad, m), navVel [N x 3] NED (m/s), tKf [K x 1] the 2 Hz epochs,
%   one per record of spfMeasLog; biasEst [N x 6] or [K x 6] optional, accelerometer then gyro
%   (m/s^2, rad/s), the value the mechanisation subtracts plus the bias state; pGyroFirst
%   optional (true when the gyro biases are states 10:12 and the accelerometer biases 13:15).
%   spfMeasLog: struct array or cell of the per-epoch spfMeas records (priorState, postState,
%   numMeas).
%
%   Position. Per sample the change of navPos (lat lon h to NED metres) minus the integral of
%   navVel over the sample is the part of the motion the velocity does not explain; summed over
%   the half epoch either side of tKf it is the correction the update put into the logged
%   position, whenever and however it was applied. Its running sum from tKf - 0.25 s, regressed
%   on the increment postState - priorState (1:3) across the epochs, is the step profile: the
%   share of the increment present in the position at tKf + tau, from which the instant of the
%   update (10, 50 and 90 % of the step) is read. The sum over the whole half epoch is regressed
%   jointly on the increments of the records k-1, k and k+1: the coefficient at k reads -1 when
%   the increment is applied in full with the error state as estimate minus truth (+1 for the
%   opposite convention) and 0 when it never reaches the logged solution; a coefficient at k-1
%   or k+1 instead says the records and the epochs are offset by one.
%   Velocity. The same, with navVel minus the line fitted to it 0.1 to 0.3 s before the epoch
%   (the line removes the vehicle's acceleration), read at tKf + tau for the profile and
%   averaged over tKf + 0.15 to 0.25 s for the change across the epoch.
%   Bias. With biasEst the sum of the accelerometer bias increments over the window is set
%   against the change of biasEst, the same -1 / 0 reading for the bias feedback.
%
%   Verdict lines name an axis whose increment does not reach the logged solution, one applied
%   in part, with the opposite sign to the other axes or over more than one epoch, records offset
%   by one epoch from the nav log, and a step spread over more than 0.05 s (a slewed feedback).
%
%   res: per-epoch t, incPos, incVel, chgPos, chgVel [K x 3] (NaN where not computed), used
%   [K x 1], the lag coefficients coefPos, coefVel [3 x 3] (rows k-1, k, k+1; columns N E D),
%   the correlations rPos, rVel [1 x 3], the profiles tau [1 x J], profPos, profVel [J x 3],
%   the step times t10, t50, t90 [2 x 3] (rows position, velocity), the bias totals incBias and
%   dBias [1 x 3] in m/s^2 (NaN without biasEst), and the verdict lines in res.verdict (cell).
%   MATLAB desktop tool, not codegen. Runs under Octave.

if iscell(spfMeasLog), spfMeasLog = [spfMeasLog{:}]; end
if nargin < 3 || isempty(win), win = [-Inf Inf]; end
K   = numel(spfMeasLog);
t   = run.t(:);
tKf = run.tKf(:);
if numel(tKf) ~= K
    error('feedbackCheck: tKf has %d epochs, spfMeasLog %d records', numel(tKf), K);
end
HALF = 0.25;                                             % half an epoch: the interval owned by each update (s)
TB1 = 0.3;  TB2 = 0.1;                                   % line fit before the epoch, tKf - TB1 to tKf - TB2 (s)
TA1 = 0.15;  TA2 = 0.25;                                 % the velocity change is read over tKf + TA1 to tKf + TA2
MG  = 9.80665e-3;                                        % m/s^2 per mg
tau = -HALF:0.01:HALF;  J = numel(tau);
dtNav  = median(diff(t));
navNed = nedFrom(run.navPos);
navVel = run.navVel;
% the position change per sample that navVel does not explain, stamped at the later sample
resPos = [NaN(1, 3); diff(navNed) - 0.5 * (navVel(1:end-1, :) + navVel(2:end, :)) .* repmat(diff(t), 1, 3)];
iAcc = 10:12;
if isfield(run, 'pGyroFirst') && ~isempty(run.pGyroFirst) && run.pGyroFirst, iAcc = 13:15; end

incPos = NaN(K, 3);  incVel = NaN(K, 3);  incAcc = NaN(K, 3);
chgPos = NaN(K, 3);  chgVel = NaN(K, 3);
cumPos = NaN(K, J, 3);  devVel = NaN(K, J, 3);
used = false(K, 1);  inWin = false(K, 1);
for k = 1:K
    s = spfMeasLog(k);
    dx = s.postState(:) - s.priorState(:);
    if all(isfinite(dx)) && numel(dx) >= max(iAcc)
        incAcc(k, :) = dx(iAcc)';                        % every finite epoch: the bias total needs them all
    end
    inWin(k) = tKf(k) >= win(1) && tKf(k) <= win(2);
    if ~inWin(k) || ~all(isfinite(dx(1:6))), continue; end
    iI = find(t >= tKf(k) - HALF & t < tKf(k) + HALF);  % the samples of this update's half epochs
    iB = find(t >= tKf(k) - TB1 & t <= tKf(k) - TB2);   % the line before
    if numel(iI) < 10 || numel(iB) < 4 || iI(1) < 2, continue; end
    incPos(k, :) = dx(1:3)';
    incVel(k, :) = dx(4:6)';
    r = resPos(iI, :);  r(~isfinite(r)) = 0;
    chgPos(k, :) = sum(r, 1);
    tI = t(iI) - tKf(k);
    cs = cumsum(r, 1);
    for j = 1:J
        n = find(tI <= tau(j), 1, 'last');
        if ~isempty(n), cumPos(k, j, :) = cs(n, :); else, cumPos(k, j, :) = 0; end
    end
    A  = [ones(numel(iB), 1), t(iB) - tKf(k)];
    c  = A \ navVel(iB, :);                              % the line before: value and slope at the epoch
    iV = find(t >= tKf(k) - HALF & t <= tKf(k) + HALF + 1e-9);
    tV = t(iV) - tKf(k);
    dv = navVel(iV, :) - (repmat(c(1, :), numel(iV), 1) + tV * c(2, :));
    devVel(k, :, :) = interp1(tV, dv, tau, 'linear', NaN);
    chgVel(k, :) = mean(dv(tV >= TA1 & tV <= TA2, :), 1);
    used(k) = true;
end
nUsed = nnz(used);

[coefPos, rPos] = lagFit(incPos, chgPos, used);
[coefVel, rVel] = lagFit(incVel, chgVel, used);
profPos = profileOf(incPos, cumPos, used);              % [J x 3] slope of the running sum on the increment
profVel = profileOf(incVel, devVel, used);
tPos = stepTimes(tau, profPos, sum(coefPos, 1));
tVel = stepTimes(tau, profVel, sum(coefVel, 1));

% accelerometer bias: the sum of the increments over the window against the change of biasEst
incBias = NaN(1, 3);  dBias = NaN(1, 3);
kk = find(used);
if isfield(run, 'biasEst') && ~isempty(run.biasEst) && size(run.biasEst, 1) > 1 && ~isempty(kk)
    B = run.biasEst(:, 1:3);
    if size(B, 1) == K                                   % per epoch: the value after the update of epoch k
        k0 = max(kk(1) - 1, 1);
        dBias   = B(kk(end), :) - B(k0, :);
        incBias = nansumRows(incAcc(k0+1:kk(end), :));
    elseif size(B, 1) == numel(t)                        % per nav sample: before the first and after the last update
        i0 = find(t <= tKf(kk(1)) - TB2, 1, 'last');
        i1 = find(t >= tKf(kk(end)) + HALF, 1, 'first');
        if ~isempty(i0) && ~isempty(i1)
            dBias   = B(i1, :) - B(i0, :);
            incBias = nansumRows(incAcc(kk(1):kk(end), :));
        end
    end
end

% ---- printout ----
tUsed = tKf(used);
fprintf('feedbackCheck: %d updates between %.1f and %.1f s (%d epochs in the window, %d without states or samples), nav samples every %.3f s\n', ...
    nUsed, min([tUsed; Inf]), max([tUsed; -Inf]), nnz(inWin), nnz(inWin) - nUsed, dtNav);
if dtNav > 0.05
    fprintf('   WARNING: nav samples every %.3f s, the profile needs a log of at least 20 Hz\n', dtNav);
end
fprintf('   position: the change of navPos that navVel does not explain, summed over tKf - %.2f to tKf + %.2f s\n', HALF, HALF);
fprintf('   velocity: navVel minus its line %.1f to %.1f s before the epoch, averaged over tKf + %.2f to %.2f s\n', TB1, TB2, TA1, TA2);
fprintf('                                           N           E           D\n');
printBlock('position increment x+ - x- (1:3)', 'change of navPos across the epoch', incPos(used, :), chgPos(used, :), coefPos, rPos, 'm');
printBlock('velocity increment x+ - x- (4:6)', 'change of navVel across the epoch', incVel(used, :), chgVel(used, :), coefVel, rVel, 'm/s');
printProfile(tau, profPos, profVel, tPos, tVel);
if all(isfinite(dBias))
    ratio = dBias ./ incBias;
    ratio(abs(incBias) < 1e-9) = NaN;
    fprintf('   accelerometer bias (%d:%d) over the window\n', iAcc(1), iAcc(end));
    fprintf('      sum of the increments        %+11.3f %+11.3f %+11.3f   mg\n', incBias / MG);
    fprintf('      change of biasEst            %+11.3f %+11.3f %+11.3f   mg\n', dBias / MG);
    fprintf('      ratio change / sum           %+11.2f %+11.2f %+11.2f       (-1 fed back in full, estimate minus truth; 0 never)\n', ratio);
end
fprintf('   reading: a coefficient of -1 at k is the increment applied in full at its own epoch (error state estimate minus\n');
fprintf('   truth), +1 the opposite convention, 0 an increment that never reaches the logged solution; a coefficient at k-1 or\n');
fprintf('   k+1 instead means the records and the epochs are offset by one\n');

% ---- verdicts ----
verdict = {};
rmsPos = sqrt(mean(incPos(used, :).^2, 1));
rmsVel = sqrt(mean(incVel(used, :).^2, 1));
if nUsed < 10
    verdict{end+1} = sprintf('only %d usable updates: no reading', nUsed);
else
    verdict = [verdict, axisVerdicts('position', coefPos, rmsPos, 1e-4)];
    verdict = [verdict, axisVerdicts('velocity', coefVel, rmsVel, 1e-5)];
    verdict = [verdict, stepVerdicts('position', tPos, coefPos, rmsPos, 1e-4)];
    verdict = [verdict, stepVerdicts('velocity', tVel, coefVel, rmsVel, 1e-5)];
end
if all(isfinite(dBias)) && any(abs(incBias) > 1e-5)
    big = abs(incBias) > 1e-5;
    rb = dBias(big) ./ incBias(big);
    if all(abs(rb) < 0.3)
        verdict{end+1} = 'the accelerometer bias increments do not reach biasEst over the window';
    elseif any(abs(rb) < 0.3)
        verdict{end+1} = sprintf('the accelerometer bias increments reach biasEst on some axes only (ratios %s)', sprintf('%+.2f ', rb));
    end
end
for i = 1:numel(verdict)
    fprintf('   VERDICT: %s\n', verdict{i});
end

res = struct('t', tKf, 'incPos', incPos, 'incVel', incVel, 'chgPos', chgPos, 'chgVel', chgVel, 'used', used, ...
    'coefPos', coefPos, 'coefVel', coefVel, 'rPos', rPos, 'rVel', rVel, 'tau', tau, 'profPos', profPos, 'profVel', profVel, ...
    't10', [tPos(1, :); tVel(1, :)], 't50', [tPos(2, :); tVel(2, :)], 't90', [tPos(3, :); tVel(3, :)], ...
    'incBias', incBias, 'dBias', dBias, 'nUsed', nUsed);
res.verdict = verdict;
end

%------------------------------------------------------------------------
function printBlock(nameInc, nameChg, inc, chg, coef, r, unit)
mInc = mean(inc, 1);
same = NaN(1, 3);
for a = 1:3
    nz = inc(:, a) ~= 0;
    if any(nz), same(a) = 100 * mean(sign(inc(nz, a)) == sign(mInc(a))); end
end
fprintf('   %s\n', nameInc);
fprintf('      mean                         %+11.4f %+11.4f %+11.4f   %s\n', mInc, unit);
fprintf('      median                       %+11.4f %+11.4f %+11.4f   %s\n', median(inc, 1), unit);
fprintf('      with the sign of the mean    %10.0f%% %10.0f%% %10.0f%%   of the epochs\n', same);
fprintf('   %s\n', nameChg);
fprintf('      mean                         %+11.4f %+11.4f %+11.4f   %s\n', mean(chg, 1), unit);
fprintf('      on the increment of record k-1 %+9.3f %+11.3f %+11.3f\n', coef(1, :));
fprintf('                             k   %+11.3f %+11.3f %+11.3f\n', coef(2, :));
fprintf('                             k+1 %+11.3f %+11.3f %+11.3f\n', coef(3, :));
fprintf('      correlation with record k    %+11.2f %+11.2f %+11.2f\n', r);
end

function printProfile(tau, profPos, profVel, tPos, tVel)
show = [-0.2 -0.1 -0.05 -0.02 -0.01 0 0.01 0.02 0.05 0.1 0.2 0.25];
j = zeros(size(show));
for i = 1:numel(show), [~, j(i)] = min(abs(tau - show(i))); end
fprintf('   step profile: share of the increment in the logged solution at tKf + tau (regression across the epochs)\n');
fprintf('      tau (s)      %s\n', sprintf('%6.2f', show));
ax = 'NED';
for a = 1:3
    fprintf('      position %s   %s\n', ax(a), sprintf('%6.2f', -profPos(j, a)));
end
for a = 1:3
    fprintf('      velocity %s   %s\n', ax(a), sprintf('%6.2f', -profVel(j, a)));
end
fprintf('      (sign flipped: 1.00 is the increment subtracted in full)\n');
fprintf('      step at 10 / 50 / 90 %% of its size, position   %s\n', stepText(tPos));
fprintf('                                     velocity   %s\n', stepText(tVel));
end

function [s] = stepText(tS)
ax = 'NED';
s = '';
for a = 1:3
    if all(isfinite(tS(:, a)))
        s = [s, sprintf('%s %+.2f %+.2f %+.2f s   ', ax(a), tS(1, a), tS(2, a), tS(3, a))];
    else
        s = [s, sprintf('%s    -     -     -     ', ax(a))];
    end
end
end

function [coef, r] = lagFit(inc, chg, used)
% chg(k) = coef(1) inc(k-1) + coef(2) inc(k) + coef(3) inc(k+1), per axis, over the epochs whose three
% records are usable; r the plain correlation of chg(k) with inc(k)
K = size(inc, 1);
coef = NaN(3, 3);  r = NaN(1, 3);
k = find(used);  k = k(k >= 2 & k <= K - 1);
k = k(used(k - 1) & used(k + 1));
for a = 1:3
    y = chg(k, a);
    X = [inc(k - 1, a), inc(k, a), inc(k + 1, a)];
    ok = all(isfinite([X, y]), 2);
    if nnz(ok) > 10 && all(sum(X(ok, :).^2, 1) > 0)
        coef(:, a) = X(ok, :) \ y(ok);
        dx = X(ok, 2) - mean(X(ok, 2));  dy = y(ok) - mean(y(ok));
        den = sqrt(sum(dx.^2) * sum(dy.^2));
        if den > 0, r(a) = sum(dx .* dy) / den; end
    end
end
end

function [prof] = profileOf(inc, cum, used)
% slope of cum(k, j, a) on inc(k, a) across the epochs, per tau and axis
J = size(cum, 2);
prof = NaN(J, 3);
for a = 1:3
    x = inc(used, a);
    sxx = sum(x.^2);
    if sxx > 0
        for j = 1:J
            c = cum(used, j, a);
            ok = isfinite(c);
            prof(j, a) = sum(x(ok) .* c(ok)) / sum(x(ok).^2);
        end
    end
end
end

function [tS] = stepTimes(tau, prof, total)
% tau at which the profile first reaches 10, 50 and 90 % of the total coefficient, per axis
tS = NaN(3, 3);
for a = 1:3
    if isfinite(total(a)) && abs(total(a)) > 0.1
        p = prof(:, a) / total(a);
        lev = [0.1 0.5 0.9];
        for i = 1:3
            j = find(p >= lev(i), 1, 'first');
            if ~isempty(j), tS(i, a) = tau(j); end
        end
    end
end
end

function [v] = axisVerdicts(what, coef, rmsInc, tiny)
ax = 'NED';
v = {};
total = sum(coef, 1);
live = rmsInc > tiny & isfinite(total);
for a = 1:3
    if ~live(a)
        continue;
    elseif abs(total(a)) < 0.3
        v{end+1} = sprintf('the %s %s increment does not reach the logged solution (coefficient %+.2f): the mechanisation keeps its %s %s through the update', ...
            ax(a), what, total(a), ax(a), what);
    elseif abs(coef(2, a)) < 0.7 * abs(total(a)) && (abs(coef(1, a)) >= 0.5 * abs(total(a)) || abs(coef(3, a)) >= 0.5 * abs(total(a)))
        if abs(coef(1, a)) > abs(coef(3, a)), lag = 'k-1: the records are one epoch ahead of the nav log'; else, lag = 'k+1: the records are one epoch behind the nav log'; end
        v{end+1} = sprintf('the %s %s change matches the increment of record %s (coefficients %+.2f %+.2f %+.2f)', ...
            ax(a), what, lag, coef(:, a));
    elseif abs(coef(2, a)) < 0.7 * abs(total(a))
        v{end+1} = sprintf('the %s %s increment is applied over more than one epoch (coefficients %+.2f %+.2f %+.2f at k-1, k, k+1)', ...
            ax(a), what, coef(:, a));
    elseif abs(total(a)) < 0.7 || abs(total(a)) > 1.3
        v{end+1} = sprintf('the %s %s increment is applied in part (coefficient %+.2f): a feedback gain below one', ax(a), what, total(a));
    end
end
full = live & abs(total) >= 0.7 & abs(total) <= 1.3;
if nnz(full) >= 2 && any(sign(total(full)) ~= sign(total(find(full, 1))))
    v{end+1} = sprintf('the %s coefficients differ in sign across the axes (%s): one axis is applied with the opposite sign', what, sprintf('%+.2f ', total));
end
end

function [v] = stepVerdicts(what, tS, coef, rmsInc, tiny)
ax = 'NED';
v = {};
total = sum(coef, 1);
for a = 1:3
    if rmsInc(a) > tiny && all(isfinite(tS(:, a))) && abs(total(a)) >= 0.3 && tS(3, a) - tS(1, a) > 0.05
        v{end+1} = sprintf('the %s %s step is spread over %.2f s (10 to 90 %% between tKf %+.2f and %+.2f s): a slewed feedback, not a jump', ...
            ax(a), what, tS(3, a) - tS(1, a), tS(1, a), tS(3, a));
    end
end
end

function [e] = nedFrom(p)
% [lat lon h] rows to NED metres from the first row (WGS84 radii at each row)
L  = p(:, 1);  h = p(:, 3);
a  = 6378137;  e2 = 6.69437999014e-3;
RN = a * (1 - e2) ./ (1 - e2 * sin(L).^2).^1.5;
RE = a ./ sqrt(1 - e2 * sin(L).^2);
e = [(p(:, 1) - p(1, 1)) .* (RN + h), ...
     (p(:, 2) - p(1, 2)) .* (RE + h) .* cos(L), ...
    -(p(:, 3) - p(1, 3))];
end

function [s] = nansumRows(X)
X(~isfinite(X)) = 0;
s = sum(X, 1);
end
