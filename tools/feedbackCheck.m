function [res] = feedbackCheck(run, spfMeasLog, win)
%FEEDBACKCHECK  Does the KF update reach the mechanisation? The jump of the logged nav solution
%   across each 2 Hz update against the increment the filter computed, per NED axis.
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
%   Per epoch with at least one row, the nav samples 0.1 to 0.3 s before the epoch and 0.1 to
%   0.3 s after it are each fitted with a line in time; the two lines evaluated at the epoch
%   differ by the jump the update put into the logged solution (the 0.1 s left out either side
%   covers the latency of the update inside the host step; the lines remove the velocity and
%   the acceleration of the vehicle). The jump is regressed on the increment postState minus
%   priorState per axis, position (1:3, lat lon h converted to NED metres) and velocity (4:6):
%   a slope of -1 says the increment is applied in full with the error state as estimate minus
%   truth (+1 for the opposite convention), a slope of 0 that it never reaches the logged
%   solution. The position validates the method on every axis: its slopes must read -1 before
%   the velocity slopes mean anything. With biasEst the sum of the accelerometer bias
%   increments over the window is set against the change of biasEst, the same -1 / 0 reading
%   for the bias feedback.
%
%   Verdict lines name an axis whose velocity or position increment does not reach the logged
%   solution (slope within 0.3 of zero while the increments are not zero), one applied with the
%   sign opposite to the other axes, one applied in part, and a log without a jump on any axis
%   (the update is applied outside the 0.1 to 0.3 s windows, or tKf is not the update instant).
%
%   res: per-epoch t, incPos, incVel, jumpPos, jumpVel [K x 3] (NaN where not computed), used
%   [K x 1], the per-axis slopes slopePos, slopeVel with their standard errors sePos, seVel and
%   correlations rPos, rVel, the bias totals incBias and dBias [1 x 3] in m/s^2 (NaN without
%   biasEst), and the verdict lines in res.verdict (cell).
%   MATLAB desktop tool, not codegen. Runs under Octave.

if iscell(spfMeasLog), spfMeasLog = [spfMeasLog{:}]; end
if nargin < 3 || isempty(win), win = [-Inf Inf]; end
K   = numel(spfMeasLog);
t   = run.t(:);
tKf = run.tKf(:);
if numel(tKf) ~= K
    error('feedbackCheck: tKf has %d epochs, spfMeasLog %d records', numel(tKf), K);
end
T1 = 0.1;  T2 = 0.3;                                     % fit windows either side of the epoch (s)
MG = 9.80665e-3;                                         % m/s^2 per mg
dtNav  = median(diff(t));
navNed = nedFrom(run.navPos);                            % NED metres from the first sample
navVel = run.navVel;
iAcc = 10:12;
if isfield(run, 'pGyroFirst') && ~isempty(run.pGyroFirst) && run.pGyroFirst, iAcc = 13:15; end

incPos = NaN(K, 3);  incVel = NaN(K, 3);  incAcc = NaN(K, 3);
jumpPos = NaN(K, 3);  jumpVel = NaN(K, 3);
used = false(K, 1);  inWin = false(K, 1);
for k = 1:K
    s = spfMeasLog(k);
    dx = s.postState(:) - s.priorState(:);
    if all(isfinite(dx)) && numel(dx) >= max(iAcc)
        incAcc(k, :) = dx(iAcc)';                        % every finite epoch: the bias total needs them all
    end
    inWin(k) = tKf(k) >= win(1) && tKf(k) <= win(2);
    if ~inWin(k) || double(s.numMeas) < 1 || ~all(isfinite(dx(1:6))), continue; end
    iB = t >= tKf(k) - T2 & t <= tKf(k) - T1;
    iA = t >= tKf(k) + T1 & t <= tKf(k) + T2;
    if nnz(iB) < 4 || nnz(iA) < 4, continue; end
    incPos(k, :) = dx(1:3)';
    incVel(k, :) = dx(4:6)';
    jumpPos(k, :) = lineAt(t(iA) - tKf(k), navNed(iA, :)) - lineAt(t(iB) - tKf(k), navNed(iB, :));
    jumpVel(k, :) = lineAt(t(iA) - tKf(k), navVel(iA, :)) - lineAt(t(iB) - tKf(k), navVel(iB, :));
    used(k) = true;
end
nUsed = nnz(used);

[slopeVel, seVel, rVel] = slopeOf(incVel(used, :), jumpVel(used, :));
[slopePos, sePos, rPos] = slopeOf(incPos(used, :), jumpPos(used, :));

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
        i0 = find(t <= tKf(kk(1)) - T1, 1, 'last');
        i1 = find(t >= tKf(kk(end)) + T2, 1, 'first');
        if ~isempty(i0) && ~isempty(i1)
            dBias   = B(i1, :) - B(i0, :);
            incBias = nansumRows(incAcc(kk(1):kk(end), :));
        end
    end
end

% ---- printout ----
ax = 'NED';
tUsed = tKf(used);
fprintf('feedbackCheck: %d updates between %.1f and %.1f s (%d epochs in the window, %d without rows or samples), nav samples every %.3f s\n', ...
    nUsed, min([tUsed; Inf]), max([tUsed; -Inf]), nnz(inWin), nnz(inWin) - nUsed, dtNav);
fprintf('   line fits %.1f to %.1f s before and after each epoch; jump = line after minus line before, both at the epoch\n', T1, T2);
if dtNav > 0.05
    fprintf('   WARNING: nav samples every %.3f s, the fits need a log of at least 20 Hz\n', dtNav);
end
fprintf('                                           N           E           D\n');
printBlock('velocity increment x+ - x- (4:6)', 'navVel jump across the update', incVel(used, :), jumpVel(used, :), slopeVel, seVel, rVel, 'm/s');
printBlock('position increment x+ - x- (1:3)', 'navPos jump across the update (NED)', incPos(used, :), jumpPos(used, :), slopePos, sePos, rPos, 'm');
if all(isfinite(dBias))
    ratio = dBias ./ incBias;
    ratio(abs(incBias) < 1e-9) = NaN;
    fprintf('   accelerometer bias (%d:%d) over the window\n', iAcc(1), iAcc(end));
    fprintf('      sum of the increments        %+11.3f %+11.3f %+11.3f   mg\n', incBias / MG);
    fprintf('      change of biasEst            %+11.3f %+11.3f %+11.3f   mg\n', dBias / MG);
    fprintf('      ratio change / sum           %+11.2f %+11.2f %+11.2f       (-1 fed back in full, estimate minus truth; 0 never)\n', ratio);
end
fprintf('   reading: slope -1 the increment is applied in full (error state estimate minus truth), +1 the opposite\n');
fprintf('   convention, 0 the increment never reaches the logged solution\n');

% ---- verdicts ----
verdict = {};
rmsVel = sqrt(mean(incVel(used, :).^2, 1));
rmsPos = sqrt(mean(incPos(used, :).^2, 1));
if nUsed < 10
    verdict{end+1} = sprintf('only %d usable updates: no reading', nUsed);
elseif all(abs([slopeVel, slopePos]) < 0.3 | ~isfinite([slopeVel, slopePos]))
    verdict{end+1} = 'no jump on any axis within 0.3 s of tKf: the update is applied outside the fit windows, or tKf is not the update instant';
else
    verdict = [verdict, axisVerdicts('velocity', slopeVel, rmsVel, 1e-5)];
    verdict = [verdict, axisVerdicts('position', slopePos, rmsPos, 1e-4)];
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

res = struct('t', tKf, 'incPos', incPos, 'incVel', incVel, 'jumpPos', jumpPos, 'jumpVel', jumpVel, 'used', used, ...
    'slopePos', slopePos, 'sePos', sePos, 'rPos', rPos, 'slopeVel', slopeVel, 'seVel', seVel, 'rVel', rVel, ...
    'incBias', incBias, 'dBias', dBias, 'nUsed', nUsed);
res.verdict = verdict;
end

%------------------------------------------------------------------------
function printBlock(nameInc, nameJump, inc, jump, slope, se, r, unit)
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
fprintf('   %s\n', nameJump);
fprintf('      mean                         %+11.4f %+11.4f %+11.4f   %s\n', mean(jump, 1), unit);
fprintf('      slope on the increment       %+11.3f %+11.3f %+11.3f\n', slope);
fprintf('      standard error of the slope  %11.3f %11.3f %11.3f\n', se);
fprintf('      correlation                  %+11.2f %+11.2f %+11.2f\n', r);
end

function [v] = axisVerdicts(what, slope, rmsInc, tiny)
ax = 'NED';
v = {};
live = rmsInc > tiny & isfinite(slope);
full = live & abs(slope) >= 0.7 & abs(slope) <= 1.3;
for a = 1:3
    if ~live(a)
        continue;
    elseif abs(slope(a)) < 0.3
        v{end+1} = sprintf('the %s %s increment does not reach the logged solution (slope %+.2f): the mechanisation keeps its %s %s through the update', ...
            ax(a), what, slope(a), ax(a), what);
    elseif ~full(a)
        v{end+1} = sprintf('the %s %s increment is applied in part (slope %+.2f): a feedback gain below one, or part of the update outside the fit windows', ...
            ax(a), what, slope(a));
    end
end
if nnz(full) >= 2 && any(sign(slope(full)) ~= sign(slope(find(full, 1))))
    v{end+1} = sprintf('the %s slopes differ in sign across the axes (%s): one axis is applied with the opposite sign', what, sprintf('%+.2f ', slope));
end
end

function [slope, se, r] = slopeOf(x, y)
% y = slope x per column, through the origin: the slope, its standard error, the correlation
slope = NaN(1, 3);  se = NaN(1, 3);  r = NaN(1, 3);
for a = 1:3
    xi = x(:, a);  yi = y(:, a);
    ok = isfinite(xi) & isfinite(yi);
    xi = xi(ok);  yi = yi(ok);
    sxx = sum(xi.^2);
    if sxx > 0 && numel(xi) > 2
        slope(a) = sum(xi .* yi) / sxx;
        resid = yi - slope(a) * xi;
        se(a) = sqrt(sum(resid.^2) / (numel(xi) - 1) / sxx);
        dx = xi - mean(xi);  dy = yi - mean(yi);
        den = sqrt(sum(dx.^2) * sum(dy.^2));
        if den > 0, r(a) = sum(dx .* dy) / den; end
    end
end
end

function [a] = lineAt(tau, X)
% intercept at tau = 0 of the least-squares line through the rows of X
A = [ones(numel(tau), 1), tau(:)];
c = A \ X;
a = c(1, :);
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
