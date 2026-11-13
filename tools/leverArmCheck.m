function [res] = leverArmCheck(spfMeasLog, tKf, run, tWin)
%LEVERARMCHECK  Do the rate-row innovations grow with the yaw rate? The lever-arm compensation read from the records.
%   res = leverArmCheck(spfMeasLog, tKf, run)         run: the runCompare struct, t and navAtt ([roll pitch yaw], rad)
%   res = leverArmCheck(spfMeasLog, tKf, run, tWin)   tWin [t1 t2] in s on tKf; default the whole log
%   An antenna r metres from the IMU moves at omega x r in a turn, and the velocity rows carry that motion. The host
%   compensates it with its lever arm; what is left in the rate-row innovations in the turns is the error of that
%   arm, in proportion to the yaw rate. Per epoch the tool projects the rate rows onto the horizontal velocity along
%   and across the heading (the row's velocity columns against the heading of navAtt) and takes the mean innovation
%   of each, then regresses the two on the yaw rate (the derivative of the unwrapped yaw, smoothed over 1 s) over the
%   epochs that turn faster than 2 deg/s: the slope is the effective arm error in metres, across track for an
%   along-body arm and along track for a cross-body arm; the intercept is the straight-leg residual. The straight
%   legs (|yaw rate| below 0.5 deg/s) give the rms the turns are compared with, and the turning epochs the rms with
%   the regression removed. A slope beyond its own two sigma and 0.1 m is an arm error of that size, whether or not
%   the turning rms shows it (a 0.3 m arm at 15 deg/s is 0.08 m/s, under the row noise, and only the regression
%   sees it); a turning rms above the straight one without a slope is not the arm; neither clears the arm.
%   Result: t, yawRate, innAlong, innAcross (per epoch, NaN without rate rows), slope [1 x 2], slopeSig, intercept,
%   rmsStraight [1 x 2], rmsTurn, rmsTurnResid, nTurn, nStraight.
%   MATLAB desktop tool, not codegen. Runs under Octave.
if nargin < 4, tWin = []; end
R2D = 180 / pi;
tKf = tKf(:);  K = numel(tKf);
if numel(spfMeasLog) < K, K = numel(spfMeasLog);  tKf = tKf(1:K); end
% the heading and the yaw rate on the epochs
yaw = unwrap(run.navAtt(:, 3));
yawK = interp1(run.t(:), yaw, tKf, 'linear', NaN);
dtN = median(diff(run.t));  w = max(1, round(1 / dtN));
rate = [diff(yaw); 0] / dtN;
rate = filter(ones(w, 1) / w, 1, rate);                         % smoothed over 1 s
rateK = interp1(run.t(:), rate, tKf, 'linear', NaN);
% the innovations along and across the heading
innAlong = NaN(K, 1);  innAcross = NaN(K, 1);
for k = 1:K
    s = spfMeasLog(k);  m = double(s.numMeas);
    if m < 1 || ~isfinite(yawK(k)), continue; end
    H = s.obsMatrix(1:m, :);  y = s.innovation(1:m);
    normV = sqrt(sum(H(:, 4:6) .^ 2, 2));  normP = sqrt(sum(H(:, 1:3) .^ 2, 2));
    rows = find(normV > 0.5 & normP < 0.1);
    if isempty(rows), continue; end
    uAlong = [cos(yawK(k)); sin(yawK(k)); 0];  uAcross = [-sin(yawK(k)); cos(yawK(k)); 0];
    cA = H(rows, 4:6) * uAlong;  cX = H(rows, 4:6) * uAcross;
    wa = abs(cA) > 0.3;  wx = abs(cX) > 0.3;
    if any(wa), innAlong(k) = mean(y(rows(wa)) ./ cA(wa)); end
    if any(wx), innAcross(k) = mean(y(rows(wx)) ./ cX(wx)); end
end
use = isfinite(rateK) & (isfinite(innAlong) | isfinite(innAcross));
if ~isempty(tWin), use = use & tKf >= tWin(1) & tKf <= tWin(2); end
turn = use & abs(rateK) > 2 / R2D;  straight = use & abs(rateK) < 0.5 / R2D;
res = struct('t', tKf, 'yawRate', rateK, 'innAlong', innAlong, 'innAcross', innAcross, 'slope', NaN(1, 2), 'slopeSig', NaN(1, 2), ...
    'intercept', NaN(1, 2), 'rmsStraight', NaN(1, 2), 'rmsTurn', NaN(1, 2), 'rmsTurnResid', NaN(1, 2), 'nTurn', nnz(turn), 'nStraight', nnz(straight));
fprintf('leverArmCheck: %d epochs, %d turning faster than 2 deg/s, %d straight below 0.5 deg/s; yaw rate up to %.1f deg/s\n', ...
    nnz(use), nnz(turn), nnz(straight), max(abs(rateK(use))) * R2D);
names = {'along track', 'across track'};
inn = [innAlong, innAcross];
fprintf('   %-13s %10s %10s %12s %12s %12s %10s\n', 'projection', 'straight', 'turn rms', 'slope [m]', 'sigma [m]', 'intercept', 'turn resid');
fprintf('   %-13s %10s %10s %12s %12s %12s %10s\n', '', 'rms [m/s]', '[m/s]', '', '', '[m/s]', 'rms [m/s]');
for a = 1:2
    v = inn(:, a);
    s = straight & isfinite(v);  tr = turn & isfinite(v);
    if nnz(s) >= 5, res.rmsStraight(a) = sqrt(mean(v(s) .^ 2)); end
    if nnz(tr) >= 10
        res.rmsTurn(a) = sqrt(mean(v(tr) .^ 2));
        A = [rateK(tr), ones(nnz(tr), 1)];
        c = A \ v(tr);
        r = v(tr) - A * c;
        cov = (r' * r / (nnz(tr) - 2)) * inv(A' * A);
        res.slope(a) = c(1);  res.slopeSig(a) = sqrt(cov(1, 1));  res.intercept(a) = c(2);
        res.rmsTurnResid(a) = sqrt(mean(r .^ 2));
    end
    fprintf('   %-13s %10.3f %10.3f %12.3f %12.3f %12.3f %10.3f\n', names{a}, res.rmsStraight(a), res.rmsTurn(a), res.slope(a), ...
        res.slopeSig(a), res.intercept(a), res.rmsTurnResid(a));
end
fprintf('   (slope: the innovation per rad/s of yaw rate, the effective arm error in metres: an along-body arm error shows across\n');
fprintf('   track, a cross-body one along track; intercept: the straight-leg residual the turns share)\n');
sig = abs(res.slope) > 2 * res.slopeSig & abs(res.slope) > 0.1;
grow = res.rmsTurn > 1.3 * res.rmsStraight;
if any(sig)
    a = find(sig);
    fprintf('   VERDICT: the rate-row innovations grow with the yaw rate, %s: an effective lever-arm error of %s m\n', ...
        strjoin(names(a), ' and '), strjoin(arrayfun(@(x) sprintf('%.2f', x), res.slope(a), 'UniformOutput', false), ' and '));
elseif any(grow)
    fprintf('   VERDICT: the innovations are larger in the turns but not in proportion to the yaw rate: not the arm; a heading that\n');
    fprintf('            changes in the turns or a latency on the rate rows\n');
else
    fprintf('   VERDICT: the rate-row innovations in the turns are no larger than on the straight legs: the lever arm is compensated\n');
end
end
