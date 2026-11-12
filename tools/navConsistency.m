function [res] = navConsistency(run, tKf, pAtEpochs, aidedWin)
%NAVCONSISTENCY  Does the host's P describe its own errors? Each state group against the sigma the filter carried.
%   res = navConsistency(run, tKf, pAtEpochs)
%   res = navConsistency(run, tKf, pAtEpochs, aidedWin)
%   run         the runCompare struct: t, navPos [lat lon h] (rad, rad, m), tTrue, truePos; and, for the groups that
%               need them, navVel and trueVel (NED, m/s, trueVel on tTrue), navAtt and trueAtt ([roll pitch yaw], rad),
%               biasEst [N x 6] or [K x 6] (accelerometer then gyro, the value the mechanisation subtracts) and
%               biasTrue [1 x 6]; pGyroFirst true when the host's P orders the gyro biases 10:12 before the
%               accelerometer biases 13:15 (default true)
%   tKf         [K x 1] the 2 Hz epochs (s)
%   pAtEpochs   the host's P per epoch as coastDriftCheck takes it: a cell of n x n (or of diagonals), [n x n x K], or
%               the diagonal [K x n]
%   aidedWin    [W x 2] windows [t1 t2] to read (s); [] or omitted: the whole run
%   Per group it prints the rms of the error, its mean over the window, the rms of the sigma the filter carried,
%   their ratio, the same ratio with the mean error removed, and the share of epochs beyond 3 sigma: a consistent
%   P reads a ratio near 1 and a share near 0.3 %; a ratio of 3 is a P three times wider than the error, the
%   pessimistic case; a ratio of 0.3 is an overconfident P, and its share beyond 3 sigma is large. A constant
%   offset between the nav and the truth (a reference-point or datum offset, a GNSS latency at steady speed) is
%   not something a P can cover: the mean-removed ratio reads the P against the part of the error that varies.
%   Groups: position N, E, D; velocity N, E, D; tilt (the horizontal pair, rotation-invariant); heading;
%   accelerometer bias; gyro bias. A group whose truth is not in the run is skipped. The tilt and bias groups
%   compare magnitudes (the mean removed per component), so the host's sign conventions do not matter.
%   Result: res.group (names), res.rmsErr, res.meanErr (the norm of the mean per component), res.rmsErrC (mean
%   removed), res.rmsSig, res.ratio, res.ratioC (mean removed), res.beyond3 (per group), res.K (epochs read).
%   MATLAB desktop tool, not codegen. Runs under Octave.
if nargin < 4, aidedWin = []; end
tKf = tKf(:);  K = numel(tKf);
gyroFirst = ~isfield(run, 'pGyroFirst') || isempty(run.pGyroFirst) || run.pGyroFirst;
if gyroFirst, gyroIdx = 10:12;  accIdx = 13:15; else, accIdx = 10:12;  gyroIdx = 13:15; end
[sig, n] = sigmaAt(pAtEpochs, K);
use = true(K, 1);
if ~isempty(aidedWin)
    use = false(K, 1);
    for w = 1:size(aidedWin, 1), use = use | (tKf >= aidedWin(w, 1) & tKf <= aidedWin(w, 2)); end
end
G = struct('names', {{}}, 'rmsErr', [], 'meanErr', [], 'rmsErrC', [], 'rmsSig', [], 'ratio', [], 'ratioC', [], 'beyond3', []);
% position
e = navErrAt(run, tKf);
ax = {'N', 'E', 'D'};
for a = 1:3
    G = add(G, ['position ' ax{a}], e(:, a), sig(:, a), use);
end
% velocity
if isfield(run, 'navVel') && ~isempty(run.navVel) && isfield(run, 'trueVel') && ~isempty(run.trueVel)
    nv = interpRows(run.t, run.navVel, tKf);  tv = interpRows(run.tTrue, run.trueVel, tKf);
    for a = 1:3
        G = add(G, ['velocity ' ax{a}], nv(:, a) - tv(:, a), sig(:, 3 + a), use);
    end
end
% tilt and heading
if isfield(run, 'navAtt') && ~isempty(run.navAtt) && isfield(run, 'trueAtt') && ~isempty(run.trueAtt)
    na = interpRows(run.t, run.navAtt, tKf);  ta = interpRows(run.tTrue, run.trueAtt, tKf);
    d = na - ta;  d = atan2(sin(d), cos(d));
    tiltSig = sqrt(sig(:, 7) .^ 2 + sig(:, 8) .^ 2);
    G = add(G, 'tilt (roll, pitch)', d(:, 1:2), tiltSig, use);
    G = add(G, 'heading', d(:, 3), sig(:, 9), use);
end
% biases
if isfield(run, 'biasEst') && ~isempty(run.biasEst) && isfield(run, 'biasTrue') && ~isempty(run.biasTrue)
    b = run.biasEst;
    if size(b, 1) == numel(run.t), b = interpRows(run.t, b, tKf);
    elseif size(b, 1) == 1, b = repmat(b, K, 1);
    elseif size(b, 1) ~= K, b = interpRows(linspace(run.t(1), run.t(end), size(b, 1))', b, tKf); end
    bt = run.biasTrue;
    if size(bt, 1) == 1, bt = repmat(bt, K, 1); elseif size(bt, 1) ~= K, bt = interpRows(linspace(run.t(1), run.t(end), size(bt, 1))', bt, tKf); end
    db = b - bt;
    accSig = sqrt(sum(sig(:, accIdx) .^ 2, 2));  gyrSig = sqrt(sum(sig(:, gyroIdx) .^ 2, 2));
    G = add(G, 'accelerometer bias', db(:, 1:3), accSig, use);
    G = add(G, 'gyro bias', db(:, 4:6), gyrSig, use);
end
fprintf('navConsistency: %d epochs read of %d, %d states in P\n', nnz(use), K, n);
fprintf('   %-22s %12s %12s %12s %8s %10s %10s\n', 'group', 'rms error', 'mean error', 'rms sigma', 'ratio', 'mean out', 'beyond 3s');
for g = 1:numel(G.names)
    fprintf('   %-22s %12.4g %12.4g %12.4g %8.2f %10.2f %9.2f%%\n', G.names{g}, G.rmsErr(g), G.meanErr(g), G.rmsSig(g), ...
        G.ratio(g), G.ratioC(g), 100 * G.beyond3(g));
end
fprintf('   (ratio = sigma / error: near 1 consistent, 3 a P three times the error, 0.3 overconfident; mean out: the same with the\n');
fprintf('   mean error over the window removed, a constant offset being no business of P; a consistent Gaussian P has 0.3 %% beyond\n');
fprintf('   3 sigma; position in m, velocity m/s, angles rad, biases in the units of biasEst; the mean of a multi-axis group is its norm)\n');
res = struct('group', {G.names}, 'rmsErr', G.rmsErr, 'meanErr', G.meanErr, 'rmsErrC', G.rmsErrC, 'rmsSig', G.rmsSig, ...
    'ratio', G.ratio, 'ratioC', G.ratioC, 'beyond3', G.beyond3, 'K', nnz(use));
end

function [G] = add(G, name, err, sg, use)
% err [K x c]: the error's components (one for a signed axis, two or three for a magnitude group); sg [K x 1]
ok = use & all(isfinite(err), 2) & isfinite(sg) & sg > 0;
if nnz(ok) < 2, return; end
mag = sqrt(sum(err(ok, :) .^ 2, 2));
mu = mean(err(ok, :), 1);
magC = sqrt(sum((err(ok, :) - mu) .^ 2, 2));
re = sqrt(mean(mag .^ 2));  reC = sqrt(mean(magC .^ 2));  rs = sqrt(mean(sg(ok) .^ 2));
G.names{end + 1} = name;  G.rmsErr(end + 1) = re;  G.meanErr(end + 1) = norm(mu);  G.rmsErrC(end + 1) = reC;
G.rmsSig(end + 1) = rs;  G.ratio(end + 1) = rs / re;  G.ratioC(end + 1) = rs / reC;
G.beyond3(end + 1) = mean(mag > 3 * sg(ok));
end

function [sig, n] = sigmaAt(p, K)
% the sigma of every state per epoch from the forms the P may take
if iscell(p)
    n = round(sqrt(numel(p{1})));  if numel(p{1}) ~= n * n, n = numel(p{1}); end
    sig = NaN(K, n);
    for k = 1:min(K, numel(p))
        x = p{k};
        if numel(x) == n * n, sig(k, :) = sqrt(max(diag(reshape(x, n, n)), 0))'; else, sig(k, :) = sqrt(max(x(:), 0))'; end
    end
elseif ndims(p) == 3
    n = size(p, 1);  sig = NaN(K, n);
    for k = 1:min(K, size(p, 3)), sig(k, :) = sqrt(max(diag(p(:, :, k)), 0))'; end
else
    n = size(p, 2);  sig = NaN(K, n);
    m = min(K, size(p, 1));  sig(1:m, :) = sqrt(max(p(1:m, :), 0));
end
end

function [y] = interpRows(t, x, tq)
y = NaN(numel(tq), size(x, 2));
for c = 1:size(x, 2), y(:, c) = interp1(t(:), x(:, c), tq(:), 'linear', NaN); end
end

function [e] = navErrAt(run, t)
% nav minus truth in NED (m) at the times t
a = 6378137;  e2 = 6.69437999014e-3;
nav = interpRows(run.t, run.navPos, t);  tr = interpRows(run.tTrue, run.truePos, t);
lat = tr(:, 1);  s2 = sin(lat) .^ 2;
Rm = a * (1 - e2) ./ (1 - e2 * s2) .^ 1.5;  Rn = a ./ sqrt(1 - e2 * s2);
e = [(nav(:, 1) - tr(:, 1)) .* (Rm + tr(:, 3)), (nav(:, 2) - tr(:, 2)) .* (Rn + tr(:, 3)) .* cos(lat), -(nav(:, 3) - tr(:, 3))];
end
