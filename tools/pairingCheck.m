function [cost, shifts] = pairingCheck(tEpoch, pr, satPos, tRef, posRef, shifts, sysId)
%PAIRINGCHECK  Which instant of the reference the pseudoranges of each epoch describe.
%   [cost, shifts] = pairingCheck(tEpoch, pr, satPos, tRef, posRef)
%   [cost, shifts] = pairingCheck(tEpoch, loggedData, [], tRef, posRef, shifts)
%   [cost, shifts] = pairingCheck(tEpoch, loggedData, 'order', tRef, posRef, shifts)
%     loggedData: the host's per-epoch log (prLog / satPosLog / sysIdLog, one block per
%     constellation), unpacked by pairingLogUnpack; satPos and sysId then come from it.
%     'order' pairs the used ranges with the live position rows by order instead of by
%     slot (see pairingLogUnpack). When the residual at the best shift is far above code
%     noise a slot diagnostic says, on a few epochs, which position slot each range fits
%     with one receiver clock per system: ranges fitting other slots mean the log is
%     compacted on one side, no fit at all means the positions are not ECEF at the epoch
%   For each candidate shift s, the reference position at tEpoch + s predicts
%   the ranges to the satellites. The measured pseudoranges minus those
%   ranges, with the per-epoch mean (the receiver clock) removed, leave a
%   residual whose spread over the satellites is smallest at the instant the
%   observables were really measured. Prints the median residual std per
%   shift; cost(i) is that median for shifts(i). Filter-free.
%     tEpoch  [K x 1]      host time the host assigns to each epoch (the pulse
%                          whose interpolated state the call compares with)
%     pr      [K x M]      pseudoranges as used in h(x): satellite clock and
%                          atmosphere corrected the way the host does, NaN
%                          where a satellite is not used (m)
%     satPos  [K x M x 3]  satellite ECEF positions as used in h(x), i.e. at
%                          transmit time in the receive-time frame (m)
%     tRef    [N x 1]      reference time on the host counter (s)
%     posRef  [N x 3]      reference [lat lon h] (rad, rad, m)
%     shifts  optional, default -1.5:0.1:0.5 (s)
%     sysId   optional [K x M] constellation id per satellite: the clock mean
%                          is removed per system (inter-system bias)
%   A minimum at s = -0.5 s means the observables compared with the pulse
%   state are one 2 Hz epoch older than the host believes; at 0 they are the
%   pulse's own epoch. The receiver clock removed per system at the best shift
%   is printed, with a warning when the residual there is far above code noise
%   (satellite clock or atmosphere still in the ranges, or the wrong satellite
%   positions) and a note when it is far below it at zero shift (a range built
%   from the nav solution, not a measurement).
%   MATLAB desktop tool, not codegen. Runs under Octave.

unpacked = isstruct(pr) || (iscell(pr) && ~isempty(pr) && isstruct(pr{1}));
if unpacked
    pairing = 'slot';
    if ischar(satPos) && ~isempty(satPos), pairing = satPos; end
    [pr, satPos, sysId] = pairingLogUnpack(pr, pairing);
end
if nargin < 6 || isempty(shifts), shifts = -1.5:0.1:0.5; end
if ~unpacked && (nargin < 7 || isempty(sysId)), sysId = ones(size(pr)); end
K = size(pr, 1);
if numel(tEpoch) ~= K
    error('pairingCheck:epochs', 'tEpoch has %d entries, the ranges %d epochs', numel(tEpoch), K);
end
refEcef = llaToEcef(posRef);
tEpoch  = tEpoch(:);
cost = NaN(size(shifts));
for i = 1:numel(shifts)
    s = residualAtShift(tEpoch, pr, satPos, tRef, refEcef, sysId, shifts(i));
    cost(i) = median(s(isfinite(s)));
end
[cBest, iBest] = min(cost);
fprintf('residual std over the satellites, median over the epochs, per shift of the reference instant:\n');
fprintf('   shift %+5.2f s: %6.2f m\n', [shifts(:)'; cost(:)']);
fprintf('smallest at %+.2f s (%.2f m): the observables describe the vehicle %.2f s %s the instant the host assigns to them\n', ...
    shifts(iBest), cBest, abs(shifts(iBest)), ifelse(shifts(iBest) < 0, 'before', 'after'));
if iBest == 1 || iBest == numel(shifts)
    fprintf('WARNING: at the edge of the scan, widen shifts\n');
end
[~, clk, ids] = residualAtShift(tEpoch, pr, satPos, tRef, refEcef, sysId, shifts(iBest));
if ~isempty(ids)
    fprintf('receiver clock removed per system at the best shift, median over the epochs:');
    for j = 1:numel(ids), fprintf(' system %d %+.0f m', ids(j), median(clk(isfinite(clk(:, j)), j))); end
    fprintf('\n');
end
if cBest > 100
    fprintf('WARNING: a residual of %.0f m at the best shift is far above code noise: the satellite clock or the atmosphere are still in\n', cBest);
    fprintf('         the pseudoranges, or the satellite positions are not the ones of h(x) (transmit time, earth rotation), or the\n');
    fprintf('         slots of the ranges and of the positions do not match\n');
    slotDiagnostic(tEpoch, pr, satPos, tRef, refEcef, sysId);
elseif abs(shifts(iBest)) < 0.015 && cBest < 0.3
    fprintf('NOTE: a minimum at zero shift with a residual of %.2f m is what a range built from the nav solution gives; a measured\n', cBest);
    fprintf('      pseudorange carries code noise of a metre or more. Check that the logged range is the measurement, not the prediction\n');
end
end

%------------------------------------------------------------------------
function slotDiagnostic(tEpoch, pr, satPos, tRef, refEcef, sysId)
% on a few epochs at zero shift: for each system, the receiver clock c and the position slot that each range
% fits with it, from the offsets D(i, j) = pr_i - |p_j - ref| (the correct pairing shares one c within TOL)
TOL = 50;
M  = size(pr, 2);
ep = find(sum(isfinite(pr), 2) >= 4);
if isempty(ep), return; end
ep = ep(unique(round(linspace(1, numel(ep), 5))));
ids = unique(sysId(isfinite(pr)))';
nFit = 0;  nTot = 0;  moved = 0;
fprintf('slot diagnostic at zero shift on %d epochs: the position slot each range fits with one receiver clock per system (within %d m)\n', numel(ep), TOL);
for k = ep'
    p0 = interp1(tRef(:), refEcef, tEpoch(k), 'linear', NaN);
    if ~all(isfinite(p0)), continue; end
    sat = reshape(satPos(k, :, :), M, 3);
    rho = sqrt(sum((sat - repmat(p0, M, 1)).^2, 2));
    for id = ids
        block = find(sysId(k, :) == id);
        ri = block(isfinite(pr(k, block)));                        % ranges of this system
        pj = block(rho(block)' > 1e6);                             % position rows that hold a satellite
        if numel(ri) < 4 || isempty(pj), continue; end
        D = repmat(pr(k, ri)', 1, numel(pj)) - repmat(rho(pj)', numel(ri), 1);   % [ranges x positions]
        best = 0;  jBest = ones(numel(ri), 1);  cBest = NaN;
        for c = D(1, :)
            [dmin, jmin] = min(abs(D - c), [], 2);
            n = nnz(dmin < TOL);
            if n > best, best = n;  jBest = jmin;  cBest = c; end
        end
        ok = abs(D(sub2ind(size(D), (1:numel(ri))', jBest)) - cBest) < TOL;
        rs = ri(ok) - block(1) + 1;  ps = pj(jBest(ok)) - block(1) + 1;        % slots within the system block
        fprintf('   epoch %d, system %d: %d of %d ranges fit a position (clock %+.0f m); range slots %s -> position slots %s\n', ...
            k, id, best, numel(ri), cBest, mat2str(rs), mat2str(ps));
        nFit = nFit + best;  nTot = nTot + numel(ri);  moved = moved + nnz(rs ~= ps);
    end
end
if nTot == 0, return; end
if nFit >= 0.8 * nTot && moved > 0
    fprintf('   the ranges fit positions in OTHER slots: one side of the log is compacted or ordered differently. If the used positions\n');
    fprintf('   are the live rows in order, rerun with pairingCheck(tEpoch, loggedData, ''order'', ...); otherwise log the positions in\n');
    fprintf('   the slots of the ranges\n');
elseif nFit >= 0.8 * nTot
    fprintf('   the slots pair; the residual comes from the values: satellite clock or atmosphere still in the ranges, or positions\n');
    fprintf('   not at the measurement epoch\n');
else
    fprintf('   no slot assignment fits most ranges: the positions are not ECEF at the measurement epoch (inertial frame, ephemeris\n');
    fprintf('   reference time, another unit), or the ranges are not pseudoranges. Compare satPosLog.position with the satellite\n');
    fprintf('   position the range-row builder uses for the same epoch\n');
end
end

function [s, clk, ids] = residualAtShift(tEpoch, pr, satPos, tRef, refEcef, sysId, shift)
% per epoch: the residual std over the satellites after the per-system mean (the receiver clock)
% is removed, with the reference at tEpoch + shift; clk [K x nSys] the means removed, ids the systems
K = size(pr, 1);  M = size(pr, 2);
ids = unique(sysId(isfinite(pr)))';
p = NaN(K, 3);
for c = 1:3, p(:, c) = interp1(tRef(:), refEcef(:, c), tEpoch + shift, 'linear', NaN); end
s   = NaN(K, 1);
clk = NaN(K, numel(ids));
for k = 1:K
    if ~all(isfinite(p(k, :))), continue; end
    sat = reshape(satPos(k, :, :), M, 3);
    rho = sqrt(sum((sat - repmat(p(k, :), M, 1)).^2, 2));
    r   = pr(k, :)' - rho;
    ok  = isfinite(r);
    if nnz(ok) < 5, continue; end
    g = sysId(k, :)';
    for j = 1:numel(ids)
        m = ok & g == ids(j);
        if any(m)
            clk(k, j) = mean(r(m));
            r(m) = r(m) - clk(k, j);
        end
    end
    s(k) = sqrt(mean(r(ok).^2));
end
end

function [s] = ifelse(c, a, b)
if c, s = a; else, s = b; end
end

function [e] = llaToEcef(p)
a = 6378137;  e2 = 6.69437999014e-3;
L = p(:, 1);  l = p(:, 2);  h = p(:, 3);
N = a ./ sqrt(1 - e2 * sin(L).^2);
e = [(N + h) .* cos(L) .* cos(l), (N + h) .* cos(L) .* sin(l), (N * (1 - e2) + h) .* sin(L)];
end
