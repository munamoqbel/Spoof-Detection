function [cost, shifts] = pairingCheck(tEpoch, pr, satPos, tRef, posRef, shifts, sysId)
%PAIRINGCHECK  Which instant of the reference the pseudoranges of each epoch describe.
%   [cost, shifts] = pairingCheck(tEpoch, pr, satPos, tRef, posRef)
%   [cost, shifts] = pairingCheck(tEpoch, loggedData, [], tRef, posRef, shifts)
%     loggedData: the host's per-epoch log (prLog / satPosLog / sysIdLog, one block per
%     constellation), unpacked by pairingLogUnpack; satPos and sysId then come from it
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
if unpacked, [pr, satPos, sysId] = pairingLogUnpack(pr); end
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
elseif abs(shifts(iBest)) < 0.015 && cBest < 0.3
    fprintf('NOTE: a minimum at zero shift with a residual of %.2f m is what a range built from the nav solution gives; a measured\n', cBest);
    fprintf('      pseudorange carries code noise of a metre or more. Check that the logged range is the measurement, not the prediction\n');
end
end

%------------------------------------------------------------------------
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
