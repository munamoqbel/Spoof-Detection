function [cost, shifts] = pairingCheck(tEpoch, pr, satPos, tRef, posRef, shifts, sysId)
%PAIRINGCHECK  Which instant of the reference the pseudoranges of each epoch describe.
%   [cost, shifts] = pairingCheck(tEpoch, pr, satPos, tRef, posRef)
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
%   pulse's own epoch. MATLAB desktop tool, not codegen. Runs under Octave.

if nargin < 6 || isempty(shifts), shifts = -1.5:0.1:0.5; end
if nargin < 7 || isempty(sysId), sysId = ones(size(pr)); end
K = size(pr, 1);  M = size(pr, 2);
refEcef = llaToEcef(posRef);
tEpoch  = tEpoch(:);
cost = NaN(size(shifts));
for i = 1:numel(shifts)
    p = NaN(K, 3);
    for c = 1:3, p(:, c) = interp1(tRef(:), refEcef(:, c), tEpoch + shifts(i), 'linear', NaN); end
    s = NaN(K, 1);
    for k = 1:K
        if ~all(isfinite(p(k, :))), continue; end
        sat = reshape(satPos(k, :, :), M, 3);
        rho = sqrt(sum((sat - repmat(p(k, :), M, 1)).^2, 2));
        r   = pr(k, :)' - rho;
        ok  = isfinite(r);
        if nnz(ok) < 5, continue; end
        g = sysId(k, :)';
        for id = unique(g(ok))'
            m = ok & g == id;
            r(m) = r(m) - mean(r(m));
        end
        s(k) = sqrt(mean(r(ok).^2));
    end
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
