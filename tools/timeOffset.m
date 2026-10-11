function [tau, stdBest, std0] = timeOffset(tA, posA, tB, posB, tauGrid)
%TIMEOFFSET  Time shift between two position records of the same vehicle.
%   tau = timeOffset(tA, posA, tB, posB)
%   finds tau such that A(t) matches B(t + tau): the shift of record B's time
%   base that minimises the horizontal difference between the two records.
%   The cost is a standard deviation (mean removed), so a constant offset
%   between the records, a lever arm or a GNSS bias, does not bias tau.
%     tA, tB    [N x 1], [M x 1]  time of each record (s), same clock units
%     posA, posB [N x 3], [M x 3] [lat lon h] in rad, rad, m
%     tauGrid   optional, default -1.5:0.01:1.5 s
%   Returns tau and the horizontal difference std at tau and at zero shift.
%
%   Filter-free use: A = the receiver's own position output (same time tags
%   as its observables), B = the reference. A tau of -0.45 s then means the
%   receiver record at tag t shows the vehicle where the reference puts it
%   0.45 s earlier, i.e. the GNSS tags run 0.45 s late against the reference,
%   and the host must apply each GNSS epoch 0.45 s earlier in IMU time.
%   The same call with A = the nav solution reproduces the position lag of
%   runCompare. MATLAB desktop tool, not codegen. Runs under Octave.

if nargin < 5 || isempty(tauGrid), tauGrid = -1.5:0.01:1.5; end
tA = tA(:);  tB = tB(:);
c  = NaN(size(tauGrid));
for i = 1:numel(tauGrid)
    c(i) = horizStd(posA, interpRows(tB, posB, tA + tauGrid(i)));
end
[stdBest, iBest] = min(c);
tau  = tauGrid(iBest);
[~, i0] = min(abs(tauGrid));
std0 = c(i0);
fprintf('A(t) = B(t %+.2f s): horizontal difference std %.2f m (%.2f m at zero shift)\n', tau, stdBest, std0);
if tau < 0
    fprintf('record A trails record B by %.2f s\n', -tau);
elseif tau > 0
    fprintf('record A leads record B by %.2f s\n', tau);
else
    fprintf('the two records are on the same time base\n');
end
if abs(tau) >= tauGrid(end) - 0.011
    fprintf('WARNING: at the edge of the scan, widen tauGrid\n');
end
end

function [v] = horizStd(p, q)
% pooled std of the horizontal NED difference p - q over the finite rows, mean removed
L  = p(:, 1);  h = p(:, 3);
a  = 6378137;  e2 = 6.69437999014e-3;
RN = a * (1 - e2) ./ (1 - e2 * sin(L).^2).^1.5;
RE = a ./ sqrt(1 - e2 * sin(L).^2);
e  = [(p(:, 1) - q(:, 1)) .* (RN + h), (p(:, 2) - q(:, 2)) .* (RE + h) .* cos(L)];
ok = all(isfinite(e), 2);
e  = e(ok, :);
e  = e - repmat(mean(e, 1), size(e, 1), 1);
v  = sqrt(mean(sum(e.^2, 2)));
end

function [y] = interpRows(tt, Y, tq)
y = NaN(numel(tq), size(Y, 2));
for k = 1:size(Y, 2)
    y(:, k) = interp1(tt(:), Y(:, k), tq(:), 'linear', NaN);
end
end
