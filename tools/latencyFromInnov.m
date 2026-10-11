function [dlt, r2] = latencyFromInnov(innov, los, vel)
%LATENCYFROMINNOV  Measurement latency per epoch, from the pseudorange innovations.
%   [dlt, r2] = latencyFromInnov(innov, los, vel)
%   A pseudorange measured dlt seconds before the state it is compared with
%   carries an innovation of dlt * (u . v) per satellite (u: unit vector from
%   the receiver to the satellite, v: receiver velocity), on top of a common
%   clock term. Per epoch the function solves innov = c + dlt * (u . v) over
%   the satellites and returns dlt (s) and r2, the fraction of the innovation
%   variance that the term explains.
%     innov [K x M]      innovations z - h(x-) BEFORE the update, NaN where a
%                        satellite is not used (m)
%     los   [K x M x 3]  unit vectors receiver -> satellite
%     vel   [K x 3]      receiver (nav) velocity at the instant of the state
%                        used in h(x), same frame as los (NED or ECEF)
%   Prints the median and the 16..84% spread of dlt over the epochs with at
%   least 5 satellites. A clean pairing gives |dlt| of a few ms; 0.5 s means
%   the observables are one 2 Hz epoch older than the state they are compared
%   with. Split dlt by epoch type (mod(tow, 1) == 0 against 0.5) to see
%   whether only one of the two epochs is paired wrongly. If the host defines
%   the innovation as h(x) - z, or u towards the receiver, dlt changes sign;
%   the magnitude stands. Doppler innovations work the same way with the
%   acceleration in place of the velocity.
%   MATLAB desktop tool, not codegen. Runs under Octave.

K   = size(innov, 1);
dlt = NaN(K, 1);
r2  = NaN(K, 1);
for k = 1:K
    u  = reshape(los(k, :, :), size(los, 2), 3);
    uv = u * vel(k, :)';
    y  = innov(k, :)';
    ok = isfinite(y) & isfinite(uv);
    if nnz(ok) < 5, continue; end
    A  = [ones(nnz(ok), 1), uv(ok)];
    x  = A \ y(ok);
    res = y(ok) - A * x;
    dlt(k) = x(2);
    yc = y(ok) - mean(y(ok));
    if sum(yc.^2) > 0, r2(k) = 1 - sum(res.^2) / sum(yc.^2); end
end
ok = isfinite(dlt);
if ~any(ok)
    fprintf('latencyFromInnov: no epoch with 5 or more satellites\n');
    return
end
d  = sort(dlt(ok));
n  = numel(d);
q  = @(f) d(max(1, min(n, round(f * n))));
fprintf('latency from the pseudorange innovations: median %+.3f s (16..84%%: %+.3f .. %+.3f s) over %d epochs;\n', ...
    q(0.5), q(0.16), q(0.84), n);
fprintf('it explains a median %.0f%% of the innovation variance per epoch\n', 100 * median(r2(ok)));
end
