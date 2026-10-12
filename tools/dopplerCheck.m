function [res] = dopplerCheck(spfMeasLog, tKf)
%DOPPLERCHECK  What the range-rate (Doppler) rows do in the host KF update, read from the
%   spfMeas records saved per 2 Hz epoch.
%   res = dopplerCheck(spfMeasLog)          spfMeasLog: struct array or cell of the per-epoch spfMeas
%   res = dopplerCheck(spfMeasLog, tKf)     tKf: epoch times for res.t (default 1:K)
%
%   spfMeas (STRUCT_SPF.setKfMeas / kfMeasFromUpdate) holds, per epoch, the innovation y, the
%   innovation covariance S = H P- H' + R, the rows H, the measurement noise R, numMeas, and the
%   state before and after the update. Nothing else is needed: H P- H' is S - R.
%
%   Rows are classified from H with the host's state layout (position 1:3, velocity 4:6, NED
%   error states in m and m/s): a row whose velocity partial has unit size is a range-rate row,
%   a row whose position partial has unit size and whose velocity partial is zero a range row.
%   A row with a single position entry (pressure height) is left out of the pairing below.
%
%   Printed per row type (medians over the epochs): rows per epoch, sqrt(R), sqrt(H P- H'),
%   innovation rms, the gain fraction H P- H' / S (what share of each innovation the filter
%   applies: 0 means the row is ignored), and the NIS y' S^-1 y / m over the rows of the type
%   (near 1 when R and P are honest). Then the velocity correction postState - priorState (4:6) actually applied per
%   epoch, and the sign/frame consistency of the two derivatives: every range-rate row is paired
%   with the range row of the same line of sight (closest position partial) and the dot product
%   of the two unit partials is +1 when the velocity partial of the range-rate row uses the
%   sign convention and frame of the position partial of the range row, -1 for the opposite
%   sign, short of 1 in size for another frame.
%
%   Verdict lines are printed when a figure is out of range:
%     gain fraction (range rate) < 0.1      R far above H P- H': the Doppler rows are switched off in effect
%     innovation rms < 0.3 sqrt(H P- H')    the measured range rate is built from the nav solution itself
%     NIS median > 3                        innovations larger than S allows: R too small, or biased rows (a wrong-sign H)
%     pair dot < -0.99                      the velocity partial has the opposite sign to the position partial
%     |pair dot| < 0.99                     the two partials are not in the same frame
%
%   res: per-epoch arrays t, nPr, nRr, sR, sHph, innRms, gain, nis (each [K x 2], column 1 range
%   rows, column 2 range-rate rows), dvRms (velocity correction norm), pairDot (mean over the
%   pairs of the epoch), and the summary medians in res.summary.
%   MATLAB desktop tool, not codegen. Runs under Octave.

if iscell(spfMeasLog), spfMeasLog = [spfMeasLog{:}]; end
K = numel(spfMeasLog);
if nargin < 2 || isempty(tKf), tKf = (1:K)'; end
t = tKf(:);

nPr = zeros(K, 1);  nRr = zeros(K, 1);
sR = NaN(K, 2);  sHph = NaN(K, 2);  innRms = NaN(K, 2);  gain = NaN(K, 2);  nis = NaN(K, 2);
dvRms = NaN(K, 1);  pairDot = NaN(K, 1);  nPair = zeros(K, 1);
for k = 1:K
    s = spfMeasLog(k);
    m = double(s.numMeas);
    if m == 0, continue; end
    H = s.obsMatrix(1:m, :);
    y = s.innovation(1:m);
    rD = diag(s.measNoiseCov(1:m, 1:m));
    Sm = s.innovationCov(1:m, 1:m);
    sD = diag(Sm);
    hph = max(sD - rD, 0);                                   % H P- H' per row
    normP = sqrt(sum(H(:, 1:3).^2, 2));
    normV = sqrt(sum(H(:, 4:6).^2, 2));
    isRr = normV > 0.5 & normV > normP;
    isPr = normP > 0.5 & normV < 0.1;
    nPr(k) = nnz(isPr);  nRr(k) = nnz(isRr);
    sel = {isPr, isRr};
    for c = 1:2
        i = sel{c};
        if ~any(i), continue; end
        sR(k, c)     = sqrt(median(rD(i)));
        sHph(k, c)   = sqrt(median(hph(i)));
        innRms(k, c) = sqrt(mean(y(i).^2));
        gain(k, c)   = median(hph(i) ./ max(sD(i), eps));
        nis(k, c)    = (y(i)' * (Sm(i, i) \ y(i))) / nnz(i);     % y' S^-1 y / m over the rows of the type: ~1 when R and P are honest
    end
    dx = s.postState(:) - s.priorState(:);
    dvRms(k) = norm(dx(4:6));
    % pair every range-rate row with the range row of the closest line of sight: a consistent host
    % gives exactly the same unit vector (dot +1), an opposite sign -1, another frame leaves the
    % best match short of 1 (a chance alignment with another satellite is never exact)
    los = isPr & sum(abs(H(:, 1:3)) > 1e-9, 2) >= 2;         % a satellite line of sight, not the pressure row
    if any(los) && any(isRr)
        uP = H(los, 1:3) ./ repmat(normP(los), 1, 3);
        uV = H(isRr, 4:6) ./ repmat(normV(isRr), 1, 3);
        d  = uV * uP';                                       % [nRr x nLos] dot products
        [~, jj] = max(abs(d), [], 2);
        dBest = d(sub2ind(size(d), (1:size(d, 1))', jj));
        pairDot(k) = mean(dBest);
        nPair(k) = numel(dBest);
    end
end

use = nRr > 0;
sm = struct('epochs', K, 'epochsRr', nnz(use), 'nPr', median(nPr(nPr > 0)), 'nRr', median(nRr(use)), ...
    'sR', nanmedian(sR), 'sHph', nanmedian(sHph), 'innRms', sqrt(nanmean(innRms.^2)), 'gain', nanmedian(gain), ...
    'nis', nanmedian(nis), 'dvRms', sqrt(nanmean(dvRms.^2)), 'pairDot', nanmedian(pairDot), 'nPair', median(nPair(nPair > 0)));
if isempty(sm.nPr), sm.nPr = 0; end
if isempty(sm.nRr), sm.nRr = 0; end
if isempty(sm.nPair), sm.nPair = 0; end

fprintf('dopplerCheck: %d epochs, %d with range rows, %d with range-rate rows\n', K, nnz(nPr > 0), nnz(use));
fprintf('                                   range rows     range-rate rows\n');
fprintf('   rows per epoch (median)      %8.0f        %8.0f\n', sm.nPr, sm.nRr);
fprintf('   sqrt(R)                      %8.3f m      %8.3f m/s\n', sm.sR);
fprintf('   sqrt(H P- H'')                %8.3f m      %8.3f m/s     (what the filter expects from its own state)\n', sm.sHph);
fprintf('   innovation rms               %8.3f m      %8.3f m/s\n', sm.innRms);
fprintf('   gain fraction H P- H'' / S    %8.3f        %8.3f         (share of the innovation applied; 0 = row ignored)\n', sm.gain);
fprintf('   NIS y'' S^-1 y / m (median)   %8.3f        %8.3f         (~1 when R and P are honest)\n', sm.nis);
fprintf('   velocity correction applied per epoch |x+ - x-|(4:6): %.3f m/s rms\n', sm.dvRms);
if isfinite(sm.pairDot)
    fprintf('   derivative consistency: range-rate velocity partial . range position partial = %+.3f (median, %.0f pairs per epoch)\n', sm.pairDot, sm.nPair);
    fprintf('   (+1 same sign convention and frame, -1 opposite sign, anything short of 1 in size another frame or mixed rows)\n');
end
if nnz(use) == 0
    fprintf('   VERDICT: no range-rate rows found in H (no row with a unit velocity partial): the Doppler is not in the update\n');
else
    if sm.gain(2) < 0.1
        fprintf('   VERDICT: the range-rate gain fraction is %.3f: R (%.3f m/s) is far above what the filter expects from its own state\n', sm.gain(2), sm.sR(2));
        fprintf('            (%.3f m/s); the Doppler rows are switched off in effect. Set R from the receiver Doppler noise, ~0.05-0.1 m/s\n', sm.sHph(2));
    end
    if sm.innRms(2) < 0.3 * sm.sHph(2)
        fprintf('   VERDICT: the range-rate innovations (%.3f m/s rms) are far below the velocity uncertainty the filter carries (%.3f m/s):\n', sm.innRms(2), sm.sHph(2));
        fprintf('            the measured range rate is built from the nav solution itself (self-referential), not from the receiver Doppler\n');
    end
    if sm.nis(2) > 3
        fprintf('   VERDICT: the range-rate NIS is %.1f: innovations larger than R + H P- H'' allows. R too small, or biased rows (a wrong-sign\n', sm.nis(2));
        fprintf('            or wrong-frame velocity partial, a lever arm or a clock-drift term missing from the predicted range rate)\n');
    end
    if isfinite(sm.pairDot) && sm.pairDot < -0.99
        fprintf('   VERDICT: the velocity partial of the range-rate rows has the OPPOSITE sign to the position partial of the range rows\n');
    elseif isfinite(sm.pairDot) && abs(sm.pairDot) < 0.99
        fprintf('   VERDICT: the velocity partial of the range-rate rows is not in the frame of the position partial of the range rows\n');
        fprintf('            (or the two row types do not describe the same satellites)\n');
    end
end
if nnz(nPr > 0) > 0 && sm.nis(1) > 3
    fprintf('   VERDICT: the range NIS is %.1f: the pseudorange innovations are larger than R + H P- H'' allows\n', sm.nis(1));
end

res = struct('t', t, 'nPr', nPr, 'nRr', nRr, 'sR', sR, 'sHph', sHph, 'innRms', innRms, 'gain', gain, 'nis', nis, ...
    'dvRms', dvRms, 'pairDot', pairDot, 'nPair', nPair, 'summary', sm);
end

%------------------------------------------------------------------------
function [m] = nanmedian(x)
m = NaN(1, size(x, 2));
for c = 1:size(x, 2)
    v = x(isfinite(x(:, c)), c);
    if ~isempty(v), m(c) = median(v); end
end
end

function [m] = nanmean(x)
m = NaN(1, size(x, 2));
for c = 1:size(x, 2)
    v = x(isfinite(x(:, c)), c);
    if ~isempty(v), m(c) = mean(v); end
end
end
