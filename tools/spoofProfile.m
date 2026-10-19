function [prof] = spoofProfile(spfMeasLog, tKf, sysId)
%SPOOFPROFILE  What a spoofed run did to the measurements, from the saved spfMeas records.
%   prof = spoofProfile(spfMeasLog, tKf)          per epoch the range rows split into two groups where a gap in
%                                                 their innovations is wider than the group spread; the group
%                                                 farther from zero is the offset group. Prints onset, end, size
%                                                 and rows of the offset, and its course in time
%   prof = spoofProfile(spfMeasLog, tKf, sysId)   sysId [K x MAX_MEAS]: the system of each row as the host built
%                                                 the rows (0 or NaN = none). Adds the mean innovation per system
%                                                 and names the system that carries the offset
%
%   The innovation is measured against the host's prediction. In NOMINAL mode the filter follows a slow spoof and
%   the innovation shows only its rate times the filter lag, so the offset read here is a lower bound; in COAST
%   (host INS-only) the prediction is the coast and the innovation is the spoof offset plus the INS drift, so the
%   offset is read whole. A step or a large push shows at full size at its first epoch either way.
%   Range rows are those with a unit line of sight in the position columns of H and nothing in the velocity
%   columns (the pressure row, a single position entry, is left out).
%
%   Result: t, nRange, nOff (rows in the offset group), offset (its mean innovation minus the clean group's),
%   meanOff, meanClean, spreadClean, present (offset group of 2 rows or more, 10 m or more away), onset and tEnd
%   (first and last epoch of 3 consecutive present epochs), sizeOnset, sizeMax, sysIds and perSys [K x nSys]
%   (mean range innovation per system), sysOff (the system holding most of the offset rows at onset).

if iscell(spfMeasLog), spfMeasLog = [spfMeasLog{:}]; end
K = numel(spfMeasLog);
if nargin < 2 || isempty(tKf), tKf = 0.5 * (0:K-1)'; end
tKf = tKf(:);
if nargin < 3, sysId = []; end
haveSys = ~isempty(sysId) && size(sysId, 1) == K;
if ~isempty(sysId) && ~haveSys
    fprintf('   sysId has %d rows, the log %d epochs: not used\n', size(sysId, 1), K);
end
MIN_GAP = 10;                                                % m: a gap narrower than this is no second group

nRange = zeros(K, 1);  nOff = zeros(K, 1);  meanOff = NaN(K, 1);  meanClean = NaN(K, 1);  spreadClean = NaN(K, 1);
sysIds = [];
if haveSys, sysIds = unique(sysId(isfinite(sysId) & sysId > 0))'; end
perSys = NaN(K, numel(sysIds));  nPerSys = zeros(K, numel(sysIds));  offSys = zeros(K, numel(sysIds));
for k = 1:K
    s = spfMeasLog(k);
    m = double(s.numMeas);
    if m < 2, continue; end
    H = s.obsMatrix(1:m, :);
    y = s.innovation(1:m);
    normP = sqrt(sum(H(:, 1:3).^2, 2));
    normV = sqrt(sum(H(:, 4:6).^2, 2));
    isPr  = normP > 0.5 & normV < 0.1 & sum(abs(H(:, 1:3)) > 1e-9, 2) >= 2;
    idx = find(isPr & isfinite(y));
    nRange(k) = numel(idx);
    if haveSys
        for j = 1:numel(sysIds)
            r = idx(sysId(k, idx) == sysIds(j));
            nPerSys(k, j) = numel(r);
            if ~isempty(r), perSys(k, j) = mean(y(r)); end
        end
    end
    if numel(idx) < 3, continue; end
    yR = y(idx);
    [ys, order] = sort(yR);
    gaps = diff(ys);
    [gMax, iG] = max(gaps);
    lo = ys(1:iG);  hi = ys(iG+1:end);
    madLo = 1.4826 * median(abs(lo - median(lo)));  madHi = 1.4826 * median(abs(hi - median(hi)));   % spread inside each group
    if gMax > max(MIN_GAP, 5 * max(madLo, madHi))
        if abs(mean(hi)) >= abs(mean(lo)), offG = hi;  cleanG = lo;  offRows = idx(order(iG+1:end));
        else, offG = lo;  cleanG = hi;  offRows = idx(order(1:iG)); end
        nOff(k) = numel(offG);  meanOff(k) = mean(offG);  meanClean(k) = mean(cleanG);
        if numel(cleanG) > 1, spreadClean(k) = std(cleanG); else, spreadClean(k) = 0; end
        if haveSys
            for j = 1:numel(sysIds), offSys(k, j) = nnz(sysId(k, offRows) == sysIds(j)); end
        end
    else
        meanClean(k) = mean(yR);  spreadClean(k) = std(yR);
    end
end
offset  = meanOff - meanClean;
present = nOff >= 2 & abs(offset) >= MIN_GAP;
run3 = conv(double(present), ones(3, 1), 'valid') == 3;      % 3 consecutive present epochs
onset = NaN;  tEnd = NaN;  kOn = NaN;  kEnd = NaN;
if any(run3)
    kOn  = find(run3, 1);  onset = tKf(kOn);
    kEnd = find(run3, 1, 'last') + 2;  tEnd = tKf(kEnd);
end
sizeOnset = NaN;  sizeMax = NaN;  sysOff = NaN;
if isfinite(kOn)
    sizeOnset = mean(offset(kOn:min(kOn + 2, K)));
    seg = kOn:kEnd;
    [~, iMax] = max(abs(offset(seg)));  sizeMax = offset(seg(iMax));
    if haveSys
        cnt = sum(offSys(kOn:min(kOn + 2, K), :), 1);
        [~, j] = max(cnt);  if cnt(j) > 0, sysOff = sysIds(j); end
    end
end

fprintf('spoofProfile: %d epochs; range rows per epoch median %d; an offset group is 2 rows or more, %d m or more from the rest\n', ...
    K, round(median(nRange(nRange > 0))), MIN_GAP);
if ~isfinite(onset)
    fprintf('   no offset group on 3 consecutive epochs: nothing in the range innovations that a group of rows carries apart\n');
else
    fprintf('   onset        %10.2f s (epoch %d)\n', onset, kOn);
    fprintf('   end          %10.2f s (epoch %d), %.1f s long\n', tEnd, kEnd, tEnd - onset);
    fprintf('   size at onset %+9.1f m (mean of the first 3 epochs), largest %+.1f m\n', sizeOnset, sizeMax);
    fprintf('   rows in the offset group at onset %d of %d range rows\n', nOff(kOn), nRange(kOn));
    if haveSys && isfinite(sysOff)
        fprintf('   the offset rows belong to system %d\n', sysOff);
    elseif haveSys
        fprintf('   the offset rows belong to no single system\n');
    end
    seg = kOn:kEnd;
    gapsIn = seg(~present(seg));
    if ~isempty(gapsIn)
        fprintf('   %d of the %d epochs between onset and end show no offset group (the rows were absent or the host left them out)\n', numel(gapsIn), numel(seg));
    end
    fprintf('   course of the offset: t [s]   offset [m]   offset rows / range rows   clean spread [m]');
    if haveSys, fprintf('   mean innovation per system [m]: %s', sprintf('%d ', sysIds)); end
    fprintf('\n');
    pick = unique(round(linspace(seg(1), seg(end), min(12, numel(seg)))));
    for k = pick
        fprintf('   %21.2f   %+10.1f   %3d / %3d                   %6.1f', tKf(k), offset(k), nOff(k), nRange(k), spreadClean(k));
        if haveSys, fprintf('   %s', sprintf('%+8.1f ', perSys(k, :))); end
        fprintf('\n');
    end
end
prof = struct('t', tKf, 'nRange', nRange, 'nOff', nOff, 'offset', offset, 'meanOff', meanOff, 'meanClean', meanClean, ...
    'spreadClean', spreadClean, 'present', present, 'onset', onset, 'tEnd', tEnd, 'sizeOnset', sizeOnset, 'sizeMax', sizeMax, ...
    'sysIds', sysIds, 'perSys', perSys, 'nPerSys', nPerSys, 'sysOff', sysOff);
end
