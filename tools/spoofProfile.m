function [prof] = spoofProfile(spfMeasLog, tKf, sysId, loggedClean)
%SPOOFPROFILE  What a spoofed run did to the measurements, from the saved spfMeas records.
%   prof = spoofProfile(spfMeasLog, tKf)          per epoch the range rows split into two groups where a gap in
%                                                 their innovations is wider than the group spread; the group
%                                                 farther from zero is the offset group. Prints onset, end, size
%                                                 and rows of the offset, and its course in time
%   prof = spoofProfile(spfMeasLog, tKf, sysId)   sysId [K x MAX_MEAS]: the system of each row as the host built
%                                                 the rows (0 or NaN = none). Adds the mean innovation per system
%                                                 and names the system that carries the offset
%   prof = spoofProfile(spfMeasLog, tKf, loggedData)   the host's pairing log instead (prLog / satPosLog per epoch,
%                                                 as pairingCheck takes it): the used range slots in slot order give
%                                                 the system of the record's range rows in the order the host built
%                                                 them; the tool counts the epochs where the two counts differ and
%                                                 warns when the host's row order cannot be the slot order
%   prof = spoofProfile(spfMeasLog, tKf, loggedData, loggedClean)   the clean twin's pairing log as well: the
%                                                 spoofed ranges minus the clean ones, per slot and epoch, are the
%                                                 injection itself, read without any filter in the way. Valid when
%                                                 the spoofed file is the clean file plus an injection, epoch for
%                                                 epoch; the tool prints the injected offset per system and its
%                                                 onset, end and size next to the innovation-based reading
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
%   (mean range innovation per system), sysOff (the system holding most of the offset rows at onset), and with the
%   clean twin: injSys, injOnset, injEnd, injSize and injMax per system, nSlotChanged.

if iscell(spfMeasLog), spfMeasLog = [spfMeasLog{:}]; end
K = numel(spfMeasLog);
if nargin < 2 || isempty(tKf), tKf = 0.5 * (0:K-1)'; end
tKf = tKf(:);
if nargin < 3, sysId = []; end
if nargin < 4, loggedClean = []; end
fromLog = isstruct(sysId) || iscell(sysId);
prSlot = [];  sysSlot = [];
if fromLog
    [prSlot, ~, sysSlot] = pairingLogUnpack(sysId);
    haveSys = size(prSlot, 1) == K;
    if ~haveSys, fprintf('   the pairing log has %d epochs, the spfMeas records %d: not used\n', size(prSlot, 1), K); end
else
    haveSys = ~isempty(sysId) && size(sysId, 1) == K;
    if ~isempty(sysId) && ~haveSys
        fprintf('   sysId has %d rows, the log %d epochs: not used\n', size(sysId, 1), K);
    end
end
nSysMismatch = 0;
MIN_GAP = 10;                                                % m: a gap narrower than this is no second group

nRange = zeros(K, 1);  nOff = zeros(K, 1);  meanOff = NaN(K, 1);  meanClean = NaN(K, 1);  spreadClean = NaN(K, 1);
sysIds = [];
if haveSys && fromLog
    sysIds = unique(sysSlot(isfinite(sysSlot) & sysSlot > 0))';
elseif haveSys
    sysIds = unique(sysId(isfinite(sysId) & sysId > 0))';
end
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
    rowSys = NaN(1, numel(idx));                              % the system of each range row of this epoch
    if haveSys && fromLog
        usedSys = sysSlot(k, isfinite(prSlot(k, :)));         % the used range slots in slot order
        if numel(usedSys) ~= numel(idx), nSysMismatch = nSysMismatch + 1; end
        nn = min(numel(usedSys), numel(idx));
        rowSys(1:nn) = usedSys(1:nn);
    elseif haveSys
        rowSys = sysId(k, idx);
    end
    if haveSys
        for j = 1:numel(sysIds)
            r = idx(rowSys == sysIds(j));
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
            offSysRow = rowSys(ismember(idx, offRows));
            for j = 1:numel(sysIds), offSys(k, j) = nnz(offSysRow == sysIds(j)); end
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
if haveSys && fromLog
    nWith = nnz(nRange > 0);
    fprintf('   row systems from the pairing log: the used range slots in slot order; their count matched the record''s range rows on %d of %d epochs\n', ...
        nWith - nSysMismatch, nWith);
    if nSysMismatch > 0.1 * max(nWith, 1)
        fprintf('   WARNING: the host does not build its range rows in slot order (or the pairing log and the records differ): the per-system\n');
        fprintf('            numbers and the named system are unreliable; give sysId per row instead\n');
    end
end
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
%% the injection itself: spoofed ranges minus the clean twin's, per slot
injSys = [];  injOnset = NaN;  injEnd = NaN;  injSize = NaN(1, 0);  injMax = NaN(1, 0);  nSlotChanged = NaN;
if ~isempty(loggedClean)
    if ~fromLog
        fprintf('   the clean twin is used only with the spoofed pairing log as the third argument\n');
    else
        [prC, ~, sysC] = pairingLogUnpack(loggedClean);
        if ~isequal(size(prC), size(prSlot))
            fprintf('   the clean twin''s pairing log is %d x %d, the spoofed one %d x %d: not compared\n', size(prC), size(prSlot));
        else
            both = isfinite(prC) & isfinite(prSlot);
            dPr  = prSlot - prC;  dPr(~both) = NaN;
            nSlotChanged = nnz(xor(isfinite(prC), isfinite(prSlot)));
            injSys = unique(sysSlot(both & sysSlot > 0))';
            injMean = NaN(K, numel(injSys));  injRows = zeros(K, numel(injSys));
            for j = 1:numel(injSys)
                m = both & sysSlot == injSys(j);
                for k = 1:K
                    if any(m(k, :)), injMean(k, j) = mean(dPr(k, m(k, :)));  injRows(k, j) = nnz(m(k, :)); end
                end
            end
            injPresent = any(abs(injMean) >= 1, 2);                 % an injection of 1 m or more on some system
            run3i = conv(double(injPresent), ones(3, 1), 'valid') == 3;
            if any(run3i)
                kI = find(run3i, 1);  kE = find(run3i, 1, 'last') + 2;
                injOnset = tKf(kI);  injEnd = tKf(kE);
                injSize = mean(injMean(kI:min(kI + 2, K), :), 1);
                [~, iM] = max(abs(injMean(kI:kE, :)), [], 1);
                injMax = NaN(1, numel(injSys));
                for j = 1:numel(injSys), injMax(j) = injMean(kI - 1 + iM(j), j); end
                fprintf('   injection (spoofed ranges minus the clean twin''s, per used slot): onset %.2f s (epoch %d), end %.2f s (epoch %d), %.1f s long\n', ...
                    injOnset, kI, injEnd, kE, injEnd - injOnset);
                fprintf('   per system: %s\n', strjoin(arrayfun(@(j) sprintf('system %d: %+.1f m at onset, largest %+.1f m, %d slots', ...
                    injSys(j), injSize(j), injMax(j), round(median(injRows(kI:kE, j)))), 1:numel(injSys), 'UniformOutput', false), '; '));
                if nSlotChanged > 0
                    fprintf('   %d slot-epochs are used in one file and not the other (satellites the spoof added or removed)\n', nSlotChanged);
                end
                pick = unique(round(linspace(kI, kE, min(12, kE - kI + 1))));
                fprintf('   course of the injection: t [s]   mean offset per system [m]: %s\n', sprintf('%d ', injSys));
                for k = pick
                    fprintf('   %25.2f   %s\n', tKf(k), sprintf('%+10.1f ', injMean(k, :)));
                end
            else
                fprintf('   injection (spoofed ranges minus the clean twin''s): none of 1 m or more on 3 consecutive epochs; the two files hold the same ranges\n');
            end
        end
    end
end
prof = struct('t', tKf, 'nRange', nRange, 'nOff', nOff, 'offset', offset, 'meanOff', meanOff, 'meanClean', meanClean, ...
    'spreadClean', spreadClean, 'present', present, 'onset', onset, 'tEnd', tEnd, 'sizeOnset', sizeOnset, 'sizeMax', sizeMax, ...
    'sysIds', sysIds, 'perSys', perSys, 'nPerSys', nPerSys, 'sysOff', sysOff, 'nSysMismatch', nSysMismatch, ...
    'injSys', injSys, 'injOnset', injOnset, 'injEnd', injEnd, 'injSize', injSize, 'injMax', injMax, 'nSlotChanged', nSlotChanged);
end
