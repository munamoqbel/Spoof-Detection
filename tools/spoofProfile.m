function [prof] = spoofProfile(spfMeasLog, tKf, sysId, tRef, posRef, spfLog)
%SPOOFPROFILE  What a spoofed run did to the measurements, from the saved spfMeas records and the pairing log.
%   prof = spoofProfile(spfMeasLog, tKf)          per epoch the range rows split into two groups where a gap in
%                                                 their innovations is wider than the group spread; the group
%                                                 farther from zero is the offset group. Prints onset, end, size
%                                                 and rows of the offset, and its course in time
%   prof = spoofProfile(spfMeasLog, tKf, sysId)   sysId [K x MAX_MEAS]: the system of each row as the host built
%                                                 the rows (0 or NaN = none). Adds the mean innovation per system,
%                                                 names the system that carries the offset, prints the course of
%                                                 the mean innovation per system with the host's row set epoch by
%                                                 epoch over the stretch where some system is 10 m or more off,
%                                                 and fits the position and clock the off system's rows ask for
%   prof = spoofProfile(spfMeasLog, tKf, loggedData)   the host's pairing log instead (prLog / satPosLog per epoch,
%                                                 as pairingCheck takes it): the used range slots in slot order give
%                                                 the system of the record's range rows in the order the host built
%                                                 them; the tool counts the epochs where the two counts differ and
%                                                 warns when the host's row order cannot be the slot order
%   prof = spoofProfile(spfMeasLog, tKf, loggedData, tRef, posRef)   the reference position as well, [lat lon h]
%                                                 (rad, rad, m) on tRef, one row for a stationary run: the profile
%                                                 with no filter in the way. Each used range minus the geometric
%                                                 range from the reference leaves the receiver clock, the system's
%                                                 inter-system bias and the injection. The biases are run medians,
%                                                 the clock is read on the epochs where two systems agree within
%                                                 5 m and interpolated across the others, and what is left per
%                                                 system is the injection: onset, end, size, system and course,
%                                                 and per epoch the position and clock the off system's ranges
%                                                 describe (a least-squares fix on them), relative to the reference
%   prof = spoofProfile(..., spfLog)              the gate log from SPF_logAppend as well: the timeline around the
%                                                 spoof carries the gate mode under the row set, and the events
%
%   The innovation is measured against the host's prediction. In NOMINAL mode the filter follows a slow spoof and
%   the innovation shows only its rate times the filter lag, so the offset read here is a lower bound; in COAST
%   (host INS-only) the prediction is the coast and the innovation is the spoof offset plus the INS drift, so the
%   offset is read whole. A step or a large push shows at full size at its first epoch either way. The two-group
%   split needs rows of both kinds in the same record: an epoch whose rows are all spoofed shows no split, only
%   the per-system course and the reference-based profile see it. A spoof that moves the position rather than the
%   clock gives every satellite its own offset: the per-system mean then depends on which satellites are in the
%   record, and the fitted position and clock are the right summary.
%   Range rows are those with a unit line of sight in the position columns of H and nothing in the velocity
%   columns (the pressure row, a single position entry, is left out).
%
%   Result: t, nRange, nOff (rows in the offset group), offset (its mean innovation minus the clean group's),
%   meanOff, meanClean, spreadClean, present (offset group of 2 rows or more, 10 m or more away), onset and tEnd
%   (first and last epoch of 3 consecutive present epochs), sizeOnset, sizeMax, sysIds, perSys and nPerSys
%   [K x nSys] (mean range innovation per system and its rows), sysOff (the system holding most of the offset
%   rows at onset), nSysMismatch, offSysInn (the system off in the innovations), fitInn [K x 6] (the position NED,
%   clock, rms and rows the off system's innovations ask for); with the reference: inj [K x nSys] (the offset per
%   system against the reference, clock and biases removed), nSlots [K x nSys], spread [K x nSys] (within the
%   system), clock [K x 1], clockGap (distance in s to the nearest epoch the clock was read on), isb [1 x nSys],
%   injSys, injOnset, injEnd, injSize, injMax, nBadSlot, fitRef [K x 6] (position NED, clock, rms, rows of the
%   fix on the off system's ranges, relative to the reference), fitSign (1 when the two fits agree, -1 when they
%   agree with the sign reversed, 0 otherwise), nRecNoRows and nLogNoSlots (epochs of the timeline window where
%   the record has no range row but the log has used slots, and the reverse), stretchesInn and stretchesRef
%   ([kStart kEnd nOff] per stretch; off epochs closer than 20 epochs form one stretch, the tables and the
%   timeline follow the strongest).
%   MATLAB desktop tool, not codegen. Runs under Octave.

if iscell(spfMeasLog), spfMeasLog = [spfMeasLog{:}]; end
K = numel(spfMeasLog);
if nargin < 2 || isempty(tKf), tKf = 0.5 * (0:K-1)'; end
tKf = tKf(:);
if nargin < 3, sysId = []; end
if nargin < 4, tRef = []; end
if nargin < 5, posRef = []; end
if nargin < 6, spfLog = []; end
fromLog = isstruct(sysId) || iscell(sysId);
prSlot = [];  satSlot = [];  sysSlot = [];
nRejSlot = [];
if fromLog
    [prSlot, satSlot, sysSlot, nRejSlot] = pairingLogUnpack(sysId);
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
BAD_RES = 1e5;                                               % m: a range this far from its slot's geometry is not its satellite
STRETCH_GAP = 20;                                            % epochs: off epochs closer than this belong to one stretch
gapS = STRETCH_GAP * medianFinite(diff(tKf));

nRange = zeros(K, 1);  nOff = zeros(K, 1);  meanOff = NaN(K, 1);  meanClean = NaN(K, 1);  spreadClean = NaN(K, 1);
sysIds = [];
if haveSys && fromLog
    sysIds = unique(sysSlot(isfinite(sysSlot) & sysSlot > 0))';
elseif haveSys
    sysIds = unique(sysId(isfinite(sysId) & sysId > 0))';
end
nSys = numel(sysIds);
perSys = NaN(K, nSys);  nPerSys = zeros(K, nSys);  offSys = zeros(K, nSys);
nUsedSys = zeros(K, nSys);                                   % the log's used slots per system
if haveSys && fromLog
    for j = 1:nSys, nUsedSys(:, j) = sum(isfinite(prSlot) & sysSlot == sysIds(j), 2); end
end
rowsOf = cell(K, nSys);                                      % the record's range rows of each system
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
        for j = 1:nSys
            r = idx(rowSys == sysIds(j));
            rowsOf{k, j} = r;
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
            for j = 1:nSys, offSys(k, j) = nnz(offSysRow == sysIds(j)); end
        end
    else
        meanClean(k) = mean(yR);  spreadClean(k) = std(yR);
    end
end
offset  = meanOff - meanClean;
present = nOff >= 2 & abs(offset) >= MIN_GAP;
[kOn, kEnd] = runOf3(present);
onset = NaN;  tEnd = NaN;
if isfinite(kOn), onset = tKf(kOn);  tEnd = tKf(kEnd); end
sizeOnset = NaN;  sizeMax = NaN;  sysOff = NaN;
if isfinite(kOn)
    sizeOnset = median(offset(kOn:min(kOn + 2, K)));
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
    fprintf('   size at onset %+9.1f m (median of the first 3 epochs), largest %+.1f m\n', sizeOnset, sizeMax);
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
%% the mean innovation per system, epoch by epoch, where some system is off: the host's row set is visible here
offSysInn = NaN;  fitInn = NaN(K, 6);  offInn = false(K, 1);  stInn = zeros(0, 3);  kA = NaN;  kB = NaN;
if haveSys
    offInn = any(abs(perSys) >= MIN_GAP, 2);
    [stInn, iS] = stretches(offInn, max(abs(perSys), [], 2), STRETCH_GAP);
    if ~isempty(stInn)
        printStretches(tKf, stInn, iS, offInn, perSys, sysIds, MIN_GAP, gapS, 'innovation');
        kA = stInn(iS, 1);  kB = stInn(iS, 2);
        fprintf('   mean range innovation per system over the strongest stretch, %.2f s (epoch %d) to %.2f s (epoch %d)\n', tKf(kA), kA, tKf(kB), kB);
        fprintf('   (the first 8 epochs in a row, then spread to the end; a system without rows is a gap in the host''s row set):\n');
        printCourse(tKf, perSys, nPerSys, sysIds, kA, kB, {}, {});
        rowSets(nPerSys, perSys, sysIds, kA, kB, 'innovation');
        cnt = sum(abs(perSys(kA:kB, :)) >= MIN_GAP, 1);
        [~, jS] = max(cnt);  offSysInn = sysIds(jS);
        for k = 1:K
            fitInn(k, :) = fitRows(spfMeasLog(k), rowsOf{k, jS});
        end
        fprintf('   the position and clock the rows of system %d ask for, relative to the host''s prediction (least squares on the\n', offSysInn);
        fprintf('   record''s rows of that system, 4 or more; rms = what the rows leave after the fit, 0 when exactly determined):\n');
        printFit(tKf, fitInn, kA, kB);
    end
end
%% the profile against the reference position: no filter in the way
inj = [];  nSlots = [];  spread = [];  clock = NaN(K, 1);  clockGap = NaN(K, 1);  isb = [];
injSys = NaN;  injOnset = NaN;  injEnd = NaN;  injSize = NaN;  injMax = NaN;  nBadSlot = 0;  fitRef = NaN(K, 6);  offRef = false(K, 1);
stRef = zeros(0, 3);  kI = NaN;  kE = NaN;
if ~isempty(posRef)
    if ~(haveSys && fromLog)
        fprintf('   the reference position is used only with the pairing log as the third argument\n');
    elseif nSys == 0
        fprintf('   the pairing log holds no used range slot: no reference-based profile\n');
    else
        refK = refAtEpochs(tRef, posRef, tKf);
        M = size(prSlot, 2);
        res = NaN(K, M);  los = NaN(K, M, 3);
        for k = 1:K
            if ~all(isfinite(refK(k, :))), continue; end
            sat = reshape(satSlot(k, :, :), M, 3);
            d = sat - repmat(refK(k, :), M, 1);
            rho = sqrt(sum(d.^2, 2));
            res(k, :) = prSlot(k, :) - rho';
            los(k, :, :) = reshape(d ./ repmat(max(rho, 1), 1, 3), [1 M 3]);
        end
        bad = abs(res) > BAD_RES;
        nBadSlot = nnz(bad);
        if nBadSlot > 0
            [kb, ~] = find(bad);
            fprintf('   %d slot-epochs with a range more than %.0f km from the slot''s geometry dropped (first at %.2f s, epoch %d): there the range and the\n', ...
                nBadSlot, BAD_RES / 1e3, tKf(min(kb)), min(kb));
            fprintf('   position in the slot are not of the same satellite or instant\n');
            kb = unique(kb);
            nLeft = nnz(nRange(kb) < sum(isfinite(prSlot(kb, :)), 2));
            fprintf('   on %d of these %d epochs the record holds fewer range rows than the log has used slots: the host left a row out\n', nLeft, numel(kb));
            res(bad) = NaN;
        end
        mS = NaN(K, nSys);  nSlots = zeros(K, nSys);  spread = NaN(K, nSys);
        for j = 1:nSys
            cols = sysSlot == sysIds(j) & isfinite(res);
            for k = 1:K
                v = res(k, cols(k, :));
                if ~isempty(v)
                    mS(k, j) = median(v);  nSlots(k, j) = numel(v);
                    if numel(v) > 1, spread(k, j) = std(v); else, spread(k, j) = 0; end
                end
            end
        end
        % inter-system biases relative to the system present most often, as run medians (a spoof on fewer than
        % half of the common epochs leaves them untouched)
        [~, jRef] = max(sum(nSlots > 0, 1));
        isb = NaN(1, nSys);  isb(jRef) = 0;
        for j = 1:nSys
            if j ~= jRef, isb(j) = medianFinite(mS(:, j) - mS(:, jRef)); end
        end
        a = mS - repmat(isb, K, 1);
        % the receiver clock on the epochs where two systems agree, interpolated across the others
        agree = false(K, 1);
        for k = 1:K
            v = a(k, isfinite(a(k, :)));
            if numel(v) >= 2
                close = abs(v - median(v)) < MIN_GAP / 2;
                if nnz(close) >= 2, agree(k) = true;  clock(k) = mean(v(close)); end
            end
        end
        kAg = find(agree);
        fitLine = numel(kAg) < 2;
        if ~fitLine
            for k = 1:K
                if agree(k), clockGap(k) = 0;  continue; end
                [g, i] = min(abs(tKf(kAg) - tKf(k)));  clockGap(k) = g;
                if k < kAg(1) || k > kAg(end), clock(k) = clock(kAg(i));
                else, clock(k) = interp1(tKf(kAg), clock(kAg), tKf(k), 'linear'); end
            end
        else
            aMed = NaN(K, 1);
            for k = 1:K, aMed(k) = medianFinite(a(k, :)); end
            ok = isfinite(aMed);
            if nnz(ok) >= 2, p = polyfit(tKf(ok), aMed(ok), 1);  clock = polyval(p, tKf); end
        end
        inj = a - repmat(clock, 1, nSys);
        fprintf('   against the reference position: range minus geometric range per used slot, inter-system bias relative to system %d:%s\n', ...
            sysIds(jRef), sprintf(' system %d %+.1f m;', [sysIds; isb]));
        if fitLine
            fprintf('   fewer than 2 epochs where two systems agree within %.0f m: the clock is a straight line fitted over the run, and an offset\n', MIN_GAP / 2);
            fprintf('   common to all systems is partly absorbed in it\n');
        else
            fprintf('   the clock read on %d of %d epochs with ranges (two systems within %.0f m), interpolated across the others\n', ...
                numel(kAg), nnz(any(nSlots > 0, 2)), MIN_GAP / 2);
        end
        offRef = any(abs(inj) >= MIN_GAP, 2);
        [stRef, iR] = stretches(offRef, max(abs(inj), [], 2), STRETCH_GAP);
        if isempty(stRef)
            fprintf('   no stretch with a system %d m or more off the reference on 3 epochs\n', MIN_GAP);
        else
            printStretches(tKf, stRef, iR, offRef, inj, sysIds, MIN_GAP, gapS, 'offset');
            kI = stRef(iR, 1);  kE = stRef(iR, 2);
            injOnset = tKf(kI);  injEnd = tKf(kE);
            seg = kI:kE;
            cnt = sum(abs(inj(seg, :)) >= MIN_GAP, 1);
            [~, jS] = max(cnt);  injSys = sysIds(jS);
            injSize = medianFinite(inj(kI:min(kI + 2, K), jS));
            v = inj(seg, jS);  [~, iM] = max(abs(v));  injMax = v(iM);
            fprintf('   onset        %10.2f s (epoch %d), the strongest stretch\n', injOnset, kI);
            fprintf('   end          %10.2f s (epoch %d), %.1f s long\n', injEnd, kE, injEnd - injOnset);
            fprintf('   size at onset %+9.1f m on system %d (median of the first 3 epochs), largest %+.1f m\n', injSize, injSys, injMax);
            others = setdiff(1:nSys, jS);
            if ~isempty(others)
                fprintf('   the other systems in this stretch: median offset%s\n', sprintf(' system %d %+.1f m;', [sysIds(others); medianRows(inj(seg, others))]));
            end
            nI = nnz(clockGap(seg) > 0.3 & isfinite(inj(seg, jS)));
            if nI > 0
                fprintf('   on %d of the %d epochs the clock is interpolated (gaps up to %.1f s): the offset there rests on the clock staying on its line\n', ...
                    nI, numel(seg), max(clockGap(seg)));
            end
            fprintf('   within-system spread of system %d in this stretch: median %.1f m (at code noise all its satellites carry the same offset;\n', ...
                injSys, medianFinite(spread(seg, jS)));
            fprintf('   far above it only some do, or the spoof moves the position rather than the clock)\n');
            fprintf('   course (the first 8 epochs in a row, then spread to the end; i = clock interpolated):\n');
            printCourse(tKf, inj, nSlots, sysIds, kI, kE, {'clock [m]', sprintf('spread %d [m]', injSys)}, ...
                {@(k) sprintf('%+11.1f%s', clock(k), ifelse(clockGap(k) > 0.3, 'i', ' ')), @(k) sprintf('%11.1f', spread(k, jS))});
            rowSets(nSlots, inj, sysIds, kI, kE, 'offset');
            % the position and clock the off system's ranges describe: a least-squares fix on them per epoch
            C = nedFrame(posRef(1, :));
            for k = 1:K
                cols = find(sysSlot(k, :) == sysIds(jS) & isfinite(res(k, :)));
                if numel(cols) < 4, continue; end
                u = reshape(los(k, cols, :), numel(cols), 3);
                A = [-u, ones(numel(cols), 1)];
                r = res(k, cols)' - clock(k);
                x = A \ r;
                e = A * x - r;
                fitRef(k, :) = [(C * x(1:3))', x(4), sqrt(mean(e.^2)), numel(cols)];
            end
            fprintf('   the position and clock the ranges of system %d describe, relative to the reference (least squares on its used slots,\n', injSys);
            fprintf('   4 or more; receiver clock removed; rms = what the ranges leave after the fit, 0 when exactly determined):\n');
            printFit(tKf, fitRef, kI, kE);
            f = fitRef(seg, :);  ok = isfinite(f(:, 1));
            if any(ok)
                fprintf('   over the stretch: median position %+.1f N %+.1f E %+.1f D m (horizontal %.1f m), clock %+.1f m, rms %.1f m; the position\n', ...
                    median(f(ok, 1)), median(f(ok, 2)), median(f(ok, 3)), sqrt(median(f(ok, 1))^2 + median(f(ok, 2))^2), median(f(ok, 4)), median(f(ok, 5)));
                fprintf('   moves within %.1f N %.1f E %.1f D m (10th to 90th percentile) on %d fitted epochs\n', ...
                    pctRange(f(ok, 1)), pctRange(f(ok, 2)), pctRange(f(ok, 3)), nnz(ok));
            end
        end
    end
end
%% the two fits side by side: same sign, opposite sign (the host's position columns carry +LOS), or apart
fitSign = 0;
inSeg = false(K, 1);
if isfinite(kI), inSeg(kI:kE) = true; end
both = isfinite(fitInn(:, 1)) & isfinite(fitRef(:, 1)) & offRef & inSeg;
if nnz(both) >= 3
    hI = median(fitInn(both, 1:2), 1);  hR = median(fitRef(both, 1:2), 1);
    dSame = norm(hI - hR);  dOpp = norm(hI + hR);  sz = norm(hR);
    if sz < MIN_GAP
        fitSign = 1;
        fprintf('   both fits put the offset in the clock: position within %.0f m of the reference, clock %+.0f m in the innovations and %+.0f m\n', ...
            max(norm(hI), sz), median(fitInn(both, 4)), median(fitRef(both, 4)));
        fprintf('   in the ranges, on %d off epochs\n', nnz(both));
    elseif dOpp < 0.2 * sz && dOpp < dSame
        fitSign = -1;
        fprintf('   the two fits agree within %.0f m horizontal with the sign reversed: the position columns of the host''s H carry +LOS (its\n', dOpp);
        fprintf('   position error state is estimate minus truth), so the innovation fit reads as minus the position the ranges describe\n');
    elseif dSame < 0.2 * sz
        fitSign = 1;
        fprintf('   the two fits agree within %.0f m horizontal on %d off epochs\n', dSame, nnz(both));
    else
        fprintf('   the two fits differ by %.0f m horizontal (reversed: %.0f m) on %d off epochs: the host''s prediction sat that far from the\n', dSame, dOpp, nnz(both));
        fprintf('   reference, or the record''s rows are not the log''s slots\n');
    end
end
%% the timeline around the spoof: the row set of each epoch, and the gate mode under it
nRecNoRows = NaN;  nLogNoSlots = NaN;
kF = min([kA, kI]);  kL = max([kB, kE]);
if haveSys && isfinite(kF)
    k1 = max(1, kF - 20);  k2 = min(K, kL + 20);
    fprintf('   timeline around the strongest stretch, %.2f s (epoch %d) to %.2f s (epoch %d), one character per epoch, 50 per line: the record''s row set (rows)\n', tKf(k1), k1, tKf(k2), k2);
    fprintf('   (. = no range rows, A = all %d systems, a digit = that system alone, x = some systems)', nSys);
    if fromLog, fprintf(', the log''s used slots under it where they differ (log)'); end
    modes = [];
    if ~isempty(spfLog) && isstruct(spfLog) && isfield(spfLog, 'mode')
        nL = numel(spfLog.mode);
        if isfield(spfLog, 'count'), nL = min(nL, double(spfLog.count)); end
        if nL == K
            modes = spfLog.mode(1:K);
            fprintf(', the gate mode under it (N nominal, C coast, P probation)');
        else
            fprintf('; the gate log has %d epochs, the records %d: no mode line', nL, K);
        end
    end
    fprintf(':\n');
    mixes = {};  mixCnt = [];
    for k0 = k1:50:k2
        kk = k0:min(k0 + 49, k2);
        line = blanks(numel(kk));  logLine = blanks(numel(kk));
        for i = 1:numel(kk)
            on = nPerSys(kk(i), :) > 0;
            if ~any(on), line(i) = '.';
            elseif all(on), line(i) = 'A';
            elseif nnz(on) == 1, line(i) = char('0' + min(sysIds(on), 9));
            else
                line(i) = 'x';
                name = strtrim(sprintf('%d ', sysIds(on)));
                m = find(strcmp(mixes, name));
                if isempty(m), mixes{end+1} = name;  mixCnt(end+1) = 1; else, mixCnt(m) = mixCnt(m) + 1; end %#ok<AGROW>
            end
            onL = nUsedSys(kk(i), :) > 0;
            if ~any(onL), logLine(i) = '.';
            elseif all(onL), logLine(i) = 'A';
            elseif nnz(onL) == 1, logLine(i) = char('0' + min(sysIds(onL), 9));
            else, logLine(i) = 'x'; end
        end
        fprintf('   %9.2f %-5s %s\n', tKf(k0), 'rows', grouped(line));
        if fromLog && ~strcmp(line, logLine), fprintf('   %9s %-5s %s\n', '', 'log', grouped(logLine)); end
        if ~isempty(modes)
            ml = blanks(numel(kk));
            for i = 1:numel(kk)
                switch round(modes(kk(i)))
                    case 1, ml(i) = 'N';
                    case 2, ml(i) = 'C';
                    case 3, ml(i) = 'P';
                    otherwise, ml(i) = '?';
                end
            end
            fprintf('   %9s %-5s %s\n', '', 'mode', grouped(ml));
        end
    end
    if nSys == 1, fprintf('   (one system: A and 1 are the same)\n'); end
    if fromLog
        w = k1:k2;
        nRecNoRows = nnz(nRange(w) == 0 & any(nUsedSys(w, :) > 0, 2));
        nLogNoSlots = nnz(nRange(w) > 0 & ~any(nUsedSys(w, :) > 0, 2));
        fprintf('   in this window the record holds no range row on %d epochs where the log has used slots (the host built none), and rows on %d\n', nRecNoRows, nLogNoSlots);
        fprintf('   epochs where the log has no used slot\n');
        nRejW = sum(nRejSlot(w, :) > 0, 1);
        if any(nRejW)
            fprintf('   the log holds ranges the host did not use (satValid false):%s\n', sprintf(' system %d on %d epochs;', [1:numel(nRejW); nRejW]));
        else
            fprintf('   every range in the log was used in this window: no slot with a range and satValid false\n');
        end
    end
    if ~isempty(mixes)
        fprintf('   x stands for:');
        for m = 1:numel(mixes), fprintf(' systems %s on %d epochs;', mixes{m}, mixCnt(m)); end
        fprintf('\n');
    end
    if ~isempty(modes)
        ev = {'evLatch', 'latch'; 'evProbStart', 'probation start'; 'evVeto', 'probation vetoed'; 'evHandback', 'hand-back'};
        parts = {};
        for i = 1:size(ev, 1)
            if isfield(spfLog, ev{i, 1})
                for e = spfLog.(ev{i, 1})(:)'
                    if e >= 1 && e <= K, parts{end+1} = sprintf('%s at %.2f s (epoch %d)', ev{i, 2}, tKf(e), e); end %#ok<AGROW>
                end
            end
        end
        if ~isempty(parts), fprintf('   gate events: %s\n', strjoin(parts, '; ')); end
    end
end
prof = struct('t', tKf, 'nRange', nRange, 'nOff', nOff, 'offset', offset, 'meanOff', meanOff, 'meanClean', meanClean, ...
    'spreadClean', spreadClean, 'present', present, 'onset', onset, 'tEnd', tEnd, 'sizeOnset', sizeOnset, 'sizeMax', sizeMax, ...
    'sysIds', sysIds, 'perSys', perSys, 'nPerSys', nPerSys, 'sysOff', sysOff, 'nSysMismatch', nSysMismatch, ...
    'offSysInn', offSysInn, 'fitInn', fitInn, 'fitSign', fitSign, 'nRecNoRows', nRecNoRows, 'nLogNoSlots', nLogNoSlots, ...
    'stretchesInn', stInn, 'stretchesRef', stRef, ...
    'inj', inj, 'nSlots', nSlots, 'spread', spread, 'clock', clock, 'clockGap', clockGap, 'isb', isb, ...
    'injSys', injSys, 'injOnset', injOnset, 'injEnd', injEnd, 'injSize', injSize, 'injMax', injMax, 'nBadSlot', nBadSlot, 'fitRef', fitRef);
end

%------------------------------------------------------------------------
function [kOn, kEnd] = runOf3(flag)
% first and last epoch of the stretch between the first and the last run of 3 consecutive true entries
kOn = NaN;  kEnd = NaN;
if numel(flag) < 3, return; end
run3 = conv(double(flag(:)), ones(3, 1), 'valid') == 3;
if any(run3), kOn = find(run3, 1);  kEnd = find(run3, 1, 'last') + 2; end
end

function [st, iStrong] = stretches(flag, strength, maxGap)
% runs of true entries of flag, merged when maxGap or fewer false entries separate them; st rows [kStart kEnd nOff]
% for the stretches with 3 off epochs or more; iStrong = the row whose off epochs have the largest median strength
st = zeros(0, 3);  iStrong = [];  med = [];
k = find(flag(:))';
if isempty(k), return; end
breaks = find(diff(k) > maxGap + 1);
starts = k([1, breaks + 1]);  ends = k([breaks, numel(k)]);
for i = 1:numel(starts)
    seg = starts(i):ends(i);
    n = nnz(flag(seg));
    if n >= 3
        st(end+1, :) = [starts(i), ends(i), n]; %#ok<AGROW>
        v = strength(seg);  v = v(flag(seg));
        med(end+1) = medianFinite(v); %#ok<AGROW>
    end
end
if ~isempty(st), [~, iStrong] = max(med); end
end

function printStretches(tKf, st, iS, flag, vals, sysIds, minGap, gapS, what)
% one line per stretch: when, how many off epochs, the median value per system over them, the strongest marked
fprintf('   stretches with a system %d m or more off (off epochs closer than %.0f s belong to one stretch):\n', minGap, gapS);
for i = 1:size(st, 1)
    seg = st(i, 1):st(i, 2);  off = seg(flag(seg));
    parts = '';
    for j = 1:numel(sysIds)
        v = vals(off, j);  v = v(isfinite(v) & abs(v) >= minGap);
        if ~isempty(v), parts = [parts, sprintf(' system %d %+.1f m on %d;', sysIds(j), median(v), numel(v))]; end %#ok<AGROW>
    end
    mark = '';
    if i == iS, mark = '   <- the tables below'; end
    fprintf('   %d) %.2f to %.2f s (epochs %d to %d), %d off epochs, median %s:%s%s\n', ...
        i, tKf(st(i, 1)), tKf(st(i, 2)), st(i, 1), st(i, 2), st(i, 3), what, parts, mark);
end
end

function [f] = fitRows(s, rows)
% the position (in the convention of H's position columns) and clock the given rows of one record ask for:
% least squares of the innovation on [H(:, 1:3), 1]; f = [pos(1:3), clock, rms, rows], NaN with fewer than 4 rows
f = NaN(1, 6);
if numel(rows) < 4, return; end
A = [s.obsMatrix(rows, 1:3), ones(numel(rows), 1)];
y = s.innovation(rows);
if rank(A) < 4, return; end
x = A \ y;
e = A * x - y;
f = [x', sqrt(mean(e.^2)), numel(rows)];
end

function printFit(tKf, fit, kA, kB)
pick = unique([kA:min(kA + 7, kB), round(linspace(kA, kB, 8))]);
fprintf('   %10s  %10s  %10s  %10s  %10s  %8s  %4s\n', 't [s]', 'N [m]', 'E [m]', 'D [m]', 'clock [m]', 'rms [m]', 'rows');
for k = pick
    if isfinite(fit(k, 1))
        fprintf('   %10.2f  %+10.1f  %+10.1f  %+10.1f  %+10.1f  %8.1f  %4d\n', tKf(k), fit(k, 1:5), fit(k, 6));
    else
        fprintf('   %10.2f  %10s  %10s  %10s  %10s  %8s  %4s\n', tKf(k), '-', '-', '-', '-', '-', '<4');
    end
end
end

function printCourse(tKf, vals, nRows, sysIds, kA, kB, extraHdr, extraFcn)
% a table over kA:kB: t, per system the value with its rows, then the extra columns
pick = unique([kA:min(kA + 7, kB), round(linspace(kA, kB, 8))]);
fprintf('   %10s', 't [s]');
for j = 1:numel(sysIds), fprintf('  %18s', sprintf('system %d [m] (rows)', sysIds(j))); end
for i = 1:numel(extraHdr), fprintf('  %12s', extraHdr{i}); end
fprintf('\n');
for k = pick
    fprintf('   %10.2f', tKf(k));
    for j = 1:numel(sysIds)
        if nRows(k, j) > 0, fprintf('  %+13.1f (%2d)', vals(k, j), nRows(k, j));
        else, fprintf('  %13s (%2d)', '-', 0); end
    end
    for i = 1:numel(extraFcn), fprintf('  %12s', extraFcn{i}(k)); end
    fprintf('\n');
end
end

function rowSets(nRows, vals, sysIds, kA, kB, what)
% which systems have rows on the epochs kA:kB, how often, and the median value per system on those epochs
seg = kA:kB;
sig = double(nRows(seg, :) > 0);
[u, ~, ic] = unique(sig, 'rows');
cnt = accumarray(ic(:), 1);
[cnt, order] = sort(cnt, 'descend');
fprintf('   row sets in this stretch (%d epochs):', numel(seg));
for i = 1:numel(cnt)
    rows = seg(ic == order(i));
    if ~any(u(order(i), :))
        fprintf(' no range rows on %d;', cnt(i));
    else
        js = find(u(order(i), :));
        fprintf(' systems %s on %d (median %s%s);', strtrim(sprintf('%d ', sysIds(js))), cnt(i), what, ...
            sprintf(' %d: %+.1f', [sysIds(js); medianRows(vals(rows, js))]));
    end
end
fprintf('\n');
end

function [s] = grouped(line)
% the characters in groups of 10
s = '';
for i = 1:10:numel(line), s = [s, line(i:min(i + 9, numel(line))), ' ']; end %#ok<AGROW>
s = strtrim(s);
end

function [refK] = refAtEpochs(tRef, posRef, tKf)
% the reference position in ECEF at each epoch: one row for a stationary run, else interpolated on tRef
K = numel(tKf);
if size(posRef, 2) ~= 3 || max(abs(posRef(:, 1))) > pi / 2 || max(abs(posRef(:, 2))) > pi || max(abs(posRef(:, 3))) > 2e4
    error('spoofProfile:posRef', 'posRef must be the reference [lat lon h] in radians, radians and metres; got %d x %d', size(posRef, 1), size(posRef, 2));
end
e = llaToEcef(posRef);
if size(posRef, 1) == 1
    refK = repmat(e, K, 1);
else
    if numel(tRef) ~= size(posRef, 1)
        error('spoofProfile:tRef', 'tRef has %d entries, posRef %d rows', numel(tRef), size(posRef, 1));
    end
    refK = NaN(K, 3);
    for c = 1:3, refK(:, c) = interp1(tRef(:), e(:, c), tKf, 'linear', NaN); end
end
end

function [e] = llaToEcef(p)
a = 6378137;  e2 = 6.69437999014e-3;
L = p(:, 1);  l = p(:, 2);  h = p(:, 3);
N = a ./ sqrt(1 - e2 * sin(L).^2);
e = [(N + h) .* cos(L) .* cos(l), (N + h) .* cos(L) .* sin(l), (N * (1 - e2) + h) .* sin(L)];
end

function [C] = nedFrame(p)
% rows: the north, east and down unit vectors in ECEF at [lat lon h]
L = p(1);  l = p(2);
C = [-sin(L) * cos(l), -sin(L) * sin(l), cos(L); -sin(l), cos(l), 0; -cos(L) * cos(l), -cos(L) * sin(l), -sin(L)];
end

function [r] = pctRange(x)
x = sort(x(isfinite(x)));
if numel(x) < 2, r = 0; return; end
r = x(max(1, round(0.9 * numel(x)))) - x(max(1, round(0.1 * numel(x))));
end

function [m] = medianRows(X)
% the median of each column over its finite entries, NaN for a column without any
m = NaN(1, size(X, 2));
for j = 1:size(X, 2), m(j) = medianFinite(X(:, j)); end
end

function [m] = medianFinite(x)
x = x(isfinite(x));
if isempty(x), m = NaN; else, m = median(x); end
end

function [s] = ifelse(c, a, b)
if c, s = a; else, s = b; end
end
