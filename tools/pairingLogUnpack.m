function [pr, satPos, sysId] = pairingLogUnpack(loggedData, pairing)
%PAIRINGLOGUNPACK  The arrays pairingCheck needs, from the host's per-epoch log of the range rows.
%   [pr, satPos, sysId] = pairingLogUnpack(loggedData)
%   [pr, satPos, sysId] = pairingLogUnpack(loggedData, 'order')
%     pairing 'slot' (default): the range in slot j pairs with the position in slot j
%     pairing 'order': the i-th used range of a system pairs with the i-th live position row
%     (norm above 1000 km) of that system, for a log whose ranges are compacted while the
%     positions stay in channel slots; both land in the range's slot
%   loggedData [K x 1] struct array (or cell of structs), one record per 2 Hz epoch, with
%     prLog      one field per constellation, e.g. gpsRange, gloRange, galRange, each a struct
%                whose geomRange [S x 1] is the corrected pseudorange per satellite slot
%                (satellite clock and atmosphere removed, receiver clock left in), 0 where the
%                slot is not used
%     satPosLog  one field per constellation in the SAME order, e.g. navstarStatus,
%                glonassStatus, galileoStatus, each a struct whose position [S x 3] is the
%                satellite ECEF position of the slot as used in h(x) (m); a satValid [S x 1]
%                field, when present, also gates the slot
%     sysIdLog   optional [m x 1] constellation id per range row of the epoch; only checked
%                against the used slots per system
%   The k-th field of prLog pairs with the k-th field of satPosLog as system k. A slot is used
%   when its geomRange is finite and not 0 (and satValid where given); unused slots are NaN in
%   pr, an epoch without GNSS a row of NaN. A constellation switched off (all zeros) contributes
%   no slot. Positions given in km are scaled to metres, with a note.
%     pr     [K x M]      M = the slots of all systems
%     satPos [K x M x 3]  ECEF (m)
%     sysId  [K x M]
%   MATLAB desktop tool, not codegen. Runs under Octave.

if nargin < 2 || isempty(pairing), pairing = 'slot'; end
byOrder = strcmpi(pairing, 'order');
if iscell(loggedData), loggedData = [loggedData{:}]; end
K = numel(loggedData);
nShort = 0;
k0 = 1;
while k0 < K && ~(isstruct(loggedData(k0).prLog) && isstruct(loggedData(k0).satPosLog))
    k0 = k0 + 1;
end
prFields = fieldnames(loggedData(k0).prLog);
spFields = fieldnames(loggedData(k0).satPosLog);
nSys = min(numel(prFields), numel(spFields));
if numel(prFields) ~= numel(spFields)
    fprintf('pairingLogUnpack: WARNING prLog has %d constellations, satPosLog %d; the first %d are paired in order\n', ...
        numel(prFields), numel(spFields), nSys);
end
slots = zeros(1, nSys);
for s = 1:nSys, slots(s) = numel(loggedData(k0).prLog.(prFields{s}).geomRange); end
first = cumsum([0 slots(1:end-1)]);
M = sum(slots);
pr = NaN(K, M);  satPos = NaN(K, M, 3);  sysId = NaN(K, M);
nUsed = zeros(K, nSys);  nIdLog = NaN(K, nSys);
hasId = isfield(loggedData, 'sysIdLog');
for k = 1:K
    if ~(isstruct(loggedData(k).prLog) && isstruct(loggedData(k).satPosLog)), continue; end
    for s = 1:nSys
        cols = first(s) + (1:slots(s));
        g  = double(loggedData(k).prLog.(prFields{s}).geomRange(:));
        st = loggedData(k).satPosLog.(spFields{s});
        p  = double(st.position);
        if numel(g) ~= slots(s) || size(p, 1) ~= slots(s) || size(p, 2) ~= 3
            error('pairingLogUnpack:slots', 'epoch %d, system %d: geomRange has %d slots and position is %d x %d; expected %d x 1 and %d x 3', ...
                k, s, numel(g), size(p, 1), size(p, 2), slots(s), slots(s));
        end
        used = isfinite(g) & g ~= 0;
        if isfield(st, 'satValid') && numel(st.satValid) == slots(s)
            used = used & logical(st.satValid(:));
        end
        g(~used) = NaN;
        if byOrder
            live = find(sqrt(sum(p.^2, 2)) > 1e6 & all(isfinite(p), 2));   % position rows that hold a satellite
            iu   = find(used);
            n    = min(numel(iu), numel(live));
            if n < numel(iu), nShort = nShort + 1; end
            q = NaN(slots(s), 3);
            q(iu(1:n), :) = p(live(1:n), :);
            g(iu(n+1:end)) = NaN;
            p = q;
        end
        pr(k, cols) = g';
        satPos(k, cols, :) = reshape(p, [1 slots(s) 3]);
        sysId(k, cols) = s;
        nUsed(k, s) = nnz(isfinite(g));
    end
    if hasId && ~isempty(loggedData(k).sysIdLog)
        ids = double(loggedData(k).sysIdLog(:));
        for s = 1:nSys, nIdLog(k, s) = nnz(ids == s); end
    end
end

% satellite positions in km would sit 2.0e4..2.7e4 from the Earth's centre
rad = sqrt(sum(satPos.^2, 3));
rad = rad(isfinite(pr));
if ~isempty(rad) && median(rad) < 1e5
    satPos = satPos * 1e3;
    fprintf('pairingLogUnpack: NOTE satellite positions read as km (median radius %.0f): scaled to metres\n', median(rad));
end

fprintf('pairingLogUnpack: %d epochs, %d with GNSS; systems:', K, nnz(any(isfinite(pr), 2)));
for s = 1:nSys
    if any(nUsed(:, s))
        fprintf(' %d = %s / %s (%d slots, median %d used)', s, prFields{s}, spFields{s}, slots(s), round(median(nUsed(nUsed(:, s) > 0, s))));
    else
        fprintf(' %d = %s / %s (switched off)', s, prFields{s}, spFields{s});
    end
end
fprintf('\n');
if byOrder
    fprintf('pairingLogUnpack: ranges paired with the live position rows by order');
    if nShort > 0, fprintf('; %d epochs had fewer live positions than ranges (the extra ranges dropped)', nShort); end
    fprintf('\n');
end
chk = all(isfinite(nIdLog), 2);
bad = find(chk & any(nIdLog ~= nUsed, 2));
if ~isempty(bad)
    fprintf('pairingLogUnpack: WARNING %d of %d epochs with sysIdLog count the rows per system differently from the used slots\n', numel(bad), nnz(chk));
    fprintf('   (first at epoch %d: sysIdLog %s, slots %s): the host and the log disagree on which satellites are used\n', ...
        bad(1), mat2str(nIdLog(bad(1), :)), mat2str(nUsed(bad(1), :)));
end
end
