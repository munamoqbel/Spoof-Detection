function [prRate, satVel, have] = pairingLogRates(loggedData)
%PAIRINGLOGRATES  The range rate per slot and the satellite velocity per slot from the host's pairing log.
%   [prRate, satVel, have] = pairingLogRates(loggedData)
%   loggedData as pairingLogUnpack takes it: prLog with a geomRangeRate [S x 1] per constellation and satPosLog
%   with a velocity [S x 3] per constellation. The slot layout is that of pairingLogUnpack (the systems' slots in
%   order). prRate [K x M], satVel [K x M x 3] in m/s; velocities given in km/s are scaled to m/s, with a note.
%   have is false when the log carries no geomRangeRate or no velocity field, or no finite rate at all.
%   MATLAB desktop tool, not codegen. Runs under Octave.
if iscell(loggedData), loggedData = [loggedData{:}]; end
K = numel(loggedData);  prRate = [];  satVel = [];  have = false;
k0 = 1;
while k0 < K && ~(isstruct(loggedData(k0).prLog) && isstruct(loggedData(k0).satPosLog)), k0 = k0 + 1; end
prFields = fieldnames(loggedData(k0).prLog);  spFields = fieldnames(loggedData(k0).satPosLog);
nSys = min(numel(prFields), numel(spFields));
if ~isfield(loggedData(k0).prLog.(prFields{1}), 'geomRangeRate') || ~isfield(loggedData(k0).satPosLog.(spFields{1}), 'velocity'), return; end
slots = zeros(1, nSys);
for s = 1:nSys, slots(s) = numel(loggedData(k0).prLog.(prFields{s}).geomRange); end
first = cumsum([0 slots(1:end-1)]);  M = sum(slots);
prRate = NaN(K, M);  satVel = NaN(K, M, 3);
for k = 1:K
    if ~(isstruct(loggedData(k).prLog) && isstruct(loggedData(k).satPosLog)), continue; end
    for s = 1:nSys
        cols = first(s) + (1:slots(s));
        pf = loggedData(k).prLog.(prFields{s});  st = loggedData(k).satPosLog.(spFields{s});
        if ~isfield(pf, 'geomRangeRate') || ~isfield(st, 'velocity'), continue; end
        r = double(pf.geomRangeRate(:));  v = double(st.velocity);
        if numel(r) ~= slots(s) || size(v, 1) ~= slots(s) || size(v, 2) ~= 3, continue; end
        prRate(k, cols) = r';  satVel(k, cols, :) = reshape(v, [1 slots(s) 3]);
    end
end
spd = sqrt(sum(satVel.^2, 3));  spd = spd(isfinite(spd) & spd > 0);
if isempty(spd), return; end
if median(spd) < 10
    satVel = satVel * 1e3;
    fprintf('   (satellite velocities read as km/s, median %.2f: scaled to m/s)\n', median(spd));
end
have = any(isfinite(prRate(:)));
end
