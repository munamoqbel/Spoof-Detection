function [res] = spoofInjectCheck(loggedSpoofed, loggedClean, tKf, tRef, posRef, scn)
%SPOOFINJECTCHECK  What the injection put into the host's ranges: the spoofed run's pairing log against its clean twin's.
%   res = spoofInjectCheck(loggedSpoofed, loggedClean, tKf, [], posRef)   posRef [lat lon h] (rad, rad, m), one row
%                                                 for a stationary receiver; tKf the epoch times of the spoofed run
%   res = spoofInjectCheck(loggedSpoofed, loggedClean, tKf, tRef, posRef)   a moving receiver: posRef [N x 3] on tRef
%   res = spoofInjectCheck(..., scn)              the scenario struct of spoofScenario as well: what the ranges carry
%                                                 is printed next to what was intended, epoch by epoch
%   The two logs are the host's pairing logs (prLog / satPosLog per epoch, as pairingCheck takes them) of the spoofed
%   run and of its clean twin, the same inputs replayed without the injection. Epochs are matched on the satellite
%   positions in the slots (the same inputs give the same sky), so a shift of the time base between the two runs does
%   not matter. On a matched epoch each used slot's range in the spoofed log minus the same slot's range in the clean
%   log is the injection and nothing else: no receiver error, no reference error, no filter in the way. A
%   least-squares fix on these differences per system ([-u, 1]: position offset in NED at the reference and clock
%   offset) reads the injected vector back exactly, and the range rates the same way (velocity offset and clock
%   drift). The reference only gives the lines of sight, so a few hundred metres of error in it do not matter.
%   Prints: the epochs matched; per system the epochs with a changed slot, the first and the last; the course of the
%   fitted offset over the changed epochs; the medians over them; with scn what was intended, the difference, the
%   epochs where something was intended but nothing changed, and the epochs with a change while nothing was intended
%   (the logs are not twins, or the host's ranges depend on its own state).
%   Result: t, kClean (the clean epoch matched to each spoofed epoch, NaN when none), sysIds, nUsed and nChanged
%   [K x nSys], fit [K x 6 x nSys] (N, E, D, clock, rms, slots), fitRate [K x 6 x nSys] (vN, vE, vD, drift, rms,
%   slots), intended [K x 8] (dPosNed, dClk, dVelNed, dClkRate from scn, NaN without one), tFirst and tLast [1 x nSys]
%   (first and last changed epoch per system), medFit and medRate [nSys x 4] (medians over the changed epochs),
%   nMissing [1 x nSys] (epochs with an offset intended and no changed slot), nOutside [1 x nSys] (epochs with a
%   changed slot and nothing intended), shift (median of kClean minus the spoofed index).
%   A slot counts as changed when its range differs by more than 1 cm or its rate by more than 1 mm/s.
%   MATLAB desktop tool, not codegen. Runs under Octave.

if nargin < 4, tRef = []; end
if nargin < 5, posRef = []; end
if nargin < 6, scn = []; end
TOL_CHANGE = 1e-2;                                           % m: a range that differs by more than this changed (a correction that
                                                             % depends on the host's own state moves a range by millimetres, not this)
TOL_RATE   = 1e-3;                                           % m/s: the same for a range rate
TOL_MATCH  = 1.0;                                            % m: satellite positions within this are the same sky
SEARCH     = 50;                                             % epochs: how far from the same index the twin is searched
[prS, satS, sysS] = pairingLogUnpack(loggedSpoofed);
[prC, satC] = pairingLogUnpack(loggedClean);
K = size(prS, 1);  KC = size(prC, 1);  M = size(prS, 2);
if size(prC, 2) ~= M
    error('spoofInjectCheck:slots', 'the spoofed log has %d slots, the clean %d: not the same host', M, size(prC, 2));
end
tKf = tKf(:);
if numel(tKf) ~= K
    error('spoofInjectCheck:tKf', 'tKf has %d entries, the spoofed log %d epochs', numel(tKf), K);
end
if isempty(posRef)
    error('spoofInjectCheck:posRef', 'posRef is needed for the lines of sight: [lat lon h] in radians, radians, metres');
end
sysIds = unique(sysS(isfinite(prS) & isfinite(sysS) & sysS > 0))';
nSys = numel(sysIds);
%% the clean epoch of each spoofed epoch: the same satellite positions in the slots
kClean = NaN(K, 1);
for k = 1:K
    used = find(isfinite(prS(k, :)));
    if numel(used) < 4, continue; end
    sk = reshape(satS(k, used, :), [1 numel(used) 3]);
    cand = max(1, k - SEARCH):min(KC, k + SEARCH);
    if isempty(cand), continue; end
    cs = satC(cand, used, :);
    d = max(max(abs(cs - repmat(sk, [numel(cand) 1 1])), [], 3), [], 2);
    nCommon = sum(all(isfinite(cs), 3), 2);
    hit = cand(d < TOL_MATCH & nCommon >= 4);
    if ~isempty(hit), [~, i] = min(abs(hit - k));  kClean(k) = hit(i); end
end
matched = isfinite(kClean);
withRanges = sum(isfinite(prS), 2) >= 4;
shift = NaN;
if any(matched), shift = median(kClean(matched) - find(matched)); end
fprintf('spoofInjectCheck: %d spoofed epochs, %d clean; %d of the %d epochs with 4 used slots or more matched on the satellite positions\n', ...
    K, KC, nnz(matched), nnz(withRanges));
fprintf('   (within %.0f m); the clean epoch sits %+d epochs from the spoofed index (median); a slot counts as changed when its range differs\n', ...
    TOL_MATCH, round(shift));
fprintf('   by more than %.0f mm or its rate by more than %.1f mm/s\n', TOL_CHANGE * 1e3, TOL_RATE * 1e3);
if nnz(matched) < nnz(withRanges)
    fprintf('   %d epochs with ranges have no clean twin within %d epochs: not the same inputs there (a different sky)\n', ...
        nnz(withRanges & ~matched), SEARCH);
end
%% the differences and their fits
refK = refAtEpochs(tRef, posRef, tKf);
C = nedFrame(mean(posRef, 1));
dPr = NaN(K, M);
dPr(matched, :) = prS(matched, :) - prC(kClean(matched), :);
[rrS, ~, haveS] = pairingLogRates(loggedSpoofed);
[rrC, ~, haveC] = pairingLogRates(loggedClean);
haveRate = haveS && haveC && isequal(size(rrS), size(prS)) && isequal(size(rrC), size(prC));
dRr = NaN(K, M);
if haveRate, dRr(matched, :) = rrS(matched, :) - rrC(kClean(matched), :); end
nUsed = zeros(K, nSys);  nChanged = zeros(K, nSys);  fit = NaN(K, 6, nSys);  fitRate = NaN(K, 6, nSys);
for k = find(matched)'
    if ~all(isfinite(refK(k, :))), continue; end
    for j = 1:nSys
        cols = find(sysS(k, :) == sysIds(j) & isfinite(dPr(k, :)));
        nUsed(k, j) = numel(cols);
        nChanged(k, j) = nnz(abs(dPr(k, cols)) > TOL_CHANGE | abs(dRr(k, cols)) > TOL_RATE);
        if numel(cols) < 4, continue; end
        sat = reshape(satS(k, cols, :), numel(cols), 3);
        d = sat - repmat(refK(k, :), numel(cols), 1);
        u = d ./ repmat(sqrt(sum(d.^2, 2)), 1, 3);
        A = [-u, ones(numel(cols), 1)];
        fit(k, :, j) = solve(A, dPr(k, cols)', C);
        if haveRate && all(isfinite(dRr(k, cols))), fitRate(k, :, j) = solve(A, dRr(k, cols)', C); end
    end
end
changedEp = nChanged > 0;
%% what was intended
intended = NaN(K, 8);
if ~isempty(scn)
    for k = 1:K
        inj = spoofScenario(tKf(k), scn);
        intended(k, :) = [inj.dPosNed(:)', inj.dClk, inj.dVelNed(:)', inj.dClkRate];
    end
end
somethingIntended = any(abs(intended) > 0, 2);
%% per system: when, how many slots, the course and the medians
tFirst = NaN(1, nSys);  tLast = NaN(1, nSys);  medFit = NaN(nSys, 4);  medRate = NaN(nSys, 4);
nMissing = NaN(1, nSys);  nOutside = NaN(1, nSys);
fprintf('\n--- per system: the epochs with a changed slot ---------------------------------------------------------------\n');
for j = 1:nSys
    ep = find(changedEp(:, j));
    if isempty(ep)
        fprintf('   system %d: no slot ever changed on %d matched epochs with used slots\n', sysIds(j), nnz(nUsed(:, j) > 0));
    else
        tFirst(j) = tKf(ep(1));  tLast(j) = tKf(ep(end));
        fprintf('   system %d: changed slots on %d epochs, first at %.2f s (epoch %d), last at %.2f s (epoch %d); on those epochs %d of %d used\n', ...
            sysIds(j), numel(ep), tFirst(j), ep(1), tLast(j), ep(end), round(median(nChanged(ep, j))), round(median(nUsed(ep, j))));
        f = fit(ep, :, j);  ok = isfinite(f(:, 1));
        if any(ok), medFit(j, :) = median(f(ok, 1:4), 1); end
        fr = fitRate(ep, :, j);  okR = isfinite(fr(:, 1));
        if any(okR), medRate(j, :) = median(fr(okR, 1:4), 1); end
        fprintf('   the offset the differences describe, medians over those epochs: N %+.2f E %+.2f D %+.2f clock %+.2f m', medFit(j, :));
        if any(okR), fprintf(', rates vN %+.3f vE %+.3f vD %+.3f drift %+.3f m/s', medRate(j, :)); end
        fprintf('\n');
        partial = ep(nChanged(ep, j) < nUsed(ep, j));
        if ~isempty(partial)
            fprintf('   on %d of these epochs only some of the system''s used slots changed (first at %.2f s): the injection reached part of the\n', ...
                numel(partial), tKf(partial(1)));
            fprintf('   constellation, or the clean twin differs on a slot\n');
        end
    end
end
anyChanged = any(changedEp, 2);
if any(anyChanged)
    fprintf('\n--- the course over the changed epochs (the first 8 in a row, then spread to the end) ----------------------------\n');
    ep = find(anyChanged);
    pick = unique([ep(1:min(8, numel(ep)))', ep(round(linspace(1, numel(ep), min(8, numel(ep)))))']);
    fprintf('   %10s', 't [s]');
    for j = 1:nSys, fprintf('  | sys %d chg/used %8s %8s %8s %8s %6s', sysIds(j), 'N [m]', 'E [m]', 'D [m]', 'clk [m]', 'rms'); end
    if ~isempty(scn), fprintf('  | intended N E D clock [m]'); end
    fprintf('\n');
    for k = pick
        fprintf('   %10.2f', tKf(k));
        for j = 1:nSys
            if isfinite(fit(k, 1, j))
                fprintf('  |       %2d/%2d   %+8.2f %+8.2f %+8.2f %+8.2f %6.2f', nChanged(k, j), nUsed(k, j), fit(k, 1:5, j));
            else
                fprintf('  |       %2d/%2d   %8s %8s %8s %8s %6s', nChanged(k, j), nUsed(k, j), '-', '-', '-', '-', '-');
            end
        end
        if ~isempty(scn), fprintf('  | %+7.2f %+7.2f %+7.2f %+7.2f', intended(k, 1:4)); end
        fprintf('\n');
    end
    if haveRate && any(any(isfinite(fitRate(:, 1, :))))
        fprintf('   the rates on the same epochs:\n');
        fprintf('   %10s', 't [s]');
        for j = 1:nSys, fprintf('  | sys %d %8s %8s %8s %8s %6s', sysIds(j), 'vN', 'vE', 'vD', 'drift', 'rms'); end
        if ~isempty(scn), fprintf('  | intended vN vE vD drift [m/s]'); end
        fprintf('\n');
        for k = pick
            fprintf('   %10.2f', tKf(k));
            for j = 1:nSys
                if isfinite(fitRate(k, 1, j)), fprintf('  |       %+8.3f %+8.3f %+8.3f %+8.3f %6.3f', fitRate(k, 1:5, j));
                else, fprintf('  |       %8s %8s %8s %8s %6s', '-', '-', '-', '-', '-'); end
            end
            if ~isempty(scn), fprintf('  | %+7.3f %+7.3f %+7.3f %+7.3f', intended(k, 5:8)); end
            fprintf('\n');
        end
    end
else
    fprintf('\n   no slot of any system changed on any matched epoch: the two logs carry the same ranges, nothing was injected\n');
end
%% against the scenario
if ~isempty(scn)
    fprintf('\n--- against the scenario -----------------------------------------------------------------------------------------\n');
    fprintf('   the scenario intends an offset on %d epochs, %.2f to %.2f s\n', nnz(somethingIntended), ...
        tKf(find(somethingIntended, 1)), tKf(find(somethingIntended, 1, 'last')));
    for j = 1:nSys
        nMissing(j) = nnz(somethingIntended & matched & ~changedEp(:, j) & nUsed(:, j) > 0);
        nOutside(j) = nnz(~somethingIntended & changedEp(:, j));
        both = somethingIntended & changedEp(:, j) & isfinite(fit(:, 1, j));
        if any(both)
            dv = fit(both, 1:4, j) - intended(both, 1:4);
            fprintf('   system %d: read minus intended on %d epochs, medians N %+.2f E %+.2f D %+.2f clock %+.2f m, largest |difference| %.2f m', ...
                sysIds(j), nnz(both), median(dv, 1), max(abs(dv(:))));
            bothR = both & isfinite(fitRate(:, 1, j));
            if any(bothR)
                dr = fitRate(bothR, 1:4, j) - intended(bothR, 5:8);
                fprintf('; rates medians vN %+.3f vE %+.3f vD %+.3f drift %+.3f m/s, largest %.3f m/s', median(dr, 1), max(abs(dr(:))));
            end
            fprintf('\n');
            if max(abs(dv(:))) < 0.05
                fprintf('   VERDICT: system %d carries the injection as intended\n', sysIds(j));
            else
                fprintf('   VERDICT: system %d carries an offset other than the intended one: check the frame of latLon (radians), the units of satPos and\n', sysIds(j));
                fprintf('            posTrueEcef (metres), and whether the host changes the ranges after the injection\n');
            end
        end
        if nMissing(j) > 0
            fprintf('   system %d: nothing changed on %d of the epochs with an offset intended (first at %.2f s)', sysIds(j), nMissing(j), ...
                tKf(find(somethingIntended & matched & ~changedEp(:, j) & nUsed(:, j) > 0, 1)));
            if ~any(both), fprintf(': no injection on this system'); end
            fprintf('\n');
        end
        if nOutside(j) > 0
            fprintf('   system %d: a slot changed on %d epochs where nothing was intended (first at %.2f s): the logs are not twins there, or the\n', ...
                sysIds(j), nOutside(j), tKf(find(~somethingIntended & changedEp(:, j), 1)));
            fprintf('            host''s ranges depend on its own state\n');
        end
    end
end
res = struct('t', tKf, 'kClean', kClean, 'sysIds', sysIds, 'nUsed', nUsed, 'nChanged', nChanged, 'fit', fit, 'fitRate', fitRate, ...
    'intended', intended, 'tFirst', tFirst, 'tLast', tLast, 'medFit', medFit, 'medRate', medRate, 'nMissing', nMissing, ...
    'nOutside', nOutside, 'shift', shift);
end

%------------------------------------------------------------------------
function [f] = solve(A, r, C)
% least squares of r on A = [-u, 1]; f = [NED offset, fourth term, rms, rows]
x = A \ r;
e = A * x - r;
f = [(C * x(1:3))', x(4), sqrt(mean(e.^2)), numel(r)];
end

function [refK] = refAtEpochs(tRef, posRef, tKf)
% the reference position in ECEF at each epoch: one row for a stationary run, else interpolated on tRef
K = numel(tKf);
if size(posRef, 2) ~= 3 || max(abs(posRef(:, 1))) > pi / 2 || max(abs(posRef(:, 2))) > pi || max(abs(posRef(:, 3))) > 2e4
    error('spoofInjectCheck:posRef', 'posRef must be the reference [lat lon h] in radians, radians and metres; got %d x %d', size(posRef, 1), size(posRef, 2));
end
e = llaToEcef(posRef);
if size(posRef, 1) == 1
    refK = repmat(e, K, 1);
else
    if numel(tRef) ~= size(posRef, 1)
        error('spoofInjectCheck:tRef', 'tRef has %d entries, posRef %d rows', numel(tRef), size(posRef, 1));
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
