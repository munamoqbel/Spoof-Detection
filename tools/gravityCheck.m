function [res] = gravityCheck(run, gHost)
%GRAVITYCHECK  The normal gravity at a run's site, and the vertical acceleration error a host carries when it uses
%   another value, in mg.
%   res = gravityCheck(run)            run: the runCompare struct; the site is the mean latitude and height of
%                                      truePos (or navPos), and the span of the run gives the spread
%   res = gravityCheck(run, gHost)     gHost: the value the host's mechanisation uses at the site (m/s^2), read from
%                                      its gravity function or constant; its difference is printed as well
%   res = gravityCheck([lat h], gHost) the site given directly (rad, m)
%
%   Normal gravity on the WGS84 ellipsoid (Somigliana) with the free-air height terms:
%     g(phi, h) = 9.7803253359 (1 + 0.00193185265241 sin^2 phi) / sqrt(1 - 0.00669437999014 sin^2 phi)
%                 - (3.0877e-6 - 4.4e-9 sin^2 phi) h + 7.2e-13 h^2                               (m/s^2)
%   A host that carries a constant (9.80665, 9.81, 9.8), the ellipsoid value without the height terms, or the
%   Somigliana closed form with its constant built as gp/ge - 1 (the semi-major axis in both places of
%   k = b gp / (a ge) - 1, which gives 0.0053 in place of 0.0019), differs from it by an amount the mechanisation
%   integrates as a vertical velocity error between the updates. feedbackCheck
%   measures that amount as the constant of its vertical budget (the acceleration error left with the bias right);
%   the table here gives, for each candidate, the Down acceleration error it carries at the site, positive when the
%   host's gravity is larger than the normal one (the nav then descends too fast). Match the candidate to the
%   measured constant.
%   res: lat, h, gNormal, the candidates with their error in mg, and the spread over the run.
%   MATLAB desktop tool, not codegen. Runs under Octave.

MG = 9.80665e-3;
if isstruct(run)
    if isfield(run, 'truePos') && ~isempty(run.truePos), P = run.truePos; src = 'truePos';
    elseif isfield(run, 'navPos') && ~isempty(run.navPos), P = run.navPos; src = 'navPos';
    else, error('gravityCheck: the run has neither truePos nor navPos');
    end
    ok = all(isfinite(P(:, [1 3])), 2);
    lat = mean(P(ok, 1));  h = mean(P(ok, 3));
    latSpan = [min(P(ok, 1)), max(P(ok, 1))];  hSpan = [min(P(ok, 3)), max(P(ok, 3))];
else
    lat = run(1);  h = run(2);  src = 'given';
    latSpan = [lat lat];  hSpan = [h h];
end
gN   = normalGravity(lat, h);
gSea = normalGravity(lat, 0);
gLo  = normalGravity(latSpan(1), hSpan(2));                   % the smallest normal gravity the run meets
gHi  = normalGravity(latSpan(2), hSpan(1));                   % the largest
gBad = somiglianaBadK(lat, h);                                % the closed form with k = gp/ge - 1 (a/a in place of b/a)
cand = {'9.80665 (standard gravity, a constant)', 9.80665; '9.81 (a constant)', 9.81; '9.8 (a constant)', 9.8; ...
        'ellipsoid value at the site without the height terms', gSea; '9.780327 (the equator, a constant)', 9.7803253359; ...
        'Somigliana with k = gp/ge - 1 (a/a in place of b/a in k)', gBad};
errMg = (cell2mat(cand(:, 2)) - gN) / MG;

fprintf('gravityCheck: site from %s, latitude %.4f deg, height %.0f m (the run spans %.4f to %.4f deg, %.0f to %.0f m)\n', ...
    src, lat * 180 / pi, h, latSpan * 180 / pi, hSpan);
fprintf('   normal gravity (WGS84, with the height terms) at the site: %.6f m/s^2; over the run it ranges %.6f to %.6f, a spread of %.2f mg\n', ...
    gN, gLo, gHi, (gHi - gLo) / MG);
fprintf('   a host that uses                                            carries a Down acceleration error of\n');
for i = 1:size(cand, 1)
    fprintf('      %-54s %+7.2f mg\n', cand{i, 1}, errMg(i));
end
if nargin >= 2 && ~isempty(gHost)
    eH = (gHost - gN) / MG;
    fprintf('      the host''s own value %.6f m/s^2                       %+7.2f mg\n', gHost, eH);
else
    eH = NaN;
end
fprintf('   (positive: the host''s gravity is larger than the normal one and the nav descends too fast; compare with the constant\n');
fprintf('   feedbackCheck prints as the acceleration error left with the bias right)\n');
res = struct('lat', lat, 'h', h, 'gNormal', gN, 'gSea', gSea, 'candidates', {cand(:, 1)}, 'candidateErrMg', errMg, ...
    'gHost', [], 'hostErrMg', eH, 'spreadMg', (gHi - gLo) / MG);
if nargin >= 2, res.gHost = gHost; end
end

%------------------------------------------------------------------------
function [g] = somiglianaBadK(lat, h)
% the closed form with the constant built from the semi-major axis in both places, k = gp / ge - 1
s2 = sin(lat)^2;
k = 9.8321849378 / 9.7803253359 - 1;
g = 9.7803253359 * (1 + k * s2) / sqrt(1 - 0.00669437999014 * s2) - (3.0877e-6 - 4.4e-9 * s2) * h + 7.2e-13 * h^2;
end

function [g] = normalGravity(lat, h)
s2 = sin(lat)^2;
g = 9.7803253359 * (1 + 0.00193185265241 * s2) / sqrt(1 - 0.00669437999014 * s2) - (3.0877e-6 - 4.4e-9 * s2) * h + 7.2e-13 * h^2;
end
