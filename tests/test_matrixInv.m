%% test_matrixInv.m
% Contract check of matrixInv ([inv, invInvalid] = matrixInv(A)) as used
% by SPF_cpiMonitor and SPF_revalidation. Run it with the SIMULATION's own
% matrixInv on the path (the repo file is only a stand-in).
%
% Required (PASS/FAIL):
%   1. zero-padded input: the live block is inverted, the padding stays
%      zero, and products with zero-padded vectors equal the live results
%   2. full-size input (numMeas = MAX_MEAS) and a single live row
%   3. accuracy on a realistic S (code + carrier rows, 5 orders of
%      magnitude between variances): A*inv(A) = I on the live block
%   4. non-finite input is flagged (invInvalid) and the result is finite
%   5. a rank-deficient live block (duplicated row, zero noise) returns a
%      finite result and, if not flagged, a valid pseudo-inverse
% Informational (printed, no verdict):
%   6. an indefinite matrix is NOT detectable by an SVD inverse (the gate
%      guards it with the live-row variance check and q >= 0)
%   7. time per call (the CPI calls it up to 3 x WINDOW_LENGTH times per epoch)
% Run from the repo root (MATLAB, or Octave with tools/octave_shim).

mMax = double(CST_spfParam.MAX_MEAS);
n    = CST_gnssHybrid.NO_STATES;
rng(11);
fprintf('matrixInv contract (MAX_MEAS = %d) ...\n', mMax);

% realistic S: 6 satellites x (code, carrier) rows, sigma 0.36 m / 0.003 m
nSv = 6; m = 2 * nSv;
los = randn(nSv, 3); los = los ./ vecnorm(los, 2, 2);
H = zeros(m, n);
for sv = 1:nSv
    H(2*sv-1, 1:3) = -los(sv, :); H(2*sv, 1:3) = -los(sv, :);
end
R = diag(repmat([0.36^2 0.003^2], 1, nSv));
P = 1e-2 * eye(n);
S = H * P * H' + R; S = (S + S') / 2;                             % live block, PD

padded = zeros(mMax); padded(1:m, 1:m) = S;                        % zero-padded as the gate builds it

%% 1. padding
[Ainv, bad1] = matrixInv(padded);
padZero  = all(all(Ainv(m+1:end, :) == 0)) && all(all(Ainv(:, m+1:end) == 0));
liveErr  = norm(Ainv(1:m, 1:m) * S - eye(m)) / norm(eye(m));
g = zeros(mMax, 1); g(1:m) = randn(m, 1);
f = zeros(mMax, 1); f(1:m) = H(:, 3);
qFull = f' * (Ainv * g);  qLive = H(:, 3)' * (S \ g(1:m));
ok1 = ~bad1 && padZero && liveErr < 1e-8 && abs(qFull - qLive) < 1e-9 * max(1, abs(qLive));
fprintf('  1. zero padding: invalid = %d, padding zero = %d, |inv*S - I| = %.1e, f''S^-1g full vs live diff = %.1e -> %s\n', ...
    bad1, padZero, liveErr, abs(qFull - qLive), pf(ok1));

%% 2. full-size and single row
Afull = randn(mMax); Afull = Afull * Afull' + mMax * eye(mMax);
[AinvF, bad2a] = matrixInv(Afull);
errF = norm(AinvF * Afull - eye(mMax)) / norm(eye(mMax));
one = zeros(mMax); one(1, 1) = 4.0;
[Ainv1, bad2b] = matrixInv(one);
ok2 = ~bad2a && errF < 1e-8 && ~bad2b && abs(Ainv1(1, 1) - 0.25) < 1e-12 && all(all(Ainv1(2:end, :) == 0));
fprintf('  2. full size (%d rows): invalid = %d, |inv*A - I| = %.1e | single row: inv(4) = %.4f -> %s\n', ...
    mMax, bad2a, errF, Ainv1(1, 1), pf(ok2));

%% 3. accuracy on the realistic S (carrier vs code: 5 orders of magnitude)
condS = cond(S);
ok3 = liveErr < 1e-8;
fprintf('  3. realistic S: cond = %.1e, |inv*S - I| = %.1e (need < 1e-8) -> %s\n', condS, liveErr, pf(ok3));

%% 4. non-finite input
bad = padded; bad(2, 3) = NaN; bad(3, 2) = NaN;
[AinvN, badN] = matrixInv(bad);
bad = padded; bad(1, 1) = Inf;
[AinvI, badI] = matrixInv(bad);
ok4 = badN && badI && all(isfinite(AinvN(:))) && all(isfinite(AinvI(:)));
fprintf('  4. NaN input: invalid = %d, result finite = %d | Inf input: invalid = %d, result finite = %d -> %s\n', ...
    badN, all(isfinite(AinvN(:))), badI, all(isfinite(AinvI(:))), pf(ok4));

%% 5. rank-deficient live block (duplicated row with zero noise)
Hd = H; Hd(4, :) = Hd(3, :);
Rd = R; Rd(3, 3) = 0; Rd(4, 4) = 0;
Sd = Hd * P * Hd' + Rd; Sd = (Sd + Sd') / 2;
paddedD = zeros(mMax); paddedD(1:m, 1:m) = Sd;
[AinvD, badD] = matrixInv(paddedD);
finiteD = all(isfinite(AinvD(:)));
if badD
    ok5 = finiteD;                                                  % flagged: fine, gate drops the epoch
    note = 'flagged (gate drops the epoch)';
else
    pinvErr = norm(Sd * AinvD(1:m, 1:m) * Sd - Sd) / norm(Sd);      % pseudo-inverse property A*A+*A = A
    ok5 = finiteD && pinvErr < 1e-8;
    note = sprintf('pseudo-inverse, |A A+ A - A|/|A| = %.1e', pinvErr);
end
fprintf('  5. rank-deficient S (rank %d of %d): invalid = %d, finite = %d, %s -> %s\n', ...
    rank(Sd), m, badD, finiteD, note, pf(ok5));

%% 6. informational: indefinite matrix
Sind = S; Sind(1, 1) = -Sind(1, 1);
paddedI = zeros(mMax); paddedI(1:m, 1:m) = Sind;
[~, badInd] = matrixInv(paddedI);
fprintf('  6. (info) indefinite S: invalid = %d. An SVD inverse cannot detect this; the gate rejects it\n', badInd);
fprintf('     through the live-row variance check and q >= 0 in re-validation.\n');

%% 7. informational: timing
nCall = 200; tic;
for k = 1:nCall, [~, ~] = matrixInv(padded); end
tCall = toc / nCall;
fprintf('  7. (info) %.3f ms per call; CPI worst case %d calls per epoch = %.1f ms\n', ...
    1e3 * tCall, 3 * double(CST_spfParam.WINDOW_LENGTH), 1e3 * tCall * 3 * double(CST_spfParam.WINDOW_LENGTH));

if ok1 && ok2 && ok3 && ok4 && ok5
    fprintf('=== test_matrixInv PASS ===\n');
else
    fprintf('=== test_matrixInv FAIL ===\n');
end

function s = pf(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end
