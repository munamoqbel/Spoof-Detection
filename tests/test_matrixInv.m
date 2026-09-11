%% test_matrixInv.m
% Contract check of matrixInv ([inv, invInvalid] = matrixInv(A)) as used
% by SPF_cpiMonitor and SPF_revalidation. Run it with the SIMULATION's own
% matrixInv on the path (the repo file is only a stand-in).
%
% The gate always calls matrixInv with a full MAX_MEAS x MAX_MEAS matrix
% whose rows/cols beyond numMeas are zero EXCEPT the diagonal, which is
% filled with the largest live variance: A = blkdiag(S, p*I). That matrix
% is non-singular, so a function that flags singular inputs is fine.
%
% Required (PASS/FAIL):
%   1. padded input: not flagged, the live block is inverted, the off-block
%      coupling is zero, and products with zero-padded vectors equal the
%      live-block results
%   2. full-size input (numMeas = MAX_MEAS) and a single live row
%   3. accuracy on a realistic S (code + carrier rows, 5 orders of
%      magnitude between variances): A*inv(A) = I on the live block
%   4. non-finite input is flagged (the gate then ignores the result)
%   5. a rank-deficient live block (duplicated row, zero noise) is either
%      flagged (gate drops the epoch) or returns a valid finite pseudo-inverse
% Informational (printed, no verdict):
%   6. whether an indefinite matrix is flagged (the gate also guards it
%      with the live-row variance check and q >= 0 in re-validation)
%   7. time per call (the CPI calls it up to 3 x WINDOW_LENGTH times per epoch)
% Run from the repo root (MATLAB, or Octave with tools/octave_shim).

mMax = double(CST_spfParam.MAX_MEAS);
n    = CST_gnssHybrid.NO_STATES;
rng(11);
% exactly what the gate passes: zero padding, padding diagonal = largest live variance
padAsGate = @(S, mMax) blkdiag(S, max(max(diag(S)), 1.0) * eye(mMax - size(S, 1)));
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

padded = padAsGate(S, mMax);                                       % exactly what the gate passes

%% 1. padding
[Ainv, bad1] = matrixInv(padded);
offBlockZero = all(all(abs(Ainv(m+1:end, 1:m)) < 1e-12)) && all(all(abs(Ainv(1:m, m+1:end)) < 1e-12));
liveErr  = norm(Ainv(1:m, 1:m) * S - eye(m)) / norm(eye(m));
g = zeros(mMax, 1); g(1:m) = randn(m, 1);
f = zeros(mMax, 1); f(1:m) = H(:, 3);
qFull = f' * (Ainv * g);  qLive = H(:, 3)' * (S \ g(1:m));
ok1 = ~bad1 && offBlockZero && liveErr < 1e-8 && abs(qFull - qLive) < 1e-9 * max(1, abs(qLive));
fprintf('  1. padded input: invalid = %d, off-block zero = %d, |inv*S - I| = %.1e, f''S^-1g full vs live diff = %.1e -> %s\n', ...
    bad1, offBlockZero, liveErr, abs(qFull - qLive), pf(ok1));

%% 2. full-size and single row
Afull = randn(mMax); Afull = Afull * Afull' + mMax * eye(mMax);
[AinvF, bad2a] = matrixInv(Afull);
errF = norm(AinvF * Afull - eye(mMax)) / norm(eye(mMax));
one = padAsGate(4.0, mMax);
[Ainv1, bad2b] = matrixInv(one);
ok2 = ~bad2a && errF < 1e-8 && ~bad2b && abs(Ainv1(1, 1) - 0.25) < 1e-12 && all(abs(Ainv1(2:end, 1)) < 1e-12);
fprintf('  2. full size (%d rows): invalid = %d, |inv*A - I| = %.1e | single row: invalid = %d, inv(4) = %.4f -> %s\n', ...
    mMax, bad2a, errF, bad2b, Ainv1(1, 1), pf(ok2));

%% 3. accuracy on the realistic S (carrier vs code: 5 orders of magnitude)
condS = cond(S);
ok3 = ~bad1 && liveErr < 1e-8;
fprintf('  3. realistic S: cond = %.1e, |inv*S - I| = %.1e (need < 1e-8) -> %s\n', condS, liveErr, pf(ok3));

%% 4. non-finite input (flag required; the gate ignores the result when flagged)
bad = padded; bad(2, 3) = NaN; bad(3, 2) = NaN;
[~, badN] = matrixInv(bad);
bad = padded; bad(1, 1) = Inf;
[~, badI] = matrixInv(bad);
ok4 = badN && badI;
fprintf('  4. NaN input: invalid = %d | Inf input: invalid = %d -> %s\n', badN, badI, pf(ok4));

%% 5. rank-deficient live block (duplicated row with zero noise)
Hd = H; Hd(4, :) = Hd(3, :);
Rd = R; Rd(3, 3) = 0; Rd(4, 4) = 0;
Sd = Hd * P * Hd' + Rd; Sd = (Sd + Sd') / 2;
[AinvD, badD] = matrixInv(padAsGate(Sd, mMax));
if badD
    ok5 = true;                                                     % flagged: gate drops the epoch
    note = 'flagged (gate drops the epoch)';
else
    finiteD = all(isfinite(AinvD(:)));
    pinvErr = norm(Sd * AinvD(1:m, 1:m) * Sd - Sd) / norm(Sd);      % pseudo-inverse property A*A+*A = A
    ok5 = finiteD && pinvErr < 1e-8;
    note = sprintf('not flagged: finite = %d, |A A+ A - A|/|A| = %.1e', finiteD, pinvErr);
end
fprintf('  5. rank-deficient S (rank %d of %d): invalid = %d, %s -> %s\n', rank(Sd), m, badD, note, pf(ok5));

%% 6. informational: indefinite matrix
Sind = S; Sind(1, 1) = -Sind(1, 1);
[~, badInd] = matrixInv(padAsGate(Sind, mMax));
fprintf('  6. (info) indefinite S: invalid = %d (the gate also rejects it through the live-row\n', badInd);
fprintf('     variance check and q >= 0 in re-validation)\n');

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
