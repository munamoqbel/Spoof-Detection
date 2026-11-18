%% test_updateStepCheck.m
% Synthetic spfMeas records from three hosts with the same rows, P- and innovations: one whose
% state takes the whole Kalman step, one whose gain was built with an R scaled per row type
% (range-rate rows x 10, range rows x 4) while the record carries the unscaled R, and one that
% puts a factor per state between K y and the state. Checks the slopes, the row-type
% coefficients and the verdict of tools/updateStepCheck.
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'tools'));
randn('seed', 7);  rand('seed', 7);

n = 22;  nSat = 14;  mMax = 31;  K = 240;
sig = [0.6 0.6 1.0, 0.013 0.013 0.025, 0.4e-3 0.4e-3 0.34e-3, 5e-8 5e-8 5e-8, 1e-3 1e-3 1e-3, 3 0.02 3 0.02 1 3 0.02];
Pm = diag(sig.^2);
Pm(1, 4) = 0.5 * sig(1) * sig(4);  Pm(4, 1) = Pm(1, 4);           % position-velocity cross term
Pm(4, 9) = 0.6 * sig(4) * sig(9);  Pm(9, 4) = Pm(4, 9);           % velocity-heading cross term: the heading moves through the rate rows
Pm(5, 9) = -0.4 * sig(5) * sig(9); Pm(9, 5) = Pm(5, 9);
el = 0.15 + 1.2 * rand(nSat, 1);  az = 2 * pi * rand(nSat, 1);
u  = [cos(el) .* cos(az), cos(el) .* sin(az), -sin(el)];
H = zeros(mMax, n);  Rd = zeros(mMax, 1);
for j = 1:nSat
    sysCol = 16 + 2 * (j > nSat / 2);                             % GPS clock 16, GLONASS clock 18
    H(j, 1:3) = -u(j, :);               H(j, sysCol) = 1;          Rd(j) = 5^2;
    H(nSat + j, 4:6) = -u(j, :);        H(nSat + j, sysCol + 1) = 1;  Rd(nSat + j) = 0.1^2;
end
H(2 * nSat + 1, 3) = 1;  Rd(2 * nSat + 1) = 1.5^2;
m = 2 * nSat + 1;
isRate = (1:m)' > nSat & (1:m)' <= 2 * nSat;
isRange = (1:m)' <= nSat;
R = diag(Rd(1:m));
Hm = H(1:m, :);

function [log] = makeHost(mode, K, n, mMax, m, Hm, R, Pm, isRate, isRange)
    % mode: 'whole' | 'R' | 'perState'
    cR = ones(m, 1);
    if strcmp(mode, 'R'), cR(isRate) = 10;  cR(isRange) = 4; end
    Rused = diag(cR) * R;
    F = ones(n, 1);
    if strcmp(mode, 'perState'), F(1:3) = 0.25;  F(4:6) = 0.10;  F(7:9) = 0.14; end
    S = Hm * Pm * Hm' + Rused;  S = (S + S') / 2;
    Kg = Pm * Hm' / S;
    I = eye(n);
    Pp = (I - Kg * Hm) * Pm * (I - Kg * Hm)' + Kg * Rused * Kg';     % Joseph form with the R the gain used
    Pp = (Pp + Pp') / 2;
    Sl = Hm * Pm * Hm' + R;  Sl = (Sl + Sl') / 2;                  % what the record carries: the host's prior P with the handed-over R
    log = cell(K, 1);
    xm = zeros(n, 1);
    for k = 1:K
        err = zeros(n, 1);                                          % a slowly varying truth error the rows see
        err(1:3) = [3 * cos(k / 15); -2 * sin(k / 20); 1];
        err(4:6) = [0.3 * sin(k / 25); 0.4 * cos(k / 30); 0.05];
        y = zeros(mMax, 1);
        y(1:m) = Hm * err + sqrt(diag(R)) .* randn(m, 1);
        dx = Kg * y(1:m);
        xp = xm + F .* dx;
        numMeas = uint8(m);
        if mod(k, 60) == 0, numMeas = uint8(0); end
        Hpad = zeros(mMax, n);  Hpad(1:m, :) = Hm;
        Rpad = zeros(mMax);  Rpad(1:m, 1:m) = R;
        Spad = zeros(mMax);  Spad(1:m, 1:m) = Sl;
        log{k} = struct('innovation', y, 'innovationCov', Spad, 'obsMatrix', Hpad, 'measNoiseCov', Rpad, ...
            'measNoiseCovHost', Rpad, 'innovationCovHost', Spad, 'numMeas', numMeas, 'numMeasClamped', false, ...
            'priorState', xm, 'postState', xp, 'postCov', Pp);
        xm = 0.9 * xp;                                              % the split: a tenth leaves the state per epoch
    end
end

ok = true;  chk = struct();
tK = 0.5 * (1:K)';
% A whole update
logA = makeHost('whole', K, n, mMax, m, Hm, R, Pm, isRate, isRange);
txtA = evalc('resA = updateStepCheck(logA, tK);');
fprintf('%s', txtA);
gA = resA.signal;
chk.verdictA = strcmp(resA.verdict, 'whole') && ~isempty(strfind(txtA, 'receives the Kalman step'));
chk.slopeA   = all(abs(resA.slope(gA) - 1) < 1e-6) && all(resA.fit(gA) > 1 - 1e-6);
chk.groupsA  = resA.signal(1) && resA.signal(2) && resA.signal(4) && resA.signal(7) && resA.signal(8) && ...
    strcmp(resA.groups{7, 1}, 'state 16') && size(resA.groups, 1) == 6 + 7;         % position, velocity, heading (cross term), clock bias and drift
chk.validA   = nnz(resA.valid) == K - 4 && numel(resA.t) == K;
% B gain built with a scaled R
logB = makeHost('R', K, n, mMax, m, Hm, R, Pm, isRate, isRange);
txtB = evalc('resB = updateStepCheck([logB{:}], tK, [0 1e9]);');          % struct array input, window
fprintf('%s', txtB);
gB = resB.signal;
aRate = resB.aRow(gB, 2);  aRange = resB.aRow(gB, 1);
chk.verdictB = strcmp(resB.verdict, 'R') && ~isempty(strfind(txtB, 'range-rate rows x 10.0')) && ...
    abs(resB.cRow(1) - 4) < 0.15 && abs(resB.cRow(2) - 10) < 0.05 && abs(resB.cRow(3) - 1) < 0.02;
chk.coefB    = all(abs(aRate(isfinite(aRate)) - 0.1) < 5e-3) && all(abs(aRange(isfinite(aRange)) - 0.25) < 1e-2) && ...
    nnz(isfinite(aRate)) >= 3 && nnz(isfinite(aRange)) >= 2 && all(resB.fitRow(gB) > 0.95);
chk.slopeB   = abs(resB.slope(2) - 0.1) < 0.01 && abs(resB.slope(4) - 0.1) < 0.01 && resB.slope(1) > 0.2 && resB.slope(1) < 0.7 && ...
    abs(resB.slope(7) - 0.25) < 0.01 && abs(resB.slope(8) - 0.1) < 0.01;            % clock bias 16 a quarter, clock drift 17 a tenth
% C a factor per state
logC = makeHost('perState', K, n, mMax, m, Hm, R, Pm, isRate, isRange);
txtC = evalc('resC = updateStepCheck(logC, tK);');
chk.verdictC = strcmp(resC.verdict, 'per state') && ~isempty(strfind(txtC, 'factor per state'));
chk.slopeC   = abs(resC.slope(1) - 0.25) < 1e-6 && abs(resC.slope(2) - 0.10) < 1e-6 && abs(resC.slope(4) - 0.14) < 1e-6 && ...
    abs(resC.slope(7) - 1) < 1e-6 && abs(resC.slope(8) - 1) < 1e-6 && all(resC.fit(resC.signal) > 1 - 1e-6);
aC = resC.aRow(2, :);  aC = aC(isfinite(aC));
chk.coefC    = all(abs(aC - 0.10) < 1e-3);
% D window with no rows
txtD = evalc('resD = updateStepCheck(logA, tK, [1e6 2e6]);');
chk.emptyD   = strcmp(resD.verdict, 'no signal') && ~isempty(strfind(txtD, 'no epoch'));

names = fieldnames(chk);
for i = 1:numel(names)
    if ~chk.(names{i}), fprintf('   FAILED sub-check: %s\n', names{i}); ok = false; end
end
fprintf('A %s, B %s (c_rate %.1f, c_range %.1f), C %s (slopes %.2f %.2f %.2f)\n', resA.verdict, resB.verdict, ...
    resB.cRow(2), resB.cRow(1), resC.verdict, resC.slope(1), resC.slope(2), resC.slope(4));
if ok
    fprintf('=== test_updateStepCheck PASS ===\n');
else
    fprintf('=== test_updateStepCheck FAIL ===\n');
end
