%% test_dopplerCheck.m
% Synthetic spfMeas records (the struct STRUCT_SPF.setKfMeas builds, padded to the host layout)
% for a healthy Doppler, an inflated R, a self-referential range rate, an opposite-sign and a
% wrong-frame velocity partial. Checks the figures and the verdict lines of tools/dopplerCheck.
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'tools'));
randn('seed', 11);  rand('seed', 11);

K = 300;  nSat = 8;  n = 60;  mMax = 51;
el = 0.2 + 0.9 * rand(nSat, 1);  az = 2 * pi * rand(nSat, 1);
u  = [cos(el) .* cos(az), cos(el) .* sin(az), -sin(el)];           % unit line of sight, NED (up is -D)
pVar = zeros(n, 1);  pVar(1:3) = 4^2;  pVar(4:6) = 0.15^2;  pVar(16) = 10^2;  pVar(17) = 0.01^2;  pVar(pVar == 0) = 1e-6;
Pm = diag(pVar);                                                   % prior covariance

function [log] = makeLog(K, nSat, n, mMax, u, Pm, sigRr, rrMode, pressure)
    % rrMode: 'ok' | 'self' | 'oppSign' | 'frame'
    m = 2 * nSat + pressure;
    log = cell(K, 1);
    for k = 1:K
        H = zeros(mMax, n);  R = zeros(mMax, mMax);
        for j = 1:nSat
            H(j, 1:3) = -u(j, :);          H(j, 16) = 1;          R(j, j) = 3^2;
            v = -u(j, :);
            if strcmp(rrMode, 'oppSign'), v = u(j, :); end
            if strcmp(rrMode, 'frame'),   v = -[u(j, 2), -u(j, 1), u(j, 3)]; end      % rotated 90 deg about D
            H(nSat + j, 4:6) = v;          H(nSat + j, 17) = 1;   R(nSat + j, nSat + j) = sigRr^2;
        end
        if pressure, H(m, 3) = 1;  R(m, m) = 2^2; end
        S = H * Pm * H' + R;  S = (S + S') / 2;
        y = zeros(mMax, 1);
        Sm = S(1:m, 1:m);
        y(1:m) = chol(Sm)' * randn(m, 1);
        if strcmp(rrMode, 'self'), y(nSat + (1:nSat)) = 0.001 * randn(nSat, 1); end
        xm = zeros(n, 1);
        Kg = Pm * H(1:m, :)' / Sm;
        xp = xm + Kg * y(1:m);
        numMeas = uint8(m);
        if mod(k, 50) == 0, numMeas = uint8(0); end                     % an epoch without GNSS
        log{k} = struct('innovation', y, 'innovationCov', S, 'obsMatrix', H, 'measNoiseCov', R, ...
            'numMeas', numMeas, 'numMeasClamped', false, 'priorState', xm, 'postState', xp, 'postCov', Pm);
    end
end

ok = true;
% A healthy: R = 0.05 m/s on the range rate, a pressure row, cell input, epoch times
logA = makeLog(K, nSat, n, mMax, u, Pm, 0.05, 'ok', 1);
txtA = evalc('resA = dopplerCheck(logA, 0.5 * (1:K));');
fprintf('%s', txtA);
sA = resA.summary;
hphRr = 0.15^2 + 0.01^2;
chk = struct();
chk.rowsA   = sA.nPr == nSat + 1 && sA.nRr == nSat && sA.nPair == nSat && sA.epochsRr == K - 6;
chk.sigmaA  = abs(sA.sR(2) - 0.05) < 1e-9 && abs(sA.sR(1) - 3) < 1e-9 && abs(sA.sHph(2) - sqrt(hphRr)) < 1e-6;
chk.gainA   = abs(sA.gain(2) - hphRr / (hphRr + 0.05^2)) < 1e-6 && abs(sA.gain(1) - 116 / 125) < 1e-6;    % position 16 + clock 100 over + R 9
chk.nisA    = abs(sA.nis(2) - 0.92) < 0.12 && abs(sA.nis(1) - 0.93) < 0.12;     % median of chi2(m)/m, m = 8 and 9
chk.innA    = abs(sA.innRms(2) - sqrt(hphRr + 0.05^2)) < 0.01;
chk.pairA   = sA.pairDot > 0.99 && abs(sA.pairDot - 1) < 1e-9;
chk.dvA     = sA.dvRms > 0.05 && sA.dvRms < 0.3;
chk.timeA   = abs(resA.t(2) - 1) < 1e-12 && numel(resA.t) == K;
chk.quietA  = isempty(strfind(txtA, 'VERDICT'));
% B inflated R: 1 m/s on the range rate
logB = makeLog(K, nSat, n, mMax, u, Pm, 1.0, 'ok', 0);
txtB = evalc('resB = dopplerCheck([logB{:}]);');                   % struct array input
chk.gainB   = resB.summary.gain(2) < 0.03 && ~isempty(strfind(txtB, 'switched off')) && abs(resB.summary.nis(2) - 0.92) < 0.12;
% C self-referential range rate: innovations near zero
logC = makeLog(K, nSat, n, mMax, u, Pm, 0.05, 'self', 0);
txtC = evalc('resC = dopplerCheck(logC);');
chk.selfC   = resC.summary.innRms(2) < 0.01 && ~isempty(strfind(txtC, 'self-referential'));
% D opposite sign on the velocity partial
logD = makeLog(K, nSat, n, mMax, u, Pm, 0.05, 'oppSign', 0);
txtD = evalc('resD = dopplerCheck(logD);');
chk.signD   = resD.summary.pairDot < -0.99 && ~isempty(strfind(txtD, 'OPPOSITE sign'));
% E velocity partial in another frame
logE = makeLog(K, nSat, n, mMax, u, Pm, 0.05, 'frame', 0);
txtE = evalc('resE = dopplerCheck(logE);');
chk.frameE  = abs(resE.summary.pairDot) < 0.99 && ~isempty(strfind(txtE, 'not in the frame'));
names = fieldnames(chk);
for i = 1:numel(names)
    if ~chk.(names{i}), fprintf('   FAILED sub-check: %s\n', names{i}); ok = false; end
end
fprintf('B gain %.3f, C innovation rms %.4f m/s, D pair dot %+.3f, E pair dot %+.3f\n', ...
    resB.summary.gain(2), resC.summary.innRms(2), resD.summary.pairDot, resE.summary.pairDot);
if ok, fprintf('=== test_dopplerCheck PASS ===\n'); else, fprintf('=== test_dopplerCheck FAIL ===\n'); end
