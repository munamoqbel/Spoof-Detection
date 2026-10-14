%% test_designConstants.m
% tools/designConstants on synthetic spfMeas records: s2 per axis against the direct Eq. 20, the window
% length against solve_N_min called directly, the constants at the common N, a cell input, the arming epoch
% from a P diagonal whose protection level settles, and the unreachable case (code-noise prior) reporting
% the spoofer error it would take. Run from the repo root (Octave with the shim).
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'tools'));
randn('seed', 5);  rand('seed', 5);

K = 120;  nSat = 7;  n = CST_gnssHybrid.NO_STATES;  mMax = double(CST_spfParam.MAX_MEAS);
prm = kujur_params();
nAx = 3;  pFaPlus = prm.P_FA_plus / nAx;  pFaMinus = prm.P_FA_minus / nAx;
sigmaT = 0.1;

function [log, s2ref] = makeLog(K, nSat, n, mMax, sigPos, sigR, sigP)
    pVar = 1e-6 * ones(n, 1);  pVar(1:3) = sigPos^2;  pVar(16) = 8^2;
    Pm = diag(pVar);
    log = cell(K, 1);  s2ref = zeros(K, 3);
    for k = 1:K
        el = 0.15 + 0.9 * rand(nSat, 1);  az = 2 * pi * rand(nSat, 1);
        u  = [cos(el) .* cos(az), cos(el) .* sin(az), -sin(el)];
        m  = nSat + 1;                                          % range rows and a pressure row
        H = zeros(mMax, n);  R = zeros(mMax, mMax);
        for j = 1:nSat, H(j, 1:3) = -u(j, :);  H(j, 16) = 1;  R(j, j) = sigR^2; end
        H(m, 3) = 1;  R(m, m) = sigP^2;
        S = H * Pm * H' + R;  S = (S + S') / 2;
        y = zeros(mMax, 1);  y(1:m) = chol(S(1:m, 1:m))' * randn(m, 1);
        numMeas = uint8(m);
        if k == 7, numMeas = uint8(2); end                      % too few rows: left out of the design
        log{k} = struct('innovation', y, 'innovationCov', S, 'obsMatrix', H, 'measNoiseCov', R, ...
            'numMeas', numMeas, 'numMeasClamped', false, 'priorState', zeros(n, 1), 'postState', zeros(n, 1), 'postCov', Pm);
        for a = 1:3
            f = H(1:m, a);  s2ref(k, a) = f' * (S(1:m, 1:m) \ f);
        end
    end
    s2ref(7, :) = NaN;
end

% A. a centimetre prior with carrier-level rows: the requirement is reachable
[logA, s2A] = makeLog(K, nSat, n, mMax, 0.01, 0.02, 0.02);
nRef = zeros(1, 3);
for a = 1:3
    [nRef(a), ~, ~, ~, ~] = solve_N_min(min(s2A(:, a)), sigmaT, pFaPlus, pFaMinus, prm.P_MD_plus, prm.P_MD_minus);
end
nC = max(nRef) + 1;
sigPosTrace = 2 + 58 * exp(-(0:K-1)' / 15);                    % position sigma 60 m down to 2 m
pDiag = [repmat(sigPosTrace.^2, 1, 3), 1e-4 * ones(K, n - 3)];
kMd = norminv(1 - prm.P_MD_minus);
inside = double(kMd * sigPosTrace < CST_spfParam.ARM_PL_MAX);
w = double(CST_spfParam.ARM_EPOCHS);
armRef = find(conv(inside, ones(w, 1), 'valid') == w, 1);
desA  = designConstants(logA, sigmaT, pDiag);                   % cell input
desA2 = designConstants([logA{:}], sigmaT);                     % struct-array input, no P
okA = all(abs(desA.s2min - min(s2A, [], 1)) < 1e-9 * max(1, max(s2A(:)))) && all(isnan(desA.s2(7, :))) ...
   && all(desA.nAxis == nRef) && desA.N == nC ...
   && abs(desA.tN - gaminv(1 - pFaPlus, nC / 2, 2)) < 1e-9 ...
   && abs(desA.kFa - norminv(1 - (pFaMinus / nC) / 2)) < 1e-9 && abs(desA.kMd - kMd) < 1e-9 ...
   && all(desA.pmdAtN <= desA.pmdAxis + 1e-15) && desA.armEpoch == armRef && isnan(desA2.armEpoch) && desA2.N == nC ...
   && all(desA.nominalBias == double(CST_spfParam.NOMINAL_BIAS)) ...
   && all(abs(desA.cpiOffsetMin - (desA.nominalBias + sqrt(desA.tN / desA.N) ./ sqrt(desA.s2med))) < 1e-9);

% B. a metre-level prior with code rows: not reachable, the spoofer error it would take is reported
[logB, s2B] = makeLog(K, nSat, n, mMax, 3.0, 2.5, 4.0);
desB = designConstants(logB, sigmaT);
sTreq = desB.sigmaTreq(1, 2);
pmdAtReq = compute_PMD_eq38(240, desB.s2min(1), sTreq, pFaPlus);
okB = all(isnan(desB.nAxis)) && isnan(desB.N) && all(isfinite(desB.sigmaTreq(:))) && all(desB.sigmaTreq(:, 1) >= desB.sigmaTreq(:, 2)) ...
   && abs(pmdAtReq - prm.P_MD_plus) < 0.02 * prm.P_MD_plus && abs(desB.kMd - kMd) < 1e-9 ...
   && all(abs(desB.cpiOffsetMin - (desB.nominalBias + sqrt(desB.current.tN / desB.current.N) ./ sqrt(desB.s2med))) < 1e-9);

if okA && okB, fprintf('=== test_designConstants PASS ===\n'); else, fprintf('=== test_designConstants FAIL (A %d, B %d) ===\n', okA, okB); end
