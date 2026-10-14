function [des] = designConstants(spfMeasLog, sigmaT, pDiag)
%DESIGNCONSTANTS  Offline design of the gate's monitor constants from the host's own update records.
%   des = designConstants(spfMeasLog)                 spfMeasLog: struct array or cell of the per-epoch spfMeas
%   des = designConstants(spfMeasLog, sigmaT)         sigmaT: spoofer tracking-error std (m), default 0.10
%   des = designConstants(spfMeasLog, sigmaT, pDiag)  pDiag [K x n] P diagonal per epoch (position states 1:3):
%                                                     also reports from which epoch the protection level allows arming
%
%   This is QUICKSTART step 2 on recorded data. Per monitored axis a (CST_spfParam.MONITORED_AXES), on every
%   epoch with at least REVAL_MIN_MEAS rows:
%       s2(k) = f' S^-1 f,  f = H(:, a)                                   (Eq. 20 of Kujur 2024)
%   s2min(a) is the minimum over the run, the worst sky. Omega = s2min * sigmaT^2 is the spoofer's tracking
%   jitter against the filter's own prior uncertainty along the axis, and the missed-detection requirement
%   is reachable only when Omega is of order one or more: with a prior position sigma of metres, as a loose
%   Q gives, a 0.1 m spoofer is invisible to the CPI at any window length. For a reachable axis solve_N_min
%   gives the window length N(a) with the false-alert allocations of kujur_params split across the axes; the
%   common window is N = max N(a) + 1 (one epoch of margin), and at that N:
%       T_N  = gaminv(1 - P_FA_plus / nAxes, N/2, 2)                       (Eq. 35, CPI_THRESHOLD)
%       k_FA = norminv(1 - (P_FA_minus / nAxes / N) / 2)                   (K_FALSE_ALERT)
%       k_MD = norminv(1 - P_MD_minus)                                     (K_MISSED_DETECTION)
%   For an axis that is not reachable within MAX_WINDOW epochs the tool prints instead the spoofer tracking
%   error sigma_t from which the requirement is met at the current WINDOW_LENGTH and at MAX_WINDOW.
%   The values are printed next to the ones hard-coded in CST_spfParam. Rerun after any change of Q, R, the
%   measurement set or the budgets: S carries the host's P and R, so these constants belong to one tuning.
%   With CST_spfParam.NOMINAL_BIAS set (B metres per axis: the nominal GNSS excursion the host tolerates) the
%   tool states what that costs: the separation test fires beyond B + k_FA sigma_SS, the CPI counts only the
%   axis innovation beyond B, so the sigma_t design applies to a spoofer whose offset exceeds B. In both cases
%   it prints the constant offset the CPI reports within one window on the median sky,
%       d_min = B + sqrt(T_N / N) / sqrt(s2med)                                    (N ((d - B) / sigma_a)^2 = T_N)
%   MATLAB desktop tool, not codegen (gaminv, norminv). Runs under Octave with tools/octave_shim.

MAX_WINDOW = 240;                                            % epochs (2 min at 2 Hz): longer is not a monitor

if iscell(spfMeasLog), spfMeasLog = [spfMeasLog{:}]; end
if nargin < 2 || isempty(sigmaT), sigmaT = 0.10; end
if nargin < 3, pDiag = []; end
K = numel(spfMeasLog);
prm = kujur_params();
axes = double(CST_spfParam.MONITORED_AXES(:))';
nAx = numel(axes);
pFaPlus  = prm.P_FA_plus  / nAx;                           % CPI false-alert allocation per axis
pFaMinus = prm.P_FA_minus / nAx;                           % SS  false-alert allocation per axis (total over the window)
minRows  = double(CST_spfParam.REVAL_MIN_MEAS);
cur = struct('N', double(CST_spfParam.WINDOW_LENGTH), 'tN', CST_spfParam.CPI_THRESHOLD, ...
    'kFa', CST_spfParam.K_FALSE_ALERT, 'kMd', CST_spfParam.K_MISSED_DETECTION);

%% Eq. 20 on every usable epoch
s2 = NaN(K, nAx);
for k = 1:K
    s = spfMeasLog(k);
    m = double(s.numMeas);
    if m >= minRows
        Sm = s.innovationCov(1:m, 1:m);
        H  = s.obsMatrix(1:m, :);
        if all(isfinite(Sm(:))) && all(isfinite(H(:)))
            for a = 1:nAx
                f = H(:, axes(a));
                s2(k, a) = f' * (Sm \ f);
            end
        end
    end
end
used = all(isfinite(s2), 2) & all(s2 > 0, 2);
if nnz(used) < 10
    error('designConstants:rows', 'only %d epochs with %d or more rows and a usable S: nothing to design on', nnz(used), minRows);
end
kUsed = find(used);
[s2min, iMin] = min(s2(used, :), [], 1);
s2med = median(s2(used, :), 1);
fprintf('designConstants: %d epochs, %d usable (>= %d rows, finite S); sigma_t %.2f m; budgets P_FA+ %.1e, P_FA- %.1e, P_MD+ %.1e, P_MD- %.1e over %d axes\n', ...
    K, nnz(used), minRows, sigmaT, prm.P_FA_plus, prm.P_FA_minus, prm.P_MD_plus, prm.P_MD_minus, nAx);

%% window length per axis where the requirement is reachable, else the spoofer error it would take
nAxis = NaN(1, nAx);  pmdAxis = NaN(1, nAx);  sigmaTreq = NaN(nAx, 2);
for a = 1:nAx
    omega = s2min(a) * sigmaT^2;
    reachable = compute_PMD_eq38(MAX_WINDOW, s2min(a), sigmaT, pFaPlus) <= prm.P_MD_plus;
    fprintf('   axis %d (state %d): s2 min %.4g at epoch %d, median %.4g, Omega %.3g', a, axes(a), s2min(a), kUsed(iMin(a)), s2med(a), omega);
    if reachable
        fprintf('\n');
        [nAxis(a), ~, ~, ~, pmdAxis(a)] = solve_N_min(s2min(a), sigmaT, pFaPlus, pFaMinus, prm.P_MD_plus, prm.P_MD_minus);
    else
        sigmaTreq(a, 1) = sigmaTrequired(cur.N, s2min(a), pFaPlus, prm.P_MD_plus);
        sigmaTreq(a, 2) = sigmaTrequired(MAX_WINDOW, s2min(a), pFaPlus, prm.P_MD_plus);
        fprintf(': P_MD+ not reachable within %d epochs at sigma_t %.2f m; it is met from sigma_t %.2f m at N = %d, %.2f m at N = %d\n', ...
            MAX_WINDOW, sigmaT, sigmaTreq(a, 1), cur.N, sigmaTreq(a, 2), MAX_WINDOW);
    end
end

%% common window and the constants at it
N = NaN;  tN = NaN;  kFa = NaN;  pmdAtN = NaN(1, nAx);
kMd = norminv(1 - prm.P_MD_minus);
if any(isfinite(nAxis))
    N   = max(nAxis(isfinite(nAxis))) + 1;                   % one epoch of margin, as in the quickstart
    tN  = gaminv(1 - pFaPlus, N / 2, 2);
    kFa = norminv(1 - (pFaMinus / N) / 2);
    for a = 1:nAx
        pmdAtN(a) = compute_PMD_eq38(N, s2min(a), sigmaT, pFaPlus);
    end
    fprintf('                            from this run      CST_spfParam now\n');
    fprintf('   WINDOW_LENGTH            %8d           %8d\n', N, cur.N);
    fprintf('   CPI_THRESHOLD   T_N      %14.6f     %14.6f\n', tN, cur.tN);
    fprintf('   K_FALSE_ALERT   k_FA     %14.6f     %14.6f\n', kFa, cur.kFa);
    fprintf('   K_MISSED_DETECTION k_MD  %14.6f     %14.6f\n', kMd, cur.kMd);
    fprintf('   PMD achieved at the common N per axis: %s(requirement %.1e; an axis above it is not reachable)\n', sprintf('%.2e ', pmdAtN), prm.P_MD_plus);
    if N ~= cur.N
        fprintf('   NOTE: the window length changes; PROBATION_LENGTH follows it (N + 8 in the harness design) and the\n');
        fprintf('         monitor buffers are sized from WINDOW_LENGTH: update CST_spfParam (root and deploy), rerun runAllTests\n');
    end
else
    fprintf('   no axis reaches P_MD+ at sigma_t %.2f m with this filter: the prior uncertainty along the axes (1/sqrt(s2) =\n', sigmaT);
    fprintf('   %s m) is far above the spoofer jitter. A tighter Q lowers it; otherwise the CPI design must assume the\n', sprintf('%.2f ', 1 ./ sqrt(s2min)));
    fprintf('   larger sigma_t printed above. K_MISSED_DETECTION k_MD %.6f (CST_spfParam %.6f) does not depend on it.\n', kMd, cur.kMd);
end

%% nominal bias allowance and the constant offset the CPI reports within one window
nominalBias = double(CST_spfParam.NOMINAL_BIAS(:))';
bAxis = nominalBias(axes);
nUse = N;  tUse = tN;
if ~isfinite(N), nUse = cur.N;  tUse = cur.tN; end           % no common window designed: at the current one
cpiOffsetMin = bAxis + sqrt(tUse / nUse) ./ sqrt(s2med);     % sigma_a = 1 / sqrt(s2) on the median sky
if any(bAxis > 0)
    fprintf('   NOMINAL_BIAS %s m on the monitored axes: the separation test fires beyond B + k_FA sigma_SS and the\n', mat2str(bAxis));
    fprintf('   CPI counts only the axis innovation beyond B. The sigma_t design above applies to a spoofer whose offset\n');
    fprintf('   exceeds B; tracking jitter alone below B is invisible to the CPI.\n');
end
fprintf('   constant offset the CPI reports within one window of %d epochs (median sky, B included): %s m\n', nUse, sprintf('%.2f ', cpiOffsetMin));

%% arming: first epoch from which the protection level stays inside ARM_PL_MAX for ARM_EPOCHS epochs
armEpoch = NaN;
if ~isempty(pDiag)
    pl = kMd * sqrt(max(max(pDiag(:, 1:3), [], 2), 0));
    inside = double(pl(:) < CST_spfParam.ARM_PL_MAX);
    w = double(CST_spfParam.ARM_EPOCHS);
    runLen = conv(inside, ones(w, 1), 'valid');              % qualifying epochs in a window of w
    k = find(runLen == w, 1);
    if ~isempty(k), armEpoch = k; end
    if isfinite(armEpoch)
        fprintf('   protection level K_MD * sigma_pos inside ARM_PL_MAX (%.0f m) from epoch %d (%.1f s at 2 Hz): the monitors arm at\n', ...
            CST_spfParam.ARM_PL_MAX, armEpoch, armEpoch / 2);
        fprintf('   epoch %d with ARM_EPOCHS = %d; a larger ARM_EPOCHS only delays that (the host runs 120 for its bias convergence)\n', ...
            armEpoch + w - 1, w);
    else
        fprintf('   protection level never stays inside ARM_PL_MAX (%.0f m) for %d epochs on this run\n', CST_spfParam.ARM_PL_MAX, w);
    end
end

des = struct('axes', axes, 's2', s2, 's2min', s2min, 's2med', s2med, 'nAxis', nAxis, 'pmdAxis', pmdAxis, ...
    'sigmaTreq', sigmaTreq, 'N', N, 'tN', tN, 'kFa', kFa, 'kMd', kMd, 'pmdAtN', pmdAtN, 'sigmaT', sigmaT, ...
    'nominalBias', bAxis, 'cpiOffsetMin', cpiOffsetMin, 'armEpoch', armEpoch, 'current', cur);
end

%------------------------------------------------------------------------
function [sT] = sigmaTrequired(N, s2, pFaPlus, pMdPlus)
% smallest spoofer tracking error at which P_MD(N) meets pMdPlus (P_MD falls with sigma_t); Inf if not by 1 km
lo = 1e-3;  hi = 1e3;
if compute_PMD_eq38(N, s2, hi, pFaPlus) > pMdPlus
    sT = Inf;
else
    for it = 1:60
        mid = sqrt(lo * hi);
        if compute_PMD_eq38(N, s2, mid, pFaPlus) > pMdPlus, lo = mid; else, hi = mid; end
    end
    sT = hi;
end
end
