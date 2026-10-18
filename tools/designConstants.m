function [des] = designConstants(spfMeasLog, sigmaT, pDiag, nominalBias)
%DESIGNCONSTANTS  Offline design of the gate's monitor constants from the host's own update records.
%   des = designConstants(spfMeasLog)                     spfMeasLog: struct array or cell of the per-epoch spfMeas
%   des = designConstants(spfMeasLog, sigmaT)             sigmaT: spoofer tracking-error std (m), default 0.10
%   des = designConstants(spfMeasLog, sigmaT, pDiag)      pDiag: P diagonal per epoch, position states 1:3, as
%                                                         [K x n], [n x K] or a cell per epoch (diagonal or full P).
%                                                         Selects the epochs the gate would be armed on and reports
%                                                         from which epoch it arms
%   des = designConstants(spfMeasLog, sigmaT, pDiag, B)   B: the host's NOMINAL_BIAS per NED axis (m); default the
%                                                         CST_spfParam on the path
%
%   QUICKSTART step 2 on recorded data. Per monitored axis a (CST_spfParam.MONITORED_AXES), on every epoch with
%   at least REVAL_MIN_MEAS rows:
%       s2(k) = f' S^-1 f,  f = H(:, a)                                   (Eq. 20 of Kujur 2024)
%   1/sqrt(s2) is the filter's prior uncertainty along the axis as the rows of that epoch see it. The design
%   sky is the 5th percentile of s2 over the armed epochs: enough rows and a protection level K_MD sigma_pos
%   inside ARM_PL_MAX (with pDiag; all usable epochs without it). The minimum over the whole run is printed
%   for information only: it sits on the first epoch after an outage, where P is legitimately large and the
%   gate is disarmed.
%   Omega = s2 sigmaT^2 is the spoofer's tracking jitter against that uncertainty. The missed-detection
%   requirement P_MD+ (kujur_params, split across the axes) is reachable only when Omega is of order one:
%   with a prior of metres, as a code-based filter has, a 0.1 m spoofer is invisible to the CPI at any window
%   length, and the tool prints the sigma_t it would take at the current window and at MAX_WINDOW.
%   Reachable: N(a) per axis from Eq. 38, common window N = max N(a) + 1. Not reachable: the current
%   WINDOW_LENGTH is kept. At the window used:
%       T_N  = gaminv(1 - P_FA_plus / nAxes, N/2, 2)                       (Eq. 35, CPI_THRESHOLD)
%       k_FA = norminv(1 - (P_FA_minus / nAxes / N) / 2)                   (K_FALSE_ALERT)
%       k_MD = norminv(1 - P_MD_minus)                                     (K_MISSED_DETECTION)
%   These depend on the budgets and N only, so on a code-based filter they come out equal to the published
%   ones; what the data sets is the offset the CPI can see, printed per axis on the median armed sky:
%       d50 = B + sqrt(T_N / N) / sqrt(s2med)        the constant offset at which q reaches T_N within one window
%       dMD = B + sqrt(lambda / (N s2med))           the offset missed with probability P_MD+ (normal
%                                                    approximation of the non-central chi-square)
%   The separation test fires beyond B + k_FA sigma_SS and the CPI counts only the axis innovation beyond B.
%   Rerun after any change of Q, R, the measurement set or the budgets: S carries the host's P and R.
%   MATLAB desktop tool, not codegen (gaminv, norminv). Runs under Octave with tools/octave_shim.
%
%   Result: axes, s2 [K x nAxes], used and armed (logical K x 1), s2min and s2minEpoch (whole run), s2minArmed,
%   s2des (design sky), s2med (median armed sky), nAxis, pmdAxis, sigmaTreq [nAxes x 2] (at the current window
%   and at MAX_WINDOW, when not reachable), N (designed common window, NaN when not reachable), nUse, tN, kFa,
%   kMd (the window used and the constants at it), pmdAtN, sigmaT, nominalBias, cpiOffsetMin (d50),
%   cpiOffsetMd (dMD), armEpoch, current (the CST_spfParam values).

MAX_WINDOW = 240;                                            % epochs (2 min at 2 Hz): longer is not a monitor
DESIGN_PCT = 5;                                              % percentile of the armed epochs taken as the design sky

if iscell(spfMeasLog), spfMeasLog = [spfMeasLog{:}]; end
if nargin < 2 || isempty(sigmaT), sigmaT = 0.10; end
if nargin < 3, pDiag = []; end
if nargin < 4 || isempty(nominalBias)
    nominalBias = double(CST_spfParam.NOMINAL_BIAS(:))';  biasSrc = 'CST_spfParam on the path';
else
    nominalBias = double(nominalBias(:))';  biasSrc = 'argument';
end
K = numel(spfMeasLog);
prm = kujur_params();
axes = double(CST_spfParam.MONITORED_AXES(:))';
nAx = numel(axes);
pFaPlus  = prm.P_FA_plus  / nAx;                           % CPI false-alert allocation per axis
pFaMinus = prm.P_FA_minus / nAx;                           % SS  false-alert allocation per axis (total over the window)
minRows  = double(CST_spfParam.REVAL_MIN_MEAS);
cur = struct('N', double(CST_spfParam.WINDOW_LENGTH), 'tN', CST_spfParam.CPI_THRESHOLD, ...
    'kFa', CST_spfParam.K_FALSE_ALERT, 'kMd', CST_spfParam.K_MISSED_DETECTION);
kMd = norminv(1 - prm.P_MD_minus);

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

%% the epochs the gate would be armed on: enough rows and a protection level inside ARM_PL_MAX
pl = NaN(K, 1);
if ~isempty(pDiag)
    D = pDiagRows(pDiag, K);
    if size(D, 1) == K && size(D, 2) >= 3
        pl = kMd * sqrt(max(max(D(:, 1:3), [], 2), 0));
    else
        fprintf('   pDiag has %d rows, the log %d epochs: not used\n', size(D, 1), K);
    end
end
armNote = '';
if any(isfinite(pl))
    armed = used & pl < CST_spfParam.ARM_PL_MAX;
    if nnz(armed) < 10
        armed = used;  armNote = ' (fewer than 10 armed epochs: all usable epochs taken)';
    end
else
    armed = used;  armNote = ' (no pDiag: all usable epochs taken as armed)';
end
kUsed = find(used);  kArmed = find(armed);
[s2min, iMin] = min(s2(used, :), [], 1);  s2minEpoch = kUsed(iMin)';
[s2minArmed, iMinA] = min(s2(armed, :), [], 1);  s2minArmedEpoch = kArmed(iMinA)';
s2des = pctRows(s2(armed, :), DESIGN_PCT);
s2med = median(s2(armed, :), 1);

fprintf('designConstants: %d epochs, %d usable (>= %d rows, finite S), %d armed (K_MD sigma_pos inside ARM_PL_MAX %.0f m)%s\n', ...
    K, nnz(used), minRows, nnz(armed), CST_spfParam.ARM_PL_MAX, armNote);
fprintf('   sigma_t %.2f m; budgets P_FA+ %.1e  P_FA- %.1e  P_MD+ %.1e  P_MD- %.1e, split over %d axes; NOMINAL_BIAS %s m (%s)\n', ...
    sigmaT, prm.P_FA_plus, prm.P_FA_minus, prm.P_MD_plus, prm.P_MD_minus, nAx, mat2str(nominalBias(axes)), biasSrc);
fprintf('1. prior uncertainty along the monitored axes, 1/sqrt(s2) [m], s2 = f'' S^-1 f (Eq. 20)\n');
fprintf('   axis  state   whole-run min (epoch)   armed min (epoch)   design sky (%dth pct)   median armed sky\n', DESIGN_PCT);
for a = 1:nAx
    fprintf('   %-4d  %-5d   %10.2f (%5d)      %10.2f (%5d)      %10.2f              %10.2f\n', a, axes(a), ...
        1 / sqrt(s2min(a)), s2minEpoch(a), 1 / sqrt(s2minArmed(a)), s2minArmedEpoch(a), 1 / sqrt(s2des(a)), 1 / sqrt(s2med(a)));
end

%% CPI missed detection per axis on the design sky
nAxis = NaN(1, nAx);  pmdAxis = NaN(1, nAx);  sigmaTreq = NaN(nAx, 2);  omega = s2des * sigmaT^2;
fprintf('2. CPI missed detection at sigma_t = %.2f m on the design sky (Eq. 38, P_MD+ per axis %.1e)\n', sigmaT, prm.P_MD_plus);
fprintf('   axis   Omega       reachable within %d epochs   N(a)   P_MD at N(a)   sigma_t needed at N = %d / %d [m]\n', MAX_WINDOW, cur.N, MAX_WINDOW);
for a = 1:nAx
    [nA, pA] = windowFor(s2des(a), sigmaT, pFaPlus, prm.P_MD_plus, MAX_WINDOW);
    if isfinite(nA)
        nAxis(a) = nA;  pmdAxis(a) = pA;
        fprintf('   %-4d   %.3g    yes                           %4d   %.2e       -\n', a, omega(a), nA, pA);
    else
        sigmaTreq(a, 1) = sigmaTrequired(cur.N, s2des(a), pFaPlus, prm.P_MD_plus);
        sigmaTreq(a, 2) = sigmaTrequired(MAX_WINDOW, s2des(a), pFaPlus, prm.P_MD_plus);
        fprintf('   %-4d   %.3g    no                               -   -              %.2f / %.2f\n', a, omega(a), sigmaTreq(a, 1), sigmaTreq(a, 2));
    end
end

%% the window used and the constants at it
N = NaN;
if any(isfinite(nAxis))
    N = max(nAxis(isfinite(nAxis))) + 1;                     % one epoch of margin, as in the quickstart
    nUse = N;  why = 'designed, max N(a) + 1';
else
    nUse = cur.N;  why = 'the current one: P_MD+ not reachable on this filter';
end
tN  = gaminv(1 - pFaPlus, nUse / 2, 2);
kFa = norminv(1 - (pFaMinus / nUse) / 2);
pmdAtN = NaN(1, nAx);
for a = 1:nAx, pmdAtN(a) = compute_PMD_eq38(nUse, s2des(a), sigmaT, pFaPlus); end
fprintf('3. constants at the window used, N = %d (%s)\n', nUse, why);
fprintf('                              from this run     CST_spfParam now\n');
fprintf('   WINDOW_LENGTH        N   %14d    %14d\n', nUse, cur.N);
fprintf('   CPI_THRESHOLD        T_N %14.6f    %14.6f\n', tN, cur.tN);
fprintf('   K_FALSE_ALERT        k_FA%14.6f    %14.6f\n', kFa, cur.kFa);
fprintf('   K_MISSED_DETECTION   k_MD%14.6f    %14.6f\n', kMd, cur.kMd);
fprintf('   P_MD+ at this N per axis: %s(requirement %.1e)\n', sprintf('%.2e ', pmdAtN), prm.P_MD_plus);

%% the offsets the CPI sees on the median armed sky, within one window
bAxis = nominalBias(axes);
cpiOffsetMin = bAxis + sqrt(tN / nUse) ./ sqrt(s2med);      % q reaches T_N: sigma_a = 1 / sqrt(s2)
lambdaMd = lambdaFor(nUse, tN, prm.P_MD_plus);
cpiOffsetMd  = bAxis + sqrt(lambdaMd ./ (nUse * s2med));
fprintf('4. constant offset the CPI sees on the median armed sky within one window of %d epochs, B included [m]\n', nUse);
fprintf('   axis    B     d50 (q reaches T_N)   dMD (missed with P_MD+, normal approximation)\n');
for a = 1:nAx
    fprintf('   %-4d  %5.2f   %10.2f            %10.2f\n', a, bAxis(a), cpiOffsetMin(a), cpiOffsetMd(a));
end
if any(bAxis > 0)
    fprintf('   the separation test fires beyond B + k_FA sigma_SS; the CPI counts only the axis innovation beyond B\n');
end

%% arming: first epoch from which the protection level stays inside ARM_PL_MAX for ARM_EPOCHS epochs
armEpoch = NaN;
if any(isfinite(pl))
    inside = double(pl(:) < CST_spfParam.ARM_PL_MAX);
    w = double(CST_spfParam.ARM_EPOCHS);
    runLen = conv(inside, ones(w, 1), 'valid');              % qualifying epochs in a window of w
    k = find(runLen == w, 1);
    if ~isempty(k), armEpoch = k; end
    if isfinite(armEpoch)
        fprintf('5. arming: K_MD sigma_pos inside ARM_PL_MAX (%.0f m) from epoch %d (%.1f s at 2 Hz); the monitors arm at epoch %d with\n', ...
            CST_spfParam.ARM_PL_MAX, armEpoch, armEpoch / 2, armEpoch + w - 1);
        fprintf('   ARM_EPOCHS = %d; a larger ARM_EPOCHS only delays that (the host runs 120 for its bias convergence)\n', w);
    else
        fprintf('5. arming: the protection level never stays inside ARM_PL_MAX (%.0f m) for %d epochs on this run\n', CST_spfParam.ARM_PL_MAX, w);
    end
end

%% verdict
if isfinite(N)
    fprintf('VERDICT: window N = %d designed at sigma_t %.2f m. Set WINDOW_LENGTH, CPI_THRESHOLD and K_FALSE_ALERT from the\n', N, sigmaT);
    fprintf('         "from this run" column in CST_spfParam (root and deploy)');
    if N ~= cur.N
        fprintf('; PROBATION_LENGTH follows it (N + 8 in the harness\n         design) and the monitor buffers are sized from WINDOW_LENGTH: rerun runAllTests and deployEquivalence\n');
    else
        fprintf('\n');
    end
else
    fprintf('VERDICT: P_MD+ at sigma_t %.2f m is out of reach on this filter (a prior of %s m along the axes against a %.2f m\n', ...
        sigmaT, strtrim(sprintf('%.1f ', 1 ./ sqrt(s2des))), sigmaT);
    fprintf('         jitter). The window stays at %d and T_N, k_FA, k_MD depend on the budgets only: nothing to change in\n', nUse);
    fprintf('         CST_spfParam. The CPI catches constant offsets from about %s m per axis within one window (B included);\n', strtrim(sprintf('%.1f ', cpiOffsetMin)));
    fprintf('         the separation test takes the rest over the anchor horizon\n');
end

des = struct('axes', axes, 's2', s2, 'used', used, 'armed', armed, 's2min', s2min, 's2minEpoch', s2minEpoch, ...
    's2minArmed', s2minArmed, 's2des', s2des, 's2med', s2med, 'nAxis', nAxis, 'pmdAxis', pmdAxis, ...
    'sigmaTreq', sigmaTreq, 'N', N, 'nUse', nUse, 'tN', tN, 'kFa', kFa, 'kMd', kMd, 'pmdAtN', pmdAtN, 'sigmaT', sigmaT, ...
    'nominalBias', bAxis, 'cpiOffsetMin', cpiOffsetMin, 'cpiOffsetMd', cpiOffsetMd, 'armEpoch', armEpoch, 'current', cur);
end

%------------------------------------------------------------------------
function [nMin, pMd] = windowFor(s2, sigmaT, pFaPlus, pMdPlus, maxWindow)
% smallest window at which Eq. 38 meets pMdPlus; NaN when not within maxWindow (the search of solve_N_min)
nMin = NaN;  pMd = NaN;
if compute_PMD_eq38(maxWindow, s2, sigmaT, pFaPlus) > pMdPlus, return; end
n = 1;  p = 1;
while p > pMdPlus && n < maxWindow
    n = n + 1;
    p = compute_PMD_eq38(n, s2, sigmaT, pFaPlus);
end
nMin = n;  pMd = p;
end

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

function [lambda] = lambdaFor(N, tN, pMdPlus)
% non-centrality at which a chi-square of N degrees of freedom falls below tN with probability pMdPlus, with the
% normal approximation of the non-central chi-square (mean N + lambda, variance 2 (N + 2 lambda))
z = norminv(1 - pMdPlus);
lo = 0;  hi = 1e8;
for it = 1:200
    mid = (lo + hi) / 2;
    if N + mid - z * sqrt(2 * (N + 2 * mid)) < tN, lo = mid; else, hi = mid; end
end
lambda = hi;
end

function [v] = pctRows(x, pct)
% the pct-th percentile of each column (nearest rank), over the finite entries
v = NaN(1, size(x, 2));
for c = 1:size(x, 2)
    s = sort(x(isfinite(x(:, c)), c));
    if ~isempty(s), v(c) = s(max(1, min(numel(s), ceil(pct / 100 * numel(s))))); end
end
end

function [D] = pDiagRows(pDiag, K)
% the P diagonal per epoch as [K x n] from [K x n], [n x K], a cell per epoch (diagonal or full P) or a K x n cell of scalars
if iscell(pDiag) && all(cellfun(@numel, pDiag(:)) == 1) && numel(pDiag) > max(size(pDiag))
    pDiag = cell2mat(pDiag);
end
if iscell(pDiag)
    n = 0;
    for k = 1:numel(pDiag), n = max(n, numel(diag(pDiag{k}))); end
    D = NaN(numel(pDiag), max(n, 3));
    for k = 1:numel(pDiag)
        v = pDiag{k};
        if isempty(v), continue; end
        if size(v, 1) == size(v, 2) && size(v, 1) > 1, v = diag(v); end
        v = v(:)';
        D(k, 1:numel(v)) = v;
    end
elseif size(pDiag, 1) ~= K && size(pDiag, 2) == K
    D = pDiag';
else
    D = pDiag;
end
end
