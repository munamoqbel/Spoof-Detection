function [res] = covPropCheck(spfMeasLog, tKf, tOut, tReacq, run, pDiag, propTelLog)
%COVPROPCHECK  Does the host's P keep its cross terms from one 2 Hz epoch to the next?
%   res = covPropCheck(spfMeasLog, tKf, tOut, tReacq)
%   res = covPropCheck(spfMeasLog, tKf, tOut, tReacq, run)         run: the runCompare struct (t, navAtt, navVel) for
%                                                                   the attitude and acceleration per epoch; without it
%                                                                   the vehicle is taken level with f = -g
%   res = covPropCheck(spfMeasLog, tKf, tOut, tReacq, run, pDiag)  pDiag [K x n]: the diagonal given to runCompare,
%                                                                   checked against diag(postCov)
%   res = covPropCheck(..., pDiag, propTelLog)   propTelLog: the propTel the gate received per epoch
%                                                (STRUCT_SPF.setPropTel: accumPhi, accumQ over the interval that ends
%                                                at the epoch) as a struct array or cell of per-epoch structs, or as
%                                                one struct whose accumPhi and accumQ hold the epochs along the third
%                                                dimension ([n x n x K]) or in a cell. Checks the accumulation and its
%                                                reset, and repeats the one-step prediction with the host's own
%                                                accumPhi and accumQ
%
%   The per-epoch spfMeas (STRUCT_SPF.setKfMeas) holds postCov, the host's P+ after the update, and the rows the host
%   applied (obsMatrix, measNoiseCov, numMeas). For every epoch the previous P+ is propagated over the interval with
%   Phi = expm(F dt) of the 15 INS error states (position <- velocity, velocity <- tilt and accelerometer bias,
%   tilt <- gyro bias; identity for the other states), without Q, the epoch's rows are applied as a standard update,
%   and the result is compared with the logged P+ of the epoch block by block (Frobenius norms, so the sign
%   conventions of the tilt and bias states do not matter). Over one epoch Q adds a few percent to the diagonals and
%   nothing of note to the cross terms, so the cross-term blocks are the test:
%     ratio near 1   the block survives the host's 2 Hz path (100 Hz propagation, updateKF, setKF)
%     ratio near 0   the propagation rebuilds the block and the stored P+ has lost it: removed in updateKF or setKF
%   A second chain propagates the P+ of tOut through the whole outage with the same Phi and rows, to show what a P that
%   keeps its cross terms carries at each epoch next to the logged one.
%
%   State order assumed: position 1:3, velocity 4:6, attitude 7:9, accelerometer bias 10:12, gyro bias 13:15 (NED,
%   metres, m/s, rad); set run.pGyroFirst = true when the host orders the gyro biases before the accelerometer biases.
%   tOut / tReacq: last aided epoch and first epoch back (s), [] for a run without an outage.
%
%   Result: res.t, res.ratio (K x 4: pos-vel, vel-tilt, vel-accel bias, tilt-gyro bias, logged / predicted),
%   res.sigLogged / res.sigOneStep / res.sigChain (K x 4: sigma of position N, E and velocity N, E), res.inOutage,
%   res.medianAided / res.medianOutage (1 x 4), res.hphRatio (S - R against H P- H', median over the aided epochs),
%   res.pDiagMatch (true / false / NaN when pDiag is not given), and with propTelLog: res.accInterval (K x 1, the
%   interval accumPhi(1,4) holds), res.accCoupling (K x 3, vel-tilt, vel-accel bias, tilt-gyro bias blocks of accumPhi
%   relative to the error dynamics) and res.ratioHost (K x 4, the cross blocks of P+ against the prediction with the
%   host's accumPhi and accumQ). res.drops lists the epochs inside the outage where the logged horizontal sigma falls
%   by more than 30% from one epoch to the next: [epoch t tau sigmaBefore sigmaAfter sigmaOneStep numMeas]; when the
%   one-step prediction with the epoch's own rows reproduces the drop it is a measurement update inside the outage,
%   otherwise P was replaced outside the update (a host reset).
%
%   accumPhi(1,4) is the position-from-velocity entry, i.e. the time the accumulation covers: equal to the epoch
%   interval when the accumulator is reset right after the gate read it, 0 when it is reset before, larger when a
%   reset is missed.

if iscell(spfMeasLog), spfMeasLog = [spfMeasLog{:}]; end
K = numel(spfMeasLog);
if nargin < 2 || isempty(tKf), tKf = 0.5 * (0:K-1)'; end
tKf = tKf(:);
if numel(tKf) ~= K
    error('covPropCheck:tKf', 'tKf has %d entries, the log %d epochs', numel(tKf), K);
end
if nargin < 3, tOut = []; end
if nargin < 4, tReacq = []; end
if nargin < 5, run = []; end
if nargin < 6, pDiag = []; end
if nargin < 7, propTelLog = []; end
propTelLog = propTelRows(propTelLog);
G = 9.80665;
n = size(spfMeasLog(1).postCov, 1);
if n < 15 || size(spfMeasLog(1).postCov, 2) ~= n
    error('covPropCheck:postCov', 'postCov must be n x n with n >= 15; got %d x %d', size(spfMeasLog(1).postCov));
end
iPos = 1:3;  iVel = 4:6;  iAtt = 7:9;  iBa = 10:12;  iBg = 13:15;
if ~isempty(run) && isfield(run, 'pGyroFirst') && run.pGyroFirst, iBa = 13:15;  iBg = 10:12; end
haveNav = ~isempty(run) && isfield(run, 'navAtt') && isfield(run, 't') && isfield(run, 'navVel');
hasOut  = ~isempty(tOut) && ~isempty(tReacq);
if hasOut, inOutage = tKf > tOut & tKf < tReacq; else, inOutage = false(K, 1); end
blocks = {iPos(1:2), iVel(1:2); iVel(1:2), iAtt(1:2); iVel(1:2), iBa; iAtt, iBg};
blockName = {'pos-vel', 'vel-tilt', 'vel-accel bias', 'tilt-gyro bias'};

ratio = NaN(K, 4);  sigLogged = NaN(K, 4);  sigOneStep = NaN(K, 4);  sigChain = NaN(K, 4);  hphRatio = NaN(K, 1);
Pc = [];  kOut = [];
if hasOut, kOut = find(tKf <= tOut, 1, 'last'); end
if ~isempty(kOut), Pc = spfMeasLog(kOut).postCov; end
for k = 1:K
    L = spfMeasLog(k).postCov;
    sigLogged(k, :) = sqrt(max([L(1,1), L(2,2), L(4,4), L(5,5)], 0));
    if k == 1, continue; end
    dt = tKf(k) - tKf(k-1);
    Phi = intervalPhi(n, dt, tKf(k-1), run, haveNav, G, iPos, iVel, iAtt, iBa, iBg);
    Pm = Phi * spfMeasLog(k-1).postCov * Phi';
    [Pp, hr] = applyRows(Pm, spfMeasLog(k));
    if ~inOutage(k), hphRatio(k) = hr; end
    sigOneStep(k, :) = sqrt(max([Pp(1,1), Pp(2,2), Pp(4,4), Pp(5,5)], 0));
    for b = 1:4
        pb = norm(Pp(blocks{b, 1}, blocks{b, 2}), 'fro');
        lb = norm(L(blocks{b, 1}, blocks{b, 2}), 'fro');
        if pb > 1e-300, ratio(k, b) = lb / pb; end
    end
    if inOutage(k) && ~isempty(Pc)
        Pc = Phi * Pc * Phi';
        Pc = applyRows(Pc, spfMeasLog(k));
        sigChain(k, :) = sqrt(max([Pc(1,1), Pc(2,2), Pc(4,4), Pc(5,5)], 0));
    end
end
aided = ~inOutage;  aided(1) = false;
medianAided  = medianRows(ratio(aided, :));
medianOutage = medianRows(ratio(inOutage, :));

if haveNav, navTxt = 'C and f from the nav solution'; else, navTxt = 'level vehicle, f = -g (no run given)'; end
fprintf('covPropCheck: %d epochs, %d states, Phi over each interval from %s, no Q\n', K, n, navTxt);
if hasOut
    fprintf('   outage: %d epochs from tOut = %.2f s to tReacq = %.2f s; %d aided epochs\n', nnz(inOutage), tOut, tReacq, nnz(aided));
else
    fprintf('   no outage given: all %d epochs counted as aided\n', nnz(aided));
end
fprintf('   cross block of P+, logged / predicted from the previous P+ (Frobenius norm), median over the epochs:\n');
fprintf('   block              aided    outage   reading\n');
reading = cell(1, 4);
for b = 1:4
    reading{b} = readRatio(medianAided(b), medianOutage(b));
    fprintf('   %-16s %8s %8s   %s\n', blockName{b}, fmtRatio(medianAided(b)), fmtRatio(medianOutage(b)), reading{b});
end
if any(isfinite(hphRatio))
    fprintf('   S - R against H P- H'' with P- the propagated previous P+, median over the aided rows: %.2f (1 = the host updates\n', ...
        medianRows(hphRatio(isfinite(hphRatio))));
    fprintf('   the propagated previous P+; well below 1 = the P- the host used is smaller than that)\n');
end
if hasOut && nnz(inOutage) >= 2
    fprintf('   outage epochs: sigma logged / one-step predicted / chain from tOut (chain = a P that keeps every cross term)\n');
    fprintf('   tau [s]      pos N sigma [m]              vel N sigma [m/s]            vel-tilt   vel-ba    tilt-bg   (logged / predicted)\n');
    tauSel = [0.5 1 2 5 10 25 50 100 150 200 300 400 600];
    ko = find(inOutage);
    tauO = tKf(ko) - tKf(kOut);
    shown = [];
    for ts = [tauSel, tauO(end)]
        [~, j] = min(abs(tauO - ts));
        if ts > tauO(end) + 1 || any(shown == j), continue; end
        shown(end+1) = j; %#ok<AGROW>
        kk = ko(j);
        fprintf('   %7.1f   %7.2f / %7.2f / %8.2f    %7.3f / %7.3f / %7.3f     %s   %s   %s\n', tauO(j), ...
            sigLogged(kk, 1), sigOneStep(kk, 1), sigChain(kk, 1), sigLogged(kk, 3), sigOneStep(kk, 3), sigChain(kk, 3), ...
            fmtRatio(ratio(kk, 2)), fmtRatio(ratio(kk, 3)), fmtRatio(ratio(kk, 4)));
    end
end
drops = [];
if hasOut
    hor = sqrt(sigLogged(:, 1).^2 + sigLogged(:, 2).^2);
    horPred = sqrt(sigOneStep(:, 1).^2 + sigOneStep(:, 2).^2);
    for k = find(inOutage)'
        if k < 2 || ~isfinite(hor(k)) || ~isfinite(hor(k-1)), continue; end
        if hor(k) < 0.7 * hor(k-1)
            drops(end+1, :) = [k, tKf(k), tKf(k) - tKf(kOut), hor(k-1), hor(k), horPred(k), double(spfMeasLog(k).numMeas)]; %#ok<AGROW>
        end
    end
    if ~isempty(drops)
        fprintf('   the logged horizontal sigma drops inside the outage (by more than 30%% from one epoch to the next):\n');
        fprintf('   epoch        t [s]    tau [s]   sigma before -> after [m]   one-step with the epoch''s rows   numMeas   reading\n');
        for j = 1:min(size(drops, 1), 8)
            d = drops(j, :);
            if abs(d(6) - d(5)) <= 0.2 * max(d(4) - d(5), 1e-9) && d(7) >= 4 && d(1) == max(find(inOutage))
                rd = sprintf('the %d rows of the epoch explain the drop and it is the last outage epoch: this is the first GNSS epoch back; give tReacq at or before t = %.2f', d(7), d(2));
            elseif abs(d(6) - d(5)) <= 0.2 * max(d(4) - d(5), 1e-9)
                rd = sprintf('the %d row(s) of the epoch explain the drop: a measurement update inside the outage', d(7));
            else
                rd = 'the rows do not explain it: P was replaced outside the update (a host reset or re-initialisation)';
            end
            fprintf('   %5d   %10.2f  %8.1f   %10.2f -> %-8.2f   %10.2f                      %4d     %s\n', d(1), d(2), d(3), d(4), d(5), d(6), d(7), rd);
        end
        if size(drops, 1) > 8, fprintf('   (%d more)\n', size(drops, 1) - 8); end
        fprintf('   the coast sigma runCompare reads at the end of the outage is the value after the last drop, not what the propagation built\n');
    end
end
lost = medianOutage <= 0.2 | (~hasOut & medianAided <= 0.2);
kept = medianOutage >= 0.7 | (~hasOut & medianAided >= 0.7);
if any(lost(2:4))
    fprintf('   VERDICT: the %s block(s) are rebuilt by the propagation and gone from the stored P+ at every epoch: the host''s 2 Hz\n', strjoin(blockName(lost), ', '));
    fprintf('            path (updateKF or setKF) removes them. Nothing in Q or Phi can repair that; find the line that writes P there\n');
elseif all(kept)
    fprintf('   VERDICT: every cross block survives the 2 Hz path and the logged P+ follows its own propagation. A coast sigma far below\n');
    fprintf('            the implied one then comes from a drop listed above, or from a pDiag that is not postCov (see the pDiag line)\n');
else
    fprintf('   VERDICT: cross blocks partly kept (ratios between 0.2 and 0.7): the host shrinks them each epoch rather than zeroing them\n');
end
pDiagMatch = NaN;
if ~isempty(pDiag)
    pDiag = pDiagRows(pDiag, K, n);
    nd = min(n, size(pDiag, 2));
    if size(pDiag, 1) ~= K
        fprintf('   pDiag has %d rows, the log %d epochs: not compared\n', size(pDiag, 1), K);
    else
        rel = zeros(K, 1);
        for k = 1:K
            d = diag(spfMeasLog(k).postCov)';
            rel(k) = max(abs(pDiag(k, 1:nd) - d(1:nd)) ./ max(abs(d(1:nd)), 1e-12));
        end
        [relMax, kMax] = max(rel);
        pDiagMatch = relMax < 1e-6;
        if pDiagMatch
            fprintf('   pDiag given to runCompare is diag(postCov): its growth law describes this P+\n');
        else
            fprintf('   pDiag given to runCompare is NOT diag(postCov): largest relative difference %.3g at epoch %d (t = %.2f s). Its growth\n', relMax, kMax, tKf(kMax));
            fprintf('   law describes another object; redo runCompare with pDiag = the diagonal of postCov\n');
        end
    end
end
accInterval = NaN(K, 1);  accCoupling = NaN(K, 3);  ratioHost = NaN(K, 4);  sigHost = NaN(K, 4);
if ~isempty(propTelLog)
    if numel(propTelLog) ~= K
        fprintf('   propTelLog has %d entries, the log %d epochs: not used\n', numel(propTelLog), K);
    else
        for k = 2:K
            PhiH = propTelLog(k).accumPhi;  QH = propTelLog(k).accumQ;
            if isempty(PhiH) || size(PhiH, 1) < 15, continue; end
            dt = tKf(k) - tKf(k-1);
            PhiE = intervalPhi(n, dt, tKf(k-1), run, haveNav, G, iPos, iVel, iAtt, iBa, iBg);
            accInterval(k) = PhiH(iPos(1), iVel(1));
            accCoupling(k, :) = [norm(PhiH(iVel, iAtt), 'fro') / max(norm(PhiE(iVel, iAtt), 'fro'), 1e-300), ...
                                 norm(PhiH(iVel, iBa),  'fro') / max(norm(PhiE(iVel, iBa),  'fro'), 1e-300), ...
                                 norm(PhiH(iAtt, iBg),  'fro') / max(norm(PhiE(iAtt, iBg),  'fro'), 1e-300)];
            if isempty(QH) || ~isequal(size(QH), size(PhiH)), QH = zeros(size(PhiH)); end
            Pm = PhiH * spfMeasLog(k-1).postCov * PhiH' + QH;
            Pp = applyRows(Pm, spfMeasLog(k));
            L  = spfMeasLog(k).postCov;
            sigHost(k, :) = sqrt(max([Pp(1,1), Pp(2,2), Pp(4,4), Pp(5,5)], 0));
            for b = 1:4
                pb = norm(Pp(blocks{b, 1}, blocks{b, 2}), 'fro');
                lb = norm(L(blocks{b, 1}, blocks{b, 2}), 'fro');
                if pb > 1e-300, ratioHost(k, b) = lb / pb; end
            end
        end
        dtE = [NaN; diff(tKf)];
        okI = isfinite(accInterval) & isfinite(dtE);
        fprintf('   host accumPhi (propTel per epoch):\n');
        if any(okI)
            rI = accInterval(okI) ./ dtE(okI);
            off = okI & abs(accInterval ./ dtE - 1) > 0.05;
            fprintf('   interval it covers, accumPhi(1,4) / epoch interval: median %.3f, min %.3f, max %.3f; %d of %d epochs off by more than 5%%', ...
                median(rI), min(rI), max(rI), nnz(off), nnz(okI));
            if any(off), fprintf(' (first at epoch %d, t = %.2f s)', find(off, 1), tKf(find(off, 1))); end
            fprintf('\n');
            if median(rI) < 0.05
                fprintf('   VERDICT: the accumulator is reset before the gate reads it (accumPhi = I at the call): the gate propagates nothing\n');
            elseif median(rI) > 1.5
                fprintf('   VERDICT: the accumulator is not reset after the gate reads it: it covers %.1f epochs instead of one\n', median(rI));
            elseif abs(median(rI) - 1) > 0.05
                fprintf('   VERDICT: the accumulator covers %.0f%% of the epoch interval: the reset is one or more ticks off the gate call\n', 100 * median(rI));
            else
                fprintf('   the accumulator covers one epoch interval: reset right after the gate reads it\n');
            end
            cm = medianRows(accCoupling(okI, :));
            fprintf('   coupling blocks of accumPhi relative to the error dynamics (1 = same size): vel-tilt %.2f, vel-accel bias %.2f, tilt-gyro bias %.2f\n', cm);
        else
            fprintf('   (no usable accumPhi in propTelLog)\n');
        end
        okH = any(isfinite(ratioHost), 2);
        if any(okH)
            mh = medianRows(ratioHost(okH, :));
            fprintf('   P+ against the prediction with the host''s accumPhi and accumQ, logged / predicted, median: pos-vel %s, vel-tilt %s, vel-accel bias %s, tilt-gyro bias %s\n', ...
                fmtRatio(mh(1)), fmtRatio(mh(2)), fmtRatio(mh(3)), fmtRatio(mh(4)));
            sh = sigLogged(okH, :) ./ sigHost(okH, :);
            fprintf('   position and velocity sigma, logged / predicted with accumQ, median: pos N %.3f E %.3f, vel N %.3f E %.3f (1 = the 100 Hz P follows the accumulators)\n', ...
                median(sh(isfinite(sh(:, 1)), 1)), median(sh(isfinite(sh(:, 2)), 2)), median(sh(isfinite(sh(:, 3)), 3)), median(sh(isfinite(sh(:, 4)), 4)));
        end
    end
end
res = struct('t', tKf, 'ratio', ratio, 'sigLogged', sigLogged, 'sigOneStep', sigOneStep, 'sigChain', sigChain, ...
    'inOutage', inOutage, 'medianAided', medianAided, 'medianOutage', medianOutage, 'hphRatio', hphRatio, 'pDiagMatch', pDiagMatch, ...
    'accInterval', accInterval, 'accCoupling', accCoupling, 'ratioHost', ratioHost, 'drops', drops);
end

%------------------------------------------------------------------------
function [Phi] = intervalPhi(n, dt, t0, run, haveNav, G, iPos, iVel, iAtt, iBa, iBg)
% Phi = expm(F dt) of the INS error states; C and f at the start of the interval from the nav solution when given
C = eye(3);  f = [0; 0; -G];
if haveNav
    tt = run.t(:);
    e  = interp1(tt, run.navAtt, t0, 'linear', 'extrap');
    C  = dcmFromEuler(e);
    i0 = find(tt <= t0, 1, 'last');  i0 = max(min(i0, numel(tt) - 5), 6);
    a  = (run.navVel(i0 + 5, :) - run.navVel(i0 - 5, :))' / (tt(i0 + 5) - tt(i0 - 5));
    f  = a - [0; 0; G];
end
F = zeros(n);
F(iPos, iVel) = eye(3);
F(iVel, iAtt) = [0, -f(3), f(2); f(3), 0, -f(1); -f(2), f(1), 0];
F(iVel, iBa)  = C;
F(iAtt, iBg)  = C;
Phi = expm(F * dt);
end

function [Pp, hr] = applyRows(Pm, s)
% the epoch's rows applied to a propagated P as a standard update; hr = median over the rows of
% (S_logged - R) / (H Pm H')
Pp = Pm;  hr = NaN;
m = double(s.numMeas);
if m <= 0, return; end
H = s.obsMatrix(1:m, :);
R = s.measNoiseCov(1:m, 1:m);
hph = H * Pm * H';
S = hph + R;
Kg = (Pm * H') / S;
Pp = Pm - Kg * H * Pm;
Pp = (Pp + Pp') / 2;
if isfield(s, 'innovationCov') && ~isempty(s.innovationCov)
    d = diag(hph);
    ok = d > 1e-12;
    if any(ok)
        hr = median((diag(s.innovationCov(1:m, 1:m)) - diag(R)) ./ max(d, 1e-12));
    end
end
end

function [L] = propTelRows(propTelLog)
% propTel per epoch as a struct array, from a struct array or cell of per-epoch structs, or from one struct whose
% accumPhi / accumQ hold the epochs along the third dimension or in a cell
if iscell(propTelLog), propTelLog = [propTelLog{:}]; end
L = propTelLog;
if isempty(propTelLog) || ~isstruct(propTelLog) || numel(propTelLog) ~= 1 || ~isfield(propTelLog, 'accumPhi'), return; end
A = propTelLog.accumPhi;
Qa = [];
if isfield(propTelLog, 'accumQ'), Qa = propTelLog.accumQ; end
if iscell(A)
    Kp = numel(A);
    L = repmat(struct('accumPhi', [], 'accumQ', []), Kp, 1);
    for k = 1:Kp
        L(k).accumPhi = A{k};
        if iscell(Qa) && numel(Qa) >= k, L(k).accumQ = Qa{k}; end
    end
elseif ndims(A) == 3 && size(A, 3) > 1
    if size(A, 1) ~= size(A, 2) && size(A, 2) == size(A, 3)
        A = permute(A, [2 3 1]);                               % epochs first: [K x n x n]
        if ~isempty(Qa) && ndims(Qa) == 3, Qa = permute(Qa, [2 3 1]); end
    end
    Kp = size(A, 3);
    L = repmat(struct('accumPhi', [], 'accumQ', []), Kp, 1);
    for k = 1:Kp
        L(k).accumPhi = A(:, :, k);
        if ~isempty(Qa) && ndims(Qa) == 3 && size(Qa, 3) >= k, L(k).accumQ = Qa(:, :, k); end
    end
end
end

function [D] = pDiagRows(pDiag, K, n)
% the P diagonal per epoch as [K x n] from the forms runCompare accepts: [K x n], [n x K], or a cell holding per
% epoch the diagonal or the full P
if iscell(pDiag) && all(cellfun(@numel, pDiag(:)) == 1) && numel(pDiag) > max(size(pDiag))
    pDiag = cell2mat(pDiag);                                  % a K x n (or n x K) cell of scalars
end
if iscell(pDiag)
    D = NaN(numel(pDiag), n);
    for k = 1:numel(pDiag)
        v = pDiag{k};
        if isempty(v), continue; end
        if size(v, 1) == size(v, 2) && size(v, 1) > 1, v = diag(v); end
        v = v(:)';
        m = min(n, numel(v));
        D(k, 1:m) = v(1:m);
    end
elseif size(pDiag, 1) ~= K && size(pDiag, 2) == K
    D = pDiag';
else
    D = pDiag;
end
end

function [m] = medianRows(x)
if isempty(x), m = NaN(1, size(x, 2)); return; end
m = NaN(1, size(x, 2));
for c = 1:size(x, 2)
    v = x(isfinite(x(:, c)), c);
    if ~isempty(v), m(c) = median(v); end
end
end

function [s] = fmtRatio(r)
if isfinite(r), s = sprintf('%.2f', r); else, s = 'n/a'; end
end

function [s] = readRatio(a, o)
r = o;  if ~isfinite(r), r = a; end
if ~isfinite(r), s = 'n/a';
elseif r >= 0.7, s = 'kept';
elseif r <= 0.2, s = 'removed each epoch';
else, s = 'partly kept'; end
end

function [C] = dcmFromEuler(e)
% body -> nav, ZYX: C = Rz(yaw) Ry(pitch) Rx(roll)
cr = cos(e(1)); sr = sin(e(1)); cp = cos(e(2)); sp = sin(e(2)); cy = cos(e(3)); sy = sin(e(3));
C = [cy*cp, cy*sp*sr - sy*cr, cy*sp*cr + sy*sr;
     sy*cp, sy*sp*sr + cy*cr, sy*sp*cr - cy*sr;
     -sp,   cp*sr,            cp*cr];
end
