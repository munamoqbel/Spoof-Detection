function [res] = reacqChecks(filePre, fileFirst, fileAcc, fileRef)
%REACQCHECKS  Re-acquisition diagnostics on four saved kfUpdate epochs.
%   reacqChecks('kf_01674.69.mat', 'kf_01857.27.mat', 'kf_01857.77.mat', 'kf_01875.00.mat')
%
%   filePre    last accepted update before the outage
%   fileFirst  first 2 Hz call after the outage (one row)
%   fileAcc    first accepted multi-row update after the outage (the wrong-way step)
%   fileRef    a healthy epoch later on
%
%   Each file holds the variables of the capture list (names below). Missing
%   variables are reported and the checks that need them are skipped.
%   Prints checks 1-8; returns them in res. MATLAB desktop tool, not codegen.
%
%   Expected variables per file
%     covarIn measH R measZ nonLinZ innov K xUpdate xPrior xPost covar
%     rowMap        [m x 4]: type (1 pseudorange, 2 range-rate, 3 pressure), constellation, PRN, accept
%     mechPos       [lat lon h] rad, rad, m      refPos     [lat lon h]
%     mechVel       [vN vE vD] m/s               refVel     [vN vE vD]
%     mechAtt       [roll pitch yaw] rad         refAtt     [roll pitch yaw] rad (optional)
%     satPos        [k x 3] satellite ECEF position (m) per pseudorange row, same order as the rows
%     rhoMeas       [k x 1] corrected pseudorange (m), rangeRateMeas [k x 1] (m/s), refPosEcef [3 x 1]
%                   (optional: enable section 0, the measurement-quality check against the truth)
%     tRx           receiver time tag (s), optional (else parsed from the file name)
%     statePropagated covPropagated accumPhi accumQ   (fileFirst only)
%   State layout: 1-3 pos NED (m), 4-6 vel NED, 7-9 attitude, 16/18/21 clock bias, 17/19/22 drift.

POS = 1:3;  VEL = 4:6;  ATT = 7:9;
G   = 9.80665;

E = {loadEpoch(filePre), loadEpoch(fileFirst), loadEpoch(fileAcc), loadEpoch(fileRef)};
names = {'PRE (last before outage)', 'FIRST (first call back)', 'ACC (first accepted update)', 'REF (healthy)'};
res = struct();

fprintf('\n================ re-acquisition checks ================\n');
for k = 1:4
    fprintf('%-30s t = %10.2f  file: %s\n', names{k}, E{k}.t, E{k}.file);
end

% ---------------------------------------------------------------- per-epoch checks
for k = [3 4 1 2]
    S = E{k};
    fprintf('\n\n############ %s  (t = %.2f) ############\n', names{k}, S.t);
    if isempty(S.measH)
        fprintf('  no measH in this file: skipped\n');
        continue
    end
    [prRows, rrRows, baroRows] = rowsOf(S);
    n = numel(S.xPrior);
    fprintf('  rows: %d pseudorange, %d range-rate, %d pressure, states %d\n', ...
        numel(prRows), numel(rrRows), numel(baroRows), n);

    % ---- 0. measurement quality against the truth (needs satPos, rhoMeas, refPosEcef) ----
    fprintf('\n-- 0. pseudoranges against the truth (r1 = rho_meas - |sat - ref|, clock = median per constellation) --\n');
    if ~isempty(S.satPos) && ~isempty(S.rhoMeas) && numel(S.refPosEcef) == 3 && size(S.satPos, 1) == numel(S.rhoMeas)
        kSat   = numel(S.rhoMeas);
        rhoRef = sqrt(sum((S.satPos - repmat(S.refPosEcef', kSat, 1)) .^ 2, 2));
        r1     = S.rhoMeas - rhoRef;
        if ~isempty(S.rowMap) && kSat == numel(prRows)
            con = S.rowMap(prRows, 2);
            prn = S.rowMap(prRows, 3);
        else
            con = ones(kSat, 1);
            prn = (1:kSat)';
        end
        d1 = r1;
        for c = unique(con)'
            idx = (con == c);
            clk = median(r1(idx));
            d1(idx) = r1(idx) - clk;
            fprintf('  constellation %d: clock (median r1) = %9.3f m   d1 rms = %7.3f m   max |d1| = %7.3f m\n', ...
                c, clk, rms(d1(idx)), max(abs(d1(idx))));
        end
        fprintf('  PRN   con   rho_ref [km]     r1 [m]      d1 [m]\n');
        for i = 1:kSat
            fprintf('  %3d   %2d   %10.1f   %10.3f   %9.3f\n', prn(i), con(i), rhoRef(i) / 1e3, r1(i), d1(i));
        end
        if any(rhoRef < 1.8e7 | rhoRef > 2.7e7)
            fprintf('  WARNING: some rho_ref outside 18000-27000 km: frame or unit mismatch between satPos and refPosEcef\n');
        end
        if ~isempty(S.rangeRateMeas) && numel(S.rangeRateMeas) == kSat && kSat >= 3
            A = [S.rangeRateMeas(:), ones(kSat, 1)];
            cf = A \ d1;
            fprintf('  d1 vs range rate: slope = %.4f s (a time offset between ranges and satellite positions), intercept = %.2f m\n', cf(1), cf(2));
        end
    else
        fprintf('  satPos / rhoMeas / refPosEcef missing or inconsistent sizes: skipped\n');
    end

    % ---- 1. reproduce the host update -------------------------------------
    fprintf('\n-- 1. reproduce the host update --\n');
    if ~isempty(S.K) && ~isempty(S.innov)
        xHat = S.xPrior + S.K * S.innov;
        fprintf('  max |xPrior + K*innov - xPost|          = %.3e\n', max(abs(xHat - S.xPost)));
    end
    if ~isempty(S.measZ) && ~isempty(S.nonLinZ) && ~isempty(S.innov)
        fprintf('  max |(measZ - nonLinZ) - innov|          = %.3e\n', max(abs(S.measZ - S.nonLinZ - S.innov)));
    end
    if ~isempty(S.covarIn) && ~isempty(S.R)
        Sm = S.measH * S.covarIn * S.measH' + S.R;
        Kc = S.covarIn * S.measH' / Sm;
        if ~isempty(S.K)
            fprintf('  max |P H''(HPH''+R)^-1 - K|  (rel.)       = %.3e\n', max(abs(Kc(:) - S.K(:))) / max(abs(S.K(:))));
        end
        xTb = S.xPrior + Kc * S.innov;
        fprintf('  textbook xHat(1:6) = %s\n', vec2str(xTb(1:6)));
        fprintf('  host     xPost(1:6) = %s\n', vec2str(S.xPost(1:6)));
        res.(fld(k)).xTextbook = xTb;
    end

    % ---- 2. h(x) against H ----------------------------------------------
    fprintf('\n-- 2. nonLinZ (h(x)) against H*xPrior --\n');
    if ~isempty(S.nonLinZ)
        Hx  = S.measH * S.xPrior;
        H2  = S.measH;  H2(prRows, 1:2) = 2 * H2(prRows, 1:2);   % what the doubled builder would give
        Hx2 = H2 * S.xPrior;
        fprintf('  row   nonLinZ      H*xPrior     H2*xPrior   (H2 = north/east doubled on pseudorange rows)\n');
        for r = [prRows(1:min(4, end)) rrRows(1:min(2, end))]
            fprintf('  %3d  %10.3f  %10.3f  %10.3f\n', r, S.nonLinZ(r), Hx(r), Hx2(r));
        end
        fprintf('  rms |nonLinZ - H*xPrior|  pseudorange rows = %.3f   range-rate rows = %.4f\n', ...
            rms(S.nonLinZ(prRows) - Hx(prRows)), rms(S.nonLinZ(rrRows) - Hx(rrRows)));
        fprintf('  rms |nonLinZ - H2*xPrior| pseudorange rows = %.3f\n', rms(S.nonLinZ(prRows) - Hx2(prRows)));
    end

    % ---- 3. innovation against the truth, all rows -------------------------
    fprintf('\n-- 3. innovation against the true state error --\n');
    [xTrue, haveAtt] = trueState(S, n, POS, VEL, ATT);
    if any(isfinite(xTrue(POS)))
        for sgn = [1 -1]
            xt = sgn * xTrue;  xt(~isfinite(xt)) = 0;
            pred = S.measH * (xt - S.xPrior);
            resid = removeClock(S.innov - pred, S, prRows, rrRows);
            if sgn == 1, lab = 'x = mech - truth'; else, lab = 'x = truth - mech'; end
            fprintf('  %-18s rms residual: pseudorange %9.3f m   range-rate %8.4f m/s\n', ...
                lab, rms(resid(prRows)), rms(resid(rrRows)));
        end
        fprintf('  (per-constellation median removed = clock; the convention with the small residual is the host''s)\n');
    else
        fprintf('  mechPos / refPos missing: skipped\n');
    end

    % ---- 4. update split by rows -------------------------------------------
    fprintf('\n-- 4. update split by row type, x(1:6) --\n');
    if ~isempty(S.covarIn) && ~isempty(S.R)
        sets = {prRows, rrRows, [prRows rrRows baroRows]};
        labs = {'pseudorange only', 'range-rate only', 'all rows'};
        for j = 1:3
            r = sets{j};
            if isempty(r), continue; end
            Hs = S.measH(r, :);
            Ks = S.covarIn * Hs' / (Hs * S.covarIn * Hs' + S.R(r, r));
            xs = S.xPrior + Ks * S.innov(r);
            fprintf('  %-18s %s\n', labs{j}, vec2str(xs(1:6)));
        end
        fprintf('  %-18s %s\n', 'prior', vec2str(S.xPrior(1:6)));
        if any(isfinite(xTrue(POS)))
            fprintf('  %-18s %s\n', 'truth (mech-truth)', vec2str(xTrue(1:6)));
        end
    end

    % ---- 5. coupling ---------------------------------------------------------
    fprintf('\n-- 5. coupling in covarIn --\n');
    if ~isempty(S.covarIn)
        P = S.covarIn;
        d = sqrt(diag(P(1:9, 1:9)));
        C = P(1:9, 1:9) ./ (d * d');
        fprintf('  sigma(1:9) = %s\n', vec2str(d'));
        fprintf('  correlation pos-vel  N %.2f  E %.2f  D %.2f\n', C(1, 4), C(2, 5), C(3, 6));
        fprintf('  correlation pos-att  N-pitch %.2f  E-roll %.2f\n', C(1, 8), C(2, 7));
        fprintf('  lever P(pos,vel)/P(vel,vel) [s]:  N %.1f  E %.1f  D %.1f\n', ...
            P(1, 4) / P(4, 4), P(2, 5) / P(5, 5), P(3, 6) / P(6, 6));
    end

    % ---- 6. prior against truth ----------------------------------------------
    fprintf('\n-- 6. prior state against the true error (mech - truth) --\n');
    if any(isfinite(xTrue(POS)))
        fprintf('  state      prior        truth       ratio\n');
        lab9 = {'pN', 'pE', 'pD', 'vN', 'vE', 'vD', 'att1', 'att2', 'att3'};
        for i = 1:9
            if isfinite(xTrue(i))
                fprintf('  %-5s %12.4f %12.4f %9.2f\n', lab9{i}, S.xPrior(i), xTrue(i), S.xPrior(i) / xTrue(i));
            else
                fprintf('  %-5s %12.4f        (no reference)\n', lab9{i}, S.xPrior(i));
            end
        end
        if ~haveAtt
            fprintf('  attitude truth not available (mechAtt/refAtt missing)\n');
        end
        res.(fld(k)).xTrue = xTrue;
    end
end

% ---------------------------------------------------------------- 7. propagation across the gap
fprintf('\n\n############ 7. propagation across the gap (PRE -> FIRST) ############\n');
Sp = E{1}; Sf = E{2};
if ~isempty(Sf.accumPhi) && ~isempty(Sp.xPost) && ~isempty(Sf.statePropagated)
    xProp = Sf.accumPhi * Sp.xPost;
    fprintf('  state      accumPhi*xPost(PRE)   statePropagated(FIRST)   diff\n');
    lab9 = {'pN', 'pE', 'pD', 'vN', 'vE', 'vD', 'att1', 'att2', 'att3'};
    for i = 1:9
        fprintf('  %-5s %18.4f %22.4f %12.4f\n', lab9{i}, xProp(i), Sf.statePropagated(i), xProp(i) - Sf.statePropagated(i));
    end
    if ~isempty(Sp.covar) && ~isempty(Sf.covPropagated) && ~isempty(Sf.accumQ)
        Pprop = Sf.accumPhi * Sp.covar * Sf.accumPhi' + Sf.accumQ;
        fprintf('  sqrt diag P(1:6): propagated by hand %s\n', vec2str(sqrt(diag(Pprop(1:6, 1:6)))'));
        fprintf('                    host               %s\n', vec2str(sqrt(diag(Sf.covPropagated(1:6, 1:6)))'));
    end
else
    fprintf('  accumPhi / statePropagated / xPost(PRE) missing: skipped\n');
end

% ---------------------------------------------------------------- 8. attitude -> velocity -> position chain
fprintf('\n\n############ 8. attitude residual -> velocity -> position over the gap ############\n');
Sa = E{3};
dt = Sa.t - Sp.t;
if isfinite(dt) && dt > 0 && ~isempty(Sa.xPrior)
    att = Sa.xPrior(ATT);
    fprintf('  gap length dt = %.1f s\n', dt);
    fprintf('  prior attitude residual [rad] = %s   (%.3f, %.3f, %.3f deg)\n', vec2str(att'), att * 180 / pi);
    fprintf('  g*att*dt      [m/s] = %s    prior vel  = %s\n', vec2str((G * att * dt)'), vec2str(Sa.xPrior(VEL)'));
    fprintf('  0.5*g*att*dt^2 [m] = %s    prior pos  = %s\n', vec2str((0.5 * G * att * dt^2)'), vec2str(Sa.xPrior(POS)'));
    fprintf('  (roll tilts east velocity, pitch tilts north velocity: compare att(1) with vE/pE, att(2) with vN/pN, signs per host convention)\n');
    xt = trueState(Sa, numel(Sa.xPrior), POS, VEL, ATT);
    if any(isfinite(xt(POS)))
        drift = xt(POS);
        fprintf('  actual drift (mech-truth) [m] = %s\n', vec2str(drift'));
        fprintf('  tilt implied by the actual drift, 2*drift/(g*dt^2) [deg] = %s\n', vec2str((2 * drift / (G * dt^2) * 180 / pi)'));
    end
    if any(isfinite(xt(ATT)))
        fprintf('  true attitude error (mech-ref) [deg] = %s\n', vec2str((xt(ATT) * 180 / pi)'));
    end
else
    fprintf('  need the time tags of PRE and ACC (tRx or file name) and xPrior(ACC): skipped\n');
end
fprintf('\n=========================================================\n');

end

% ============================================================================
function [S] = loadEpoch(file)
D = load(file);
want = {'covarIn', 'measH', 'R', 'measZ', 'nonLinZ', 'innov', 'K', 'xUpdate', 'xPrior', 'xPost', 'covar', ...
        'rowMap', 'mechPos', 'mechVel', 'mechAtt', 'refPos', 'refPosEcef', 'refVel', 'refAtt', ...
        'statePropagated', 'covPropagated', 'accumPhi', 'accumQ', 'tRx', ...
        'satPos', 'satVel', 'rhoMeas', 'rangeRateMeas'};
S = struct();
missing = {};
for i = 1:numel(want)
    if isfield(D, want{i})
        S.(want{i}) = double(D.(want{i}));
    else
        S.(want{i}) = [];
        missing{end + 1} = want{i}; %#ok<AGROW>
    end
end
for f = {'xPrior', 'xPost', 'measZ', 'nonLinZ', 'innov', 'xUpdate', 'statePropagated', 'mechPos', 'mechVel', 'mechAtt', 'refPos', 'refVel', 'refAtt', 'refPosEcef', 'rhoMeas', 'rangeRateMeas'}
    S.(f{1}) = S.(f{1})(:);
end
S.file = file;
if ~isempty(S.tRx)
    S.t = S.tRx(1);
else
    tok = regexp(file, '(\d+\.\d+)', 'tokens', 'once');
    if isempty(tok), S.t = NaN; else, S.t = str2double(tok{1}); end
end
if ~isempty(missing)
    fprintf('  [%s] missing: %s\n', file, strjoin(missing, ' '));
end
end

function [prRows, rrRows, baroRows] = rowsOf(S)
m = size(S.measH, 1);
if ~isempty(S.rowMap) && size(S.rowMap, 2) >= 1
    typ = S.rowMap(:, 1);
    prRows = find(typ == 1)';  rrRows = find(typ == 2)';  baroRows = find(typ == 3)';
else
    % infer from H: pseudorange rows have unit position entries, range-rate rows unit velocity entries
    hp = sqrt(sum(S.measH(:, 1:3) .^ 2, 2));
    hv = sqrt(sum(S.measH(:, 4:6) .^ 2, 2));
    prRows   = find(hp > 0.5 & hv < 0.5)';
    rrRows   = find(hv > 0.5)';
    baroRows = setdiff(1:m, [prRows rrRows]);
end
end

function [xTrue, haveAtt] = trueState(S, n, POS, VEL, ATT)
% true error of the mechanisation in the state space: mech - truth (NED metres, m/s, rad)
xTrue = NaN(n, 1);
haveAtt = false;
if numel(S.mechPos) == 3 && numel(S.refPos) == 3
    L  = S.mechPos(1);  h = S.mechPos(3);
    a  = 6378137;  e2 = 6.69437999014e-3;
    RN = a * (1 - e2) / (1 - e2 * sin(L)^2)^1.5;
    RE = a / sqrt(1 - e2 * sin(L)^2);
    xTrue(POS(1)) = (S.mechPos(1) - S.refPos(1)) * (RN + h);
    xTrue(POS(2)) = (S.mechPos(2) - S.refPos(2)) * (RE + h) * cos(L);
    xTrue(POS(3)) = -(S.mechPos(3) - S.refPos(3));
end
if numel(S.mechVel) == 3 && numel(S.refVel) == 3
    xTrue(VEL) = S.mechVel - S.refVel;
end
if numel(S.mechAtt) == 3 && numel(S.refAtt) == 3
    d = S.mechAtt - S.refAtt;
    d(3) = atan2(sin(d(3)), cos(d(3)));
    xTrue(ATT) = d;
    haveAtt = true;
end
end

function [r] = removeClock(r, S, prRows, rrRows)
% subtract the per-constellation median over pseudorange rows and over range-rate rows
if ~isempty(S.rowMap) && size(S.rowMap, 2) >= 2
    con = S.rowMap(:, 2);
else
    con = ones(numel(r), 1);
end
for rows = {prRows, rrRows}
    rr = rows{1};
    for c = unique(con(rr))'
        idx = rr(con(rr) == c);
        r(idx) = r(idx) - median(r(idx));
    end
end
end

function [s] = vec2str(v)
s = sprintf('%11.4f', v);
end

function [r] = rms(v)
if isempty(v), r = NaN; else, r = sqrt(mean(v(:) .^ 2)); end
end

function [f] = fld(k)
names = {'pre', 'first', 'acc', 'ref'};
f = names{k};
end
