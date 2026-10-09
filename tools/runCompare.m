function [res] = runCompare(runA, runB, tOut, tReacq, aidedWin)
%RUNCOMPARE  Decompose an INS-only outage drift from the logs of one or two runs.
%   res = runCompare(runA, runB, tOut, tReacq, aidedWin)
%   res = runCompare(runA, [],   tOut, tReacq, aidedWin)      one run only
%
%   Post-processing only: nothing is re-run. Each run is a struct of the
%   arrays the host already logs (nav solution, truth, bias estimates):
%     t        [N x 1]  host time of the nav samples (s)
%     navPos   [N x 3]  nav  [lat lon h]   rad, rad, m
%     navVel   [N x 3]  nav  [vN vE vD]    m/s
%     navAtt   [N x 3]  nav  [roll pitch yaw] rad
%     biasEst  [N x 6]  IMU bias the mechanisation SUBTRACTS from the raw
%                       IMU: accel (m/s^2) then gyro (rad/s), body axes.
%                       With a split host this is the fed-back value plus
%                       the bias state. May be [K x 6] on the 2 Hz epochs
%                       tKf instead (held between updates), or a single row
%                       [1 x 6] holding the value at tOut.
%     tKf      [K x 1]  2 Hz update epochs (s), for biasEst or pDiag logged
%                       per update. Optional: without it the K records are
%                       spread evenly over the nav time span (valid when the
%                       2 Hz function logged one record per call)
%     tTrue    [M x 1]  truth time (s), may be a denser or offset time base
%     truePos  [M x 3]  truth [lat lon h],   trueAtt [M x 3] [roll pitch yaw]
%     trueVel  [M x 3]  truth NED velocity, optional: omitted, empty, or not
%                       matching d(truePos)/dt (another frame, e.g. ECEF) it
%                       is replaced by the velocity differenced from truePos
%     biasTrue [1 x 6]  simulated IMU biases (or [M x 6] if they vary)
%     pDiag    optional sigma source: [N x 22] diagonal of P per nav sample,
%                       or [K x 22] per 2 Hz epoch (K = numel(tKf)), or a cell
%                       {K x 1} holding per epoch either the 22 x 22 P or its
%                       diagonal, or a single 22 x 22 P (or its diagonal) taken
%                       at tOut. A field named pDig is accepted as pDiag.
%     name     char, optional
%   tOut     last GNSS-aided epoch before the outage (s)
%   tReacq   first GNSS epoch after the outage (s)
%   aidedWin [K x 2]  GNSS-aided windows [t1 t2] used for the truth-lag
%            estimate and the aided-accuracy figures, e.g. [100 760; 900 1600].
%            Optional: omitted or [], the 300 s before tOut are used.
%
%   What it prints, per run
%     0. time alignment. (a) truth arrays: trueVel and the yaw of trueAtt
%        are checked against the differenced truePos (a logged array shifted
%        in time, in another frame or with another convention shows here);
%        a shifted velocity is re-aligned, a shifted yaw only when it then
%        matches the course to a few mrad (sideslip dynamics also shift the
%        yaw against the course). (b) nav arrays: navVel against the
%        differenced navPos (a mechanised solution must agree between the
%        2 Hz jumps). (c) the lag of the nav solution behind the truth, found
%        separately on the position, the velocity and the yaw. The attitude
%        lag is the reference instant for all errors: a gyro-integrated
%        attitude cannot trail the truth, so its lag is the time-base offset
%        of the logged truth; what the position and velocity trail beyond it
%        is error of the nav solution (a position error proportional to the
%        velocity and a velocity error proportional to the acceleration are
%        the signature of GNSS measurements applied late). The regression of
%        the errors on the truth velocity / acceleration quantifies that.
%     1. aided accuracy: RMS position, velocity and attitude error, bias error
%     2. state error at tOut: velocity, tilt, bias, with sigma from pDiag
%     3. outage trace: the slope over the first 10 s (must match the tOut
%        velocity error, or the truth or the lag is wrong), then the fit
%        err = a + b x + c x^2 + d x^3 per NED axis, x from tOut; the terms
%        b T, c T^2, d T^3 at T = tReacq - tOut say whether the drift is
%        velocity (t), bias/tilt (t^2) or gyro (t^3)
%     4. state budget: the same drift predicted from the tOut errors,
%        including the yaw error acting on the manoeuvres of the outage
%   and, with two runs, the difference of the tOut state errors in sigmas.
%
%   Conventions used in section 4 (if a predicted term has the right size
%   and the wrong sign, the host uses the opposite convention for that
%   quantity; the magnitudes stand):
%     tilt psi from C_nav * C_true' = I - [psi x], i.e. for small errors
%     psi = -C_b^n * (navEuler - trueEuler); specific-force error
%     (g psi_E, -g psi_N, 0), so a nose-up pitch error drifts south;
%     accel-bias error enters as -(biasEst - biasTrue) on the corrected IMU;
%     yaw error psi_D rotates the integrated specific force, i.e. the part
%     of the true displacement not explained by the initial velocity:
%     -psi x (dp - v0 T) with psi = [0 0 psi_D]. The tilt and yaw rows use
%     the exact rotation (C_nav C_true' - I), the gyro row psi_dot = C db.
%   Result: res.run1 (and res.run2) with tau (reference lag), tauP, tauV,
%   tauA, sV, sA, sNV, kP, r2P, kV, r2V, velSource, rmsPos, rmsVel, meanVel,
%   rmsAtt, tOut, T, posErrOut, dv, slope10, attErr, psi, db, sig, drift,
%   coef, velTerm, biasTerm, tiltTerm, gyroTerm, yawTerm and the traces
%   posErrTrace, velErrTrace, attErrTrace, biasErrTrace on t. sig is the
%   square root of the P diagonal (22 x 1) at the last update <= tOut; it
%   is NaN, and the sigma columns print n/a, when pDiag is missing or its
%   records cannot be matched to tKf (a NOTE says so).
%   MATLAB desktop tool, not codegen. Runs under Octave.

G = 9.80665;
if nargin < 5 || isempty(aidedWin), aidedWin = [tOut - 300, tOut]; end
runs = {runA};
if ~isempty(runB), runs{end+1} = runB; end
res = struct();
for k = 1:numel(runs)
    r = runs{k};
    if ~isstruct(r) || numel(r) ~= 1
        error('runCompare:runStruct', ['run %d must be a single struct; it is a %d x %d %s. A cell (pDiag) passed ' ...
            'to struct() makes a struct array: assign the fields one by one, or write struct(..., ''pDiag'', {pCell})'], ...
            k, size(r, 1), size(r, 2), class(r));
    end
    if ~isfield(r, 'name') || isempty(r.name), r.name = sprintf('run %d', k); end
    fprintf('\n==================== %s ====================\n', r.name);
    if ~isfield(r, 'pDiag') && isfield(r, 'pDig')                                  % common short name
        r.pDiag = r.pDig;  r = rmfield(r, 'pDig');
        fprintf('   NOTE: field pDig used as pDiag\n');
    end
    extra = setdiff(fieldnames(r), {'t', 'navPos', 'navVel', 'navAtt', 'biasEst', 'tKf', 'tTrue', 'truePos', 'trueVel', ...
        'trueAtt', 'biasTrue', 'pDiag', 'name'});
    if ~isempty(extra)
        fprintf('   NOTE: run fields the tool does not use (misspelt?): %s\n', strjoin(extra(:)', ', '));
    end
    res.(sprintf('run%d', k)) = analyseRun(r, tOut, tReacq, aidedWin, G);
end
if numel(runs) == 2
    compareRuns(res.run1, res.run2, runs{1}.name, runs{2}.name);
end
end

%------------------------------------------------------------------------
function [o] = analyseRun(r, tOut, tReacq, aidedWin, G)
t = r.t(:);
N = numel(t);
yawT = r.trueAtt(:, 3);  okY = isfinite(yawT);  yawT(okY) = unwrap(yawT(okY));  r.trueAtt(:, 3) = yawT;   % continuous across +-pi for interpolation
rA = NaN;

%% quantities logged per 2 Hz update: bias held between updates, P taken at the last update <= tOut
P = [];
if isfield(r, 'pDiag') && ~isempty(r.pDiag), P = pDiagMatrix(r.pDiag); end
if size(r.biasEst, 1) == 1, r.biasEst = repmat(r.biasEst, N, 1); end          % value at tOut only
if size(P, 1) == 1, P = repmat(P, N, 1); end                                   % value at tOut only
nKf = 0;
if isfield(r, 'tKf') && ~isempty(r.tKf), nKf = numel(r.tKf); end
if ~isempty(P) && size(P, 1) ~= N && size(P, 1) ~= nKf && any(size(P, 2) == [N, nKf, size(r.biasEst, 1)])
    P = P';                                                                    % given as states x records
end
K = max(size(r.biasEst, 1) * (size(r.biasEst, 1) ~= N), size(P, 1) * (size(P, 1) ~= N));
if K > 0 && (~isfield(r, 'tKf') || isempty(r.tKf))
    r.tKf = t(1) + (t(end) - t(1)) * (0:K-1)' / (K - 1);
    fprintf('   NOTE: tKf not given; %d 2 Hz records spread evenly over %.2f..%.2f s, %.3f s apart (expect 0.5)\n', ...
        K, t(1), t(end), (t(end) - t(1)) / (K - 1));
end
if isfield(r, 'tKf') && ~isempty(r.tKf)
    tKf = r.tKf(:);
    if size(r.biasEst, 1) == numel(tKf) && numel(tKf) ~= N
        r.biasEst = interpRows(tKf, r.biasEst, t, 'previous');
    end
    if ~isempty(P) && size(P, 1) == numel(tKf)
        iKf = find(tKf <= tOut, 1, 'last');
        if isempty(iKf)
            fprintf('   NOTE: no tKf epoch at or before tOut = %.2f (tKf runs %.2f..%.2f): sigmas not available\n', tOut, tKf(1), tKf(end));
            P = [];
        else
            P = repmat(P(iKf, :), N, 1);
        end
    elseif ~isempty(P) && size(P, 1) ~= N
        fprintf('   NOTE: pDiag has %d records, tKf %d epochs, the nav %d samples: sigmas not available (log pDiag and tKf at the same place)\n', ...
            size(P, 1), numel(tKf), N);
        P = [];
    end
end
if size(r.biasEst, 1) ~= N
    error('runCompare:bias', 'biasEst has %d rows; give it per nav sample (%d) or per 2 Hz epoch with tKf', size(r.biasEst, 1), N);
end
inAided = false(N, 1);
for w = 1:size(aidedWin, 1)
    inAided = inAided | (t >= aidedWin(w, 1) & t <= aidedWin(w, 2));
end

%% 0. time alignment: nav(t) is compared with the truth at t + lag
tauGrid = -1.5:0.01:1.5;
tT = r.tTrue(:);
tA = t(inAided);
inT = false(numel(tT), 1);
for w = 1:size(aidedWin, 1)
    inT = inT | (tT >= aidedWin(w, 1) - 2 & tT <= aidedWin(w, 2) + 2);
end
fprintf('0. time alignment\n');

%% (a) truth arrays against the differenced truth position
velPos = velFromPos(tT, r.truePos);                           % NED velocity differenced from the truth position
velSource = 'position';
sV = 0;  rV = NaN;
if isfield(r, 'trueVel') && ~isempty(r.trueVel)
    [sV, ~, ~] = lagScan(tauGrid, @(s) stdRows(interpRows(tT, r.trueVel, tT(inT) + s) - velPos(inT, :), 1:3, 0));
    rV = rmsRows(interpRows(tT, r.trueVel, tT(inT) + sV) - velPos(inT, :), 1:3);
    if rV < 1 && abs(sV) < 1.49
        velSource = 'logged';
        fprintf('   truth arrays: trueVel(t %+.2f s) = d(truePos)/dt(t), residual %.3f m/s rms', sV, rV);
        if rV > 0.1
            fprintf('\n   (truePos carries ~%.0f mm of sample-to-sample jitter; trueVel is the smooth source and gives the course and the acceleration)\n', ...
                rV * 0.02 / sqrt(2) * 1e3);
        else
            fprintf('\n');
        end
    elseif rV < 1
        sV = 0;
        fprintf('   WARNING: trueVel matches d(truePos)/dt only at the edge of the scan (%+.2f s): not trusted, the velocity\n', sV);
        fprintf('            differenced from truePos is used instead\n');
    else
        sV = 0;
        fprintf('   WARNING: trueVel differs from d(truePos)/dt by %.2f m/s rms at its best shift (another frame or unit?);\n', rV);
        fprintf('            the velocity differenced from truePos is used instead\n');
    end
end
if strcmp(velSource, 'logged'), trueVel = interpRows(tT, r.trueVel, tT + sV); else, trueVel = velPos; end
velRef = trueVel;                                             % velocity truth for the course, the acceleration and the moving mask
sA = 0;  rA = NaN;  rA0 = NaN;
spdT = hypot(velRef(:, 1), velRef(:, 2));
movT = inT & spdT > 3 & all(isfinite(velRef), 2);
if any(movT)
    cogT = atan2(velRef(movT, 2), velRef(movT, 1));
    [sA, ~, ~] = lagScan(tauGrid, @(s) yawStd(interpRows(tT, r.trueAtt(:, 3), tT(movT) + s), cogT));
    rA0 = yawRms(interpRows(tT, r.trueAtt(:, 3), tT(movT)), cogT);
    rA  = yawRms(interpRows(tT, r.trueAtt(:, 3), tT(movT) + sA), cogT);
    if abs(sA) >= 1.49 || abs(sA) <= 0.05
        sA = 0;
        fprintf('   truth arrays: trueAtt yaw = truth course over ground to %.1f mrad rms at zero shift (sideslip, or a convention if large)\n', rA0 * 1e3);
    elseif rA < 0.005
        fprintf('   truth arrays: trueAtt yaw(t %+.2f s) = truth course over ground(t), residual %.1f mrad rms: a logging offset, re-aligned\n', sA, rA * 1e3);
    else
        fprintf('   truth arrays: trueAtt yaw(t %+.2f s) = truth course over ground(t) best, %.1f mrad rms left (%.1f at zero shift): yaw and\n', sA, rA * 1e3, rA0 * 1e3);
        fprintf('   course differ by sideslip-like dynamics, so the shift is not applied (a logging offset would leave a few mrad); the attitude\n');
        fprintf('   lag below is uncertain by this much\n');
        sA = 0;
    end
end
trueAtt = interpRows(tT, r.trueAtt, tT + sA);
if abs(sV) > 0.05 || abs(sA) > 0.05
    fprintf('   (the shifted truth arrays are re-aligned to the truth position before anything below)\n');
end

%% (b) nav arrays against the differenced nav position
velNav = velFromPos(t, r.navPos);                             % d(navPos)/dt on the nav time base
[sNV, ~, ~] = lagScan(tauGrid, @(s) stdRows(interpRows(t, r.navVel, tA + s) - velNav(inAided, :), 1:2, 0.1));
rNV = stdRows(interpRows(t, r.navVel, tA + sNV) - velNav(inAided, :), 1:2, 0.1);
spdN = hypot(velNav(:, 1), velNav(:, 2));
movN = inAided & spdN > 3 & all(isfinite(velNav), 2);
rNA = NaN;
if any(movN)
    rNA = yawRmsTrim(r.navAtt(movN, 3), atan2(velNav(movN, 2), velNav(movN, 1)), 0.1);
end
fprintf('   nav arrays:   navVel(t %+.2f s) = d(navPos)/dt(t), residual %.3f m/s std (2 Hz jumps excluded); navAtt yaw minus the nav course %.1f mrad rms\n', ...
    sNV, rNV, rNA * 1e3);
if abs(sNV) > 0.05
    fprintf('   WARNING: a mechanised solution has navVel = d(navPos)/dt; a shift means the two nav arrays were logged at different points of the host step\n');
end

%% (c) lag of the nav solution behind the truth, per quantity
navPosA = r.navPos(inAided, :);
[tauP, pBest, p0] = lagScan(tauGrid, @(s) stdRows(nedDelta(navPosA, interpRows(tT, r.truePos, tA + s)), 1:2, 0));
[tauV, vBest, v0] = lagScan(tauGrid, @(s) stdRows(r.navVel(inAided, :) - interpRows(tT, trueVel, tA + s), 1:2, 0));
tauVp = NaN;  vpBest = NaN;
if strcmp(velSource, 'logged')                                % the nav against the other velocity truth, for comparison
    [tauVp, vpBest, ~] = lagScan(tauGrid, @(s) stdRows(r.navVel(inAided, :) - interpRows(tT, velPos, tA + s), 1:2, 0));
end
vqP = interpRows(tT, velRef, tA + tauP);
mov = all(isfinite(vqP), 2) & hypot(vqP(:, 1), vqP(:, 2)) > 3;      % moving aided samples: turns make the yaw lag visible
tauA = tauP;  aBest = NaN;  a0 = NaN;  aAtP = NaN;
if any(mov)
    yawNav = r.navAtt(inAided, 3);  yawNav = yawNav(mov);
    tM = tA(mov);
    yawCost = @(s) yawStd(yawNav, interpRows(tT, trueAtt(:, 3), tM + s));
    [tauA, aBest, a0] = lagScan(tauGrid, yawCost);
    aAtP = yawCost(tauP);
end
haveA = any(mov) && abs(tauA) < 1.49 && ~(min(rA, rA0) > 0.3);    % an attitude lag from a yaw consistent with the course
fprintf('   lag of the nav solution behind the truth (nav(t) = truth(t + lag)): position %+.2f s, velocity %+.2f s, attitude %+.2f s\n', tauP, tauV, tauA);
fprintf('   (aided error std at its own lag / at zero lag: position %.2f / %.2f m, velocity %.3f / %.3f m/s, yaw %.1f / %.1f mrad)\n', ...
    pBest, p0, vBest, v0, aBest * 1e3, a0 * 1e3);
if strcmp(velSource, 'logged')
    fprintf('   nav velocity error std against trueVel %.3f m/s (lag %+.2f s), against d(truePos)/dt %.3f m/s (lag %+.2f s)\n', vBest, tauV, vpBest, tauVp);
end
if ~haveA
    tau = tauP;
    fprintf('   reference instant for the errors below: the position lag (no usable attitude lag: no moving samples, truth yaw off the\n');
    fprintf('   course, or at the scan edge)\n');
elseif abs(tauA - tauP) <= 0.05
    tau = tauP;
    fprintf('   position and attitude lags agree: the logged truth is %+.2f s off the nav time base and the nav position does not trail\n', tauP);
    fprintf('   the truth; this lag is the reference instant for the errors below\n');
elseif aBest < 0.9 * aAtP
    tau = tauA;
    fprintf('   reference instant for the errors below: the attitude lag. A gyro-integrated attitude cannot trail the truth, so its lag is the\n');
    fprintf('   time-base offset of the logged truth; what the position and velocity trail beyond it is error of the nav solution.\n');
else
    tau = tauP;
    fprintf('   WARNING: the yaw error does not prefer its own lag over the position lag (std %.1f vs %.1f mrad): position lag used as the\n', aBest * 1e3, aAtP * 1e3);
    fprintf('            reference instant for the errors below\n');
end
edgeName = {'position lag', 'velocity lag', 'attitude lag', 'trueVel shift', 'navVel shift'};
edgeVal  = [tauP, tauV, tauA, sV, sNV];
for i = find(abs(edgeVal) >= 1.49)
    fprintf('   WARNING: %s at the edge of the scan (%+.2f s): extend tauGrid or check the time bases\n', edgeName{i}, edgeVal(i));
end
if haveA, tauK = tauA;  kName = 'attitude lag';  else, tauK = tauP;  kName = 'position lag';  end
tqK = tA + tauK;
eP  = nedDelta(navPosA, interpRows(tT, r.truePos, tqK));
eV  = r.navVel(inAided, :) - interpRows(tT, trueVel, tqK);
vT  = interpRows(tT, velRef, tqK);
aT  = interpRows(tT, diffRows(tT, velRef, 5), tqK);              % truth acceleration, 0.1 s stencil
[kP, r2P] = fitK(eP(:, 1:2), vT(:, 1:2));
[kV, r2V] = fitK(eV(:, 1:2), aT(:, 1:2));
fprintf('   at the %s: pos err = %+.3f s x truth velocity (%.0f%% of the horizontal error variance), vel err = %+.3f s x truth acceleration (%.0f%%)\n', ...
    kName, kP, 100 * r2P, kV, 100 * r2V);
if haveA && (abs(tauP - tauA) > 0.05 || abs(tauV - tauA) > 0.05)
    fprintf('   NOTE: relative to the attitude lag the nav position trails by %+.2f s and the velocity by %+.2f s. With consistent truth and\n', tauA - tauP, tauA - tauV);
    fprintf('         nav arrays (above) that is not a logging artefact: the GNSS measurements are applied late relative to the IMU\n');
    fprintf('         (pseudoranges ~%.2f s, Doppler ~%.2f s), or the truth used to build them is taken earlier than the nav epoch.\n', tauA - tauP, tauA - tauV);
end

%% error traces at the reference instant
tq      = t + tau;
posErr  = posErrNed(r, tq, true(N, 1));
velErr  = r.navVel - interpRows(tT, trueVel, tq);
attTrue = interpRows(tT, trueAtt, tq);
attErr  = r.navAtt - attTrue;
attErr(:, 3) = atan2(sin(attErr(:, 3)), cos(attErr(:, 3)));
if size(r.biasTrue, 1) > 1
    biasTrue = interpRows(tT, r.biasTrue, tq);
else
    biasTrue = repmat(r.biasTrue(:)', N, 1);
end
biasErr = r.biasEst - biasTrue;

%% 1. aided accuracy
ok = inAided & all(isfinite(posErr), 2) & all(isfinite(velErr), 2);
rmsPos = sqrt(mean(posErr(ok, :).^2, 1));
rmsVel = sqrt(mean(velErr(ok, :).^2, 1));
meanPos = mean(posErr(ok, :), 1);
meanVel = mean(velErr(ok, :), 1);
okA = ok & all(isfinite(attErr), 2);
rmsAtt = sqrt(mean(attErr(okA, :).^2, 1));
fprintf('1. aided windows: pos rms N %.2f E %.2f D %.2f m (mean N %.2f E %.2f D %.2f)\n', rmsPos, meanPos);
fprintf('   vel rms N %.3f E %.3f D %.3f m/s (mean N %+.3f E %+.3f D %+.3f), att rms roll %.2f pitch %.2f yaw %.2f mrad\n', ...
    rmsVel, meanVel, rmsAtt * 1e3);
lastWin = inAided & t <= tOut & t >= tOut - 60;
if any(lastWin)
    bm = mean(biasErr(lastWin, :), 1);
    fprintf('   bias error, mean over the last 60 s before tOut: accel %+.3f %+.3f %+.3f mg, gyro %+.3f %+.3f %+.3f deg/h\n', ...
        bm(1:3) / (G * 1e-3), bm(4:6) * 180 / pi * 3600);
end

%% 2. state error at tOut
iOut = find(t <= tOut, 1, 'last');
iEnd = find(t < tReacq, 1, 'last');            % last INS-only sample before the first update
T    = t(iEnd) - t(iOut);
dv   = velErr(iOut, :)';
db   = biasErr(iOut, :)';
Cnav  = dcmFromEuler(r.navAtt(iOut, :));
Ctrue = dcmFromEuler(attTrue(iOut, :));
dC  = Cnav * Ctrue';
psi = [dC(2,3) - dC(3,2); dC(3,1) - dC(1,3); dC(1,2) - dC(2,1)] / 2;   % I - [psi x]
sig = NaN(22, 1);
if ~isempty(P) && size(P, 1) == N
    sig = sqrt(max(P(iOut, :), 0))';
end
if ~any(isfinite(sig))
    fprintf('   (no sigmas: pDiag missing or not matched to the epochs, see the NOTEs above)\n');
end
fprintf('2. state error at t = %.2f (truth lag applied), outage T = %.1f s to t = %.2f\n', t(iOut), T, t(iEnd));
fprintf('   pos err   N %+8.2f E %+8.2f D %+8.2f m\n', posErr(iOut, :));
fprintf('   vel err   N %+8.3f E %+8.3f D %+8.3f m/s     sigma %s\n', dv, sigStr(sig(4:6)));
fprintf('   euler err roll %+7.3f pitch %+7.3f yaw %+7.3f mrad   sigma(mrad) %s\n', attErr(iOut, :) * 1e3, sigStr(sig(7:9) * 1e3));
fprintf('   tilt psi  N %+7.3f E %+7.3f D %+7.3f mrad (from the DCMs)\n', psi * 1e3);
fprintf('   accel bias err %+7.3f %+7.3f %+7.3f mg      sigma(mg) %s\n', db(1:3) / (G * 1e-3), sigStr(sig(10:12) / (G * 1e-3)));
fprintf('   gyro  bias err %+7.3f %+7.3f %+7.3f deg/h   sigma(deg/h) %s\n', db(4:6) * 180 / pi * 3600, sigStr(sig(13:15) * 180 / pi * 3600));

%% 3. outage trace fit
sel = (t >= t(iOut)) & (t <= t(iEnd)) & all(isfinite(posErr), 2);
x   = t(sel) - t(iOut);
drift = posErr(iEnd, :) - posErr(iOut, :);
i10 = min(find(t <= t(iOut) + 10, 1, 'last'), iEnd);
slope10 = (posErr(i10, :) - posErr(iOut, :)) / (t(i10) - t(iOut));
[jump, iJ] = max(sqrt(sum(diff(posErr(iOut:i10, :), 1, 1).^2, 2)));
fprintf('3. outage drift observed: N %+8.2f E %+8.2f D %+8.2f m\n', drift);
aErr = Cnav * (-db(1:3)) + G * [psi(2); -psi(1); 0];        % bias + tilt acceleration error of section 2
dvImplied = slope10 - 0.5 * aErr' * (t(i10) - t(iOut));
fprintf('   slope over the first %.0f s of the outage: N %+.3f E %+.3f D %+.3f m/s = vel err at tOut + %.0f s x the bias+tilt acceleration\n', ...
    t(i10) - t(iOut), slope10, 0.5 * (t(i10) - t(iOut)));
fprintf('   implied vel err at tOut: N %+.3f E %+.3f D %+.3f m/s; a gap of more than ~0.1 m/s to section 2 means the truth velocity or its lag is wrong\n', ...
    dvImplied);
if jump > 0.5
    fprintf('   WARNING: the position error jumps by %.2f m at t = %.2f inside those %.0f s (a solution switch?); the slope above includes it\n', ...
        jump, t(iOut + iJ), t(i10) - t(iOut));
end
fprintf('   fit err = a + b x + c x^2 + d x^3 over %d samples, contributions at T = %.0f s:\n', nnz(sel), T);
fprintf('   axis      b T (vel)    c T^2 (bias/tilt)    d T^3 (gyro)    fit rms   b [m/s]\n');
coef = zeros(3, 4);
axName = 'NED';
for ax = 1:3
    p = polyfit(x, posErr(sel, ax), 3);              % [d c b a]
    coef(ax, :) = p;
    fitRms = sqrt(mean((polyval(p, x) - posErr(sel, ax)).^2));
    fprintf('   %s     %+9.2f     %+9.2f            %+9.2f       %6.2f    %+.3f\n', ...
        axName(ax), p(3) * T, p(2) * T^2, p(1) * T^3, fitRms, p(3));
end
fprintf('   (b should agree with the vel err of section 2; a large gap means the lag or the logged solution is wrong)\n');

%% 4. budget from the tOut state
velTerm  = dv * T;
biasTerm = 0.5 * (Cnav * (-db(1:3))) * T^2;                 % -(est - true) is the residual on the corrected IMU
tiltTerm = 0.5 * (dC - eye(3)) * [0; 0; -G] * T^2;          % exact: (C_nav C_true' - I) acting on the gravity reaction; -psi x f to first order
psiDot   = Cnav * db(4:6);                                  % the nav integrates omega_true - (est - true): psi_dot = +C db for C_nav = (I - [psi x]) C_true
gyroTerm = G * [psiDot(2); -psiDot(1); 0] * T^3 / 6;
dpTrue   = nedDelta(interpRows(tT, r.truePos, tq(iEnd)), interpRows(tT, r.truePos, tq(iOut)))';
v0       = interpRows(tT, trueVel, tq(iOut))';
u        = dpTrue - v0 * T;                                 % displacement from the integrated specific force
yawTerm  = (dC - eye(3)) * u;                               % exact; [psi_D u_E, -psi_D u_N, 0] to first order
total    = velTerm + biasTerm + tiltTerm + gyroTerm + yawTerm;
fprintf('4. budget from the tOut state (m):   N         E         D\n');
fprintf('   velocity  dv T           %+8.2f  %+8.2f  %+8.2f\n', velTerm);
fprintf('   accel bias 1/2 C db T^2  %+8.2f  %+8.2f  %+8.2f\n', biasTerm);
fprintf('   tilt   (dC-I) g T^2/2    %+8.2f  %+8.2f  %+8.2f\n', tiltTerm);
fprintf('   gyro bias g psidot T^3/6 %+8.2f  %+8.2f  %+8.2f\n', gyroTerm);
fprintf('   yaw   (dC-I)(dp - v0 T)  %+8.2f  %+8.2f  %+8.2f   (manoeuvre displacement |dp - v0 T| = %.0f m)\n', yawTerm, norm(u(1:2)));
fprintf('   sum                      %+8.2f  %+8.2f  %+8.2f\n', total);
fprintf('   observed                 %+8.2f  %+8.2f  %+8.2f\n', drift);
bt = r.biasTrue(1, 1:3);
fprintf('   scale: true accel bias uncompensated over T = %.0f m, applied with the wrong sign = %.0f m\n', ...
    0.5 * norm(bt) * T^2, norm(bt) * T^2);

o = struct('name', r.name, 'tau', tau, 'tauP', tauP, 'tauV', tauV, 'tauVp', tauVp, 'tauA', tauA, 'sV', sV, 'sA', sA, 'sNV', sNV, ...
    'kP', kP, 'r2P', r2P, 'kV', kV, 'r2V', r2V, 'velSource', velSource, 'rmsPos', rmsPos, ...
    'rmsVel', rmsVel, 'meanVel', meanVel, 'rmsAtt', rmsAtt, 'tOut', t(iOut), 'T', T, ...
    'posErrOut', posErr(iOut, :), 'dv', dv, 'slope10', slope10', 'dvImplied', dvImplied', 'attErr', attErr(iOut, :)', 'psi', psi, 'db', db, 'sig', sig, ...
    'drift', drift, 'coef', coef, 'velTerm', velTerm, 'biasTerm', biasTerm, 'tiltTerm', tiltTerm, ...
    'gyroTerm', gyroTerm, 'yawTerm', yawTerm, 'posErrTrace', posErr, 'velErrTrace', velErr, 'attErrTrace', attErr, 'biasErrTrace', biasErr, 't', t);
end

%------------------------------------------------------------------------
function compareRuns(a, b, nameA, nameB)
fprintf('\n==================== %s minus %s at tOut ====================\n', nameB, nameA);
sig = b.sig;  sig(~isfinite(sig)) = a.sig(~isfinite(sig));
d = [b.dv - a.dv; b.psi - a.psi; b.db - a.db];
s = [sig(4:6); sig(7:9); sig(10:15)];
labels = {'vel N', 'vel E', 'vel D', 'psi N', 'psi E', 'psi D', 'ba X', 'ba Y', 'ba Z', 'bg X', 'bg Y', 'bg Z'};
units  = [1 1 1, 1e3 1e3 1e3, 1/(9.80665e-3) 1/(9.80665e-3) 1/(9.80665e-3), 180/pi*3600 180/pi*3600 180/pi*3600];
unitNames = {'m/s', 'm/s', 'm/s', 'mrad', 'mrad', 'mrad', 'mg', 'mg', 'mg', 'deg/h', 'deg/h', 'deg/h'};
fprintf('   quantity   difference        in sigmas of %s\n', nameB);
for i = 1:12
    if isfinite(s(i)) && s(i) > 0
        fprintf('   %-7s  %+10.4f %-6s  %+7.2f\n', labels{i}, d(i) * units(i), unitNames{i}, d(i) / s(i));
    else
        fprintf('   %-7s  %+10.4f %-6s      n/a\n', labels{i}, d(i) * units(i), unitNames{i});
    end
end
fprintf('   drift observed   %s: N %+7.2f E %+7.2f | %s: N %+7.2f E %+7.2f\n', nameA, a.drift(1:2), nameB, b.drift(1:2));
fprintf('   drift difference %s minus %s, predicted from the state differences (truth offsets cancel):\n', nameB, nameA);
fprintf('                           N         E\n');
fprintf('   velocity         %+8.2f  %+8.2f\n', b.velTerm(1:2) - a.velTerm(1:2));
fprintf('   accel bias       %+8.2f  %+8.2f\n', b.biasTerm(1:2) - a.biasTerm(1:2));
fprintf('   tilt             %+8.2f  %+8.2f\n', b.tiltTerm(1:2) - a.tiltTerm(1:2));
fprintf('   gyro bias        %+8.2f  %+8.2f\n', b.gyroTerm(1:2) - a.gyroTerm(1:2));
fprintf('   yaw              %+8.2f  %+8.2f\n', b.yawTerm(1:2) - a.yawTerm(1:2));
fprintf('   sum              %+8.2f  %+8.2f\n', (b.velTerm(1:2) + b.biasTerm(1:2) + b.tiltTerm(1:2) + b.gyroTerm(1:2) + b.yawTerm(1:2)) ...
    - (a.velTerm(1:2) + a.biasTerm(1:2) + a.tiltTerm(1:2) + a.gyroTerm(1:2) + a.yawTerm(1:2)));
fprintf('   observed         %+8.2f  %+8.2f\n', b.drift(1:2) - a.drift(1:2));
fprintf('   aided pos rms    %s: N %5.2f E %5.2f | %s: N %5.2f E %5.2f m\n', nameA, a.rmsPos(1:2), nameB, b.rmsPos(1:2));
fprintf('   aided vel rms    %s: N %5.3f E %5.3f | %s: N %5.3f E %5.3f m/s\n', nameA, a.rmsVel(1:2), nameB, b.rmsVel(1:2));
fprintf('   A difference of a few sigmas or less on every row is two draws of the same filter; many sigmas on the\n');
fprintf('   bias or tilt rows means the two filters disagree beyond their own uncertainty, i.e. at least one is inconsistent.\n');
end

%------------------------------------------------------------------------
function [e] = posErrNed(r, tq, mask)
% nav(mask) minus truth interpolated at tq, in NED metres
e = nedDelta(r.navPos(mask, :), interpRows(r.tTrue, r.truePos, tq));
end

function [e] = nedDelta(p, q)
% p minus q, both [lat lon h] rows, in NED metres (WGS84 radii at p)
L  = p(:, 1);  h = p(:, 3);
a  = 6378137;  e2 = 6.69437999014e-3;
RN = a * (1 - e2) ./ (1 - e2 * sin(L).^2).^1.5;
RE = a ./ sqrt(1 - e2 * sin(L).^2);
e = [(p(:, 1) - q(:, 1)) .* (RN + h), ...
     (p(:, 2) - q(:, 2)) .* (RE + h) .* cos(L), ...
    -(p(:, 3) - q(:, 3))];
end

function [tau, best, at0] = lagScan(tauGrid, cost)
% the lag on tauGrid that minimises cost(lag); best = that minimum, at0 = cost at zero lag.
% The costs are standard deviations (mean removed), so a constant offset between nav and
% truth (lever arm, GNSS bias, a real filter error) does not bias the lag estimate.
c = NaN(size(tauGrid));
for i = 1:numel(tauGrid), c(i) = cost(tauGrid(i)); end
[best, iBest] = min(c);
tau = tauGrid(iBest);
[~, i0] = min(abs(tauGrid));
at0 = c(i0);
end

function [v] = rmsRows(e, cols)
% rms of the row norms over the finite rows, using the given columns
ok = all(isfinite(e), 2);
v  = sqrt(mean(sum(e(ok, cols).^2, 2)));
end

function [v] = stdRows(e, cols, trim)
% pooled standard deviation (mean of each column removed) over the finite rows, using the given columns;
% trim > 0 drops that fraction of the rows with the largest norm first (the 2 Hz jumps of a nav solution)
ok = all(isfinite(e), 2);
e  = e(ok, cols);
e  = e - repmat(mean(e, 1), size(e, 1), 1);
if trim > 0 && size(e, 1) > 10
    n  = sqrt(sum(e.^2, 2));
    ns = sort(n);
    e  = e(n <= ns(max(1, floor(numel(ns) * (1 - trim)))), :);
    e  = e - repmat(mean(e, 1), size(e, 1), 1);
end
v  = sqrt(mean(sum(e.^2, 2)));
end

function [v] = yawRmsTrim(a, b, trim)
% rms of the wrapped difference a - b after dropping the fraction trim with the largest |difference|
d  = a(:) - b(:);
d  = atan2(sin(d), cos(d));
d  = abs(d(isfinite(d)));
ds = sort(d);
d  = d(d <= ds(max(1, floor(numel(ds) * (1 - trim)))));
v  = sqrt(mean(d.^2));
end

function [k, r2] = fitK(e, x)
% least-squares gain k in e = k x, pooled over the columns, means removed; r2 = fraction of the variance of e explained
ok = all(isfinite([e x]), 2);
e  = e(ok, :);  x = x(ok, :);
e  = e - repmat(mean(e, 1), size(e, 1), 1);
x  = x - repmat(mean(x, 1), size(x, 1), 1);
k  = sum(sum(e .* x)) / sum(sum(x .* x));
r2 = 1 - sum(sum((e - k * x).^2)) / sum(sum(e.^2));
end

function [d] = diffRows(tt, Y, h)
% central difference of each column with a stencil of h samples on each side (one-sided at the ends)
tt = tt(:);
M  = numel(tt);
i0 = max((1:M) - h, 1);  i1 = min((1:M) + h, M);
d  = (Y(i1, :) - Y(i0, :)) ./ repmat(tt(i1) - tt(i0), 1, size(Y, 2));
end

function [v] = yawRms(a, b)
% rms of the wrapped difference a - b over the finite samples
d  = a(:) - b(:);
d  = atan2(sin(d), cos(d));
v  = sqrt(mean(d(isfinite(d)).^2));
end

function [v] = yawStd(a, b)
% standard deviation of the wrapped difference a - b over the finite samples
d  = a(:) - b(:);
d  = atan2(sin(d), cos(d));
d  = d(isfinite(d));
v  = sqrt(mean((d - mean(d)).^2));
end

function [v] = velFromPos(tt, pos)
% NED velocity by central differences of truth [lat lon h] (one-sided at the ends)
tt = tt(:);
M  = numel(tt);
i0 = [1, 1:M-1];  i1 = [2:M, M];
dt = tt(i1) - tt(i0);
L  = pos(:, 1);  h = pos(:, 3);
a  = 6378137;  e2 = 6.69437999014e-3;
RN = a * (1 - e2) ./ (1 - e2 * sin(L).^2).^1.5;
RE = a ./ sqrt(1 - e2 * sin(L).^2);
v = [(pos(i1, 1) - pos(i0, 1)) ./ dt .* (RN + h), ...
     (pos(i1, 2) - pos(i0, 2)) ./ dt .* (RE + h) .* cos(L), ...
    -(pos(i1, 3) - pos(i0, 3)) ./ dt];
end

function [y] = interpRows(tt, Y, tq, method)
if nargin < 4, method = 'linear'; end
y = NaN(numel(tq), size(Y, 2));
for c = 1:size(Y, 2)
    y(:, c) = interp1(tt(:), Y(:, c), tq(:), method, NaN);
end
end

function [P] = pDiagMatrix(p)
% per-epoch P diagonal as a K x n matrix, from a matrix, a cell of P / diag(P), or one P / diag(P)
if ~iscell(p)
    if size(p, 1) == size(p, 2) && size(p, 1) > 1, P = diag(p)';       % one full P at tOut
    elseif isvector(p),                           P = p(:)';           % one diagonal at tOut
    else,                                         P = p;               % K x n or N x n
    end
    return
end
if all(cellfun(@(x) isnumeric(x) && isscalar(x), p(:)))                   % num2cell of a K x n matrix
    P = cell2mat(p);
    if isvector(P), P = P(:)'; end
    return
end
K = numel(p);
P = NaN(K, numel(diag(p{1})));
for k = 1:K
    v = p{k};
    if isvector(v), P(k, :) = v(:)'; else, P(k, :) = diag(v)'; end
end
end

function [C] = dcmFromEuler(e)
% body -> nav, ZYX: C = Rz(yaw) Ry(pitch) Rx(roll)
cr = cos(e(1)); sr = sin(e(1)); cp = cos(e(2)); sp = sin(e(2)); cy = cos(e(3)); sy = sin(e(3));
C = [cy*cp, cy*sp*sr - sy*cr, cy*sp*cr + sy*sr;
     sy*cp, sy*sp*sr + cy*cr, sy*sp*cr - cy*sr;
     -sp,   cp*sr,            cp*cr];
end

function [s] = sigStr(v)
if all(isfinite(v))
    s = sprintf('%.3g ', v);
else
    s = 'n/a';
end
end
