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
%                       tKf instead (held between updates).
%     tKf      [K x 1]  2 Hz update epochs (s), needed when biasEst or
%                       pDiag are logged per update rather than per nav sample
%     tTrue    [M x 1]  truth time (s), may be a denser or offset time base
%     truePos  [M x 3]  truth [lat lon h],   trueVel [M x 3], trueAtt [M x 3]
%     biasTrue [1 x 6]  simulated IMU biases (or [M x 6] if they vary)
%     pDiag    optional sigma source: [N x 22] diagonal of P per nav sample,
%                       or [K x 22] per 2 Hz epoch, or a cell {K x 1} holding
%                       per epoch either the 22 x 22 P or its diagonal
%     name     char, optional
%   tOut     last GNSS-aided epoch before the outage (s)
%   tReacq   first GNSS epoch after the outage (s)
%   aidedWin [K x 2]  GNSS-aided windows [t1 t2] used for the truth-lag
%            estimate and the aided-accuracy figures, e.g. [100 760; 900 1600].
%            Optional: omitted or [], the 300 s before tOut are used.
%
%   What it prints, per run
%     0. truth lag: the shift of the truth time base that minimises the
%        aided horizontal position error (the call-time sampling artefact)
%     1. aided accuracy: RMS position and velocity error, bias error
%     2. state error at tOut: velocity, tilt, bias, with sigma from pDiag
%     3. outage trace fit: err = a + b x + c x^2 + d x^3 per NED axis, x from
%        tOut; the contributions b T, c T^2, d T^3 at T = tReacq - tOut say
%        whether the drift is velocity (t), bias/tilt (t^2) or gyro (t^3)
%     4. state budget: the same drift predicted from the tOut errors
%   and, with two runs, the difference of the tOut state errors in sigmas.
%
%   Conventions used in section 4 (if a predicted term has the right size
%   and the wrong sign, the host uses the opposite convention for that
%   quantity; the magnitudes stand):
%     tilt psi from C_nav * C_true' = I - [psi x], i.e. for small errors
%     psi = -C_b^n * (navEuler - trueEuler); specific-force error
%     (g psi_E, -g psi_N, 0), so a nose-up pitch error drifts south;
%     accel-bias error enters as -(biasEst - biasTrue) on the corrected IMU.
%   MATLAB desktop tool, not codegen. Runs under Octave.

G = 9.80665;
if nargin < 5 || isempty(aidedWin), aidedWin = [tOut - 300, tOut]; end
runs = {runA};
if ~isempty(runB), runs{end+1} = runB; end
res = struct();
for k = 1:numel(runs)
    r = runs{k};
    if ~isfield(r, 'name') || isempty(r.name), r.name = sprintf('run %d', k); end
    fprintf('\n==================== %s ====================\n', r.name);
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

%% quantities logged per 2 Hz update: bias held between updates, P taken at the last update <= tOut
P = [];
if isfield(r, 'pDiag') && ~isempty(r.pDiag), P = pDiagMatrix(r.pDiag); end
if isfield(r, 'tKf') && ~isempty(r.tKf)
    tKf = r.tKf(:);
    if size(r.biasEst, 1) == numel(tKf) && numel(tKf) ~= N
        r.biasEst = interpRows(tKf, r.biasEst, t, 'previous');
    end
    if ~isempty(P) && size(P, 1) == numel(tKf)
        iKf = find(tKf <= tOut, 1, 'last');
        P = repmat(P(iKf, :), N, 1);
    end
end
if size(r.biasEst, 1) ~= N
    error('runCompare:bias', 'biasEst has %d rows; give it per nav sample (%d) or per 2 Hz epoch with tKf', size(r.biasEst, 1), N);
end
inAided = false(N, 1);
for w = 1:size(aidedWin, 1)
    inAided = inAided | (t >= aidedWin(w, 1) & t <= aidedWin(w, 2));
end

%% 0. truth lag: nav(t) is compared with truth(t + tau)
tauGrid = -1.5:0.01:1.5;
rmsH = NaN(size(tauGrid));
for i = 1:numel(tauGrid)
    e = posErrNed(r, t(inAided) + tauGrid(i), inAided);
    ok = all(isfinite(e), 2);
    rmsH(i) = sqrt(mean(sum(e(ok, 1:2).^2, 2)));
end
[rmsBest, iBest] = min(rmsH);
tau = tauGrid(iBest);
rms0 = rmsH(tauGrid == 0);
fprintf('0. truth lag: nav(t) matches truth(t %+.2f s); aided horizontal rms %.2f m (was %.2f m at zero lag)\n', ...
    tau, rmsBest, rms0);
if abs(tau) >= 1.49
    fprintf('   WARNING: lag at the edge of the scan, extend tauGrid or check the time bases\n');
end

%% error traces with the lag applied
tq      = t + tau;
posErr  = posErrNed(r, tq, true(N, 1));
velErr  = r.navVel - interpRows(r.tTrue, r.trueVel, tq);
attTrue = interpRows(r.tTrue, r.trueAtt, tq);
attErr  = r.navAtt - attTrue;
attErr(:, 3) = atan2(sin(attErr(:, 3)), cos(attErr(:, 3)));
if size(r.biasTrue, 1) > 1
    biasTrue = interpRows(r.tTrue, r.biasTrue, tq);
else
    biasTrue = repmat(r.biasTrue(:)', N, 1);
end
biasErr = r.biasEst - biasTrue;

%% 1. aided accuracy
ok = inAided & all(isfinite(posErr), 2) & all(isfinite(velErr), 2);
rmsPos = sqrt(mean(posErr(ok, :).^2, 1));
rmsVel = sqrt(mean(velErr(ok, :).^2, 1));
meanPos = mean(posErr(ok, :), 1);
fprintf('1. aided windows: pos rms N %.2f E %.2f D %.2f m (mean N %.2f E %.2f D %.2f), vel rms N %.3f E %.3f D %.3f m/s\n', ...
    rmsPos, meanPos, rmsVel);
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
fprintf('3. outage drift observed: N %+8.2f E %+8.2f D %+8.2f m\n', drift);
fprintf('   fit err = a + b x + c x^2 + d x^3 over %d samples, contributions at T = %.0f s:\n', nnz(sel), T);
fprintf('   axis      b T (vel)    c T^2 (bias/tilt)    d T^3 (gyro)    fit rms   b [m/s]\n');
coef = zeros(3, 4);
for ax = 1:3
    p = polyfit(x, posErr(sel, ax), 3);              % [d c b a]
    coef(ax, :) = p;
    fitRms = sqrt(mean((polyval(p, x) - posErr(sel, ax)).^2));
    fprintf('   %s     %+9.2f     %+9.2f            %+9.2f       %6.2f    %+.3f\n', ...
        'NED'(ax), p(3) * T, p(2) * T^2, p(1) * T^3, fitRms, p(3));
end
fprintf('   (b should agree with the vel err of section 2; a large gap means the lag or the logged solution is wrong)\n');

%% 4. budget from the tOut state
velTerm  = dv * T;
biasTerm = 0.5 * (Cnav * (-db(1:3))) * T^2;                 % -(est - true) is the residual on the corrected IMU
dfTilt   = G * [psi(2); -psi(1); 0];
tiltTerm = 0.5 * dfTilt * T^2;
psiDot   = Cnav * (-db(4:6));
gyroTerm = G * [psiDot(2); -psiDot(1); 0] * T^3 / 6;
total    = velTerm + biasTerm + tiltTerm + gyroTerm;
fprintf('4. budget from the tOut state (m):   N         E         D\n');
fprintf('   velocity  dv T           %+8.2f  %+8.2f  %+8.2f\n', velTerm);
fprintf('   accel bias 1/2 C db T^2  %+8.2f  %+8.2f  %+8.2f\n', biasTerm);
fprintf('   tilt      1/2 g psi T^2  %+8.2f  %+8.2f  %+8.2f\n', tiltTerm);
fprintf('   gyro bias g psidot T^3/6 %+8.2f  %+8.2f  %+8.2f\n', gyroTerm);
fprintf('   sum                      %+8.2f  %+8.2f  %+8.2f\n', total);
fprintf('   observed                 %+8.2f  %+8.2f  %+8.2f\n', drift);
bt = r.biasTrue(1, 1:3);
fprintf('   scale: true accel bias uncompensated over T = %.0f m, applied with the wrong sign = %.0f m\n', ...
    0.5 * norm(bt) * T^2, norm(bt) * T^2);

o = struct('name', r.name, 'tau', tau, 'rmsPos', rmsPos, 'rmsVel', rmsVel, 'tOut', t(iOut), 'T', T, ...
    'posErrOut', posErr(iOut, :), 'dv', dv, 'attErr', attErr(iOut, :)', 'psi', psi, 'db', db, 'sig', sig, ...
    'drift', drift, 'coef', coef, 'velTerm', velTerm, 'biasTerm', biasTerm, 'tiltTerm', tiltTerm, ...
    'gyroTerm', gyroTerm, 'posErrTrace', posErr, 'velErrTrace', velErr, 'attErrTrace', attErr, 'biasErrTrace', biasErr, 't', t);
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
fprintf('   aided pos rms    %s: N %5.2f E %5.2f | %s: N %5.2f E %5.2f m\n', nameA, a.rmsPos(1:2), nameB, b.rmsPos(1:2));
fprintf('   aided vel rms    %s: N %5.3f E %5.3f | %s: N %5.3f E %5.3f m/s\n', nameA, a.rmsVel(1:2), nameB, b.rmsVel(1:2));
fprintf('   A difference of a few sigmas or less on every row is two draws of the same filter; many sigmas on the\n');
fprintf('   bias or tilt rows means the two filters disagree beyond their own uncertainty, i.e. at least one is inconsistent.\n');
end

%------------------------------------------------------------------------
function [e] = posErrNed(r, tq, mask)
% nav(mask) minus truth interpolated at tq, in NED metres
navPos = r.navPos(mask, :);
tp = interpRows(r.tTrue, r.truePos, tq);
L  = navPos(:, 1);  h = navPos(:, 3);
a  = 6378137;  e2 = 6.69437999014e-3;
RN = a * (1 - e2) ./ (1 - e2 * sin(L).^2).^1.5;
RE = a ./ sqrt(1 - e2 * sin(L).^2);
e = [(navPos(:, 1) - tp(:, 1)) .* (RN + h), ...
     (navPos(:, 2) - tp(:, 2)) .* (RE + h) .* cos(L), ...
    -(navPos(:, 3) - tp(:, 3))];
end

function [y] = interpRows(tt, Y, tq, method)
if nargin < 4, method = 'linear'; end
y = NaN(numel(tq), size(Y, 2));
for c = 1:size(Y, 2)
    y(:, c) = interp1(tt(:), Y(:, c), tq(:), method, NaN);
end
end

function [P] = pDiagMatrix(p)
% per-epoch P diagonal as a K x n matrix, from a matrix or a cell of P / diag(P)
if ~iscell(p), P = p; return; end
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
