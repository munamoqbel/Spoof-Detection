function [res] = coastDriftCheck(run, tKf, tOut, tReacq, propTelLog, pAtOut, imu, runC)
%COASTDRIFTCHECK  An INS-only coast against the host's P and the datasheet Q, both along the host's own Phi.
%   res = coastDriftCheck(run, tKf, tOut, tReacq, propTelLog, pAtOut)
%   res = coastDriftCheck(run, tKf, tOut, tReacq, propTelLog, pAtOut, imu)
%   res = coastDriftCheck(run, tKf, tOut, tReacq, propTelLog, pAtOut, imu, runC)
%   run         the run struct of runCompare for the coast run: t, navPos [lat lon h] (rad, rad, m), tTrue, truePos
%   tKf         [K x 1] the 2 Hz epochs (s)
%   tOut        the last aided epoch (s); tReacq the first epoch back, or the end of the coast
%   propTelLog  accumPhi and accumQ per epoch, in any form covPropCheck takes: a struct array or cell of per-epoch
%               structs, the host's per-epoch log with a propTel field, or one struct with the epochs along the third
%               dimension of accumPhi and accumQ
%   pAtOut      the host's P at tOut, n x n; or the P of every epoch as a cell, [n x n x K] or the diagonal [K x n],
%               from which the last epoch <= tOut is taken (a diagonal starts the prediction without the cross terms
%               between velocity, tilt and bias, which the printout notes) and whose later epochs give the logged sigma
%   imu         the datasheet in its measured form: arw (deg/sqrt(h)), vrw (m/s/sqrt(h)), gyroBias (deg/h, the in-run
%               sigma), accBias (micro-g), tauBias (s), gyroWalk (deg/h per sqrt(h)), accWalk (micro-g per sqrt(h)),
%               gyroFirst (true: gyro biases on states 10:12, accelerometer on 13:15). Default, or []: the vehicle IMU,
%               0.05, 0.0147, 0.03, 4, 300, 0.016, 10, gyro first. Optional: floorStates and floorSigmas, the gate's
%               COAST_SIGMA_FLOOR lists, which add the column the gate itself carries (the host's P at tOut floored,
%               the datasheet Q); headingState, the heading's index (default 9).
%   runC        the same struct for the errors-off run: its coast error is the mechanisation's own
%
%   From tOut the host's P is propagated with the host's accumPhi three ways: with the host's accumQ (the P the host
%   carries), with the datasheet Q on the inertial states (velocity from the VRW, tilt from the ARW, the bias states
%   from the in-run stability and the walk, the host's accumQ on the other states) and with the datasheet Q from a
%   zero P (what the sensor alone adds), and with the gate's floors on the P at tOut when they are given. The table
%   gives, at 5, 10, 30, 60, 90 and 120 s into the coast and at its end, the measured nav error (N, E, D,
%   horizontal), the logged horizontal sigma when pAtOut holds every epoch, the predicted horizontal sigmas and the
%   errors-off run's horizontal error; then the ratios at 120 s or at the end: the host's sigma over the measured
%   error (a P three times the drift reads 3), the datasheet's (near 1 when the datasheet and the host's errors at
%   tOut explain the coast; below 1 the coast carries something the datasheet does not, and runC says whether it is
%   the mechanisation) and the gate's. The heading is the state a wheeled vehicle's host cannot observe on straight
%   legs and its P understates: the tool also solves for the heading sigma at tOut that would explain the coast at
%   120 s with the datasheet Q along this Phi (the heading reaches the position only where the host's Phi couples
%   it, in the manoeuvres of the coast), which is the heading floor the gate needs, read from one coast.
%   Result: t (the coast epochs), T (seconds from tOut), err [m x 3], errH, sigLogged, sigHost, sigDs, sigQ, sigGate
%   (NaN without floors), errC, ratioHost, ratioDs, ratioGate, headingSigmaHost, headingFloor, headingImplied (rad;
%   0 when the datasheet already covers the coast, Inf when 0.5 rad does not), tEnd, note.
%   MATLAB desktop tool, not codegen. Runs under Octave.
G = 9.80665;  R2D = 180 / pi;
if nargin < 7 || isempty(imu), imu = struct(); end
def = struct('arw', 0.05, 'vrw', 0.0147, 'gyroBias', 0.03, 'accBias', 4, 'tauBias', 300, 'gyroWalk', 0.016, 'accWalk', 10, 'gyroFirst', true, ...
    'floorStates', zeros(1, 0), 'floorSigmas', zeros(1, 0), 'headingState', 9);
f = fieldnames(def);
for i = 1:numel(f), if ~isfield(imu, f{i}) || isempty(imu.(f{i})), imu.(f{i}) = def.(f{i}); end, end
if nargin < 8, runC = []; end
tKf = tKf(:);  K = numel(tKf);
L = propTelRows(propTelLog);
if isempty(L) || ~isfield(L, 'accumPhi'), error('coastDriftCheck: propTelLog carries no accumPhi'); end
if numel(L) ~= K
    fprintf('coastDriftCheck: propTelLog has %d entries, tKf %d epochs: the first %d are used\n', numel(L), K, min(numel(L), K));
    K = min(numel(L), K);  tKf = tKf(1:K);
end
n = size(L(1).accumPhi, 1);
kOut = find(tKf <= tOut, 1, 'last');
if isempty(tReacq), kEnd = K; else, kEnd = find(tKf <= tReacq, 1, 'last'); end
if isempty(kOut) || kEnd <= kOut, error('coastDriftCheck: no coast epochs between tOut and tReacq'); end
[P0, sigLoggedAll, note] = pAtEpochs(pAtOut, kOut, K, n);

%% the datasheet Q on the inertial states, continuous densities in SI
if imu.gyroFirst, gyroIdx = 10:12;  accIdx = 13:15; else, accIdx = 10:12;  gyroIdx = 13:15; end
velIdx = 4:6;  attIdx = 7:9;  inertial = 1:15;
qVel  = (imu.vrw / 60) ^ 2;                                            % (m/s)^2 per s
qAtt  = (imu.arw / R2D / 60) ^ 2;                                      % rad^2 per s
sG = imu.gyroBias / R2D / 3600;  kG = imu.gyroWalk / R2D / 3600 / 60;  % rad/s; rad/s per sqrt(s)
sA = imu.accBias * G * 1e-6;     kA = imu.accWalk * G * 1e-6 / 60;     % m/s^2; m/s^2 per sqrt(s)
qGyro = 2 * sG ^ 2 / imu.tauBias + kG ^ 2;                             % the Markov bias as a walk over a coast, plus the walk
qAcc  = 2 * sA ^ 2 / imu.tauBias + kA ^ 2;
Qc = zeros(n);
Qc(velIdx, velIdx) = qVel * eye(3);  Qc(attIdx, attIdx) = qAtt * eye(3);
Qc(gyroIdx, gyroIdx) = qGyro * eye(3);  Qc(accIdx, accIdx) = qAcc * eye(3);

%% the predictions along the host's accumPhi
m = kEnd - kOut;
prop = struct('L', L, 'kOut', kOut, 'kEnd', kEnd, 'tKf', tKf, 'Qc', Qc, 'inertial', inertial, 'n', n);
[sigHost, sigHostD] = propagateSigma(prop, P0, true);
[sigDs, sigDsD] = propagateSigma(prop, P0, false);
sigQ = propagateSigma(prop, zeros(n), false);
hasFloor = ~isempty(imu.floorStates);
P0f = P0;
if hasFloor
    for k = 1:numel(imu.floorStates)
        idx = imu.floorStates(k);
        if P0f(idx, idx) < imu.floorSigmas(k)^2, P0f(idx, idx) = imu.floorSigmas(k)^2; end
    end
    sigGate = propagateSigma(prop, P0f, false);
else
    sigGate = NaN(m, 1);
end
t = tKf(kOut + 1:kEnd);  T = t - tKf(kOut);
sigLogged = sigLoggedAll(kOut + 1:kEnd);

%% the measured coast
err = navErrAt(run, t);  errH = sqrt(sum(err(:, 1:2) .^ 2, 2));
e0 = navErrAt(run, tKf(kOut));
errC = NaN(m, 1);
if ~isempty(runC), eC = navErrAt(runC, t);  errC = sqrt(sum(eC(:, 1:2) .^ 2, 2)); end

%% the table
fprintf('coastDriftCheck: coast from %.1f s to %.1f s (%d epochs, %.1f s); the nav error at tOut N %+.2f E %+.2f D %+.2f m\n', ...
    tKf(kOut), tKf(kEnd), m, T(end), e0(1), e0(2), e0(3));
fprintf('   datasheet: ARW %g deg/sqrt(h), VRW %g m/s/sqrt(h), bias sigma %g deg/h and %g micro-g with tau %g s, walks %g deg/h/sqrt(h) and %g micro-g/sqrt(h)\n', ...
    imu.arw, imu.vrw, imu.gyroBias, imu.accBias, imu.tauBias, imu.gyroWalk, imu.accWalk);
if ~isempty(note), fprintf('   NOTE: %s\n', note); end
if hasFloor
    fprintf('   the gate''s floors on the P at tOut: states %s at sigma %s\n', mat2str(imu.floorStates), mat2str(imu.floorSigmas, 3));
end
fprintf('   %7s %9s %9s %9s %9s | %9s %9s %9s %9s %9s | %9s\n', 'T [s]', 'err N', 'err E', 'err D', 'err H', 'logged', 'host P', 'datasheet', 'Q only', 'gate', 'run C H');
fprintf('   %7s %9s %9s %9s %9s | %9s %9s %9s %9s %9s | %9s\n', '', '[m]', '[m]', '[m]', '[m]', 'sig H', 'sig H', 'sig H', 'sig H', 'sig H', '[m]');
pick = [];
for x = [5 10 30 60 90 120]
    j = find(T >= x, 1);
    if ~isempty(j), pick(end + 1) = j; end %#ok<AGROW>
end
pick = unique([pick, m]);
for j = pick
    fprintf('   %7.1f %9.2f %9.2f %9.2f %9.2f | %9s %9.2f %9.2f %9.2f %9s | %9s\n', T(j), err(j, 1), err(j, 2), err(j, 3), errH(j), ...
        num2(sigLogged(j)), sigHost(j), sigDs(j), sigQ(j), num2(sigGate(j)), num2(errC(j)));
end
jv = find(T >= 120, 1);  if isempty(jv), jv = m; end
ratioHost = sigHost(jv) / errH(jv);  ratioDs = sigDs(jv) / errH(jv);  ratioGate = sigGate(jv) / errH(jv);
fprintf('   at %.0f s: the host''s sigma is %.2f times the measured horizontal error, the datasheet''s %.2f times', T(jv), ratioHost, ratioDs);
if hasFloor, fprintf(', the gate''s %.2f times', ratioGate); end
fprintf('; the sensor alone gives %.2f m\n', sigQ(jv));
if isfinite(errC(jv)), fprintf('   the errors-off run drifted %.2f m by then: the mechanisation''s own share of the coast\n', errC(jv)); end
% the heading sigma at tOut that would explain the coast: the floor the gate needs, from this one coast
hs = imu.headingState;
headingSigmaHost = sqrt(max(P0(hs, hs), 0));
headingFloor = 0;
if hasFloor
    iH = find(imu.floorStates == hs, 1);
    if ~isempty(iH), headingFloor = imu.floorSigmas(iH); end
end
headingImplied = impliedHeading(prop, P0f, hs, jv, errH(jv));
fprintf('   heading at tOut: the host''s P carries %.1f mrad', headingSigmaHost * 1e3);
if hasFloor, fprintf(', the floor in force %.1f mrad', headingFloor * 1e3); end
if headingImplied == 0 && hasFloor
    fprintf('; with the floors in force the datasheet Q covers the coast at %.0f s: the floors are enough for this coast\n', T(jv));
elseif headingImplied == 0
    fprintf('; the datasheet Q already covers the coast at %.0f s: no heading floor needed from this coast\n', T(jv));
elseif isinf(headingImplied)
    fprintf('; no heading sigma up to 500 mrad explains the coast at %.0f s along this Phi: the gap is not the heading\n', T(jv));
    fprintf('           (the velocity, the tilt or the biases at tOut, or the mechanisation; run C and runCompare section 2 say which)\n');
else
    fprintf('; the heading sigma that explains the coast at %.0f s along this Phi is %.1f mrad\n', T(jv), headingImplied * 1e3);
    fprintf('           (the heading floor this coast asks for, read from one coast: a sample with its scatter, not an rms)\n');
end
if ratioDs >= 0.7 && ratioDs <= 2
    fprintf('   VERDICT: the datasheet Q along the host''s Phi explains the coast; the host''s P overstates it %.1f times\n', ratioHost);
elseif ratioDs < 0.7
    if hasFloor && ratioGate >= 0.7 && ratioGate <= 2
        fprintf('   VERDICT: the coast carries %.1f times what the datasheet and the host''s tOut state predict; with the gate''s floors it reads\n', 1 / ratioDs);
        fprintf('           %.2f of the measured error: the floors carry the state the host''s P understates\n', ratioGate);
    else
        fprintf('   VERDICT: the coast carries %.1f times what the datasheet and the tOut state predict: an error the datasheet does not describe\n', 1 / ratioDs);
        fprintf('           (the mechanisation if run C shows it, the host''s state at tOut otherwise; the heading line above says whether a floor covers it)\n');
    end
else
    fprintf('   VERDICT: even the datasheet Q overstates this coast %.1f times: the host''s P at tOut carries more than the real errors there\n', ratioDs);
end
res = struct('t', t, 'T', T, 'err', err, 'errH', errH, 'sigLogged', sigLogged, 'sigHost', sigHost, 'sigDs', sigDs, 'sigQ', sigQ, ...
    'sigGate', sigGate, 'sigHostD', sigHostD, 'sigDsD', sigDsD, 'errC', errC, 'ratioHost', ratioHost, 'ratioDs', ratioDs, 'ratioGate', ratioGate, ...
    'headingSigmaHost', headingSigmaHost, 'headingFloor', headingFloor, 'headingImplied', headingImplied, 'tEnd', tKf(kEnd), 'note', note, 'Qc', Qc);
end

function [sigH, sigD] = propagateSigma(prop, P0, useHostQ)
% the horizontal and vertical position sigma along the coast from P0 with the host's accumPhi, and the host's accumQ or
% the datasheet Q on the inertial states (the host's Q on the others)
m = prop.kEnd - prop.kOut;  n = prop.n;
sigH = NaN(m, 1);  sigD = NaN(m, 1);
P = (P0 + P0') / 2;
for j = 1:m
    k = prop.kOut + j;
    Phi = prop.L(k).accumPhi;
    if isfield(prop.L, 'accumQ') && ~isempty(prop.L(k).accumQ), Qh = prop.L(k).accumQ; else, Qh = zeros(n); end
    if useHostQ
        Q = Qh;
    else
        dt = Phi(1, 4);                                                % the interval the accumulation covers
        if ~(dt > 0) || ~isfinite(dt), dt = prop.tKf(k) - prop.tKf(k - 1); end
        Q = dt / 2 * (Phi * prop.Qc * Phi' + prop.Qc);
        Qo = Qh;  Qo(prop.inertial, :) = 0;  Qo(:, prop.inertial) = 0;
        Q = Q + Qo;
    end
    P = Phi * P * Phi' + Q;  P = (P + P') / 2;
    sigH(j) = sqrt(P(1, 1) + P(2, 2));  sigD(j) = sqrt(P(3, 3));
end
end

function [s] = impliedHeading(prop, P0, hs, jv, target)
% the heading sigma at tOut at which the datasheet prediction meets the measured horizontal error at epoch jv:
% 0 when the prediction already covers it, Inf when 0.5 rad does not, else by bisection on the (hs, hs) entry
sig0 = propagateSigma(prop, P0, false);
if sig0(jv) >= target, s = 0; return; end
sMax = 0.5;
if predictAt(prop, P0, hs, sMax, jv) < target, s = Inf; return; end
lo = sqrt(max(P0(hs, hs), 0));  hi = sMax;
for it = 1:40
    mid = 0.5 * (lo + hi);
    if predictAt(prop, P0, hs, mid, jv) < target, lo = mid; else, hi = mid; end
end
s = 0.5 * (lo + hi);
end

function [v] = predictAt(prop, P0, hs, s, jv)
P = P0;  P(hs, hs) = max(P(hs, hs), s^2);
sig = propagateSigma(prop, P, false);
v = sig(jv);
end

function [P0, sigLogged, note] = pAtEpochs(pAtOut, kOut, K, n)
% the host's P at tOut and the logged horizontal position sigma per epoch, from the forms pAtOut may take
note = '';  sigLogged = NaN(K, 1);
if iscell(pAtOut)
    P0 = full22(pAtOut{min(kOut, numel(pAtOut))}, n);
    if numel(pAtOut{1}) ~= n * n, note = 'pAtOut holds diagonals: the prediction starts without the cross terms between velocity, tilt and bias'; end
    for k = 1:min(K, numel(pAtOut)), Pk = full22(pAtOut{k}, n);  sigLogged(k) = sqrt(Pk(1, 1) + Pk(2, 2)); end
elseif ndims(pAtOut) == 3
    P0 = pAtOut(:, :, min(kOut, size(pAtOut, 3)));
    for k = 1:min(K, size(pAtOut, 3)), sigLogged(k) = sqrt(pAtOut(1, 1, k) + pAtOut(2, 2, k)); end
elseif isequal(size(pAtOut), [n n])
    P0 = pAtOut;
elseif size(pAtOut, 2) == n && size(pAtOut, 1) >= kOut
    P0 = diag(pAtOut(kOut, :));
    note = 'pAtOut is a diagonal per epoch: the prediction starts without the cross terms between velocity, tilt and bias';
    for k = 1:min(K, size(pAtOut, 1)), sigLogged(k) = sqrt(pAtOut(k, 1) + pAtOut(k, 2)); end
else
    error('coastDriftCheck: pAtOut is neither an n x n P, a cell or [n x n x K] of P, nor a [K x n] diagonal');
end
P0 = (P0 + P0') / 2;
end

function [P] = full22(x, n)
if numel(x) == n * n, P = reshape(x, n, n); else, P = diag(x(:)); end
end

function [e] = navErrAt(run, t)
% nav minus truth in NED (m) at the times t: the nav interpolated on its own time base, the truth on its own
a = 6378137;  e2 = 6.69437999014e-3;
nav = NaN(numel(t), 3);  tr = NaN(numel(t), 3);
for c = 1:3
    nav(:, c) = interp1(run.t(:), run.navPos(:, c), t(:), 'linear', NaN);
    tr(:, c) = interp1(run.tTrue(:), run.truePos(:, c), t(:), 'linear', NaN);
end
lat = tr(:, 1);  s2 = sin(lat) .^ 2;
Rm = a * (1 - e2) ./ (1 - e2 * s2) .^ 1.5;  Rn = a ./ sqrt(1 - e2 * s2);
e = [(nav(:, 1) - tr(:, 1)) .* (Rm + tr(:, 3)), (nav(:, 2) - tr(:, 2)) .* (Rn + tr(:, 3)) .* cos(lat), -(nav(:, 3) - tr(:, 3))];
end

function [L] = propTelRows(propTelLog)
% propTel per epoch as a struct array, from a struct array or cell of per-epoch structs, from the host's per-epoch
% log with a propTel field, or from one struct whose accumPhi / accumQ hold the epochs along the third dimension or
% in a cell (as covPropCheck reads them)
if iscell(propTelLog), propTelLog = [propTelLog{:}]; end
if isstruct(propTelLog) && isfield(propTelLog, 'propTel') && ~isfield(propTelLog, 'accumPhi')
    propTelLog = [propTelLog.propTel];
end
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
        A = permute(A, [2 3 1]);
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

function [P] = sym(P)
P = (P + P') / 2;
end

function [s] = num2(x)
if isnan(x), s = '-'; else, s = sprintf('%.2f', x); end
end
