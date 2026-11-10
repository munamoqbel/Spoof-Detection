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
%               0.05, 0.0147, 0.03, 4, 300, 0.016, 10, gyro first.
%   runC        the same struct for the errors-off run: its coast error is the mechanisation's own
%
%   From tOut the host's P is propagated with the host's accumPhi three ways: with the host's accumQ (the P the host
%   carries), with the datasheet Q on the inertial states (velocity from the VRW, tilt from the ARW, the bias states
%   from the in-run stability and the walk, the host's accumQ on the other states) and with the datasheet Q from a
%   zero P (what the sensor alone adds). The table gives, at 5, 10, 30, 60, 90 and 120 s into the coast and at its
%   end, the measured nav error (N, E, D, horizontal), the logged horizontal sigma when pAtOut holds every epoch, the
%   three predicted horizontal sigmas and the errors-off run's horizontal error; then the ratios at the end of the
%   coast or at 120 s: the host's sigma over the measured error (a P three times the drift reads 3) and the
%   datasheet's (near 1 when the datasheet and the host's errors at tOut explain the coast; below 1 the coast
%   carries something the datasheet does not, and runC says whether it is the mechanisation).
%   Result: t (the coast epochs), T (seconds from tOut), err [m x 3], errH, sigLogged, sigHost, sigDs, sigQ, errC,
%   ratioHost, ratioDs, tEnd, note. MATLAB desktop tool, not codegen. Runs under Octave.
G = 9.80665;  R2D = 180 / pi;
if nargin < 7 || isempty(imu), imu = struct(); end
def = struct('arw', 0.05, 'vrw', 0.0147, 'gyroBias', 0.03, 'accBias', 4, 'tauBias', 300, 'gyroWalk', 0.016, 'accWalk', 10, 'gyroFirst', true);
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

%% the three predictions along the host's accumPhi
PH = P0;  PD = P0;  PQ = zeros(n);
m = kEnd - kOut;
sigHost = NaN(m, 1);  sigDs = NaN(m, 1);  sigQ = NaN(m, 1);  sigHostD = NaN(m, 1);  sigDsD = NaN(m, 1);
for j = 1:m
    k = kOut + j;
    Phi = L(k).accumPhi;
    if isfield(L, 'accumQ') && ~isempty(L(k).accumQ), Qh = L(k).accumQ; else, Qh = zeros(n); end
    dt = Phi(1, 4);                                                    % the interval the accumulation covers
    if ~(dt > 0) || ~isfinite(dt), dt = tKf(k) - tKf(k - 1); end
    Qds = dt / 2 * (Phi * Qc * Phi' + Qc);                             % the datasheet Q over the interval
    Qo = Qh;  Qo(inertial, :) = 0;  Qo(:, inertial) = 0;               % the host's Q on the states the IMU does not drive
    Qds = Qds + Qo;
    PH = sym(Phi * PH * Phi' + Qh);
    PD = sym(Phi * PD * Phi' + Qds);
    PQ = sym(Phi * PQ * Phi' + Qds);
    sigHost(j) = sqrt(PH(1, 1) + PH(2, 2));  sigHostD(j) = sqrt(PH(3, 3));
    sigDs(j) = sqrt(PD(1, 1) + PD(2, 2));    sigDsD(j) = sqrt(PD(3, 3));
    sigQ(j) = sqrt(PQ(1, 1) + PQ(2, 2));
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
fprintf('   %7s %9s %9s %9s %9s | %9s %9s %9s %9s | %9s\n', 'T [s]', 'err N', 'err E', 'err D', 'err H', 'logged', 'host P', 'datasheet', 'Q only', 'run C H');
fprintf('   %7s %9s %9s %9s %9s | %9s %9s %9s %9s | %9s\n', '', '[m]', '[m]', '[m]', '[m]', 'sig H', 'sig H', 'sig H', 'sig H', '[m]');
pick = [];
for x = [5 10 30 60 90 120]
    j = find(T >= x, 1);
    if ~isempty(j), pick(end + 1) = j; end %#ok<AGROW>
end
pick = unique([pick, m]);
for j = pick
    fprintf('   %7.1f %9.2f %9.2f %9.2f %9.2f | %9s %9.2f %9.2f %9.2f | %9s\n', T(j), err(j, 1), err(j, 2), err(j, 3), errH(j), ...
        num2(sigLogged(j)), sigHost(j), sigDs(j), sigQ(j), num2(errC(j)));
end
jv = find(T >= 120, 1);  if isempty(jv), jv = m; end
ratioHost = sigHost(jv) / errH(jv);  ratioDs = sigDs(jv) / errH(jv);
fprintf('   at %.0f s: the host''s sigma is %.2f times the measured horizontal error, the datasheet''s %.2f times; the sensor alone gives %.2f m\n', ...
    T(jv), ratioHost, ratioDs, sigQ(jv));
if isfinite(errC(jv)), fprintf('   the errors-off run drifted %.2f m by then: the mechanisation''s own share of the coast\n', errC(jv)); end
if ratioDs >= 0.7 && ratioDs <= 2
    fprintf('   VERDICT: the datasheet Q along the host''s Phi explains the coast; the host''s P overstates it %.1f times\n', ratioHost);
elseif ratioDs < 0.7
    fprintf('   VERDICT: the coast carries %.1f times what the datasheet and the tOut state predict: an error the datasheet does not describe\n', 1 / ratioDs);
    fprintf('           (the mechanisation if run C shows it, the host''s state at tOut otherwise)\n');
else
    fprintf('   VERDICT: even the datasheet Q overstates this coast %.1f times: the host''s P at tOut carries more than the real errors there\n', ratioDs);
end
res = struct('t', t, 'T', T, 'err', err, 'errH', errH, 'sigLogged', sigLogged, 'sigHost', sigHost, 'sigDs', sigDs, 'sigQ', sigQ, ...
    'sigHostD', sigHostD, 'sigDsD', sigDsD, 'errC', errC, 'ratioHost', ratioHost, 'ratioDs', ratioDs, 'tEnd', tKf(kEnd), 'note', note, 'Qc', Qc);
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
