function test_coastDriftCheck()
%% test_coastDriftCheck.m
% tools/coastDriftCheck on a synthetic 22-state INS error model: Phi from the error dynamics over 0.5 s (position from
% velocity, velocity from tilt through gravity and from the accelerometer bias, tilt from the gyro bias), a host Q nine
% times the datasheet's, and a Monte Carlo of coasts driven by the datasheet Q from a known P at tOut. The datasheet
% sigma the tool propagates must match the ensemble rms of the coasts at 120 s within 10 %, the host sigma must sit
% well above it, the measured column must reproduce the one realisation given as the run, the errors-off run must
% read zero, and a diagonal P must be accepted with its note. Run from anywhere:
% octave --no-gui --quiet --eval test_coastDriftCheck
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'tools'));
rng(4);
n = 22;  dt = 0.5;  G = 9.80665;
F = zeros(n);
F(1:3, 4:6) = eye(3);                           % position from velocity
F(4, 8) = G;  F(5, 7) = -G;                     % velocity from tilt: (g psi_E, -g psi_N, 0)
F(4:6, 13:15) = eye(3);                         % velocity from the accelerometer bias (level vehicle)
F(7:9, 10:12) = eye(3);                         % tilt from the gyro bias
Phi = expm(F * dt);
Phi(1, 4) = dt;                                 % the interval the accumulation covers, as the host's accumPhi holds it
imu = struct('arw', 0.05, 'vrw', 0.0147, 'gyroBias', 0.03, 'accBias', 4, 'tauBias', 300, 'gyroWalk', 0.016, 'accWalk', 10, 'gyroFirst', true);
R2D = 180 / pi;
Qc = zeros(n);
Qc(4:6, 4:6) = (imu.vrw / 60)^2 * eye(3);  Qc(7:9, 7:9) = (imu.arw / R2D / 60)^2 * eye(3);
Qc(10:12, 10:12) = (2 * (imu.gyroBias / R2D / 3600)^2 / imu.tauBias + (imu.gyroWalk / R2D / 3600 / 60)^2) * eye(3);
Qc(13:15, 13:15) = (2 * (imu.accBias * G * 1e-6)^2 / imu.tauBias + (imu.accWalk * G * 1e-6 / 60)^2) * eye(3);
Qds = dt / 2 * (Phi * Qc * Phi' + Qc);  Qds = (Qds + Qds') / 2;
Qhost = 9 * Qds;                                % the host's Q: three times the datasheet in sigma
% the state at tOut: 1 m position, 0.05 m/s velocity, 0.1 mrad tilt, 0.05 deg/h and 30 micro-g bias
sig0 = zeros(n, 1);  sig0(1:3) = 1;  sig0(4:6) = 0.05;  sig0(7:9) = 1e-4;  sig0(10:12) = 0.05 / R2D / 3600;  sig0(13:15) = 30 * G * 1e-6;
sig0(16:22) = 1e-3;
P0 = diag(sig0 .^ 2);
% the coast: 600 s aided (not simulated), then 240 epochs (120 s) of coast, 400 realisations
K = 1200 + 240;  tKf = 0.5 * (1:K)';  tOut = 600;  kOut = 1200;
[V, D] = eig(Qds);  Sq = V * sqrt(max(D, 0));            % Qds^(1/2)
nMc = 400;  m = 240;
posEnd = zeros(nMc, 2);  one = zeros(m, 3);
for r = 1:nMc
    x = sig0 .* randn(n, 1);
    for j = 1:m
        x = Phi * x + Sq * randn(n, 1);
        if r == 1, one(j, :) = x(1:3)'; end
    end
    posEnd(r, :) = x(1:2)';
end
rmsH = sqrt(mean(sum(posEnd .^ 2, 2)));
% the run struct of the first realisation: a stationary truth, the nav offset by the simulated error
lat0 = 0.43;  lon0 = 0.81;  h0 = 43;
a = 6378137;  e2 = 6.69437999014e-3;
Rm = a * (1 - e2) / (1 - e2 * sin(lat0)^2)^1.5;  Rn = a / sqrt(1 - e2 * sin(lat0)^2);
t = tKf;  errN = zeros(K, 1);  errE = zeros(K, 1);  errD = zeros(K, 1);
errN(kOut + 1:end) = one(:, 1);  errE(kOut + 1:end) = one(:, 2);  errD(kOut + 1:end) = one(:, 3);
run = struct('t', t, 'navPos', [lat0 + errN / (Rm + h0), lon0 + errE / ((Rn + h0) * cos(lat0)), h0 - errD], ...
    'tTrue', [0; t(end) + 1], 'truePos', [lat0 lon0 h0; lat0 lon0 h0]);
runC = struct('t', t, 'navPos', repmat([lat0 lon0 h0], K, 1), 'tTrue', [0; t(end) + 1], 'truePos', [lat0 lon0 h0; lat0 lon0 h0]);
propTelLog = repmat(struct('accumPhi', Phi, 'accumQ', Qhost), K, 1);
fprintf('--- the full P at tOut ---\n');
r1 = coastDriftCheck(run, tKf, tOut, [], propTelLog, P0, imu, runC);
jEnd = numel(r1.T);
okDs = abs(r1.sigDs(jEnd) / rmsH - 1) < 0.1;
okHost = r1.sigHost(jEnd) > 1.4 * r1.sigDs(jEnd) && r1.ratioHost > r1.ratioDs;   % the tOut state is common to both; the host's Q is 3 times the datasheet's
okMeas = max(max(abs(r1.err - one))) < 1e-3 && all(r1.errC < 1e-6);
okQ = r1.sigQ(jEnd) < r1.sigDs(jEnd) && r1.sigQ(jEnd) > 0.3 * r1.sigDs(jEnd) && abs(r1.sigQ(jEnd) - 7.1) < 1.0;   % the ARW alone: 5 m per axis, 7.1 m horizontal at 120 s
fprintf('   datasheet sigma at the end %.2f m against the ensemble rms %.2f m; host sigma %.2f m; Q alone %.2f m\n', r1.sigDs(jEnd), rmsH, r1.sigHost(jEnd), r1.sigQ(jEnd));
fprintf('--- the diagonal of P per epoch ---\n');
pDiag = repmat(sig0' .^ 2, K, 1);
r2 = coastDriftCheck(run, tKf, tOut, tKf(end), propTelLog, pDiag, imu);
okDiag = ~isempty(r2.note) && abs(r2.sigDs(end) - r1.sigDs(jEnd)) < 0.2 * r1.sigDs(jEnd) && all(isfinite(r2.sigLogged)) && all(isnan(r2.errC));
fprintf('--- propTel as one struct with the epochs along the third dimension ---\n');
P3 = struct('accumPhi', repmat(Phi, [1 1 K]), 'accumQ', repmat(Qhost, [1 1 K]));
evalc('r3 = coastDriftCheck(run, tKf, tOut, [], P3, P0, imu);');
okForm = max(abs(r3.sigDs - r1.sigDs)) < 1e-9;
fprintf('datasheet sigma within 10 %% of the ensemble: %s | host sigma above it: %s | measured column and run C: %s | Q alone below the full prediction: %s | diagonal P with its note: %s | third-dimension propTel: %s\n', ...
    pf(okDs), pf(okHost), pf(okMeas), pf(okQ), pf(okDiag), pf(okForm));
if okDs && okHost && okMeas && okQ && okDiag && okForm, fprintf('=== test_coastDriftCheck PASS ===\n');
else, fprintf('=== test_coastDriftCheck FAIL ===\n'); end
end

function s = pf(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end
