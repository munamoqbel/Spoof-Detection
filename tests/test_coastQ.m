function test_coastQ()
%% test_coastQ.m
% The gate's own coast Q and the covariance floor. (1) SPF_coastQ on a Phi with the error dynamics (position from
% velocity, velocity from tilt and accelerometer bias, tilt from gyro bias, over 0.5 s) and a host Q nine times the
% datasheet: the inertial block must equal dt/2 (Phi Qc Phi' + Qc) built here from the same datasheet, the position
% block carry no noise of its own, the clock and the other states keep the host's entries, the cross terms between
% inertial and other states be zero, the result symmetric; with a Phi without the position-from-velocity entry the
% interval must fall back to COAST_Q_INTERVAL. (2) SPF_floorCov raises the listed diagonals that sit below their floor
% and leaves the rest, keeping the matrix positive semi-definite. (3) In the repository configuration (the flag off,
% no floors) SPF_coastPropTel and SPF_coastStartCov pass their inputs unchanged, so the test references stand.
% Run from the repo root (MATLAB, or Octave with tools/octave_shim).
n = CST_gnssHybrid.NO_STATES;  dt = 0.5;  G = 9.80665;  R2D = 180 / pi;
F = zeros(n);
F(1:3, 4:6) = eye(3);  F(4, 8) = G;  F(5, 7) = -G;
F(4:6, CST_spfParam.COAST_ACC_BIAS_STATES) = eye(3);  F(7:9, CST_spfParam.COAST_GYRO_BIAS_STATES) = eye(3);
F(16, 17) = 1;                                                         % clock bias from clock drift
Phi = expm(F * dt);  Phi(1, 4) = dt;
% the datasheet densities as the gate derives them
qVel = (CST_spfParam.IMU_VRW / 60)^2;  qAtt = (CST_spfParam.IMU_ARW / R2D / 60)^2;
sG = CST_spfParam.IMU_GYRO_BIAS_SIGMA / R2D / 3600;  sA = CST_spfParam.IMU_ACC_BIAS_SIGMA * G * 1e-6;
kG = CST_spfParam.IMU_GYRO_WALK / R2D / 3600 / 60;    kA = CST_spfParam.IMU_ACC_WALK * G * 1e-6 / 60;
Qc = zeros(n);
Qc(4:6, 4:6) = qVel * eye(3);  Qc(7:9, 7:9) = qAtt * eye(3);
Qc(CST_spfParam.COAST_GYRO_BIAS_STATES, CST_spfParam.COAST_GYRO_BIAS_STATES) = (2 * sG^2 / CST_spfParam.IMU_BIAS_TAU + kG^2) * eye(3);
Qc(CST_spfParam.COAST_ACC_BIAS_STATES, CST_spfParam.COAST_ACC_BIAS_STATES) = (2 * sA^2 / CST_spfParam.IMU_BIAS_TAU + kA^2) * eye(3);
Qds = dt / 2 * (Phi * Qc * Phi' + Qc);
Qhost = 9 * Qds + 1e-3 * eye(n);  Qhost(16, 16) = 0.1;  Qhost(17, 17) = 0.01;  Qhost(16, 17) = 0.02;  Qhost(17, 16) = 0.02;
Qhost(1:3, 1:3) = Qhost(1:3, 1:3) + 0.04 * eye(3);                      % the host's position noise
Qg = SPF_coastQ(Phi, Qhost);
inertial = false(n, 1);  inertial(1:15) = true;
diffIn = max(max(abs(Qg(inertial, inertial) - Qds(inertial, inertial))));
scaleIn = max(max(abs(Qds(inertial, inertial))));
okBlock = diffIn < 1e-12 * scaleIn;
okOther = isequal(Qg(~inertial, ~inertial), Qhost(~inertial, ~inertial)) && all(all(Qg(inertial, ~inertial) == 0)) && all(all(Qg(~inertial, inertial) == 0));
okSym = isequal(Qg, Qg') && min(eig(Qg(1:15, 1:15))) > -1e-18;
okPos = Qg(1, 1) < 1e-6 && Qg(1, 1) > 0;                               % the position block: only what the velocity noise integrates to
% the interval from Phi: half the interval halves the per-step density
Phi2 = Phi;  Phi2(1, 4) = dt / 2;
Qg2 = SPF_coastQ(Phi2, Qhost);
okDt = abs(Qg2(5, 5) / Qg(5, 5) - 0.5) < 1e-9;
% no position-from-velocity entry: the constant interval
Phi3 = Phi;  Phi3(1, 4) = 0;
Qg3 = SPF_coastQ(Phi3, Qhost);
Qds3 = CST_spfParam.COAST_Q_INTERVAL / 2 * (Phi3 * Qc * Phi3' + Qc);
okFallback = max(max(abs(Qg3(inertial, inertial) - Qds3(inertial, inertial)))) < 1e-12 * scaleIn;
fprintf('coast Q: inertial block %s, other states kept %s, symmetric and PSD %s, position from the velocity only %s, interval from Phi %s, fallback %s\n', ...
    pf(okBlock), pf(okOther), pf(okSym), pf(okPos), pf(okDt), pf(okFallback));
fprintf('   at 120 s the velocity block alone gives %.2f m of position sigma per axis from the datasheet, %.2f m from this test''s host Q\n', ...
    sqrt(Qds(4, 4) / dt) * 120^1.5 / sqrt(3), sqrt(Qhost(4, 4) / dt) * 120^1.5 / sqrt(3));
%% the floor
P = diag([1 4 9 0.01 0.01 0.01 1e-8 1e-8 1e-8 ones(1, n - 9) * 0.5]);  P(1, 2) = 0.5;  P(2, 1) = 0.5;
Pf = SPF_floorCov(P, [1 2 9 16], [2.0 1.0 0.01 0.1]);
okFloor = Pf(1, 1) == 4 && Pf(2, 2) == 4 && Pf(9, 9) == 1e-4 && Pf(16, 16) == 0.5 && Pf(3, 3) == 9 && Pf(1, 2) == 0.5 ...
    && min(eig(Pf(1:3, 1:3))) > 0 && isequal(Pf(4:8, 4:8), P(4:8, 4:8));
%% the repository configuration: the flag off, no floors
propTel = STRUCT_SPF.setPropTel(Phi, Qhost);
propTelOut = SPF_coastPropTel(propTel);
okRepo = ~CST_spfParam.COAST_Q_DATASHEET && isequal(propTelOut.accumQ, Qhost) && isequal(propTelOut.accumPhi, Phi) ...
    && isempty(CST_spfParam.COAST_SIGMA_FLOOR_STATES) && isequal(SPF_coastStartCov(P), P);
fprintf('floor: %s | repository configuration passes the host''s matrices unchanged: %s\n', pf(okFloor), pf(okRepo));
if okBlock && okOther && okSym && okPos && okDt && okFallback && okFloor && okRepo, fprintf('=== test_coastQ PASS ===\n');
else, fprintf('=== test_coastQ FAIL ===\n'); end
end

function s = pf(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end
