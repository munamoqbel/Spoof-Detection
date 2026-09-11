%% test_error_handlers.m
% Numerical exception handlers of the gate (audit of every division,
% solve, sqrt, index and NaN path on the runtime path):
%   1. SPF_cpiMonitor: a non-PD S (duplicated row, zero noise) is dropped with
%      solveFault and a finite q; a PD S still alarms on a bias, and the
%      Cholesky form equals the explicit f'(S\g) form.
%   2. SPF_revalidation: an indefinite residual covariance cannot pass and is
%      flagged; a PD one reproduces r'(S\r) and the REVAL_MIN_MEAS floor.
%   3. numMeas > MAX_MEAS is clamped in both kfMeas builders and flagged.
%   4. a non-finite input epoch resets the gate (inputFault), nothing
%      enters the persistent state, and the gate runs on cleanly.
%   5. the startup anchor is stamped with epoch 1 (a latch onto it reports
%      eventAnchorEpoch = 1, never the 'none' sentinel 0).
% Run from the repo root (MATLAB, or Octave with tools/octave_shim).

n    = CST_gnssHybrid.NO_STATES;
mMax = double(CST_spfParam.MAX_MEAS);
N    = double(CST_spfParam.WINDOW_LENGTH);
rng(7);

%% 1. SPF_cpiMonitor with a singular / PD innovation covariance
fprintf('SPF_cpiMonitor exception handler (non-PD S) ...\n');
m = 4;  P = eye(n);
H = zeros(m, n); H(:, 1:3) = randn(m, 3); H(4, :) = H(3, :);      % rows 3 and 4 identical
innBuf = zeros(mMax, N); sBuf = zeros(mMax, mMax, N); hBuf = zeros(mMax, n, N);
numBuf = uint8(m) * ones(N, 1, 'uint8');
bias = zeros(n, 1); bias(1) = 5.0;                                 % 5 m along axis 1
for k = 1:N
    Rsing = diag([1 1 0 0]);                                       % zero noise on the duplicated rows
    S = H * P * H' + Rsing;
    innBuf(1:m, k) = H * bias + 0.1 * randn(m, 1);
    sBuf(1:m, 1:m, k) = S;  hBuf(1:m, :, k) = H;
end
[alarmSing, qSing, ~, faultSing] = SPF_cpiMonitor(innBuf, sBuf, hBuf, numBuf, 1);
okSing = faultSing && isfinite(qSing) && (qSing == 0) && ~alarmSing;
fprintf('  singular S: solveFault = %d, q = %g (finite, epochs dropped), alarm = %d -> %s\n', ...
    faultSing, qSing, alarmSing, pf(okSing));

qRef = 0;
for k = 1:N
    S = H * P * H' + eye(m);  sBuf(1:m, 1:m, k) = S;
    g = innBuf(1:m, k); f = H(:, 1);
    qRef = qRef + (f' * (S \ g))^2 / (f' * (S \ f));               % explicit Eq. 17/20/29 form
end
[alarmPD, qPD, ~, faultPD] = SPF_cpiMonitor(innBuf, sBuf, hBuf, numBuf, 1);
okPD = alarmPD && ~faultPD && abs(qPD - qRef) < 1e-9 * max(1, qRef);
fprintf('  PD S, 5 m bias: alarm = %d, |q_chol - q_explicit| = %.1e -> %s\n\n', alarmPD, abs(qPD - qRef), pf(okPD));
ok1 = okSing && okPD;

%% 2. SPF_revalidation with an indefinite / PD residual covariance
fprintf('SPF_revalidation exception handler (non-PD residual covariance) ...\n');
m = 5;
Hr = zeros(mMax, n); Hr(1:m, 1:3) = randn(m, 3);
Rr = zeros(mMax); Rr(1:m, 1:m) = eye(m);
inn = zeros(mMax, 1); inn(1:m) = 0.5 * randn(m, 1);
Pc = eye(n); Pc(1, 1) = -50;                                       % indefinite along an observed direction
[passBad, qBad, faultBad] = SPF_revalidation(inn, Hr, Rr, m, zeros(n, 1), Pc);
okBad = ~passBad && faultBad && (qBad == 0);
fprintf('  indefinite P_C: passed = %d, solveFault = %d, q = %g -> %s\n', passBad, faultBad, qBad, pf(okBad));
Pc = 4 * eye(n);
Sr = Hr(1:m, :) * Pc * Hr(1:m, :)' + eye(m);
qRefR = inn(1:m)' * (Sr \ inn(1:m));
[passGood, qGood, faultGood] = SPF_revalidation(inn, Hr, Rr, m, zeros(n, 1), Pc);
okGood = passGood && ~faultGood && abs(qGood - qRefR) < 1e-9 * max(1, qRefR);
[passFew, ~, ~] = SPF_revalidation(inn, Hr, Rr, 1, zeros(n, 1), Pc);    % pressure-only style epoch
okFew = ~passFew;
fprintf('  PD P_C: passed = %d, |q_chol - q_explicit| = %.1e -> %s | numMeas = 1 cannot pass -> %s\n\n', ...
    passGood, abs(qGood - qRefR), pf(okGood), pf(okFew));
ok2 = okBad && okGood && okFew;

%% 3. numMeas beyond MAX_MEAS
fprintf('numMeas > MAX_MEAS clamp ...\n');
mBig = mMax + 1;
yBig = randn(mBig, 1); Hbig = randn(mBig, n); Rbig = eye(mBig);
kfMeasBig = STRUCT_SPF.kfMeasFromUpdate(yBig, Hbig, Rbig, mBig, zeros(n, 1), eye(n), zeros(n, 1), eye(n));
okFrom = (kfMeasBig.numMeas == CST_spfParam.MAX_MEAS) && kfMeasBig.numMeasClamped ...
    && isequal(size(kfMeasBig.innovation), [mMax 1]) && isequal(size(kfMeasBig.obsMatrix), [mMax n]) ...
    && isequal(size(kfMeasBig.innovationCov), [mMax mMax]);
kfMeasSet = STRUCT_SPF.setKfMeas(zeros(mMax, 1), eye(mMax), zeros(mMax, n), eye(mMax), 35, zeros(n, 1), zeros(n, 1), eye(n));
okSet = (kfMeasSet.numMeas == CST_spfParam.MAX_MEAS) && kfMeasSet.numMeasClamped;
propTelI = STRUCT_SPF.zeroPropTel;
tel = SPF_gate(kfMeasBig, propTelI, true, true);           % must run, and report the clamp
okRun = tel.info.numMeasClamped && (tel.info.mode == CST_spfMode.NOMINAL);
ok3 = okFrom && okSet && okRun;
fprintf('  kfMeasFromUpdate(31 rows) -> numMeas %d, clamped %d, fixed layout %d | setKfMeas(35) -> %d | gate runs, flags %d -> %s\n\n', ...
    kfMeasBig.numMeas, kfMeasBig.numMeasClamped, okFrom, kfMeasSet.numMeas, tel.info.numMeasClamped, pf(ok3));

%% 4./5. non-finite input epoch, and the startup anchor stamp
fprintf('Non-finite input and startup anchor ...\n');
prm_boot = kujur_params();
[~, ~, H_all, ~, ~, Phi, Q, ~, z_all, V, P0] = generate_test_data(prm_boot, 5, 0.0, 0.0, [0 0 1]);
propTel = STRUCT_SPF.setPropTel(Phi, Q);
mm = size(z_all, 1);
kf_x = zeros(n, 1); kf_P = P0;
faultSeen = false; alarmsAfter = 0; finiteOK = true; modeOK = true;
for k = 1:60
    [kf_x, kf_P, y, ~, ~, xp, xpP] = kalman_update_step(kf_x, kf_P, z_all(:, k), H_all(:, :, k), propTel, V);
    if k == 21, y(3) = NaN; end                                    % one corrupted epoch
    kfMeas = STRUCT_SPF.kfMeasFromUpdate(y, H_all(:, :, k), V, mm, xp, xpP, kf_x, kf_P);
    tel = SPF_gate(kfMeas, propTel, true, k == 1);
    if k == 21
        faultSeen = tel.info.inputFault && (tel.info.mode == CST_spfMode.NOMINAL) && ~tel.nav.applyCorrection;
    elseif k > 21
        alarmsAfter = alarmsAfter + (tel.info.ssAlarm || tel.info.cpiAlarm);
        finiteOK = finiteOK && isfinite(tel.info.maxProtectionLevel) && all(isfinite(tel.nav.sigmaPosition)) && ~tel.info.inputFault;
        modeOK = modeOK && (tel.info.mode == CST_spfMode.NOMINAL);
    end
end
ok4 = faultSeen && (alarmsAfter == 0) && finiteOK && modeOK;
fprintf('  NaN at epoch 21: inputFault = %d | epochs 22-60: alarms %d, finite %d, NOMINAL %d -> %s\n', ...
    faultSeen, alarmsAfter, finiteOK, modeOK, pf(ok4));

kf_x = zeros(n, 1); kf_P = P0; anchorEpoch = -1; latched = false;
for k = 1:6
    z = z_all(:, k); if k == 4, z = z + 50.0; end                  % gross spoof on all rows at epoch 4
    [kf_x, kf_P, y, ~, ~, xp, xpP] = kalman_update_step(kf_x, kf_P, z, H_all(:, :, k), propTel, V);
    kfMeas = STRUCT_SPF.kfMeasFromUpdate(y, H_all(:, :, k), V, mm, xp, xpP, kf_x, kf_P);
    tel = SPF_gate(kfMeas, propTel, true, k == 1);
    if tel.info.eventLatched && ~latched
        latched = true; anchorEpoch = double(tel.info.eventAnchorEpoch);
    end
end
ok5 = latched && (anchorEpoch == 1);
fprintf('  50 m spoof at epoch 4 before any clean close: latched = %d, eventAnchorEpoch = %d (expect 1) -> %s\n\n', ...
    latched, anchorEpoch, pf(ok5));

if ok1 && ok2 && ok3 && ok4 && ok5
    fprintf('=== test_error_handlers PASS ===\n');
else
    fprintf('=== test_error_handlers FAIL ===\n');
end

function s = pf(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end
