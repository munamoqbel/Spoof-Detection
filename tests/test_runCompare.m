%% test_runCompare.m
% Synthetic check of tools/runCompare.m (desktop post-processing tool, not
% codegen). A 300 s truth trajectory with speed and heading changes, two nav
% runs built from it with known errors, and the tool must read back:
%   - the truth arrays logged at different instants (trueVel shifted -0.3 s,
%     trueAtt +0.3 s against truePos) and the nav time-base lag (0.5 s)
%   - the velocity, attitude and bias errors at the outage start, the cubic
%     fit coefficients of the outage trace and the yaw-coupling budget term
%   - 2 Hz bias / P records with and without their epoch times, single
%     values at tOut, a truth velocity in ECEF (fallback to d(truePos)/dt),
%     no truth velocity at all, a truth yaw in degrees (exposed, not fixed)
%   - GNSS measurements applied 0.43 s late: the position and velocity trail
%     the attitude by 0.43 s and the regressions read -0.43 s
% Run from anywhere: octave --no-gui --quiet --eval test_runCompare
% (about a minute: the lag scans run on 30k samples). Prints PASS/FAIL lines.
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'tools'));
G = 9.80665;

tTrue = (0:0.01:300)';
lat0 = 0.98; lon0 = 0.2; h0 = 100;
a = 6378137; e2 = 6.69437999014e-3;
spd = 20 + 8 * sin(2 * pi * tTrue / 25);                            % accelerating (2 m/s^2 peaks): velocity lag observable
yaw = (170 + 20 * sin(2 * pi * tTrue / 60) + 90 * (tTrue >= 150)) * pi / 180;  % turning, crossing +-180 deg: attitude lag observable, wrap exercised
vN = spd .* cos(yaw); vE = spd .* sin(yaw); vD = zeros(size(tTrue));
lat = zeros(size(tTrue)); lon = lat; lat(1) = lat0; lon(1) = lon0;
for i = 2:numel(tTrue)
    L = lat(i-1);
    RN = a * (1 - e2) / (1 - e2 * sin(L)^2)^1.5;  RE = a / sqrt(1 - e2 * sin(L)^2);
    lat(i) = lat(i-1) + 0.5 * (vN(i-1) + vN(i)) * 0.01 / (RN + h0);        % trapezoidal: no half-sample lag
    lon(i) = lon(i-1) + 0.5 * (vE(i-1) + vE(i)) * 0.01 / ((RE + h0) * cos(L));
end
truePos = [lat lon h0 * ones(size(tTrue))];
trueVel = [vN vE vD];
trueAtt = [zeros(size(tTrue)) zeros(size(tTrue)) yaw];                 % yaw continuous here; wrapped where logged
biasTrue = [0.3e-3*G, -0.2e-3*G, 0.1e-3*G, 1/3600*pi/180, -0.5/3600*pi/180, 2/3600*pi/180];

function [C] = dcmT(e)
    cr = cos(e(1)); sr = sin(e(1)); cp = cos(e(2)); sp = sin(e(2)); cy = cos(e(3)); sy = sin(e(3));
    C = [cy*cp, cy*sp*sr - sy*cr, cy*sp*cr + sy*sr; sy*cp, sy*sp*sr + cy*cr, sy*sp*cr - cy*sr; -sp, cp*sr, cp*cr];
end

function [p] = mechCoast(Cnav, Ctrue, db, T)
    % stationary level truth with DCM Ctrue; the nav starts at Cnav, integrates the corrected gyro (0 - db_g)
    % and the corrected specific force (f_b - db_a), with gravity compensation, for T s at 100 Hz
    g = 9.80665;  dt = 0.01;
    fb = Ctrue' * [0; 0; -g];
    w = -db(4:6);  th = norm(w) * dt;  k = w / max(norm(w), eps);
    K = [0 -k(3) k(2); k(3) 0 -k(1); -k(2) k(1) 0];
    R = eye(3) + sin(th) * K + (1 - cos(th)) * K * K;             % body rotation per step
    C = Cnav;  v = zeros(3, 1);  p = zeros(3, 1);
    for i = 1:round(T / dt)
        C = C * R;
        a = C * (fb - db(1:3)) + [0; 0; g];
        v = v + a * dt;
        p = p + v * dt;
    end
end

function [M] = pDiagCat(c)
    M = zeros(numel(c), 22);
    for k = 1:numel(c), M(k, :) = diag(c{k})'; end
end

function [e] = nedDeltaT(p, q)
    a = 6378137; e2 = 6.69437999014e-3;
    L = p(:,1); h = p(:,3);
    RN = a * (1 - e2) ./ (1 - e2 * sin(L).^2).^1.5;  RE = a ./ sqrt(1 - e2 * sin(L).^2);
    e = [(p(:,1) - q(:,1)) .* (RN + h), (p(:,2) - q(:,2)) .* (RE + h) .* cos(L), -(p(:,3) - q(:,3))];
end

function [r, yawExp] = makeRun(name, lags, dv, acc, jerk, attErr, dBias, yawCouple, tTrue, truePos, trueVel, trueAtt, biasTrue)
    % lags(1): the nav at time t is the truth at t + lags(1). lags(2), lags(3): the logged truth velocity
    % and attitude arrays are shifted so that the tool must find lags(2) and lags(3) for them.
    a = 6378137; e2 = 6.69437999014e-3;
    t = (1:0.01:298)';                      % nav samples inside the truth span even with the lag
    tp = [interp1(tTrue, truePos(:,1), t + lags(1)) interp1(tTrue, truePos(:,2), t + lags(1)) interp1(tTrue, truePos(:,3), t + lags(1))];
    tv = [interp1(tTrue, trueVel(:,1), t + lags(1)) interp1(tTrue, trueVel(:,2), t + lags(1)) interp1(tTrue, trueVel(:,3), t + lags(1))];
    ta = [interp1(tTrue, trueAtt(:,1), t + lags(1)) interp1(tTrue, trueAtt(:,2), t + lags(1)) interp1(tTrue, trueAtt(:,3), t + lags(1))];
    rng_state = 42; randn('seed', rng_state);
    w = randn(numel(t) + 4000, 3);                                  % smooth aided position error (two 10 s moving averages)
    w = filter(ones(1, 1000) / 1000, 1, w);  w = filter(ones(1, 1000) / 1000, 1, w);
    eNed = w(4001:end, :);  eNed = eNed / std(eNed(:)) * 0.1;        % 0.1 m, slow: its derivative stays ~0.03 m/s
    vErr = ([eNed(2:end, :); eNed(end, :)] - [eNed(1, :); eNed(1:end-1, :)]) / 0.02;   % its derivative: the nav is consistent
    out = (t >= 100) & (t < 250);
    x = t(out) - 100;
    eNed(out, :) = x * dv + 0.5 * x.^2 * acc + x.^3 * jerk / 6;
    vErr(out, :) = repmat(dv, nnz(out), 1) + x * acc + 0.5 * x.^2 * jerk;
    yawExp = [0 0 0];
    if yawCouple                             % yaw error rotates the integrated specific force: -psi x (dp - v0 x)
        iOut = find(out, 1);
        psiD = -attErr(3);
        u = nedDeltaT(tp(out, :), repmat(tp(iOut, :), nnz(out), 1)) - x * tv(iOut, :);
        eNed(out, :) = eNed(out, :) + [psiD * u(:, 2), -psiD * u(:, 1), zeros(nnz(out), 1)];
        vErr(out, :) = vErr(out, :) + [psiD * (tv(out, 2) - tv(iOut, 2)), -psiD * (tv(out, 1) - tv(iOut, 1)), zeros(nnz(out), 1)];
        yawExp = [psiD * u(end, 2), -psiD * u(end, 1), 0];
    end
    if numel(lags) > 3                       % lags(4): GNSS measurements applied late: pos err = -lat v, vel err = -lat a
        aTrue = ([trueVel(2:end, :); trueVel(end, :)] - [trueVel(1, :); trueVel(1:end-1, :)]) / 0.02;
        ta4 = [interp1(tTrue, aTrue(:,1), t + lags(1)) interp1(tTrue, aTrue(:,2), t + lags(1)) interp1(tTrue, aTrue(:,3), t + lags(1))];
        eNed = eNed - lags(4) * tv;
        vErr = vErr - lags(4) * ta4;
    end
    L = tp(:,1); h = tp(:,3);
    RN = a * (1 - e2) ./ (1 - e2 * sin(L).^2).^1.5;  RE = a ./ sqrt(1 - e2 * sin(L).^2);
    r.name = name;  r.t = t;
    r.navPos = [tp(:,1) + eNed(:,1) ./ (RN + h), tp(:,2) + eNed(:,2) ./ ((RE + h) .* cos(L)), tp(:,3) - eNed(:,3)];
    r.navVel = tv + vErr;
    r.navAtt = ta + repmat(attErr, numel(t), 1);
    r.navAtt(:, 3) = atan2(sin(r.navAtt(:, 3)), cos(r.navAtt(:, 3)));          % logged wrapped to +-pi
    r.biasEst = repmat(biasTrue + dBias, numel(t), 1);
    r.tTrue = tTrue; r.truePos = truePos; r.biasTrue = biasTrue;
    sV = lags(1) - lags(2);  sA = lags(1) - lags(3);      % logged(tTrue) = truth(tTrue + s): nav(t) = logged(t + lags(1) - s)
    r.trueVel = [interp1(tTrue, trueVel(:,1), tTrue + sV) interp1(tTrue, trueVel(:,2), tTrue + sV) interp1(tTrue, trueVel(:,3), tTrue + sV)];
    r.trueAtt = [interp1(tTrue, trueAtt(:,1), tTrue + sA) interp1(tTrue, trueAtt(:,2), tTrue + sA) interp1(tTrue, trueAtt(:,3), tTrue + sA)];
    r.trueAtt(:, 3) = atan2(sin(r.trueAtt(:, 3)), cos(r.trueAtt(:, 3)));        % logged wrapped to +-pi
    r.pDiag = repmat([ones(1,3) 0.05^2*ones(1,3) 1e-6*ones(1,3) (0.2e-3*9.80665)^2*ones(1,3) (1/3600*pi/180)^2*ones(1,3) ones(1,7)], numel(t), 1);
end

dvA = [0.10 -0.05 0.02]; jerkA = [0 2e-5 0];
attErrA = [0.5e-3 -0.3e-3 1e-3];  dBiasA = [-0.6e-3*G 0.1e-3*G 0 0.3/3600*pi/180 0 -1/3600*pi/180];
yaw100 = interp1(tTrue, trueAtt(:,3), 100.5);                                  % truth heading at the nav instant of tOut
CtrueA0 = dcmT([0 0 yaw100]);  CnavA0 = dcmT([0 0 yaw100] + attErrA);
dCA0 = CnavA0 * CtrueA0';
psiA0 = [dCA0(2,3) - dCA0(3,2); dCA0(3,1) - dCA0(1,3); dCA0(1,2) - dCA0(2,1)] / 2;
accA = (CnavA0 * (-dBiasA(1:3)') + G * [psiA0(2); -psiA0(1); 0])';           % so the trace matches the state errors
dvB = [0.30  0.10 0.00]; accB = [0.006 -0.003 0]; jerkB = [0 0 0];
lagsA = [0.5 0.2 0.8];                     % position, velocity and attitude truth logged at different instants
lagsB = [0.5 0.5 0.5];
runA = makeRun('A', lagsA, dvA, accA, jerkA, attErrA, dBiasA, false, tTrue, truePos, trueVel, trueAtt, biasTrue);
[runB, yawExpB] = makeRun('B', lagsB, dvB, accB, jerkB, [0.2e-3  0.4e-3 -2e-3], [ 0.2e-3*G 0.3e-3*G 0 1.0/3600*pi/180 0  0.5/3600*pi/180], true, tTrue, truePos, trueVel, trueAtt, biasTrue);

% run B logged per 2 Hz epoch: bias on tKf, P as a cell of 22 x 22 matrices
runB.tKf = (1:0.5:298)';
runB.biasEst = repmat(runB.biasEst(1, :), numel(runB.tKf), 1);
pRow = runB.pDiag(1, :);
runB.pDiag = cell(numel(runB.tKf), 1);
for k = 1:numel(runB.tKf), runB.pDiag{k} = diag(pRow); end
res = runCompare(runA, runB, 100, 250, [10 95; 255 295]);

T = res.run1.T;
expDriftA = dvA * T + 0.5 * accA * T^2 + jerkA * T^3 / 6;
chk = struct();
chk.truthShifts  = abs(res.run1.sV + 0.3) < 0.011 && abs(res.run1.sA - 0.3) < 0.011;                 % logged truth arrays shifted
chk.lagsA        = abs(res.run1.tauP - 0.5) < 0.011 && abs(res.run1.tauV - 0.5) < 0.011 && abs(res.run1.tauA - 0.5) < 0.011 && res.run1.tau == res.run1.tauA;
chk.lagsB        = abs(res.run2.sV) < 0.011 && abs(res.run2.sA) < 0.011 && abs(res.run2.tau - 0.5) < 0.011 && abs(res.run2.tauP - 0.5) < 0.011;
chk.navConsistent = abs(res.run1.sNV) < 0.011 && abs(res.run2.sNV) < 0.011 && abs(res.run1.kP) < 0.03 && abs(res.run1.kV) < 0.03;   % no latency
chk.dvA          = all(abs(res.run1.dv' - dvA) < 0.02) && all(abs(res.run1.slope10' - dvA) < 0.08);
chk.dvImpliedA   = all(abs(res.run1.dvImplied' - dvA) < 0.03);                                       % secant corrected for the bias+tilt growth
chk.driftA       = all(abs(res.run1.drift - expDriftA) < 1.0);
chk.fitB         = all(abs(res.run1.coef(:,3)' - dvA) < 0.02);
chk.fitC         = all(abs(res.run1.coef(:,2)' - 0.5*accA) < 2e-4);
chk.fitD         = all(abs(res.run1.coef(:,1)' - jerkA/6) < 5e-7);
chk.dvB          = all(abs(res.run2.dv' - dvB) < 0.02);
chk.attitudeA    = abs(res.run1.attErr(3) - 1e-3) < 1e-5 && abs(res.run1.psi(3) + 1e-3) < 2e-5;
chk.biasA        = abs(res.run1.db(1) + 0.6e-3*G) < 1e-9;
chk.biasSigmaB   = abs(res.run2.db(1) - 0.2e-3*G) < 1e-9 && abs(res.run2.sig(4) - 0.05) < 1e-12;    % 2 Hz bias and cell P
chk.yawTermB     = all(abs(res.run2.yawTerm(1:2)' - yawExpB(1:2)) < 0.1) && norm(yawExpB) > 1;       % yaw coupling term (horizontal)
% bias, tilt and gyro rows against a mechanised stationary coast from the same tOut state (signs and conventions)
iO = find(runA.t <= 100, 1, 'last');
CnavA  = dcmT(runA.navAtt(iO, :));
CtrueA = dcmT(interp1(tTrue, trueAtt, 100 + res.run1.tau));
pMech  = mechCoast(CnavA, CtrueA, res.run1.db, res.run1.T);
rows   = res.run1.biasTerm + res.run1.tiltTerm + res.run1.gyroTerm;
fprintf('mechanised coast N %.2f E %.2f vs budget rows N %.2f E %.2f m\n', pMech(1:2), rows(1:2));
chk.mechanised   = all(abs(pMech(1:2) - rows(1:2)) < 0.01 * norm(rows(1:2)) + 0.5);
names = fieldnames(chk);  ok = true;
for i = 1:numel(names)
    if ~chk.(names{i}), fprintf('   FAILED sub-check: %s\n', names{i}); ok = false; end
end
fprintf('\nexpected drift A: N %.2f E %.2f D %.2f;  yaw term B expected N %.2f E %.2f, got N %.2f E %.2f\n', ...
    expDriftA, yawExpB(1:2), res.run2.yawTerm(1:2));
if ok, fprintf('=== runCompare synthetic check PASS ===\n'); else, fprintf('=== runCompare synthetic check FAIL ===\n'); end

% single-run mode
res1 = runCompare(runA, [], 100, 250);
% 2 Hz records without tKf (evenly spread), and a single P at tOut
runC = runA;  runC.name = 'C no tKf';
runC.biasEst = repmat(runA.biasEst(1, :), 595, 1);          % 595 records over 297 s = 0.5 s apart
runC.pDiag = diag(runA.pDiag(1, :));                        % one 22 x 22 P at tOut
resC = runCompare(runC, [], 100, 250);
% a single bias row at tOut
runD = runA;  runD.name = 'D bias row';  runD.biasEst = runA.biasEst(1, :);  runD.pDiag = runA.pDiag(1, :);
resD = runCompare(runD, [], 100, 250);
okC = abs(resC.run1.db(1) + 0.6e-3*G) < 1e-9 && abs(resC.run1.sig(4) - 0.05) < 1e-12 ...
   && abs(resD.run1.db(1) + 0.6e-3*G) < 1e-9 && abs(resD.run1.sig(10)/(9.80665e-3) - 0.2) < 1e-9;
if okC, fprintf('=== no-tKf / single-value check PASS ===\n'); else, fprintf('=== no-tKf / single-value check FAIL ===\n'); end
% truth velocity logged in ECEF (wrong frame): the tool must warn and fall back to d(truePos)/dt
runE = runA;  runE.name = 'E ecef trueVel';
L = runA.truePos(:,1); l = runA.truePos(:,2);
runE.trueVel = [(-sin(L).*cos(l)).*trueVel(:,1) + (-sin(l)).*trueVel(:,2) + (-cos(L).*cos(l)).*trueVel(:,3), ...
                (-sin(L).*sin(l)).*trueVel(:,1) + ( cos(l)).*trueVel(:,2) + (-cos(L).*sin(l)).*trueVel(:,3), ...
                ( cos(L)).*trueVel(:,1)                                   + (-sin(L)).*trueVel(:,3)];
resE = runCompare(runE, [], 100, 250);
% no truth velocity at all
runF = runA;  runF.name = 'F no trueVel';  runF = rmfield(runF, 'trueVel');
resF = runCompare(runF, [], 100, 250);
% truth yaw in degrees: the yaw-vs-course line must expose it
runG = runA;  runG.name = 'G yaw in degrees';  runG.trueAtt(:,3) = runG.trueAtt(:,3) * 180 / pi;
resG = runCompare(runG, [], 100, 250);
% GNSS measurements applied 0.43 s late: position and velocity trail the attitude, the regressions read the latency
runH = makeRun('H latency 0.43 s', [0.5 0.5 0.5 0.43], dvA, accA, jerkA, [0.5e-3 -0.3e-3 1e-3], [-0.6e-3*G 0.1e-3*G 0 0.3/3600*pi/180 0 -1/3600*pi/180], false, tTrue, truePos, trueVel, trueAtt, biasTrue);
resH = runCompare(runH, [], 100, 250, [10 95; 255 295]);
okH = abs(resH.run1.tauA - 0.5) < 0.011 && abs(resH.run1.tauP - 0.07) < 0.011 && abs(resH.run1.tauV - 0.07) < 0.011 ...
   && abs(resH.run1.kP + 0.43) < 0.02 && resH.run1.r2P > 0.95 && abs(resH.run1.kV + 0.43) < 0.03 && abs(resH.run1.sNV) < 0.011;
if okH, fprintf('=== latency check PASS ===\n'); else, fprintf('=== latency check FAIL ===\n'); end
% pDiag given as states x records, with the wrong record count, and under a misspelt field name
runP = runB;  runP.name = 'P transposed';  runP.pDiag = pDiagCat(runB.pDiag)';           % 22 x K
resP = runCompare(runP, [], 100, 250, [10 95; 255 295]);
runQ = runB;  runQ.name = 'Q pDiag count mismatch';  runQ.pDiag = [pDiagCat(runB.pDiag); ones(3, 22)];
resQ = runCompare(runQ, [], 100, 250, [10 95; 255 295]);
runR = runB;  runR.name = 'R misspelt pDig';  runR.pDig = runB.pDiag;  runR = rmfield(runR, 'pDiag');
resR = runCompare(runR, [], 100, 250, [10 95; 255 295]);
okP = abs(resP.run1.sig(4) - 0.05) < 1e-12 && ~any(isfinite(resQ.run1.sig)) && ~any(isfinite(resR.run1.sig));
if okP, fprintf('=== pDiag shape / mismatch / misspelt check PASS ===\n'); else, fprintf('=== pDiag shape / mismatch / misspelt check FAIL ===\n'); end
okG = resG.run1.tau == resG.run1.tauP && abs(resG.run1.tauP - 0.5) < 0.011;
if okG, fprintf('=== unusable attitude fallback check PASS ===\n'); else, fprintf('=== unusable attitude fallback check FAIL ===\n'); end
okE = strcmp(resE.run1.velSource, 'position') && all(abs(resE.run1.dv' - dvA) < 0.02) ...
   && strcmp(resF.run1.velSource, 'position') && all(abs(resF.run1.dv' - dvA) < 0.02) ...
   && strcmp(res1.run1.velSource, 'logged');
if okE, fprintf('=== truth velocity fallback check PASS ===\n'); else, fprintf('=== truth velocity fallback check FAIL ===\n'); end
if abs(res1.run1.tau - 0.5) < 0.011 && all(abs(res1.run1.dv' - dvA) < 0.02), fprintf('=== default window check PASS ===\n'); else, fprintf('=== default window check FAIL ===\n'); end
