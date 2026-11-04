function test_measNoiseCheck()
%% test_measNoiseCheck.m
% tools/measNoiseCheck on synthetic spfMeas records (8 range rows at 5 m, 8 rate rows, a pressure row): on epochs 30 to
% 45 the host scales five range rows to 15 m and leaves two rows out; the old interface (one R) must read the gate's
% R as scaled and the rows as missing; the new one (setKfMeasBoth layout, the gate's R nominal, the host's scaled with
% the left-out rows at 1e6) must read the gate's R at the nominal and the host's as scaled with two rows left out.
% Run from anywhere: octave --no-gui --quiet --eval test_measNoiseCheck
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'tools'));
randn('seed', 4);  rand('seed', 4);
K = 60;  n = 22;  mMax = 51;  tKf = 0.5 * (0:K-1)';
el = 0.2 + 0.8 * rand(8, 1);  az = 2 * pi * rand(8, 1);
u = [cos(el) .* cos(az), cos(el) .* sin(az), -sin(el)];
H = zeros(17, n);  Rd = zeros(17, 1);
for j = 1:8, H(j, 1:3) = -u(j, :);  H(j, 16) = 1;  Rd(j) = 25; end
for j = 1:8, H(8 + j, 4:6) = -u(j, :);  H(8 + j, 17) = 1;  Rd(8 + j) = 0.01; end
H(17, 3) = 1;  Rd(17) = 2.25;
P = eye(n);
old = cell(K, 1);  both = cell(K, 1);
for k = 1:K
    y = 0.5 * randn(17, 1);  R = diag(Rd);  Rh = R;  rows = 1:17;
    if k >= 30 && k <= 45
        y(1:8) = y(1:8) + 20 * H(1:8, 1:3) * [1 0 0]';
        Rh(1, 1) = 1e6;  Rh(2, 2) = 1e6;                      % the host leaves rows 1 and 2 out
        for j = 3:7, Rh(j, j) = 225; end                      % and scales five rows to 15 m
        rows = 3:17;                                          % the old interface: two rows gone, R scaled
        R = Rh;
    end
    m = numel(rows);
    Hk = zeros(mMax, n);  Hk(1:m, :) = H(rows, :);  yk = zeros(mMax, 1);  yk(1:m) = y(rows);
    Rk = zeros(mMax);  Rk(1:m, 1:m) = R(rows, rows);  Sk = zeros(mMax);  Sk(1:m, 1:m) = H(rows, :) * P * H(rows, :)' + R(rows, rows);
    old{k} = struct('innovation', yk, 'innovationCov', Sk, 'obsMatrix', Hk, 'measNoiseCov', Rk, 'numMeas', m, ...
        'priorState', zeros(n, 1), 'postState', zeros(n, 1), 'postCov', P);
    Hb = zeros(mMax, n);  Hb(1:17, :) = H;  yb = zeros(mMax, 1);  yb(1:17) = y;
    Rn = zeros(mMax);  Rn(1:17, 1:17) = diag(Rd);  Sn = zeros(mMax);  Sn(1:17, 1:17) = H * P * H' + diag(Rd);
    Rhk = zeros(mMax);  Rhk(1:17, 1:17) = Rh;  Shk = zeros(mMax);  Shk(1:17, 1:17) = H * P * H' + Rh;
    both{k} = struct('innovation', yb, 'innovationCov', Sn, 'obsMatrix', Hb, 'measNoiseCov', Rn, 'numMeas', 17, ...
        'priorState', zeros(n, 1), 'postState', zeros(n, 1), 'postCov', P, 'innovationCovHost', Shk, 'measNoiseCovHost', Rhk);
end
fprintf('--- the old interface: one R, scaled, a row left out ---\n');
r1 = measNoiseCheck(old, tKf, [tKf(28) tKf(48)]);
fprintf('--- the new interface: the gate''s R nominal, the host''s beside ---\n');
r2 = measNoiseCheck(both, tKf, [tKf(28) tKf(48)]);
ok1 = r1.gateScaled && ~r1.hasHost && r1.nScaled(35) == 5 && abs(r1.maxSqrtR(35) - 15) < 1e-9 && r1.nRange(35) == 6 && r1.nRange(20) == 8 ...
   && r1.nScaled(20) == 0 && abs(r1.rNominal(1) - 25) < 1e-9 && r1.nMissing == 16 && isnan(r1.nScaledHost(35));
ok2 = ~r2.gateScaled && r2.hasHost && r2.hostScaled && r2.nScaled(35) == 0 && abs(r2.maxSqrtR(35) - 5) < 1e-9 && r2.nRange(35) == 8 ...
   && r2.nScaledHost(35) == 5 && r2.nRejectedHost(35) == 2 && abs(r2.maxSqrtRHost(35) - 15) < 1e-9 && r2.nRejectedHost(20) == 0 && r2.nMissing == 0;
if ok1 && ok2, fprintf('=== test_measNoiseCheck PASS ===\n');
else, fprintf('ok1 %d ok2 %d\n=== test_measNoiseCheck FAIL ===\n', ok1, ok2); end
end
