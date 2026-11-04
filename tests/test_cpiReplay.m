function test_cpiReplay()
%% test_cpiReplay.m
% tools/cpiReplay on synthetic spfMeas records (8 range rows, 8 rate rows, a pressure row; S = H P H' + R with a 1 m
% position sigma): a 12 m north offset in the range rows on epochs 30 to 45 must be read back as dHat 12 m with the
% sigma of the projection, xi = (12 - B)/sigma, the window sums must equal SPF_cpiMonitor's and the gate log's when
% the log holds them, and must not when the log holds other numbers. Run from the repo root (the stand-in matrixInv):
% octave --no-gui --quiet --eval test_cpiReplay
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'tools'));  addpath(fullfile(here, '..', 'deploy'));  addpath(fullfile(here, '..'));
randn('seed', 3);  rand('seed', 3);
K = 60;  n = 22;  mMax = 51;  tKf = 0.5 * (0:K-1)';
P = zeros(n);  P(1:3, 1:3) = eye(3) * 1.0;  P(4:6, 4:6) = eye(3) * 0.01;  P(16, 16) = 4;  P(17, 17) = 0.01;
el = 0.2 + 0.8 * rand(8, 1);  az = 2 * pi * rand(8, 1);
u = [cos(el) .* cos(az), cos(el) .* sin(az), -sin(el)];
H = zeros(17, n);  Rd = zeros(17, 1);
for j = 1:8, H(j, 1:3) = -u(j, :);  H(j, 16) = 1;  Rd(j) = 25; end
for j = 1:8, H(8 + j, 4:6) = -u(j, :);  H(8 + j, 17) = 1;  Rd(8 + j) = 0.01; end
H(17, 3) = 1;  Rd(17) = 2.25;
R = diag(Rd);  S = H * P * H' + R;
log = cell(K, 1);
for k = 1:K
    d = [0 0 0];
    if k >= 30 && k <= 45, d = [12 0 0]; end
    y = H(:, 1:3) * d' + 0.1 * randn(17, 1);
    Hk = zeros(mMax, n);  Hk(1:17, :) = H;  yk = zeros(mMax, 1);  yk(1:17) = y;
    Sk = zeros(mMax);  Sk(1:17, 1:17) = S;  Rk = zeros(mMax);  Rk(1:17, 1:17) = R;
    log{k} = struct('innovation', yk, 'innovationCov', Sk, 'obsMatrix', Hk, 'measNoiseCov', Rk, 'numMeas', 17, ...
        'priorState', zeros(n, 1), 'postState', zeros(n, 1), 'postCov', P);
end
B = [5 5 8];
% a gate log whose cpiRatio is SPF_cpiMonitor on these records, and one that holds other numbers
spfLog = SPF_logInit(K);  spfLog.count = K;
N = double(CST_spfParam.WINDOW_LENGTH);  T = CST_spfParam.CPI_THRESHOLD;
for k = N:K
    innB = zeros(mMax, N);  covB = zeros(mMax, mMax, N);  obsB = zeros(mMax, n, N);  numB = zeros(N, 1);
    for i = 1:N, r = log{k - N + i};  innB(:, i) = r.innovation;  covB(:, :, i) = r.innovationCov;  obsB(:, :, i) = r.obsMatrix;  numB(i) = r.numMeas; end
    for a = 1:3, [~, q] = SPF_cpiMonitor(innB, covB, obsB, numB, a, B(a));  spfLog.cpiRatio(k, a) = q / T; end
end
spfOther = spfLog;  spfOther.cpiRatio = 0.5 * spfLog.cpiRatio;
fprintf('--- 12 m north on epochs 30 to 45, B = [5 5 8], the log from the monitor itself ---\n');
r1 = cpiReplay(log, tKf, spfLog, [tKf(25) tKf(50)], B);
fprintf('--- the same against a log holding half the numbers, B from CST_spfParam ---\n');
r2 = cpiReplay(log, tKf, spfOther, [tKf(25) tKf(50)]);
f = H(:, 1);  s2 = f' * (S \ f);  sigExp = 1 / sqrt(s2);  sigRexp = 1 / sqrt(f' * (R \ f));
xiExp = (12 - 5) / sigExp;
qExp = sum(r1.xi(36:45, 1).^2) / T;
ok1 = abs(r1.dHat(35, 1) - 12) < 0.3 && isfinite(r1.dHat(35, 2)) && abs(r1.sig(35, 1) - sigExp) < 1e-9 && abs(r1.sigR(35, 1) - sigRexp) < 1e-9 ...   % E reads the leak of the N step through the sky, not 0
   && abs(r1.xi(35, 1) - xiExp) < 0.1 && abs(r1.qReplay(45, 1) - qExp) < 1e-6 && abs(r1.qReplay(45, 1) - spfLog.cpiRatio(45, 1)) < 1e-9 ...
   && r1.qReplay(25, 1) == 0 && r1.xi(20, 1) == 0 && r1.matches == 1 && abs(r1.sigP(35, 1) - 1) < 1e-9 && isequal(r1.window, [tKf(25) tKf(50)]) ...
   && abs(r1.dHat(20, 1)) < 0.3;
ok2 = r2.matches == 0 && abs(r2.xi(35, 1) - 12 / sigExp) < 0.1 && abs(r2.dHat(35, 1) - 12) < 0.3;   % B = 0 from the deploy constants
if ok1 && ok2, fprintf('=== test_cpiReplay PASS ===\n');
else, fprintf('ok1 %d ok2 %d\n=== test_cpiReplay FAIL ===\n', ok1, ok2); end
end
