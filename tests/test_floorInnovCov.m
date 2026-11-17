function test_floorInnovCov()
%% test_floorInnovCov.m
% SPF_floorInnovCov: the host's S with the gate's covariance carried in, S + H dP H' on the live rows. (1) With the
% gate's P equal to the host's and no clock floor the matrix comes back as it went, padding included. (2) With 2.5 m
% floors on the position the live block must equal S + H diag(dP) H', dP the floor variance over the host's, the
% padding untouched and the result symmetric. (3) A clock floor of 30 m on a clock sigma of 2 m must add the clock
% term; a clock state beyond the layout is ignored. (4) Through SPF_cpiMonitor: 30 m on every range row, a clock jump,
% crosses the threshold on the host's S and stays below a twentieth of it with the floors in; a 20 m north step still
% crosses with them. (5) Rows beyond numMeas are not touched. (6) In the repository configuration, no floors and no
% clock floor, the function passes S unchanged, so the test references stand.
% Run from the repo root (MATLAB, or Octave with tools/octave_shim).
randn('seed', 5);  rand('seed', 5);
n = CST_gnssHybrid.NO_STATES;  mMax = double(CST_spfParam.MAX_MEAS);
P = zeros(n);  P(1:3, 1:3) = eye(3) * 1.0;  P(4:6, 4:6) = eye(3) * 0.01;  P(16, 16) = 4;  P(17, 17) = 0.01;
el = 0.2 + 0.8 * rand(8, 1);  az = 2 * pi * rand(8, 1);
u = [cos(el) .* cos(az), cos(el) .* sin(az), -sin(el)];
m = 17;  H = zeros(mMax, n);  Rd = zeros(mMax, 1);
for j = 1:8, H(j, 1:3) = -u(j, :);  H(j, 16) = 1;  Rd(j) = 25; end
for j = 1:8, H(8 + j, 4:6) = -u(j, :);  H(8 + j, 17) = 1;  Rd(8 + j) = 0.01; end
H(17, 3) = 1;  Rd(17) = 25;
R = diag(Rd);  S = H * P * H' + R;  S = (S + S') / 2;                        % the padding rows and columns are zero
live = 1:m;
%% (1) inert
S1 = SPF_floorInnovCov(S, H, uint8(m), P, P, zeros(1, 0), 0.0);
ok1 = isequal(S1, S);
%% (2) position floors
Pg = SPF_floorCov(P, [1 2 3], [2.5 2.5 2.5]);
S2 = SPF_floorInnovCov(S, H, uint8(m), P, Pg, zeros(1, 0), 0.0);
dP = zeros(n, 1);  dP(1:3) = 2.5^2 - 1;
S2exp = S;  S2exp(live, live) = S(live, live) + H(live, :) * diag(dP) * H(live, :)';
ok2 = max(max(abs(S2 - S2exp))) < 1e-9 && isequal(S2(m+1:end, :), S(m+1:end, :)) && isequal(S2(:, m+1:end), S(:, m+1:end)) ...
   && max(max(abs(S2 - S2'))) < 1e-12;
%% (3) the clock floor, and a clock state outside the layout
S3 = SPF_floorInnovCov(S, H, uint8(m), P, P, [16 n + 3], 30.0);
dPc = zeros(n, 1);  dPc(16) = 30^2 - 4;
S3exp = S;  S3exp(live, live) = S(live, live) + H(live, :) * diag(dPc) * H(live, :)';
ok3 = max(max(abs(S3 - S3exp))) < 1e-9;
%% (4) through the CPI: a clock jump and a position step
N = double(CST_spfParam.WINDOW_LENGTH);  T = CST_spfParam.CPI_THRESHOLD;
innC = zeros(mMax, N);  innP = zeros(mMax, N);  covH = zeros(mMax, mMax, N);  covF = zeros(mMax, mMax, N);
obsB = zeros(mMax, n, N);  numB = m * ones(N, 1);
SF = SPF_floorInnovCov(S, H, uint8(m), P, Pg, 16, 30.0);                   % the floors and the clock floor together
for i = 1:N
    noise = 0.1 * randn(m, 1);
    yC = zeros(mMax, 1);  yC(1:m) = noise;  yC(1:8) = yC(1:8) + 30;        % 30 m on every range row: a clock jump
    yP = zeros(mMax, 1);  yP(1:m) = noise + H(1:m, 1:3) * [20; 0; 0];      % a 20 m north step
    innC(:, i) = yC;  innP(:, i) = yP;  covH(:, :, i) = S;  covF(:, :, i) = SF;  obsB(:, :, i) = H;
end
qC = zeros(1, 3);  qCF = zeros(1, 3);  qP = zeros(1, 3);  qPF = zeros(1, 3);
for a = 1:3
    [~, qC(a)]  = SPF_cpiMonitor(innC, covH, obsB, numB, a, 0.0);
    [~, qCF(a)] = SPF_cpiMonitor(innC, covF, obsB, numB, a, 0.0);
    [~, qP(a)]  = SPF_cpiMonitor(innP, covH, obsB, numB, a, 0.0);
    [~, qPF(a)] = SPF_cpiMonitor(innP, covF, obsB, numB, a, 0.0);
end
fprintf('clock jump of 30 m: q/T_N %.3f on the host''s S, %.3f with the floors in; 20 m north step: %.3f and %.3f\n', ...
    max(qC) / T, max(qCF) / T, qP(1) / T, qPF(1) / T);
ok4 = max(qC) / T > 1 && max(qCF) / T < 0.05 && qP(1) / T > 1 && qPF(1) / T > 1;
%% (5) rows beyond numMeas untouched
S5 = SPF_floorInnovCov(S, H, uint8(8), P, Pg, 16, 30.0);
ok5 = max(max(abs(S5(1:8, 1:8) - SF(1:8, 1:8)))) < 1e-9 && isequal(S5(9:end, :), S(9:end, :)) && isequal(S5(:, 9:end), S(:, 9:end));
%% (6) the repository configuration
S6 = SPF_floorInnovCov(S, H, uint8(m), P, SPF_coastStartCov(P), CST_spfParam.CPI_CLOCK_STATES, CST_spfParam.CPI_CLOCK_SIGMA);
ok6 = isequal(S6, S);
if ok1 && ok2 && ok3 && ok4 && ok5 && ok6, fprintf('=== test_floorInnovCov PASS ===\n');
else, fprintf('ok1 %d ok2 %d ok3 %d ok4 %d ok5 %d ok6 %d\n=== test_floorInnovCov FAIL ===\n', ok1, ok2, ok3, ok4, ok5, ok6); end
end
