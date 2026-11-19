function test_rowScreen()
%% test_rowScreen.m
% SPF_rowScreen: the host's rows screened against their own fit, not against the filter. Rows: 8 GPS ranges, 6 GLONASS
% ranges, the same 14 range rates, a pressure row beyond numRows. (1) Noise only: nothing dropped, kept rows at the
% nominal R, the pressure row and the padding untouched, the fit near zero. (2) A 60 m single range row is dropped to
% sigma^2 times the drop factor, every other row nominal. (3) The filter wrong: 20 m north on every range row and
% 0.5 m/s north on every range rate, five sigma against the filter on the rows that see it, nothing dropped, every row
% nominal, the fit reading the shift and the velocity. (4) A 30 m jump on one constellation's ranges, a clock, nothing
% dropped. (5) Two outliers dropped one after the other, worst first, the smaller one on a 10 degree satellite whose
% leverage would hide it from a plain residual test. (6) A 1 m/s single range-rate row dropped while
% the consistent 0.5 m/s rows are kept. (7) No redundancy: four rows of one constellation, one at 60 m, nothing can be
% decided and nothing is dropped. (8) The override: the input comes back as it went. (9) Rows the host's setup had
% already inflated against the filter come back nominal when they agree with the fit: the screen supersedes the test.
% Run from the repo root (MATLAB, or Octave with tools/octave_shim).
randn('seed', 9);  rand('seed', 9);
n = CST_gnssHybrid.NO_STATES;  mMax = double(CST_spfParam.MAX_MEAS);
nG = 8;  nR = 6;  nSat = nG + nR;  numRows = 2 * nSat;
el = 0.15 + 1.2 * rand(nSat, 1);  az = 2 * pi * rand(nSat, 1);
u = [cos(el) .* cos(az), cos(el) .* sin(az), -sin(el)];
H = zeros(mMax, n);
for j = 1:nSat
    clk = 16 + 2 * double(j > nG);                                      % GPS clock 16 / drift 17, GLONASS 18 / 19
    H(j, 1:3) = -u(j, :);           H(j, clk) = 1;
    H(nSat + j, 4:6) = -u(j, :);    H(nSat + j, clk + 1) = 1;
end
H(numRows + 1, 3) = 1;                                                  % the pressure row, beyond numRows
sigR = 5.0;  sigD = 0.1;  K = 3.0;  DROP = 1e6;
rNom = zeros(mMax, 1);  rNom(1:nSat) = sigR^2;  rNom(nSat + 1:numRows) = sigD^2;  rNom(numRows + 1) = 25;
isRange = (1:mMax)' <= nSat;  isRate = (1:mMax)' > nSat & (1:mMax)' <= numRows;
noise = zeros(mMax, 1);  noise(1:nSat) = sigR * 0.5 * randn(nSat, 1);  noise(nSat + 1:numRows) = sigD * 0.5 * randn(nSat, 1);
clockCols = [16 18 21];  driftCols = [17 19 22];
nominalOut = @(r) all(abs(r(isRange) - sigR^2) < 1e-12) && all(abs(r(isRate) - sigD^2) < 1e-12);
%% (1) noise only
[r1, nd1, f1, fit1] = SPF_rowScreen(H, noise, rNom, uint8(numRows), false, clockCols, driftCols, [sigR sigD], K, DROP);
ok1 = isequal(nd1, [0 0]) && ~any(f1) && nominalOut(r1) && isequal(r1(numRows + 1:end), rNom(numRows + 1:end)) ...
    && max(abs(fit1(1, 1:2))) < 3 && max(abs(fit1(2, 1:2))) < 0.06;              % horizontal: the vertical is weak with two clocks
%% (2) a single 60 m range row
z2 = noise;  z2(3) = z2(3) + 60;
[r2, nd2, f2] = SPF_rowScreen(H, z2, rNom, uint8(numRows), false, clockCols, driftCols, [sigR sigD], K, DROP);
ok2 = isequal(nd2, [1 0]) && f2(3) && nnz(f2) == 1 && abs(r2(3) - sigR^2 * DROP) < 1e-6 ...
    && all(abs(r2([1 2 4:nSat]) - sigR^2) < 1e-12) && all(abs(r2(isRate) - sigD^2) < 1e-12);
%% (3) the filter wrong: a consistent shift on every row
z3 = noise + H(:, 1:3) * [20; 0; 0] + H(:, 4:6) * [0.5; 0; 0];
sigmaSeen = max(abs(z3(isRate))) / sigD;
[r3, nd3, f3, fit3] = SPF_rowScreen(H, z3, rNom, uint8(numRows), false, clockCols, driftCols, [sigR sigD], K, DROP);
ok3 = isequal(nd3, [0 0]) && ~any(f3) && nominalOut(r3) && abs(fit3(1, 1) - 20) < 3 && abs(fit3(2, 1) - 0.5) < 0.06 ...
    && sigmaSeen > 4;
fprintf('filter wrong by 20 m and 0.5 m/s: the rows reach %.1f sigma against it, none dropped, fit %.1f m and %.2f m/s north\n', ...
    sigmaSeen, fit3(1, 1), fit3(2, 1));
%% (4) a clock jump on one constellation
z4 = noise;  z4(nG + 1:nSat) = z4(nG + 1:nSat) + 30;
[r4, nd4, f4, fit4] = SPF_rowScreen(H, z4, rNom, uint8(numRows), false, clockCols, driftCols, [sigR sigD], K, DROP);
ok4 = isequal(nd4, [0 0]) && ~any(f4) && nominalOut(r4) && max(abs(fit4(1, :))) < 3;
%% (5) two outliers, worst first
z5 = noise;  z5(2) = z5(2) + 40;  z5(nG + 2) = z5(nG + 2) - 80;
[r5, nd5, f5] = SPF_rowScreen(H, z5, rNom, uint8(numRows), false, clockCols, driftCols, [sigR sigD], K, DROP);
ok5 = isequal(nd5, [2 0]) && f5(2) && f5(nG + 2) && nnz(f5) == 2 && abs(r5(2) - sigR^2 * DROP) < 1e-6 ...
    && abs(r5(nG + 2) - sigR^2 * DROP) < 1e-6;
%% (6) a single 1 m/s range-rate row among consistent 0.5 m/s rows
z6 = noise + H(:, 4:6) * [0.5; 0; 0];  z6(nSat + 5) = z6(nSat + 5) + 1.0;
[r6, nd6, f6, fit6] = SPF_rowScreen(H, z6, rNom, uint8(numRows), false, clockCols, driftCols, [sigR sigD], K, DROP);
ok6 = isequal(nd6, [0 1]) && f6(nSat + 5) && nnz(f6) == 1 && abs(r6(nSat + 5) - sigD^2 * DROP) < 1e-9 ...
    && abs(fit6(2, 1) - 0.5) < 0.06 && all(abs(r6(isRange) - sigR^2) < 1e-12);
%% (7) no redundancy: four rows of one constellation
H7 = zeros(mMax, n);  H7(1:4, :) = H(1:4, :);  z7 = zeros(mMax, 1);  z7(1:4) = noise(1:4);  z7(1) = z7(1) + 60;
[r7, nd7, f7] = SPF_rowScreen(H7, z7, rNom, uint8(4), false, clockCols, driftCols, [sigR sigD], K, DROP);
ok7 = isequal(nd7, [0 0]) && ~any(f7) && all(abs(r7(1:4) - sigR^2) < 1e-12) && isequal(r7(5:end), rNom(5:end));
%% (8) the override
rIn8 = rNom;  rIn8(1:numRows) = 225;
[r8, nd8, f8] = SPF_rowScreen(H, z2, rIn8, uint8(numRows), true, clockCols, driftCols, [sigR sigD], K, DROP);
ok8 = isequal(r8, rIn8) && isequal(nd8, [0 0]) && ~any(f8);
%% (9) the screen supersedes the host's test against the filter
rIn9 = rNom;  big = abs(z3(1:numRows)) > [sigR * ones(nSat, 1); sigD * ones(nSat, 1)];
rIn9(big) = z3(big).^2;                                                 % what the host's rule would have set
[r9, nd9, f9] = SPF_rowScreen(H, z3, rIn9, uint8(numRows), false, clockCols, driftCols, [sigR sigD], K, DROP);
ok9 = isequal(nd9, [0 0]) && ~any(f9) && nominalOut(r9) && nnz(big) > 10;
if ok1 && ok2 && ok3 && ok4 && ok5 && ok6 && ok7 && ok8 && ok9, fprintf('=== test_rowScreen PASS ===\n');
else, fprintf('ok1 %d ok2 %d ok3 %d ok4 %d ok5 %d ok6 %d ok7 %d ok8 %d ok9 %d\n=== test_rowScreen FAIL ===\n', ...
        ok1, ok2, ok3, ok4, ok5, ok6, ok7, ok8, ok9); end
end
