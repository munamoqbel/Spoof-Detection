function test_imuModelCheck()
%% test_imuModelCheck.m
% tools/imuModelCheck on synthetic error streams with known terms, 2000 s at 100 Hz, as an IMU model outputs them
% one at a time or summed: white noise at the datasheet densities, a first-order Gauss-Markov bias at the datasheet
% stability with a 30 s correlation time, a constant turn-on bias at the repeatability. On the sum the random walks
% must read within 10 % at 1 s and the constant within a factor 1.5; on the bias term alone the floor must be
% reached and read within a factor 1.5 of the Markov plateau, 0.62 sigma; white noise alone must read no floor; rates in
% place of increments must give the same readings. Run from anywhere: octave --no-gui --quiet --eval test_imuModelCheck
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'tools'));
rng(1);
dt = 0.01;  K = 200000;
arw = 0.05 * pi / 180 / 60;            % rad/sqrt(s)
vrw = 0.0147 / 60;                     % m/s/sqrt(s)
gBias = 0.48 * pi / 180 / 3600;        % rad/s, the stability
aBias = 10 * 9.80665e-6;               % m/s^2
gRep = 0.5 * pi / 180 / 3600;          % rad/s, the turn-on bias of this run
aRep = 300 * 9.80665e-6;
tauC = 30;  a = exp(-dt / tauC);
gm = @(s) filter(1, [1 -a], s * sqrt(1 - a^2) * randn(K, 3));
gWhite = arw / sqrt(dt) * randn(K, 3);  gMarkov = gm(gBias);
aWhite = vrw / sqrt(dt) * randn(K, 3);  aMarkov = gm(aBias);
fprintf('--- the terms summed ---\n');
r = imuModelCheck(struct('dt', dt, 'gyro', (gWhite + gMarkov + gRep) * dt, 'acc', (aWhite + aMarkov + aRep) * dt));
okRw = all(abs(r.gyro.rwRatio - 1) < 0.1) && all(abs(r.acc.rwRatio - 1) < 0.1);
okBias = all(r.gyro.bias / 0.5 > 0.5 & r.gyro.bias / 0.5 < 1.5) && all(r.acc.bias / 300 > 0.5 & r.acc.bias / 300 < 1.5);
fprintf('--- the bias term alone ---\n');
b = imuModelCheck(struct('dt', dt, 'gyro', gMarkov * dt, 'acc', aMarkov * dt));
plateau = @(fl, s) all(fl / (0.62 * s) > 1 / 1.5 & fl / (0.62 * s) < 1.5);
okFloor = all(b.gyro.floorReached) && all(b.acc.floorReached) && plateau(b.gyro.floor, 0.48) && plateau(b.acc.floor, 10);
fprintf('--- white noise alone ---\n');
w = imuModelCheck(struct('dt', dt, 'gyro', gWhite * dt));
okWhite = all(abs(w.gyro.rwRatio - 1) < 0.1) && ~any(w.gyro.floorReached) && all(abs(w.gyro.bias) < 0.25);   % the mean of 2000 s of white noise: 0.067 deg/h one sigma
fprintf('--- rates in place of increments ---\n');
q = imuModelCheck(struct('dt', dt, 'gyro', gWhite + gMarkov + gRep, 'isRate', true));
okRate = max(abs(q.gyro.allan(:) - r.gyro.allan(:))) < 1e-9 && max(abs(q.gyro.rwAt1s - r.gyro.rwAt1s)) < 1e-12;
fprintf('random walks within 10 %%: %s | turn-on bias within 1.5: %s | the bias term reads its plateau: %s | white noise alone reads no floor: %s | rates equal increments: %s\n', ...
    pf(okRw), pf(okBias), pf(okFloor), pf(okWhite), pf(okRate));
if okRw && okBias && okFloor && okWhite && okRate, fprintf('=== test_imuModelCheck PASS ===\n');
else, fprintf('=== test_imuModelCheck FAIL ===\n'); end
end

function s = pf(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end
