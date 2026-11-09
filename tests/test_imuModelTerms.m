function test_imuModelTerms()
%% test_imuModelTerms.m
% tools/imuModelTerms on a synthetic record in the host's layout: the model's terms per sample as rates (white noise
% at the datasheet densities, first-order Markov bias stability at the datasheet sigma, rate random walks, constant
% turn-on biases, scale factors in ppm, a g-sensitivity matrix, the sums), the sample times, and the model's input
% and output beside them, the output composed as input (1 + SF) + drift f + error, times the interval. The tool must
% read the noise densities within 10 % with no floor, the stability sigmas within 25 %, the walks' slope near +1/2,
% the constants exactly, the sums as the sum of the terms, and name the right composition at lag 0, with the input
% taken at the sample or integrated over the interval; the struct-array layout (one element per sample) must read
% the same as the struct of arrays. Run from anywhere:
% octave --no-gui --quiet --eval test_imuModelTerms
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'tools'));
rng(2);
G = 9.80665;  R2D = 180 / pi;
dt = 0.009;  K = 100000;  recTime = dt * (1:K)';
arw = 0.06 / R2D / 60;  vrw = 0.0147 / 60;                      % rad/sqrt(s), m/s/sqrt(s)
gS = 0.48 / R2D / 3600;  aS = 10 * G * 1e-6;                     % rad/s, m/s^2
tauC = 30;  a = exp(-dt / tauC);
gm = @(s) filter(1, [1 -a], s * sqrt(1 - a^2) * randn(K, 3));
qG = 0.02 / R2D / 3600;  qA = 1 * G * 1e-6;                      % rate random walk densities per sqrt(s)
rec = struct();
rec.rndNoiseGyro = arw / sqrt(dt) * randn(K, 3);  rec.rndNoiseAcc = vrw / sqrt(dt) * randn(K, 3);
rec.Bias_Stability_Gyro = gm(gS);  rec.Bias_Stability_Acc = gm(aS);
rec.rndWalkGyro = cumsum(qG * sqrt(dt) * randn(K, 3), 1);  rec.rndWalkAcc = cumsum(qA * sqrt(dt) * randn(K, 3), 1);
biasG = [1.38 -0.67 1.52] / R2D / 3600;  biasA = [255 96 -387] * G * 1e-6;
rec.BiasGyro = repmat(biasG, K, 1);  rec.BiasAcc = repmat(biasA, K, 1);
sfG = [-61.4574 -205.6508 -50.1267];  sfA = [43.8806 -301.6522 -96.7046];
rec.SFGyro = repmat(sfG, K, 1);  rec.SFAcc = repmat(sfA, K, 1);
drift = 1e-8 * [-0.2 1.39 -1.19; -0.12 1.4 0.71; 1.47 0.66 1.61];
rec.GyroAccDriftMatrix = drift;
rec.errorGyro = rec.rndNoiseGyro + rec.Bias_Stability_Gyro + rec.rndWalkGyro + rec.BiasGyro;
rec.errorAcc = rec.rndNoiseAcc + rec.Bias_Stability_Acc + rec.rndWalkAcc + rec.BiasAcc;
io.fibb = repmat([0 0 -G], K, 1) + 0.05 * randn(K, 3);          % the reference's own noise rides on the input
io.wibb = repmat([0 5.3e-5 -5.0e-5], K, 1) + 1e-3 * randn(K, 3);
io.dAng = (io.wibb .* (1 + sfG * 1e-6) + (drift * io.fibb')' + rec.errorGyro) * dt;
io.dVel = (io.fibb .* (1 + sfA * 1e-6) + rec.errorAcc) * dt;
fprintf('--- the struct of arrays ---\n');
r = imuModelTerms(rec, recTime, io);
okNoise = all(abs(r.noise.gyro.rwRatio - 1) < 0.1) && all(abs(r.noise.acc.rwRatio - 1) < 0.1) ...
    && ~any(r.noise.gyro.floorReached) && ~any(r.noise.acc.floorReached);
okStab = all(abs(r.stability.gyroSigma / 0.48 - 1) < 0.25) && all(abs(r.stability.accSigma / 10 - 1) < 0.25);
[slG, ~] = slopeOf(r.walk.gyro);  [slA, ~] = slopeOf(r.walk.acc);
okWalk = all(slG > 0.3 & slG < 0.7) && all(slA > 0.3 & slA < 0.7);
rmsOf = @(x) sqrt(mean(x .^ 2, 1));
okConst = max(abs(r.constants.biasGyro - biasG)) < 1e-12 && max(abs(r.constants.biasAcc - biasA)) < 1e-12 ...
    && isequal(r.constants.sfGyro, sfG) && isequal(r.constants.sfAcc, sfA) && isequal(r.constants.drift, drift);
okSum = r.sumGyroResidual < 1e-12 && r.sumAccResidual < 1e-12;
okOut = ~isempty(strfind(r.outGyro.best, 'at the sample (1 + SF) + drift')) && ~isempty(strfind(r.outGyro.best, 'lag +0')) ...
    && ~isempty(strfind(r.outAcc.best, 'at the sample (1 + SF)')) && ~isempty(strfind(r.outAcc.best, 'lag +0')) ...
    && all(squeeze(r.outGyro.rms(2, 3, :))' < 1e-8 * rmsOf(io.dAng)) && all(squeeze(r.outAcc.rms(2, 3, :))' < 1e-8 * rmsOf(io.dVel));
% the same model composing its output from the input integrated over the interval
trap = @(x) [x(1, :); (x(1:end - 1, :) + x(2:end, :)) / 2];
ioT = io;
ioT.dAng = (trap(io.wibb) .* (1 + sfG * 1e-6) + (drift * trap(io.fibb)')' + rec.errorGyro) * dt;
ioT.dVel = (trap(io.fibb) .* (1 + sfA * 1e-6) + rec.errorAcc) * dt;
fprintf('--- the output composed over the interval ---\n');
evalc('rT = imuModelTerms(rec, recTime, ioT);');
okTrap = ~isempty(strfind(rT.outGyro.best, 'over the interval (1 + SF) + drift')) && ~isempty(strfind(rT.outGyro.best, 'lag +0')) ...
    && ~isempty(strfind(rT.outAcc.best, 'over the interval (1 + SF)')) && ~isempty(strfind(rT.outAcc.best, 'lag +0')) ...
    && all(squeeze(rT.outGyro.rms(4, 3, :))' < 1e-8 * rmsOf(ioT.dAng)) && all(squeeze(rT.outAcc.rms(4, 3, :))' < 1e-8 * rmsOf(ioT.dVel));
fprintf('   best: %s | %s\n', rT.outGyro.best, rT.outAcc.best);
% the struct-array layout, one element per sample, on the first samples
n = 3000;
names = fieldnames(rec);  args = {};
for i = 1:numel(names)
    x = rec.(names{i});
    if strcmp(names{i}, 'GyroAccDriftMatrix'), args(end+1:end+2) = {names{i}, repmat({drift}, n, 1)};
    else, args(end+1:end+2) = {names{i}, num2cell(x(1:n, :), 2)}; end
end
recArr = struct(args{:});
fprintf('--- the struct array, %d samples ---\n', n);
ioN = struct('fibb', io.fibb(1:n, :), 'wibb', io.wibb(1:n, :), 'dAng', io.dAng(1:n, :), 'dVel', io.dVel(1:n, :));
s = imuModelTerms(recArr, recTime(1:n), ioN);
okArr = numel(recArr) == n && isequal(s.constants.drift, drift) && isequal(s.constants.sfAcc, sfA) ...
    && max(abs(s.constants.biasGyro - biasG)) < 1e-12 && ~isempty(strfind(s.outGyro.best, 'lag +0')) ...
    && abs(s.noise.gyro.sigmaSample(1) - std(rec.rndNoiseGyro(1:n, 1)) * dt) < 1e-12;
fprintf('noise densities, no floor: %s | stability sigmas: %s | walks rise: %s | constants: %s | sums: %s | composition at lag 0: %s | over the interval: %s | struct array reads the same: %s\n', ...
    pf(okNoise), pf(okStab), pf(okWalk), pf(okConst), pf(okSum), pf(okOut), pf(okTrap), pf(okArr));
if okNoise && okStab && okWalk && okConst && okSum && okOut && okTrap && okArr, fprintf('=== test_imuModelTerms PASS ===\n');
else, fprintf('=== test_imuModelTerms FAIL ===\n'); end
end

function [sl, q] = slopeOf(c)
% the slope of the Allan deviation from 1 s up, as the tool reads it
sl = NaN(1, 3);  q = NaN(1, 3);
w = c.tau >= min(1, c.tau(end) / 10);
for a = 1:3
    p = polyfit(log10(c.tau(w)), log10(c.allan(w, a)), 1);  sl(a) = p(1);
end
end

function s = pf(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end
