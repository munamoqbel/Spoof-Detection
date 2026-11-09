function [res] = imuModelTerms(recErrorIMU, recTime, io, spec, verbose)
%IMUMODELTERMS  Read an IMU model's recorded error terms, one at a time, against the datasheet.
%   res = imuModelTerms(recErrorIMU, recTime)       recErrorIMU: the per-sample record of the model's terms, a struct
%                                                   array (one element per sample) or a struct of [K x 3] arrays, as
%                                                   rates: rndNoiseAcc (m/s^2), rndNoiseGyro (rad/s), Bias_Stability_Acc,
%                                                   Bias_Stability_Gyro, rndWalkAcc, rndWalkGyro, BiasAcc, BiasGyro
%                                                   (constants per run), SFAcc, SFGyro (ppm), GyroAccDriftMatrix
%                                                   (rad/s per m/s^2, 3 x 3), errorAcc, errorGyro (the sums the model
%                                                   adds). A term that is missing is skipped. recTime: the sample times
%                                                   (s); the interval is their median difference.
%   res = imuModelTerms(recErrorIMU, recTime, io)   io: the model's input and output beside the terms, for the sum check:
%                                                   io.fibb, io.wibb [K x 3] (the input specific force m/s^2 and angular
%                                                   rate rad/s) and io.dVel, io.dAng [K x 3] (the output delta velocity
%                                                   m/s and delta angle rad). [] to skip.
%   res = imuModelTerms(..., spec, verbose)         spec: the datasheet as for imuModelCheck; verbose true prints the
%                                                   full imuModelCheck report of every term (default: the summary only).
%   What each term is read for: the noise terms for the random walk at 1 s against the ARW and VRW, and that they
%   carry no floor; the bias-stability terms for their standard deviation over the record against the datasheet
%   stability, and for their Allan floor; the random-walk terms for the slope of their Allan deviation (+1/2) and the
%   density it implies; the constant biases against the repeatability (one sigma per run); the scale factors and the
%   g-sensitivity matrix as constants, for the record; the sums against the sum of the terms. With io, the output is
%   compared with the input and the sums sample by sample, with and without the scale factor and the g-sensitivity
%   on the input, at lags of -2 to 2 samples: the smallest residual says how the model composes its output and
%   whether the recorded terms are all that it adds. MATLAB desktop tool, not codegen. Runs under Octave.
if nargin < 3, io = []; end
if nargin < 4 || isempty(spec), spec = struct('arw', 0.05, 'vrw', 0.0147, 'gyroBias', 0.48, 'accBias', 10); end
if nargin < 5 || isempty(verbose), verbose = false; end
G = 9.80665;  R2D = 180 / pi;
t = recTime(:);
d = diff(t);  dt = median(d);
K = numel(t);
fprintf('imuModelTerms: %d samples, %.1f s, interval %.4g s (%.4g Hz); the interval spreads from %.4g to %.4g s\n', ...
    K, t(end) - t(1), dt, 1 / dt, min(d), max(d));
fprintf('   datasheet: ARW %g deg/sqrt(h), VRW %g m/s/sqrt(h), gyro bias stability %g deg/h, accelerometer %g micro-g\n', ...
    spec.arw, spec.vrw, spec.gyroBias, spec.accBias);
res = struct('dt', dt, 'K', K, 'spec', spec);
term = @(name) termOf(recErrorIMU, name, K);

%% 1. the white noise terms: the random walks at 1 s, no floor
[gN, aN] = deal(term('rndNoiseGyro'), term('rndNoiseAcc'));
if ~isempty(gN) || ~isempty(aN)
    c = check(struct('dt', dt, 'gyro', gN, 'acc', aN, 'isRate', true), spec, verbose, 'the noise terms');
    res.noise = c;
    fprintf('1. noise terms\n');
    if isfield(c, 'gyro')
        fprintf('   gyro  ARW at 1 s %s deg/sqrt(h), datasheet %g: ratios %s; floor reached on any axis: %s\n', ...
            v3(c.gyro.rwAt1s, '%.4g'), spec.arw, v3(c.gyro.rwRatio, '%.2f'), yesNo(any(c.gyro.floorReached)));
    end
    if isfield(c, 'acc')
        fprintf('   accel VRW at 1 s %s m/s/sqrt(h), datasheet %g: ratios %s; floor reached on any axis: %s\n', ...
            v3(c.acc.rwAt1s, '%.4g'), spec.vrw, v3(c.acc.rwRatio, '%.2f'), yesNo(any(c.acc.floorReached)));
    end
    fprintf('   (a noise term reads its density at 1 s and no floor; a floor here means the term is coloured)\n');
end

%% 2. the bias stability terms: their sigma over the record, their Allan floor
[gS, aS] = deal(term('Bias_Stability_Gyro'), term('Bias_Stability_Acc'));
if ~isempty(gS) || ~isempty(aS)
    c = check(struct('dt', dt, 'gyro', gS, 'acc', aS, 'isRate', true), spec, verbose, 'the bias stability terms');
    res.stability = c;
    fprintf('2. bias stability terms\n');
    if ~isempty(gS)
        sg = std(gS, 0, 1) * R2D * 3600;
        fprintf('   gyro  sigma over the record %s deg/h, datasheet %g: ratios %s; Allan floor %s deg/h at tau %s s, reached: %s\n', ...
            v3(sg, '%.3g'), spec.gyroBias, v3(sg / spec.gyroBias, '%.2f'), v3(c.gyro.floor, '%.3g'), v3(c.gyro.floorTau, '%.3g'), ...
            yesNo(all(c.gyro.floorReached)));
        res.stability.gyroSigma = sg;
    end
    if ~isempty(aS)
        sa = std(aS, 0, 1) / (G * 1e-6);
        fprintf('   accel sigma over the record %s micro-g, datasheet %g: ratios %s; Allan floor %s micro-g at tau %s s, reached: %s\n', ...
            v3(sa, '%.3g'), spec.accBias, v3(sa / spec.accBias, '%.2f'), v3(c.acc.floor, '%.3g'), v3(c.acc.floorTau, '%.3g'), ...
            yesNo(all(c.acc.floorReached)));
        res.stability.accSigma = sa;
    end
    fprintf('   (the sigma is the direct reading of a Markov or flicker bias; the floor reads 0.62 to 0.66 of it)\n');
end

%% 3. the random walk terms: the Allan slope and the density
[gW, aW] = deal(term('rndWalkGyro'), term('rndWalkAcc'));
if ~isempty(gW) || ~isempty(aW)
    c = check(struct('dt', dt, 'gyro', gW, 'acc', aW, 'isRate', true), spec, verbose, 'the random walk terms');
    res.walk = c;
    fprintf('3. random walk terms (no datasheet line; a rate random walk rises with the slope +1/2)\n');
    if ~isempty(gW)
        [sl, q] = walkOf(c.gyro.tau, c.gyro.allan / (R2D * 3600));
        fprintf('   gyro  Allan slope from 1 s up %s; density %s deg/h per sqrt(h); the term spans %s deg/h over the record\n', ...
            v3(sl, '%+.2f'), v3(q * R2D * 3600 * 60, '%.3g'), v3((max(gW, [], 1) - min(gW, [], 1)) * R2D * 3600, '%.3g'));
    end
    if ~isempty(aW)
        [sl, q] = walkOf(c.acc.tau, c.acc.allan * G * 1e-6);
        fprintf('   accel Allan slope from 1 s up %s; density %s micro-g per sqrt(h); the term spans %s micro-g over the record\n', ...
            v3(sl, '%+.2f'), v3(q / (G * 1e-6) * 60, '%.3g'), v3((max(aW, [], 1) - min(aW, [], 1)) / (G * 1e-6), '%.3g'));
    end
end

%% 4. the constants of this run
fprintf('4. constants of this run (the datasheet repeatability is one sigma per run)\n');
bG = constOf(term('BiasGyro'));  bA = constOf(term('BiasAcc'));
if ~isempty(bG), fprintf('   gyro turn-on bias  %s deg/h (repeatability 0.5 deg/h one sigma on the vehicle IMU)\n', v3(bG * R2D * 3600, '%.3g')); end
if ~isempty(bA), fprintf('   accel turn-on bias %s micro-g (repeatability 300 micro-g one sigma on the vehicle IMU)\n', v3(bA / (G * 1e-6), '%.3g')); end
sG = constOf(term('SFGyro'));  sA = constOf(term('SFAcc'));
if ~isempty(sG), fprintf('   gyro scale factor  %s ppm\n', v3(sG, '%.4g')); end
if ~isempty(sA), fprintf('   accel scale factor %s ppm\n', v3(sA, '%.4g')); end
drift = driftOf(recErrorIMU);
if ~isempty(drift)
    fprintf('   g-sensitivity matrix [deg/h per g]:\n');
    for i = 1:3, fprintf('      %s\n', v3(drift(i, :) * G * R2D * 3600, '%10.4g')); end
end
res.constants = struct('biasGyro', bG, 'biasAcc', bA, 'sfGyro', sG, 'sfAcc', sA, 'drift', drift);

%% 5. the sums against the terms
[gE, aE] = deal(term('errorGyro'), term('errorAcc'));
fprintf('5. the sums the model adds\n');
if ~isempty(gE) && ~isempty(gN)
    parts = gN + opt(gS, K) + opt(gW, K) + opt(term('BiasGyro'), K);
    r = gE - parts;
    fprintf('   errorGyro minus (noise + stability + walk + bias): rms %s rad/s against the sum''s rms %s\n', ...
        v3(rms(r), '%.3g'), v3(rms(gE), '%.3g'));
    res.sumGyroResidual = rms(r);
end
if ~isempty(aE) && ~isempty(aN)
    parts = aN + opt(aS, K) + opt(aW, K) + opt(term('BiasAcc'), K);
    r = aE - parts;
    fprintf('   errorAcc  minus (noise + stability + walk + bias): rms %s m/s^2 against the sum''s rms %s\n', ...
        v3(rms(r), '%.3g'), v3(rms(aE), '%.3g'));
    res.sumAccResidual = rms(r);
end
if ~isempty(gE) || ~isempty(aE)
    c = check(struct('dt', dt, 'gyro', gE, 'acc', aE, 'isRate', true), spec, verbose, 'the sums');
    res.sum = c;
    if ~isempty(gE), fprintf('   on errorGyro: ARW at 1 s ratios %s, floor %s deg/h reached: %s\n', v3(c.gyro.rwRatio, '%.2f'), v3(c.gyro.floor, '%.3g'), yesNo(all(c.gyro.floorReached))); end
    if ~isempty(aE), fprintf('   on errorAcc:  VRW at 1 s ratios %s, floor %s micro-g reached: %s\n', v3(c.acc.rwRatio, '%.2f'), v3(c.acc.floor, '%.3g'), yesNo(all(c.acc.floorReached))); end
end

%% 6. the output against the input and the sums
if ~isempty(io)
    fprintf('6. the output against the input plus the sums, lags of -2 to 2 samples (the smallest residual tells how the model composes its output)\n');
    if isfield(io, 'dAng') && ~isempty(gE)
        base = io.wibb + gE;
        withSf = io.wibb .* (1 + sG(:)' * 1e-6) + gE;
        if ~isempty(drift), withSf = withSf + (drift * io.fibb')'; end
        res.outGyro = compose(io.dAng, dt, {base, withSf}, {'input + errorGyro', 'input (1 + SF) + drift * f + errorGyro'}, 'delta angle [rad]');
    end
    if isfield(io, 'dVel') && ~isempty(aE)
        base = io.fibb + aE;
        withSf = io.fibb .* (1 + sA(:)' * 1e-6) + aE;
        res.outAcc = compose(io.dVel, dt, {base, withSf}, {'input + errorAcc', 'input (1 + SF) + errorAcc'}, 'delta velocity [m/s]');
    end
end
end

function [x] = termOf(rec, name, K)
% the term as [K x n], from a struct array (one element per sample) or a struct of arrays; [] when absent
x = [];
if ~isfield(rec, name), return; end
if numel(rec) > 1
    x = cat(1, rec.(name));
    if size(x, 1) ~= K && size(x, 1) == 3 * K, x = reshape(x', [], K)'; end   % 3 x 3 per element stacked
else
    x = rec.(name);
    if size(x, 1) ~= K && size(x, 2) == K, x = x'; end
end
end

function [m] = driftOf(rec)
m = [];
if ~isfield(rec, 'GyroAccDriftMatrix'), return; end
if numel(rec) > 1, m = rec(1).GyroAccDriftMatrix; else, m = rec.GyroAccDriftMatrix; end
if numel(m) ~= 9, m = []; else, m = reshape(m, 3, 3); end
end

function [c] = constOf(x)
% the constant of a term that should not move: its median (exact on a constant); a warning when it moves
c = [];
if isempty(x), return; end
c = median(x, 1);
sp = max(x, [], 1) - min(x, [], 1);
if any(sp > 1e-9 * max(abs(c), 1e-30))
    fprintf('   NOTE: a constant term moves over the record by %s\n', v3(sp, '%.3g'));
end
end

function [y] = opt(x, K)
if isempty(x), y = zeros(K, 3); else, y = x; end
end

function [c] = check(rec, spec, verbose, what)
if verbose
    fprintf('--- imuModelCheck on %s ---\n', what);
    c = imuModelCheck(rec, spec);
else
    evalc('c = imuModelCheck(rec, spec);');                                      % the report captured, not printed
end
end

function [sl, q] = walkOf(tau, allan)
% the slope of the Allan deviation from 1 s up (a walk term alone rises as +1/2 throughout) and the rate random
% walk density q (ad = q sqrt(tau / 3)) at the end of the record
sl = NaN(1, 3);  q = NaN(1, 3);
w = tau >= min(1, tau(end) / 10);
for a = 1:3
    ok = w & isfinite(allan(:, a)) & allan(:, a) > 0;
    if nnz(ok) >= 3
        p = polyfit(log10(tau(ok)), log10(allan(ok, a)), 1);  sl(a) = p(1);
    end
    q(a) = allan(end, a) * sqrt(3 / tau(end));
end
end

function [r] = compose(out, dt, models, names, what)
% the output against each model of it, at lags of -2 to 2 samples: rms residual per axis, the first and last two
% samples left out
K = size(out, 1);  lags = -2:2;
r = struct('names', {names}, 'lags', lags, 'rms', NaN(numel(models), numel(lags), 3), 'best', '');
best = Inf;
for i = 1:numel(models)
    for j = 1:numel(lags)
        L = lags(j);  k = 3:K - 2;
        res = out(k, :) - models{i}(k + L, :) * dt;
        r.rms(i, j, :) = rms(res);
        if sum(r.rms(i, j, :)) < best, best = sum(r.rms(i, j, :));  r.best = sprintf('%s at lag %+d', names{i}, L);  r.bestMean = mean(res, 1); end
    end
end
fprintf('   %s: rms of output minus model * dt per axis\n', what);
for i = 1:numel(models)
    for j = 1:numel(lags)
        fprintf('      %-45s lag %+d: %s\n', names{i}, lags(j), v3(squeeze(r.rms(i, j, :))', '%.3g'));
    end
end
fprintf('      output rms %s; best: %s, mean residual %s\n', v3(rms(out), '%.3g'), r.best, v3(r.bestMean, '%.3g'));
end

function [v] = rms(x)
v = sqrt(mean(x .^ 2, 1));
end

function [s] = v3(v, fmt)
s = strtrim(sprintf([fmt ' '], v));
end

function [s] = yesNo(b)
if b, s = 'yes'; else, s = 'no'; end
end
