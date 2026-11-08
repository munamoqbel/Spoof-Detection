function [res] = imuModelCheck(rec, spec)
%IMUMODELCHECK  Allan deviation of an IMU model's injected error, read against the datasheet.
%   res = imuModelCheck(rec)        rec.dt: the sample interval (s); rec.gyro: [K x 3] delta angles per sample (rad);
%                                   rec.acc: [K x 3] delta velocities per sample (m/s); either may be left out. The
%                                   streams are the model's injected error: its output with the error terms on minus
%                                   its output with them off on the same input, or the model's own error terms, one
%                                   at a time or summed. Rates instead of increments (rad/s, m/s^2): rec.isRate true.
%   res = imuModelCheck(rec, spec)  spec: the datasheet: arw (deg/sqrt(h)), vrw (m/s/sqrt(h)), gyroBias (deg/h),
%                                   accBias (micro-g); default the vehicle IMU: 0.05, 0.0147, 0.48, 10.
%   Per axis: the constant part (this run's turn-on bias), the per-sample sigma against the one the datasheet's
%   white noise gives at this rate, the overlapping Allan deviation from the sample interval to a tenth of the
%   record, the random walk read at 1 s (white noise: the slope of -1/2) and the floor: the deviation at the
%   averaging time between 10 s and a twentieth of the record where the curve is flattest, which is the plateau of
%   a bias-stability term (a flicker floor reads 0.664 times the instability, a first-order Markov bias peaks at
%   0.62 times its sigma). The floor counts as reached when the local slope there is within 0.3 of flat and the
%   deviation stands above what the white noise read at 1 s leaves at that averaging time by three standard
%   errors; otherwise the white noise still covers the bias terms at this record length, the floor printed is an
%   upper bound, and a longer record or the bias term on its own reads them. Each reading is
%   also printed as a ratio to the datasheet: a model that emits the datasheet reads ratios near 1. For the vehicle
%   IMU the expected readings are 3.0 deg/h and 25 micro-g at 1 s and a floor near 0.32 deg/h and 6.6 micro-g.
%   Result per channel (res.gyro, res.acc): bias and biasSpec [1 x 3, deg/h or micro-g], sigmaSample and
%   sigmaSampleSpec (per sample), tau and allan [nTau x 3], rwAt1s, rwSpec and rwRatio, floor, floorTau,
%   floorSlope, floorReached (per axis), floorSpec (0.664 times the datasheet instability) and floorRatio.
%   MATLAB desktop tool, not codegen. Runs under Octave.
if nargin < 2 || isempty(spec), spec = struct('arw', 0.05, 'vrw', 0.0147, 'gyroBias', 0.48, 'accBias', 10); end
isRate = isfield(rec, 'isRate') && rec.isRate;
res = struct('dt', rec.dt, 'spec', spec);
if isfield(rec, 'gyro') && ~isempty(rec.gyro), res.gyro = channel(rec.gyro, rec.dt, isRate, true, spec); end
if isfield(rec, 'acc') && ~isempty(rec.acc), res.acc = channel(rec.acc, rec.dt, isRate, false, spec); end
end

function [c] = channel(x, dt, isRate, isGyro, spec)
FLAT = 0.3;                                                                      % |slope| below this: a floor
K = size(x, 1);
if isRate, inc = x * dt; else, inc = x; end                                     % increments per sample
if isGyro
    dens = spec.arw * pi / 180 / 60;             % rad/sqrt(s): the datasheet white noise
    biasSpec = spec.gyroBias * pi / 180 / 3600;  % rad/s
    toRw = 180 / pi * 60;  toRate = 180 / pi * 3600;
    unitRw = 'deg/sqrt(h)';  unitRate = 'deg/h';  rwName = 'ARW';  name = 'gyro';  unitInc = 'rad';
else
    dens = spec.vrw / 60;                        % m/s/sqrt(s)
    biasSpec = spec.accBias * 9.80665e-6;        % m/s^2
    toRw = 60;  toRate = 1 / 9.80665e-6;
    unitRw = 'm/s/sqrt(h)';  unitRate = 'micro-g';  rwName = 'VRW';  name = 'accelerometer';  unitInc = 'm/s';
end
bias = mean(inc, 1) / dt;                                                        % the constant part, rate units
sigmaSample = std(inc, 0, 1);
sigmaSpec = dens * sqrt(dt);
mList = unique(round(logspace(0, log10(max(floor(K / 10), 1)), 40)));
mList = mList(mList >= 1 & 2 * mList < K);
tau = mList(:) * dt;
allan = allanDev(inc, dt, mList);
m1 = round(1 / dt);
rw = NaN(1, 3);
if m1 >= 1 && 2 * m1 < K, rw = allanDev(inc, dt, m1); end
% the floor: the flattest point of the log-log curve between 10 s and a twentieth of the record
fl = NaN(1, 3);  flTau = NaN(1, 3);  flSlope = NaN(1, 3);  reached = false(1, 3);
centre = find(tau >= 10 & tau <= K * dt / 20);
if numel(centre) >= 3
    for a = 1:3
        slope = NaN(numel(centre), 1);
        for j = 1:numel(centre)
            w = tau >= tau(centre(j)) / sqrt(10) & tau <= tau(centre(j)) * sqrt(10) & isfinite(allan(:, a)) & allan(:, a) > 0;
            if nnz(w) >= 3
                p = polyfit(log10(tau(w)), log10(allan(w, a)), 1);
                slope(j) = p(1);
            end
        end
        [~, j] = min(abs(slope));
        if isfinite(slope(j))
            fl(a) = allan(centre(j), a);  flTau(a) = tau(centre(j));  flSlope(a) = slope(j);
            whiteLeft = rw(a) / sqrt(flTau(a));                                  % what the 1 s white noise leaves here
            relErr = 1 / sqrt(2 * (K / mList(centre(j)) - 1));                   % standard error of the deviation here
            reached(a) = abs(slope(j)) < FLAT && fl(a) > whiteLeft * (1 + 3 * relErr);
        end
    end
end
floorSpec = 0.664 * biasSpec;
c = struct('name', name, 'bias', bias * toRate, 'biasUnit', unitRate, 'biasSpec', biasSpec * toRate, ...
    'sigmaSample', sigmaSample, 'sigmaSampleSpec', sigmaSpec, 'tau', tau, 'allan', allan * toRate, 'allanUnit', unitRate, ...
    'rwAt1s', rw * toRw, 'rwUnit', unitRw, 'rwSpec', dens * toRw, 'rwRatio', rw / dens, ...
    'floor', fl * toRate, 'floorTau', flTau, 'floorSlope', flSlope, 'floorReached', reached, ...
    'floorSpec', floorSpec * toRate, 'floorRatio', fl / floorSpec);
fprintf('%s: %d samples at %g Hz, %.0f s\n', name, K, 1 / dt, K * dt);
fprintf('   constant part, this run''s turn-on bias [%s]: %s (the datasheet repeatability is one sigma per run)\n', ...
    unitRate, v3(c.bias, '%.3g'));
fprintf('   per-sample sigma [%s]: %s; the datasheet white noise gives %.3g at this rate: ratios %s\n', ...
    unitInc, v3(sigmaSample, '%.3g'), sigmaSpec, v3(sigmaSample / sigmaSpec, '%.2f'));
fprintf('   Allan deviation [%s]:\n   %10s %12s %12s %12s\n', unitRate, 'tau [s]', 'X', 'Y', 'Z');
pick = unique([1:3:numel(tau), numel(tau)]);
for i = pick
    fprintf('   %10.3g %12.4g %12.4g %12.4g\n', tau(i), c.allan(i, 1), c.allan(i, 2), c.allan(i, 3));
end
fprintf('   %s at 1 s [%s]: %s; datasheet %.4g: ratios %s\n', rwName, unitRw, v3(c.rwAt1s, '%.4g'), c.rwSpec, v3(c.rwRatio, '%.2f'));
fprintf('   floor, the flattest point beyond 10 s [%s]: %s at tau %s s, local slope %s; a flicker floor at the datasheet\n', ...
    unitRate, v3(c.floor, '%.3g'), v3(flTau, '%.3g'), v3(flSlope, '%+.2f'));
fprintf('   instability %.3g would read %.3g, a Markov bias of that sigma %.3g: ratios to the flicker floor %s\n', ...
    c.biasSpec, c.floorSpec, 0.62 * c.biasSpec, v3(c.floorRatio, '%.2f'));
if all(reached)
    fprintf('   VERDICT: white noise %s the datasheet; the floor %s the datasheet''s\n', times(c.rwRatio), times(c.floorRatio));
else
    fprintf('   VERDICT: white noise %s the datasheet; no floor reached on axis %s within the record: the white noise\n', ...
        times(c.rwRatio), mat2str(find(~reached)));
    fprintf('           still covers the bias terms there, the floor above is an upper bound; a longer record or the bias\n');
    fprintf('           term on its own reads them\n');
end
end

function [ad] = allanDev(inc, dt, mList)
% the overlapping Allan deviation of the increment streams [K x 3] at the cluster sizes mList, rate units
K = size(inc, 1);
theta = [zeros(1, size(inc, 2)); cumsum(inc, 1)];
ad = NaN(numel(mList), size(inc, 2));
for i = 1:numel(mList)
    m = mList(i);
    if 2 * m < K
        d = theta(1 + 2 * m:end, :) - 2 * theta(1 + m:end - m, :) + theta(1:end - 2 * m, :);
        ad(i, :) = sqrt(sum(d .^ 2, 1) / (2 * (m * dt) ^ 2 * (K + 1 - 2 * m)));
    end
end
end

function [s] = v3(v, fmt)
s = strtrim(sprintf([fmt ' '], v));
end

function [s] = times(r)
r = max(r);
if r > 2, s = sprintf('%.1f times', r);
elseif r < 0.5, s = sprintf('%.2f of', r);
else, s = 'near'; end
end
