function test_spoofInject()
%% test_spoofInject.m
% tools/spoofInject and tools/spoofScenario on one constellation of a stationary receiver: six used slots of
% fifteen, satellites 2.0e7 m away moving at 3 km/s, a 1500 m receiver clock and a 0.3 m/s drift. A step, a
% ramp, a clock step and a clock ramp are injected and read back with a least-squares fix on the ranges
% (position and clock) and on the rates (velocity and drift); off the interval the arrays are untouched, and
% so are the unused slots. A moving true receiver is checked too. Run from anywhere:
% octave --no-gui --quiet --eval test_spoofInject
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'tools'));
randn('seed', 5);  rand('seed', 5);

lat0 = 0.90;  lon0 = 0.10;  h0 = 50;
a = 6378137;  e2 = 6.69437999014e-3;
Nr = a / sqrt(1 - e2 * sin(lat0)^2);
refEcef = [(Nr + h0) * cos(lat0) * cos(lon0), (Nr + h0) * cos(lat0) * sin(lon0), (Nr * (1 - e2) + h0) * sin(lat0)];
C = [-sin(lat0) * cos(lon0), -sin(lat0) * sin(lon0), cos(lat0); -sin(lon0), cos(lon0), 0; ...
     -cos(lat0) * cos(lon0), -cos(lat0) * sin(lon0), -sin(lat0)];                      % rows: N, E, D in ECEF
S = 15;  used = [1 3 4 7 9 12];  RHO = 2.0e7;  CLK = 1500;  DRIFT = 0.3;
up = refEcef / norm(refEcef);
dirs = randn(S, 3);  dirs = dirs - (dirs * up') * up;  dirs = dirs ./ repmat(sqrt(sum(dirs.^2, 2)), 1, 3);
elv = 0.2 + 0.6 * rand(S, 1);
dirs = repmat(cos(elv), 1, 3) .* dirs + repmat(sin(elv), 1, 3) .* repmat(up, S, 1);
satPos = zeros(S, 3);  satVel = zeros(S, 3);  range0 = zeros(S, 1);  rate0 = zeros(S, 1);
for j = used
    satPos(j, :) = refEcef + RHO * dirs(j, :);
    v = randn(1, 3);  satVel(j, :) = 3000 * v / norm(v);
    range0(j) = RHO + CLK;
    rate0(j)  = dirs(j, :) * satVel(j, :)' + DRIFT;
end
u = dirs(used, :);  A = [-u, ones(numel(used), 1)];
fixPos  = @(r)  A \ (r(used) - RHO - CLK);                                   % [dp ecef; clock]
fixRate = @(rr) A \ (rr(used) - sum(u .* satVel(used, :), 2) - DRIFT);       % [dv ecef; drift]
ned = @(x) (C * x(1:3))';
still = [0 0 0];  ll = [lat0 lon0];

fprintf('%-12s %10s %10s %10s %10s | %10s %10s %10s %10s\n', 'case', 'N [m]', 'E [m]', 'D [m]', 'clock [m]', 'vN [m/s]', 'vE', 'vD', 'drift');
%% step 20 m north, read back
scn = struct('kind', 'step', 'tOn', 100, 'tOff', 160, 'size', 20, 'dirNed', [1 0 0]);
[r, rr] = spoofInject(range0, rate0, satPos, satVel, refEcef, still, ll, spoofScenario(130, scn));
xp = fixPos(r);  xv = fixRate(rr);  pStep = ned(xp);  vStep = ned(xv);
fprintf('%-12s %+10.3f %+10.3f %+10.3f %+10.3f | %+10.4f %+10.4f %+10.4f %+10.4f\n', 'step 20 N', pStep, xp(4), vStep, xv(4));
okStep = max(abs(pStep - [20 0 0])) < 1e-3 && abs(xp(4)) < 1e-3 && max(abs(vStep)) < 1e-2 && abs(xv(4)) < 1e-2 ...
    && r(2) == 0 && rr(2) == 0 && r(15) == 0;
%% ramp 1 m/s east, 10 s in
scn = struct('kind', 'ramp', 'tOn', 100, 'tOff', 160, 'size', 1.0, 'dirNed', [0 1 0]);
[r, rr] = spoofInject(range0, rate0, satPos, satVel, refEcef, still, ll, spoofScenario(110, scn));
xp = fixPos(r);  xv = fixRate(rr);  pRamp = ned(xp);  vRamp = ned(xv);
fprintf('%-12s %+10.3f %+10.3f %+10.3f %+10.3f | %+10.4f %+10.4f %+10.4f %+10.4f\n', 'ramp 1 E', pRamp, xp(4), vRamp, xv(4));
okRamp = max(abs(pRamp - [0 10 0])) < 1e-3 && abs(xp(4)) < 1e-3 && max(abs(vRamp - [0 1 0])) < 1e-2 && abs(xv(4)) < 1e-2;
%% clock step 30 m
scn = struct('kind', 'clock', 'tOn', 100, 'tOff', 160, 'size', 30, 'dirNed', [1 0 0]);
[r, rr] = spoofInject(range0, rate0, satPos, satVel, refEcef, still, ll, spoofScenario(130, scn));
xp = fixPos(r);  xv = fixRate(rr);  pClk = ned(xp);  vClk = ned(xv);
fprintf('%-12s %+10.3f %+10.3f %+10.3f %+10.3f | %+10.4f %+10.4f %+10.4f %+10.4f\n', 'clock 30', pClk, xp(4), vClk, xv(4));
okClk = max(abs(pClk)) < 1e-3 && abs(xp(4) - 30) < 1e-3 && max(abs(vClk)) < 1e-6 && abs(xv(4)) < 1e-6;
%% clock ramp 0.5 m/s, 20 s in
scn = struct('kind', 'clockRamp', 'tOn', 100, 'tOff', 160, 'size', 0.5, 'dirNed', [1 0 0]);
[r, rr] = spoofInject(range0, rate0, satPos, satVel, refEcef, still, ll, spoofScenario(120, scn));
xp = fixPos(r);  xv = fixRate(rr);  pCr = ned(xp);  vCr = ned(xv);
fprintf('%-12s %+10.3f %+10.3f %+10.3f %+10.3f | %+10.4f %+10.4f %+10.4f %+10.4f\n', 'clockRamp .5', pCr, xp(4), vCr, xv(4));
okCr = max(abs(pCr)) < 1e-3 && abs(xp(4) - 10) < 1e-3 && max(abs(vCr)) < 1e-6 && abs(xv(4) - 0.5) < 1e-6;
%% off the interval: nothing changes
[r, rr] = spoofInject(range0, rate0, satPos, satVel, refEcef, still, ll, spoofScenario(50, scn));
okOff = isequal(r, range0) && isequal(rr, rate0);
%% a true receiver moving 10 m/s north with a 20 m step: the rates change only by the line-of-sight term
velTrue = (C' * [10 0 0]')';
rateMov = rate0;
for j = used, rateMov(j) = dirs(j, :) * (satVel(j, :) - velTrue)' + DRIFT; end
scn = struct('kind', 'step', 'tOn', 100, 'tOff', 160, 'size', 20, 'dirNed', [1 0 0]);
[r, rr] = spoofInject(range0, rateMov, satPos, satVel, refEcef, velTrue, ll, spoofScenario(130, scn));
xp = fixPos(r);  dv = ned(A \ (rr(used) - rateMov(used)));
fprintf('%-12s %+10.3f %+10.3f %+10.3f %+10.3f | %+10.4f %+10.4f %+10.4f\n', 'moving+step', ned(xp), xp(4), dv);
okMov = max(abs(ned(xp) - [20 0 0])) < 1e-3 && max(abs(dv)) < 1e-2;
%% a bad kind errors
okErr = false;
try
    spoofScenario(130, struct('kind', 'jump', 'tOn', 100, 'tOff', 160, 'size', 1, 'dirNed', [1 0 0]));
catch
    okErr = true;
end
if okStep && okRamp && okClk && okCr && okOff && okMov && okErr
    fprintf('=== test_spoofInject PASS ===\n');
else
    fprintf('step %d ramp %d clock %d clockRamp %d off %d moving %d err %d\n=== test_spoofInject FAIL ===\n', ...
        okStep, okRamp, okClk, okCr, okOff, okMov, okErr);
end
end
