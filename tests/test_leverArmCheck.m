function test_leverArmCheck()
%% test_leverArmCheck.m
% tools/leverArmCheck on synthetic records: a vehicle at 15 m/s driving straight legs and 15 deg/s turns, GNSS
% velocity rows along N, E, D with 0.1 m/s noise; the innovations carry omega x (arm error) in the turns. With a
% 0.5 m along-body arm error the across-track slope must read 0.5 m within 0.1 and the verdict name the arm; with
% no arm error the slopes must sit within two sigma and the verdict clear it. Run from anywhere:
% octave --no-gui --quiet --eval test_leverArmCheck
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'tools'));
rng(6);
dt = 0.01;  T = 600;  t = (0:dt:T - dt)';  N = numel(t);
R2D = 180 / pi;
% the yaw: 60 s straight, 12 s turns at 15 deg/s, alternating sense
rate = zeros(N, 1);
for k = 1:8
    t0 = 60 + (k - 1) * 72;  sgn = (-1) ^ k;
    rate(t >= t0 & t < t0 + 12) = sgn * 15 / R2D;
end
yaw = cumsum(rate) * dt;
run = struct('t', t, 'navAtt', [zeros(N, 2), yaw]);
tKf = (0.5:0.5:T - 0.5)';  K = numel(tKf);
maxMeas = 51;  n = 22;
function [log] = records(armErr)
    log = repmat(struct('numMeas', uint8(0), 'obsMatrix', zeros(maxMeas, n), 'innovation', zeros(maxMeas, 1)), K, 1);
    for k = 1:K
        i = find(t >= tKf(k), 1);
        w = [0; 0; rate(i)];
        C = [cos(yaw(i)) -sin(yaw(i)) 0; sin(yaw(i)) cos(yaw(i)) 0; 0 0 1];   % body to NED
        dv = C * cross(w, armErr);                                           % the uncompensated antenna motion in NED
        H = zeros(maxMeas, n);  H(1:3, 4:6) = eye(3);
        y = zeros(maxMeas, 1);  y(1:3) = dv + 0.1 * randn(3, 1);
        log(k).numMeas = uint8(3);  log(k).obsMatrix = H;  log(k).innovation = y;
    end
end
fprintf('--- a 0.5 m along-body arm error ---\n');
r1 = leverArmCheck(records([0.5; 0; 0]), tKf, run);
ok1 = abs(r1.slope(2) - 0.5) < 0.1 && abs(r1.slope(1)) < 0.1 && r1.rmsTurn(2) > 1.3 * r1.rmsStraight(2) && r1.nTurn > 100;
fprintf('--- no arm error ---\n');
r0 = leverArmCheck(records([0; 0; 0]), tKf, run);
ok0 = all(abs(r0.slope) < 2 * r0.slopeSig + 0.05) && all(r0.rmsTurn < 1.3 * r0.rmsStraight);
fprintf('--- a 0.3 m cross-body arm error, window on the second half ---\n');
r2 = leverArmCheck(records([0; 0.3; 0]), tKf, run, [300 600]);
ok2 = abs(r2.slope(1) + 0.3) < 0.1 && abs(r2.slope(2)) < 0.1 && r2.nTurn < r1.nTurn && abs(r2.slope(1)) > 2 * r2.slopeSig(1);   % named by the regression alone
fprintf('along-body arm read across track: %s | no arm clears: %s | cross-body arm read along track, windowed: %s\n', pf(ok1), pf(ok0), pf(ok2));
if ok1 && ok0 && ok2, fprintf('=== test_leverArmCheck PASS ===\n'); else, fprintf('=== test_leverArmCheck FAIL ===\n'); end
end

function s = pf(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end
