function test_covPropCheck()
%% test_covPropCheck.m
% Synthetic check of tools/covPropCheck.m: a log whose P+ follows the standard propagation and update must read every
% cross block as kept (ratio 1) with the chain equal to the logged sigma; a log whose P+ has the (pos, vel) x (att,
% bias) cross terms zeroed after every update must read vel-tilt, vel-accel bias and tilt-gyro bias as removed each
% epoch with pos-vel and tilt-gyro bias kept; pDiag equal to diag(postCov) is recognised, a scaled one is not.
% Run from anywhere: octave --no-gui --quiet --eval test_covPropCheck
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'tools'));
rand('seed', 11);  randn('seed', 11);
G = 9.80665;  n = 22;  K = 200;  tKf = 0.5 * (0:K-1)';
tOut = tKf(100);  tReacq = tKf(161);
sig0 = [2 2 2, 0.03 0.03 0.04, 0.5e-3 0.5e-3 2e-3, 7e-4*[1 1 1], 2e-6*[1 1 1], 3, 0.1, ones(1, 5)];
F = zeros(n);  F(1:3,4:6) = eye(3);  F(4:6,7:9) = [0 G 0; -G 0 0; 0 0 0];  F(4:6,10:12) = eye(3);  F(7:9,13:15) = eye(3);
Phi = expm(F * 0.5);
logA = cell(K, 1);  logB = cell(K, 1);
PA = diag(sig0.^2);  PB = PA;
for k = 1:K
    inGap = tKf(k) > tOut & tKf(k) < tReacq;
    if inGap
        H = zeros(1, n);  H(3) = 1;  R = 25;                        % the barometer row
    else
        u = randn(6, 3);  u = u ./ repmat(sqrt(sum(u.^2, 2)), 1, 3);
        H = zeros(6, n);  H(:, 1:3) = u;  H(:, 16) = 1;  R = 4 * eye(6);    % six ranges and the clock
    end
    for c = 1:2
        if c == 1, P = PA; else, P = PB; end
        if k > 1, P = Phi * P * Phi'; end
        S = H * P * H' + R;  Kg = (P * H') / S;
        P = P - Kg * H * P;  P = (P + P') / 2;
        if c == 2, P(1:6, 7:15) = 0;  P(7:15, 1:6) = 0; end          % the fault: cross terms dropped after the update
        m = size(H, 1);
        rec = struct('innovation', zeros(m, 1), 'innovationCov', S, 'obsMatrix', H, 'measNoiseCov', R, ...
            'numMeas', m, 'priorState', zeros(n, 1), 'postState', zeros(n, 1), 'postCov', P);
        if c == 1, PA = P;  logA{k} = rec; else, PB = P;  logB{k} = rec; end
    end
end
pdA = cell2mat(cellfun(@(s) diag(s.postCov)', logA, 'UniformOutput', false));
% propTel logs: right (one interval, reset after use), reset before use (identity), reset missed (two intervals)
telOk = cell(K, 1);  telEarly = cell(K, 1);  telMissed = cell(K, 1);
for k = 1:K
    telOk{k}     = struct('accumPhi', Phi, 'accumQ', zeros(n));
    telEarly{k}  = struct('accumPhi', eye(n), 'accumQ', zeros(n));
    telMissed{k} = struct('accumPhi', Phi * Phi, 'accumQ', zeros(n));
end
fprintf('--- A: consistent host ---\n');
resA = covPropCheck(logA, tKf, tOut, tReacq, [], pdA, telOk);
fprintf('--- A with the accumulator reset before the gate reads it ---\n');
resE = covPropCheck(logA, tKf, tOut, tReacq, [], [], telEarly);
fprintf('--- A with a missed reset ---\n');
resM = covPropCheck(logA, tKf, tOut, tReacq, [], [], telMissed);
fprintf('--- B: cross terms dropped after every update ---\n');
resB = covPropCheck(logB, tKf, tOut, tReacq, [], 1.5 * pdA);
fprintf('--- A without an outage ---\n');
resN = covPropCheck(logA, tKf, [], []);
okA = all(abs(resA.medianAided - 1) < 1e-6) && all(abs(resA.medianOutage - 1) < 1e-6) && resA.pDiagMatch == true ...
   && max(abs(resA.sigChain(resA.inOutage, :) - resA.sigLogged(resA.inOutage, :))) < 1e-9 ...
   && abs(resA.hphRatio(2) - 1) < 1e-9;
okB = all(resB.medianOutage([1 4]) > 0.7) && all(resB.medianOutage(2:3) < 0.2) && resB.pDiagMatch == false ...
   && resB.sigChain(160, 1) > 2 * resB.sigLogged(160, 1);           % the fault leaves pos-vel and tilt-gyro bias intact
okN = all(abs(resN.medianAided - 1) < 1e-6) && ~any(resN.inOutage) && isnan(resN.pDiagMatch);
okT = abs(median(resA.accInterval(2:end)) - 0.5) < 1e-9 && all(abs(resA.accCoupling(2, :) - 1) < 1e-9) ...
   && all(all(abs(resA.ratioHost(2:end, :) - 1) < 1e-6)) ...
   && abs(median(resE.accInterval(2:end))) < 1e-12 && abs(median(resM.accInterval(2:end)) - 1.0) < 1e-9 ...
   && ~any(isfinite(resN.accInterval));
if okA && okB && okN && okT
    fprintf('=== test_covPropCheck PASS ===\n');
else
    fprintf('okA %d okB %d okN %d okT %d\n', okA, okB, okN, okT);
    fprintf('=== test_covPropCheck FAIL ===\n');
end
end
