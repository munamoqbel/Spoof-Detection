function [res] = feedbackCheck(run, spfMeasLog, win)
%FEEDBACKCHECK  Does the KF update reach the mechanisation, when, and in full? The change of the
%   logged nav solution around each 2 Hz epoch against the increment the filter computed and the
%   state it carries, per NED axis, with the instant and the profile of the change read from the
%   log itself, and the host's bookkeeping of the state between updates read from the records.
%   res = feedbackCheck(run, spfMeasLog)         every epoch
%   res = feedbackCheck(run, spfMeasLog, win)    epochs with win(1) <= tKf <= win(2)
%
%   run: the struct runCompare takes: t [N x 1] host time of the nav samples (100 Hz), navPos
%   [N x 3] [lat lon h] (rad, rad, m), navVel [N x 3] NED (m/s), tKf [K x 1] the 2 Hz epochs,
%   one per record of spfMeasLog; biasEst [N x 6] or [K x 6] optional, accelerometer then gyro
%   (m/s^2, rad/s), the value the mechanisation subtracts plus the bias state; pGyroFirst
%   optional (true when the gyro biases are states 10:12 and the accelerometer biases 13:15).
%   spfMeasLog: struct array or cell of the per-epoch spfMeas records (priorState, postState,
%   numMeas), host layout: position 1:3, velocity 4:6, attitude 7:9, biases 10:15, clock 16:17,
%   pressure altitude 20.
%
%   Position. Per sample the change of navPos (lat lon h to NED metres) minus the integral of
%   navVel over the sample is the part of the motion the velocity does not explain; summed over
%   the half epoch either side of tKf it is the correction the update put into the logged
%   position, whenever and however it was applied. Its running sum from tKf - 0.25 s, regressed
%   on the increment postState - priorState (1:3) across the epochs, is the step profile: the
%   share of the increment present in the position at tKf + tau, from which the instant of the
%   update (10, 50 and 90 % of the step) is read. The sum over the half epochs is regressed
%   jointly on the increments of the records k-1, k and k+1, and on the state x+ of record k
%   (position, velocity, pressure altitude): a host that applies GAIN x+ shows -GAIN on x+ with
%   a correlation near -1; on the increments it shows -GAIN at k and -GAIN (1-GAIN) at k-1.
%   Velocity. The same, with navVel against a line fitted 0.3 to 0.05 s before the step the
%   position profile found and a line 0.05 to 0.25 s after it, both evaluated at the step.
%   Bookkeeping. x- of record k+1 regressed on x+ of record k: the share of each state kept
%   after the feedback (1 - GAIN for a split feedback, 0 for a full one, 1 for a state never fed
%   back), the velocity residual integrated into the position residual (0.5 s times the kept
%   share when the residual is propagated with Phi), and the position residual leaking into the
%   velocity residual (0 when Phi is applied the right way round).
%   Bias. With biasEst the sum of the accelerometer bias increments over the window is set
%   against the change of biasEst, the same -1 / 0 reading for the bias feedback.
%
%   Verdict lines name an axis whose increment does not reach the logged solution, one applied
%   in part, with the opposite sign to the other axes or over more than one epoch (a split
%   feedback), records offset by one epoch from the nav log, a step spread over more than
%   0.05 s (a slewed feedback), a logged solution that receives less of the state than the
%   bookkeeping removes (part of each correction lost between setKF and the mechanisation), and
%   a residual propagated with Phi transposed or without its coupling.
%
%   res: per-epoch t, incPos, incVel, chgPos, chgVel [K x 3] (NaN where not computed), used
%   [K x 1], the lag coefficients coefPos, coefVel [3 x 3] (rows k-1, k, k+1; columns N E D)
%   with their standard errors sePos, seVel [1 x 3], the correlations rPos, rVel [1 x 3], the
%   coefficients on the state gainPos [3 x 3] (rows x+ position, x+ velocity, x+ pressure) and
%   gainVel [1 x 3] with the correlations rxPos, rxVel, the profiles tau [1 x J], profPos,
%   profVel [J x 3], the step times t10, t50, t90 [2 x 3] (rows position, velocity), tStep, the
%   bookkeeping keep (struct per state group), aPV, aVP [1 x 3], the state means xpMean [1 x n],
%   the bias totals incBias and dBias [1 x 3] in m/s^2 (NaN without biasEst), and the verdict
%   lines in res.verdict (cell).
%   MATLAB desktop tool, not codegen. Runs under Octave.

if iscell(spfMeasLog), spfMeasLog = [spfMeasLog{:}]; end
if nargin < 3 || isempty(win), win = [-Inf Inf]; end
K   = numel(spfMeasLog);
t   = run.t(:);
tKf = run.tKf(:);
if numel(tKf) ~= K
    error('feedbackCheck: tKf has %d epochs, spfMeasLog %d records', numel(tKf), K);
end
HALF = 0.25;                                             % half an epoch: the interval owned by each update (s)
LATE = 0.45;                                             % the profile is followed to tKf + LATE, before the next update (s)
TB1 = 0.2;  TB2 = 0.05;                                  % line before the step, tStep - TB1 to tStep - TB2 (s), clear of the previous epoch's late changes
TA1 = 0.05;  TA2 = 0.22;                                 % line after the step, tStep + TA1 to tStep + TA2 (s)
TL1 = 0.30;  TL2 = 0.45;                                 % the late read of the velocity, tStep + TL1 to tStep + TL2 (s)
MG  = 9.80665e-3;                                        % m/s^2 per mg
tau = -HALF:0.01:LATE;  J = numel(tau);
dtNav  = median(diff(t));
tol    = dtNav / 4;                                      % sample-time jitter allowed at the interval ends
navNed = nedFrom(run.navPos);
navVel = run.navVel;
% the position change per sample that navVel does not explain, stamped at the later sample
resPos = [NaN(1, 3); diff(navNed) - 0.5 * (navVel(1:end-1, :) + navVel(2:end, :)) .* repmat(diff(t), 1, 3)];
iAcc = 10:12;  iGyr = 13:15;
if isfield(run, 'pGyroFirst') && ~isempty(run.pGyroFirst) && run.pGyroFirst, iAcc = 13:15;  iGyr = 10:12; end
n = numel(spfMeasLog(1).postState);
iPrs = 20;
if n < 20, iPrs = []; end

% ---- the records ----
XM = NaN(K, n);  XP = NaN(K, n);
hasCov = isfield(spfMeasLog, 'postCov');
cVB = NaN(K, 3, 3);  cVA = NaN(K, 3, 3);  cPV = NaN(K, 3);       % |correlation| velocity x accel bias, velocity x attitude, position-velocity
for k = 1:K
    s = spfMeasLog(k);
    if numel(s.postState) == n && numel(s.priorState) == n
        XM(k, :) = s.priorState(:)';
        XP(k, :) = s.postState(:)';
    end
    if hasCov && isequal(size(s.postCov), [n n])
        P = s.postCov;
        d = sqrt(diag(P));
        if all(isfinite(d)) && all(d > 0)
            C = P ./ (d * d');
            cVB(k, :, :) = abs(C(4:6, iAcc));
            cVA(k, :, :) = abs(C(4:6, 7:9));
            cPV(k, :) = abs([C(1, 4), C(2, 5), C(3, 6)]);
        end
    end
end
INC = XP - XM;
inWin = tKf >= win(1) & tKf <= win(2);
okRec = inWin & all(isfinite(XP(:, 1:6)), 2) & all(isfinite(XM(:, 1:6)), 2);

% ---- pass 1: the position, model-free ----
incPos = NaN(K, 3);  incVel = NaN(K, 3);
chgPos = NaN(K, 3);  chgVel = NaN(K, 3);
cumPos = NaN(K, J, 3);  devVel = NaN(K, J, 3);
used = false(K, 1);
for k = 1:K
    if ~okRec(k), continue; end
    iI = find(t >= tKf(k) - HALF - tol & t < tKf(k) + HALF - tol);   % the samples of this update's half epochs
    iL = find(t >= tKf(k) - HALF - tol & t <= tKf(k) + LATE + tol);  % followed into the next half epoch, before its update
    if numel(iI) < 10 || iI(1) < 2, continue; end
    incPos(k, :) = INC(k, 1:3);
    incVel(k, :) = INC(k, 4:6);
    r = resPos(iI, :);  r(~isfinite(r)) = 0;
    chgPos(k, :) = sum(r, 1);
    r = resPos(iL, :);  r(~isfinite(r)) = 0;
    tI = t(iL) - tKf(k);
    cs = cumsum(r, 1);
    for j = 1:J
        m = find(tI <= tau(j) + tol, 1, 'last');
        if ~isempty(m), cumPos(k, j, :) = cs(m, :); else, cumPos(k, j, :) = 0; end
    end
    used(k) = true;
end
nUsed = nnz(used);
[coefPos, sePos, rPos] = lagFit(incPos, chgPos, used);
profPos = profileOf(incPos, cumPos, used);              % [J x 3] slope of the running sum on the increment
tPos = stepTimes(tau, profPos, sum(coefPos, 1));
t50 = tPos(2, isfinite(tPos(2, :)));
tStep = 0;                                               % the instant of the update, from the position
if ~isempty(t50), tStep = median(t50); end

% ---- pass 2: the velocity, lines either side of the step ----
chgVelLate = NaN(K, 3);
for k = find(used)'
    iB = find(t >= tKf(k) + tStep - TB1 & t <= tKf(k) + tStep - TB2);
    iA = find(t >= tKf(k) + tStep + TA1 & t <= tKf(k) + tStep + TA2);
    iLa = find(t >= tKf(k) + tStep + TL1 & t <= tKf(k) + tStep + TL2);
    iV = find(t >= tKf(k) - HALF - tol & t <= tKf(k) + LATE + tol);
    if numel(iB) < 4 || numel(iA) < 4
        used(k) = false;  chgPos(k, :) = NaN;  continue;
    end
    cB = [ones(numel(iB), 1), t(iB) - tKf(k)] \ navVel(iB, :);   % line before: value at tKf and slope
    cA = [ones(numel(iA), 1), t(iA) - tKf(k)] \ navVel(iA, :);   % line after
    chgVel(k, :) = (cA(1, :) + tStep * cA(2, :)) - (cB(1, :) + tStep * cB(2, :));
    if numel(iLa) >= 4
        cL = [ones(numel(iLa), 1), t(iLa) - tKf(k)] \ navVel(iLa, :);   % line late in the epoch
        chgVelLate(k, :) = (cL(1, :) + tStep * cL(2, :)) - (cB(1, :) + tStep * cB(2, :));
    end
    tV = t(iV) - tKf(k);
    dv = navVel(iV, :) - (repmat(cB(1, :), numel(iV), 1) + tV * cB(2, :));
    devVel(k, :, :) = interp1(tV, dv, tau, 'linear', NaN);
end
nUsed = nnz(used);
[coefPos, sePos, rPos] = lagFit(incPos, chgPos, used);
[coefVel, seVel, rVel] = lagFit(incVel, chgVel, used);
profPos = profileOf(incPos, cumPos, used);
profVel = profileOf(incVel, devVel, used);
tPos = stepTimes(tau, profPos, sum(coefPos, 1));
tVel = stepTimes(tau, profVel, sum(coefVel, 1));

% ---- the change against the state x+ of record k ----
gainPos = NaN(3, 3);  rxPos = NaN(1, 3);  sxPos = NaN(1, 3);
gainVel = NaN(1, 3);  rxVel = NaN(1, 3);  sxVel = NaN(1, 3);
gainVelLate = NaN(1, 3);  sxVelLate = NaN(1, 3);
for a = 1:3
    X = [XP(used, a), XP(used, a + 3)];
    if ~isempty(iPrs) && a == 3, X = [X, XP(used, iPrs)]; end
    [c, rr, ss] = fitThrough(X, chgPos(used, a));
    gainPos(1:numel(c), a) = c;  rxPos(a) = rr;  sxPos(a) = ss;
    [c, rr, ss] = fitThrough(XP(used, a + 3), chgVel(used, a));
    gainVel(a) = c;  rxVel(a) = rr;  sxVel(a) = ss;
    [c, ~, ss] = fitThrough(XP(used, a + 3), chgVelLate(used, a));
    gainVelLate(a) = c;  sxVelLate(a) = ss;
end

% ---- bookkeeping: x- of record k+1 on x+ of record k ----
kk = find(okRec(1:K-1) & okRec(2:K));                   % consecutive usable pairs
groups = {'position', 1:3; 'velocity', 4:6; 'attitude', 7:9; 'accel bias', iAcc; 'gyro bias', iGyr; 'clock', 16:17; 'pressure alt', iPrs};
keep = struct();
keepVal = NaN(1, size(groups, 1));
for g = 1:size(groups, 1)
    idx = groups{g, 2};
    idx = idx(idx <= n);
    kv = NaN(1, numel(idx));
    for i = 1:numel(idx)
        kv(i) = fitThrough(XP(kk, idx(i)), XM(kk + 1, idx(i)));
    end
    kv = kv(isfinite(kv));
    if ~isempty(kv), keepVal(g) = median(kv); end
end
aPV = NaN(1, 3);  aVP = NaN(1, 3);  aPP = NaN(1, 3);  aVV = NaN(1, 3);
for a = 1:3
    c = fitThrough([XP(kk, a), XP(kk, a + 3)], XM(kk + 1, a));          % position residual on position and velocity
    aPP(a) = c(1);  aPV(a) = c(2);
    c = fitThrough([XP(kk, a + 3), XP(kk, a)], XM(kk + 1, a + 3));      % velocity residual on velocity and position
    aVV(a) = c(1);  aVP(a) = c(2);
end
keepVal(1) = median(aPP);  keepVal(2) = median(aVV);     % own-state share from the joint fit, free of the coupling
xpMean = mean(XP(okRec, :), 1);
xpRms  = sqrt(mean(XP(okRec, :).^2, 1));
incAccU = INC(okRec, iAcc);                              % the accelerometer bias increments, for the chase / tracking reading
accSame = NaN(1, 3);
for a = 1:3
    nz = incAccU(:, a) ~= 0 & isfinite(incAccU(:, a));
    if any(nz), accSame(a) = 100 * mean(sign(incAccU(nz, a)) == sign(mean(incAccU(nz, a)))); end
end
medVB = NaN(3, 3);  medVA = NaN(3, 3);  medPV = NaN(1, 3);
if hasCov
    for i = 1:3
        for j = 1:3
            v = cVB(okRec, i, j);  v = v(isfinite(v));  if ~isempty(v), medVB(i, j) = median(v); end
            v = cVA(okRec, i, j);  v = v(isfinite(v));  if ~isempty(v), medVA(i, j) = median(v); end
        end
        v = cPV(okRec, i);  v = v(isfinite(v));  if ~isempty(v), medPV(i) = median(v); end
    end
end
xmMean = mean(XM(okRec, :), 1);
xmRms  = sqrt(mean(XM(okRec, :).^2, 1));

% ---- accelerometer bias: the sum of the increments over the window against the change of biasEst ----
incBias = NaN(1, 3);  dBias = NaN(1, 3);
ku = find(used);
if isfield(run, 'biasEst') && ~isempty(run.biasEst) && size(run.biasEst, 1) > 1 && ~isempty(ku)
    B = run.biasEst(:, 1:3);
    incAcc = INC(:, iAcc);  incAcc(~isfinite(incAcc)) = 0;
    if size(B, 1) == K                                   % per epoch: the value after the update of epoch k
        k0 = max(ku(1) - 1, 1);
        dBias   = B(ku(end), :) - B(k0, :);
        incBias = sum(incAcc(k0+1:ku(end), :), 1);
    elseif size(B, 1) == numel(t)                        % per nav sample: before the first and after the last update
        i0 = find(t <= tKf(ku(1)) - TB2, 1, 'last');
        i1 = find(t >= tKf(ku(end)) + HALF, 1, 'first');
        if ~isempty(i0) && ~isempty(i1)
            dBias   = B(i1, :) - B(i0, :);
            incBias = sum(incAcc(ku(1):ku(end), :), 1);
        end
    end
end

% ---- printout ----
tUsed = tKf(used);
fprintf('feedbackCheck: %d updates between %.1f and %.1f s (%d epochs in the window, %d without states or samples), nav samples every %.3f s\n', ...
    nUsed, min([tUsed; Inf]), max([tUsed; -Inf]), nnz(inWin), nnz(inWin) - nUsed, dtNav);
if dtNav > 0.05
    fprintf('   WARNING: nav samples every %.3f s, the profile needs a log of at least 20 Hz\n', dtNav);
end
fprintf('   position: the change of navPos that navVel does not explain, summed over tKf - %.2f to tKf + %.2f s\n', HALF, HALF);
fprintf('   velocity: navVel lines %.2f to %.2f s before and %.2f to %.2f s after the step at tKf %+.2f s, both at the step\n', TB1, TB2, TA1, TA2, tStep);
fprintf('                                           N           E           D\n');
printBlock('position increment x+ - x- (1:3)', 'change of navPos across the epoch', incPos(used, :), chgPos(used, :), coefPos, sePos, rPos, 'm');
fprintf('      on the state x+ of record k\n');
fprintf('         x+ position                 %+11.3f %+11.3f %+11.3f       (-GAIN when the host applies GAIN x+)\n', gainPos(1, :));
fprintf('         x+ velocity                 %+11.3f %+11.3f %+11.3f       (0 unless the velocity state is applied to the position)\n', gainPos(2, :));
if ~isempty(iPrs)
    fprintf('         x+ pressure altitude              -           - %+11.3f       (0 unless the pressure state is applied to the altitude)\n', gainPos(3, 3));
end
fprintf('         correlation with x+ position %+10.2f %+11.2f %+11.2f\n', rxPos);
printBlock('velocity increment x+ - x- (4:6)', 'change of navVel across the step', incVel(used, :), chgVel(used, :), coefVel, seVel, rVel, 'm/s');
fprintf('      on the state x+ of record k\n');
fprintf('         x+ velocity                 %+11.3f %+11.3f %+11.3f       (-GAIN when the host applies GAIN x+)\n', gainVel);
fprintf('         standard error              %11.3f %11.3f %11.3f\n', sxVel);
fprintf('         correlation                 %+11.2f %+11.2f %+11.2f\n', rxVel);
fprintf('      read %.2f to %.2f s after the step instead (the correction still there late in the epoch)\n', TL1, TL2);
fprintf('         mean                        %+11.4f %+11.4f %+11.4f   m/s\n', mean(chgVelLate(used, :), 1));
fprintf('         on x+ velocity              %+11.3f %+11.3f %+11.3f       (equal to the early read when the correction stays)\n', gainVelLate);
fprintf('         standard error              %11.3f %11.3f %11.3f\n', sxVelLate);
printProfile(tau, profPos, profVel, tPos, tVel);
fprintf('   accelerometer bias increment x+ - x- (%d:%d), in mg\n', iAcc(1), iAcc(end));
fprintf('      mean                         %+11.4f %+11.4f %+11.4f   mg\n', mean(incAccU, 1) / MG);
fprintf('      median                       %+11.4f %+11.4f %+11.4f   mg\n', median(incAccU, 1) / MG);
fprintf('      with the sign of the mean    %10.0f%% %10.0f%% %10.0f%%   of the epochs (near 50 tracking a wandering bias, well above it a chase)\n', accSame);
fprintf('      x+ accel bias  mean / rms     %+8.3f /%6.3f %+8.3f /%6.3f %+8.3f /%6.3f   mg\n', [xpMean(iAcc); xpRms(iAcc)] / MG);
if hasCov
    fprintf('   coupling the covariance carries (postCov), |correlation| between the states, median over the epochs\n');
    fprintf('      position - velocity, same axis  %8.2f    %8.2f    %8.2f\n', medPV);
    fprintf('      velocity x accel bias            bX        bY        bZ        velocity x attitude   roll      pitch     yaw\n');
    rowName = {'vN', 'vE', 'vD'};
    for i = 1:3
        fprintf('         %s                        %6.2f    %6.2f    %6.2f           %s                %6.2f    %6.2f    %6.2f\n', ...
            rowName{i}, medVB(i, :), rowName{i}, medVA(i, :));
    end
    fprintf('      (a velocity row near zero on every bias and attitude column is an axis whose error the filter cannot attribute to\n');
    fprintf('      the IMU: no coupling in F for that axis)\n');
end
fprintf('   state carried between updates (records only): x- of record k+1 on x+ of record k\n');
fprintf('      share of x+ kept in x-, median over the states of the group\n');
for g = 1:size(groups, 1)
    if ~isempty(groups{g, 2}) && isfinite(keepVal(g))
        fprintf('         %-14s %6.3f\n', groups{g, 1}, keepVal(g));
    end
end
fprintf('      (0 full feedback, 1 - GAIN a split feedback, 1 a state never fed back)\n');
fprintf('      position residual on x+ position %+10.3f %+11.3f %+11.3f\n', aPP);
fprintf('                           x+ velocity %+10.3f %+11.3f %+11.3f       (0.5 s x the kept share when the residual is propagated with Phi)\n', aPV);
fprintf('      velocity residual on x+ velocity %+10.3f %+11.3f %+11.3f\n', aVV);
fprintf('                           x+ position %+10.3f %+11.3f %+11.3f       (0 with Phi the right way round; 0.5 s x the kept share with Phi transposed)\n', aVP);
fprintf('      x+ position  mean / rms       %+8.3f /%6.3f %+8.3f /%6.3f %+8.3f /%6.3f   m\n', [xpMean(1:3); xpRms(1:3)]);
fprintf('      x+ velocity  mean / rms       %+8.3f /%6.3f %+8.3f /%6.3f %+8.3f /%6.3f   m/s\n', [xpMean(4:6); xpRms(4:6)]);
fprintf('      x- position  mean / rms       %+8.3f /%6.3f %+8.3f /%6.3f %+8.3f /%6.3f   m\n', [xmMean(1:3); xmRms(1:3)]);
fprintf('      x- velocity  mean / rms       %+8.3f /%6.3f %+8.3f /%6.3f %+8.3f /%6.3f   m/s\n', [xmMean(4:6); xmRms(4:6)]);
if ~isempty(iPrs)
    fprintf('      x+ pressure altitude  mean / rms  %+8.3f /%6.3f m,  x- %+8.3f /%6.3f m\n', xpMean(iPrs), xpRms(iPrs), xmMean(iPrs), xmRms(iPrs));
end
if all(isfinite(dBias))
    ratio = dBias ./ incBias;
    ratio(abs(incBias) < 1e-9) = NaN;
    fprintf('   accelerometer bias (%d:%d) over the window\n', iAcc(1), iAcc(end));
    fprintf('      sum of the increments        %+11.3f %+11.3f %+11.3f   mg\n', incBias / MG);
    fprintf('      change of biasEst            %+11.3f %+11.3f %+11.3f   mg\n', dBias / MG);
    fprintf('      ratio change / sum           %+11.2f %+11.2f %+11.2f       (-1 fed back in full, estimate minus truth; 0 never)\n', ratio);
end
fprintf('   reading: on the increments, -1 at k is the increment applied in full at its own epoch (error state estimate minus\n');
fprintf('   truth), 0 one that never reaches the logged solution, -GAIN at k with -GAIN (1-GAIN) at k-1 a split feedback; on the\n');
fprintf('   state, -GAIN is the share of x+ the logged solution receives, to be compared with the share the bookkeeping removes\n');

% ---- verdicts ----
verdict = {};
rmsPos = sqrt(mean(incPos(used, :).^2, 1));
rmsVel = sqrt(mean(incVel(used, :).^2, 1));
if nUsed < 10
    verdict{end+1} = sprintf('only %d usable updates: no reading', nUsed);
else
    % on an axis where the bookkeeping keeps a residual and the change follows the state x+ closely, the
    % state-based reading replaces the increment-based one
    byStatePos = splitAxes(gainPos(1, :), sxPos, rxPos, keepVal(1), xpRms(1:3), 1e-3);
    byStateVel = splitAxes(gainVel, sxVel, rxVel, keepVal(2), xpRms(4:6), 1e-4);
    verdict = [verdict, axisVerdicts('position', coefPos, rmsPos, 1e-4, byStatePos)];
    verdict = [verdict, axisVerdicts('velocity', coefVel, rmsVel, 1e-5, byStateVel)];
    verdict = [verdict, stepVerdicts('position', tPos, coefPos, rmsPos, 1e-4)];
    verdict = [verdict, stepVerdicts('velocity', tVel, coefVel, rmsVel, 1e-5)];
    [vS, nS] = splitVerdicts(gainPos(1, :), keepVal(1), 'position', byStatePos);
    verdict = [verdict, vS];  notes = nS;
    [vS, nS] = splitVerdicts(gainVel, keepVal(2), 'velocity', byStateVel);
    verdict = [verdict, vS];  notes = [notes, nS];
    for i = 1:numel(notes), fprintf('   NOTE: %s\n', notes{i}); end
    ax = 'NED';
    for a = 1:3
        early = gainVel(a);  late = gainVelLate(a);
        if isfinite(early) && isfinite(late) && abs(early) > 3 * sxVel(a) && abs(early) > 0.03 && abs(late) < 0.3 * abs(early)
            verdict{end+1} = sprintf('the %s velocity correction applied at the step (%+.3f of x+) is gone by tKf %+.2f s (%+.3f of x+ left): something restores the velocity between the updates', ...
                ax(a), early, tStep + TL1, late);
        end
    end
    if isfinite(keepVal(2)) && keepVal(2) > 0.2 && any(abs(aVP) > 0.1)
        verdict{end+1} = sprintf('the velocity residual receives %s x the position residual each epoch: the residual is propagated with Phi transposed (kfState * phi in place of (phi * kfState'')'')', sprintf('%+.2f ', aVP));
    end
    if isfinite(keepVal(1)) && keepVal(1) > 0.2 && all(abs(aPV) < 0.1 * keepVal(1)) && all(abs(aVP) < 0.1)
        verdict{end+1} = 'the position residual does not receive the velocity residual: the residual is carried without the Phi coupling';
    end
    if ~isempty(iPrs) && isfinite(gainPos(3, 3)) && abs(gainPos(3, 3)) > 0.1 && xpRms(iPrs) > 1e-3
        verdict{end+1} = sprintf('the altitude moves by %+.2f x the pressure altitude state at each update, on top of the position feedback', gainPos(3, 3));
    end
end
if all(isfinite(dBias)) && any(abs(incBias) > 1e-5)
    big = abs(incBias) > 1e-5;
    rb = dBias(big) ./ incBias(big);
    if all(abs(rb) < 0.3)
        verdict{end+1} = 'the accelerometer bias increments do not reach biasEst over the window';
    elseif any(abs(rb) < 0.3)
        verdict{end+1} = sprintf('the accelerometer bias increments reach biasEst on some axes only (ratios %s)', sprintf('%+.2f ', rb));
    end
end
for i = 1:numel(verdict)
    fprintf('   VERDICT: %s\n', verdict{i});
end

keep = cell2struct(num2cell(keepVal(:)), regexprep(groups(:, 1), ' ', '_'), 1);
res = struct('t', tKf, 'incPos', incPos, 'incVel', incVel, 'chgPos', chgPos, 'chgVel', chgVel, 'used', used, ...
    'coefPos', coefPos, 'coefVel', coefVel, 'sePos', sePos, 'seVel', seVel, 'rPos', rPos, 'rVel', rVel, ...
    'gainPos', gainPos, 'gainVel', gainVel, 'gainVelLate', gainVelLate, 'rxPos', rxPos, 'rxVel', rxVel, 'chgVelLate', chgVelLate, ...
    'accSame', accSame, 'corrVelBias', medVB, 'corrVelAtt', medVA, 'corrPosVel', medPV, ...
    'tau', tau, 'profPos', profPos, 'profVel', profVel, 'tStep', tStep, ...
    't10', [tPos(1, :); tVel(1, :)], 't50', [tPos(2, :); tVel(2, :)], 't90', [tPos(3, :); tVel(3, :)], ...
    'keep', keep, 'aPP', aPP, 'aPV', aPV, 'aVV', aVV, 'aVP', aVP, 'xpMean', xpMean, 'xmMean', xmMean, ...
    'incBias', incBias, 'dBias', dBias, 'nUsed', nUsed);
res.verdict = verdict;
end

%------------------------------------------------------------------------
function printBlock(nameInc, nameChg, inc, chg, coef, se, r, unit)
mInc = mean(inc, 1);
same = NaN(1, 3);
for a = 1:3
    nz = inc(:, a) ~= 0;
    if any(nz), same(a) = 100 * mean(sign(inc(nz, a)) == sign(mInc(a))); end
end
fprintf('   %s\n', nameInc);
fprintf('      mean                         %+11.4f %+11.4f %+11.4f   %s\n', mInc, unit);
fprintf('      median                       %+11.4f %+11.4f %+11.4f   %s\n', median(inc, 1), unit);
fprintf('      with the sign of the mean    %10.0f%% %10.0f%% %10.0f%%   of the epochs\n', same);
fprintf('   %s\n', nameChg);
fprintf('      mean                         %+11.4f %+11.4f %+11.4f   %s\n', mean(chg, 1), unit);
fprintf('      on the increment of record k-1 %+9.3f %+11.3f %+11.3f\n', coef(1, :));
fprintf('                             k   %+11.3f %+11.3f %+11.3f\n', coef(2, :));
fprintf('                             k+1 %+11.3f %+11.3f %+11.3f\n', coef(3, :));
fprintf('         standard error          %11.3f %11.3f %11.3f\n', se);
fprintf('      correlation with record k    %+11.2f %+11.2f %+11.2f\n', r);
end

function printProfile(tau, profPos, profVel, tPos, tVel)
show = [-0.2 -0.1 -0.05 -0.02 -0.01 0 0.01 0.02 0.05 0.1 0.2 0.25 0.3 0.35 0.4 0.45];
j = zeros(size(show));
for i = 1:numel(show), [~, j(i)] = min(abs(tau - show(i))); end
fprintf('   step profile: share of the increment in the logged solution at tKf + tau (regression across the epochs)\n');
fprintf('      tau (s)      %s\n', sprintf('%6.2f', show));
ax = 'NED';
for a = 1:3
    fprintf('      position %s   %s\n', ax(a), sprintf('%6.2f', -profPos(j, a)));
end
for a = 1:3
    fprintf('      velocity %s   %s\n', ax(a), sprintf('%6.2f', -profVel(j, a)));
end
fprintf('      (sign flipped: 1.00 is the increment subtracted in full; before the step the velocity rows show the error the\n');
fprintf('      filter is about to correct, not a feedback)\n');
fprintf('      step at 10 / 50 / 90 %% of its size, position   %s\n', stepText(tPos));
fprintf('                                     velocity   %s\n', stepText(tVel));
end

function [s] = stepText(tS)
ax = 'NED';
s = '';
for a = 1:3
    if all(isfinite(tS(:, a)))
        s = [s, sprintf('%s %+.2f %+.2f %+.2f s   ', ax(a), tS(1, a), tS(2, a), tS(3, a))];
    else
        s = [s, sprintf('%s    -     -     -     ', ax(a))];
    end
end
end

function [c, r, se] = fitThrough(X, y)
% least squares through the origin, rows with a NaN and columns without variation left out; r the
% correlation of y with the first column, se the standard error of the first coefficient
c = NaN(size(X, 2), 1);  r = NaN;  se = NaN;
ok = all(isfinite([X, y]), 2);
live = sum(X(ok, :).^2, 1) > 0;
if nnz(ok) > 10 && any(live)
    Xo = X(ok, live);  yo = y(ok);
    c(live) = Xo \ yo;
    resid = yo - Xo * c(live);
    s2 = sum(resid.^2) / max(numel(yo) - nnz(live), 1);
    if live(1)
        V = s2 * inv(Xo' * Xo);
        se = sqrt(V(1, 1));
    end
    dx = X(ok, 1) - mean(X(ok, 1));  dy = yo - mean(yo);
    den = sqrt(sum(dx.^2) * sum(dy.^2));
    if den > 0, r = sum(dx .* dy) / den; end
end
end

function [coef, se, r] = lagFit(inc, chg, used)
% chg(k) = coef(1) inc(k-1) + coef(2) inc(k) + coef(3) inc(k+1), per axis, over the epochs whose three
% records are usable; se the largest standard error of the three; r the plain correlation with inc(k)
K = size(inc, 1);
coef = NaN(3, 3);  se = NaN(1, 3);  r = NaN(1, 3);
k = find(used);  k = k(k >= 2 & k <= K - 1);
k = k(used(k - 1) & used(k + 1));
for a = 1:3
    y = chg(k, a);
    X = [inc(k - 1, a), inc(k, a), inc(k + 1, a)];
    ok = all(isfinite([X, y]), 2);
    if nnz(ok) > 10 && all(sum(X(ok, :).^2, 1) > 0)
        Xo = X(ok, :);  yo = y(ok);
        c = Xo \ yo;
        coef(:, a) = c;
        resid = yo - Xo * c;
        s2 = sum(resid.^2) / max(numel(yo) - 3, 1);
        se(a) = sqrt(max(diag(s2 * inv(Xo' * Xo))));
        dx = Xo(:, 2) - mean(Xo(:, 2));  dy = yo - mean(yo);
        den = sqrt(sum(dx.^2) * sum(dy.^2));
        if den > 0, r(a) = sum(dx .* dy) / den; end
    end
end
end

function [prof] = profileOf(inc, cum, used)
% slope of cum(k, j, a) on inc(k, a) across the epochs, per tau and axis
J = size(cum, 2);
prof = NaN(J, 3);
for a = 1:3
    x = inc(used, a);
    if sum(x.^2) > 0
        for j = 1:J
            c = cum(used, j, a);
            ok = isfinite(c);
            if nnz(ok) > 10, prof(j, a) = sum(x(ok) .* c(ok)) / sum(x(ok).^2); end
        end
    end
end
end

function [tS] = stepTimes(tau, prof, total)
% tau at which the profile first reaches 10, 50 and 90 % of its own final value, per axis, on the axes
% whose total coefficient is not negligible
tS = NaN(3, 3);
for a = 1:3
    if isfinite(total(a)) && abs(total(a)) > 0.1 && isfinite(prof(end, a)) && abs(prof(end, a)) > 0.05
        p = prof(:, a) / prof(end, a);
        lev = [0.1 0.5 0.9];
        for i = 1:3
            j = find(p >= lev(i), 1, 'first');
            if ~isempty(j), tS(i, a) = tau(j); end
        end
    end
end
end

function [v] = axisVerdicts(what, coef, rmsInc, tiny, byState)
ax = 'NED';
v = {};
total = sum(coef, 1);
live = rmsInc > tiny & isfinite(total) & ~byState;
for a = 1:3
    if ~live(a)
        continue;
    elseif abs(total(a)) < 0.3
        v{end+1} = sprintf('the %s %s increment does not reach the logged solution (coefficients %+.2f %+.2f %+.2f at k-1, k, k+1): the mechanisation keeps its %s %s through the update', ...
            ax(a), what, coef(:, a), ax(a), what);
    elseif abs(coef(2, a)) < 0.3 * abs(total(a)) && abs(coef(1, a)) >= 0.5 * abs(total(a))
        v{end+1} = sprintf('the %s %s change matches the increment of record k-1 alone (coefficients %+.2f %+.2f %+.2f): the records are one epoch ahead of the nav log', ...
            ax(a), what, coef(:, a));
    elseif abs(coef(2, a)) < 0.3 * abs(total(a)) && abs(coef(3, a)) >= 0.5 * abs(total(a))
        v{end+1} = sprintf('the %s %s change matches the increment of record k+1 alone (coefficients %+.2f %+.2f %+.2f): the records are one epoch behind the nav log', ...
            ax(a), what, coef(:, a));
    elseif abs(coef(1, a)) >= 0.5 * abs(coef(2, a)) && abs(coef(2, a)) >= 0.05
        v{end+1} = sprintf('the %s %s increment is applied through a split feedback: %+.2f of it at its own epoch, the rest from the state over the following epochs (coefficients %+.2f %+.2f %+.2f)', ...
            ax(a), what, coef(2, a), coef(:, a));
    elseif abs(total(a)) < 0.7 || abs(total(a)) > 1.3
        v{end+1} = sprintf('the %s %s increment is applied in part (coefficient %+.2f): a feedback gain below one', ax(a), what, total(a));
    end
end
full = live & abs(total) >= 0.7 & abs(total) <= 1.3;
if nnz(full) >= 2 && any(sign(total(full)) ~= sign(total(find(full, 1))))
    v{end+1} = sprintf('the %s coefficients differ in sign across the axes (%s): one axis is applied with the opposite sign', what, sprintf('%+.2f ', total));
end
end

function [v] = stepVerdicts(what, tS, coef, rmsInc, tiny)
ax = 'NED';
v = {};
total = sum(coef, 1);
for a = 1:3
    if rmsInc(a) > tiny && all(isfinite(tS(:, a))) && abs(total(a)) >= 0.3 && tS(3, a) - tS(1, a) > 0.05
        v{end+1} = sprintf('the %s %s step is spread over %.2f s (10 to 90 %% between tKf %+.2f and %+.2f s): a slewed feedback, not a jump', ...
            ax(a), what, tS(3, a) - tS(1, a), tS(1, a), tS(3, a));
    end
end
end

function [by] = splitAxes(gain, se, rx, keepShare, xpRms, tiny)
% the axes read from the state: the bookkeeping keeps a residual (a split feedback) and the change of the
% logged solution follows x+ (a correlation of 0.8, or a coefficient three standard errors from zero)
by = false(1, 3);
if isfinite(keepShare) && keepShare >= 0.2 && keepShare <= 0.98
    by = xpRms > tiny & isfinite(gain) & ((isfinite(rx) & abs(rx) >= 0.8) | (isfinite(se) & abs(gain) > 3 * se));
end
end

function [v, note] = splitVerdicts(gain, keepShare, what, byState)
% the share of x+ the logged solution receives against the share the bookkeeping removes
ax = 'NED';
v = {};  note = {};
removed = 1 - keepShare;
for a = find(byState)
    if abs(gain(a)) < 0.05
        v{end+1} = sprintf('the logged %s %s does not move with the state x+ (%+.2f per epoch) while the bookkeeping removes %.2f of it: the %s feedback never reaches the mechanisation', ...
            ax(a), what, gain(a), removed, what);
    elseif abs(gain(a)) < 0.7 * removed
        v{end+1} = sprintf('the logged %s %s receives %+.2f of the state x+ per epoch while the bookkeeping removes %.2f of it (keeps %.2f): the difference is lost between setKF and the mechanisation every epoch', ...
            ax(a), what, gain(a), removed, keepShare);
    elseif abs(gain(a)) > 1.3 * removed
        v{end+1} = sprintf('the logged %s %s receives %+.2f of the state x+ per epoch while the bookkeeping removes only %.2f of it: the mechanisation gets more than the state gives up', ...
            ax(a), what, gain(a), removed);
    else
        note{end+1} = sprintf('%s %s: a split feedback, the logged solution receives %+.2f of x+ per epoch and the bookkeeping removes %.2f, consistent', ...
            ax(a), what, gain(a), removed);
    end
end
end

function [e] = nedFrom(p)
% [lat lon h] rows to NED metres from the first row (WGS84 radii at each row)
L  = p(:, 1);  h = p(:, 3);
a  = 6378137;  e2 = 6.69437999014e-3;
RN = a * (1 - e2) ./ (1 - e2 * sin(L).^2).^1.5;
RE = a ./ sqrt(1 - e2 * sin(L).^2);
e = [(p(:, 1) - p(1, 1)) .* (RN + h), ...
     (p(:, 2) - p(1, 2)) .* (RE + h) .* cos(L), ...
    -(p(:, 3) - p(1, 3))];
end
