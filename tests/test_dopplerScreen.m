function test_dopplerScreen()
%% test_dopplerScreen.m
% SPF_dopplerScreen, the host-side Doppler episode screen, called epoch by epoch on a synthetic run in the host's
% error-state form and held against tools/dopplerEpisodeCheck on the same records: the two must flag the same
% epochs, and the tool needs no feedback log or fit for it. The run (tests/dopplerSynthRun): east at 20 m/s with a
% braking and a gentle turn under a 5 mrad heading error and an accelerometer bias, the filter measuring the rows
% against the nav less its prior state and the nav taking a tenth of the posterior; the receiver's Doppler wrong
% three times, 0.35 m/s for 5 s, a 0.4 m/s ramp for 20 s through a second braking, and 0.5 m/s for 35 s, longer than
% the hold; a gap without rows. Checks: the flags identical to the tool's; the first two episodes caught within two
% epochs and covered; no flag in the honest manoeuvres; the gap undecided with the coast kept; the long episode held
% through its hold limit and the rows' return afterwards accepted within the grown tolerance; the increment the
% screen forms equal to the IMU's; the state and parameter constructors as documented.
% Run from the repo root (MATLAB, or Octave with tools/octave_shim).
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'tools'));

[run, log, tKf, aux] = dopplerSynthRun(struct('T', 430, 'seed', 5, 'corners', false, 'gap', [350 356], ...
    'arm', [0 0], 'truthPoint', [0 0], 'episodes', [200 205 0.35 0 0; 300 323 0 -0.4 3; 360 395 0.5 0 0]));
K = numel(tKf);

%% the host's screen, epoch by epoch: the nav velocity and the prior in, the update's step written back after it
prm = STRUCT_SPF.dopplerScreenParam;
screen = STRUCT_SPF.dopplerScreenState;
flagHost = false(K, 1);  validHost = false(K, 1);  holdHost = zeros(K, 1);  dvHost = NaN(K, 3);
for k = 1:K
    velNav = run.navVel(aux.jBefore(k), :)';
    priorVel = aux.xPri(k, :)';
    [screen, drop, info] = SPF_dopplerScreen(screen, velNav, priorVel, log{k}.obsMatrix, log{k}.innovation, ...
        log{k}.numMeas, [17 19 22], prm);
    screen.incPrev = aux.inc(k, :);                                     % the host, after its update
    flagHost(k) = drop;  validHost(k) = info.valid;  holdHost(k) = info.hold;
    if k > 1, dvHost(k, :) = info.dvIns; end
    if k == 100, rateRows100 = info.rateRows; end
end

%% the offline tool on the same records: the prior in the records, nothing fitted
runTool = rmfield(run, 'fbVel');
txt = evalc('res = dopplerEpisodeCheck(log, tKf, runTool);');
fprintf('%s', txt);

ok = true;  chk = struct();
chk.same     = isequal(flagHost, res.flag) && strcmp(res.fbModel, 'total');
ep1 = tKf >= 200 & tKf < 205;  ep2 = tKf >= 302 & tKf < 320;
chk.caught1  = find(flagHost & ep1, 1) - find(ep1, 1) <= 2 && mean(flagHost(ep1)) >= 0.7;
chk.caught2  = find(flagHost & ep2, 1) - find(ep2, 1) <= 2 && mean(flagHost(ep2)) >= 0.8;
chk.quiet    = ~any(flagHost(tKf >= 95 & tKf < 170)) && ~any(flagHost(tKf >= 245 & tKf < 270));
chk.gap      = ~any(validHost(tKf >= 350 & tKf < 356)) && ~any(flagHost(tKf >= 350 & tKf < 356)) && ...
    validHost(find(tKf >= 356, 1));
% the long episode: flagged from its start, the hold limit reached at 20 s, the episode reopened at once since the
% rows are still beyond the grown tolerance, and the rows' return accepted within it, so the return costs little
kHold = find(holdHost >= prm.HOLD - 1 & tKf > 360, 1);
chk.holdLimit = ~isempty(kHold) && tKf(kHold) > 378 && tKf(kHold) < 382 && mean(flagHost(tKf >= 361 & tKf < 395)) >= 0.9;
chk.return   = nnz(flagHost(tKf >= 396 & tKf < 430)) <= 8;
% the increment the screen forms is the IMU's: the truth's change plus the heading error's work and the bias, not
% the hand-over, so it stays within a few mm/s of the truth's change on the straight legs
dTrue = [NaN(1, 3); diff(aux.vTrue(aux.jBefore, :))];
straight = tKf > 20 & tKf < 95;
chk.increment = sqrt(mean(sum((dvHost(straight, 1:2) - dTrue(straight, 1:2)) .^ 2, 2))) < 0.01;
chk.structs  = islogical(screen.anchored) && isfield(screen, 'incPrev') && isfield(screen, 'vHatPrev') && ...
    prm.FLOOR == 0.20 && prm.ROWS_SIGN == -1 && prm.RELEASE_GROW == 10;
chk.rateRows = nnz(rateRows100) == 10 && all(find(rateRows100) == (11:20)');
names = fieldnames(chk);
for i = 1:numel(names)
    if ~chk.(names{i}), fprintf('   FAILED sub-check: %s\n', names{i}); ok = false; end
end
fprintf('host flags %d, tool flags %d, differing epochs %d; long episode hold limit at %.1f s; increment error %.4f m/s rms\n', ...
    nnz(flagHost), nnz(res.flag), nnz(flagHost ~= res.flag), tKf(max(kHold, 1)), ...
    sqrt(mean(sum((dvHost(straight, 1:2) - dTrue(straight, 1:2)) .^ 2, 2))));
if ok
    fprintf('=== test_dopplerScreen PASS ===\n');
else
    fprintf('=== test_dopplerScreen FAIL ===\n');
end
end
