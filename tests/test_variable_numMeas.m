%% test_variable_numMeas.m
% The gate must accept a measurement count that changes every epoch
% (satellites rise/set). Feeds the fed FSM (protectedNav) the harness
% scenario with 5..8 satellites visible per epoch, two ways:
%   (a) inputs padded to CST_spfParam.MAX_MEAS rows (Coder style)
%   (b) inputs at their exact size (numMeas rows)
% and checks (a) == (b) bit-for-bit, no alarm before the attack, a latch
% after attack onset, and a handback before the end.
% Run from the repo root (MATLAB, or Octave with tools/octave_shim).

prm_boot = kujur_params();
[~, ~, H_all, ~, ~, Phi, Q, scn, z_all, V, P0] = ...
    generate_test_data(prm_boot, 42, 0.10, 0.5, [1 1 1]);
spoofInfo.phiAcc = Phi; spoofInfo.qAcc = Q;
n      = prm_boot.n_states;
mMax   = double(CST_spfParam.MAX_MEAS);
N      = scn.N_total;

rng(7);
nSv    = randi([5 8], 1, N);          % visible satellites per epoch
mk     = 2 * nSv;                     % code + carrier rows

res = cell(1, 2);
for variant = 1:2
    padded = (variant == 1);
    mode = CST_spfMode.NOMINAL;
    sys  = STRUCT_SPF.setSys(mode, STRUCT_SPF.setFilter(zeros(n,1), P0), ...
        STRUCT_SPF.setTrial(zeros(n,1), P0), STRUCT_SPF.zeroMonitorPool, ...
        STRUCT_SPF.setAnchor(false, zeros(n,1), P0, uint32(0)), 0, 0, 0);
    kf_x = zeros(n, 1); kf_P = P0;
    st = zeros(1, N); xn = zeros(n, N); det = []; hb = []; al = false(1, N);
    for k = 1:N
        m = mk(k);
        z = z_all(1:m, k); H = H_all(1:m, :, k); R = V(1:m, 1:m);
        [kf_x, kf_P, y, S, ~, xp] = kalman_update_step(kf_x, kf_P, z, H, spoofInfo, R);
        if padded
            yP = zeros(mMax, 1); yP(1:m) = y;
            SP = zeros(mMax);    SP(1:m, 1:m) = S;
            HP = zeros(mMax, n); HP(1:m, :) = H;
            RP = zeros(mMax);    RP(1:m, 1:m) = R;
        else
            yP = y; SP = S; HP = H; RP = R;
        end
        [sys, tel] = protectedNav(sys, yP, SP, HP, m, xp, kf_x, kf_P, RP, spoofInfo, k);
        if tel.kfCommand.reseedKF
            kf_x = tel.kfCommand.reseedState; kf_P = tel.kfCommand.reseedCov;
        end
        st(k) = tel.info.mode; xn(:, k) = tel.nav.state;
        al(k) = tel.info.ssAlarm || tel.info.cpiAlarm;
        if tel.info.eventLatched,  det(end+1) = k; end %#ok<AGROW>
        if tel.info.eventHandback, hb(end+1)  = k; end %#ok<AGROW>
    end
    res{variant} = struct('st', st, 'xn', xn, 'det', det, 'hb', hb, 'al', al);
end

same   = isequal(res{1}.st, res{2}.st) && isequal(res{1}.xn, res{2}.xn) && ...
         isequal(res{1}.det, res{2}.det) && isequal(res{1}.hb, res{2}.hb);
quiet  = ~any(res{1}.al(1:scn.spoof_start - 1));
latch  = ~isempty(res{1}.det) && res{1}.det(1) >= scn.spoof_start && ...
         res{1}.det(1) <= scn.spoof_start + double(CST_spfParam.WINDOW_LENGTH);
back   = ~isempty(res{1}.hb) && res{1}.hb(end) > scn.spoof_end;
err    = vecnorm(res{1}.xn(1:3, :));
bounded = max(err(scn.spoof_start:end)) < 1.0;

fprintf('numMeas per epoch: min %d max %d\n', min(mk), max(mk));
fprintf('padded == exact          : %s\n', pf(same));
fprintf('no alarm before attack   : %s\n', pf(quiet));
fprintf('latch within N of onset  : %s (epoch %s)\n', pf(latch), mat2str(res{1}.det));
fprintf('handback after attack    : %s (epoch %s)\n', pf(back), mat2str(res{1}.hb));
fprintf('|3-D err| < 1 m post-onset: %s (%.3f m)\n', pf(bounded), max(err(scn.spoof_start:end)));
if same && quiet && latch && back && bounded
    fprintf('=== test_variable_numMeas PASS ===\n');
else
    fprintf('=== test_variable_numMeas FAIL ===\n');
end

function s = pf(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end
