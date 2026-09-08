%% test_variable_numMeas.m
% The gate must accept a measurement count that changes every epoch
% (satellites rise/set). Runs the host contract (docs/HOST_2HZ_WIRING.m)
% on the harness scenario with 5..8 satellites visible per epoch, two ways:
%   (a) setKfMeas with S from the host, arrays padded to MAX_MEAS (Coder style)
%   (b) kfMeasFromUpdate with y, H, R, P_bar at exact size (S formed by the gate)
% and checks (a) == (b) bit-for-bit, no alarm before the attack, a latch
% after attack onset, and a handback before the end.
% Run from the repo root (MATLAB, or Octave with tools/octave_shim).

prm_boot = kujur_params();
[~, ~, H_all, ~, ~, Phi, Q, scn, z_all, V, P0] = ...
    generate_test_data(prm_boot, 42, 0.10, 0.5, [1 1 1]);
propTel = STRUCT_SPF.setPropTel(Phi, Q);
n      = prm_boot.n_states;
mMax   = double(CST_spfParam.MAX_MEAS);
N      = scn.N_total;

rng(7);
nSv    = randi([5 8], 1, N);          % visible satellites per epoch
mk     = 2 * nSv;                     % code + carrier rows

res = cell(1, 2);
for variant = 1:2
    padded = (variant == 1);
    sys = STRUCT_SPF.zeroSys; sys.coastCov = P0;
    sys.anchor = STRUCT_SPF.setAnchor(true, zeros(n,1), P0, uint32(0));
    kf_x = zeros(n, 1); kf_P = P0; tr_x = kf_x; tr_P = P0; mode = CST_spfMode.NOMINAL;
    st = zeros(1, N); xn = zeros(n, N); det = []; hb = []; al = false(1, N);
    for k = 1:N
        m = mk(k);
        z = z_all(1:m, k); H = H_all(1:m, :, k); R = V(1:m, 1:m);
        inProb = (mode == CST_spfMode.PROBATION); kfUpd = (mode == CST_spfMode.NOMINAL);
        if inProb
            [tr_x, tr_P, y, S, ~, xp, xpP] = kalman_update_step(tr_x, tr_P, z, H, propTel, R);
            post = tr_x; postP = tr_P;
        else
            [ax, aP, y, S, ~, xp, xpP] = kalman_update_step(kf_x, kf_P, z, H, propTel, R);
            post = ax; postP = aP;
        end
        if padded      % (a) host supplies S, arrays padded to MAX_MEAS
            yP = zeros(mMax, 1); yP(1:m) = y;
            SP = zeros(mMax);    SP(1:m, 1:m) = S;
            HP = zeros(mMax, n); HP(1:m, :) = H;
            RP = zeros(mMax);    RP(1:m, 1:m) = R;
            kfMeas = STRUCT_SPF.setKfMeas(yP, SP, HP, RP, m, xp, post, postP);
        else           % (b) host supplies y, H, R, P_bar at exact size; gate forms S
            kfMeas = STRUCT_SPF.kfMeasFromUpdate(y, H, R, m, xp, xpP, post, postP);
        end
        [sys, tel] = protectedNav(sys, kfMeas, propTel, k);
        mode = tel.info.mode;
        if kfUpd, kf_x = post; kf_P = postP; else, [kf_x, kf_P] = insCoast(kf_x, kf_P, propTel); end
        if tel.nav.applyCorrection, kf_x = tel.nav.state; kf_P = tel.nav.covar; end
        if tel.kfCommand.reseedKF, tr_x = kf_x; tr_P = kf_P; end
        st(k) = tel.info.mode; xn(:, k) = kf_x;
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
