%% test_kfMeasBoth.m
% The gate reads the S and R of setKfMeas / setKfMeasBoth (innovationCov, measNoiseCov) and never the host's pair
% (innovationCovHost, measNoiseCovHost). On the harness scenario (docs/HOST_2HZ_WIRING.m contract, 5..8 satellites):
%   (a) setKfMeas with the nominal S and R (the reference)
%   (b) setKfMeasBoth with the same nominal pair and a host pair scaled nine-fold: bit-for-bit (a)
%   (c) setKfMeas handed the scaled pair, the old interface of a host that scales R: the CPI margin shrinks and the
%       latch comes no earlier
% and the four constructors give one struct layout, kfMeasFromUpdateBoth forming both S from the same P_bar.
% Run from the repo root (MATLAB, or Octave with tools/octave_shim).

prm_boot = kujur_params();
[~, ~, H_all, ~, ~, Phi, Q, scn, z_all, V, P0] = ...
    generate_test_data(prm_boot, 42, 0.10, 0.5, [1 1 1]);
propTel = STRUCT_SPF.setPropTel(Phi, Q);
n      = prm_boot.n_states;
mMax   = double(CST_spfParam.MAX_MEAS);
N      = scn.N_total;
SCALE  = 9;

rng(7);
nSv    = randi([5 8], 1, N);
mk     = 2 * nSv;

res = cell(1, 3);
for variant = 1:3
    sys = STRUCT_SPF.zeroSys; sys.coastCov = P0;
    sys.anchor = STRUCT_SPF.setAnchor(true, zeros(n,1), P0, uint32(0));
    kf_x = zeros(n, 1); kf_P = P0; tr_x = kf_x; tr_P = P0; mode = CST_spfMode.NOMINAL;
    st = zeros(1, N); xn = zeros(n, N); det = []; hb = []; cpi = zeros(1, N);
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
        yP = zeros(mMax, 1); yP(1:m) = y;
        HP = zeros(mMax, n); HP(1:m, :) = H;
        SP = zeros(mMax);    SP(1:m, 1:m) = S;
        RP = zeros(mMax);    RP(1:m, 1:m) = R;
        SH = zeros(mMax);    SH(1:m, 1:m) = H * xpP * H' + SCALE * R;   % the host's scaled pair
        RH = zeros(mMax);    RH(1:m, 1:m) = SCALE * R;
        switch variant
            case 1, kfMeas = STRUCT_SPF.setKfMeas(yP, SP, HP, RP, m, xp, post, postP);
            case 2, kfMeas = STRUCT_SPF.setKfMeasBoth(yP, SP, HP, RP, m, xp, post, postP, SH, RH);
            case 3, kfMeas = STRUCT_SPF.setKfMeas(yP, SH, HP, RH, m, xp, post, postP);
        end
        [sys, tel] = SPF_protectedNav(sys, kfMeas, propTel, k);
        mode = tel.info.mode;
        if kfUpd, kf_x = post; kf_P = postP; else, [kf_x, kf_P] = SPF_insCoast(kf_x, kf_P, propTel); end
        if tel.nav.applyCorrection, kf_x = tel.nav.state; kf_P = tel.nav.covar; end
        if tel.kfCommand.startTrial, tr_x = kf_x; tr_P = kf_P; end
        st(k) = tel.info.mode; xn(:, k) = kf_x; cpi(k) = max(tel.info.cpiRatio);
        if tel.info.eventLatched,  det(end+1) = k; end %#ok<AGROW>
        if tel.info.eventHandback, hb(end+1)  = k; end %#ok<AGROW>
    end
    res{variant} = struct('st', st, 'xn', xn, 'det', det, 'hb', hb, 'cpi', cpi);
end

same = isequal(res{1}.st, res{2}.st) && isequal(res{1}.xn, res{2}.xn) && isequal(res{1}.det, res{2}.det) ...
    && isequal(res{1}.hb, res{2}.hb) && isequal(res{1}.cpi, res{2}.cpi);
cpiA = max(res{1}.cpi(scn.spoof_start:scn.spoof_end));  cpiC = max(res{3}.cpi(scn.spoof_start:scn.spoof_end));
shrunk = cpiC < 0.5 * cpiA;
noEarlier = isempty(res{3}.det) || (~isempty(res{1}.det) && res{3}.det(1) >= res{1}.det(1));
latched = ~isempty(res{1}.det);

% one layout from every constructor, and kfMeasFromUpdateBoth's two S from the same P_bar
m = 6; y = randn(m, 1); H = randn(m, n); R = diag(1 + rand(m, 1)); Rh = 4 * R; Pb = eye(n) * 2;
kA = STRUCT_SPF.setKfMeas(zeros(mMax, 1), zeros(mMax), zeros(mMax, n), zeros(mMax), 0, zeros(n, 1), zeros(n, 1), eye(n));
kB = STRUCT_SPF.setKfMeasBoth(zeros(mMax, 1), zeros(mMax), zeros(mMax, n), zeros(mMax), 0, zeros(n, 1), zeros(n, 1), eye(n), zeros(mMax), zeros(mMax));
kC = STRUCT_SPF.kfMeasFromUpdate(y, H, R, m, zeros(n, 1), Pb, zeros(n, 1), eye(n));
kD = STRUCT_SPF.kfMeasFromUpdateBoth(y, H, R, Rh, m, zeros(n, 1), Pb, zeros(n, 1), eye(n));
kZ = STRUCT_SPF.zeroKfMeas;
layout = isequal(fieldnames(kA), fieldnames(kB)) && isequal(fieldnames(kA), fieldnames(kC)) ...
    && isequal(fieldnames(kA), fieldnames(kD)) && isequal(fieldnames(kA), fieldnames(kZ));
SD = H * Pb * H' + R;  SDh = H * Pb * H' + Rh;
formed = max(max(abs(kD.innovationCov(1:m, 1:m) - SD))) < 1e-12 && max(max(abs(kD.innovationCovHost(1:m, 1:m) - SDh))) < 1e-12 ...
    && isequal(kD.measNoiseCov(1:m, 1:m), R) && isequal(kD.measNoiseCovHost(1:m, 1:m), Rh) && ~any(any(kD.innovationCov(m+1:end, :))) ...
    && isequal(kC.innovationCovHost, kC.innovationCov) && isequal(kC.measNoiseCovHost, kC.measNoiseCov) && isequal(kD.numMeas, uint8(m));

fprintf('nominal pair with the host''s beside == nominal pair alone (bitwise): %s\n', pf(same));
fprintf('the gate handed the scaled pair: CPI margin %.3f against %.3f, latch no earlier: %s\n', cpiC, cpiA, pf(shrunk && noEarlier));
fprintf('one layout from every constructor, both S formed from P_bar: %s\n', pf(layout && formed));
if same && shrunk && noEarlier && latched && layout && formed
    fprintf('=== test_kfMeasBoth PASS ===\n');
else
    fprintf('=== test_kfMeasBoth FAIL ===\n');
end

function s = pf(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end
