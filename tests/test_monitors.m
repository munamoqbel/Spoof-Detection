%% test_monitors.m
% Unit + Monte-Carlo checks of the deployed monitor engine (supersedes the
% dual_monitor-era test_all.m, whose Test 5 cannot drive SPF_cpiMonitor since
% the window length and thresholds moved into CST_spfParam).
% Run from the repo root (MATLAB, or Octave with tools/octave_shim).

fprintf('=============================================\n');
fprintf('  Tests - Kujur 2024 dual monitor (engine)\n');
fprintf('=============================================\n\n');
pass = 0; fail = 0;

prm_boot = kujur_params();
n_states = prm_boot.n_states;
m_meas   = prm_boot.m_meas;
n_sv     = prm_boot.n_sv;
idx_z    = prm_boot.idx_z;

% ---- shared static geometry for unit tests ----
sigma_code = 0.36; sigma_carrier = 0.003;
rng(1);
los = randn(n_sv,3); los = los ./ vecnorm(los,2,2);
H_t = zeros(m_meas, n_states);
V_diag = zeros(m_meas,1);
for sv = 1:n_sv
    H_t(2*sv-1,1:3) = -los(sv,:);
    H_t(2*sv,  1:3) = -los(sv,:);
    V_diag(2*sv-1)  = sigma_code^2;
    V_diag(2*sv)    = sigma_carrier^2;
end
V   = diag(V_diag);
S_t = H_t*(1e-6*eye(n_states))*H_t' + V;  S_t = (S_t+S_t')/2;
f_t     = H_t(:, idx_z);
s2gu_t  = f_t' * (S_t \ f_t);

% Exact N(0,S) sampler for tests (S != V, so sqrt(V_diag) is NOT valid)
[U_s, D_s] = eig(S_t); D_s = max(diag(D_s),0);
L_s = U_s * diag(sqrt(D_s));

%% Test 1: Eq. 20 positive and carrier-dominated
fprintf('Test 1: sigma2_gamma_u positive & large ...\n');
if s2gu_t > 1.0/(2.0*sigma_carrier)^2
    fprintf('  PASS  s2gu=%.0f\n\n', s2gu_t); pass=pass+1;
else
    fprintf('  FAIL  s2gu=%.0f\n\n', s2gu_t); fail=fail+1;
end

%% Test 2: solve_N_min sane
fprintf('Test 2: solve_N_min ...\n');
sigma_t_test = 0.10;
[N_min_t, T_N_t, k_FA_t, k_MD_t, PMD_t] = solve_N_min( ...
    s2gu_t, sigma_t_test, prm_boot.P_FA_plus, prm_boot.P_FA_minus, ...
    prm_boot.P_MD_plus, prm_boot.P_MD_minus);
if N_min_t>=1 && N_min_t<=100 && PMD_t<=prm_boot.P_MD_plus
    fprintf('  PASS  N_min=%d PMD=%.2e\n\n', N_min_t, PMD_t); pass=pass+1;
else
    fprintf('  FAIL  N_min=%d PMD=%.2e\n\n', N_min_t, PMD_t); fail=fail+1;
end

%% Test 3: kujur_params 0-arg
fprintf('Test 3: kujur_params() ...\n');
try
    p0 = kujur_params();
    assert(p0.n_states>0 && p0.N_min>0);
    fprintf('  PASS\n\n'); pass=pass+1;
catch e
    fprintf('  FAIL  %s\n\n', e.message); fail=fail+1;
end

%% Test 4: kujur_params 4-arg
fprintf('Test 4: kujur_params(4 args) ...\n');
try
    p4 = kujur_params(N_min_t, T_N_t, k_FA_t, k_MD_t);
    assert(p4.N_min==N_min_t && abs(p4.CPI_threshold-T_N_t)<1e-12);
    fprintf('  PASS\n\n'); pass=pass+1;
catch e
    fprintf('  FAIL  %s\n\n', e.message); fail=fail+1;
end

%% Test 5: CPI false-alarm rate under H0 (gamma ~ N(0,S) exactly)
fprintf('Test 5: CPI P_FA under H0 ...\n');
n_mc = 3000; P_FA_t = 1e-3;
T_fa  = gaminv(1-P_FA_t, N_min_t/2, 2);
gw = zeros(m_meas, N_min_t); Sw = repmat(S_t, [1 1 N_min_t]);
Hw = repmat(H_t, [1 1 N_min_t]); mw = uint8(m_meas) * ones(N_min_t, 1, 'uint8');
n_alarm = 0; q_sum = 0;
for k = 1:n_mc
    for j = 1:N_min_t, gw(:,j) = L_s * randn(m_meas,1); end
    [a, q, ~] = SPF_cpiMonitor(gw, Sw, Hw, mw, idx_z, N_min_t, T_fa);
    n_alarm = n_alarm + a; q_sum = q_sum + q;
end
q_mean = q_sum / n_mc;                       % Gamma(N/2,2) mean = N
if n_alarm <= 10 && abs(q_mean - N_min_t) < 0.15*N_min_t
    fprintf('  PASS  alarms=%d/%d (expect ~%d), mean q=%.2f (expect %d)\n\n', ...
        n_alarm, n_mc, round(n_mc*P_FA_t), q_mean, N_min_t); pass=pass+1;
else
    fprintf('  FAIL  alarms=%d/%d, mean q=%.2f\n\n', n_alarm, n_mc, q_mean); fail=fail+1;
end

%% Test 6: CPI detects a 2 cm position-domain bias with the design threshold
fprintf('Test 6: CPI detection under a bias ...\n');
bias = 0.02;                                 % m along the monitored axis
for j = 1:N_min_t, gw(:,j) = L_s * randn(m_meas,1) + f_t * bias; end
[a, q, ~] = SPF_cpiMonitor(gw, Sw, Hw, mw, idx_z, N_min_t, T_N_t);
if a
    fprintf('  PASS  q=%.1f > T_N=%.1f\n\n', q, T_N_t); pass=pass+1;
else
    fprintf('  FAIL  q=%.1f <= T_N=%.1f\n\n', q, T_N_t); fail=fail+1;
end

%% Test 7: SS calibration under H0 (truth driven by the KF's own Q)
% Checks Eq. 50 directly: on a long clean run, over NON-overlapping windows,
% the normalised separation r = q_SS / sigma_SS must be N(0,1):
%   mean(r^2) ~ 1  and  P(|r| > 2) ~ 4.55 %.
% (Counting 1e-3 exceedances over overlapping windows is too clustered to
% be a useful test; that count is printed for information only.)
% NOTE: generate_test_data uses a static truth without process noise, so
% the KF is conservative there and sigma_SS overstates the spread; here
% the truth carries w ~ N(0, Q) so the KF is exactly matched.
fprintf('Test 7: SS calibration under H0 ...\n');
[Phi, Q] = build_Phi_Q(prm_boot);
propTel = STRUCT_SPF.setPropTel(Phi, Q);
sqQ = sqrt(diag(Q)); sqV = sqrt(V_diag);
rng(11);
P = eye(n_states);
for it = 1:200                                % Riccati warm-up
    Pb = Phi*P*Phi' + Q; Sk = H_t*Pb*H_t' + V; L = Pb*H_t'/Sk;
    P = (eye(n_states) - L*H_t)*Pb; P = (P+P')/2;
end
N  = double(CST_spfParam.WINDOW_LENGTH);
N_run = 600 * N;                              % 600 independent windows
x_true = zeros(n_states,1); x_hat = zeros(n_states,1);
r_all = zeros(1, 3 * (N-1) * 600); n_r = 0;
n_1e3 = 0; k_1e3 = norminv(1 - 1e-3/2);
dS = zeros(n_states, 1); cP = P; x_prev = x_hat; P_prev = P;
for k = 1:N_run
    x_true = Phi*x_true + sqQ .* randn(n_states,1);
    z  = H_t*x_true + sqV .* randn(m_meas,1);
    P_prev = P;
    xb = Phi*x_hat; Pb = Phi*P*Phi' + Q;
    Sk = H_t*Pb*H_t' + V; Sk = (Sk+Sk')/2; L = Pb*H_t'/Sk;
    x_hat = xb + L*(z - H_t*xb);
    P = (eye(n_states) - L*H_t)*Pb; P = (P+P')/2;
    inc = x_hat - Phi * x_prev;               % K*y (increment form) of this epoch
    if mod(k, N) == 1                         % open a fresh window on this epoch
        dS = zeros(n_states, 1); cP = P;
    else
        [dS, cP, r] = SPF_ssMonitor(dS, cP, inc, P, propTel, 2.0, k_MD_t);   % gate at 2 sigma
        rr = r.separation ./ max(r.sigmaSeparation, 1e-12);
        r_all(n_r+1 : n_r+3) = rr; n_r = n_r + 3;
        n_1e3 = n_1e3 + sum(abs(rr) > k_1e3);
    end
    x_prev = x_hat;
end
r_all = r_all(1:n_r);
ms  = mean(r_all.^2);
p2  = mean(abs(r_all) > 2);
ok7 = abs(ms - 1) < 0.15 && p2 > 0.025 && p2 < 0.07;
fprintf('  mean(r^2) = %.3f (expect 1) | P(|r|>2) = %.3f (expect 0.046) | 1e-3 exceedances = %d of %d (expect ~%d, info)\n', ...
    ms, p2, n_1e3, n_r, round(n_r*1e-3));
if ok7
    fprintf('  PASS\n\n'); pass=pass+1;
else
    fprintf('  FAIL\n\n'); fail=fail+1;
end

%% Test 8: SS detects a 1 m KF-vs-coast separation with the design gate
fprintf('Test 8: SS detection under a separation ...\n');
% last epoch of the Test-7 run: window opened at k-1 (cov P_prev), tested at
% k with that epoch's real increment plus 1 m on axis 3
dS = zeros(n_states, 1); cP = P_prev;
inc(3) = inc(3) + 1.0;
[~, ~, r] = SPF_ssMonitor(dS, cP, inc, P, propTel, CST_spfParam.K_FALSE_ALERT, CST_spfParam.K_MISSED_DETECTION);
if r.anyAlarm && r.alarmPerAxis(3) && ~r.alarmPerAxis(1)
    fprintf('  PASS  axis-3 alarm, PL=%.3f m\n\n', r.maxProtectionLevel); pass=pass+1;
else
    fprintf('  FAIL  alarmPerAxis=%s\n\n', mat2str(r.alarmPerAxis)); fail=fail+1;
end

%% Test 9: SPF_revalidation passes on consistent GNSS, fails on a 2 m offset
fprintf('Test 9: SPF_revalidation ...\n');
cP = P; cmp = zeros(n_states, 1);      % coast - prior = 0 (host not updating)
S_r = H_t * cP * H_t' + V; S_r = (S_r + S_r')/2;
[U_r, D_r] = eig(S_r); L_r = U_r * diag(sqrt(max(diag(D_r), 0)));
n_ok = 0; n_tr = 500;
for k = 1:n_tr
    y = L_r * randn(m_meas, 1);              % innovation vs the coast (prior = coast)
    ok = SPF_revalidation(y, H_t, V, m_meas, cmp, cP);
    n_ok = n_ok + ok;
end
y_off = L_r * randn(m_meas, 1) + H_t(:, 3) * 2.0;
[ok_off, q_off] = SPF_revalidation(y_off, H_t, V, m_meas, cmp, cP);
if n_ok >= 0.99 * n_tr && ~ok_off
    fprintf('  PASS  pass rate %.3f, offset q=%.1f rejected\n\n', n_ok/n_tr, q_off); pass=pass+1;
else
    fprintf('  FAIL  pass rate %.3f, offset ok=%d\n\n', n_ok/n_tr, ok_off); fail=fail+1;
end

%% Summary
fprintf('=============================================\n');
fprintf('  %d passed, %d failed\n', pass, fail);
if fail == 0, fprintf('  === ALL MONITOR TESTS PASS ===\n'); end
fprintf('=============================================\n');
