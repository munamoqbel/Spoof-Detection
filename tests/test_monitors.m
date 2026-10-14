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
% The engine reads N and T_N from CST_spfParam, so the test runs on the
% deployed window (N = WINDOW_LENGTH) and applies its own Monte-Carlo
% threshold (P_FA = 1e-3 of the same Gamma(N/2, 2) law) to the returned q.
fprintf('Test 5: CPI P_FA under H0 ...\n');
n_mc = 3000; P_FA_t = 1e-3;
N_w   = double(CST_spfParam.WINDOW_LENGTH);
T_fa  = gaminv(1-P_FA_t, N_w/2, 2);
gw = zeros(m_meas, N_w); Sw = repmat(S_t, [1 1 N_w]);
Hw = repmat(H_t, [1 1 N_w]); mw = uint8(m_meas) * ones(N_w, 1, 'uint8');
n_alarm = 0; q_sum = 0;
for k = 1:n_mc
    for j = 1:N_w, gw(:,j) = L_s * randn(m_meas,1); end
    [~, q, ~] = SPF_cpiMonitor(gw, Sw, Hw, mw, idx_z, 0.0);
    n_alarm = n_alarm + (q > T_fa); q_sum = q_sum + q;
end
q_mean = q_sum / n_mc;                       % Gamma(N/2,2) mean = N
if n_alarm <= 10 && abs(q_mean - N_w) < 0.15*N_w
    fprintf('  PASS  alarms=%d/%d (expect ~%d), mean q=%.2f (expect %d)\n\n', ...
        n_alarm, n_mc, round(n_mc*P_FA_t), q_mean, N_w); pass=pass+1;
else
    fprintf('  FAIL  alarms=%d/%d, mean q=%.2f\n\n', n_alarm, n_mc, q_mean); fail=fail+1;
end

%% Test 6: CPI detects a 2 cm position-domain bias with the design threshold
fprintf('Test 6: CPI detection under a bias ...\n');
bias = 0.02;                                 % m along the monitored axis
for j = 1:N_w, gw(:,j) = L_s * randn(m_meas,1) + f_t * bias; end
[a, q, ~] = SPF_cpiMonitor(gw, Sw, Hw, mw, idx_z, 0.0);
if a
    fprintf('  PASS  q=%.1f > T_N=%.1f\n\n', q, CST_spfParam.CPI_THRESHOLD); pass=pass+1;
else
    fprintf('  FAIL  q=%.1f <= T_N=%.1f\n\n', q, CST_spfParam.CPI_THRESHOLD); fail=fail+1;
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
        [dS, cP, r] = SPF_ssMonitor(dS, cP, inc, P, propTel, zeros(1, 3));
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
[~, ~, r] = SPF_ssMonitor(dS, cP, inc, P, propTel, zeros(1, 3));
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

%% Test 10: nominal bias allowance B (CST_spfParam.NOMINAL_BIAS, passed in by the pool)
% SS: the gate moves from k_FA*sigma_SS to k_FA*sigma_SS + B and the protection
% level by exactly B; a separation past the plain gate but inside the raised one
% is no alarm. CPI: the position-domain innovation counts only beyond B, so the
% 2 cm bias of Test 6 is invisible under B = 5 cm, reduced under B = 1 cm, and
% xi follows max(|xi0| - B*sqrt(s2), 0) epoch by epoch; B = 0 is Test 6 itself.
fprintf('Test 10: nominal bias allowance ...\n');
B = 1.0;
dS = zeros(n_states, 1); cP = P_prev;
inc0 = inc; inc0(3) = 0.0;                     % window just opened: the axis-3 separation is inc0(3)
[~, ~, rRef] = SPF_ssMonitor(dS, cP, inc0, P, propTel, zeros(1, 3));
thr0 = CST_spfParam.K_FALSE_ALERT * rRef.sigmaSeparation(3);
incA = inc0; incA(3) = thr0 + 0.5 * B;         % past the plain gate, inside the raised one
[~, ~, rA0] = SPF_ssMonitor(dS, cP, incA, P, propTel, zeros(1, 3));
[~, ~, rAB] = SPF_ssMonitor(dS, cP, incA, P, propTel, [0 0 B]);
incC = inc0; incC(3) = thr0 + 1.5 * B;         % past the raised gate
[~, ~, rCB] = SPF_ssMonitor(dS, cP, incC, P, propTel, [0 0 B]);
okSS = rA0.alarmPerAxis(3) && ~rAB.alarmPerAxis(3) && rCB.alarmPerAxis(3) ...
    && abs((rAB.protectionLevel(3) - rA0.protectionLevel(3)) - B) < 1e-9 ...
    && abs(rAB.separation(3) - rA0.separation(3)) < 1e-12;
[a0, q0, xi0] = SPF_cpiMonitor(gw, Sw, Hw, mw, idx_z, 0.0);        % the Test-6 window
[a5, q5, ~]   = SPF_cpiMonitor(gw, Sw, Hw, mw, idx_z, 0.05);
[a1, q1, xi1] = SPF_cpiMonitor(gw, Sw, Hw, mw, idx_z, 0.01);
xiExp = max(abs(xi0) - 0.01 * sqrt(s2gu_t), 0);
okCPI = a0 && (q0 == q) && ~a5 && (q5 == 0) && a1 && (q1 > 0) && (q1 < q0) ...
    && max(abs(xi1 - xiExp)) < 1e-6 * max(1, max(abs(xi0)));
if okSS && okCPI
    fprintf('  PASS  SS gate %.3f -> %.3f m, PL + %.3f m; CPI q %.1f (B=0), %.1f (B=1 cm), %.1f (B=5 cm)\n\n', ...
        thr0, thr0 + B, rAB.protectionLevel(3) - rA0.protectionLevel(3), q0, q1, q5); pass=pass+1;
else
    fprintf('  FAIL  SS ok=%d (alarms %d %d %d), CPI ok=%d (q %.1f %.1f %.1f)\n\n', okSS, ...
        rA0.alarmPerAxis(3), rAB.alarmPerAxis(3), rCB.alarmPerAxis(3), okCPI, q0, q1, q5); fail=fail+1;
end

%% Summary
fprintf('=============================================\n');
fprintf('  %d passed, %d failed\n', pass, fail);
if fail == 0, fprintf('  === ALL MONITOR TESTS PASS ===\n'); end
fprintf('=============================================\n');
