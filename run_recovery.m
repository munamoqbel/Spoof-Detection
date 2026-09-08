%% run_recovery.m
% Full paper architecture demo, 3-AXIS MONITOR BANK:
%   per-axis detection (u = e1, e2, e3) -> LATCH on any axis ->
%   fall back to last clean window's coast -> INS-only ->
%   GNSS re-validation -> monitored PROBATION -> handback.
% Compared against the unprotected KF on the identical measurements.
%
% OFFLINE DESIGN (per-axis, then common constants)
%   - False-alert budgets split across the 3 axes (Bonferroni):
%       P_FA_plus/3 per-axis CPI, P_FA_minus/3 per-axis SS.
%     Missed-detection budgets are per fault hypothesis: NOT split.
%   - solve_N_min per axis -> N_min(a); run ALL axes at the common
%     window length N_c = max(N_min): longer-than-solved N only lowers
%     P_MD, and T_N / k_FA are recomputed at N_c with the split budgets.
%
% Run test_all.m first (all existing tests must still pass).

clear; clc; close all;

% ---------------- scenario toggles ----------------
seed           = 42;
sigma_t_attack = 0.10;    % injected tracking-error std (m). 0 => ramp-only
ramp_rate      = 0.5;     % injected ramp (m/epoch).         0 => tracking-only
attack_dir     = [1 1 1]; % attack direction: [0 0 1]=Down, [1 0 0]=axis1,
                          % [1 1 1]=oblique diagonal (auto-normalised)
alpha_reval    = 1e-3;    % re-validation false-alarm allocation per epoch
M_dwell        = 10;      % consecutive passes to open probation
axes_mon       = [1 2 3]; % monitored axes (the bank)

fprintf('=============================================\n');
fprintf('  Kujur 2024 - Sec. 5 Recovery, 3-axis bank\n');
fprintf('=============================================\n\n');

%% Step 1 - data (one measurement stream feeds both systems)
prm_boot = kujur_params();

[~, S_all, H_all, xh_base, ~, Phi, Q, scn, z_all, V, P0] = ...
    generate_test_data(prm_boot, seed, sigma_t_attack, ramp_rate, attack_dir);

spoofInfo.phiAcc = Phi;
spoofInfo.qAcc   = Q;

N  = scn.N_total;
fs = prm_boot.fs;
t  = (0:N-1) / fs;

fprintf('Attack: dir [%.2f %.2f %.2f] | sigma_t_attack = %.2f m | ramp = %.2f m/epoch\n', ...
    scn.attack_dir(1), scn.attack_dir(2), scn.attack_dir(3), ...
    scn.sigma_t_attack, scn.ramp_rate);
fprintf('Phase 1 clean 1..%d | Phase 2 spoofed %d..%d (ramp from %d) | Phase 3 clean %d..%d\n\n', ...
    scn.clean_end, scn.spoof_start, scn.spoof_end, scn.ramp_start, scn.clean2_start, N);

%% Step 2 - offline design, per axis (DESIGN sigma_t, never the attack value)
n_ax    = numel(axes_mon);
PFAp_ax = prm_boot.P_FA_plus  / n_ax;     % CPI budget per axis
PFAm_ax = prm_boot.P_FA_minus / n_ax;     % SS  budget per axis

s2min = zeros(1, 3);
for a = axes_mon
    s2 = zeros(N, 1);
    for kk = 1:N
        f      = H_all(:, a, kk);
        s2(kk) = f' * (S_all(:,:,kk) \ f);
    end
    s2min(a) = min(s2);                   % worst geometry per axis
end

N_ax = zeros(1, 3);
for a = axes_mon
    fprintf('=== axis %d design ===\n', a);
    [N_ax(a), ~, ~, ~, ~] = solve_N_min(s2min(a), scn.sigma_t, ...
        PFAp_ax, PFAm_ax, prm_boot.P_MD_plus, prm_boot.P_MD_minus);
end

N_c  = max(N_ax(axes_mon));                       % common window length
T_Nc = gaminv(1 - PFAp_ax, N_c/2, 2);             % Eq. 35 at N_c
kFAc = norminv(1 - (PFAm_ax / N_c) / 2);          % SS gate at N_c
kMDc = norminv(1 - prm_boot.P_MD_minus);

prm          = kujur_params(N_c, T_Nc, kFAc, kMDc);
prm.mon_axes = axes_mon;

fprintf('--- Bank constants ---\n');
fprintf('Per-axis N_min = [%d %d %d] -> common N = %d (%.1f s)\n', ...
    N_ax(1), N_ax(2), N_ax(3), N_c, N_c/fs);
fprintf('T_N = %.4f | k_FA = %.4f | k_MD = %.4f  (FA budgets /%d per axis)\n', ...
    T_Nc, kFAc, kMDc, n_ax);
for a = axes_mon
    fprintf('axis %d: PMD achieved at N=%d : %.3e (req %.1e)\n', a, N_c, ...
        compute_PMD_eq38(N_c, s2min(a), scn.sigma_t, PFAp_ax), prm_boot.P_MD_plus);
end
fprintf('\n');

% re-validation constants (offline; chi2inv(p, m) == gaminv(p, m/2, 2))
rec.T_reval = gaminv(1 - alpha_reval, prm_boot.m_meas / 2, 2);
rec.M_dwell = M_dwell;
rec.M_prob  = prm.N_min + 8;
fprintf('Re-validation: T_reval = %.2f (chi2_%d), dwell = %d ep, probation = %d ep\n\n', ...
    rec.T_reval, prm_boot.m_meas, rec.M_dwell, rec.M_prob);

% ---- consistency guard: the runtime reads CST_spfParam, not prm/rec ----
% The monitors (monitorPool/cpiMonitor/ssMonitor/protectedNav) take their
% constants from CST_spfParam.m (Coder-friendly). The offline design above
% is informational; warn if the two have drifted apart.
chk = {'WINDOW_LENGTH',        double(CST_spfParam.WINDOW_LENGTH),        N_c;
       'CPI_THRESHOLD',        CST_spfParam.CPI_THRESHOLD,                T_Nc;
       'K_FALSE_ALERT',        CST_spfParam.K_FALSE_ALERT,                kFAc;
       'K_MISSED_DETECTION',   CST_spfParam.K_MISSED_DETECTION,           kMDc;
       'REVAL_THRESHOLD',      CST_spfParam.REVAL_THRESHOLD,              rec.T_reval;
       'REVAL_DWELL_REQUIRED', double(CST_spfParam.REVAL_DWELL_REQUIRED), rec.M_dwell;
       'PROBATION_LENGTH',     double(CST_spfParam.PROBATION_LENGTH),     rec.M_prob};
for c = 1:size(chk, 1)
    if abs(chk{c,2} - chk{c,3}) > 1e-9 * max(1, abs(chk{c,3}))
        warning('run_recovery:designMismatch', ...
            'CST_spfParam.%s = %.10g but offline design gives %.10g (runtime uses CST_spfParam)', ...
            chk{c,1}, chk{c,2}, chk{c,3});
    end
end

%% Step 3 - run the protected system
use_fed = true;
out = recovery_nav_sim(z_all, H_all, V, spoofInfo, P0, prm, use_fed);

%% Step 4 - report (3-D position error norm; panel 1 shows dominant axis)
ia       = scn.idx_attack;                        % dominant attack axis
err_base = vecnorm(xh_base(1:3, :));              % |3-D position error|
err_rec  = vecnorm(out.x_nav(1:3, :));
sig_rss  = vecnorm(out.sig_pos);                  % RSS 1-sigma of position

fprintf('--- Events ---\n');
if isempty(out.ev.t_detect)
    fprintf('No alarm raised (attack below detection floor?) - stayed NOMINAL.\n');
else
    pre  = 1:out.ev.t_detect(1) - 1;
    dmax = max(max(abs(out.x_nav(:, pre) - xh_base(:, pre))));
    fprintf('Pre-detection consistency (recovery KF vs baseline KF): %.2e m\n', dmax);

    for i = 1:numel(out.ev.t_detect)
        td = out.ev.t_detect(i);
        ax = find(out.alarm_axis(td, :));
        fprintf('Latch %d : alarm at epoch %d (t = %.1f s), axis/axes [%s]\n', ...
            i, td, td/fs, num2str(ax));
        ta = out.ev.t_anchor(i);
        if isnan(ta) || ta == 0
            fprintf('          NO clean anchor available - froze current state\n');
        else
            fprintf('          anchored to clean window closed at epoch %d (t = %.1f s)\n', ta, ta/fs);
        end
        if numel(out.ev.t_handback) >= i
            th = out.ev.t_handback(i);
            fprintf('          handback at epoch %d (t = %.1f s); total = %.1f s\n', ...
                th, th/fs, (th - td)/fs);
        else
            fprintf('          still in COAST/PROBATION at end of run\n');
        end
    end
end

span = scn.spoof_start:N;
fprintf('\nMax |3-D position error| from attack onset to end of run:\n');
fprintf('  unprotected KF : %8.3f m\n', max(err_base(span)));
fprintf('  recovery system: %8.3f m\n', max(err_rec(span)));
fprintf('Probation attempts: %d started, %d vetoed, %d committed\n', ...
    numel(out.ev.t_prob_start), numel(out.ev.t_prob_fail), numel(out.ev.t_handback));
n_viol = sum(err_rec > max(3*sig_rss, 1e-9));
fprintf('Integrity check: epochs with |3-D err| > 3-sigma (RSS) claim: %d\n\n', n_viol);

%% Step 5 - plots (attacked axis shown)
t_sp = scn.spoof_start / fs;
t_rp = scn.ramp_start  / fs;
t_se = scn.spoof_end   / fs;

figure('Name','Sec-5 Recovery, 3-axis bank','NumberTitle','off','Position',[50 50 1050 850]);

% Panel 1: navigation output on the attacked axis
subplot(4,1,1); hold on;
ymx = max([max(err_base), max(err_rec), 1]) * 1.1;
fill([t_sp t_se t_se t_sp], [-ymx -ymx ymx ymx], [1 .85 .85], ...
    'EdgeColor','none','FaceAlpha',.3,'HandleVisibility','off');
plot(t, xh_base(ia,:), 'r-', 'LineWidth', 1.2, 'DisplayName', 'Unprotected KF (captured)');
plot(t, out.x_nav(ia,:), 'b-', 'LineWidth', 1.4, 'DisplayName', 'Recovery system output');
yline(0, 'k:', 'HandleVisibility','off');
xline(t_rp, 'm:', 'HandleVisibility','off');
ylabel(sprintf('Axis-%d pos (m)', ia));
legend('Location','northwest'); grid on; xlim([0 t(end)]); ylim([-ymx ymx]);
title(sprintf('Navigation output vs truth (truth = 0), attack dir [%.2f %.2f %.2f], dominant axis %d', ...
    scn.attack_dir(1), scn.attack_dir(2), scn.attack_dir(3), ia));

% Panel 2: |3-D position error| on log scale + RSS 3-sigma
subplot(4,1,2);
semilogy(t, max(err_base, 1e-6), 'r-', 'LineWidth', 1.2); hold on;
semilogy(t, max(err_rec,  1e-6), 'b-', 'LineWidth', 1.4);
semilogy(t, max(3*sig_rss, 1e-6), 'b:', 'LineWidth', 1.0);
xline(t_sp,'r:','HandleVisibility','off'); xline(t_se,'r:','HandleVisibility','off');
ylabel('|3-D pos err| (m)'); grid on; xlim([0 t(end)]);
legend('Unprotected','Recovery','Recovery 3\sigma (RSS)','Location','southeast');

% Panel 3: re-validation statistic
subplot(4,1,3);
q_plot = out.q_reval;
q_plot(~out.reval_computed) = NaN;
semilogy(t, q_plot, 'k.-', 'LineWidth', 0.8, 'MarkerSize', 6); hold on;
yline(rec.T_reval, 'r-', 'T_{reval}', 'LineWidth', 1.2);
xline(t_se, 'r:', 'HandleVisibility','off');
ylabel('q_{reval}'); grid on; xlim([0 t(end)]);
title(sprintf('GNSS-vs-coast re-validation, \\chi^2_{%d}; needs %d consecutive passes', ...
    prm_boot.m_meas, rec.M_dwell));

% Panel 4: FSM state + per-axis alarm markers
subplot(4,1,4); hold on;
fill([t_sp t_se t_se t_sp], [0.8 0.8 3.4 3.4], [1 .85 .85], ...
    'EdgeColor','none','FaceAlpha',.3,'HandleVisibility','off');
stairs(t, out.state, 'k-', 'LineWidth', 1.6, 'HandleVisibility','off');
mk = {'r^', 'gs', 'bo'};
for a = 1:3
    iax = find(out.alarm_axis(:, a));
    if ~isempty(iax)
        scatter(t(iax), out.state(iax) + 0.10 + 0.08*a, 16, mk{a}, 'filled', ...
            'DisplayName', sprintf('alarm axis %d', a));
    end
end
xline(t_sp,'r:','HandleVisibility','off');
xline(t_rp,'m:','HandleVisibility','off');
xline(t_se,'r:','HandleVisibility','off');
ylim([0.8 3.4]); yticks([1 2 3]); yticklabels({'NOMINAL','COAST','PROBATION'});
xlabel('Time (s)'); grid on; xlim([0 t(end)]);
legend('Location','northeast');
title('FSM state (per-axis alarm markers)');

figure; plot(t, out.SS_PL, 'm-'); grid on; ylabel('SS PL(m)'); xlabel('Time(s)'); title('Eq.51 protection level (max over live windows/axes)');
%% set anchor age limit
alertLimt = 10;
P=P0;
for k =1:10
    P = spoofInfo.phiAcc * P * (spoofInfo.phiAcc)' + spoofInfo.qAcc;
end
horiz = zeros(600, 1);
for k=1:600
    P = spoofInfo.phiAcc * P * (spoofInfo.phiAcc)' + spoofInfo.qAcc;
    horiz(k) = 3*sqrt(P(1,1) + P(2,2));
end
maxAnchorAge = find(horiz < alertLimt, 1, 'last');
fprintf('maxAnchorAge = %d epochs (%.0f s)\n', maxAnchorAge, maxAnchorAge/2);
