function SPF_plotLog(spfLog, fs)
%SPF_PLOTLOG  Standard diagnostic figures for a gate log (SPF_logInit /
%   SPF_logAppend). fs = 2 Hz epoch rate (default 2).
%
%   Figure 1  mode timeline with events (grey = warm-up, monitors not
%             armed); numMeas
%   Figure 2  monitor margins: SS ratio and CPI ratio per axis (alarm > 1),
%             protection level, host NIS, innovation mean / spread
%   Figure 3  re-validation: qReval vs its threshold, dwell count, coast time
%   Figure 4  fault flags
%
%   How to read it (jamming / false-latch investigation):
%   - a latch inside or right after the grey warm-up band: the filter had
%     not settled (few rows, large PL). Raise ARM_EPOCHS or lower
%     ARM_PL_MAX in CST_spfParam so arming waits for convergence.
%   - NIS >> 1 while jammed and the SS/CPI ratios cross 1 at the same time:
%     the KF's R is optimistic under jamming; the monitors are reporting a
%     real model inconsistency. Fix R (or exclude the rows) on the host side
%     before touching the gate thresholds.
%   - qReval stays above its threshold long after the jamming ends, and
%     grows with coast time: coast capture (P_C too small for the real INS
%     drift). Needs the coast-time budget / Q, not the thresholds.
%   - probation opens and is vetoed repeatedly: look at which ratio vetoes
%     (SS = the trial drifts from the coast, CPI = biased innovations).
%   - innovation mean vs spread (figure 2, last panel): a large MEAN with a
%     small spread is an error common to every row, i.e. the receiver clock
%     bias; compare it with the filter's clock sigma. A large SPREAD is a
%     position error (different projection on every line of sight). At
%     re-acquisition after a gap the mean says whether the clock state was
%     re-initialised or dragged the filter for a minute.
%   - covariance honesty (host Q / initial P): at the end of a pure-INS
%     coast compare the true position error (e.g. the jump at
%     re-acquisition or at a host reset) with sigma_N/E/D in the PL panel.
%     A consistent filter gives a ratio of about 1 to 3; 20 or more means
%     Q (or the post-reset P of the bias states) is too small and the SS
%     and re-validation tests will fire on legitimate corrections.

if nargin < 2, fs = 2; end
n = spfLog.count;
if n > spfLog.numEpochs, n = spfLog.numEpochs; end
t = (0:n-1)' / fs;
ev = @(idx) t(idx(idx <= n));

%% 1. mode and measurement count
figure('Name', 'SPF gate: mode');
subplot(2,1,1);
unarmed = ~spfLog.armed(1:n);
d = diff([false; unarmed; false]);               % one grey band per contiguous unarmed run
bandStart = find(d == 1); bandEnd = find(d == -1) - 1;
for b = 1:numel(bandStart)
    area([t(bandStart(b)) t(bandEnd(b)) + 1/fs], [3.5 3.5], -0.5, ...
        'FaceColor', [0.85 0.85 0.85], 'EdgeColor', 'none'); hold on;
end
stairs(t, spfLog.mode(1:n), 'k', 'LineWidth', 1.2); hold on;
plot(ev(spfLog.evLatch), 2*ones(size(ev(spfLog.evLatch))), 'rv', 'MarkerFaceColor', 'r');
plot(ev(spfLog.evProbStart), 3*ones(size(ev(spfLog.evProbStart))), 'b^', 'MarkerFaceColor', 'b');
plot(ev(spfLog.evVeto), 2*ones(size(ev(spfLog.evVeto))), 'mx', 'MarkerSize', 8, 'LineWidth', 1.5);
plot(ev(spfLog.evHandback), ones(size(ev(spfLog.evHandback))), 'go', 'MarkerFaceColor', 'g');
yticks([1 2 3]); yticklabels({'NOMINAL', 'COAST', 'PROBATION'}); ylim([0.5 3.5]); grid on;
ylabel('mode'); title('gate mode (grey: warm-up, not armed; v latch, ^ probation, x veto, o hand-back)');
subplot(2,1,2);
stairs(t, spfLog.numMeas(1:n), 'k'); grid on; ylabel('numMeas'); xlabel('t [s]');

%% 2. monitor margins
figure('Name', 'SPF gate: monitor margins');
subplot(5,1,1);
plot(t, spfLog.ssRatio(1:n, :)); hold on; yline(1, 'r--'); grid on;
ylabel('SS |d|/(k_{FA}\sigma_{SS})'); legend('N', 'E', 'D', 'alarm', 'Location', 'northwest');
title('solution-separation margin (alarm > 1)');
subplot(5,1,2);
plot(t, spfLog.cpiRatio(1:n, :)); hold on; yline(1, 'r--'); grid on;
ylabel('CPI q/T_N'); title('CPI margin of windows closed this epoch (alarm > 1)');
subplot(5,1,3);
semilogy(t, max(spfLog.pl(1:n), 1e-2), 'k', 'LineWidth', 1.2); hold on;
semilogy(t, max(spfLog.sigmaPos(1:n, :), 1e-2)); grid on;
ylabel('[m]'); legend('PL', '\sigma_N', '\sigma_E', '\sigma_D', 'Location', 'northwest');
title('max protection level and host position sigma sqrt(P+) (compare sigma with the true drift after a coast)');
subplot(5,1,4);
semilogy(t, max(spfLog.nis(1:n), 1e-3), 'k'); hold on; yline(1, 'r--'); grid on;
ylabel('NIS'); title('host KF normalised innovation squared (~1 when R, P honest)');
subplot(5,1,5);
plot(t, spfLog.innMean(1:n), 'k'); hold on; plot(t, spfLog.innStd(1:n), 'b'); grid on;
ylabel('[m]'); xlabel('t [s]'); legend('mean (common mode: clock)', 'spread (geometry: position)', 'Location', 'northwest');
title('innovation mean and spread over the live rows');

%% 3. re-validation
figure('Name', 'SPF gate: re-validation');
subplot(3,1,1);
q = spfLog.qReval(1:n); q(~spfLog.revalDone(1:n)) = NaN;
thr = spfLog.revalThr(1:n); thr(~spfLog.revalDone(1:n)) = NaN;
semilogy(t, max(q, 1e-3), 'b'); hold on; semilogy(t, thr, 'r--'); grid on;
ylabel('q_{reval}'); legend('q', 'threshold(numMeas)', 'Location', 'northwest');
title('GNSS vs coast chi-square (COAST: decision, PROBATION: diagnostic)');
subplot(3,1,2);
stairs(t, spfLog.dwell(1:n), 'b'); hold on; yline(double(CST_spfParam.REVAL_DWELL_REQUIRED), 'r--'); grid on;
ylabel('dwell'); title('consecutive re-validation passes');
subplot(3,1,3);
plot(t, spfLog.coastEpochs(1:n) / fs, 'k'); hold on;
yline(double(CST_spfParam.COAST_BUDGET_EPOCHS) / fs, 'r--'); grid on;
ylabel('coast time [s]'); xlabel('t [s]'); title('time in COAST + PROBATION (dashed: coast budget, host may reset beyond it)');

%% 4. faults
figure('Name', 'SPF gate: faults');
stairs(t, [spfLog.inputFault(1:n) spfLog.solveFault(1:n) + 2 spfLog.clamped(1:n) + 4]); grid on;
yticks([0 1 2 3 4 5]); yticklabels({'', 'inputFault', '', 'solveFault', '', 'numMeasClamped'});
ylim([-0.5 5.5]); xlabel('t [s]'); title('fault flags');
end
