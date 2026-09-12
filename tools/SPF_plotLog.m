function SPF_plotLog(spfLog, fs)
%SPF_PLOTLOG  Standard diagnostic figures for a gate log (SPF_logInit /
%   SPF_logAppend). fs = 2 Hz epoch rate (default 2).
%
%   Figure 1  mode timeline with events; numMeas
%   Figure 2  monitor margins: SS ratio and CPI ratio per axis (alarm > 1),
%             protection level, host NIS
%   Figure 3  re-validation: qReval vs its threshold, dwell count, coast time
%   Figure 4  fault flags
%
%   How to read it (jamming / false-latch investigation):
%   - NIS >> 1 while jammed and the SS/CPI ratios cross 1 at the same time:
%     the KF's R is optimistic under jamming; the monitors are reporting a
%     real model inconsistency. Fix R (or exclude the rows) on the host side
%     before touching the gate thresholds.
%   - qReval stays above its threshold long after the jamming ends, and
%     grows with coast time: coast capture (P_C too small for the real INS
%     drift). Needs the coast-time budget / Q, not the thresholds.
%   - probation opens and is vetoed repeatedly: look at which ratio vetoes
%     (SS = the trial drifts from the coast, CPI = biased innovations).

if nargin < 2, fs = 2; end
n = spfLog.count;
if n > spfLog.numEpochs, n = spfLog.numEpochs; end
t = (0:n-1)' / fs;
ev = @(idx) t(idx(idx <= n));

%% 1. mode and measurement count
figure('Name', 'SPF gate: mode');
subplot(2,1,1);
stairs(t, spfLog.mode(1:n), 'k', 'LineWidth', 1.2); hold on;
plot(ev(spfLog.evLatch), 2*ones(size(ev(spfLog.evLatch))), 'rv', 'MarkerFaceColor', 'r');
plot(ev(spfLog.evProbStart), 3*ones(size(ev(spfLog.evProbStart))), 'b^', 'MarkerFaceColor', 'b');
plot(ev(spfLog.evVeto), 2*ones(size(ev(spfLog.evVeto))), 'mx', 'MarkerSize', 8, 'LineWidth', 1.5);
plot(ev(spfLog.evHandback), ones(size(ev(spfLog.evHandback))), 'go', 'MarkerFaceColor', 'g');
yticks([1 2 3]); yticklabels({'NOMINAL', 'COAST', 'PROBATION'}); ylim([0.5 3.5]); grid on;
ylabel('mode'); title('gate mode (v latch, ^ probation, x veto, o hand-back)');
subplot(2,1,2);
stairs(t, spfLog.numMeas(1:n), 'k'); grid on; ylabel('numMeas'); xlabel('t [s]');

%% 2. monitor margins
figure('Name', 'SPF gate: monitor margins');
subplot(4,1,1);
plot(t, spfLog.ssRatio(1:n, :)); hold on; yline(1, 'r--'); grid on;
ylabel('SS |d|/(k_{FA}\sigma_{SS})'); legend('N', 'E', 'D', 'alarm', 'Location', 'northwest');
title('solution-separation margin (alarm > 1)');
subplot(4,1,2);
plot(t, spfLog.cpiRatio(1:n, :)); hold on; yline(1, 'r--'); grid on;
ylabel('CPI q/T_N'); title('CPI margin of windows closed this epoch (alarm > 1)');
subplot(4,1,3);
plot(t, spfLog.pl(1:n), 'k'); grid on; ylabel('PL [m]'); title('max protection level');
subplot(4,1,4);
semilogy(t, max(spfLog.nis(1:n), 1e-3), 'k'); hold on; yline(1, 'r--'); grid on;
ylabel('NIS'); xlabel('t [s]'); title('host KF normalised innovation squared (~1 when R, P honest)');

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
plot(t, spfLog.coastEpochs(1:n) / fs, 'k'); grid on; ylabel('coast time [s]'); xlabel('t [s]');

%% 4. faults
figure('Name', 'SPF gate: faults');
stairs(t, [spfLog.inputFault(1:n) spfLog.solveFault(1:n) + 2 spfLog.clamped(1:n) + 4]); grid on;
yticks([0 1 2 3 4 5]); yticklabels({'', 'inputFault', '', 'solveFault', '', 'numMeasClamped'});
ylim([-0.5 5.5]); xlabel('t [s]'); title('fault flags');
end
