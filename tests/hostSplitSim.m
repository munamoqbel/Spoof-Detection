function out = hostSplitSim(z_all, H_all, V, propTel, P0, prm, gain, drainCoast, drainProbation)
%HOSTSPLITSIM  Test harness: a host with the REAL feedback split of setKF
%              (stateFB = GAIN*x+, states = x+ - GAIN*x+), following the
%              gate contract in docs/HOST_2HZ_WIRING.m.
%
% The mechanization is emulated by 'fed', the corrections it has absorbed
% (propagated with Phi); 'r' is KF.states (the residual error still in the
% mechanization). The navigation OUTPUT is fed; the filter ESTIMATE is
% fed + r, which must equal the bare-state host of recovery_nav_sim.
%
% INPUTS  as recovery_nav_sim, plus
%   gain            setKF feedback gain (one constant)
%   drainCoast      true: on COAST epochs apply the split to the extrapolated
%                   state (stateFB = GAIN*r, r = r - GAIN*r)
%   drainProbation  true: same on PROBATION epochs (WRONG: the trial shares
%                   the mechanization; included to demonstrate)
%
% OUTPUT out: estimate [n x N] (fed + r), output [n x N] (fed), residual
%   [n x N] (r), mode [1 x N], ev.* as recovery_nav_sim, latchState [n x 1]
%   (nav.state at the first latch).

n = prm.n_states; N = size(z_all, 2); m = size(z_all, 1);
Phi = propTel.accumPhi; Q = propTel.accumQ; I = eye(n);

sys = STRUCT_SPF.zeroSys;
sys.coastCov = P0;
sys.anchor   = STRUCT_SPF.setAnchor(true, zeros(n, 1), P0, uint32(0));

fed = zeros(n, 1);                 % mechanization corrections (output)
r   = zeros(n, 1);  P  = P0;       % operational KF.states / covariance
t   = zeros(n, 1);  tP = P0;       % trial
mode = CST_spfMode.NOMINAL;

out.estimate = zeros(n, N); out.output = zeros(n, N); out.residual = zeros(n, N);
out.mode = zeros(1, N); out.latchState = zeros(n, 1);
ev_detect = []; ev_prob_start = []; ev_prob_fail = []; ev_handback = [];

for k = 1:N
    inProbation = (mode == CST_spfMode.PROBATION);

    % ---- 100 Hz side: everything extrapolates ----
    fed  = Phi * fed;
    rBar = Phi * r;  PBar = Phi * P * Phi' + Q;
    if inProbation
        tBar = Phi * t;  tPBar = Phi * tP * Phi' + Q;
        xBar = tBar;     PB = tPBar;
    else
        xBar = rBar;     PB = PBar;
    end

    % ---- kfUpdate on the active filter (host innovation: z - h(mech - residual)) ----
    H = H_all(:, :, k);
    y = z_all(:, k) - H * (fed + xBar);
    S = H * PB * H' + V;  S = (S + S') / 2;
    K = PB * H' / S;
    xPlus = xBar + K * y;
    PPlus = (I - K * H) * PB;  PPlus = (PPlus + PPlus') / 2;

    kfMeas = STRUCT_SPF.kfMeasFromUpdate(y, H, V, m, xBar, PB, xPlus, PPlus);
    [sys, tel] = SPF_protectedNav(sys, kfMeas, propTel, k);
    mode = tel.info.mode;

    % ---- host: ONE setKF per epoch on the operational KF ----
    if tel.nav.applyCorrection                       % LATCH or COMMIT: normal case
        xs  = tel.nav.state;
        fed = fed + gain * xs;  r = xs - gain * xs;  P = tel.nav.covar;
        if tel.info.eventLatched && ~any(out.latchState), out.latchState = xs; end
    elseif ~inProbation && (mode == CST_spfMode.NOMINAL)   % normal update
        fed = fed + gain * xPlus;  r = xPlus - gain * xPlus;  P = PPlus;
    else                                             % operational KF coasts
        r = rBar;  P = PBar;
        if (mode == CST_spfMode.COAST && drainCoast) || ...
           (mode == CST_spfMode.PROBATION && drainProbation)
            fed = fed + gain * r;  r = r - gain * r;  % drain (covariance untouched)
        end
    end
    if inProbation, t = xPlus; tP = PPlus; end       % trial keeps its update, no split
    if tel.kfCommand.startTrial, t = r; tP = P; end

    out.estimate(:, k) = fed + r;  out.output(:, k) = fed;  out.residual(:, k) = r;
    out.mode(k) = double(mode);
    if tel.info.eventLatched,          ev_detect(end+1)     = k; end %#ok<AGROW>
    if tel.info.eventProbationStarted, ev_prob_start(end+1) = k; end %#ok<AGROW>
    if tel.info.eventProbationVetoed,  ev_prob_fail(end+1)  = k; end %#ok<AGROW>
    if tel.info.eventHandback,         ev_handback(end+1)   = k; end %#ok<AGROW>
end
out.ev.t_detect = ev_detect; out.ev.t_prob_start = ev_prob_start;
out.ev.t_prob_fail = ev_prob_fail; out.ev.t_handback = ev_handback;
end
