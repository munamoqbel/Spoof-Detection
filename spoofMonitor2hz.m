%******************************************************************************************
% DESCRIPTION:
% 2 Hz spoofing-detection GATE for the host GNSS/INS EKF (persistent state).
%
% Call ONCE per GNSS epoch from the host's 2 Hz function, right after the
% host's kfUpdate on the ACTIVE filter (operational KF in NOMINAL/COAST,
% trial copy in PROBATION). Call it on EVERY 2 Hz epoch, also when no GNSS
% rows are available (numMeas = 0): windows age, coast covariances
% propagate and increments accumulate only when the gate runs. The gate
% never touches the host filters: it watches the update increments (CPI +
% SS monitors on overlapping windows), declares the mode, tests the raw
% GNSS against the certified-clean coast, runs a probation, and hands back.
% It only ever COMMANDS the host through spoofTel:
%
%   spoofTel.info.mode              NOMINAL / COAST / PROBATION (CST_spfMode)
%   spoofTel.nav.applyCorrection    true on LATCH and COMMIT epochs:
%   spoofTel.nav.state / .covar     a complete (x+, P+) for the operational
%                                   KF -> through the host's normal setKF
%   spoofTel.kfCommand.startTrial   true when probation opens: trial <- copy
%                                   of the operational KF (states + covariance)
%   spoofTel.info.*                 alarms per axis, PL, qReval, dwell,
%                                   coastEpochs, events
%
% INPUTS:
%   - kfMeas        STRUCT_SPF.setKfMeas / kfMeasFromUpdate: GNSS rows ONLY
%                   (y, H, R of the range / range-rate rows; the pressure
%                   row excluded), numMeas = GNSS rows (0 = no GNSS), and
%                   xPrior, xPost, PPost of the ACTIVE filter's update
%   - propTel       STRUCT_SPF.setPropTel(accumPhi, accumQ): interval Phi / Q
%                   accumulated on the 100 Hz side since the previous epoch
%   - navActive     logical, true in NAVIGATION mode. False (alignment):
%                   the gate returns STRUCT_SPF.zeroTel (NOMINAL, no
%                   commands), discards its state and epoch counter, and
%                   re-initialises on the first active call afterwards
%   - resetRequest  logical, true forces a full re-init (host first call)
%
% OUTPUTS:
%   - spoofTel      STRUCT_SPF.setTel
%
% ASSUMPTIONS AND LIMITATIONS:
% - Position states are CST_spfParam.MONITORED_AXES (NED, 1:3).
% - Works in increment space: invariant to the host's feedback/gain
%   bookkeeping of KF.states (see STRUCT_SPF header).
% - Persistents survive between runs: pass resetRequest = true (or
%   'clear functions') to start fresh.
%
% REQUIREMENT TRACEABILITY:
%
%******************************************************************************************
%#codegen
function [spoofTel] = spoofMonitor2hz(kfMeas, propTel, navActive, resetRequest)

persistent sys epoch needInit

if isempty(needInit)
    needInit = true;
    sys      = STRUCT_SPF.zeroSys;
    epoch    = uint32(0);
end

if ~navActive
    % alignment (or any non-navigation mode): gate off, re-arm
    needInit = true;
    epoch    = uint32(0);
    spoofTel = STRUCT_SPF.zeroTel;
    return;
end

if resetRequest || needInit
    sys = STRUCT_SPF.zeroSys;
    sys.coastCov = kfMeas.postCov;
    % the initial solution is the fallback until the first window closes clean
    sys.anchor = STRUCT_SPF.setAnchor(true, zeros(CST_gnssHybrid.NO_STATES, 1), ...
        kfMeas.postCov, uint32(0));
    epoch    = uint32(0);
    needInit = false;
end

epoch = epoch + uint32(1);

[sys, spoofTel] = protectedNav(sys, kfMeas, propTel, epoch);

end
% ------------------------------------------------------------------------
