%******************************************************************************************
% DESCRIPTION:
% 2 Hz spoofing-detection GATE for the host GNSS/INS EKF (persistent state).
%
% Call ONCE per GNSS epoch from the host's 2 Hz function, right after the
% host's kfUpdate on the ACTIVE filter (operational KF in NOMINAL/COAST,
% trial copy in PROBATION). The gate never touches the host filters: it
% watches the update increments (CPI + SS monitors on overlapping windows),
% declares the mode, tests the raw GNSS against the certified-clean coast,
% runs a probation, and hands back. It only ever COMMANDS the host through
% spoofTel:
%
%   spoofTel.info.mode              NOMINAL / COAST / PROBATION (CST_spfMode)
%   spoofTel.nav.applyCorrection    true on LATCH and COMMIT epochs:
%   spoofTel.nav.correction         x+ = KF.states + correction, P+ = nav.covar
%   spoofTel.nav.covar              -> run the usual setKF bookkeeping on it
%   spoofTel.kfCommand.reseedKF     true when probation opens: trial <- copy
%                                   of the operational KF (cov reseedCov)
%   spoofTel.info.*                 alarms per axis, PL, qReval, dwell,
%                                   coastEpochs, events
%
% INPUTS:
%   - kfMeas        STRUCT_SPF.setKfMeas(y, S, H, R, numMeas, xPrior, xPost, PPost)
%                   from this epoch's kfUpdate; rows 1:numMeas valid, arrays
%                   may be exact-size or padded to CST_spfParam.MAX_MEAS
%   - propTel       STRUCT_SPF.setPropTel(accumPhi, accumQ): interval Phi / Q
%                   accumulated on the 100 Hz side since the previous epoch
%   - resetRequest  logical, true forces a full re-init (first call /
%                   commanded restart)
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
function [spoofTel] = spoofMonitor2hz(kfMeas, propTel, resetRequest)

persistent sys epoch

if resetRequest || isempty(sys)
    sys = STRUCT_SPF.zeroSys;
    sys.coastCov = kfMeas.postCov;
    % the initial solution is the fallback until the first window closes clean
    sys.anchor = STRUCT_SPF.setAnchor(true, zeros(CST_gnssHybrid.NO_STATES, 1), ...
        kfMeas.postCov, uint32(0));
    epoch = uint32(0);
end

epoch = epoch + uint32(1);

[sys, spoofTel] = protectedNav(sys, kfMeas, propTel, epoch);

end
% ------------------------------------------------------------------------
