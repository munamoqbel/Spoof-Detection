%***************************************************************************
% DESCRIPTION:
% The 2 Hz spoofing-detection gate of the host GNSS/INS EKF, with its
% persistent state. Called once per GNSS epoch right after the host's
% update of the ACTIVE filter (the operational KF in NOMINAL and COAST,
% the trial copy in PROBATION), on every epoch, with rows or without.
% It never touches the host filters: it watches the update increments,
% declares the mode and commands the host through spoofTel, a complete
% x+ and P+ on a latch or a hand-back and the start of a trial filter
% when probation opens.
%
% INPUTS:
%   - kfMeas        STRUCT_SPF.setKfMeas or kfMeasFromUpdate: the rows
%                   the filter used (y, S, H, R, numMeas) and x_bar, x+,
%                   P+ of the ACTIVE filter's update
%   - propTel       STRUCT_SPF.setPropTel: interval Phi and Q accumulated
%                   on the 100 Hz side since the previous epoch
%   - navActive     false during alignment: the gate returns zeroTel,
%                   drops its state and re-initialises on the first
%                   active call afterwards
%   - resetRequest  true forces a full re-initialisation
%
% OUTPUTS:
%   - spoofTel      STRUCT_SPF.setTel (.info, .kfCommand, .nav)
%
% ASSUMPTIONS AND LIMITATIONS:
% - Position states are CST_spfParam.MONITORED_AXES (NED, 1:3).
% - Works on increments x+ - x_bar: invariant to the host's split of its
%   estimate between feedback and residual states.
% - A non-finite input drops the epoch (info.inputFault) and the gate
%   re-initialises on the next good one; a re-initialised gate is not
%   armed until the warm-up of SPF_protectedNav is over.
% - Persistents survive between runs: resetRequest = true or 'clear
%   functions' starts fresh.
%
% REQUIREMENT TRACEABILITY:
% Paper Sec. 5 recovery (detection, re-validation, probation, hand-back)
% in SPF_protectedNav; the monitors (Eq. 17 to 35, 49 to 52) in
% SPF_monitorPool.
%***************************************************************************
%#codegen
function [spoofTel] = SPF_gate(kfMeas, propTel, navActive, resetRequest)

persistent sys epoch needInit

if isempty(needInit)
    needInit = true;
    sys      = STRUCT_SPF.zeroSys;
    epoch    = uint32(0);
end

if ~navActive
    % alignment: gate off, state dropped
    needInit = true;
    epoch    = uint32(0);
    spoofTel = STRUCT_SPF.zeroTel;

elseif ~SPF_inputsFinite(kfMeas, propTel)
    % non-finite input: the epoch is dropped, re-initialisation on the next good one
    needInit = true;
    epoch    = uint32(0);
    spoofTel = STRUCT_SPF.zeroTel;
    spoofTel.info.inputFault = true;

else
    if resetRequest || needInit
        sys          = STRUCT_SPF.zeroSys;      % unarmed: warm-up in SPF_protectedNav
        sys.coastCov = kfMeas.postCov;
        epoch        = uint32(0);
        needInit     = false;
    end

    epoch = epoch + uint32(1);

    [sys, spoofTel] = SPF_protectedNav(sys, kfMeas, propTel, epoch);
end

end

function [ok] = SPF_inputsFinite(kfMeas, propTel)
%#codegen
% every field the gate reads must be finite
ok = all(isfinite(kfMeas.innovation(:))) && all(isfinite(kfMeas.innovationCov(:))) && ...
     all(isfinite(kfMeas.obsMatrix(:)))  && all(isfinite(kfMeas.measNoiseCov(:)))  && ...
     all(isfinite(kfMeas.priorState(:))) && all(isfinite(kfMeas.postState(:)))     && ...
     all(isfinite(kfMeas.postCov(:)))    && all(isfinite(propTel.accumPhi(:)))     && ...
     all(isfinite(propTel.accumQ(:)));
end
%------------------------------------------------------------------------
