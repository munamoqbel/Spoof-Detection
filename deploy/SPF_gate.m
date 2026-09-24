%******************************************************************************************
% DESCRIPTION:
% 2 Hz spoofing-detection gate (persistent state). Call once per GNSS epoch after the
% host update of the ACTIVE filter. Never touches the host filters: commands via spoofTel.
%
% INPUTS:  kfMeas (STRUCT_SPF.setKfMeas), propTel (STRUCT_SPF.setPropTel),
%          navActive (false = alignment: gate off, re-init later), resetRequest (force re-init)
% OUTPUTS: spoofTel (STRUCT_SPF.setTel)
%******************************************************************************************
%#codegen
function [spoofTel] = SPF_gate(kfMeas, propTel, navActive, resetRequest)

persistent sys epoch needInit

if isempty(needInit)
    needInit = true;
    sys      = STRUCT_SPF.zeroSys;
    epoch    = uint32(0);
end

if ~navActive
    needInit = true;
    epoch    = uint32(0);
    spoofTel = STRUCT_SPF.zeroTel;

elseif ~SPF_inputsFinite(kfMeas, propTel)
    % non-finite input: drop the epoch, re-init on the next good one
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
ok = all(isfinite(kfMeas.innovation(:))) && all(isfinite(kfMeas.innovationCov(:))) && ...
     all(isfinite(kfMeas.obsMatrix(:)))  && all(isfinite(kfMeas.measNoiseCov(:)))  && ...
     all(isfinite(kfMeas.priorState(:))) && all(isfinite(kfMeas.postState(:)))     && ...
     all(isfinite(kfMeas.postCov(:)))    && all(isfinite(propTel.accumPhi(:)))     && ...
     all(isfinite(propTel.accumQ(:)));
end
%------------------------------------------------------------------------
