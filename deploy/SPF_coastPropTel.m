%******************************************************************************************
% DESCRIPTION:
% The interval matrices the gate propagates its coasts with: the host's
% accumPhi, and either the host's accumQ as it comes or the gate's own Q
% from the IMU datasheet (SPF_coastQ), by CST_spfParam.COAST_Q_DATASHEET.
%
% INPUTS:
%   propTel  STRUCT_SPF.setPropTel: the host's accumPhi, accumQ
%
% OUTPUTS:
%   propTelOut  the same struct with accumQ replaced when the flag is set
%******************************************************************************************
%#codegen
function [propTelOut] = SPF_coastPropTel(propTel)

propTelOut = propTel;
if CST_spfParam.COAST_Q_DATASHEET
    propTelOut.accumQ = SPF_coastQ(propTel.accumPhi, propTel.accumQ);
end
end
