%******************************************************************************************
% DESCRIPTION:
% The gate's own process noise for a coast: the IMU datasheet on the inertial
% states, the host's accumQ everywhere else.
%
% An INS-only coast drifts by the sensor: the velocity by the velocity random
% walk, the attitude by the angle random walk, the biases by their in-run
% stability and walk. The host's Q is a filter tuning, inflated to keep the
% filter responsive under the GNSS, and carried into a coast it overstates
% the drift several times (measured on the vehicle host: three in sigma).
% Every gate threshold is k times the coast sigma, so the gate propagates its
% coasts with the sensor's Q, built here from CST_spfParam.IMU_* along the
% host's own accumPhi: dt/2 (Phi Qc Phi' + Qc) on the inertial states, the
% position block without noise (an INS position error comes through the
% velocity and the attitude), and the host's accumQ on the states the IMU
% does not drive (clock, pressure, GNSS states). dt is the interval the
% accumulation covers, accumPhi(pos, vel), or CST_spfParam.COAST_Q_INTERVAL
% when the host's Phi carries no position-from-velocity entry.
%
% INPUTS:
%   accumPhi  [n x n]  the host's interval transition (STRUCT_SPF.setPropTel)
%   accumQ    [n x n]  the host's interval process noise
%
% OUTPUTS:
%   coastQ    [n x n]  the gate's interval process noise, symmetric
%
% ASSUMPTIONS AND LIMITATIONS:
% - CST_spfParam.COAST_*_STATES name the position, velocity, attitude, gyro
%   bias and accelerometer bias states in the host's layout.
% - Used by SPF_coastPropTel when CST_spfParam.COAST_Q_DATASHEET is true.
%******************************************************************************************
%#codegen
function [coastQ] = SPF_coastQ(accumPhi, accumQ)

numStates = CST_gnssHybrid.NO_STATES;
R2D = 180.0 / pi;
G0  = 9.80665;
posIdx  = CST_spfParam.COAST_POS_STATES;
velIdx  = CST_spfParam.COAST_VEL_STATES;
attIdx  = CST_spfParam.COAST_ATT_STATES;
gyroIdx = CST_spfParam.COAST_GYRO_BIAS_STATES;
accIdx  = CST_spfParam.COAST_ACC_BIAS_STATES;

% the continuous densities, SI per second
qVel  = (CST_spfParam.IMU_VRW / 60.0)^2;                                  % (m/s)^2 per s
qAtt  = (CST_spfParam.IMU_ARW / R2D / 60.0)^2;                            % rad^2 per s
sigG  = CST_spfParam.IMU_GYRO_BIAS_SIGMA / R2D / 3600.0;                  % rad/s
sigA  = CST_spfParam.IMU_ACC_BIAS_SIGMA * G0 * 1.0e-6;                    % m/s^2
walkG = CST_spfParam.IMU_GYRO_WALK / R2D / 3600.0 / 60.0;                 % rad/s per sqrt(s)
walkA = CST_spfParam.IMU_ACC_WALK * G0 * 1.0e-6 / 60.0;                   % m/s^2 per sqrt(s)
qGyro = 2.0 * sigG^2 / CST_spfParam.IMU_BIAS_TAU + walkG^2;               % the Markov bias as a walk over a coast, plus the walk
qAcc  = 2.0 * sigA^2 / CST_spfParam.IMU_BIAS_TAU + walkA^2;

contQ = zeros(numStates, numStates);
for k = 1:numel(velIdx)
    contQ(velIdx(k), velIdx(k)) = qVel;
end
for k = 1:numel(attIdx)
    contQ(attIdx(k), attIdx(k)) = qAtt;
end
for k = 1:numel(gyroIdx)
    contQ(gyroIdx(k), gyroIdx(k)) = qGyro;
end
for k = 1:numel(accIdx)
    contQ(accIdx(k), accIdx(k)) = qAcc;
end

% the interval the accumulation covers
dt = CST_spfParam.COAST_Q_INTERVAL;
if (numel(posIdx) >= 1) && (numel(velIdx) >= 1)
    dtPhi = accumPhi(posIdx(1), velIdx(1));
    if (dtPhi > 0.0) && isfinite(dtPhi)
        dt = dtPhi;
    end % ELSE: no position-from-velocity entry in this Phi: the constant interval
end % ELSE: no position or velocity states named: the constant interval

% the inertial states: the host's entries go, the datasheet's come
inertial = false(numStates, 1);
for k = 1:numel(posIdx)
    inertial(posIdx(k)) = true;
end
for k = 1:numel(velIdx)
    inertial(velIdx(k)) = true;
end
for k = 1:numel(attIdx)
    inertial(attIdx(k)) = true;
end
for k = 1:numel(gyroIdx)
    inertial(gyroIdx(k)) = true;
end
for k = 1:numel(accIdx)
    inertial(accIdx(k)) = true;
end
datasheetQ = (dt / 2.0) * (accumPhi * contQ * accumPhi' + contQ);
coastQ = accumQ;
for i = 1:numStates
    for j = 1:numStates
        if inertial(i) && inertial(j)
            coastQ(i, j) = datasheetQ(i, j);
        elseif inertial(i) || inertial(j)
            coastQ(i, j) = 0.0;
        end % ELSE: neither state is inertial: the host's entry stays
    end
end
coastQ = (coastQ + coastQ') / 2.0;
end
