function [Phi, Q] = build_Phi_Q(prm)
%BUILD_PHI_Q  State transition and process noise matrices.
%             MATLAB Coder compatible (fixed sizes, simple loops).
%
% In your REAL EKF, pass your own Phi_k/Q_k to dual_monitor instead -
% this function only supplies matrices for the standalone simulation.
%
% State layout (Appendix A): 1:3 pos | 4:6 vel | 7:9 att |
% 10:12 acc bias | 13:15 gyro bias | 16:17 clock | 18:end GNSS states.
%#codegen

n  = prm.n_states;
dt = prm.dt;

Phi = eye(n);
Q   = zeros(n);

% pos <- pos + vel*dt
for i = 1:3
    Phi(i, i+3) = dt;
end

% velocity driven by accelerometer VRW (Appendix F)
for i = 4:6
    Q(i,i) = prm.VRW_sq * dt;
end
% attitude driven by gyro ARW
for i = 7:9
    Q(i,i) = prm.ARW_sq * dt;
end
% IMU bias instabilities (Eq. A3)
for i = 10:12
    Q(i,i) = prm.acc_bias_sq * dt;
end
for i = 13:15
    Q(i,i) = prm.gyro_bias_sq * dt;
end
% receiver clock (Eq. A7)
Q(16,16) = 1.0e-4 * dt;
Q(17,17) = 1.0e-8 * dt;
end
