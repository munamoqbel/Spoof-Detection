function [state, covariance, innovation, innovation_cov, kalman_gain] = ...
    kalman_update_step(state, covariance, measurement, obs_matrix, ...
                       spoofInfo, meas_noise_cov)
%KALMAN_UPDATE_STEP  One complete Kalman filter epoch (time + measurement).
%
% INPUTS
%   state            [n x 1]  filter state at previous epoch (post-update)
%   covariance       [n x n]  filter covariance at previous epoch
%   measurement      [m x 1]  GNSS measurement vector z
%   obs_matrix       [m x n]  observation matrix H
%   Phi, Q           [n x n]  state transition and process noise
%   meas_noise_cov   [m x m]  measurement noise covariance V
%
% OUTPUTS
%   state            [n x 1]  updated filter state
%   covariance       [n x n]  updated filter covariance
%   innovation       [m x 1]  gamma = z - H*x_bar                (paper Eq. 3)
%   innovation_cov   [m x m]  S = H*P_bar*H' + V
%   kalman_gain      [n x m]  L = P_bar*H'/S
%
% NOTE  Linear measurement model (harness). In a real EKF, replace
%       obs_matrix*predicted_state with your nonlinear h(x).
Phi = spoofInfo.phiAcc;
Q = spoofInfo.qAcc;
% ---- time update ----
predicted_state      = Phi * state;
predicted_covariance = Phi * covariance * Phi' + Q;

% ---- innovation ----
innovation     = measurement - obs_matrix * predicted_state;
innovation_cov = obs_matrix * predicted_covariance * obs_matrix' + meas_noise_cov;
innovation_cov = (innovation_cov + innovation_cov') / 2;   % keep symmetric

% ---- measurement update ----
kalman_gain = predicted_covariance * obs_matrix' / innovation_cov;
state       = predicted_state + kalman_gain * innovation;

num_states = size(covariance, 1);
covariance = (eye(num_states) - kalman_gain * obs_matrix) * predicted_covariance;
covariance = (covariance + covariance') / 2;                % keep symmetric

end
