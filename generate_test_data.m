function [gamma_all, S_all, H_all, x_hat_all, P_hat_all, Phi, Q, scn, z_all, V, P0] = ...
    generate_test_data(prm_boot, seed, sigma_t_attack, ramp_rate, attack_dir)
%GENERATE_TEST_DATA  Realistic EKF test data: clean / spoofed / clean.
%
% WHY A GENUINE CLOSED-LOOP KF IS REQUIRED
% ---------------------------------------
% Earlier attempts drew innovations as independent random vectors and
% did x_hat = Phi*x_hat + L*gamma. That is NOT a Kalman filter:
% real innovations gamma = z - H*x_bar are anti-correlated with the
% state error and keep x_hat bounded near truth. With independent
% gamma, x_hat is an unbounded random walk, so u'(x_hat - x_C) grows
% without any spoofing and the SS monitor (whose Eq. 50 variance
% P_C - P_hat assumes true KF error dynamics, Appendix E) alarms
% constantly. This generator therefore simulates:
%   - a true (static) state,
%   - real measurements  z_k = H_k*x_true + v_k,  v ~ N(0,V),
%   - spoofing injected into the MEASUREMENT z (Eq. C16),
%   - a genuine Riccati-propagated KF producing gamma, S, x_hat, P_hat.
% Everything the monitors assume is then true by construction.
%
% SCENARIO (all indices in epochs; fs = prm_boot.fs)
%   Phase 1  clean          1 .. 80
%   Phase 2  spoofed       81 .. 180
%              WGN tracking error, sigma_t_attack, from epoch 81
%              position ramp, ramp_rate m/epoch, from epoch 120
%   Phase 3  clean again  181 .. 480
%
% ATTACK PARAMETERISATION (dual-role fix)
%   sigma_t_attack : what is actually INJECTED into z. Default 0.10 m.
%                    Set 0 for the ramp-only case.
%   ramp_rate      : injected ramp, m/epoch. Default 0.5.
%                    Set 0 for the tracking-only case.
%   scn.sigma_t    : the DESIGN value fed to solve_N_min. Stays 0.10
%                    regardless of the attack toggles. NEVER set this
%                    to zero (Omega = 0 makes Eq. 38 unsolvable).
%   The scalar randn() for the tracking error is drawn on every spoofed
%   epoch even when sigma_t_attack = 0, so the RNG stream is identical
%   across attack cases with the same seed (epoch-by-epoch comparable).
%
% RECOVERY PHYSICS (why Phase 3 is long)
%   At spoof end the KF holds a large spoofed offset. With mm-level
%   steady-state P the filter trusts its prior, so the offset decays
%   only by ~0.93x per epoch. Both monitors CORRECTLY keep alarming
%   until the residual offset falls to their mm-scale noise floor.
%   scn.recovery_buf = 150 epochs bounds this conservatively; alarms
%   inside (spoof_end, spoof_end+recovery_buf] are detection of the
%   still-corrupted KF, NOT false alarms.
%
% All arrays are FIXED SIZE (Coder-style); this file itself is host-side
% (uses randn) and is not code-generated.
%
% Inputs
%   prm_boot         kujur_params() bootstrap struct
%   seed             RNG seed (default 42)
%   sigma_t_attack   injected tracking-error std, m (default 0.10)
%   ramp_rate        injected ramp, m/epoch (default 0.5)
%
% Outputs
%   gamma_all   [m_meas x N]                 real KF innovations
%   S_all       [m_meas x m_meas x N]        H*P_bar*H' + V (per epoch)
%   H_all       [m_meas x n_states x N]
%   x_hat_all   [n_states x N]               post-update KF states
%   P_hat_all   [n_states x n_states x N]    post-update covariances
%   Phi, Q                                    constant matrices used
%   scn         scenario struct (see below)
%   z_all       [m_meas x N]  raw measurements AS THE RECEIVER SAW THEM
%                             (spoofing included) - input to recovery_nav_sim
%   V           [m_meas x m_meas]  measurement noise covariance
%   P0          [n_states x n_states] Riccati-converged initial covariance

if nargin < 2, seed = 42; end
if nargin < 3, sigma_t_attack = 0.10; end
if nargin < 4, ramp_rate = 0.5; end
if nargin < 5, attack_dir = [0 0 1]; end    % attack direction in position axes
if isscalar(attack_dir)                     % legacy: axis index 1/2/3
    tmp = zeros(3, 1);
    tmp(attack_dir) = 1;
    attack_dir = tmp;
end
attack_dir = attack_dir(:) / norm(attack_dir);   % unit vector: drag magnitude
rng(seed);                                        % along d is ramp_rate m/epoch

n_states = prm_boot.n_states;
m_meas   = prm_boot.m_meas;
n_sv     = prm_boot.n_sv;
idx_z    = prm_boot.idx_z;

% ---- scenario definition ----
scn.N_total        = 480;
scn.clean_end      = 80;
scn.spoof_start    = 81;
scn.ramp_start     = 120;
scn.spoof_end      = 180;
scn.clean2_start   = 181;
scn.recovery_buf   = 150;     % epochs after spoof_end before quiet is expected
scn.sigma_t        = 0.10;    % DESIGN tracking-error std (m) -> solve_N_min
scn.sigma_t_attack = sigma_t_attack;    % INJECTED tracking-error std (m)
scn.ramp_rate      = ramp_rate;         % INJECTED ramp (m/epoch), along attack_dir
scn.attack_dir     = attack_dir;        % unit attack direction [3x1]
[~, scn.idx_attack] = max(abs(attack_dir));    % dominant axis, for display only
scn.sigma_code     = 0.36;    % Appendix A
scn.sigma_carrier  = 0.003;

N = scn.N_total;

[Phi, Q] = build_Phi_Q(prm_boot);

% ---- measurement noise (interleaved code/carrier per SV) ----
V_diag = zeros(m_meas, 1);
for sv = 1:n_sv
    V_diag(2*sv-1) = scn.sigma_code^2;
    V_diag(2*sv)   = scn.sigma_carrier^2;
end
V     = diag(V_diag);
sqrtV = sqrt(V_diag);

% ---- LOS geometry (per epoch) ----
H_all = zeros(m_meas, n_states, N);
for k = 1:N
    los = randn(n_sv, 3);
    los = los ./ vecnorm(los, 2, 2);
    for sv = 1:n_sv
        H_all(2*sv-1, 1:3, k) = -los(sv,:);
        H_all(2*sv,   1:3, k) = -los(sv,:);
    end
end

% ---- Riccati warm-up: converge P before epoch 1 ----
P_hat = 1.0 * eye(n_states);
H_1   = H_all(:,:,1);
for it = 1:200
    P_bar = Phi * P_hat * Phi' + Q;
    S_1   = H_1 * P_bar * H_1' + V;
    L_1   = P_bar * H_1' / S_1;
    P_hat = (eye(n_states) - L_1 * H_1) * P_bar;
    P_hat = (P_hat + P_hat') / 2;
end
P0 = P_hat;                            % exported: converged initial covariance

% ---- allocate outputs ----
gamma_all = zeros(m_meas, N);
S_all     = zeros(m_meas, m_meas,    N);
x_hat_all = zeros(n_states, N);
P_hat_all = zeros(n_states, n_states, N);
z_all     = zeros(m_meas, N);

% ---- closed-loop simulation ----
x_true = zeros(n_states, 1);     % static truth (relative frame)
x_hat  = zeros(n_states, 1);     % KF starts at truth (converged filter)

for k = 1:N
    H_k = H_all(:,:,k);

    % truth propagation (static under this Phi: zero velocity)
    x_true = Phi * x_true;

    % real measurement
    z_k = H_k * x_true + sqrtV .* randn(m_meas, 1);

    % ---- spoofing injected into z (Eq. C16): z^s = z + H*u*offset ----
    if k >= scn.spoof_start && k <= scn.spoof_end
        w    = randn();                 % drawn EVERY spoofed epoch (RNG alignment)
        Hd   = H_k(:, 1:3) * scn.attack_dir;       % H*d: 3-D attack signature
        z_k = z_k + Hd * (scn.sigma_t_attack * w);
        if k >= scn.ramp_start
            ramp_m = scn.ramp_rate * (k - scn.ramp_start + 1);
            z_k    = z_k + Hd * ramp_m;
        end
    end
    % Phase 3: nothing injected - spoofer gone.

    z_all(:,k) = z_k;                 % log what the receiver reported

    % ---- genuine KF time + measurement update ----
    x_bar = Phi * x_hat;
    P_bar = Phi * P_hat * Phi' + Q;

    S_k     = H_k * P_bar * H_k' + V;
    S_k     = (S_k + S_k') / 2;
    gamma_k = z_k - H_k * x_bar;                          % REAL innovation

    L_k   = P_bar * H_k' / S_k;
    x_hat = x_bar + L_k * gamma_k;
    P_hat = (eye(n_states) - L_k * H_k) * P_bar;
    P_hat = (P_hat + P_hat') / 2;

    % ---- log ----
    gamma_all(:,k)   = gamma_k;
    S_all(:,:,k)     = S_k;
    x_hat_all(:,k)   = x_hat;
    P_hat_all(:,:,k) = P_hat;
end

end
