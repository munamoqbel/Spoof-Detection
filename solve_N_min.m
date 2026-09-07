function [N_min, CPI_threshold, k_FA, k_MD, P_MD_achieved] = ...
    solve_N_min(sigma2_gamma_u, sigma_t, P_FA_plus, P_FA_minus, P_MD_plus, P_MD_minus)
%SOLVE_N_MIN  Offline design: minimum window length from Eq. (38).
%             NOT code-generated (uses gaminv/gammainc/norminv).
%
% Equations
% ---------
%   Omega  = sigma2_gamma_u * sigma_t^2             (Eq. 26 / below Eq. 37)
%   T_N    = gaminv(1 - P_FA_plus, N/2, 2)          (Eq. 35, upper tail)
%   P_MD   = gammainc(T_N/(2*(1+Omega)), N/2)       (Eq. 38, regularised lower)
%
% Inputs
%   sigma2_gamma_u   Eq. 20 from the SAME S_k the monitor will use at
%                    runtime. Use the MINIMUM over a representative run
%                    (worst geometry -> conservative design).
%   sigma_t          worst-case spoofer tracking error std dev (m)
%   P_FA_plus        CPI false-alert allocation
%   P_FA_minus       SS false-alert allocation (TOTAL; divided by N_min here)
%   P_MD_plus        CPI missed-detection requirement
%   P_MD_minus       SS  missed-detection requirement
%
% Outputs
%   N_min, CPI_threshold(T_N), k_FA, k_MD, P_MD_achieved
%
% NOTE on the two-pass pattern: k_FA depends on P_FA_minus/N_min, and
% N_min is unknown before the first solve. This function solves N_min
% first (independent of k_FA), then computes k_FA with the solved N_min -
% so a single call is sufficient; no external iteration needed.

Omega = sigma2_gamma_u * sigma_t^2;

fprintf('--- solve_N_min ---\n');
fprintf('sigma2_gamma_u = %.4f\n', sigma2_gamma_u);
fprintf('sigma_t        = %.4f m\n', sigma_t);
fprintf('Omega          = %.4f\n', Omega);

% ---- search N (Eq. 38 depends only on N, P_FA_plus, Omega) ----
N_min = 1;
P_MD  = 1.0;
while P_MD > P_MD_plus && N_min < 10000
    N_min = N_min + 1;
    T_N   = gaminv(1 - P_FA_plus, N_min/2.0, 2.0);                       % Eq. 35
    P_MD  = gammainc(T_N / (2.0*(1.0 + Omega)), N_min/2.0, 'lower'); % Eq. 38
end
if N_min >= 10000
    error(['solve_N_min: P_MD_plus=%.1e unreachable for sigma_t=%.3f m. ', ...
           'Increase sigma_t assumption or relax P_MD_plus.'], P_MD_plus, sigma_t);
end
P_MD_achieved = P_MD;

% ---- constants with the solved N_min ----
CPI_threshold = gaminv(1 - P_FA_plus, N_min/2.0, 2.0);
P_FA_ss_win   = P_FA_minus / N_min;              % Bonferroni across N_min windows
k_FA          = norminv(1 - P_FA_ss_win/2.0);    % two-sided SS threshold
k_MD          = norminv(1 - P_MD_minus);         % one-sided PL term

fprintf('N_min         = %d epochs\n', N_min);
fprintf('CPI_threshold = %.10f\n', CPI_threshold);
fprintf('k_FA          = %.10f\n', k_FA);
fprintf('k_MD          = %.10f\n', k_MD);
fprintf('P_MD achieved = %.4e  (req: %.4e)\n', P_MD_achieved, P_MD_plus);
fprintf('-------------------\n\n');

end
