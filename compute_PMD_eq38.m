function P_MD = compute_PMD_eq38(N, sigma2_gamma_u, sigma_t, P_FA_plus)
%COMPUTE_PMD_EQ38  Exact Eq. (38). Offline analysis only (not code-generated).
%
%   sigma2_Delta = sigma2_gamma_u^2 * sigma_t^2          (Eq. 26)
%   Omega        = sigma2_Delta / sigma2_gamma_u
%                = sigma2_gamma_u * sigma_t^2
%   T_N          = gaminv(1 - P_FA_plus, N/2, 2)         (Eq. 35)
%   P_MD         = gamma(N/2, T_N/(2(1+Omega)))/Gamma(N/2)
%                = gammainc(T_N/(2*(1+Omega)), N/2, 'lower')   (Eq. 38)

Omega = sigma2_gamma_u * sigma_t^2;
T_N   = gaminv(1 - P_FA_plus, N/2.0, 2.0);
P_MD  = gammainc(T_N / (2.0*(1.0 + Omega)), N/2.0, 'lower');

end
