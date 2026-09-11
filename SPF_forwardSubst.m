%******************************************************************************************
% DESCRIPTION:
% Fixed-size forward substitution: solve L * w = b for a lower-triangular L
% with a strictly positive diagonal (the factor from SPF_cholesky). Plain
% loops over the leading numActive rows; the padding rows of w are zero
% (b is zero there by construction). Compiles with MATLAB Coder with
% variable sizing disabled, no library solve.
%
% INPUTS:
%   - L          [n x n]  lower-triangular, diag(L) > 0
%   - b          [n x 1]
%   - numActive  scalar   rows of the live block (<= n)
%
% OUTPUTS:
%   - w      [n x 1]  L \ b on the live block, zeros in the padding
%
% ASSUMPTIONS AND LIMITATIONS:
% Only call with the factor of a successful SPF_cholesky (ok = true), so
% no diagonal entry of the live block is zero.
%
% REQUIREMENT TRACEABILITY:
%
%******************************************************************************************
%#codegen
function [w] = SPF_forwardSubst(L, b, numActive)

% Define variables
n       = size(L, 1);
nActive = min(double(numActive), n);
w       = zeros(n, 1);

for rowIdx = 1:nActive
    t = b(rowIdx);
    for k = 1:(rowIdx - 1)
        t = t - L(rowIdx, k) * w(k);
    end
    w(rowIdx) = t / L(rowIdx, rowIdx);
end

end
%------------------------------------------------------------------------
