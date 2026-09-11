%******************************************************************************************
% DESCRIPTION:
% Fixed-size Cholesky factorisation A = L * L' (L lower triangular) with a
% failure flag instead of an error. Written as plain loops so that it
% compiles with MATLAB Coder with variable sizing disabled (the library
% [R, p] = chol(A) returns a variable-size R on failure) and needs no
% LAPACK on the target.
%
% Only the leading numActive x numActive block is factorised (loop bounds
% may be run-time values; array sizes stay fixed). Rows/cols beyond it are
% the padding of the MAX_MEAS layout: their off-diagonal entries must be
% zero and their diagonal positive, and L there is simply the square root
% of the diagonal (A = blkdiag(A_active, D) -> L = blkdiag(L_active, sqrt(D))).
%
% INPUTS:
%   - A          [n x n]  symmetric matrix (only the lower triangle is read)
%   - numActive  scalar   rows/cols of the live block (<= n)
%
% OUTPUTS:
%   - L      [n x n]  lower-triangular factor; zeros where the factorisation
%                     failed
%   - ok     logical  true if every pivot was finite and strictly positive
%                     (A positive definite up to rounding); false for a
%                     singular, indefinite or non-finite A
%
% ASSUMPTIONS AND LIMITATIONS:
% A zero or negative pivot stops the factorisation (ok = false); the
% caller decides what to do (drop the epoch, raise a fault flag). The
% relative pivot floor (CST_spfParam.PIVOT_REL_TOL) is applied by the
% caller on diag(L).
%
% REQUIREMENT TRACEABILITY:
%
%******************************************************************************************
%#codegen
function [L, ok] = SPF_cholesky(A, numActive)

% Define variables
n       = size(A, 1);
nActive = min(double(numActive), n);
L       = zeros(n, n);
ok      = true;

for colIdx = 1:nActive
    if (ok)
        % pivot
        pivot = A(colIdx, colIdx);
        for k = 1:(colIdx - 1)
            pivot = pivot - L(colIdx, k) * L(colIdx, k);
        end
        if (pivot > 0.0) && isfinite(pivot)
            L(colIdx, colIdx) = sqrt(pivot);
            % column below the pivot (live block only)
            for rowIdx = (colIdx + 1):nActive
                t = A(rowIdx, colIdx);
                for k = 1:(colIdx - 1)
                    t = t - L(rowIdx, k) * L(colIdx, k);
                end
                L(rowIdx, colIdx) = t / L(colIdx, colIdx);
            end
        else
            ok = false;              % Exception handler: not positive definite
        end
    end % ELSE: already failed, leave the rest zero
end

% padding block: diagonal only
for rowIdx = (nActive + 1):n
    if (ok) && (A(rowIdx, rowIdx) > 0.0) && isfinite(A(rowIdx, rowIdx))
        L(rowIdx, rowIdx) = sqrt(A(rowIdx, rowIdx));
    else
        ok = false;
    end
end

if (~ok)
    L = zeros(n, n);
end % ELSE is trivial

end
%------------------------------------------------------------------------
