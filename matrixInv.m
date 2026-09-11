%******************************************************************************************
% DESCRIPTION:
% STAND-IN for the simulation's SVD-based matrix inverse. The target
% simulation has its own matrixInv with this exact signature (the one used
% for the S inverse inside kfUpdate); this file exists only so the
% standalone harness and tests run. Keep this file OFF the path when
% running inside the simulation.
%
% CONTRACT THE GATE RELIES ON (SPF_cpiMonitor, SPF_revalidation):
%   - the input is the full fixed MAX_MEAS x MAX_MEAS matrix: the live
%     block in the leading rows/cols, zero elsewhere EXCEPT the padding
%     diagonal, which the gate fills with the largest live variance so the
%     matrix is non-singular (blkdiag(S, p*I));
%   - the result inverts the live block and has no coupling between the
%     live block and the padding, so products with zero-padded vectors
%     equal the live-block results;
%   - invInvalid is true when no usable inverse exists (non-finite,
%     singular / rank-deficient input); the gate then ignores the result.
%
% INPUTS:
%   - matrix         [n x n]
%
% OUTPUTS:
%   - matrixInverse  [n x n]  pseudo-inverse (zeros where invalid)
%   - invInvalid     logical
%
%******************************************************************************************
%#codegen
function [matrixInverse, invInvalid] = matrixInv(matrix)

n = size(matrix, 1);
matrixInverse = zeros(n, n);
invInvalid = ~all(isfinite(matrix(:)));

if (~invInvalid)
    [U, S, V] = svd(matrix);
    singularValues = diag(S);
    tolerance = n * eps(max(singularValues));
    sInv = zeros(n, 1);
    for idx = 1:n
        if (singularValues(idx) > tolerance)
            sInv(idx) = 1.0 / singularValues(idx);
        end % ELSE: zero (padding / rank deficient) -> pseudo-inverse
    end
    matrixInverse = V * diag(sInv) * U';
    invInvalid = ~all(isfinite(matrixInverse(:)));
    if (invInvalid)
        matrixInverse = zeros(n, n);
    end % ELSE is trivial
end % ELSE: non-finite input

end
%------------------------------------------------------------------------
