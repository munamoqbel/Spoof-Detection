%******************************************************************************************
% DESCRIPTION:
% The innovation covariance as the gate's monitors read it: the host's S
% with the gate's covariance carried in. S = H P H' + R comes from the host
% with its own P; the gate reads that P floored (SPF_coastStartCov) and
% carries the receiver clock bias states at CPI_CLOCK_SIGMA, so the same
% increase is added to S: S + H dP H', dP the diagonal increase of the
% gate's P over the host's on every state and, on the clock states, of the
% clock floor's variance over the host's. A shift common to every range is
% a clock, by definition: with the clock carried free the CPI projection
% reads the position shift the rows agree on, where a clock held at the
% host's P after minutes of aiding made it read a common mode as position.
%
% INPUTS:
%   - innovationCov   [MAX_MEAS x MAX_MEAS] S from the host, nominal R; rows
%                     and columns beyond numMeas are padding
%   - obsMatrix       [MAX_MEAS x n] H
%   - numMeas         rows live this epoch
%   - covHost         [n x n] the host's P as it came
%   - covGate         [n x n] the same P as the gate reads it (floored)
%   - clockStates     [1 x k] the receiver clock bias states
%                     (CST_spfParam.CPI_CLOCK_STATES; a state outside the
%                     layout is ignored)
%   - clockSigma      the clock sigma S carries, m (CST_spfParam.CPI_CLOCK_SIGMA;
%                     0 carries none)
%
% OUTPUTS:
%   - covOut          S + H dP H' on the live rows, symmetric; the padding as
%                     it came
%
% ASSUMPTIONS AND LIMITATIONS:
% dP is diagonal: the floors raise diagonals only and the clock floor is a
% diagonal too, so the cross terms of the host's P stay in S as they came.
% With covGate equal to covHost and clockSigma 0 the matrix comes back as it
% went (the repository configuration). Element loops on the fixed MAX_MEAS
% layout (Coder).
%
% REQUIREMENT TRACEABILITY:
% docs/HOST_DECISIONS.md section 14
%
%******************************************************************************************
%#codegen
function [covOut] = SPF_floorInnovCov(innovationCov, obsMatrix, numMeas, covHost, covGate, clockStates, clockSigma)

% Define variables
numStates = size(obsMatrix, 2);
rows      = min(double(numMeas), double(size(innovationCov, 1)));
deltaP    = zeros(numStates, 1);
covOut    = innovationCov;

% the increase of the gate's P over the host's, per state
for idx = 1:numStates
    delta = covGate(idx, idx) - covHost(idx, idx);
    if (delta > 0.0)
        deltaP(idx) = delta;
    end
end

% the clock floor on the receiver clock bias states
clockVar = clockSigma * clockSigma;
for k = 1:numel(clockStates)
    idx = clockStates(k);
    if (idx >= 1) && (idx <= numStates)
        delta = clockVar - covHost(idx, idx);
        if (delta > deltaP(idx))
            deltaP(idx) = delta;
        end
    end
end

% S + H dP H' on the live rows
for rowIdx = 1:rows
    for colIdx = rowIdx:rows
        term = 0.0;
        for idx = 1:numStates
            if (deltaP(idx) > 0.0)
                term = term + obsMatrix(rowIdx, idx) * deltaP(idx) * obsMatrix(colIdx, idx);
            end
        end
        covOut(rowIdx, colIdx) = innovationCov(rowIdx, colIdx) + term;
        covOut(colIdx, rowIdx) = covOut(rowIdx, colIdx);
    end
end

end

%------------------------------------------------------------------------
