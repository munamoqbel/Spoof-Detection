%******************************************************************************************
% DESCRIPTION:
% Row screen for the host's measurement setup: which rows of an epoch
% disagree with the OTHER rows of their type, judged against the rows'
% own fit, not against the filter's prediction. The ranges are fitted
% alone to a position shift and one clock per constellation, the range
% rates alone to a velocity and one clock drift per constellation; the
% worst row whose residual against its fit exceeds SCREEN_K times the
% residual's own sigma (nominal sigma times sqrt(1 - leverage), Baarda's
% w-test, so a low satellite that pulls the fit toward itself is judged
% on what it hides) is dropped and the fit repeated, while the fit keeps
% a row of redundancy. A kept row gets
% the nominal R whatever its misclosure against the filter; a dropped row
% gets the nominal R times DROP_FACTOR, which takes it out of the update
% without changing the row layout.
%
% Why the reference matters. A rule that compares each row with the
% filter's prediction cannot tell a bad row from a filter that is wrong,
% and treats both alike: it shrinks the row's weight until the
% disagreement counts as one sigma. When the INS carries an error, the
% rows that see it are exactly the rows such a rule weights down, in
% proportion to the error, so the larger the error the less the filter
% is allowed to learn about it; a heading error then survives the
% manoeuvres that should correct it. On the vehicle data (long run,
% braking from 35 to 9 m/s under 21 mrad of heading error) the honest
% Doppler rows sat at four to six sigma against the filter, were weighted
% down 25 to 42 times, and the north velocity error grew to 0.55 m/s
% against fifteen rows that pointed at it every epoch. Against the rows'
% own fit the same rows agreed to 0.02 to 0.10 m/s, within noise: they
% were right and the filter was wrong. A single satellite's multipath or
% cycle slip, by contrast, stands out against the other rows and is still
% caught here. Under a spoof all rows move together and pass: following
% consistent rows is the host's job, deciding whether they are honest is
% the gate's, and the gate reads the nominal R by design.
%
% INPUTS:
%   - obsMatrix     [MAX_MEAS x n] H as the host built it
%   - misclosure    [MAX_MEAS x 1] measured minus predicted per row (the
%                   misclosure the host tests today)
%   - measNoiseIn   [MAX_MEAS x 1] the diagonal of R as the host's setup
%                   left it
%   - numRows       the GNSS rows live this epoch (the pressure row, when
%                   it sits beyond them, is not touched)
%   - overrideFlag  true: the host's override is on, nothing is changed
%   - clockStates   [1 x k] receiver clock bias states, one per
%                   constellation, read on the range rows (a state outside
%                   the layout is ignored)
%   - driftStates   [1 x k] receiver clock drift states, read on the
%                   range-rate rows
%   - sigmaNominal  [1 x 2] nominal sigma of a range row and of a range-rate
%                   row, m and m/s (the host's SQRT_R_RANGE, SQRT_R_DOPPLER)
%   - screenK       the threshold on a row's residual against the fit, in
%                   nominal sigma (3)
%   - dropFactor    a dropped row's R is sigma^2 times this (1e6)
%
% OUTPUTS:
%   - measNoiseOut  [MAX_MEAS x 1] the diagonal of R to use: nominal on
%                   every kept range and range-rate row, nominal times
%                   DROP_FACTOR on a dropped one, measNoiseIn elsewhere
%   - numDropped    [1 x 2] rows dropped, ranges and range rates
%   - droppedFlag   [MAX_MEAS x 1] logical, true on a dropped row
%   - fitShift      [2 x 3] the fit's position shift of the range rows (m)
%                   and velocity of the range-rate rows (m/s), NED, from
%                   the last fit of each type; zeros when no fit was made
%
% ASSUMPTIONS AND LIMITATIONS:
% Row types are read from H with the host's layout: position 1:3, velocity
% 4:6 in NED metres and m/s. A range row has a unit position partial on at
% least two axes and no velocity partial; a range-rate row a unit velocity
% partial. A row of neither kind, the pressure row among them, is left as
% it came. A constellation with a single row has its own clock unknown, so
% that row can never be dropped; a type with fewer rows than unknowns plus
% one is not screened. The normal equations (at most MAX_UNKNOWNS) are
% solved with the host's matrixInv on a fixed-size padded matrix; an
% invalid inverse ends the screening of that type with every row kept.
% Element loops on the fixed MAX_MEAS layout (Coder).
%
% REQUIREMENT TRACEABILITY:
% docs/HOST_DECISIONS.md section 16
%
%******************************************************************************************
%#codegen
function [measNoiseOut, numDropped, droppedFlag, fitShift] = SPF_rowScreen(obsMatrix, misclosure, measNoiseIn, ...
    numRows, overrideFlag, clockStates, driftStates, sigmaNominal, screenK, dropFactor)

% Define variables
MAX_UNKNOWNS = 10;                                   % 3 axes + up to 7 clock columns
POS_COLS     = [1 2 3];
VEL_COLS     = [4 5 6];
UNIT_MIN     = 0.5;                                  % a unit partial reads above this
ZERO_MAX     = 0.1;                                  % and an absent one below this
LEVERAGE_MIN = 1.0e-6;                               % 1 - h below this: the row alone sets an unknown, no test
numMeasMax   = size(obsMatrix, 1);
numStates    = size(obsMatrix, 2);
rows         = min(double(numRows), double(numMeasMax));
measNoiseOut = measNoiseIn;
numDropped   = zeros(1, 2);
droppedFlag  = false(numMeasMax, 1);
fitShift     = zeros(2, 3);

if (~overrideFlag)
    for rowType = 1:2
        % the rows of this type and the columns of its unknowns
        if (rowType == 1)
            axisCols = POS_COLS;
            sysCols  = clockStates;
        else
            axisCols = VEL_COLS;
            sysCols  = driftStates;
        end
        sigma   = sigmaNominal(rowType);
        inType  = false(numMeasMax, 1);
        for rowIdx = 1:rows
            posNorm = sqrt(obsMatrix(rowIdx, 1)^2 + obsMatrix(rowIdx, 2)^2 + obsMatrix(rowIdx, 3)^2);
            velNorm = sqrt(obsMatrix(rowIdx, 4)^2 + obsMatrix(rowIdx, 5)^2 + obsMatrix(rowIdx, 6)^2);
            posEntries = double(abs(obsMatrix(rowIdx, 1)) > 1e-9) + double(abs(obsMatrix(rowIdx, 2)) > 1e-9) + ...
                double(abs(obsMatrix(rowIdx, 3)) > 1e-9);
            if (rowType == 1)
                inType(rowIdx) = (posNorm > UNIT_MIN) && (velNorm < ZERO_MAX) && (posEntries >= 2);
            else
                inType(rowIdx) = (velNorm > UNIT_MIN) && (velNorm > posNorm);
            end
        end
        % the unknown columns: the three axes, then every clock column some row of the type uses
        unkCols = zeros(1, MAX_UNKNOWNS);
        unkCols(1:3) = axisCols;
        numUnk = 3;
        for k = 1:numel(sysCols)
            col = sysCols(k);
            if (col >= 1) && (col <= numStates) && (numUnk < MAX_UNKNOWNS)
                used = false;
                for rowIdx = 1:rows
                    if inType(rowIdx) && (abs(obsMatrix(rowIdx, col)) > 1e-9)
                        used = true;
                    end
                end
                if (used)
                    numUnk = numUnk + 1;
                    unkCols(numUnk) = col;
                end
            end
        end
        % drop the worst row beyond SCREEN_K sigma of the fit, refit, while a row of redundancy remains
        kept = inType;
        numKept = 0;
        for rowIdx = 1:rows
            if kept(rowIdx)
                numKept = numKept + 1;
            end
        end
        screening = (numKept - numUnk) >= 1;
        passes = 0;
        while (screening) && (passes < rows)
            passes = passes + 1;
            % normal equations on the kept rows, padded with the identity beyond numUnk
            normalMat = eye(MAX_UNKNOWNS);
            rhs = zeros(MAX_UNKNOWNS, 1);
            for i = 1:numUnk
                for j = 1:numUnk
                    term = 0.0;
                    for rowIdx = 1:rows
                        if kept(rowIdx)
                            term = term + obsMatrix(rowIdx, unkCols(i)) * obsMatrix(rowIdx, unkCols(j));
                        end
                    end
                    normalMat(i, j) = term;
                end
                term = 0.0;
                for rowIdx = 1:rows
                    if kept(rowIdx)
                        term = term + obsMatrix(rowIdx, unkCols(i)) * misclosure(rowIdx);
                    end
                end
                rhs(i) = term;
            end
            [normalInv, invInvalid] = matrixInv(normalMat);
            if (invInvalid)
                screening = false;                                      % no usable fit: every row of the type is kept
            else
                solution = zeros(MAX_UNKNOWNS, 1);
                for i = 1:numUnk
                    term = 0.0;
                    for j = 1:numUnk
                        term = term + normalInv(i, j) * rhs(j);
                    end
                    solution(i) = term;
                end
                fitShift(rowType, 1) = solution(1);
                fitShift(rowType, 2) = solution(2);
                fitShift(rowType, 3) = solution(3);
                % the worst standardised residual against the fit: the residual over its own sigma,
                % sigma * sqrt(1 - h), h the row's leverage A (A'A)^-1 A', so that a row which pulls the
                % fit toward itself (a low satellite on the vertical) is judged on what it hides, not on
                % what it leaves (Baarda's w-test)
                worstRow = 0;
                worstNorm = 0.0;
                for rowIdx = 1:rows
                    if kept(rowIdx)
                        fitted = 0.0;
                        leverage = 0.0;
                        for i = 1:numUnk
                            fitted = fitted + obsMatrix(rowIdx, unkCols(i)) * solution(i);
                            for j = 1:numUnk
                                leverage = leverage + obsMatrix(rowIdx, unkCols(i)) * normalInv(i, j) * obsMatrix(rowIdx, unkCols(j));
                            end
                        end
                        residualVar = 1.0 - leverage;
                        if (residualVar > LEVERAGE_MIN)
                            residNorm = abs(misclosure(rowIdx) - fitted) / (sigma * sqrt(residualVar));
                        else
                            residNorm = 0.0;                              % a row that alone sets an unknown cannot be tested
                        end
                        if (residNorm > worstNorm)
                            worstNorm = residNorm;
                            worstRow = rowIdx;
                        end
                    end
                end
                if (worstRow >= 1) && (worstNorm > screenK)
                    kept(worstRow) = false;
                    droppedFlag(worstRow) = true;
                    numDropped(rowType) = numDropped(rowType) + 1;
                    numKept = numKept - 1;
                    screening = (numKept - numUnk) >= 1;              % another pass only with redundancy left
                else
                    screening = false;                                % every remaining row within SCREEN_K of the fit
                end
            end
        end
        % the R to use on the rows of this type
        for rowIdx = 1:rows
            if inType(rowIdx)
                if droppedFlag(rowIdx)
                    measNoiseOut(rowIdx) = sigma * sigma * dropFactor;
                else
                    measNoiseOut(rowIdx) = sigma * sigma;
                end
            end
        end
    end
end

end
