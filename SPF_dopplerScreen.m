%******************************************************************************************
% DESCRIPTION:
% Doppler episode screen for the host's measurement setup: are the
% range-rate rows of this epoch moving in a way the IMU did not measure?
% A receiver whose Doppler goes wrong on every row at once under strong
% transient dynamics is invisible to any consensus test, and a filter
% that follows those rows carries the error into the tilt, where a coast
% integrates it for as long as it lasts. The screen keeps a velocity coast
% anchored on the rows it has accepted and propagated with the
% mechanisation's own increment, and holds the range-rate rows out while
% they sit further from that coast than a tolerance allows.
%
% Per call, in this order:
%   increment     dvIns = velNavPre - velPrev + velCorrApplied: the change
%                 of the mechanisation velocity since the last call with
%                 the velocity correction the mechanisation applied in
%                 between added back, so the increment is the IMU's alone
%   rows          the range-rate rows (unit velocity partial) fitted alone
%                 to a velocity and one clock drift per constellation
%                 (normal equations on the host's matrixInv, as in
%                 SPF_rowScreen); vRows = velNavPre + ROWS_SIGN * dv is the
%                 velocity the rows say, ROWS_SIGN -1 for an innovation of
%                 predicted minus measured
%   coast         vCoast follows dvIns every call; each accepted epoch
%                 re-anchors it on the rows through a first-order filter
%                 (ALPHA of the difference), so the anchor is the recent
%                 honest rows and not one epoch; the tolerance
%                 T = sqrt(FLOOR^2 + (DRIFT_RATE * tAnchor)^2
%                        + (HEADING_SIGMA * speedSince)^2)
%                 grows with the time the coast has run on the IMU alone,
%                 for the tilt and bias drift and a heading error on the
%                 speed change since the anchor
%   flag          CONFIRM consecutive epochs with |vRows - vCoast| > T
%                 (horizontal) start an episode; while it lasts the coast
%                 runs on the IMU and dropRateRows is true; RELEASE
%                 consecutive epochs back within T end it, plus one more
%                 for every RELEASE_GROW epochs the rows have been out, so
%                 rows that were wrong for long must agree for longer; the
%                 coast is kept at a release, so rows that only came half
%                 way back are flagged again at once; at HOLD epochs the
%                 episode is closed with the coast kept and its tolerance
%                 still growing, so rows still beyond it start a new
%                 episode and rows within it re-anchor; a long fault is held
%                 out until the grown tolerance accepts it
%
% The flag is consumed by the measurement setup alone: the range-rate rows'
% R goes to the host's maximum while it is true, the row layout stays, the
% ranges and the pressure row update as before, and the gate reads every
% row at the nominal R as it does for the host's own reject. A wrong
% Doppler row in a turn does not help the heading, it damages it; in an
% honest turn the rows and the IMU agree and nothing is flagged.
%
% INPUTS:
%   - screenIn        the screen's state from the previous call
%                     (STRUCT_SPF.dopplerScreenState to start)
%   - velNavPre       [3 x 1] NED m/s, the mechanisation velocity this
%                     update predicts with, before this epoch's feedback
%   - velCorrApplied  [3 x 1] NED m/s, the velocity correction the
%                     mechanisation applied since the previous call (the
%                     previous epoch's hand-over as applied; zeros when
%                     none), positive when it was taken off the velocity
%   - obsMatrix       [MAX_MEAS x n] H as the host built it
%   - innovation      [MAX_MEAS x 1] the rows' innovation
%   - numRows         the GNSS rows live this epoch
%   - driftStates     [1 x k] receiver clock drift states, one per
%                     constellation (a state outside the layout is ignored)
%   - prm             the screen's constants (STRUCT_SPF.dopplerScreenParam):
%                     FLOOR (m/s), DRIFT_RATE (m/s per s), HEADING_SIGMA
%                     (rad), ALPHA, CONFIRM, RELEASE, RELEASE_GROW, HOLD
%                     (epochs), EPOCH_DT (s), MIN_RATE_ROWS, ROWS_SIGN
%
% OUTPUTS:
%   - screenOut       the state for the next call
%   - dropRateRows    true: hold the range-rate rows out of this update
%   - info            valid (enough rate rows and a coast to judge
%                     against), numRate, rateRows ([MAX_MEAS x 1] logical,
%                     the rows read as range-rate rows, for the host to
%                     put at its maximum R when dropRateRows is true), sep,
%                     tolerance, vRows [1 x 3], dvIns [1 x 3], inEpisode,
%                     hold
%
% ASSUMPTIONS AND LIMITATIONS:
% Row types are read from H with the host's layout: velocity 4:6 in NED
% m/s; a range-rate row has a unit velocity partial larger than its
% position partial. Fewer than MIN_RATE_ROWS such rows, or an invalid
% inverse of the normal equations, leave the epoch undecided: the coast
% keeps propagating, no row is dropped, an episode in progress stays as it
% is. A call without a previous velocity, or a previous call that left no
% increment, loses the coast, which the next rows anchor again. The
% horizontal plane only is judged. A slow drift of the rows under
% FLOOR / ALPHA per epoch is followed by the anchor and passes, by design.
% Element loops on the fixed MAX_MEAS layout (Coder); no persistent state,
% the host carries screenIn / screenOut.
%
% REQUIREMENT TRACEABILITY:
% docs/HOST_DECISIONS.md section 17
%
%******************************************************************************************
%#codegen
function [screenOut, dropRateRows, info] = SPF_dopplerScreen(screenIn, velNavPre, velCorrApplied, obsMatrix, ...
    innovation, numRows, driftStates, prm)

% Define variables
MAX_UNKNOWNS = 10;                                   % 3 velocity axes + up to 7 drift columns
VEL_COLS     = [4 5 6];
UNIT_MIN     = 0.5;                                  % a unit partial reads above this
numMeasMax   = size(obsMatrix, 1);
numStates    = size(obsMatrix, 2);
rows         = min(double(numRows), double(numMeasMax));
screenOut    = screenIn;
dropRateRows = false;
info         = struct('valid', false, 'numRate', 0, 'sep', 0.0, 'tolerance', 0.0, 'vRows', zeros(1, 3), ...
    'dvIns', zeros(1, 3), 'inEpisode', false, 'hold', 0, 'rateRows', false(numMeasMax, 1));
velNow       = [velNavPre(1), velNavPre(2), velNavPre(3)];
corrNow      = [velCorrApplied(1), velCorrApplied(2), velCorrApplied(3)];

% the IMU's increment since the last call: the mechanisation's change with the applied correction added back
haveIns = false;
dvIns   = zeros(1, 3);
if (screenOut.havePrev)
    dvIns   = velNow - screenOut.velPrev + corrNow;
    haveIns = true;
end % ELSE: first call, no increment yet
screenOut.velPrev  = velNow;
screenOut.havePrev = true;

% the rate rows and their velocity: normal equations on the velocity axes and the drift columns in use
inRate = false(numMeasMax, 1);
numRate = 0;
for rowIdx = 1:rows
    posNorm = sqrt(obsMatrix(rowIdx, 1)^2 + obsMatrix(rowIdx, 2)^2 + obsMatrix(rowIdx, 3)^2);
    velNorm = sqrt(obsMatrix(rowIdx, 4)^2 + obsMatrix(rowIdx, 5)^2 + obsMatrix(rowIdx, 6)^2);
    if (velNorm > UNIT_MIN) && (velNorm > posNorm)
        inRate(rowIdx) = true;
        numRate = numRate + 1;
    end % ELSE: a range row, the pressure row, or padding
end
unkCols = zeros(1, MAX_UNKNOWNS);
unkCols(1:3) = VEL_COLS;
numUnk = 3;
for k = 1:numel(driftStates)
    col = driftStates(k);
    if (col >= 1) && (col <= numStates) && (numUnk < MAX_UNKNOWNS)
        used = false;
        for rowIdx = 1:rows
            if inRate(rowIdx) && (abs(obsMatrix(rowIdx, col)) > 1e-9)
                used = true;
            end % ELSE: this row does not carry the column
        end
        if (used)
            numUnk = numUnk + 1;
            unkCols(numUnk) = col;
        end % ELSE: no rate row uses this drift column
    end % ELSE: a state outside the layout, or no room left: ignored
end
haveRows = false;
vRows = zeros(1, 3);
if (numRate >= prm.MIN_RATE_ROWS) && (numRate > numUnk)
    normalMat = eye(MAX_UNKNOWNS);
    rhs = zeros(MAX_UNKNOWNS, 1);
    for i = 1:numUnk
        for j = 1:numUnk
            term = 0.0;
            for rowIdx = 1:rows
                if inRate(rowIdx)
                    term = term + obsMatrix(rowIdx, unkCols(i)) * obsMatrix(rowIdx, unkCols(j));
                end % ELSE: not a rate row
            end
            normalMat(i, j) = term;
        end
        term = 0.0;
        for rowIdx = 1:rows
            if inRate(rowIdx)
                term = term + obsMatrix(rowIdx, unkCols(i)) * innovation(rowIdx);
            end % ELSE: not a rate row
        end
        rhs(i) = term;
    end
    [normalInv, invInvalid] = matrixInv(normalMat);
    if (~invInvalid)
        dv = zeros(1, 3);
        for i = 1:3
            term = 0.0;
            for j = 1:numUnk
                term = term + normalInv(i, j) * rhs(j);
            end
            dv(i) = term;
        end
        vRows = velNow + prm.ROWS_SIGN * dv;
        haveRows = true;
    end % ELSE: no usable fit: the epoch is not judged
end % ELSE: too few rate rows for a fit with redundancy: the epoch is not judged

% the coast follows the IMU; without an increment it is lost
if (screenOut.anchored) && (haveIns)
    screenOut.vCoast     = screenOut.vCoast + dvIns;
    screenOut.tAnchor    = screenOut.tAnchor + prm.EPOCH_DT;
    screenOut.speedSince = screenOut.speedSince + sqrt(dvIns(1)^2 + dvIns(2)^2);
elseif (screenOut.anchored)
    screenOut.anchored  = false;
    screenOut.inEpisode = false;
    screenOut.above     = 0;
end % ELSE: no coast yet

% the decision
if (haveRows)
    if (~screenOut.anchored)
        screenOut.vCoast     = vRows;                           % the first rows anchor the coast
        screenOut.anchored   = true;
        screenOut.tAnchor    = 0.0;
        screenOut.speedSince = 0.0;
    else
        sep = sqrt((vRows(1) - screenOut.vCoast(1))^2 + (vRows(2) - screenOut.vCoast(2))^2);
        tolerance = sqrt(prm.FLOOR^2 + (prm.DRIFT_RATE * screenOut.tAnchor)^2 + (prm.HEADING_SIGMA * screenOut.speedSince)^2);
        if (~screenOut.inEpisode)
            if (sep > tolerance)
                screenOut.above = screenOut.above + 1;              % not accepted: the coast keeps to the IMU
                if (screenOut.above >= prm.CONFIRM)
                    screenOut.inEpisode = true;
                    screenOut.hold  = 0;
                    screenOut.rel   = 0;
                    screenOut.above = 0;
                    dropRateRows    = true;
                end % ELSE: one epoch beyond the tolerance is not yet an episode
            else
                screenOut.above      = 0;
                screenOut.vCoast     = screenOut.vCoast + prm.ALPHA * (vRows - screenOut.vCoast);   % accepted: re-anchored, filtered
                screenOut.tAnchor    = 0.0;
                screenOut.speedSince = 0.0;
            end
        else
            screenOut.hold = screenOut.hold + 1;
            dropRateRows = true;
            if (sep <= tolerance)
                screenOut.rel = screenOut.rel + 1;
            else
                screenOut.rel = 0;
            end
            if (screenOut.rel >= prm.RELEASE + floor(screenOut.hold / prm.RELEASE_GROW))
                screenOut.inEpisode  = false;                        % released: the coast kept, the rows must keep agreeing
                screenOut.rel        = 0;
                screenOut.hold       = 0;
                screenOut.vCoast     = screenOut.vCoast + prm.ALPHA * (vRows - screenOut.vCoast);
                screenOut.tAnchor    = 0.0;
                screenOut.speedSince = 0.0;
            elseif (screenOut.hold >= prm.HOLD)
                screenOut.inEpisode  = false;                        % the hold is up: the coast and its tolerance go on growing,
                screenOut.rel        = 0;                            % rows still beyond it start a new episode, rows within it re-anchor
                screenOut.hold       = 0;
            end % ELSE: the episode goes on
        end
        info.valid     = true;
        info.sep       = sep;
        info.tolerance = tolerance;
    end
else
    dropRateRows = screenOut.inEpisode;                     % no rate rows to judge: an episode in progress stands
end
info.numRate   = numRate;
info.rateRows  = inRate;
info.vRows     = vRows;
info.dvIns     = dvIns;
info.inEpisode = screenOut.inEpisode;
info.hold      = screenOut.hold;

end
