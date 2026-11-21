%******************************************************************************************
% DESCRIPTION:
% Doppler episode screen for the host's state update: are the range-rate
% rows of this epoch moving in a way the IMU did not measure? A receiver
% whose Doppler goes wrong on every row at once under strong transient
% dynamics is invisible to any consensus test, and a filter that follows
% those rows carries the error into the tilt, where a coast integrates it
% for as long as it lasts. The screen keeps a velocity coast anchored on
% the rows it has accepted and propagated with the IMU's own increment,
% and holds the range-rate rows out while they sit further from that
% coast than a tolerance allows.
%
% The host's innovation is total-state: predicted minus measured minus H
% times the prior error state. The velocity the rows are measured against
% is therefore vHat = velNav - priorVel, the mechanisation velocity less
% the prior state's velocity, and vHat does not move at a feedback, since
% a hand-over moves the mechanisation and the state by the same amount;
% it moves only by the IMU's increment and by the update's own step. Per
% call, in this order:
%   increment     dvIns = vHat - vHatPrev + incPrev: the change of vHat
%                 since the last call with the previous update's velocity
%                 step added back (the host writes screen.incPrev after
%                 each update), which is the IMU's increment less the
%                 error growth the filter's model expects; no hand-over
%                 enters and nothing about the feedback is assumed
%   rows          the range-rate rows (unit velocity partial) fitted alone
%                 to a velocity and one clock drift per constellation
%                 (normal equations on the host's matrixInv, as in
%                 SPF_rowScreen); vRows = vHat + ROWS_SIGN * dv is the
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
%                 episode is closed and the coast is re-anchored on vHat,
%                 the filter's own velocity, which has run without the
%                 range-rate rows for the whole hold and carries the ranges'
%                 aiding, a Doppler-free second opinion the pure-IMU coast
%                 is not; its tolerance is kept as grown, so rows within it
%                 are accepted again and rows still beyond it start a new
%                 episode; a long fault is held out for as long as the
%                 filter's velocity disagrees with it
%   blind         the IMU alone carries the coast for HOLD epochs, flagged
%                 or without rows alike: after HOLD consecutive calls
%                 without rows to judge, the coast is lost and the next rows
%                 anchor it afresh, so a tunnel does not leave a coast
%                 metres per second off to flag the honest rows at its exit.
%                 The hold is time: it advances on every call while an
%                 episode is on, judged or not, and the hold limit applies on
%                 an unjudged call too, since the re-anchor needs no rows; an
%                 unjudged call neither confirms nor releases
%
% The flag is consumed by the measurement noise alone: the range-rate rows'
% R goes to the host's maximum while it is true, the row layout stays, the
% ranges and the pressure row update as before, and the gate reads every
% row at the nominal R as it does for the host's own reject. A wrong
% Doppler row in a turn does not help the heading, it damages it; in an
% honest turn the rows and the IMU agree and nothing is flagged.
%
% INPUTS:
%   - screenIn        the screen's state from the previous call
%                     (STRUCT_SPF.dopplerScreenState to start); its incPrev
%                     is the velocity step the previous update added to the
%                     state, [1 x 3] NED m/s, written by the host after that
%                     update (zeros when none was applied); it is consumed
%                     here and zeroed
%   - velNav          [3 x 1] NED m/s, the mechanisation velocity the rows
%                     were predicted with, at the measurement time
%   - priorVel        [3 x 1] NED m/s, the velocity components of the prior
%                     error state this update starts from
%   - obsMatrix       [MAX_MEAS x n] H as the host built it
%   - innovation      [MAX_MEAS x 1] the rows' innovation, predicted minus
%                     measured minus H times the prior
%   - numRows         the GNSS rows live this epoch
%   - driftStates     [1 x k] receiver clock drift states, one per
%                     constellation (a state outside the layout is ignored)
%   - prm             the screen's constants (STRUCT_SPF.dopplerScreenParam):
%                     FLOOR (m/s), DRIFT_RATE (m/s per s), HEADING_SIGMA
%                     (rad), ALPHA, CONFIRM, RELEASE, RELEASE_GROW, HOLD
%                     (epochs), EPOCH_DT (s), MIN_RATE_ROWS, ROWS_SIGN
%
% OUTPUTS:
%   - screenOut       the state for the next call (the host then writes
%                     screenOut.incPrev with this update's velocity step)
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
% is, for up to HOLD such calls in a row; then the coast is lost and the
% episode with it. A call without a previous velocity loses the coast, which the next
% rows anchor again. A step the state took outside the update (a gate
% commit) is not in incPrev and reads as an IMU increment for that epoch.
% The horizontal plane only is judged. A slow drift of the rows under
% FLOOR / ALPHA per epoch is followed by the anchor and passes, by design.
% Element loops on the fixed MAX_MEAS layout (Coder); no persistent state,
% the host carries screenIn / screenOut.
%
% REQUIREMENT TRACEABILITY:
% docs/HOST_DECISIONS.md section 17
%
%******************************************************************************************
%#codegen
function [screenOut, dropRateRows, info] = SPF_dopplerScreen(screenIn, velNav, priorVel, obsMatrix, ...
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
vHat         = [velNav(1) - priorVel(1), velNav(2) - priorVel(2), velNav(3) - priorVel(3)];

% the IMU's increment since the last call: the change of nav less prior, the previous update's step added back
haveIns = false;
dvIns   = zeros(1, 3);
if (screenOut.havePrev)
    dvIns   = vHat - screenOut.vHatPrev + screenOut.incPrev;
    haveIns = true;
end
screenOut.vHatPrev = vHat;
screenOut.incPrev  = zeros(1, 3);                    % consumed; the host writes this update's step after it
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
    end
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
            end
        end
        if (used)
            numUnk = numUnk + 1;
            unkCols(numUnk) = col;
        end
    end
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
                end
            end
            normalMat(i, j) = term;
        end
        term = 0.0;
        for rowIdx = 1:rows
            if inRate(rowIdx)
                term = term + obsMatrix(rowIdx, unkCols(i)) * innovation(rowIdx);
            end
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
        vRows = vHat + prm.ROWS_SIGN * dv;
        haveRows = true;
    end
end

% the coast follows the IMU; without an increment it is lost
if (screenOut.anchored) && (haveIns)
    screenOut.vCoast     = screenOut.vCoast + dvIns;
    screenOut.tAnchor    = screenOut.tAnchor + prm.EPOCH_DT;
    screenOut.speedSince = screenOut.speedSince + sqrt(dvIns(1)^2 + dvIns(2)^2);
elseif (screenOut.anchored)
    screenOut.anchored  = false;
    screenOut.inEpisode = false;
    screenOut.above     = 0;
end

% the decision
if (haveRows)
    screenOut.blind = 0;
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
                end
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
                screenOut.inEpisode  = false;                        % the hold is up: the coast re-anchored on the filter's own velocity,
                screenOut.rel        = 0;                            % Doppler-free through the hold, its tolerance kept as grown; rows
                screenOut.hold       = 0;                            % within it are accepted again, rows still beyond it start a new episode
                screenOut.vCoast     = vHat;
            end
        end
        info.valid     = true;
        info.sep       = sep;
        info.tolerance = tolerance;
    end
else
    dropRateRows = screenOut.inEpisode && (numRate > 0);    % too few rate rows to judge: an episode in progress stands; none: nothing to drop
    if (screenOut.inEpisode)
        screenOut.hold = screenOut.hold + 1;                % the hold is time, judged or not
        if (screenOut.hold >= prm.HOLD)
            screenOut.inEpisode  = false;                   % the hold is up on an unjudged call: the same re-anchor on the filter's velocity
            screenOut.rel        = 0;
            screenOut.hold       = 0;
            screenOut.vCoast     = vHat;
            dropRateRows         = false;
        end
    end
    if (screenOut.anchored)
        screenOut.blind = screenOut.blind + 1;
        if (screenOut.blind > prm.HOLD)
            screenOut.anchored  = false;                    % blind for longer than the IMU alone can carry it: the coast is lost,
            screenOut.inEpisode = false;                    % the next rows anchor it afresh
            screenOut.above     = 0;
            screenOut.blind     = 0;
            dropRateRows        = false;
        end
    end
end
info.numRate   = numRate;
info.rateRows  = inRate;
info.vRows     = vRows;
info.dvIns     = dvIns;
info.inEpisode = screenOut.inEpisode;
info.hold      = screenOut.hold;

end
