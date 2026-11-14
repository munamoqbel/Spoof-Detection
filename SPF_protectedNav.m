%******************************************************************************************
% DESCRIPTION:
% One epoch of the spoofing-protected navigator (the FSM behind
% SPF_gate). HOST-OWNED FILTERS, INCREMENT FORM.
%
% MODES (constants in CST_spfMode)
%   NOMINAL    host KF updated normally; monitors watch its increments;
%              clean-closed windows refresh the anchor; any alarm => LATCH.
%   COAST      host does NOT update its operational KF (100 Hz extrapolation
%              makes it the INS-only coast); the update is still computed
%              on a scratch copy to feed y, S here; re-validation counts
%              consecutive passes (an epoch with too few rows to test holds
%              the count); full dwell => PROBATION.
%              Monitor silence NEVER exits this mode (capture).
%   PROBATION  host updates a TRIAL copy; monitors watch the trial; alarm
%              => VETO; full quiet probation => COMMIT (handback).
%
% HOST CONTRACT (see docs/HOST_2HZ_WIRING.m)
%   - every epoch: kfMeas from the ACTIVE filter's update (operational KF
%     in NOMINAL/COAST, trial in PROBATION)
%   - nav.applyCorrection true (LATCH, COMMIT): pass (nav.state, nav.covar)
%     to the normal setKF as this epoch's update result of the OPERATIONAL
%     KF (LATCH: x+ of the operational KF minus the anchor separation;
%     COMMIT: the trial's x+, which by construction equals the coasting
%     KF.states plus the accumulated trial increments)
%   - kfCommand.startTrial true (probation opens): trial <- copy of the
%     operational KF (states and covariance)
%
% INPUTS:
%   - sys      STRUCT_SPF.setSys (persistent)
%   - kfMeas   STRUCT_SPF.setKfMeas (y, S, H, R, numMeas, x_bar, x+, P+)
%   - propTel  STRUCT_SPF.setPropTel (.accumPhi, .accumQ)
%   - epoch    scalar epoch counter (event stamps, anchor age)
%
% OUTPUTS:
%   - sys      updated
%   - spoofTel STRUCT_SPF.setTel (.info, .kfCommand, .nav)
%
% ASSUMPTIONS AND LIMITATIONS:
%
% REQUIREMENT TRACEABILITY:
%
%******************************************************************************************
%#codegen
function [sys, spoofTel] = SPF_protectedNav(sys, kfMeas, propTel, epoch)

info      = STRUCT_SPF.zeroInfo;
kfCommand = STRUCT_SPF.zeroCommand;
numStates = CST_gnssHybrid.NO_STATES;
applyCorrection = false;
correction      = zeros(numStates, 1);    % increment relative to the OPERATIONAL KF
navState        = kfMeas.postState;       % absolute x+ to hand to setKF (LATCH / COMMIT)

propTel = SPF_coastPropTel(propTel);                     % the gate's own Q on its coasts (CST_spfParam.COAST_Q_DATASHEET)
postCovHost    = kfMeas.postCov;                         % the host's P as it came: handed back as it came at a commit
kfMeas.postCov = SPF_coastStartCov(kfMeas.postCov);      % the host's P as the gate reads it: floored where its own
                                                         % consistency reads say it underestimates (CST_spfParam.COAST_SIGMA_FLOOR_*)
Phi = propTel.accumPhi;
Q   = propTel.accumQ;
kfIncrement = kfMeas.postState - kfMeas.priorState;      % K*y of the active filter

switch sys.mode
    % ==================================================================
    case CST_spfMode.NOMINAL
    % ==================================================================

        % ---- 1. the protected solution IS the host KF ----
        sys.coastCov = kfMeas.postCov;
        sys.probSep  = zeros(numStates, 1);

        % ---- 1b. warm-up: arm the monitors only on a converged filter ----
        %          ARM_EPOCHS qualifying epochs (enough rows and a protection
        %          level below ARM_PL_MAX). A protection level above the limit
        %          restarts the count (filter diverging or coasting); a short
        %          row dropout with the PL still inside only pauses it.
        %          The monitors run in SHADOW meanwhile: windows open and
        %          close as in NOMINAL and the margins are logged, but nothing
        %          latches and the alarm flags stay down (unarmed: no alarm,
        %          no latch, the host contract). A shadow alarm restarts the
        %          count, so the gate arms only after ARM_EPOCHS epochs on
        %          which it would not have alarmed: a transient the host's
        %          covariance does not admit (a heading still converging, a
        %          constellation joining the solution, a reacquisition)
        %          delays the arming instead of latching right after it
        %          (measured on the vehicle data: a latch at 83 s, 7 s after
        %          the arming, on the start transient of 82.6 to 84.6 s).
        %          Sticky once armed; the pool carries over into NOMINAL.
        armedBefore = sys.armed;
        if (~sys.armed)
            sigmaPosMax = 0.0;
            for idx = 1:numel(CST_spfParam.MONITORED_AXES)
                axisIdx = CST_spfParam.MONITORED_AXES(idx);
                sigmaPosMax = max(sigmaPosMax, sqrt(max(kfMeas.postCov(axisIdx, axisIdx), 0.0)));
            end
            plInside   = (CST_spfParam.K_MISSED_DETECTION * sigmaPosMax < CST_spfParam.ARM_PL_MAX);
            enoughRows = (kfMeas.numMeas >= CST_spfParam.REVAL_MIN_MEAS);
            [sys.pool, report] = SPF_monitorPool(sys.pool, kfMeas, propTel);
            info.maxProtectionLevel = report.maxProtectionLevel;
            info.ssRatio            = report.ssRatio;
            info.cpiRatio           = report.cpiRatio;
            info.solveFault         = report.solveFault;
            if (report.anyAlarm)
                sys.armCount = 0;                     % a shadow alarm: restart
                sys.pool     = STRUCT_SPF.closeAllWindows(sys.pool);
            else
                if (~plInside)
                    sys.armCount = 0;                 % restart
                elseif (enoughRows)
                    sys.armCount = sys.armCount + 1;  % qualifying epoch
                end % ELSE: dropout with the PL inside: pause, keep the count
                sys.pool = STRUCT_SPF.openWindow(sys.pool, kfMeas);
            end
            if (sys.armCount >= double(CST_spfParam.ARM_EPOCHS))
                sys.armed = true;
            end % ELSE is trivial
        end % ELSE: already armed

        if (~armedBefore)
            % a warm-up epoch, the arming one included (its pool pass is done
            % above): no latch; the anchor follows the current solution so it
            % is fresh when the monitors start; the long coasts of the return
            % path start with the arming
            sys.anchor    = STRUCT_SPF.setAnchor(true, zeros(numStates, 1), kfMeas.postCov, uint32(epoch));
            sys.longCoast = STRUCT_SPF.zeroLongCoast;
        else
            % ---- 2. monitor bank on the host increments ----
            [sys.pool, report] = SPF_monitorPool(sys.pool, kfMeas, propTel);

            % ---- 2b. keep the anchor LIVE: (host solution now) - (anchor coast
            %          now) accumulates this epoch's increment; P_C propagates ----
            if (sys.anchor.valid)
                sys.anchor.separation = Phi * sys.anchor.separation + kfIncrement;
                sys.anchor.covariance = Phi * sys.anchor.covariance * Phi' + Q;
                sys.anchor.covariance = (sys.anchor.covariance + sys.anchor.covariance') / 2;
            end % ELSE is trivial

            % ---- 2c. the long coasts of the return path, kept live the same
            %          way and restarted in turn: the memory the anchor lacks ----
            sys.longCoast = SPF_longCoastStep(sys.longCoast, Phi, Q, kfIncrement, kfMeas.postCov, uint32(epoch), true);

            info.ssAlarm            = report.ssAlarm;
            info.cpiAlarm           = report.cpiAlarm;
            info.alarmPerAxis       = report.alarmPerAxis;
            info.maxProtectionLevel = report.maxProtectionLevel;
            info.solveFault         = report.solveFault;
            info.ssRatio            = report.ssRatio;
            info.cpiRatio           = report.cpiRatio;

            % ---- 3. refresh anchor (only on alarm-free epochs, so the
            %         anchor always ends strictly before detection) ----
            if (report.cleanCloseFound) && (~report.anyAlarm)
                sys.anchor.valid      = true;
                sys.anchor.separation = report.cleanCloseSeparation;   % already at 'now'
                sys.anchor.covariance = report.cleanCloseCovar;
                sys.anchor.epoch      = uint32(epoch);
            end

            % ---- 4. latch on any alarm, else open the next window ----
            if (report.anyAlarm)
                info.eventLatched = true;
                anchorAge = uint32(epoch) - sys.anchor.epoch;
                % A return or an onset? "The GNSS is 30 m from the nav" reads
                % both ways; the inertial sensors tell them apart: in a
                % followed spoof the nav has wandered from where the INS alone
                % would have carried it, and the returning GNSS lands on that
                % INS coast within its drift; an onset lands away from both.
                % The wander is the anchor's (the nav before the event): a
                % latch that comes a few epochs late, once the KF has already
                % moved towards the returning GNSS, still reads as a return.
                if (sys.anchor.valid)
                    refSep = sys.anchor.separation;
                else
                    refSep = kfIncrement;                          % no anchor: the nav before this update
                end
                [isReturn, jLong, info] = SPF_returnCheck(sys.longCoast, kfMeas, kfIncrement, refSep, info);

                if (isReturn)
                    % back to the long coast, velocity included (the ramp's
                    % velocity goes with the position); the coast's
                    % uncertainty is that coast's own, the floor of the bound
                    applyCorrection       = true;
                    correction            = -sys.longCoast.separation(:, jLong);
                    navState              = kfMeas.postState + correction;
                    sys.coastCov          = sys.longCoast.covariance(:, :, jLong);
                    info.eventAnchorEpoch = sys.longCoast.epoch(jLong);
                    info.eventReturn      = true;
                    sys.longCoast         = STRUCT_SPF.zeroLongCoast;
                elseif (sys.anchor.valid) && (anchorAge <= CST_spfParam.MAX_ANCHOR_AGE)
                    % fall back to the certified-clean coast: the host must
                    % move its solution by -(accumulated increments since the
                    % anchor window opened) and take the coast covariance;
                    % the long coasts follow the nav through the correction
                    applyCorrection       = true;
                    correction            = -sys.anchor.separation;
                    navState              = kfMeas.postState + correction;   % operational x+ - separation
                    sys.coastCov          = sys.anchor.covariance;
                    info.eventAnchorEpoch = sys.anchor.epoch;
                    sys.longCoast         = SPF_longCoastShift(sys.longCoast, correction);
                else
                    % no usable anchor: freeze the current solution
                    info.eventAnchorEpoch = uint32(0);
                    info.anchorMissing    = true;
                end
                % the floor of the bound the re-admission tests use: the
                % coast's own sigma at the latch (3 m or so from the anchor,
                % the long coast's own on a return)
                sys.sigmaLatch = SPF_sigmaOf(sys.coastCov);

                sys.pool       = STRUCT_SPF.closeAllWindows(sys.pool);
                sys.dwellCount = 0;
                sys.coastCount = 0;
                sys.mode       = CST_spfMode.COAST;

            else
                sys.pool = STRUCT_SPF.openWindow(sys.pool, kfMeas);
            end

        end % armed

    % ==================================================================
    case CST_spfMode.COAST
    % ==================================================================

        % ---- 1. protected output = INS-only coast (host extrapolates the
        %         operational KF the same way; GNSS severed) ----
        sys.coastCount = sys.coastCount + 1;
        sys.coastCov   = Phi * sys.coastCov * Phi' + Q;
        sys.coastCov   = (sys.coastCov + sys.coastCov') / 2;

        % ---- 1b. the long coasts kept live (no update in COAST: no increment)
        %          and the return read here as well: the clean GNSS back after
        %          a followed spoof reaches the gate only once the host's spike
        %          filter has released its rows, and a followed drag can end
        %          while the gate already coasts on a later onset ----
        sys.longCoast = SPF_longCoastStep(sys.longCoast, Phi, Q, zeros(numStates, 1), kfMeas.postCov, uint32(epoch), false);
        [isReturn, jLong, info] = SPF_returnCheck(sys.longCoast, kfMeas, zeros(numStates, 1), zeros(numStates, 1), info);
        coastMinusPrior = zeros(numStates, 1);
        if (isReturn)
            % back to the long coast, velocity included; the coast's
            % uncertainty is that coast's own. The operational KF IS the coast
            % here, so the correction is relative to the prior; the host
            % applies it through its setKF as it does a latch
            applyCorrection       = true;
            correction            = -sys.longCoast.separation(:, jLong);
            navState              = kfMeas.priorState + correction;
            sys.coastCov          = sys.longCoast.covariance(:, :, jLong);
            sys.sigmaLatch        = SPF_sigmaOf(sys.coastCov);
            info.eventReturn      = true;
            info.eventAnchorEpoch = sys.longCoast.epoch(jLong);
            sys.longCoast         = STRUCT_SPF.zeroLongCoast;
            sys.dwellCount        = 0;
            sys.pool              = STRUCT_SPF.closeAllWindows(sys.pool);
            coastMinusPrior       = correction;             % this epoch's rows are read against the new coast
        end % ELSE: the coast stands
        coastCovBound = STRUCT_SPF.boundCoastCov(sys.coastCov, sys.sigmaLatch);        % the coast as the re-admission tests may use it
        coastCovOpen  = STRUCT_SPF.boundCoastCov(sys.coastCov, zeros(numStates, 1));   % the first probation window keeps the short bound

        % ---- 2. test the (untrusted) GNSS against the coast ----
        %         HOST CONTRACT: the operational KF is not updated in COAST
        %         (the scratch update only feeds y, S here), so its prior IS
        %         the coast: x_coast - x_prior = 0 and the residual is y
        %         (on a return the new coast minus the prior is the correction).
        sys.probSep = zeros(numStates, 1);
        [passed, q_value, revalFault] = SPF_revalidation(kfMeas.innovation, kfMeas.obsMatrix, ...
            kfMeas.measNoiseCov, kfMeas.numMeas, coastMinusPrior, coastCovBound);
        info.qReval = q_value;
        info.revalComputed = true;
        info.solveFault = revalFault;
        % ---- 2b. the same innovations through the CPI over the last N epochs ----
        %         The chi-square is one epoch at a time at 1e-3 on all rows; on
        %         the vehicle data it admits an offset below about 24 m even with
        %         the bounded coast covariance, while the CPI latches from 20 m.
        %         So the dwell also requires the CPI statistic of the innovations
        %         against the coast (a zero increment: no update in COAST, the SS
        %         stays silent; S = H P_C,bound H' + R) to stay below T_N: GNSS
        %         the monitors would latch on is not re-admitted, with the same B
        %         and T_N as in NOMINAL. The windows of the pool serve as the
        %         sliding window; they are closed at the latch and at the
        %         probation open, and none of them touches the anchor here.
        coastS = kfMeas.obsMatrix * coastCovBound * kfMeas.obsMatrix' + kfMeas.measNoiseCov;
        kfMeasCoast = kfMeas;
        kfMeasCoast.innovation    = kfMeas.innovation - kfMeas.obsMatrix * coastMinusPrior;
        kfMeasCoast.innovationCov = (coastS + coastS') / 2;
        kfMeasCoast.postState     = kfMeas.priorState;
        [sys.pool, report] = SPF_monitorPool(sys.pool, kfMeasCoast, propTel);
        sys.pool = STRUCT_SPF.openWindow(sys.pool, kfMeasCoast);
        info.cpiRatio   = report.cpiRatio;             % the CPI against the coast, logged through the coast
        info.solveFault = info.solveFault || report.solveFault;
        quiet = ~report.cpiAlarm;

        if (passed && quiet)
            sys.dwellCount = sys.dwellCount + 1;
        elseif (double(kfMeas.numMeas) >= double(CST_spfParam.REVAL_MIN_MEAS))
            sys.dwellCount = 0;                      % a failed test restarts the dwell
        end % ELSE: too few rows to test (no GNSS, baro only): no evidence either way, the dwell holds

        % ---- 3. dwell full: open PROBATION (never direct handback) ----
        if (sys.dwellCount >= CST_spfParam.REVAL_DWELL_REQUIRED)
            info.eventProbationStarted = true;
            sys.probationCount = 0;
            sys.dwellCount     = 0;
            sys.probSep        = zeros(numStates, 1);   % trial starts ON the coast
            sys.pool           = STRUCT_SPF.closeAllWindows(sys.pool);
            % The first probation window opens on the COAST, separation zero,
            % covariance the bounded coast covariance: the trial's first update
            % then counts in that window's separation. Without it a trial that
            % jumps onto a spoof in one update (a wide coast covariance against
            % 14 rows at 5 m absorbs 85 % at once) is seen by no window, since
            % the windows of the probation open on the trial's post-update
            % solution, and the probation only watches what the trial does after.
            kfMeasOpen         = kfMeas;
            kfMeasOpen.postCov = coastCovOpen;   % the short bound whatever opened the probation: a trial that jumps
            sys.pool           = STRUCT_SPF.openWindow(sys.pool, kfMeasOpen);   % more than about 20 m is vetoed
            sys.mode           = CST_spfMode.PROBATION;
            % order the host to start its TRIAL as a copy of the coasting KF
            kfCommand = STRUCT_SPF.setCommand(true);
        end

    % ==================================================================
    case CST_spfMode.PROBATION
    % ==================================================================

        % ---- 1. protected output STAYS on the coast; the long coasts kept live ----
        sys.coastCount = sys.coastCount + 1;
        sys.coastCov   = Phi * sys.coastCov * Phi' + Q;
        sys.coastCov   = (sys.coastCov + sys.coastCov') / 2;
        sys.longCoast  = SPF_longCoastStep(sys.longCoast, Phi, Q, zeros(numStates, 1), kfMeas.postCov, uint32(epoch), false);
        coastCovBound  = STRUCT_SPF.boundCoastCov(sys.coastCov, sys.sigmaLatch);   % the coast as the re-admission tests may use it

        % ---- 2. diagnostic: GNSS-vs-coast on the trial's innovation ----
        coastMinusPrior = -(Phi * sys.probSep);
        [~, q_value, revalFault] = SPF_revalidation(kfMeas.innovation, kfMeas.obsMatrix, ...
            kfMeas.measNoiseCov, kfMeas.numMeas, coastMinusPrior, coastCovBound);
        info.qReval = q_value;
        info.revalComputed = true;
        info.solveFault = revalFault;

        % ---- 3. trial - coast accumulates the trial's increments ----
        sys.probSep = Phi * sys.probSep + kfIncrement;

        % ---- 4. monitor bank on THE TRIAL filter ----
        [sys.pool, report] = SPF_monitorPool(sys.pool, kfMeas, propTel);

        info.ssAlarm            = report.ssAlarm;
        info.cpiAlarm           = report.cpiAlarm;
        info.alarmPerAxis       = report.alarmPerAxis;
        info.maxProtectionLevel = report.maxProtectionLevel;
        info.solveFault         = info.solveFault || report.solveFault;
        info.ssRatio            = report.ssRatio;
        info.cpiRatio           = report.cpiRatio;
        % probation windows never refresh the anchor: the trial is not
        % yet trusted, so report.cleanClose* is deliberately ignored.

        if (report.anyAlarm)
            % ---- VETO: monitors caught the trial; discard it ----
            info.eventProbationVetoed = true;
            sys.pool       = STRUCT_SPF.closeAllWindows(sys.pool);
            sys.dwellCount = 0;
            sys.probSep    = zeros(numStates, 1);   % operational KF is still the coast
            sys.mode       = CST_spfMode.COAST;
        else
            sys.pool = STRUCT_SPF.openWindow(sys.pool, kfMeas);

            sys.probationCount = sys.probationCount + 1;

            if (sys.probationCount >= CST_spfParam.PROBATION_LENGTH)
                % ---- COMMIT (handback): host solution <- trial ----
                info.eventHandback = true;
                applyCorrection    = true;
                correction         = sys.probSep;          % trial - coast (relative to the operational KF)
                navState           = kfMeas.postState;     % the trial's x+ IS the clean solution
                sys.coastCov       = postCovHost;          % the host takes its own P back, not the floored one

                % just certified by a full quiet probation: new anchor; the
                % long coasts start again on the new solution
                sys.anchor.valid      = true;
                sys.anchor.separation = zeros(numStates, 1);
                sys.anchor.covariance = kfMeas.postCov;
                sys.anchor.epoch      = uint32(epoch);
                sys.longCoast         = STRUCT_SPF.zeroLongCoast;
                sys.sigmaLatch        = zeros(numStates, 1);

                sys.probSep    = zeros(numStates, 1);
                sys.coastCount = 0;
                sys.mode       = CST_spfMode.NOMINAL;
                % pool carries over seamlessly
            end
        end
end

% ---- navigation output of this epoch ----
sigmaPosition = zeros(3, 1);
for idx = 1:3
    axisIdx = CST_spfParam.MONITORED_AXES(idx);
    sigmaPosition(idx) = sqrt(max(sys.coastCov(axisIdx, axisIdx), 0.0));
end
nav = STRUCT_SPF.setNav(applyCorrection, navState, correction, sys.coastCov, sigmaPosition);

info.mode           = sys.mode;
info.dwellCount     = sys.dwellCount;
info.coastEpochs    = sys.coastCount;    % time-in-coast (COAST + PROBATION), epochs
info.numMeasClamped = kfMeas.numMeasClamped;
info.armed          = sys.armed;
% protective-reset arbitration: while the gate holds GNSS off (COAST,
% PROBATION) or alarms this epoch, a GNSS-vs-INS discrepancy is its finding
% and the host must not reset on it; the coast budget bounds the hold:
% once the coast's own position uncertainty (K_MD * sigma_C, the PL of an
% INS-only solution) exceeds ARM_PL_MAX the gate could no longer certify a
% re-validation, so the host may reset. Grows with Phi/Q, not with time.
% COAST_BUDGET_EPOCHS is the belt behind it: a Q too small keeps sigma_C
% under the limit however long the hold, so a long time bound closes it.
info.resetInhibit        = (sys.mode ~= CST_spfMode.NOMINAL) || info.ssAlarm || info.cpiAlarm;
info.coastBudgetExceeded = (CST_spfParam.K_MISSED_DETECTION * max(sigmaPosition) ...
    > CST_spfParam.ARM_PL_MAX) ...      % coast too uncertain to re-validate: host may reset
    || (sys.coastCount > double(CST_spfParam.COAST_BUDGET_EPOCHS));   % time belt (Q too small)

spoofTel = STRUCT_SPF.setTel(info, kfCommand, nav);

end

%------------------------------------------------------------------------

function [longCoast] = SPF_longCoastStep(longCoast, Phi, Q, increment, postCov, epoch, inNominal)
% the long coasts kept live (E28, E31 in increment form) and, in NOMINAL, restarted in turn every RETURN_RESTART epochs
numStates = CST_gnssHybrid.NO_STATES;
nb = double(CST_spfParam.RETURN_BANK);
for idx = 1:nb
    if (longCoast.epoch(idx) > 0)
        longCoast.separation(:, idx) = Phi * longCoast.separation(:, idx) + increment;
        covNext = Phi * longCoast.covariance(:, :, idx) * Phi' + Q;
        longCoast.covariance(:, :, idx) = (covNext + covNext') / 2;
    end % ELSE: not started
end
if (inNominal) && (epoch >= longCoast.restartAt)
    idx = double(longCoast.next);
    longCoast.separation(:, idx)    = zeros(numStates, 1);
    longCoast.covariance(:, :, idx) = postCov;
    longCoast.epoch(idx)            = epoch;
    longCoast.next                  = uint8(mod(idx, nb) + 1);
    longCoast.restartAt             = epoch + CST_spfParam.RETURN_RESTART;
end % ELSE: no restart this epoch
end

function [longCoast] = SPF_longCoastShift(longCoast, correction)
% the nav moved by correction: every started long coast's separation (nav minus coast) moves with it
nb = double(CST_spfParam.RETURN_BANK);
for idx = 1:nb
    if (longCoast.epoch(idx) > 0)
        longCoast.separation(:, idx) = longCoast.separation(:, idx) + correction;
    end % ELSE: not started
end
end

function [isReturn, jLong, info] = SPF_returnCheck(longCoast, kfMeas, increment, refSep, info)
% the return check against the oldest long coast: the horizontal wander from it of the nav before the event (the
% separation refSep: the anchor's at a latch, zero for the coasting nav in COAST), and where the position the rows ask
% for lands relative to it; a return when the wander exceeds RETURN_WANDER and the landing is within RETURN_RADIUS
% (the GNSS agrees with the inertial sensors while the nav does not)
isReturn = false;
jLong    = 1;
found    = false;
oldest   = uint32(0);
nb = double(CST_spfParam.RETURN_BANK);
for idx = 1:nb
    if (longCoast.epoch(idx) > 0)
        if (~found) || (longCoast.epoch(idx) < oldest)
            oldest = longCoast.epoch(idx);
            jLong  = idx;
            found  = true;
        end % ELSE: a younger one
    end % ELSE: not started
end
if (found)
    [dHat, dHatOk] = SPF_positionAsked(kfMeas);
    if (dHatOk)
        wander2  = 0.0;
        landing2 = 0.0;
        for idx = 1:2                                            % the horizontal axes
            axisIdx   = CST_spfParam.MONITORED_AXES(idx);
            sepBefore = longCoast.separation(axisIdx, jLong) - increment(axisIdx);   % nav minus the long coast before this update
            sepRef    = longCoast.separation(axisIdx, jLong) - refSep(axisIdx);      % the nav before the event minus the long coast
            wander2   = wander2 + sepRef * sepRef;
            landing2  = landing2 + (dHat(idx) + sepBefore) * (dHat(idx) + sepBefore);
        end
        info.returnChecked = true;
        info.returnWander  = sqrt(wander2);
        info.returnLanding = sqrt(landing2);
        isReturn = (info.returnWander >= CST_spfParam.RETURN_WANDER) && (info.returnLanding <= CST_spfParam.RETURN_RADIUS);
    end % ELSE: too few rows or an unusable S: no reading
end % ELSE: no long coast yet
end

function [dHat, ok] = SPF_positionAsked(kfMeas)
% the position offset the rows ask for on the monitored axes, f' S^-1 gamma / f' S^-1 f per axis (the CPI
% projection with the full S); ok false with fewer than REVAL_MIN_MEAS rows or an unusable S
maxMeas = double(CST_spfParam.MAX_MEAS);
numMeas = double(kfMeas.numMeas);
dHat = zeros(1, 3);
ok   = false;
if (numMeas >= double(CST_spfParam.REVAL_MIN_MEAS))
    innovationCov = kfMeas.innovationCov;
    padValue = max(max(diag(innovationCov)), 1.0);
    varianceOk = true;
    for rowIdx = 1:maxMeas
        if (rowIdx > numMeas)
            innovationCov(rowIdx, rowIdx) = padValue;             % padding diagonal: non-singular for matrixInv
        elseif ~(innovationCov(rowIdx, rowIdx) > 0.0)
            varianceOk = false;
        end % ELSE: a live row with a positive variance
    end
    [sInverse, invInvalid] = matrixInv(innovationCov);
    if (varianceOk) && (~invInvalid)
        sInvInnovation = sInverse * kfMeas.innovation;
        ok = true;
        for idx = 1:3
            projection     = kfMeas.obsMatrix(:, CST_spfParam.MONITORED_AXES(idx));
            sInvProjection = sInverse * projection;
            gammaProjection  = 0.0;
            sigma2Projection = 0.0;
            for rowIdx = 1:maxMeas                                % element loop: deterministic rounding
                gammaProjection  = gammaProjection + projection(rowIdx) * sInvInnovation(rowIdx);
                sigma2Projection = sigma2Projection + projection(rowIdx) * sInvProjection(rowIdx);
            end
            if (sigma2Projection > 0.0) && isfinite(sigma2Projection) && isfinite(gammaProjection)
                dHat(idx) = gammaProjection / sigma2Projection;
            else
                ok = false;
            end
        end
    end % ELSE: S unusable
end % ELSE: too few rows
end

function [sigma] = SPF_sigmaOf(covariance)
% the square root of the diagonal, floored at zero
numStates = CST_gnssHybrid.NO_STATES;
sigma = zeros(numStates, 1);
for idx = 1:numStates
    sigma(idx) = sqrt(max(covariance(idx, idx), 0.0));
end
end
