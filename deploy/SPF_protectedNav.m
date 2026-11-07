%***************************************************************************
% DESCRIPTION:
% One epoch of the spoofing-protected navigator, the FSM behind
% SPF_gate, on host-owned filters in increment form.
%   NOMINAL    the host KF is updated normally; the monitors watch its
%              increments; a clean-closed window refreshes the anchor;
%              any alarm latches: the output falls back to the anchor's
%              coast and the mode goes to COAST.
%   COAST      the host does not update its operational KF, so the
%              100 Hz extrapolation is the INS-only coast; the raw GNSS
%              rows are re-validated against that coast; a dwell of
%              consecutive passes opens PROBATION. Monitor silence never
%              ends COAST.
%   PROBATION  the host updates a trial copy; the monitors watch it; an
%              alarm vetoes back to COAST; a full quiet probation
%              commits: the trial becomes the solution and the anchor.
% The host contract (docs/HOST_2HZ_WIRING.m): every epoch kfMeas from
% the ACTIVE filter's update; on nav.applyCorrection (latch, commit)
% the host passes nav.state and nav.covar to its setKF as the
% operational KF's update result; on kfCommand.startTrial it starts the
% trial filter as a copy of the operational KF.
%
% INPUTS:
%   - sys      STRUCT_SPF.setSys, the persistent FSM state
%   - kfMeas   STRUCT_SPF.setKfMeas: y, S, H, R, numMeas, x_bar, x+, P+
%   - propTel  STRUCT_SPF.setPropTel: interval Phi and Q
%   - epoch    epoch counter for the event stamps and the anchor age
%
% OUTPUTS:
%   - sys      updated
%   - spoofTel STRUCT_SPF.setTel: .info (mode, alarms, PL, qReval,
%              dwell, events, faults, armed, resetInhibit,
%              coastBudgetExceeded), .kfCommand (startTrial), .nav
%              (applyCorrection, state, correction, covar, sigmaPosition)
%
% ASSUMPTIONS AND LIMITATIONS:
% - Warm-up: after a (re)initialisation the monitors arm only after
%   ARM_EPOCHS NOMINAL epochs with REVAL_MIN_MEAS rows and a protection
%   level below ARM_PL_MAX on which the monitors, run in shadow, would
%   not have alarmed; a level above or a shadow alarm restarts the count, a row
%   dropout pauses it; arming is sticky.
% - A latch falls back to the anchor only when the anchor is at most
%   MAX_ANCHOR_AGE old; otherwise the current solution is frozen and
%   info.anchorMissing is raised.
% - resetInhibit asks the host not to apply its protective reset while
%   the gate holds GNSS off or alarms; the re-validation and the first
%   probation window use the coast covariance bounded by COAST_SIGMA_MAX
%   (STRUCT_SPF.boundCoastCov), and the dwell also needs the CPI of the
%   innovations against the coast below T_N, so a spoof the monitors
%   would latch on stays refused however long the coast; coastBudgetExceeded (coast PL
%   above ARM_PL_MAX, or COAST_BUDGET_EPOCHS passed) says the gate can
%   no longer certify a re-validation and the host may reset.
% - The return path: RETURN_BANK long INS-only coasts are kept beside the
%   anchor, restarted in turn (the oldest spans 80 to 120 s). At a latch,
%   and on every COAST epoch, the position the rows ask for is read
%   against the oldest: when the nav had wandered more than RETURN_WANDER
%   from it and the GNSS lands within RETURN_RADIUS of it, the event is
%   the end of a followed spoof, not an onset (eventReturn): the nav is
%   corrected back to the long coast, velocity included (a correction in
%   COAST goes through the host's setKF as a latch does), and the bound's
%   floor becomes that coast's own sigma, so the clean GNSS re-validates;
%   the first probation window keeps the short bound. Otherwise nothing
%   changes.
%
% REQUIREMENT TRACEABILITY:
% Paper Sec. 5 recovery: detection and latch, re-validation (Guide
% Sec. 5 chi-square gate), probation, hand-back; E28 and E31 for every
% coast propagation; Eq. 51 protection level for the warm-up and the
% budget (Implementation Guide rule 4).
%***************************************************************************
%#codegen
function [sys, spoofTel] = SPF_protectedNav(sys, kfMeas, propTel, epoch)

info      = STRUCT_SPF.zeroInfo;
kfCommand = STRUCT_SPF.zeroCommand;
numStates = CST_gnssHybrid.NO_STATES;
applyCorrection = false;
correction      = zeros(numStates, 1);    % increment relative to the operational KF
navState        = kfMeas.postState;       % x+ handed to setKF on a latch or a commit

Phi = propTel.accumPhi;
Q   = propTel.accumQ;
kfIncrement = kfMeas.postState - kfMeas.priorState;      % K y of the active filter

switch sys.mode
    % ==================================================================
    case CST_spfMode.NOMINAL
    % ==================================================================

        % 1. the protected solution is the host KF
        sys.coastCov = kfMeas.postCov;
        sys.probSep  = zeros(numStates, 1);

        % 1b. warm-up: a PL above ARM_PL_MAX restarts the count, a row dropout pauses it, and the monitors run in
        %     shadow (windows open and close as in NOMINAL, the margins logged, nothing latches and the alarm flags stay
        %     down): a shadow alarm restarts the count, so the gate arms only after ARM_EPOCHS epochs on which it
        %     would not have alarmed
        armedBefore = sys.armed;
        if (~sys.armed)
            sigmaPosMax = 0.0;
            for idx = 1:numel(CST_spfParam.MONITORED_AXES)
                axisIdx = CST_spfParam.MONITORED_AXES(idx);
                sigmaPosMax = max(sigmaPosMax, sqrt(max(kfMeas.postCov(axisIdx, axisIdx), 0.0)));
            end
            plInside   = (CST_spfParam.K_MISSED_DETECTION * sigmaPosMax < CST_spfParam.ARM_PL_MAX);   % Eq. 51 term
            enoughRows = (kfMeas.numMeas >= CST_spfParam.REVAL_MIN_MEAS);
            [sys.pool, report] = SPF_monitorPool(sys.pool, kfMeas, propTel);
            info.maxProtectionLevel = report.maxProtectionLevel;
            info.ssRatio            = report.ssRatio;
            info.cpiRatio           = report.cpiRatio;
            info.solveFault         = report.solveFault;
            if (report.anyAlarm)
                sys.armCount = 0;
                sys.pool     = STRUCT_SPF.closeAllWindows(sys.pool);
            else
                if (~plInside)
                    sys.armCount = 0;
                elseif (enoughRows)
                    sys.armCount = sys.armCount + 1;
                end
                sys.pool = STRUCT_SPF.openWindow(sys.pool, kfMeas);
            end
            if (sys.armCount >= double(CST_spfParam.ARM_EPOCHS))
                sys.armed = true;
            end
        end

        if (~armedBefore)
            % a warm-up epoch, the arming one included: no latch, the anchor follows the solution, the long coasts
            % start with the arming
            sys.anchor    = STRUCT_SPF.setAnchor(true, zeros(numStates, 1), kfMeas.postCov, uint32(epoch));
            sys.longCoast = STRUCT_SPF.zeroLongCoast;
        else
            % 2. monitor bank on the host increments
            [sys.pool, report] = SPF_monitorPool(sys.pool, kfMeas, propTel);

            % 2b. the anchor kept live (E28, E31 in increment form)
            if (sys.anchor.valid)
                sys.anchor.separation = Phi * sys.anchor.separation + kfIncrement;
                sys.anchor.covariance = Phi * sys.anchor.covariance * Phi' + Q;
                sys.anchor.covariance = (sys.anchor.covariance + sys.anchor.covariance') / 2;
            end
            % 2c. the long coasts of the return path kept live the same way, restarted in turn
            sys.longCoast = SPF_longCoastStep(sys.longCoast, Phi, Q, kfIncrement, kfMeas.postCov, uint32(epoch), true);

            info.ssAlarm            = report.ssAlarm;
            info.cpiAlarm           = report.cpiAlarm;
            info.alarmPerAxis       = report.alarmPerAxis;
            info.maxProtectionLevel = report.maxProtectionLevel;
            info.solveFault         = report.solveFault;
            info.ssRatio            = report.ssRatio;
            info.cpiRatio           = report.cpiRatio;

            % 3. an alarm-free clean close refreshes the anchor
            if (report.cleanCloseFound) && (~report.anyAlarm)
                sys.anchor.valid      = true;
                sys.anchor.separation = report.cleanCloseSeparation;
                sys.anchor.covariance = report.cleanCloseCovar;
                sys.anchor.epoch      = uint32(epoch);
            end

            % 4. any alarm latches, otherwise the next window opens
            if (report.anyAlarm)
                info.eventLatched = true;
                anchorAge = uint32(epoch) - sys.anchor.epoch;
                % a return or an onset? the GNSS lands where the inertial sensors alone would have carried the vehicle
                % over the last minutes while the nav had wandered from there: the end of a followed spoof
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
                    % back to the long coast, velocity included; the coast's uncertainty is that coast's own
                    applyCorrection       = true;
                    correction            = -sys.longCoast.separation(:, jLong);
                    navState              = kfMeas.postState + correction;
                    sys.coastCov          = sys.longCoast.covariance(:, :, jLong);
                    info.eventAnchorEpoch = sys.longCoast.epoch(jLong);
                    info.eventReturn      = true;
                    sys.longCoast         = STRUCT_SPF.zeroLongCoast;
                elseif (sys.anchor.valid) && (anchorAge <= CST_spfParam.MAX_ANCHOR_AGE)
                    % fall back to the certified-clean coast; the long coasts follow the nav
                    applyCorrection       = true;
                    correction            = -sys.anchor.separation;
                    navState              = kfMeas.postState + correction;
                    sys.coastCov          = sys.anchor.covariance;
                    info.eventAnchorEpoch = sys.anchor.epoch;
                    sys.longCoast         = SPF_longCoastShift(sys.longCoast, correction);
                else
                    % no usable anchor: the current solution is frozen
                    info.eventAnchorEpoch = uint32(0);
                    info.anchorMissing    = true;
                end
                sys.sigmaLatch = SPF_sigmaOf(sys.coastCov);   % the floor of the bound: the coast's own sigma at the latch

                sys.pool       = STRUCT_SPF.closeAllWindows(sys.pool);
                sys.dwellCount = 0;
                sys.coastCount = 0;
                sys.mode       = CST_spfMode.COAST;

            else
                sys.pool = STRUCT_SPF.openWindow(sys.pool, kfMeas);
            end

        end

    % ==================================================================
    case CST_spfMode.COAST
    % ==================================================================

        % 1. the protected output is the INS-only coast (E31)
        sys.coastCount = sys.coastCount + 1;
        sys.coastCov   = Phi * sys.coastCov * Phi' + Q;
        sys.coastCov   = (sys.coastCov + sys.coastCov') / 2;

        % 1b. the long coasts kept live (no update in COAST: no increment), and the return read here too: the clean
        %     GNSS back after a followed spoof reaches the gate only once the host's spike filter has released its rows
        sys.longCoast = SPF_longCoastStep(sys.longCoast, Phi, Q, zeros(numStates, 1), kfMeas.postCov, uint32(epoch), false);
        [isReturn, jLong, info] = SPF_returnCheck(sys.longCoast, kfMeas, zeros(numStates, 1), zeros(numStates, 1), info);
        coastMinusPrior = zeros(numStates, 1);
        if (isReturn)
            % back to the long coast, velocity included; the coast's uncertainty is that coast's own; the operational
            % KF is the coast here, so the correction is relative to the prior
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
        end
        coastCovBound = STRUCT_SPF.boundCoastCov(sys.coastCov, sys.sigmaLatch);        % the coast as the re-admission tests may use it
        coastCovOpen  = STRUCT_SPF.boundCoastCov(sys.coastCov, zeros(numStates, 1));   % the first probation window keeps the short bound

        % 2. the raw GNSS against the coast; the operational KF is not updated, so prior = coast (the new one on a return)
        sys.probSep = zeros(numStates, 1);
        [passed, q_value, revalFault] = SPF_revalidation(kfMeas.innovation, kfMeas.obsMatrix, ...
            kfMeas.measNoiseCov, kfMeas.numMeas, coastMinusPrior, coastCovBound);
        info.qReval = q_value;
        info.revalComputed = true;
        info.solveFault = revalFault;
        % 2b. the same innovations through the CPI over the last N epochs, S with the bounded coast covariance and a zero
        %     increment: GNSS the monitors would latch on is not re-admitted (B and T_N as in NOMINAL), whatever the chi-square says
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
            sys.dwellCount = 0;                      % a fail restarts the dwell; too few rows to test holds it
        end

        % 3. dwell full: PROBATION opens, the host starts the trial as a copy of the coasting KF
        if (sys.dwellCount >= CST_spfParam.REVAL_DWELL_REQUIRED)
            info.eventProbationStarted = true;
            sys.probationCount = 0;
            sys.dwellCount     = 0;
            sys.probSep        = zeros(numStates, 1);
            sys.pool           = STRUCT_SPF.closeAllWindows(sys.pool);
            % the first probation window opens on the coast with the bounded covariance, so the trial's first
            % update counts in its separation (otherwise a trial that jumps onto a spoof in one update is seen by nothing)
            kfMeasOpen         = kfMeas;
            kfMeasOpen.postCov = coastCovOpen;
            sys.pool           = STRUCT_SPF.openWindow(sys.pool, kfMeasOpen);
            sys.mode           = CST_spfMode.PROBATION;
            kfCommand = STRUCT_SPF.setCommand(true);
        end

    % ==================================================================
    case CST_spfMode.PROBATION
    % ==================================================================

        % 1. the protected output stays on the coast (E31); the long coasts kept live
        sys.coastCount = sys.coastCount + 1;
        sys.coastCov   = Phi * sys.coastCov * Phi' + Q;
        sys.coastCov   = (sys.coastCov + sys.coastCov') / 2;
        sys.longCoast  = SPF_longCoastStep(sys.longCoast, Phi, Q, zeros(numStates, 1), kfMeas.postCov, uint32(epoch), false);
        coastCovBound  = STRUCT_SPF.boundCoastCov(sys.coastCov, sys.sigmaLatch);   % the coast as the re-admission tests may use it

        % 2. diagnostic: the GNSS against the coast on the trial's innovation
        coastMinusPrior = -(Phi * sys.probSep);
        [~, q_value, revalFault] = SPF_revalidation(kfMeas.innovation, kfMeas.obsMatrix, ...
            kfMeas.measNoiseCov, kfMeas.numMeas, coastMinusPrior, coastCovBound);
        info.qReval = q_value;
        info.revalComputed = true;
        info.solveFault = revalFault;

        % 3. trial minus coast accumulates the trial's increments (E28)
        sys.probSep = Phi * sys.probSep + kfIncrement;

        % 4. monitor bank on the trial filter; it never refreshes the anchor
        [sys.pool, report] = SPF_monitorPool(sys.pool, kfMeas, propTel);

        info.ssAlarm            = report.ssAlarm;
        info.cpiAlarm           = report.cpiAlarm;
        info.alarmPerAxis       = report.alarmPerAxis;
        info.maxProtectionLevel = report.maxProtectionLevel;
        info.solveFault         = info.solveFault || report.solveFault;
        info.ssRatio            = report.ssRatio;
        info.cpiRatio           = report.cpiRatio;

        if (report.anyAlarm)
            % veto: the trial is discarded
            info.eventProbationVetoed = true;
            sys.pool       = STRUCT_SPF.closeAllWindows(sys.pool);
            sys.dwellCount = 0;
            sys.probSep    = zeros(numStates, 1);
            sys.mode       = CST_spfMode.COAST;
        else
            sys.pool = STRUCT_SPF.openWindow(sys.pool, kfMeas);

            sys.probationCount = sys.probationCount + 1;

            if (sys.probationCount >= CST_spfParam.PROBATION_LENGTH)
                % commit (hand-back): the host solution becomes the trial, which becomes the anchor
                info.eventHandback = true;
                applyCorrection    = true;
                correction         = sys.probSep;
                navState           = kfMeas.postState;
                sys.coastCov       = kfMeas.postCov;

                sys.anchor.valid      = true;
                sys.anchor.separation = zeros(numStates, 1);
                sys.anchor.covariance = kfMeas.postCov;
                sys.anchor.epoch      = uint32(epoch);
                sys.longCoast         = STRUCT_SPF.zeroLongCoast;   % the long coasts start again on the new solution
                sys.sigmaLatch        = zeros(numStates, 1);

                sys.probSep    = zeros(numStates, 1);
                sys.coastCount = 0;
                sys.mode       = CST_spfMode.NOMINAL;
            end
        end
end

% the navigation output of this epoch
sigmaPosition = zeros(3, 1);
for idx = 1:3
    axisIdx = CST_spfParam.MONITORED_AXES(idx);
    sigmaPosition(idx) = sqrt(max(sys.coastCov(axisIdx, axisIdx), 0.0));
end
nav = STRUCT_SPF.setNav(applyCorrection, navState, correction, sys.coastCov, sigmaPosition);

info.mode           = sys.mode;
info.dwellCount     = sys.dwellCount;
info.coastEpochs    = sys.coastCount;
info.numMeasClamped = kfMeas.numMeasClamped;
info.armed          = sys.armed;
% reset arbitration: the host must not reset while the gate holds GNSS off or alarms
info.resetInhibit        = (sys.mode ~= CST_spfMode.NOMINAL) || info.ssAlarm || info.cpiAlarm;
info.coastBudgetExceeded = (CST_spfParam.K_MISSED_DETECTION * max(sigmaPosition) ...
    > CST_spfParam.ARM_PL_MAX) ...      % coast PL (Eq. 51 term) beyond the limit: the host may reset
    || (sys.coastCount > double(CST_spfParam.COAST_BUDGET_EPOCHS));   % time belt (rule 4)

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
    end
end
if (inNominal) && (epoch >= longCoast.restartAt)
    idx = double(longCoast.next);
    longCoast.separation(:, idx)    = zeros(numStates, 1);
    longCoast.covariance(:, :, idx) = postCov;
    longCoast.epoch(idx)            = epoch;
    longCoast.next                  = uint8(mod(idx, nb) + 1);
    longCoast.restartAt             = epoch + CST_spfParam.RETURN_RESTART;
end
end

function [longCoast] = SPF_longCoastShift(longCoast, correction)
% the nav moved by correction: every started long coast's separation (nav minus coast) moves with it
nb = double(CST_spfParam.RETURN_BANK);
for idx = 1:nb
    if (longCoast.epoch(idx) > 0)
        longCoast.separation(:, idx) = longCoast.separation(:, idx) + correction;
    end
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
        end
    end
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
    end
end
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
        end
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
    end
end
end

function [sigma] = SPF_sigmaOf(covariance)
% the square root of the diagonal, floored at zero
numStates = CST_gnssHybrid.NO_STATES;
sigma = zeros(numStates, 1);
for idx = 1:numStates
    sigma(idx) = sqrt(max(covariance(idx, idx), 0.0));
end
end
