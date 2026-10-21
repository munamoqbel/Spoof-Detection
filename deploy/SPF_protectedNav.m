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
%   level below ARM_PL_MAX; a level above restarts the count, a row
%   dropout pauses it; arming is sticky.
% - A latch falls back to the anchor only when the anchor is at most
%   MAX_ANCHOR_AGE old; otherwise the current solution is frozen and
%   info.anchorMissing is raised.
% - resetInhibit asks the host not to apply its protective reset while
%   the gate holds GNSS off or alarms; coastBudgetExceeded (coast PL
%   above ARM_PL_MAX, or COAST_BUDGET_EPOCHS passed) says the gate can
%   no longer certify a re-validation and the host may reset.
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

        % 1b. warm-up: a PL above ARM_PL_MAX restarts the count, a row dropout pauses it
        if (~sys.armed)
            sigmaPosMax = 0.0;
            for idx = 1:numel(CST_spfParam.MONITORED_AXES)
                axisIdx = CST_spfParam.MONITORED_AXES(idx);
                sigmaPosMax = max(sigmaPosMax, sqrt(max(kfMeas.postCov(axisIdx, axisIdx), 0.0)));
            end
            plInside   = (CST_spfParam.K_MISSED_DETECTION * sigmaPosMax < CST_spfParam.ARM_PL_MAX);   % Eq. 51 term
            enoughRows = (kfMeas.numMeas >= CST_spfParam.REVAL_MIN_MEAS);
            if (~plInside)
                sys.armCount = 0;
            elseif (enoughRows)
                sys.armCount = sys.armCount + 1;
            end
            if (sys.armCount >= double(CST_spfParam.ARM_EPOCHS))
                sys.armed = true;
            end
        end

        if (~sys.armed)
            % not armed: no windows, no alarms, the anchor follows the solution
            sys.anchor = STRUCT_SPF.setAnchor(true, zeros(numStates, 1), kfMeas.postCov, uint32(epoch));
            sys.pool   = STRUCT_SPF.closeAllWindows(sys.pool);
        else
            % 2. monitor bank on the host increments
            [sys.pool, report] = SPF_monitorPool(sys.pool, kfMeas, propTel);

            % 2b. the anchor kept live (E28, E31 in increment form)
            if (sys.anchor.valid)
                sys.anchor.separation = Phi * sys.anchor.separation + kfIncrement;
                sys.anchor.covariance = Phi * sys.anchor.covariance * Phi' + Q;
                sys.anchor.covariance = (sys.anchor.covariance + sys.anchor.covariance') / 2;
            end

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

                if (sys.anchor.valid) && (anchorAge <= CST_spfParam.MAX_ANCHOR_AGE)
                    % fall back to the certified-clean coast
                    applyCorrection       = true;
                    correction            = -sys.anchor.separation;
                    navState              = kfMeas.postState + correction;
                    sys.coastCov          = sys.anchor.covariance;
                    info.eventAnchorEpoch = sys.anchor.epoch;
                else
                    % no usable anchor: the current solution is frozen
                    info.eventAnchorEpoch = uint32(0);
                    info.anchorMissing    = true;
                end

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

        % 2. the raw GNSS against the coast; the operational KF is not updated, so prior = coast
        sys.probSep = zeros(numStates, 1);
        coastMinusPrior = zeros(numStates, 1);
        [passed, q_value, revalFault] = SPF_revalidation(kfMeas.innovation, kfMeas.obsMatrix, ...
            kfMeas.measNoiseCov, kfMeas.numMeas, coastMinusPrior, sys.coastCov);
        info.qReval = q_value;
        info.revalComputed = true;
        info.solveFault = revalFault;

        if (passed)
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
            sys.mode           = CST_spfMode.PROBATION;
            kfCommand = STRUCT_SPF.setCommand(true);
        end

    % ==================================================================
    case CST_spfMode.PROBATION
    % ==================================================================

        % 1. the protected output stays on the coast (E31)
        sys.coastCount = sys.coastCount + 1;
        sys.coastCov   = Phi * sys.coastCov * Phi' + Q;
        sys.coastCov   = (sys.coastCov + sys.coastCov') / 2;

        % 2. diagnostic: the GNSS against the coast on the trial's innovation
        coastMinusPrior = -(Phi * sys.probSep);
        [~, q_value, revalFault] = SPF_revalidation(kfMeas.innovation, kfMeas.obsMatrix, ...
            kfMeas.measNoiseCov, kfMeas.numMeas, coastMinusPrior, sys.coastCov);
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
