%******************************************************************************************
% DESCRIPTION:
% One epoch of the spoofing-protected navigator FSM (NOMINAL / COAST / PROBATION),
% host-owned filters, increment form. See docs/HOST_2HZ_WIRING.m for the host contract.
%
% INPUTS:  sys (STRUCT_SPF.setSys), kfMeas (STRUCT_SPF.setKfMeas), propTel, epoch
% OUTPUTS: sys (updated), spoofTel (STRUCT_SPF.setTel)
%******************************************************************************************
%#codegen
function [sys, spoofTel] = SPF_protectedNav(sys, kfMeas, propTel, epoch)

info      = STRUCT_SPF.zeroInfo;
kfCommand = STRUCT_SPF.zeroCommand;
numStates = CST_gnssHybrid.NO_STATES;
applyCorrection = false;
correction      = zeros(numStates, 1);    % increment relative to the operational KF
navState        = kfMeas.postState;       % x+ to hand to setKF on LATCH / COMMIT

Phi = propTel.accumPhi;
Q   = propTel.accumQ;
kfIncrement = kfMeas.postState - kfMeas.priorState;      % K*y of the active filter

switch sys.mode
    % ==================================================================
    case CST_spfMode.NOMINAL
    % ==================================================================

        % 1. the protected solution is the host KF
        sys.coastCov = kfMeas.postCov;
        sys.probSep  = zeros(numStates, 1);

        % 1b. warm-up: PL above ARM_PL_MAX restarts the count, a row dropout pauses it
        if (~sys.armed)
            sigmaPosMax = 0.0;
            for idx = 1:numel(CST_spfParam.MONITORED_AXES)
                axisIdx = CST_spfParam.MONITORED_AXES(idx);
                sigmaPosMax = max(sigmaPosMax, sqrt(max(kfMeas.postCov(axisIdx, axisIdx), 0.0)));
            end
            plInside   = (CST_spfParam.K_MISSED_DETECTION * sigmaPosMax < CST_spfParam.ARM_PL_MAX);
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
            % not armed: no windows, no alarms, anchor follows the solution
            sys.anchor = STRUCT_SPF.setAnchor(true, zeros(numStates, 1), kfMeas.postCov, uint32(epoch));
            sys.pool   = STRUCT_SPF.closeAllWindows(sys.pool);
        else
            % 2. monitor bank on the host increments
            [sys.pool, report] = SPF_monitorPool(sys.pool, kfMeas, propTel);

            % 2b. keep the anchor live
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

            % 3. refresh the anchor on an alarm-free clean close
            if (report.cleanCloseFound) && (~report.anyAlarm)
                sys.anchor.valid      = true;
                sys.anchor.separation = report.cleanCloseSeparation;
                sys.anchor.covariance = report.cleanCloseCovar;
                sys.anchor.epoch      = uint32(epoch);
            end

            % 4. latch on any alarm, else open the next window
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
                    % no usable anchor: freeze the current solution
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

        % 1. protected output = INS-only coast
        sys.coastCount = sys.coastCount + 1;
        sys.coastCov   = Phi * sys.coastCov * Phi' + Q;
        sys.coastCov   = (sys.coastCov + sys.coastCov') / 2;

        % 2. test the raw GNSS against the coast (operational KF not updated: prior = coast)
        sys.probSep = zeros(numStates, 1);
        coastMinusPrior = zeros(numStates, 1);
        [passed, q_value, revalFault] = SPF_revalidation(kfMeas.innovation, kfMeas.obsMatrix, ...
            kfMeas.measNoiseCov, kfMeas.numMeas, coastMinusPrior, sys.coastCov);
        info.qReval = q_value;
        info.revalComputed = true;
        info.solveFault = revalFault;

        if (passed)
            sys.dwellCount = sys.dwellCount + 1;
        else
            sys.dwellCount = 0;
        end

        % 3. dwell full: open PROBATION
        if (sys.dwellCount >= CST_spfParam.REVAL_DWELL_REQUIRED)
            info.eventProbationStarted = true;
            sys.probationCount = 0;
            sys.dwellCount     = 0;
            sys.probSep        = zeros(numStates, 1);
            sys.pool           = STRUCT_SPF.closeAllWindows(sys.pool);
            sys.mode           = CST_spfMode.PROBATION;
            kfCommand = STRUCT_SPF.setCommand(true);          % host: trial <- copy of the coasting KF
        end

    % ==================================================================
    case CST_spfMode.PROBATION
    % ==================================================================

        % 1. protected output stays on the coast
        sys.coastCount = sys.coastCount + 1;
        sys.coastCov   = Phi * sys.coastCov * Phi' + Q;
        sys.coastCov   = (sys.coastCov + sys.coastCov') / 2;

        % 2. diagnostic: GNSS vs coast on the trial's innovation
        coastMinusPrior = -(Phi * sys.probSep);
        [~, q_value, revalFault] = SPF_revalidation(kfMeas.innovation, kfMeas.obsMatrix, ...
            kfMeas.measNoiseCov, kfMeas.numMeas, coastMinusPrior, sys.coastCov);
        info.qReval = q_value;
        info.revalComputed = true;
        info.solveFault = revalFault;

        % 3. trial - coast accumulates the trial's increments
        sys.probSep = Phi * sys.probSep + kfIncrement;

        % 4. monitor bank on the trial filter (never refreshes the anchor)
        [sys.pool, report] = SPF_monitorPool(sys.pool, kfMeas, propTel);

        info.ssAlarm            = report.ssAlarm;
        info.cpiAlarm           = report.cpiAlarm;
        info.alarmPerAxis       = report.alarmPerAxis;
        info.maxProtectionLevel = report.maxProtectionLevel;
        info.solveFault         = info.solveFault || report.solveFault;
        info.ssRatio            = report.ssRatio;
        info.cpiRatio           = report.cpiRatio;

        if (report.anyAlarm)
            % VETO: discard the trial
            info.eventProbationVetoed = true;
            sys.pool       = STRUCT_SPF.closeAllWindows(sys.pool);
            sys.dwellCount = 0;
            sys.probSep    = zeros(numStates, 1);
            sys.mode       = CST_spfMode.COAST;
        else
            sys.pool = STRUCT_SPF.openWindow(sys.pool, kfMeas);

            sys.probationCount = sys.probationCount + 1;

            if (sys.probationCount >= CST_spfParam.PROBATION_LENGTH)
                % COMMIT (handback): host solution <- trial, new anchor
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

% navigation output of this epoch
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
% protective-reset arbitration: the host must not reset while the gate holds GNSS off or alarms
info.resetInhibit        = (sys.mode ~= CST_spfMode.NOMINAL) || info.ssAlarm || info.cpiAlarm;
info.coastBudgetExceeded = (sys.coastCount > double(CST_spfParam.COAST_BUDGET_EPOCHS));

spoofTel = STRUCT_SPF.setTel(info, kfCommand, nav);

end
%------------------------------------------------------------------------
