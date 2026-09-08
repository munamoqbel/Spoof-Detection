%******************************************************************************************
% DESCRIPTION:
% One epoch of the spoofing-protected navigator (the FSM behind
% spoofMonitor2hz). HOST-OWNED FILTERS, INCREMENT FORM.
%
% MODES (constants in CST_spfMode)
%   NOMINAL    host KF updated normally; monitors watch its increments;
%              clean-closed windows refresh the anchor; any alarm => LATCH.
%   COAST      host does NOT update its operational KF (100 Hz extrapolation
%              makes it the INS-only coast); the update is still computed
%              on a scratch copy to feed y, S here; re-validation counts
%              consecutive passes; full dwell => PROBATION.
%              Monitor silence NEVER exits this mode (capture).
%   PROBATION  host updates a TRIAL copy; monitors watch the trial; alarm
%              => VETO; full quiet probation => COMMIT (handback).
%
% HOST CONTRACT (see docs/HOST_2HZ_WIRING.m)
%   - every epoch: kfMeas from the ACTIVE filter's update (operational KF
%     in NOMINAL/COAST, trial in PROBATION)
%   - nav.applyCorrection true (LATCH, COMMIT): treat
%         x+ = KF.states + nav.correction,  P+ = nav.covar
%     as this epoch's update result of the OPERATIONAL KF and run the usual
%     feedback bookkeeping on it
%   - kfCommand.reseedKF true (probation opens): trial <- copy of the
%     operational KF (covariance kfCommand.reseedCov)
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
function [sys, spoofTel] = protectedNav(sys, kfMeas, propTel, epoch)

info      = STRUCT_SPF.zeroInfo;
kfCommand = STRUCT_SPF.zeroCommand;
numStates = CST_gnssHybrid.NO_STATES;
applyCorrection = false;
correction      = zeros(numStates, 1);

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

        % ---- 2. monitor bank on the host increments ----
        [sys.pool, report] = monitorPool(sys.pool, kfMeas, propTel);

        % ---- 2b. keep the anchor LIVE: (host solution now) - (anchor coast
        %          now) accumulates this epoch's increment; P_C propagates ----
        if (sys.anchor.valid)
            sys.anchor.separation = Phi * sys.anchor.separation + kfIncrement;
            sys.anchor.covariance = Phi * sys.anchor.covariance * Phi' + Q;
            sys.anchor.covariance = (sys.anchor.covariance + sys.anchor.covariance') / 2;
        end % ELSE is trivial

        info.ssAlarm            = report.ssAlarm;
        info.cpiAlarm           = report.cpiAlarm;
        info.alarmPerAxis       = report.alarmPerAxis;
        info.maxProtectionLevel = report.maxProtectionLevel;

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

            if (sys.anchor.valid) && (anchorAge <= CST_spfParam.MAX_ANCHOR_AGE)
                % fall back to the certified-clean coast: the host must
                % move its solution by -(accumulated increments since the
                % anchor window opened) and take the coast covariance
                applyCorrection       = true;
                correction            = -sys.anchor.separation;
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

    % ==================================================================
    case CST_spfMode.COAST
    % ==================================================================

        % ---- 1. protected output = INS-only coast (host extrapolates the
        %         operational KF the same way; GNSS severed) ----
        sys.coastCount = sys.coastCount + 1;
        sys.coastCov   = Phi * sys.coastCov * Phi' + Q;
        sys.coastCov   = (sys.coastCov + sys.coastCov') / 2;

        % ---- 2. test the (untrusted) GNSS against the coast ----
        %         HOST CONTRACT: the operational KF is not updated in COAST
        %         (the scratch update only feeds y, S here), so its prior IS
        %         the coast: x_coast - x_prior = 0 and the residual is y.
        sys.probSep = zeros(numStates, 1);
        coastMinusPrior = zeros(numStates, 1);
        [passed, q_value] = revalidation(kfMeas.innovation, kfMeas.obsMatrix, ...
            kfMeas.measNoiseCov, kfMeas.numMeas, coastMinusPrior, sys.coastCov);
        info.qReval = q_value;
        info.revalComputed = true;

        if (passed)
            sys.dwellCount = sys.dwellCount + 1;
        else
            sys.dwellCount = 0;
        end

        % ---- 3. dwell full: open PROBATION (never direct handback) ----
        if (sys.dwellCount >= CST_spfParam.REVAL_DWELL_REQUIRED)
            info.eventProbationStarted = true;
            sys.probationCount = 0;
            sys.dwellCount     = 0;
            sys.probSep        = zeros(numStates, 1);   % trial starts ON the coast
            sys.pool           = STRUCT_SPF.closeAllWindows(sys.pool);
            sys.mode           = CST_spfMode.PROBATION;
            % order the host to start its TRIAL as a copy of the coasting KF
            kfCommand = STRUCT_SPF.setCommand(true, sys.coastCov);
        end

    % ==================================================================
    case CST_spfMode.PROBATION
    % ==================================================================

        % ---- 1. protected output STAYS on the coast ----
        sys.coastCount = sys.coastCount + 1;
        sys.coastCov   = Phi * sys.coastCov * Phi' + Q;
        sys.coastCov   = (sys.coastCov + sys.coastCov') / 2;

        % ---- 2. diagnostic: GNSS-vs-coast on the trial's innovation ----
        coastMinusPrior = -(Phi * sys.probSep);
        [~, q_value] = revalidation(kfMeas.innovation, kfMeas.obsMatrix, ...
            kfMeas.measNoiseCov, kfMeas.numMeas, coastMinusPrior, sys.coastCov);
        info.qReval = q_value;
        info.revalComputed = true;

        % ---- 3. trial - coast accumulates the trial's increments ----
        sys.probSep = Phi * sys.probSep + kfIncrement;

        % ---- 4. monitor bank on THE TRIAL filter ----
        [sys.pool, report] = monitorPool(sys.pool, kfMeas, propTel);

        info.ssAlarm            = report.ssAlarm;
        info.cpiAlarm           = report.cpiAlarm;
        info.alarmPerAxis       = report.alarmPerAxis;
        info.maxProtectionLevel = report.maxProtectionLevel;
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
                correction         = sys.probSep;          % trial - coast
                sys.coastCov       = kfMeas.postCov;

                % just certified by a full quiet probation: new anchor
                sys.anchor.valid      = true;
                sys.anchor.separation = zeros(numStates, 1);
                sys.anchor.covariance = kfMeas.postCov;
                sys.anchor.epoch      = uint32(epoch);

                sys.probSep    = zeros(numStates, 1);
                sys.coastCount = 0;
                sys.mode       = CST_spfMode.NOMINAL;
                % pool carries over seamlessly
            end
        end
end

% ---- navigation output of this epoch ----
sigmaPosition = zeros(3, 1);
for axisIdx = 1:3
    sigmaPosition(axisIdx) = sqrt(max(sys.coastCov(axisIdx, axisIdx), 0.0));
end
nav = STRUCT_SPF.setNav(applyCorrection, correction, sys.coastCov, sigmaPosition);

info.mode        = sys.mode;
info.dwellCount  = sys.dwellCount;
info.coastEpochs = sys.coastCount;    % time-in-coast (COAST + PROBATION), epochs

spoofTel = STRUCT_SPF.setTel(info, kfCommand, nav);

end

%------------------------------------------------------------------------
