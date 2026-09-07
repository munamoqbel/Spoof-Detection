function [sys, spoofTel] = protected_nav_step(sys, measurement, ...
    obs_matrix, meas_noise_cov, spoofInfo, epoch)
%PROTECTED_NAV_STEP  One epoch of the spoofing-protected navigator (FSM).
%
% MODES (constants in NavMode.m)
%   NOMINAL    KF + monitors; clean-closed windows refresh the anchor;
%              any alarm => LATCH (sticky).
%   COAST      output = INS-only; revalidation_test counts consecutive
%              passes; full dwell => PROBATION.
%              Monitor silence NEVER exits this mode (capture).
%   PROBATION  trial filter takes GNSS under monitor watch; output stays
%              on the coast; alarm => VETO; full quiet probation =>
%              COMMIT (handback).
%
% STRUCTS (defined centrally): sys = NavState, cfg = KujurConfig,
%   info = NavInfo, pool report = PoolReport.
% LEAF FUNCTIONS: kalman_update_step, ins_coast_step,
%   revalidation_test, monitor_pool_step / _open_window / ,._clear.
%
% INPUTS
%   sys              persistent system state (NavState.create)
%   cfg              configuration (KujurConfig.create)
%   measurement      [m x 1]  raw GNSS z at this epoch
%   obs_matrix       [m x n]  H at this epoch
%   meas_noise_cov   [m x m]  V
%   Phi, Q           [n x n]  INTERVAL matrices (Phi_acc, Q_acc from
%                             the 100 Hz side, or the harness constants)
%   epoch            scalar   epoch counter (event stamps, anchor loop)
%
% OUTPUTS
%   sys              updated
%   nav.state            [n x 1]  navigation output of this epoch
%   nav.sigma_position   [3 x 1]  1-sigma of the position states
%   info                 telemetry (see NavInfo.m)

info = STRUCT_SPF.zeroInfo;
kfCommand = STRUCT_SPF.zeroCommand;

switch sys.mode

    % ==================================================================
    case CST_spfMode.NOMINAL
    % ==================================================================

        % ---- 1. normal Kalman filter epoch ----
        [sys.filter.state, sys.filter.covariance, ...
         innovation, innovation_cov, ~] = kalman_update_step( ...
            sys.filter.state, sys.filter.covariance, ...
            measurement, obs_matrix, spoofInfo, meas_noise_cov);

        % ---- 2. monitor bank on the filter ----
        [sys.pool, report] = monitorPool(sys.pool, ...
            innovation, innovation_cov, obs_matrix, ...
            sys.filter.state, sys.filter.covariance, spoofInfo);

        info.ssAlarm            = report.ssAlarm;
        info.cpiAlarm           = report.cpiAlarm;
        info.alarmPerAxis       = report.alarmPerAxis;
        info.maxProtectionLevel = report.maxProtectionLevel;

        % ---- 3. refresh anchor (only on alarm-free epochs, so the
        %         anchor always ends strictly before detection) ----
        if report.cleanCloseFound && ~report.anyAlarm
            sys.anchor.valid      = true;
            sys.anchor.state      = report.cleanCloseState;
            sys.anchor.covariance = report.cleanCloseCovar;
            sys.anchor.epoch      = epoch;
        end

        % ---- 4. latch on any alarm, else open the next window ----
        if (report.anyAlarm)

            info.eventLatched = true;

            if (sys.anchor.valid)
                % fall back to the certified-clean coast, brought
                % forward from its close epoch to 'now' INS-only
                fallback_state = sys.anchor.state;
                fallback_cov   = sys.anchor.covariance;
                for j = sys.anchor.epoch + 1 : epoch
                    [fallback_state, fallback_cov] = ...
                        insCoast(fallback_state, fallback_cov, ...
                            spoofInfo);
                end
                sys.filter.state      = fallback_state;
                sys.filter.covariance = fallback_cov;
                info.eventAnchorEpoch = sys.anchor.epoch;
            else
                % degenerate: alarm before any clean window ever closed
                warning('protected_nav_step:noAnchor', ...
                    ['Alarm at epoch %d before any clean window ', ...
                     'closed; freezing current state.'], epoch);
            end

            sys.pool       = STRUCT_SPF.closeAllWindows(sys.pool);
            sys.dwellCount = 0;
            sys.mode       = CST_spfMode.COAST;

        else
            sys.pool = STRUCT_SPF.openWindow(sys.pool, ...
                sys.filter.state, sys.filter.covariance, ...
                innovation, innovation_cov, obs_matrix);
        end

    % ==================================================================
    case CST_spfMode.COAST
    % ==================================================================

        % ---- 1. output = INS-only propagation (GNSS severed) ----
        [sys.filter.state, sys.filter.covariance] = insCoast( ...
            sys.filter.state, sys.filter.covariance, spoofInfo);

        % ---- 2. test the (untrusted) GNSS against the coast ----
        [passed, q_value] = revalidation(measurement, obs_matrix, ...
            meas_noise_cov, sys.filter.state, sys.filter.covariance);
        info.qReval = q_value;

        if (passed)
            sys.dwellCount = sys.dwellCount + 1;
        else
            sys.dwellCount = 0;
        end

        % ---- 3. dwell full: open PROBATION (never direct handback) ----
        if (sys.dwellCount >= CST_spfParam.REVAL_DWELL_REQUIRED)
            info.eventProbationStarted = true;
            sys.trial.state       = sys.filter.state;
            sys.trial.covariance  = sys.filter.covariance;
            sys.probationCount    = 0;
            sys.dwellCount        = 0;
            sys.pool              = STRUCT_SPF.closeAllWindows(sys.pool);
            sys.mode              = CST_spfMode.PROBATION;
        end

    % ==================================================================
    case CST_spfMode.PROBATION
    % ==================================================================

        % ---- 1. output STAYS on the coast ----
        [sys.filter.state, sys.filter.covariance] = insCoast( ...
            sys.filter.state, sys.filter.covariance, spoofInfo);

        % ---- 2. trial filter takes the GNSS update ----
        [sys.trial.state, sys.trial.covariance, ...
         trial_innovation, trial_innovation_cov, ~] = ...
            kalman_update_step(sys.trial.state, sys.trial.covariance, ...
                measurement, obs_matrix, spoofInfo, meas_noise_cov);

        % ---- 3. diagnostic: GNSS-vs-coast, log only ----
        [~, q_value] = revalidation(measurement, obs_matrix, ...
            meas_noise_cov, sys.filter.state, sys.filter.covariance);
        info.qReval = q_value;

        % ---- 4. monitor bank on THE TRIAL filter ----
        [sys.pool, report] = monitorPool(sys.pool, ...
            trial_innovation, trial_innovation_cov, obs_matrix, ...
            sys.trial.state, sys.trial.covariance, spoofInfo);

        info.ssAlarm            = report.ssAlarm;
        info.cpiAlarm           = report.cpiAlarm;
        info.alarmPerAxis       = report.alarmPerAxis;
        info.maxProtectionLevel = report.maxProtectionLevel;
        % probation windows never refresh the anchor: the trial is not
        % yet trusted, so report.clean_close_* is deliberately ignored.

        if (report.anyAlarm)

            % ---- VETO: monitors caught the trial; discard it ----
            info.eventProbationVetoed = true;
            sys.pool       = STRUCT_SPF.closeAllWindows(sys.pool);
            sys.dwellCount = 0;
            sys.mode       = CST_spfMode.COAST;

        else

            sys.pool = STRUCT_SPF.openWindow(sys.pool, ...
                sys.trial.state, sys.trial.covariance, ...
                trial_innovation, trial_innovation_cov, obs_matrix);

            sys.probationCount = sys.probationCount + 1;

            if sys.probationCount >= CST_spfParam.PROBATION_LENGTH
                % ---- COMMIT (handback) ----
                info.eventHandback    = true;
                sys.filter.state      = sys.trial.state;
                sys.filter.covariance = sys.trial.covariance;

                % just certified by a full quiet probation: new anchor
                sys.anchor.valid      = true;
                sys.anchor.state      = sys.filter.state;
                sys.anchor.covariance = sys.filter.covariance;
                sys.anchor.epoch      = epoch;

                sys.mode = CST_spfMode.NOMINAL;
                % pool carries over seamlessly
            end
        end

end

% ---- navigation output of this epoch ----
nav = STRUCT_SPF.setNav(sys.filter.state, sys.filter.covariance, zeros(3, 1));
for axis_idx = 1:3
    nav.sigmaPosition(axis_idx) = ...
        sqrt(max(sys.filter.covariance(axis_idx, axis_idx), 0.0));
end

info.mode       = sys.mode;
info.dwellCount = sys.dwellCount;

spoofTel = STRUCT_SPF.setTel(info, kfCommand, nav);
end
