%******************************************************************************************
% DESCRIPTION:
%
%
% INPUTS:
%
% OUTPUTS:
%
% ASSUMPTIONS AND LIMITATIONS:
%
%
% REQUIREMENT TRACEABILITY:
%   <traces to requirements (comma-separated list or one trace per line)>
%
%******************************************************************************************
%#codegen
classdef STRUCT_SPF
    methods(Static)
        function [pool] = setMonitorPool(active, hadAlarm, windowAge, ...
                coastState, coastCovariance, innovationBuffer, ...
                innovationCovBuffer, obsMatrixBuffer)

            % Define structure
            pool = struct( ...
                'active',              logical(active), ...
                'hadAlarm',            logical(hadAlarm), ...
                'windowAge',           windowAge, ...
                'coastState',          coastState, ...
                'coastCovariance',     coastCovariance, ...
                'innovationBuffer',    innovationBuffer, ...
                'innovationCovBuffer', innovationCovBuffer,...
                'obsMatrixBuffer',     obsMatrixBuffer);
        end

        function [pool] = zeroMonitorPool(num_meas)

            % Define variables
            windowLength = CST_spfParam.WINDOW_LENGTH;
            numStates = CST_gnssHybrid.NO_STATES;

            active   = false(windowLength, 1);
            hadAlarm = false(windowLength, 1);
            windowAge = zeros(windowLength, 1);

            coastState      = zeros(numStates, windowLength);
            coastCovariance = zeros(numStates, numStates, windowLength);

            innovationBuffer = zeros(num_meas, windowLength, windowLength);
            innovationCovBuffer = zeros(num_meas, num_meas, windowLength, ...
                windowLength);
            obsMatrixBuffer = zeros(num_meas, numStates, windowLength, ...
                windowLength);

            % Set the output
            pool = STRUCT_SPF.setMonitorPool(active, hadAlarm, windowAge, ...
                coastState, coastCovariance, innovationBuffer, ...
                innovationCovBuffer, obsMatrixBuffer);
        end

        function [ssMonitor] = setMonitorReport(alarmPerAxis, anyAlarm, ...
                maxProtectionLevel, ssAlarm, cpiAlarm, cleanCloseFound, ...
                cleanCloseState, cleanCloseCovar)

            % Define structure
            ssMonitor = struct( ...
                'alarmPerAxis',       logical(alarmPerAxis), ...
                'anyAlarm',           logical(anyAlarm), ...
                'maxProtectionLevel', maxProtectionLevel, ...
                'ssAlarm',            logical(ssAlarm), ...
                'cpiAlarm',           logical(cpiAlarm), ...
                'cleanCloseFound',    logical(cleanCloseFound), ...
                'cleanCloseState',    cleanCloseState, ...
                'cleanCloseCovar',    cleanCloseCovar);
        end

        function [ssMonitor] = zeroMonitorReport

            % Init values
            axisAlarm = zeros(1, 3, 'logical');
            alarm = false;
            scalar = 0.0;
            zeroState = zeros(CST_gnssHybrid.NO_STATES, 1);
            zeroCovar = zeros(CST_gnssHybrid.NO_STATES, CST_gnssHybrid.NO_STATES);

            ssMonitor = STRUCT_SPF.setMonitorReport(axisAlarm, alarm, ...
                scalar, alarm, alarm, alarm, zeroState, zeroCovar);
        end

        function [poolOut] = closeAllWindows(poolIn)

            % Set the output
            poolOut = poolIn;
            poolOut.active(:)    = false;
            poolOut.hadAlarm(:)  = false;
            poolOut.windowAge(:) = 0;

        end

        function [poolOut] = openWindow(poolIn, kf_state, kf_covariance, ...
                innovation, innovation_cov, obs_matrix)

            % Set the output
            poolOut = poolIn;

            % ---- find the first free slot ----
            slotFound = false;
            freeSlot = 1;
            windowLength = CST_spfParam.WINDOW_LENGTH;

            for wIndx = 1:windowLength
                if (~slotFound) && (~poolIn.active(wIndx))
                    freeSlot = wIndx;
                    slotFound = true;
                end
            end

            % ---- initialise the window ----
            if (slotFound)
                poolOut.active(freeSlot)    = true;
                poolOut.hadAlarm(freeSlot)  = false;
                poolOut.windowAge(freeSlot) = 1;

                poolOut.coastState(:, freeSlot)         = kf_state;         % E24
                poolOut.coastCovariance(:, :, freeSlot) = kf_covariance;

                poolOut.innovationBuffer(:, 1, freeSlot)       = innovation;
                poolOut.innovationCovBuffer(:, :, 1, freeSlot) = innovation_cov;
                poolOut.obsMatrixBuffer(:, :, 1, freeSlot)     = obs_matrix;
            end

        end
        function [ssMonitor] = setSSmonitor(alarmPerAxis, anyAlarm, ...
                separation, sigmaSeparation, protectionLevel, ...
                maxProtectionLevel)

            % Define structure
            ssMonitor = struct( ...
                'alarmPerAxis',       logical(alarmPerAxis), ...
                'anyAlarm',           logical(anyAlarm), ...
                'separation',         separation, ...
                'sigmaSeparation',    sigmaSeparation, ...
                'protectionLevel',    protectionLevel, ...
                'maxProtectionLevel', maxProtectionLevel);
        end

        function [ssMonitor] = zeroSSmonitor

            % Init values
            axisAlarm = zeros(1, 3, 'logical');
            alarm = false;
            zeroVector = zeros(1, 3);
            scalar = 0.0;

            ssMonitor = STRUCT_SPF.setSSmonitor(axisAlarm, alarm, ...
                zeroVector, zeroVector, zeroVector, scalar);
        end

        function [filter] = setFilter(state, covariance)

            % Define structure
            filter = struct( ...
                'state',      state, ...
                'covariance', covariance);
        end

        function [filter] = zeroFilter

            % Init values
            state = zeros(CST_gnssHybrid.NO_STATES, 1);
            covariance = zeros(CST_gnssHybrid.NO_STATES, CST_gnssHybrid.NO_STATES);

            filter = STRUCT_SPF.setFilter(state, covariance);

        end

        function [trial] = setTrial(state, covariance)

            % Define structure
            trial = struct( ...
                'state',      state, ...
                'covariance', covariance);
        end

        function [trial] = zeroTrial

            % Init values
            state = zeros(CST_gnssHybrid.NO_STATES, 1);
            covariance = zeros(CST_gnssHybrid.NO_STATES, CST_gnssHybrid.NO_STATES);

            trial = STRUCT_SPF.setTrial(state, covariance);

        end

        function [anchor] = setAnchor(valid, state, covariance, epoch)

            % Define structure
            anchor = struct( ...
                'valid',      valid, ...
                'state',      state, ...
                'covariance', covariance, ...
                'epoch',      epoch);
        end

        function [anchor] = zeroAnchor

            % Init values
            valid = false;
            state = zeros(CST_gnssHybrid.NO_STATES, 1);
            covariance = zeros(CST_gnssHybrid.NO_STATES, CST_gnssHybrid.NO_STATES);
            epoch = 0.0;

            anchor = STRUCT_SPF.setAnchor(valid, state, covariance, epoch);

        end

        function [sys] = setSys(mode, filter, trial, pool, anchor, ...
                dwellCount, probationCount)

            % Define structure
            sys = struct( ...
                'mode',           mode, ...
                'filter',         filter, ...
                'trial',          trial, ...
                'pool',           pool, ...
                'anchor',         anchor, ...
                'dwellCount',     dwellCount, ...
                'probationCount', probationCount);
        end

        function [sys] = zeroSys(numMeas)

            % Init values
            mode = CST_spfMode.NOMINAL;
            filter = STRUCT_SPF.zeroFilter;
            trial = STRUCT_SPF.zeroTrial;
            pool = STRUCT_SPF.zeroMonitorPool(numMeas);
            anchor = STRUCT_SPF.zeroAnchor;
            dwellCount = 0.0;
            probationCount = 0.0;

            sys = STRUCT_SPF.setSys(mode, filter, trial, pool, ...
                anchor, dwellCount, probationCount);
        end

        function [altXCheck] = setAltXCheck(suspect, altDiff)

            % Define structure
            altXCheck = struct( ...
                'suspect', suspect, ...
                'altDiff', altDiff);

        end

        function [altXCheck] = zeroAltXCheck

            % Init values
            suspect = false;
            altDiff = 0.0;

            altXCheck = STRUCT_SPF.setAltXCheck(suspect, altDiff);
        end

        function [spoofTel] = setTel(info, kfCommand, nav)

            % Define structure
            spoofTel = struct( ...
                'info',      info, ...
                'kfCommand', kfCommand, ...
                'nav',       nav);
        end

        function [spoofTel] = zeroTel

            % Init values
            info = STRUCT_SPF.zeroInfo;
            kfCommand = STRUCT_SPF.zeroCommand;
            nav = STRUCT_SPF.zeroNav;

            % Define structure
            spoofTel = STRUCT_SPF.setTel(info, kfCommand, nav);
        end

        function [info] = setInfo(mode, ssAlarm, cpiAlarm, ...
                alarmPerAxis, maxProtectionLevel, qReval, revalComputed,...
                dwellCount, eventLatched, eventAnchorEpoch, eventProbationStarted,...
                eventProbationVetoed, eventHandback, anchorMissing)

            % Define structure
            info = struct( ...
                'mode',                  mode, ...
                'ssAlarm',               logical(ssAlarm), ...
                'cpiAlarm',              logical(cpiAlarm), ...
                'alarmPerAxis',          logical(alarmPerAxis), ...
                'maxProtectionLevel',    maxProtectionLevel, ...
                'qReval',                qReval, ...
                'revalComputed',         logical(revalComputed), ...
                'dwellCount',            dwellCount, ...
                'eventLatched',          logical(eventLatched), ...
                'eventAnchorEpoch',      eventAnchorEpoch, ...
                'eventProbationStarted', logical(eventProbationStarted), ...
                'eventProbationVetoed',  logical(eventProbationVetoed), ...
                'eventHandback',         logical(eventHandback), ...
                'anchorMissing',         logical(anchorMissing));
        end

        function [info] = zeroInfo

            % Init values
            mode                  = 1;
            ssAlarm               = false;
            cpiAlarm              = false;
            alarmPerAxis          = false(1, 3);
            maxProtectionLevel    = 0.0;
            qReval                = 0.0;
            revalComputed         = false;
            dwellCount            = 0;
            eventLatched          = false;
            eventAnchorEpoch      = 0;
            eventProbationStarted = false;
            eventProbationVetoed  = false;
            eventHandback         = false;
            anchorMissing         = false;

            % Define structure
            info = STRUCT_SPF.setInfo(mode, ssAlarm, cpiAlarm, ...
                alarmPerAxis, maxProtectionLevel, qReval, revalComputed,...
                dwellCount, eventLatched, eventAnchorEpoch, eventProbationStarted,...
                eventProbationVetoed, eventHandback, anchorMissing);
        end

        function [command] = setCommand(reseedKF, reseedState, reseedCov)

            % Define structure
            command = struct( ...
                'reseedKF',    logical(reseedKF), ...
                'reseedState', reseedState, ...
                'reseedCov',   reseedCov);
        end

        function [command] = zeroCommand

            % Init values
            reseedKF = false;
            reseedState = zeros(CST_gnssHybrid.NO_STATES, 1);
            reseedCov = zeros(CST_gnssHybrid.NO_STATES, CST_gnssHybrid.NO_STATES);

            % Define structure
            command = STRUCT_SPF.setCommand(reseedKF, reseedState, reseedCov);
        end

        function [nav] = setNav(state, covar, sigmaPosition)

            % Define structure
            nav = struct( ...
                'state',         state, ...
                'covar',         covar, ...
                'sigmaPosition', sigmaPosition);
        end

        function [nav] = zeroNav

            % Init values
            state = zeros(CST_gnssHybrid.NO_STATES, 1);
            covar = eye(CST_gnssHybrid.NO_STATES);
            sigmaPosition = zeros(3, 1);

            % Define structure
            nav = STRUCT_SPF.setNav(state, covar, sigmaPosition);
        end

    end
end

%------------------------------------------------------------------------------------------
