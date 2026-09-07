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

% !!! TRANSCRIPTION CONTINUES - lines 206 to end (zeroAnchor body, setSys,
% zeroInfo, setCommand, zeroCommand, setNav, setTel, ...) not yet provided.
